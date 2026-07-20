
overall_folder = 'data/';

load('data/11-06-2025_190150_Gimbal.mat')
load('data/2025_06_11_18_59_42_Alta_data_struc_10_59_41.mat')
focal_length = 35; %mm

Alta_data_struc.DTime.TimeZone = 'UTC';
Alta_Dtime = Alta_data_struc.DTime;
IMU_Dtime = IMU_struc.DTime;
flen_mm = focal_length;
pixp_microns = 3.45;

load('tidal_data_june2025.mat') 
dtime_future.TimeZone = 'UTC';
buoy_ellh = interp1(dtime_future,z,Alta_Dtime);
alt_ellh_diff = Alta_data_struc.vehicle_gps_position.alt_ellipsoid - buoy_ellh;

alt = interp1(Alta_Dtime,alt_ellh_diff,IMU_Dtime);
freeboard = alt - 0.3;

s = load('dolp_theta_vecs.mat');
DOLP_vec = s.DOLP_full;
theta_vec = s.theta_full;
ind_max = find(DOLP_vec==max(DOLP_vec),1,'first');
DOLP_full = linspace(0,1,10000)';
theta_full = interp1(DOLP_vec(1:ind_max),theta_vec(1:ind_max),DOLP_full,'pchip');
fprintf('preparation for slope fields completed.')

%% Calculating Slope Fields, Georectification
magnetic_north = -14.4;


% 120 meter altitude
frame_start = 1;
frame_end = frame_start+89;
frame_offset = 14699;

mean_heading = mean(IMU_struc.heading(frame_start:frame_end));
% pre-allocate
m_per_px = zeros(frame_end-frame_start+1,1);
frame_extrema_SNWE = zeros(frame_end-frame_start+1,4);
frame_struc_Sx_transformed = repmat(struct(),frame_end-frame_start+1,1)';
frame_struc_Sy_transformed = repmat(struct(),frame_end-frame_start+1,1)';
% obtain mean DOLP to find an empirical gain
counter = 0;
sum_DOLP = 0;
for frame_num = frame_start:frame_end
    counter = counter+1;
%frame_num
    frame_raw = frame_raw_array(:,:,frame_num);

    % [~,S1,S2] = Compute_StokesVecs_by_Conv_Demodul(double(frame_raw),'4x4');

    [~,S1,S2] = Compute_StokesVecs_by_KernelAveraging(frame_raw, '4x4');

    DOLP = sqrt(S1.^2+S2.^2);
    ORI = 0.5*atan2(S2,S1)*180/pi;
    sum_DOLP = DOLP + sum_DOLP;
end

% Introduce Empirical Gain
incidence_angle_vec = IMU_struc.pitch(frame_start+frame_offset:frame_end+frame_offset)+90;
incidence_mean = mean(incidence_angle_vec);
avg_DOLP = sum_DOLP/counter;
mean_DOLP = mean(avg_DOLP,2);
[aov_h,~] = get_aov(2048,2448,pixp_microns,focal_length);
theta = linspace(incidence_mean-aov_h/2,incidence_mean+aov_h/2,length(S1(:,1)));
emp_gain = median(DOLP_full(theta_full> theta(1) & theta_full < theta(end)))/median(mean_DOLP);
fprintf([' emp gain is ' num2str(emp_gain)])


counter = 0;
%disp([file_folder ' and ' file_prefix ' and ' num2str(frame_num)]
for frame_num = frame_start:frame_end

    counter = counter+1;
%disp([file_folder ' and ' file_prefix ' and ' num2str(frame_num)])

    frame_raw = read_pyxis_raw_imagingsource(file_folder,file_prefix,frame_num);

    % [~,S1,S2] = Compute_StokesVecs_by_Conv_Demodul(double(frame_raw),'4x4');

    [~,S1,S2] = Compute_StokesVecs_by_KernelAveraging(frame_raw, '4x4');

    DOLP = sqrt(S1.^2+S2.^2);
    ORI = 0.5*atan2(S2,S1)*180/pi;

    DOLP_int = floor(emp_gain*DOLP*10000);
    DOLP_int(DOLP_int<1) = 1;
    DOLP_int(DOLP_int>10000) = 10000;
    AOI = theta_full(DOLP_int);

    Sx = sind(ORI).*tand(AOI);
    Sy = cosd(ORI).*tand(AOI);

    Sx = Sx - mean(Sx,'all','omitnan');
    Sy = Sy - mean(Sy,'all','omitnan');

    [aov_h,~] = get_aov(2048,2448,pixp_microns,focal_length);

    pitch = IMU_struc.pitch(frame_num+frame_offset);
    roll = IMU_struc.roll(frame_num+frame_offset);
    heading = IMU_struc.heading(frame_num+frame_offset) -mean_heading;
    
    [Sx_out,m_per_px(counter),frame_extrema_SN_WE] = rectifier_deluxe(Sx,aov_h,freeboard(frame_num+frame_offset),pitch,roll,heading);

    frame_extrema_SNWE(counter,:) = [frame_extrema_SN_WE(1,:) frame_extrema_SN_WE(2,:)];

    [Sy_out,~,~] = rectifier_deluxe(Sy,aov_h,freeboard(frame_num+frame_offset),pitch,roll,heading);

    frame_struc_Sx_transformed(counter).Sx = single(Sx_out);
    frame_struc_Sy_transformed(counter).Sy = single(Sy_out);

    % imagesc(sqrt(Sy_out.^2+Sx_out.^2))
    % shading('flat')
    % colormap('gray')
    % pause(0.0001)


end

m_per_px_mean = mean(m_per_px);
disp('slope fields calculated, m_per_px_mean found.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero Padding of the initial rectified images 
WE_mean = mean(frame_extrema_SNWE(:,3:4),'all');
SN_mean = mean(frame_extrema_SNWE(:,1:2),'all');

subframes_modified = length(frame_struc_Sy_transformed);
for k = 1:subframes_modified
    [vert_pixels(k), horiz_pixels(k)] = size(frame_struc_Sy_transformed(k).Sy);
end

ref_S = min(frame_extrema_SNWE(:,1));
ref_W = min(frame_extrema_SNWE(:,3));
offset_S = floor((frame_extrema_SNWE(:,1)-ref_S)./m_per_px)+1;
offset_W = (floor((frame_extrema_SNWE(:,3)-ref_W)./m_per_px))+1;

num_cols = max(horiz_pixels) + max(offset_W);
num_rows = max(vert_pixels) + max(offset_S);

slopefield_stack_x = single(zeros(num_rows,num_cols,subframes_modified));
slopefield_stack_y = slopefield_stack_x;

for j = 1:subframes_modified

    % Xcomp_transformed = reprojected_frame_struc(j).Sx_reprojected;
    % Ycomp_transformed = reprojected_frame_struc(j).Sy_reprojected;
    Xcomp_transformed = frame_struc_Sx_transformed(j).Sx;
    Ycomp_transformed = frame_struc_Sy_transformed(j).Sy;

    Xcomp_padded = padarray(Xcomp_transformed,[offset_S(j) offset_W(j)],NaN,'pre');
    [s1,s2] = size(Xcomp_padded);
    Xcomp_padded = padarray(Xcomp_padded,[num_rows-s1 num_cols-s2],NaN,'post');

    Ycomp_padded = padarray(Ycomp_transformed,[offset_S(j) offset_W(j)],NaN,'pre');
    Ycomp_padded = padarray(Ycomp_padded,[num_rows-s1 num_cols-s2],NaN,'post');

    slopefield_stack_x(:,:,j) = single(Xcomp_padded);
    slopefield_stack_y(:,:,j) = single(Ycomp_padded);

end

[s1,s2] = size(slopefield_stack_y(:,:,1));
N_spect = 2^floor(log(min([s1 s2]))/log(2));
center_rows = floor(s1/2)-N_spect/2:floor(s1/2)+N_spect/2-1;
center_cols = center_rows;
slopefield_stack_x_test = double(slopefield_stack_x(center_rows,center_cols,:));
slopefield_stack_y_test = double(slopefield_stack_y(center_rows,center_cols,:));
size(slopefield_stack_y_test)
slopestack_x_test_name = ['slopestack_x_test' flight_day '_Rec' num2str(acq_num) '_frame'  num2str(frame_start) '.mat'];
slopestack_y_test_name = ['slopestack_y_test' flight_day '_Rec' num2str(acq_num) '_frame'  num2str(frame_start) '.mat'];

save(slopestack_x_test_name, "slopefield_stack_x_test", "-v7.3")
save(slopestack_y_test_name, "slopefield_stack_y_test", "-v7.3")

%% Obtain initial spectrum 
%slopefield_stack_x_test(isnan(slopefield_stack_x_test)) = 0;
%slopefield_stack_y_test(isnan(slopefield_stack_y_test)) = 0;

%fps = 30;
%w = circular_tukey(0*slopefield_stack_x_test+1,0.2);
%[dirspect,~] = compute_slope_spectrum(w.*slopefield_stack_x_test,w.*slopefield_stack_y_test,m_per_px_mean,fps,N_spect);

%Skf_initial = dirspect.Skf;
%nframes = size(slopefield_stack_y_test,3);
%df = fps/nframes; %fs/Nt
dk =2*pi/(m_per_px_mean*length(center_rows)); % 2*k_nyq/2048;
%spectra_struc = compute_sub_spectra(Skf_initial,dk,df,mean_heading + magnetic_north,6);

%figure(1);clf;
%loglog(spectra_struc.k,spectra_struc.F_k,'-','linewidth',2)
%hold on
%loglog(spectra_struc.k,spectra_struc.k.^-2.*spectra_struc.S_k,'-','linewidth',2)

%S_k = spectra_struc.S_k;
% Compare F_f and F_f obtained by F_k and dispersion
%k = spectra_struc.k;
%f = spectra_struc.f;
%F_k = spectra_struc.F_k;
%F_f = spectra_struc.F_f;
%f_disp = sqrt(9.81*k)/(2*pi);
%c = 2*pi*f_disp./k;
%cg = c/2;
%F_f_disp = F_k.*k./cg/(2*pi);
%figure(101);clf;
%loglog(f,F_f,'-',f_disp,F_f_disp,'-','linewidth',2)

%disp('initial spectrum calculated.')
%% Find peak wavenumber
%k_mean = (sum(k.*F_k,'all','omitnan')*dk)./(sum(F_k,'all','omitnan')*dk);
%peak_val = findpeaks(F_k,k);
%k_peak = k(F_k == peak_val(2));
%disp(['peak wavenumber = ' num2str(k_peak)])

disp('filtering preparation completed.')
%% Filtering approach main (Section: Automated)

min_kvec = [4*dk;6*dk;8*dk;12*dk];
t = (0:size(slopefield_stack_y_test,3)-2)./30;
Cpeak = nan(length(t),length(min_kvec));
Thetapeak = nan(length(t),length(min_kvec));

for k_inds = 1:length(min_kvec)

angle = 0:5:355;
dkmin = min_kvec(k_inds);
dkmax = dkmin + 2;


mean_speed_vec = zeros(72,1);
total_shift_stack = zeros(length(t),length(angle));

for loop = 1:length(angle)

[total_shift] = find_speed_from_a_slice(slopefield_stack_y_test,dk,angle(loop),dkmin,dkmax);
total_shift_stack(:,loop) = total_shift;

end

speed_at_every_direction = 30*m_per_px_mean*(total_shift_stack);

figure(k_inds);clf;
imagesc(angle,t,speed_at_every_direction);view([0 -90])
xticks(0:10:360)
xlabel('\theta')
ylabel('t [s]')
title('{\Delta}U(k,\theta,t) at 7 < k < 9 rad/m')
colorbar

fitted_smoothed_speed = NaN*speed_at_every_direction;
% without difference
Cp = nan(length(t),1);
Theta_p = nan(length(t),1);
for i=1:length(t)
f = fit(angle',inpaint_nans(speed_at_every_direction(i,:)'),'sin1');
%Cp = f.a1;
%Theta_p = -f.c1/f.b1;
fitted_smoothed_speed(i,:) = f(angle);
Cp(i) = max(fitted_smoothed_speed(i,:));
Theta_p_vals = angle(fitted_smoothed_speed(i,:)==Cp(i));
Theta_p(i) = mean(Theta_p_vals);
end

figure(k_inds+10);clf;
imagesc(angle,t,fitted_smoothed_speed);view([0 -90])
xticks(0:10:360)
xlabel('\theta')
ylabel('t [s]')
title('{\Delta}U(k,\theta,t) at 7 < k < 9 rad/m')
colorbar

% store the Cpeak values and their corresponding angles as 
Cpeak(:,k_inds) =Cp;
Thetapeak(:,k_inds) = Theta_p;
end

peak_names = [flight_day '_Rec' num2str(acq_num) '_frame'  num2str(frame_start) '.mat'];

save(['cpeak' peak_names], "Cpeak", "-v7.3")
save(['theta_peak' peak_names], "Thetapeak", "-v7.3")

fprintf('filtering section completed.')
%%
Cpeak_east = Cpeak.*sind(Thetapeak);
Cpeak_north = Cpeak.*cosd(Thetapeak);

t = (0:size(slopefield_stack_y_test,3)-2)./30;
figure(30);clf;
imagesc(min_kvec,t,Cpeak_north);view([0 -90])
xlabel('k [rad/m')
ylabel('t [s]')
title('Cpeak north')
colorbar


figure(4);clf;
imagesc(min_kvec,t,Cpeak_east);view([0 -90])
xlabel('k [rad/m')
ylabel('t [s]')
title('Cpeak East')
colorbar
% FROM HERE
% Peak speet thru North or East should always be similar except the times when jitter occurs
% In order to find out when and how much of jitter happened, you should detrend Cpeak North and East
% These detrended values should be similar along all wavenumbers, so you can average across wavenumbers
detrended_Cpeak_north = Cpeak_north - movmean(Cpeak_north,60);
detrended_Cpeak_east = Cpeak_east - movmean(Cpeak_east,60);
figure(66);clf;
for l=1:length(min_kvec)
plot(t,detrended_Cpeak_north(:,l),'-','LineWidth',2)
hold on 
end
hold on
d_jitter_north = 1.0*mean(detrended_Cpeak_north,2); % m/s
d_jitter_east = 1.0*mean(detrended_Cpeak_east,2); % m/s
plot(t,d_jitter_north,'-k','LineWidth',2) 
title('Detrended Cpeak north and jitter')
saveas(gcf,'Detrended_Cpeak_north_and_jitter_June11.png')

figure(77);clf;
imagesc(min_kvec,t,1./d_jitter_north.*detrended_Cpeak_north);view([0 -90])
xlabel('k [rad/m')
ylabel('t [s]')
title('{\Delta}U_{jitter} at 7 < k < 9 rad/m')
colorbar

dy_net = 1/30*[0; cumsum(d_jitter_north)];
dx_net = 1/30*[0; cumsum(d_jitter_east)];
d_jitter_north = [0;d_jitter_north];
d_jitter_east = [0;d_jitter_east];
disp('d_jitter values found.')
%%
% if no subpixel consideration
frame_extrema_SNWE_for_padded = repmat([0 size(slopefield_stack_y,1) 0 size(slopefield_stack_y,2)],size(slopefield_stack_y,3),1);
dy_pix_shift = round(dy_net./m_per_px_mean);
dx_pix_shift = round(dx_net./m_per_px_mean);

frame_extrema_SNWE_for_padded(:,1:2) = frame_extrema_SNWE_for_padded(:,1:2) + dy_pix_shift;
frame_extrema_SNWE_for_padded(:,3:4) = frame_extrema_SNWE_for_padded(:,3:4)  - dx_pix_shift;
 
ref_S = min(frame_extrema_SNWE_for_padded(:,1));
ref_W = min(frame_extrema_SNWE_for_padded(:,3));
offset_S = floor((frame_extrema_SNWE_for_padded(:,1)-ref_S))+1;
offset_W = floor((frame_extrema_SNWE_for_padded(:,3)-ref_W))+1;

num_cols = size(slopefield_stack_y,2) + max(offset_W);
num_rows = size(slopefield_stack_y,1) + max(offset_S);

subframes_modified = size(slopefield_stack_y,3);
slopefield_stack_x_new = single(zeros(num_rows,num_cols,subframes_modified));
slopefield_stack_y_new = slopefield_stack_x_new;

for j = 1:subframes_modified

    % Xcomp_transformed = reprojected_frame_struc(j).Sx_reprojected;
    % Ycomp_transformed = reprojected_frame_struc(j).Sy_reprojected;
    Xcomp_transformed = slopefield_stack_x(:,:,j);
    Ycomp_transformed = slopefield_stack_y(:,:,j);

    Xcomp_padded = padarray(Xcomp_transformed,[offset_S(j) offset_W(j)],NaN,'pre');
    [s1,s2] = size(Xcomp_padded);
    Xcomp_padded = padarray(Xcomp_padded,[num_rows-s1 num_cols-s2],NaN,'post');

    Ycomp_padded = padarray(Ycomp_transformed,[offset_S(j) offset_W(j)],NaN,'pre');
    Ycomp_padded = padarray(Ycomp_padded,[num_rows-s1 num_cols-s2],NaN,'post');

    slopefield_stack_x_new(:,:,j) = single(Xcomp_padded);
    slopefield_stack_y_new(:,:,j) = single(Ycomp_padded);

end

flight_day(regexp(flight_day,'-'))=[];
slope_stack_y_name = ['slope_stack_y_' flight_day '_Rec' num2str(acq_num) '_frame'  num2str(frame_start) '.mat'];
save(slope_stack_y_name, "slopefield_stack_y_new", "-v7.3")

slope_stack_x_name = ['slope_stack_x_' flight_day '_Rec' num2str(acq_num) '_frame'  num2str(frame_start) '.mat'];
save(slope_stack_x_name, "slopefield_stack_x_new", "-v7.3")

m_per_px_name = ['m_per_px_' flight_day '_Rec' num2str(acq_num) '_frame'  num2str(frame_start) '.mat'];
save(m_per_px_name, "m_per_px", "-v7.3") 

fprintf('new slopefield stacks are saved.')
  
disp('script ended.')
