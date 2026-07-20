% This function calculates the total shift [in pixels] between 2 consecutive wave slope
% field matrices. The slopefield is filtered out so that only a certain
% wavenumber band and a certain direction of waves are retained.

% WARNING: The input variable 'angle' must be a multipler of 5 due to the
% angle discretization. 
% Jan, 2025 

function [total_shift] = find_speed_from_a_slice(slopefield_stack_y_test,dk,angle,dkmin,dkmax)

[s1,s2] = size(slopefield_stack_y_test(:,:,1));

% Radius Matrix: create wavenumber radius values
min_dim = min([s1 s2]);
[x,y] = meshgrid((1:s2)*dk,(1:s1)*dk);
x = x - mean(mean(x(1,:)));
y = y - mean(mean(y(:,1)));
r = sqrt(x.^2+y.^2); %discretization gets incorrect for angle_matrix close to center, so start r at 5, even the max value is at a smaller r.


% Angle Matrix: move the origin to the center to find the angles
x = 1:s2; %horizontal
y = 1:s1; %vertical
[X,Y] = meshgrid(x,y);
X = X-mean(X(1,:));
Y = Y - mean(Y(:,1));
angles_matrix = -((atand(X./Y)));
total_angle_matrix  = [360+angles_matrix(1:s2/2,1:s2/2) angles_matrix(1:s2/2,s2/2+1:end); 360-90+rot90(angles_matrix(1:s2/2,1:s2/2),1) 90+rot90(angles_matrix(1:s2/2,s2/2+1:end),-1)];

total_angle_matrix = 5*round(total_angle_matrix./5);

%find all the ypeaks and xpeaks for a given angle(deg) and radii(rad/m):
[ROW,COL] = find(total_angle_matrix == angle & r> dkmin & r<dkmax);
if angle == 0
    [ROW,COL] = find(total_angle_matrix == angle | total_angle_matrix == 360 & r> dkmin & r<dkmax);
end
% Only keep these indices; filter out the rest
slopefield_stack_y_test(isnan(slopefield_stack_y_test)) = 0;
Ay = fftshift(fftn(slopefield_stack_y_test));
Ay_padded = 0*Ay;
Ay_padded(ROW,COL,:) = Ay(ROW,COL,:);

Sy_filt = real(ifftn(ifftshift(Ay_padded)));

% apply cross correlation 
dx = NaN*ones(size(Sy_filt,3)-1,1);
total_shift = dx;


[size1,size2] = size(imrotate(Sy_filt(:,:,1),angle));
center_of_the_center = floor(size1/100)*100/2;

for i = 1:size(Sy_filt,3)-1

    frame1 = imrotate(Sy_filt(:,:,i+1),angle);
    frame2 = imrotate(Sy_filt(:,:,i),angle);

    center_rows = (floor(size1/2)-center_of_the_center:floor(size1/2)+center_of_the_center-1) + 1;
    center_cols = (floor(size2/2)-center_of_the_center:floor(size2/2)+center_of_the_center-1) +1;

    frame1 = frame1(center_rows,center_cols);
    frame2 = frame2(center_rows,center_cols);
   
    m = frame1(:,floor(size(frame1,1)/2));
    n = frame2(:,floor(size(frame1,1)/2));
    [xlag, corr] = my_cross_corr(m, n, 1, 200);  

    if isempty(xlag(corr==max(corr)))
        disp(['no energy at the direction ', num2str(angle)])
        break
    end 
    TF = islocalmax(corr);
    TM = islocalmin(corr);
    xrange = [xlag(TF);xlag(TM)];
    total_shift(i) = max(xrange(xrange == min(abs(xrange)) | xrange == -min(abs(xrange))));
    % dstr_now = datestr(now,'mm/dd/yyyy HH:MM:SS');
    % disp([dstr_now '... DONE WITH PART ' num2str(i) '/' num2str(size(Sy_filt,3)-2)])

end


