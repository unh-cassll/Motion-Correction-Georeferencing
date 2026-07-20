% rotfix: residual in-plane rotation correction, applied on top of the
% de-jitter. Total_Dejitter_Process_Premise2025 removes the common-mode
% TRANSLATION jitter, but a pure translation shift cannot remove the leftover
% rotation about the look axis (the "roll" that appears at higher altitudes).
% This routine estimates the per-pair in-plane rotation from the slope stack,
% keeps only its fast (jitter) part, and re-rectifies each frame with the
% exact rotation homography H = K*R*K^-1 taken about the principal point,
% composed with the de-jitter pixel displacement so both come out in a single
% interpolation.

% Drop-in after the de-jitter (uses its outputs):
%   dejitter_disp = [-dx_net dy_net]/m_per_px_mean;   % [col row] shift, px
%   f_px = flen_mm/(pixp_microns*1e-3*2);             % focal length in grid px
%   stack_fixed = rotfix_dejittered_stack(slopefield_stack_y,30,f_px,dejitter_disp);

% NOTE: the rotation estimate is the limiter. A coarse correlation search
% removes most of the in-plane jitter rotation on real data; low-altitude coarse
% texture and parallax hold it back. A gradient-based angle estimate would tighten
% it further.
% Nathan Laxague, 2026

function [stack_fixed,theta_deg] = rotfix_dejittered_stack(stack,fps,f_px,dejitter_disp)

stack = double(stack);
[s1,s2,T] = size(stack);
cx = (s2-1)/2;
cy = (s1-1)/2;

% per consecutive pair, estimate the in-plane rotation about the principal
% point, then accumulate to a frame-0-relative series by summation
theta_pair = zeros(T,1);
for k = 2:T
    theta_pair(k) = find_rotation_from_a_slice(stack(:,:,k-1),stack(:,:,k),1.0);
end
theta_cum = cumsum(theta_pair);

% keep only the oscillatory (jitter) part: remove a low-order polynomial trend,
% then a raised-cosine spectral high-pass with a 0.5 Hz low cut so the slow
% wave drift is left in place. re-zeroed at the first frame.
n = length(theta_cum);
tt = (0:n-1)';
theta_detrended = theta_cum - polyval(polyfit(tt,theta_cum,2),tt);
freqs = min((0:n-1)',n-(0:n-1)')*(fps/n);
lo = 0.5;
taper = 0.1;
mask = double(freqs >= lo);
edge = freqs >= lo-taper & freqs < lo;
mask(edge) = 0.5*(1 - cos(pi*(freqs(edge)-(lo-taper))/taper));
theta_osc = real(ifft(fft(theta_detrended).*mask));
theta_osc = theta_osc - theta_osc(1);
theta_deg = theta_osc;

% re-rectify: backward-warp each frame through H = T(disp)*K*R*K^-1, i.e.
% rotation about the principal point (by theta) plus the de-jitter shift.
% sampling the input frame once avoids a second interpolation.
K = [f_px 0 cx; 0 f_px cy; 0 0 1];
[U,V] = meshgrid(0:s2-1,0:s1-1);
uvw = [U(:)'; V(:)'; ones(1,s1*s2)];
xcol = 0:s2-1;
yrow = 0:s1-1;
stack_fixed = zeros(s1,s2,T);
for k = 1:T
    th = theta_osc(k)*pi/180;
    R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
    H = K*R/K;
    H(1,3) = H(1,3) + dejitter_disp(k,1);
    H(2,3) = H(2,3) + dejitter_disp(k,2);
    src = H*uvw;
    Xq = reshape(src(1,:)./src(3,:),s1,s2);
    Yq = reshape(src(2,:)./src(3,:),s1,s2);
    stack_fixed(:,:,k) = interp2(xcol,yrow,stack(:,:,k),Xq,Yq,'linear');
end

disp('rotfix complete.')

end
