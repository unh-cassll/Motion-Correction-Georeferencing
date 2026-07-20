% This function calculates the in-plane rotation [in degrees] between 2
% consecutive wave slope field matrices, taken about the image center. The
% inter-frame translation is removed first (else wave advection is misread as
% rotation), then the second frame is rotated through a range of candidate
% angles and the one best matching the first frame is kept (zero-mean
% normalized cross-correlation over a central window), refined with a parabola.
% It is the rotation counterpart of find_speed_from_a_slice, used by rotfix.
%
% This is area-based (correlation) image registration by a search over the
% rotation parameter (Zitova & Flusser, Image Vision Comput 21:977-1000, 2003).
% The similarity measure is the zero-mean normalized cross-correlation (ZNCC)
% of Lewis (Fast Normalized Cross-Correlation, Vision Interface, 1995, 120-123);
% the zero-mean, normalized form is invariant to a linear intensity offset/gain
% between frames, which is why it is the standard robust criterion in digital
% image correlation (Pan et al., Meas Sci Technol 20:062001, 2009).

% WARNING: theta_max is the +/- search bound in degrees; keep it near the
% expected per-frame jitter rotation (a fraction of a degree at 30 fps).
% Nathan Laxague, 2026

function theta = find_rotation_from_a_slice(frame1,frame2,theta_max)

[s1,s2] = size(frame1);
cx = (s2-1)/2;
cy = (s1-1)/2;
xcol = 0:s2-1;
yrow = 0:s1-1;

% a slope field is much rougher than an intensity image; blur both frames a
% little so the correlation locks onto wave structure, not per-pixel noise
frame1(~isfinite(frame1)) = 0;
frame2(~isfinite(frame2)) = 0;
sigma = 1.5;
g = exp(-(-ceil(3*sigma):ceil(3*sigma)).^2/(2*sigma^2));
g = g/sum(g);
frame1 = conv2(g,g,frame1,'same');
frame2 = conv2(g,g,frame2,'same');

% inter-frame translation from Hanning-windowed phase correlation; folded into
% the search sampling below so frame2 is compared already advection-aligned
wy = 0.5*(1-cos(2*pi*(0:s1-1)'/(s1-1)));
wx = 0.5*(1-cos(2*pi*(0:s2-1)/(s2-1)));
cross = conj(fft2(frame1.*(wy*wx))).*fft2(frame2.*(wy*wx));
r = real(ifft2(cross./(abs(cross)+1e-10)));
[~,idx] = max(r(:));
[py,px] = ind2sub([s1,s2],idx);
dy = py-1; dx = px-1;
if dy > s1/2, dy = dy-s1; end
if dx > s2/2, dx = dx-s2; end

% central window, kept inside the rotation/shift-induced border
mrg = ceil(theta_max*pi/180*max(s1,s2)/2 + hypot(dx,dy) + 4);
rows = mrg+1:s1-mrg;
cols = mrg+1:s2-mrg;
[U,V] = meshgrid(cols-1,rows-1);
a = frame1(rows,cols);
a = a(:);

% coarse-to-fine search for the angle that best rotates frame2 onto frame1
center = 0;
span = theta_max;
best_theta = 0;
for level = 1:3
    if level == 1
        candidates = linspace(-theta_max,theta_max,25);
    else
        candidates = linspace(center-span,center+span,15);
    end
    ccs = -Inf*ones(size(candidates));
    for c = 1:length(candidates)
        th = candidates(c)*pi/180;
        Xq = cos(th)*(U-cx) - sin(th)*(V-cy) + cx + dx;
        Yq = sin(th)*(U-cx) + cos(th)*(V-cy) + cy + dy;
        b = reshape(interp2(xcol,yrow,frame2,Xq,Yq,'linear'),[],1);
        valid = isfinite(a) & isfinite(b);
        if nnz(valid) < 0.5*numel(b)
            continue
        end
        aa = a(valid) - mean(a(valid));
        bb = b(valid) - mean(b(valid));
        denom = sqrt(sum(aa.^2)*sum(bb.^2));
        if denom > 0
            ccs(c) = sum(aa.*bb)/denom;
        end
    end
    [~,ind_best] = max(ccs);
    best_theta = candidates(ind_best);
    % 3-point parabola refine on an interior, concave peak
    if ind_best > 1 && ind_best < length(candidates)
        y1 = ccs(ind_best-1); y2 = ccs(ind_best); y3 = ccs(ind_best+1);
        denom = y1 - 2*y2 + y3;
        if isfinite(denom) && denom < 0
            step = candidates(2) - candidates(1);
            best_theta = best_theta + 0.5*(y1-y3)/denom*step;
        end
    end
    center = best_theta;
    span = span*0.2;
end

theta = max(-theta_max,min(theta_max,best_theta));

end
