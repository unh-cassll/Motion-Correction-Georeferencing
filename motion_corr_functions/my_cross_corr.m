function [xlag, corr] = my_cross_corr(ts1, ts2, delta, tlagmax)

%this function excludes observations exceeding 3 standard deviations from the mean
num_of_std_dev = 3;
n_times = 3;

[x, ~] = my_stddev(ts1, num_of_std_dev,n_times);
[y, ~] = my_stddev(ts2, num_of_std_dev,n_times);
x = ts1;
y = ts2;

N = length(x); %or length(y) --Both vectors must have the same length anyway

K = tlagmax/delta;

corr = zeros(2*K + 1,1);
xlag = zeros(2*K + 1,1);

%positive lag

for k = 0:K

    Ng = 0; 
    Sx = 0; Sxx = 0; 
    Sy = 0; Syy = 0; 
    Sxy = 0;

    for n = 1: N-k
        x1 = x(n);
        y1 = y(n+k);

        if(~isnan(x1) && ~isnan(y1))
            Ng = Ng + 1;
            Sx = Sx + x1;
            Sxx = Sxx + x1*x1;
            Sy = Sy + y1;
            Syy = Syy + y1*y1;
            Sxy = Sxy + x1*y1;
        end

    end
    corr(k + K + 1,1) = (Sxy/Ng - (Sx/Ng)*(Sy/Ng))/((Sxx/Ng - (Sx/Ng)^2)^0.5 * (Syy/Ng - (Sy/Ng)^2)^0.5);
    xlag(k + K + 1,1) = delta * k;
end

%plot(xlag((K+1) : end, 1),corr((K+1) : end, 1));
    

%negative lag

for k = 1:K

    Ng = 0; 
    Sx = 0; Sxx = 0; 
    Sy = 0; Syy = 0; 
    Sxy = 0;

    for n = 1: N-k
        x1 = x(n + k);
        y1 = y(n);

        if(~isnan(x1) && ~isnan(y1))
            Ng = Ng + 1;
            Sx = Sx + x1;
            Sxx = Sxx + x1*x1;
            Sy = Sy + y1;
            Syy = Syy + y1*y1;
            Sxy = Sxy + x1*y1;
        end

    end
    corr(K + 1 - k,1) = (Sxy/Ng - (Sx/Ng)*(Sy/Ng))/((Sxx/Ng - (Sx/Ng)^2)^0.5 * (Syy/Ng - (Sy/Ng)^2)^0.5);
    xlag(K + 1 - k,1) = -delta * k;
end

end