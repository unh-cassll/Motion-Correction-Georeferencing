
function [ts_new, y_mean, y_var, y_std,  N_good] = my_stddev(ts, num_of_std_dev,n_times)

y = ts;
y_mean = 0 ;
y_std = 100*max(y); %this initialization is necessary for 1st calculation of the stats and this doesnt count.

%THE FIRST LOOP IS NOT A FILTER; IT CALCULATES THE MEAN, STD etc. for the first time. 
% If you want to filter the data 3 times, repeat the loop one more.
n_times = n_times+1;

for k = 1:n_times

SS = 0; %sum of squares
S = 0;
N_good = 0;


for n = 1:length(y)
    
    if abs(y(n) - y_mean) > abs(num_of_std_dev*y_std)
        y(n) = NaN;
    
    else
        
        if (~isnan(y(n)))
        N_good = N_good +1;
        S = S + y(n);
        SS = SS + y(n) * y(n);
        end
    
    end
end
    
y_mean = S/N_good;
y_var = SS/N_good - y_mean^2;
y_std = sqrt(y_var);

end

ts_new = y;

end
