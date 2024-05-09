function [v1,v2] = hist_with_gaussian_fit(x_data,zz,lw)
%FIT_GAUSSIAN_PLOT Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    zz = -5:0.01:5;
end

if nargin < 3
    lw = 1.5;
end

[xx,yy] = fit_gaussian_dist(x_data,zz);

ax1 = histogram(x_data,'normalization','pdf','numbins',50);
hold on
ax2 = plot(xx,yy,'r','linewidth',lw);
hold off
xlim([min(xx),max(xx)])

if nargout >= 1
    v1 = ax1;
    if nargout >= 2
        v2 = ax2;
    end
end

end

