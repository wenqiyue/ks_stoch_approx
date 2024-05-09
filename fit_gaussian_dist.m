function [x,y] = fit_gaussian_dist(x_data,zz)
%FIT_DIST_WITH_GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    zz = -5:0.01:5;
end

mu_x = mean(x_data);
var_x = var(x_data);

xx = mu_x + sqrt(var_x)*zz;

yy = 1/sqrt(2*pi*var_x)*exp(-1/2*zz.^2);

x = xx;
y = yy;

end

