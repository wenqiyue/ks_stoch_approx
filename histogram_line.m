function [v1] = histogram_line(x_data,nbins)
%HISTOGRAM_LINE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    nbins = 50;
end

[heights,edges] = histcounts(x_data,'normalization','pdf','numbins',nbins);
edges_midpts = 1/2*(edges(1:end-1)+edges(2:end));

ax1 = plot(edges_midpts,heights,'linewidth',1);

if nargout > 0
    v1 = ax1;
end

end

