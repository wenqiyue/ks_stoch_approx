function [out_prod] = outer_prod(xt,yt,d)
%OUTER_PROD Summary of this function goes here
%   Detailed explanation goes here

s = size(xt);
n = length(s);

%xt_a = permute(shiftdim(xt,-1),[2:(d+1),1,(d+2):(n+1)];
%yt_tr = 

xt_a = permute(xt,[1:d,n+1,(d+1):n]);
yt_tr = permute(yt,[1:(d-1),n+1,d:n]);

out_prod = xt_a.*yt_tr;

end

