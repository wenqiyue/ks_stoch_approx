function [value] = kura_saka_alltoall(ti,N,K,l,w,phii)
%KURA_SAKA_FULL Summary of this function goes here
%   Detailed explanation goes here

[a,b] = size(phii);

if (b == 1)
    phii = phii';
end

reipsi = 1/N*sum(exp(1i*phii));
r = abs(reipsi);
psi = angle(reipsi);

f = w + K*r*sin(psi-phii-l);

value = f';

end

