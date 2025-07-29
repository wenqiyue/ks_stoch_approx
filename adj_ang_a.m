function [th_a] = adj_ang_a(th,dim)
%ADJ_ANG_A Summary of this function goes here
%   Detailed explanation goes here

[a,b] = size(th);

if nargin<2
    if isvector(th)
        if a == 1
            dim = 2;
        else
            dim = 1;
        end
    else
        dim = 1;
    end
end

if dim == 1
    th1 = th(1,:);
    dth = diff(th,[],1);
    dth_a = mod(dth+pi,2*pi)-pi;
    th_a = th1 + cumsum(cat(1,zeros(1,b),dth_a),1);
elseif dim == 2
    th1 = th(:,1);
    dth = diff(th,[],2);
    dth_a = mod(dth+pi,2*pi)-pi;
    th_a = th1 + cumsum(cat(2,zeros(a,1),dth_a),2);
end

end

