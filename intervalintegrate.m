function [tt_Delta,eta_xt_Delta] = intervalintegrate(tt,xt,Delta)
%INTERVALINTEGRATE Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    Delta = xt;
    xt = tt;
    tt = 1:length(xt);
end

Delta = round(Delta);

[a,b] = size(xt);

if a == 1
    xt = xt';
    [a,b] = size(xt);
end

if length(tt) == 1
    tt = (1:a)*tt;
end

tl = a;


tl_Delta = floor((tl+0.5)/Delta);


eta_xt = zeros(tl_Delta, b);

for Di = 1:tl_Delta
%     if Delta == 1
%         if Di < tl_Delta
%             eta_Di = xt(Di,:)*(tt(Di+1)-tt(Di));
%         else
%             eta_Di = xt(Di,:)*(tt(Di)-tt(Di-1));
%         end
%     else
%         intg_ind = (Di-1)*Delta+1:Di*Delta;
%         %eta_Di = trapz(tt(intg_ind),xt(intg_ind,:));
%         eta_Di = (xt(intg_ind).*diff([tt((Di-1)*Delta),tt(intg_ind)])
%     end

    intg_ind = (Di-1)*Delta+1:Di*Delta;
    if Di == 1
        eta_Di = sum(xt(intg_ind,:).*[tt((Di-1)*Delta+2)-tt((Di-1)*Delta+1),diff(tt(intg_ind))]');
    else
        eta_Di = sum(xt(intg_ind,:).*diff([tt((Di-1)*Delta),tt(intg_ind)])');
    end
    eta_xt(Di,:) = eta_Di;
end

if nargout > 1
    tt_Delta = tt((1:tl_Delta)*Delta);
    eta_xt_Delta = eta_xt;
else
    tt_Delta = eta_xt;
end



end

