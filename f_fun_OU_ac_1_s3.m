function [value] = f_fun_OU_ac_1_s3(p,tau,flg_var)
%F_FUN_CURVEFIT Summary of this function goes here
%   Detailed explanation goes here

[a,b] = size(tau);

if a == 1
    tau = tau';
    [a,b] = size(tau);
end

gamma = p(1);
omega = p(2);

sigma11 = p(3);
sigma12 = p(4);
sigma22 = p(5);

y11 = 1/(4*gamma*(gamma^2+omega^2)).*exp(-gamma*tau).*...
    ((2*gamma^2*(sigma11^2+sigma12^2)+2*gamma*omega*sigma12*(sigma11+sigma22)+omega^2*(sigma11^2+2*sigma12^2+sigma22^2)).*cos(omega*tau)+...
    gamma*(sigma11+sigma22)*(2*gamma*sigma12+omega*(-sigma11+sigma22))*sin(omega*tau));

y12 = 1/(4*gamma*(gamma^2+omega^2)).*exp(-gamma*tau).*...
    (gamma*(sigma11+sigma22)*(2*gamma*sigma12+omega*(-sigma11+sigma22)).*cos(omega*tau)-...
    (2*gamma^2*(sigma11^2+sigma12^2)+2*gamma*omega*sigma12*(sigma11+sigma22)+omega^2*(sigma11^2+2*sigma12^2+sigma22^2)).*sin(omega*tau));

y21 = 1/(4*gamma*(gamma^2+omega^2)).*exp(-gamma*tau).*....
    (gamma*(sigma11+sigma22)*(2*gamma*sigma12+omega*(-sigma11+sigma22)).*cos(omega*tau)+...
    (2*gamma^2*(sigma12^2+sigma22^2)-2*gamma*omega*sigma12*(sigma11+sigma22)+omega^2*(sigma11^2+2*sigma12^2+sigma22^2)).*sin(omega*tau));

y22 = 1/(4*gamma*(gamma^2+omega^2))*exp(-gamma*tau).*...
   ((2*gamma^2*(sigma12^2+sigma22^2)-2*gamma*omega*sigma12*(sigma11+sigma22)+omega^2*(sigma11^2+2*sigma12^2+sigma22^2)).*cos(omega*tau)-...
   gamma*(sigma11+sigma22)*(2*gamma*sigma12+omega*(-sigma11+sigma22)).*sin(omega*tau));


switch flg_var
    case 's'
        y = y11;
    case 'c'
        y = y22;
    case '3'
        
    case 'a'
        y = [y11,y21,y12,y22];
end


value = y;


% if nargin < 3 || strcmp(flg_var,'y')
%     value = y;
% elseif strcmp(flg_out,'y_diff')
%     
%     yd11 = y_data(:,1,1);
%     yd21 = y_data(:,2,1);
%     yd12 = y_data(:,1,2);
%     yd22 = y_data(:,2,2);
%     
%     yd = [yd11,yd12,yd21,yd22];
%     
%     if 
% 
%     %value = y-yd;
%     value = [y-yd; Lambda*(y(1,:)-yd(1,:))];
%     %value = y(:,4) - yd(:,4);
%     %value = [y(:,4)-yd(:,4); 100*(y(1,4)-yd(1,4))];
% else
%     error('Invalid flg_out input');
% end

end

