function [t,xt,dxt,xt_unadj] = rk_4_mod_2pi(frhs,tt,x0,mod_ind,n)
%RK_4 Summary of this function goes here
%   Detailed explanation goes here


if nargin<5
    n=1;
end

if nargin<4
    mod_ind = 1:length(x0);
end


if nargout < 3
    flg_dxt = false;
else
    flg_dxt = true;
end

if nargout < 4
    flg_xt_unadj = false;
else
    flg_xt_unadj = true;
end

[a,b]=size(x0);

if a==1
    x0=x0';
end

tl=ceil(length(tt)/n);

xt=zeros(tl,length(x0));
if flg_dxt == true
    dxt=zeros(tl,length(x0));
end

if flg_xt_unadj == true
    xt_unadj = zeros(tl,length(x0));
end

xi=x0;

xi_unadj = xi;

xi(mod_ind)=adjust_angles(xi(mod_ind));

i_n=0;

for i=1:length(tt)-1
    
    ti=tt(i);    
    dti=tt(i+1)-tt(i);
    
    dxi=frhs(ti,xi);
    
    if mod(i-1,n)==0
        i_n=i_n+1;
        xt(i_n,:)=xi;
        if flg_dxt == true
            dxt(i_n,:)=dxi;
        end
        if flg_xt_unadj == true
            xt_unadj(i_n,:)=xi_unadj;
        end
    end
       
    f1=frhs(ti,xi);
    f2=frhs(ti+dti/2,xi+f1*dti/2);
    f3=frhs(ti+dti/2,xi+f2*dti/2);
    f4=frhs(ti+dti,xi+f3*dti);
    
    xi=xi+1/6*(f1+2*f2+2*f3+f4)*dti;
    xi(mod_ind) = adjust_angles(xi(mod_ind));
    xi_unadj = xi_unadj + 1/6*(f1+2*f2+2*f3+f4)*dti; 
end

if mod(length(tt)-1,n)==0
    i_n=i_n+1;
    xt(i_n,:)=xi;
    dxi=frhs(ti(end),xi);
    if flg_dxt == true
        dxt(i_n,:)=dxi;
    end
    if flg_xt_unadj == true
        xt_unadj(i_n,:)=xi_unadj;
    end
end

t=tt(1:n:end);

end