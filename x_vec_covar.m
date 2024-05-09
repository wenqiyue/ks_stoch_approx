function value = x_vec_covar(tauNt,xt,t1_ind,flg_out)
%X_MSQ Calculate auto-covariance of a function from its time series
%   x_covar(tauNt,xt)
%   tauNt : array of no. of steps at which covariance is calculated
%   xt: the time series of the function

if nargin == 1
    xt = tauNt;
end

if size(xt,1) < size(xt,2)
    xt = xt';
end

xl = size(xt,1);
n_x = size(xt,2);

if nargin < 3 || t1_ind > ceil(xl/2)
    t1_ind = ceil(xl/2);
end

if nargin < 4
    flg_out = 'm';
end

% removing mean from the function
xt = xt - mean(xt,1);

if nargin == 1
    %tauNt = 0:1:floor(xl/2);
    tauNt = 0:1:99;
end

tauNl = length(tauNt);

% if max(tauNt) > xl - 0.5
%     error('range of tau larger than length of time series')
% end

%t1_ind = ceil(xl/2);

x_vec_covar_t = zeros(tauNl,n_x,n_x);

for tauNi = 1:tauNl
    
    tauN = tauNt(tauNi);
    
    %x_vec_covar_tauNi = mean((xt((1+tauN):(t1_ind+tauN),:) .* xt(1:t1_ind,:)),1);
    x_vec_covar_tauNi = mean(outer_prod(xt(1:t1_ind,:),xt((1+tauN):(t1_ind+tauN),:),2));
    
    x_vec_covar_t(tauNi,:,:) = x_vec_covar_tauNi;
    
end

switch flg_out
    case 'm'
        value = x_vec_covar_t;
    case 'v'
        value = reshape(x_vec_covar_t,[tauNl,n_x^2,1]);
end

end

