%% Simulate the approximade SDE for the KS system




%% clear variables
clearvars

%% location of setting and data files
%cd('D:\MATLAB\scraps\ks_3')
%fileloc = 'D:\MATLAB\scraps\ks_3\';
fileloc = '.\';

%% Parameter setting

setting_param_sgn = '1_a';
%setting_param_sgn = '10_c_l_0';
%setting_param_sgn = '10_a_l_pi_8';

%setting_N_sgn = '';
%setting_N_sgn = '_2N';
setting_N_sgn = '_4N';

setting_sgn = [setting_param_sgn, setting_N_sgn];

setting_sgn_pa = setting_sgn;
%setting_sgn_pa = '1_a_32N';

%% Loading parameters
%param_1 = load(['D:\MATLAB\scraps\ks_3\params_ks_sim_',setting_sgn,'.mat']);
param_1 = load([fileloc,'params_ks_sim_',setting_sgn,'.mat']);

param = param_1;

N = param.N;
K = param.K;
A = param.A;
l = param.l;
w = param.w;

%% Loading statistics
%data_5 = load(['D:\MATLAB\scraps\ks_3\data_ks_sim_',setting_sgn,'_a','.mat']);
data_5 = load([fileloc,'data_ks_sim_',setting_sgn,'_a','.mat']);

Nl_ind = data_5.Nl_ind;
Nr_ind = data_5.Nr_ind;

om = data_5.omega_c;

mu_s = data_5.mu_S_t;
mu_c = data_5.mu_C_t;

mu_x = [mu_s/N,mu_c/N];

Nc_l = length(Nl_ind);
Nr_l = length(Nr_ind);

%% Loading fitted parameters for the OU process
sgn_fit_fun = 's3';

sgn_tau1 = '_0d5';

sgn_taumax = '_2d5';

sgn_L = '';
%sgn_L = '_n'; 


%pa_fit_1 = load(['..\','fit_result_','ks_sim_Delta_',setting_sgn,'_autocov_OU_','a3','_L','.mat']);
%pa_fit_1 = load(['..\','fit_result_','ks_sim_Delta_',setting_sgn,'_autocov_OU_',sgn_fit_fun,'_L','_taumax',sgn_taumax,'.mat']);
pa_fit_1 = load([fileloc,'fit_result_','ks_sim_Delta_',setting_sgn_pa,'_autocov_OU_',sgn_fit_fun,'_L',sgn_L,'_tau1',sgn_tau1,'_taumax',sgn_taumax,'.mat']);
p_a = pa_fit_1.p_a;


switch sgn_fit_fun
    case 'ga2_s3'
        gamma1 = p_a(1);
        gamma2 = p_a(2);
        omega = p_a(3);
        Sig = [p_a(4),p_a(5);p_a(5),p_a(6)];
        SS = Sig*Sig';
        a1 = SS(1,1);
        a2 = SS(1,2);
        a3 = SS(2,2);
        %p_a_1 = [gamma1,gamma2,omega,a1,a2,a3];
    case 's3'
        gamma1 = p_a(1);
        gamma2 = p_a(1);
        omega = p_a(2);
        Sig = [p_a(3),p_a(4);p_a(4),p_a(5)];
        SS = Sig*Sig';
        a1 = SS(1,1);
        a2 = SS(1,2);
        a3 = SS(2,2);
    case 's1'
        gamma1 = p_a(1);
        gamma2 = p_a(1);
        omega = p_a(2);
        Sig = [p_a(3),0;0,p_a(3)];
        SS = Sig*Sig';
        a1 = SS(1,1);
        a2 = SS(1,2);
        a3 = SS(2,2);
    case 'a3'
        gamma1 = p_a(1);
        gamma2 = p_a(1);
        omega = p_a(2);
        a1 = p_a(3);
        a2 = p_a(4);
        a3 = p_a(5);
end

p_a_1 = [gamma1,gamma2,omega,a1,a2,a3];

%p_a(1) = 0.8;
    

%% compute the corresponding auto-correlation, and integrate it

flg_sig = 's3';
f_fun = @(p,tt,flg_var) f_fun_OU_ac_1_s3(p,tt,flg_var);

tt = 0:0.01:100;

f_fun_a = f_fun(p_a,tt,'a');


a = trapz(tt,f_fun_a,1);

disp(trapz(tt,f_fun_a,1) - f_fun_OU_ac_1_s3_int(p_a))
%%
figure
plot(tt,f_fun_a(:,1))


%% Setting up for simulation of the approximate SDE

% time-step size
dt = 0.01;

% time
t0 = 0;
tmax = 1e4;

% forming a vector for the time-steps
tt = t0:dt:tmax;
%tl = length(tt);
tl = floor((tmax-t0+0.5*dt)/dt)+1;

% transient
flg_transient = true;
t1 = 1e3;


%%% dividing into intervals of time to simulate within each interval

% no. of time-steps within each interval
Delta_tt = 1e4;

% total no. of intervals
tl_Delta = floor((tl-1+0.5)/Delta_tt);


% setting up storage for simulation of the SDE
%zt = zeros(tl,1);
%xt = zeros(tl,Nc_l);


zt_traj = zeros(tl,1);
yt_traj = zeros(tl,2);

% ind_xt_a = [1:5,10,50,100,112:116];
% ind_xt_a_l = length(ind_xt_a);

ind_xt_a = Nl_ind(1:300:end);
ind_xt_a_l = length(ind_xt_a);

xt_a_traj = zeros(tl,ind_xt_a_l);

xt_a_Delta_mean = zeros(tl_Delta,Nc_l);
xt_a_2_Delta_mean = zeros(tl_Delta,Nc_l);

% xt_a_1_Delta_mean = zeros(tl_Delta,Nc_l-1);
% xt_a_1_2_Delta_mean = zeros(tl_Delta,Nc_l-1);


%% Simulate the SDE in intervals of time-steps

% setting rng
seed = 0;
rng(seed);

% initial conditions
x0 = rand(1,Nc_l)*2*pi-pi;
%y0 = [mu_s/N,mu_c/N];
y0 = [0,0];

% transient
if flg_transient == true
    [~,xt1,yt1] = euma_ks_sde_OU_1([0,t1],param_1,data_5,p_a_1,x0,y0);
    x_t1 = xt1(end,:);
    y_t1 = yt1(end,:);
else
    x_t1 = x0;
    y_t1 = y0;
end

x_ti = adjust_angles(x_t1-angle(mean(exp(1i*x_t1))));
y_ti = y_t1;

% simulate the approximate SDE with Euler-Maruyama algorithm
%[~,~,yt,zt] = euma_ks_sde_OU_1(tt,param_1,data_5,p_a,x_ti,y_ti);

tic
for ti_Delta = 1:tl_Delta
    
    tt_ti_Delta_ind = (ti_Delta-1)*Delta_tt+1:ti_Delta*Delta_tt+1;
    %tt_ti_Delta = tt(tt_ti_Delta_ind);
    tt_ti_Delta = ((ti_Delta-1)*Delta_tt:ti_Delta*Delta_tt)*dt;
    
    [~,xt,yt,zt] = euma_ks_sde_OU_1(tt_ti_Delta,param_1,data_5,p_a_1,x_ti,y_ti);
    
    zt_traj(tt_ti_Delta_ind,:) = zt;
    yt_traj(tt_ti_Delta_ind,:) = yt;
    
    reipsic_t = mean(exp(1i*xt),2);
    rc_t = abs(reipsic_t);
    psic_t = angle(reipsic_t);
    
    %psic_1_t = angle(mean(exp(1i*xt(:,1:115)),2));
    %psic_1_t = psic_t;
    
    %ang_adj_diff = (xt(1,:)-psic_1_t(1)) - adjust_angles(xt(1,:)-psic_1_t(1));
    
    %xt_a_Delta_mean(ti_Delta,:) = mean(adj_ang_a(xt(1:end-1,:) - psic_1_t(1:end-1)) - ang_adj_diff,1);
    %xt_a_2_Delta_mean(ti_Delta,:) = mean((adj_ang_a(xt(1:end-1,:) - psic_1_t(1:end-1)) - ang_adj_diff).^2,1);
    
    xt_a_Delta_mean(ti_Delta,:) = mean(adjust_angles(xt(1:end-1,:)-psic_t(1:end-1)),1);
    xt_a_2_Delta_mean(ti_Delta,:) = mean(adjust_angles(xt(1:end-1,:)-psic_t(1:end-1)).^2,1);
    
%     xt_a_1_Delta_mean(ti_Delta,:) = mean(adjust_angles(xt(1:end-1,1:115) - psic_1_t(1:end-1)),1);
%     xt_a_1_2_Delta_mean(ti_Delta,:) = mean(adjust_angles(xt(1:end-1,1:115) - psic_1_t(1:end-1)).^2,1);
    
    xt_a_traj(tt_ti_Delta_ind,:) = adjust_angles(xt(:,ind_xt_a) - psic_t);
    
    x_ti = xt(end,:);
    y_ti = yt(end,:);
end
t = toc

%% calculating some statistics
xt_a_Delta_mu = mean(xt_a_Delta_mean,1);
xt_a_Delta_var = 1/(tl-2)*(tl-1)*(mean(xt_a_2_Delta_mean,1) - mean(xt_a_Delta_mean,1).^2);

xt_a_1_Delta_mu = mean(xt_a_1_Delta_mean,1);
xt_a_1_Delta_var = 1/(tl-2)*(tl-1)*(mean(xt_a_1_2_Delta_mean,1) - mean(xt_a_1_Delta_mean,1).^2);

%% plotting some results
figure
plot(tt_ti_Delta,xt)
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)

%%
% save(['sde_approx_1_a_1_','ks_sim_Delta_',setting_sgn,'_autocov_OU','_a3','_L','_tau1',sgn_tau1,'_taumax',sgn_taumax,'.mat'],...
%     'p_a','t0','dt','tmax','flg_transient','t1','Delta_tt','zt_traj','yt_traj')
save(['sde_approx_1_a_1_','ks_sim_Delta_',setting_sgn,'_autocov_OU_','pa_',setting_sgn_pa,'_',sgn_fit_fun,'_L',sgn_L,'_tau1',sgn_tau1,'_taumax',sgn_taumax,'_t1_',int2str(t1),'_tmax_',int2str(tmax),'.mat'],...
    'p_a','p_a_1','t0','dt','tmax','tl','flg_transient','t1','Delta_tt','tl_Delta','zt_traj','yt_traj',...
    'xt_a_Delta_mean','xt_a_2_Delta_mean','xt_a_Delta_mu','xt_a_Delta_var','ind_xt_a','xt_a_traj') 

%%

% tt = t0:dt:tmax;
% 
% %
% figure
% plot(tt,abs(zt_traj))
% 
% 
% figure
% histogram(abs(zt_traj),'normalization','pdf','numbins',50)

