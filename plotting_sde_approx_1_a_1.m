%% Plotting and comparing results from the approxmiate SDE



%% clear variables
clearvars

%% location of setting, data
fileloc = '.\';
%fileloc = 'D:\MATLAB\scraps\ks_3\';

%% parameter setting
setting_param_sgn = '1_a';

setting_N_sgn = '_4N';

setting_sgn = [setting_param_sgn,setting_N_sgn];

setting_sgn_pa = setting_sgn;
%setting_sgn_pa = '1_a_32N';


sgn_fit_fun = 's3';

sgn_tau1 = '_0d5';

sgn_taumax = '_2d5';

sgn_L = '';
%sgn_L = '_n';

t1 = 1e3;

tmax = 1e4;

%% load parameters
param_1 = load([fileloc,'params_ks_sim_',setting_sgn,'.mat']);

N = param_1.N;
K = param_1.K;

l = param_1.l;
%% load data

data_5 = load([fileloc,'data_ks_sim_',setting_sgn,'_a','.mat']);

Nc_l = length(data_5.Nl_ind);
Nr_l = length(data_5.Nr_ind);

mu_s = data_5.mu_S_t;
mu_c = data_5.mu_C_t;

mu_x = [mu_s/N,mu_c/N];

%% load values of fitted parameters
pa_fit_1 = load([fileloc,'fit_result_','ks_sim_Delta_',setting_sgn,'_autocov_OU_',sgn_fit_fun,'_L',sgn_L,'_tau1',sgn_tau1,'_taumax',sgn_taumax,'.mat']);
p_a = pa_fit_1.p_a;

%% load results from simulation of the SDE
%r1_sdeq  q = load(['sde_approx_1_a_1_','ks_sim_Delta_',setting_sgn,'_autocov_OU','_a3','_L','_taumax',sgn_taumax,'.mat']);
r1_sde_OU = load(['sde_approx_1_a_1_','ks_sim_Delta_',setting_sgn,'_autocov_OU_','pa_',setting_sgn_pa,'_',sgn_fit_fun,'_L',sgn_L,'_tau1',sgn_tau1,'_taumax',sgn_taumax,'_t1_',int2str(t1),'_tmax_',int2str(tmax),'.mat']);

%r1_sde_BM = load('sde_approx_BM_1_ks_sim_Delta_1_a_4N_t1_1000_tmax_1000.mat');
%r1_sde_BM = load('sde_approx_BM_1_a_ks_sim_Delta_1_a_4N_t1_1000_tmax_1000.mat');
%r1_sde_BM = load('sde_approx_BM_1_a_int_ks_sim_Delta_1_a_4N_t1_1000_tmax_10000.mat');

r1_sde_BM = load('sde_approx_BM_with_correction_1_a_int_ks_sim_Delta_1_a_4N_t1_1000_tmax_10000.mat');

r1 = r1_sde_OU;
%r1 = r1_sde_BM;

%% load results from direct simulation of the full KS system
r1_sim = load([fileloc,'ks_sim_Delta_',setting_sgn,'_1','.mat']);
%r1_sim = load([fileloc,'ks_sim_Delta_',setting_sgn,'_10','.mat']);

%%
%r1_sim_10 = load([fileloc,'ks_sim_',setting_sgn,'_10','.mat']);

%% load statistics from direct simulation results
data_10 = load([fileloc,'data_ks_sim_Delta_',setting_sgn,'_15','.mat']);

%% calculate the mean-field variables from SDE result with OU

% result from SDE
yt_sde_OU = r1_sde_OU.yt_traj;
zt_sde_OU = r1_sde_OU.zt_traj;

% x_t = (\xi_t, \eta_t)
x_t_sde_OU = mu_x + 1/sqrt(N)*yt_sde_OU;

% mean-field variables for the synchronized cluster, from SDE
rc_t_sde_OU = abs(zt_sde_OU);
psic_t_sde_OU = angle(zt_sde_OU);

%% calculate the mean-field variables from SDE result with BM

% result from SDE
yt_sde_BM = r1_sde_BM.yt_traj;
zt_sde_BM = r1_sde_BM.zt_traj;

% x_t = (\xi_t, \eta_t)
x_t_sde_BM = mu_x + 1/sqrt(N)*yt_sde_BM;

% mean-field variables for the synchronized cluster, from SDE
rc_t_sde_BM = abs(zt_sde_BM);
psic_t_sde_BM = angle(zt_sde_BM);

%% Plotting some results from SDE & compare with result from full simulation
%% histogram of r_c for SDE with OU
figure
histogram(r1_sim.rc_t_traj,'normalization','pdf','numbins',50)
hold on
histogram(rc_t_sde_OU,'normalization','pdf','numbins',50)
set(gca,'fontsize',15)
xlabel('$r_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(r_c)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'interpreter','latex','fontsize',15,'location','nw')

%disp('mean of rc')
fprintf('mean of rc: KS: %.4g,   SDE: %.4g\n',mean(r1_sim.rc_t_traj),mean(rc_t_sde_OU))
fprintf('Var of rc: KS: %.4g,   SDE: %.4g\n',var(r1_sim.rc_t_traj),var(rc_t_sde_OU))
fprintf('skewness of rc: KS: %.4e,   SDE: %.4e\n',skewness(r1_sim.rc_t_traj),skewness(rc_t_sde_OU))
fprintf('kurtosis of rc: KS: %.4e,   SDE: %.4e\n',kurtosis(r1_sim.rc_t_traj),kurtosis(rc_t_sde_OU))



figure
hold on
histogram_line(r1_sim.rc_t_traj,[],'-',2.4)
histogram_line(rc_t_sde_OU,[],'-.',2.4)
hold off
ylim([0,8000])
set(gca,'fontsize',15)
xlabel('$r_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(r_c)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','nw')


%% histogram of r_c for SDE with BM
figure
histogram(r1_sim.rc_t_traj,'normalization','pdf','numbins',50)
hold on
histogram(rc_t_sde_BM,'normalization','pdf','numbins',50)
set(gca,'fontsize',15)
xlabel('$r_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(r_c)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'interpreter','latex','fontsize',15,'location','nw')

%disp('mean of rc')
fprintf('mean of rc: KS: %.4g,   SDE: %.4g\n',mean(r1_sim.rc_t_traj),mean(rc_t_sde_BM))
fprintf('Var of rc: KS: %.4g,   SDE: %.4g\n',var(r1_sim.rc_t_traj),var(rc_t_sde_BM))
fprintf('skewness of rc: KS: %.4e,   SDE: %.4e\n',skewness(r1_sim.rc_t_traj),skewness(rc_t_sde_BM))


figure
hold on
histogram_line(r1_sim.rc_t_traj,50,'-',2.4)
histogram_line(rc_t_sde_BM,50,'-.',2.4)
hold off
ylim([0,170])
set(gca,'fontsize',15)
xlabel('$r_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(r_c)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (BM)'},'interpreter','latex','fontsize',15,'location','nw')


%% for SDE with OU, calculate r by combining SDE solution giving rc, and S,C

% mean-field variable for the rogue oscillators, from OU process
reipsir_t_sde_OU = 1/Nr_l*N*(x_t_sde_OU(:,2)+1i*x_t_sde_OU(:,1))*exp(-1i*l);

% estimate mean-field variables from SDE by combining cluster and rogue
reipsit_sde_OU = Nc_l/N*rc_t_sde_OU + Nr_l/N*reipsir_t_sde_OU;
rt_sde_OU = abs(reipsit_sde_OU);
psit_sde_OU = angle(reipsit_sde_OU);



% estimate mean-field variables where synchronized cluster appear only as
% mean
rc_t_sde_OU_mean = mean(rc_t_sde_OU);
reipsit_sde_OU_rc_mean = Nc_l/N*rc_t_sde_OU_mean + Nr_l/N*reipsir_t_sde_OU;
rt_sde_OU_rc_mean = abs(reipsit_sde_OU_rc_mean);
psit_sde_OU_rc_mean = angle(reipsit_sde_OU_rc_mean);


% estimate mean-field variables
a = Nc_l/N*rc_t_sde_OU_mean + mean(Nr_l/N*reipsir_t_sde_OU,1);
b = Nr_l/N*reipsir_t_sde_OU - mean(Nr_l/N*reipsir_t_sde_OU,1);
rt_sde_OU_rc_mean_approx = abs(a)+ real(b*exp(-1i*angle(a)));

a_1 = Nc_l/N*rc_t_sde_OU_mean + (mu_x(2)+1i*mu_x(1))*exp(-1i*l);
b_1 = 1/sqrt(N)*(yt_sde_OU(:,2)+1i*yt_sde_OU(:,1))*exp(-1i*l);
z_1 = abs(a_1)+real(b_1*exp(-1i*angle(a_1)));

%
dz = 1/sqrt(N)*real(((diff(yt_sde_OU(:,2))+1i*diff(yt_sde_OU(:,1)))*exp(-1i*l)*exp(-1i*angle(a))));
dz1 = diff(real(b*exp(-1i*angle(a))));

dz_1 = 1/sqrt(N)*real(((diff(yt_sde_OU(:,2))+1i*diff(yt_sde_OU(:,1)))*exp(-1i*l)*exp(-1i*angle(a_1))));

z_t = cumsum([abs(a);dz1]);

zt_1 = cumsum([abs(a_1);dz_1]);

phi = angle(a_1);

ga = p_a(1);
om = p_a(2);
Sig = [p_a(3),p_a(4);p_a(4),p_a(5)];
SS = Sig*Sig';
a1 = SS(1,1);
a2 = SS(1,2);
a3 = SS(2,2);


Gamma = [-ga,0;0,-ga];
Omega = [0,om;-om,0];

t0 = r1_sde_OU.t0;
dt = r1_sde_OU.dt;
tmax = r1_sde_OU.tmax;

tt = t0:dt:tmax;
tl = length(tt);

dr = zeros(tl,1);

yt = zeros(tl,2);
y0 = [0;0];

yi = y0;

for ni = 1:tl
    yt(ni,:) = yi;
    
    dyi = (Gamma+Omega)*yi*dt + mvnrnd([0,0],SS,1)'*sqrt(dt);
    
    yi = yi + dyi;
    
    dr(ni,:) = 1/sqrt(N)*(dyi(2)*cos(phi+l)+dyi(1)*sin(phi+l));
    
end

rt = abs(a_1) + cumsum([0;dr(1:end-1,:)]);

figure
histogram_line(rt_sde_OU_rc_mean_approx,[],'-.',2.4)
hold on
histogram_line(rt,[],'-.',2.4)


%% comparing histogram of r for whole system from SDE and from direct sim
figure
histogram(r1_sim.rt_traj,'normalization','pdf','numbins',50)
hold on
histogram(rt_sde_OU,'normalization','pdf','numbins',50)
set(gca,'fontsize',15)
xlabel('$r$','interpreter','latex','fontsize',20)
ylabel('$\rho(r)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'fontsize',15,'location','nw')

fprintf('mean of r: KS: %.4g,   SDE: %.4g\n',mean(r1_sim.rt_traj),mean(rt_sde_OU))
fprintf('Var of r: KS: %.4e,   SDE: %.4e\n',var(r1_sim.rt_traj),var(rt_sde_OU))
fprintf('skewness of r: KS: %.4e,   SDE: %.4e\n',skewness(r1_sim.rt_traj),skewness(rt_sde_OU))
fprintf('kurtosis of r: KS: %.4e,   SDE: %.4e\n',kurtosis(r1_sim.rt_traj),kurtosis(rt_sde_OU))

figure
hold on
histogram_line(r1_sim.rt_traj,[],'-',2.4)
histogram_line(rt_sde_OU,[],'-.',2.4)
histogram_line(rt_sde_OU_rc_mean,[],'-.',2.4)
histogram_line(rt_sde_OU_rc_mean_approx,[],'-.',2.4)
hold off
%xlim([0.6,0.9])
ylim([0,400])
set(gca,'fontsize',15)
xlabel('$r$','interpreter','latex','fontsize',20)
ylabel('$\rho(r)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','nw')


%% for SDE with BM, calculate r by combining SDE solution giving rc, and S,C


% mean-field variable for the rogue oscillators, from BM process
reipsir_t_sde_BM = 1/Nr_l*N*(x_t_sde_BM(:,2)+1i*x_t_sde_BM(:,1))*exp(-1i*l);

% estimate mean-field variables from SDE by combining cluster and rogue
reipsit_sde_BM = Nc_l/N*rc_t_sde_BM + Nr_l/N*reipsir_t_sde_BM;
rt_sde_BM = abs(reipsit_sde_BM);
psit_sde_BM = angle(reipsit_sde_BM);

% estimate mean-field variables where synchronized cluster appear only as
% mean
rc_t_sde_BM_mean = mean(rc_t_sde_BM);
reipsit_sde_BM_rc_mean = Nc_l/N*rc_t_sde_BM_mean + Nr_l/N*reipsir_t_sde_BM;
rt_sde_BM_rc_mean = abs(reipsit_sde_BM_rc_mean);
psit_sde_BM_rc_mean = angle(reipsit_sde_BM_rc_mean);


%% comparing histogram of r for whole system from SDE and from direct sim
figure
histogram(r1_sim.rt_traj,'normalization','pdf','numbins',50)
hold on
histogram(rt_sde_BM,'normalization','pdf','numbins',50)
set(gca,'fontsize',15)
xlabel('$r$','interpreter','latex','fontsize',20)
ylabel('$\rho(r)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'fontsize',15,'location','nw')

fprintf('mean of r: KS: %.4g,   SDE: %.4g\n',mean(r1_sim.rt_traj),mean(rt_sde_BM))
fprintf('Var of r: KS: %.4e,   SDE: %.4e\n',var(r1_sim.rt_traj),var(rt_sde_BM))
fprintf('skewness of r: KS: %.4e,   SDE: %.4e\n',skewness(r1_sim.rt_traj),skewness(rt_sde_BM))

figure
hold on
histogram_line(r1_sim.rt_traj,[],'-',2.4)
histogram_line(rt_sde_BM,[],'-.',2.4)
%histogram_line(rt_sde_BM_rc_mean,[],'-.',2.4)
hold off
%xlim([0.6,0.9])
ylim([0,25])
set(gca,'fontsize',15)
xlabel('$r$','interpreter','latex','fontsize',20)
ylabel('$\rho(r)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (BM)'},'interpreter','latex','fontsize',15,'location','nw')



%% comparing histogram of S(t), C(t) from SDE and from direct simulation

%% S
% histogram of S
figure
histogram(r1_sim.S_t_traj,'normalization','pdf','numbins',50)
hold on
histogram(x_t_sde_OU(:,1),'normalization','pdf','numbins',50)
set(gca,'fontsize',15)
xlabel('$S$','interpreter','latex','fontsize',20)
ylabel('$\rho(S)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'fontsize',15,'location','nw')

disp('mean of S')
disp([mean(r1_sim.S_t_traj),mean(x_t_sde_OU(:,1))])
disp('variance of S')
disp([var(r1_sim.S_t_traj),var(x_t_sde_OU(:,1))])

%% C
% histogram of C
figure
histogram(r1_sim.C_t_traj,'normalization','pdf','numbins',50)
hold on
histogram(x_t_sde_OU(:,2),'normalization','pdf','numbins',50)
xlabel('$C$','interpreter','latex','fontsize',20)
ylabel('$\rho(C)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'fontsize',15)

disp('mean of C')
disp([mean(r1_sim.C_t_traj),mean(x_t_sde_OU(:,2))])

disp('variance of C')
disp([var(r1_sim.C_t_traj),var(x_t_sde_OU(:,2))])


%% Compare xi, zeta from SDE, analy and from KS

%% some data and analytical result and full KS
%  plotting labels
plot_label_range_1 = {'\xi','\zeta'};

yt_cov_analy = f_fun_OU_ac_1_s3(p_a,0,'a');
yt_var_analy = yt_cov_analy([1,4]);

yt_traj_sim = sqrt(r1_sim.N)*([r1_sim.S_t_traj,r1_sim.C_t_traj]);
yt_traj_sim = yt_traj_sim - mean(yt_traj_sim,1);

%% check convergence
ind_plot = 1;
plot_label = plot_label_range_1{ind_plot};

mu_x = 0; var_x = yt_var_analy(ind_plot);
zz = -5:0.01:5;
xx = mu_x + sqrt(var_x)*zz;
yy = 1/sqrt(2*pi*var_x)*exp(-1/2*zz.^2);


figure
hold on
histogram_line(r1.yt_traj(1:floor(r1.tl/2),ind_plot))
histogram_line(r1.yt_traj(ceil(r1.tl/2):r1.tl,ind_plot))
plot(xx,yy,'--','color',[0.5,0.5,0.5])
hold off
set(gca,'fontsize',15)
xlabel(sprintf('$%s$',plot_label),'interpreter','latex','fontsize',20)
ylabel(sprintf('$\\rho(%s)$',plot_label),'interpreter','latex','fontsize',20)
legend({'$t=0:5e4$','$t=5e4:1e5$'},'interpreter','latex','fontsize',15)

% legend({sprintf('$t = %d \\textrm{ to } %d$',0,tmax/2),...
%     sprintf('$t = %d \\textrm{ to } %d$',tmax/2,tmax)},'interpreter','latex','fontsize',10)

figure
hold on
histogram(r1.yt_traj(1:floor(r1.tl/2),ind_plot),'normalization','pdf','numbins',50)
histogram(r1.yt_traj(ceil(r1.tl/2):r1.tl,ind_plot),'normalization','pdf','numbins',50)
plot(xx,yy,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold off
set(gca,'fontsize',15)
xlabel(sprintf('$%s$',plot_label),'interpreter','latex','fontsize',20)
ylabel(sprintf('$\\rho(%s)$',plot_label),'interpreter','latex','fontsize',20)
legend({'$t=0:5e4$','$t=5e4:1e5$'},'interpreter','latex','fontsize',15)


%% histogram of xi, zeta, compare with analy and KS

ind_plot = 2;
plot_label = plot_label_range_1{ind_plot};

mu_x = 0; var_x = yt_var_analy(ind_plot);
zz = -5:0.01:5;
xx = mu_x + sqrt(var_x)*zz;
yy = 1/sqrt(2*pi*var_x)*exp(-1/2*zz.^2);


figure
hold on
histogram_line(yt_traj_sim(:,ind_plot),[],'-',2)
histogram_line(r1.yt_traj(:,ind_plot),[],'-.',2)
plot(xx,yy,':','linewidth',2)
hold off
set(gca,'fontsize',15)
xlabel(sprintf('$%s$',plot_label),'interpreter','latex','fontsize',20)
ylabel(sprintf('$\\rho(%s)$',plot_label),'interpreter','latex','fontsize',20)
%legend({'reduced SDE','analy'},'interpreter','latex','fontsize',15)
legend({'KS model','reduced SDE','analy'},'interpreter','latex','fontsize',15)


% figure
% hold on
% histogram(yt_traj(:,ind_plot),'normalization','pdf','numbins',50)
% plot(xx,yy,'r-')



%%% calculate the auto-correlation function for solution from SDE

%% selecting data

%% time-steps for calculating the autocovariance function
tauN0 = 0; dtauN = 10; tauNmax = 3000;
tauNt = tauN0:dtauN:tauNmax;
tauNl = length(tauNt);

dt = r1.dt;

tt_tauNt = tauNt*dt;

%tauN1 = 32e4;
tauN1 = 5e6;
t_tauN1 = tauN1*dt;

%%
tauNt_1 = tauNt/5;
tauN1_1 = tauN1/5;

%tt_tauNt_1 = tt_tauNt*r1_sim.dt;
tt_tauNt_1 = tauNt_1*r1_sim.dt;

%% calculating the auto-covariance

%%
% auto-covariance of OU component of the SDE
yt_sde_OU_covar = x_vec_covar(tauNt,yt_sde_OU,tauN1,'v');

% auto-covaraince from analytical result of an actual 2d OU process
yt_autocov_analy = f_fun_OU_ac_1_s3(p_a,tt_tauNt,'a');

% auto-covariance of S, C from direct simulation trajectory
%yt_autocov_sim = x_vec_covar(tauNt,sqrt(r1_sim.N)*[r1_sim.S_t_traj,r1_sim.C_t_traj],tauN1,'v');
yt_autocov_sim = x_vec_covar(tauNt_1,yt_traj_sim,tauN1_1,'v');

%x_t_ve
    
%x_t_


%% auto-covariance of rc_t from SDE and from KS
rc_t_sim_vec_covar = x_vec_covar(tauNt_1,r1_sim.rc_t_traj,tauN1_1);
rc_t_sde_OU_vec_covar = x_vec_covar(tauNt,rc_t_sde_OU,tauN1);
rc_t_sde_BM_vec_covar = x_vec_covar(tauNt,rc_t_sde_BM,tauN1);


%% auto-covariance of rt from SDE and from KS
rt_sim_vec_covar = x_vec_covar(tauNt_1,r1_sim.rt_traj,tauN1_1);
rt_sde_OU_vec_covar = x_vec_covar(tauNt,rt_sde_OU,tauN1);
rt_sde_BM_vec_covar = x_vec_covar(tauNt,rt_sde_BM,tauN1);


%% plotting the auto-covariance

%% auto-covariance of OU component of SDE

plot_label_range = {'\xi\xi','\zeta\xi','\xi\zeta','\zeta\zeta'};

ind_plot = 4;
plot_label = plot_label_range{ind_plot};
figure
hold on
plot(tt_tauNt_1,yt_autocov_sim(:,ind_plot),'linewidth',1)
plot(tt_tauNt,yt_sde_OU_covar(:,ind_plot),'-.','linewidth',1)
plot(tt_tauNt,yt_autocov_analy(:,ind_plot),'--','color',[0.5,0.5,0.5],'linewidth',1)
xlim([0,10])
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel(sprintf('$R_{%s}(\\tau)$',plot_label),'interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE','analy'},'interpreter','latex','fontsize',15)
%legend({'reduced SDE','analy'},'interpreter','latex','fontsize',15)

%% autocovariance of rc_t
figure
plot(tt_tauNt_1,rc_t_sim_vec_covar,'linewidth',1)
hold on
plot(tt_tauNt,rc_t_sde_OU_vec_covar,'-.','linewidth',1)
plot(tt_tauNt,rc_t_sde_BM_vec_covar,'-.','linewidth',1)
xlim([0,20])
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel('$R_{r_c}(\tau)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)','reduced SDE (BM)'},'fontsize',15)


%% autocovariance of rt
figure
plot(tt_tauNt_1,rt_sim_vec_covar,'linewidth',1)
hold on
plot(tt_tauNt,rt_sde_OU_vec_covar,'-.','linewidth',1)
plot(tt_tauNt,rt_sde_BM_vec_covar,'-.','linewidth',1)
xlim([0,10])
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel('$R_r(\tau)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'fontsize',15)

%%
%figure
%plot()



%% plotting the time series of some variables

tt_sde_OU = r1_sde_OU.t0:r1_sde_OU.dt:r1_sde_OU.tmax;
tt_sde_BM = r1_sde_BM.t0:r1_sde_BM.dt:r1_sde_BM.tmax;

%% trajectory of rc_t
figure
plot(r1_sim.tt,r1_sim.rc_t_traj,'-','linewidth',2.4)
hold on
plot(tt_sde_OU,rc_t_sde_OU,'-','linewidth',2.4)
%plot(tt_sde_BM,rc_t_sde_BM,'-','linewidth',2.4)
xlim([0,20])
ylim([0.91,0.95])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$r_c(t)$','interpreter','latex','fontsize',20)
%legend({'KS model','reduced SDE (OU)','reduced SDE (BM)'},'interpreter','latex','fontsize',15,'location','nw')
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','nw')


%% trajectory of rt
figure
plot(r1_sim.tt,r1_sim.rt_traj,'b','linewidth',0.5)
hold on
plot(tt_sde_OU,rt_sde_OU,'r','linewidth',0.5)
plot(tt_sde_BM,rt_sde_BM,'m','linewidth',0.5)
%plot(tt,rt_sde_rc_mean,'.-','linewidth',1)
xlim([0,10])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$r(t)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'fontsize',15)

%%
figure
i_ind = 10;
phit_i_a = adj_ang_a(r1_sim_10.phit(:,i_ind)-r1_sim_10.omega_*r1_sim_10.tt');
plot(r1_sim_10.tt,phit_t_a-phit_i_a(1))



%%
figure
plot(data_10)

%% calculate statistics of individual synchronized oscillators from SDE

phit_1 = xt_sde_OU(1:end-1,1:115);
%phit_1 = xt(:,1:end);

psic_1_t_sde = angle(mean(exp(1i*phit_1),2));

% ang_adj_diff = (phit_1(1,:)-psic_1_t_sde(1)) - adjust_angles(phit_1(1,:)-psic_1_t_sde(1));
% phit_c_a_sde = adj_ang_a(phit_1 - ang_adj_diff-psic_1_t_sde);

phit_c_a_1_sde = adjust_angles(phit_1 - psic_1_t_sde);

phit_c_a_1_sde_mu = mean(phit_c_a_1_sde);
phit_c_a_1_sde_var = var(phit_c_a_1_sde);

%%
phit_c_a_sde = adjust_angles(xt(1:end-1,:) - psic_t_sde(1:end-1));
phit_c_a_sde_mu = mean(phit_c_a_sde);
phit_c_a_sde_var = var(phit_c_a_sde);

%%


%%
figure
plot(phit_c_a_sde_mu,'.')
hold on
plot(data_10.phit_c_a_Delta_mu,'.')

xlim([1,115])

figure
plot(phit_c_a_sde_mu - data_10.phit_c_a_Delta_mu(1:116),'.')
xlim([1,115])

figure
plot(tt,adj_ang_a(angle(zt)))

figure
plot(phit_c_a_sde_var,'.')
hold on
plot(data_10.phit_c_a_Delta_var,'.')
xlim([1,115])

figure
plot(phit_c_a_sde_var - data_10.phit_c_a_Delta_var,'.')


%% plotting mean and variance for phit_c_a from SDE and from full KS

% mean
figure
hold on
plot(data_10.phit_c_a_Delta_mu,'o','markersize',10,'linewidth',1)
plot(r1_sde_OU.xt_a_Delta_mu,'x','markersize',10,'linewidth',1)
hold off
xlim([1,116])
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
ylabel('$\overline{\theta_i - \psi_c}$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','nw')

% difference in the mean betweem SDE and full KS
figure
plot(r1_sde_OU.xt_a_Delta_mu - data_10.phit_c_a_Delta_mu,'o')
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)

% variance
figure
hold on
plot(data_10.phit_c_a_Delta_var,'o','markersize',10,'linewidth',1)
plot(r1_sde_OU.xt_a_Delta_var,'x','markersize',10,'linewidth',1)
hold off
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
ylabel('$\textrm{Var}(\theta_i-\psi_c)$','interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','nw')
xlim([1,115])

% difference in the variance between sde and full KS
figure
hold on
plot(r1_sde_OU.xt_a_Delta_var - data_10.phit_c_a_Delta_var,'o','markersize',7,'linewidth',1)
plot(data_5.Nl_ind,zeros(1,length(data_5.Nl_ind)),'--','color',[0.5,0.5,0.5],'linewidth',2)
hold off
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
ylabel('Var(SDE) - Var(full KS)','interpreter','latex','fontsize',20)
xlim([1,115])



%%  plotting for phit_c_a_1
figure
hold on
plot(data_10.phit_c_a_1_Delta_mu,'o')
plot(xt_a_1_Delta_mu,'x')
hold off

figure
plot(xt_a_1_Delta_mu - data_10.phit_c_a_1_Delta_mu,'.')


figure
hold on
plot(data_10.phit_c_a_1_Delta_var,'o')
plot(xt_a_1_Delta_var,'x')
hold off

figure
plot(xt_a_1_Delta_var - data_10.phit_c_a_1_Delta_var,'.')

%% histogram of phases of some synchronized oscillators

% figure
% histogram(phit_c_a_sde(:,50),'normalization','pdf','numbins',50)

ind_xt_a = r1_sde_OU.ind_xt_a;

%ind_i_plot_range = 1:length(ind_xt_a);
ind_i_plot_range = [1:3,6:8,10:12];

for ind_i = ind_i_plot_range

    ind_phi_i = ind_xt_a(ind_i);

    figure
    hold on
    histogram_line(r1_sim.phit_c_a_traj(:,ind_i),50,'-',2.4)
    histogram_line(r1_sde_OU.xt_a_traj(:,ind_i),50,'-.',2.4)
    hold off
    set(gca,'fontsize',15)
    xlabel('$\theta_i-\psi_c$','interpreter','latex','fontsize',20)
    ylabel('$\rho(\theta_i-\psi_c)$','interpreter','latex','fontsize',20)
    title(sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
    legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','ne')

end

%%

ind_i = 1;
ind_phi_i = ind_xt_a(ind_i);

figure
hold on
histogram_line(r1_sim.phit_c_a_traj(:,ind_i),50,'-',2.4)
histogram_line(r1_sde_OU.xt_a_traj(:,ind_i),50,'-.',2.4)
ylim([0,25])
hold off
set(gca,'fontsize',15)
xlabel('$\theta_i-\psi_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(\theta_i-\psi_c)$','interpreter','latex','fontsize',20)
%title(sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
text(-1.085,22.5,sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','ne')

%%
ind_i = 7;
ind_phi_i = ind_xt_a(ind_i);

figure
hold on
histogram_line(r1_sim.phit_c_a_traj(:,ind_i),50,'-',2.4)
histogram_line(r1_sde_OU.xt_a_traj(:,ind_i),50,'-.',2.4)
ylim([0,120])
hold off
set(gca,'fontsize',15)
xlabel('$\theta_i-\psi_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(\theta_i-\psi_c)$','interpreter','latex','fontsize',20)
%title(sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
text(-0.097,109,sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','ne')

%%
ind_i = 12;
ind_phi_i = ind_xt_a(ind_i);

figure
hold on
histogram_line(r1_sim.phit_c_a_traj(:,ind_i),50,'-',2.4)
histogram_line(r1_sde_OU.xt_a_traj(:,ind_i),50,'-.',2.4)
ylim([0,15])
hold off
set(gca,'fontsize',15)
xlabel('$\theta_i-\psi_c$','interpreter','latex','fontsize',20)
ylabel('$\rho(\theta_i-\psi_c)$','interpreter','latex','fontsize',20)
%title(sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
text(0.625,13.5,sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE (OU)'},'interpreter','latex','fontsize',15,'location','ne')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ind_i = 8;
ind_phi_i = ind_xt_a(ind_i);
figure
hold on
histogram_line(r1_sim.phit_c_a_traj(:,ind_i),50,'-',2)
histogram_line(r1_sde_OU.xt_a_traj(:,ind_i),50,'-.',2)
hold off
set(gca,'fontsize',15)
xlabel('$\theta_i$','interpreter','latex','fontsize',20)
ylabel('$\rho(\theta_i)$','interpreter','latex','fontsize',20)
title(sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
%legend({'KS model','reduced SDE'},'fontsize',15,'location','nw')

%% 
ind_i = 13;
ind_phi_i = ind_xt_a(ind_i);
figure
hold on
histogram(r1_sim.phit_c_a_traj(:,ind_i),'normalization','pdf')
histogram(r1_sde_OU.xt_a_traj(:,ind_i),'normalization','pdf')
hold off
xlim([0.6,1.25])
set(gca,'fontsize',15)
xlabel('$\hat{\theta}_i$','interpreter','latex','fontsize',20)
ylabel('$\rho(\hat{\theta}_i)$','interpreter','latex','fontsize',20)
title(sprintf('$i = %d$',ind_phi_i),'interpreter','latex','fontsize',20)
legend({'KS model','reduced SDE'},'interpreter','latex','fontsize',15)




%%
figure


%% calculate statistics of phases of synchronized oscillators from data from
%  SDE

xt_a_Delta_mu = mean(xt_a_Delta_mean,1);