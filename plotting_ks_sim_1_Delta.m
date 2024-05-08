%% Plotting some results from the time series of simulation data

%% clear variables
clearvars

%% Setting for the parameters of the system

% setting for parameters
%setting_param_sgn = '1_a';
setting_param_sgn = '1_a_K_7';   
%setting_param_sgn = '1_a_K_1d75';
%setting_param_sgn = '1_a_K_1d25';
%setting_param_sgn = '1_a_K_0d1';
%setting_param_sgn = '10_c_l_0';

% setting for no. of oscillators
%setting_N_sgn = '_1N';
%setting_N_sgn = '_2N';
setting_N_sgn = '_64N';

% setting for random draws of intrinsic frequencies
setting_rand_sgn  = '';
%setting_rand_sgn  = '_rand_3'; 

%
setting_sgn = [setting_param_sgn, setting_N_sgn];

%% loading parameters 

param_1 = load(['D:\MATLAB\scraps\ks_3\params_ks_sim_',setting_sgn,'.mat']);
%param_1 = load([fileloc,'params_ks_sim_',setting_sgn,'.mat']);

param = param_1;

N = param.N;
K = param.K;
A = param.A;
l = param.l;
w = param.w;

%% Loading statistics
data_5 = load(['D:\MATLAB\scraps\ks_3\data_ks_sim_',setting_sgn,'_a','.mat']);
%data_5 = load([fileloc,'data_ks_sim_',setting_sgn,'_a','.mat']);

%data_5 = load(['D:\MATLAB\scraps\ks_3\data_ks_sim_',setting_sgn,'.mat']);

Nl_ind = data_5.Nl_ind;
Nr_ind = data_5.Nr_ind;

om = data_5.omega_c;

mu_s = data_5.mu_S_t;
mu_c = data_5.mu_C_t;

mu_x = [mu_s/N,mu_c/N];

Nc_l = length(Nl_ind);
Nr_l = length(Nr_ind);

%%
data_10 = load(['D:\MATLAB\scraps\ks_3\','data_ks_sim_Delta_',setting_sgn,'_15','.mat']);
omega_c = data_10.omega_c;

%% Loading result from simulation

%r1 = load(['ks_sim_',setting_sgn,'_10','.mat']);
r1 = load(['ks_sim_Delta_',setting_sgn,'_1','.mat']);


%% Load results from corresponding thermodynamic limit
%r1_thmd = load(['ks_thermodyn_num_sol_1_a.mat']);
r1_thmd = load(['ks_thmd_num_sol_',setting_param_sgn,'.mat']);

C_t_mean_sol = r1_thmd.C_t_mean_sol;
S_t_mean_sol = r1_thmd.S_t_mean_sol;

r_sol = r1_thmd.r_sol;
Om_sol = r1_thmd.Om_sol;

psic_1_sol = r1_thmd.psic_1_sol;



%% 

% histogram of r
figure
histogram(r1.rt_traj,'normalization','pdf','numbins',50)
set(gca,'fontsize',15)
xlabel('$r$','interpreter','latex','fontsize',20)
ylabel('$\rho(r)$','interpreter','latex','fontsize',20)

figure
histogram_line(r1.rt_traj,50)


% histogram of r_c
figure
histogram(r1.rc_t_traj,'normalization','pdf','numbins',50)
xlabel('$r_c$','interpreter','latex')

% 


% time-series of S_t
figure
plot(r1.tt,r1.S_t_traj)
hold on


% time series of C_t
figure
plot(r1.tt,r1.C_t_traj)
xlim([0,100])


%% calculate the mean values of S and C

rbar = mean(r1.rt_traj);
%omega_c = om;

psic_1_mean = mean(adj_ang_a(r1.psic_t_traj-r1.psit_traj));
%mean(sin(adj_ang_a(r1.psit_traj-r1.psic_t_traj)))
%sin(mean(adj_ang_a(r1.psit_traj-r1.psic_t_traj)))

%mean(cos(adj_ang_a(r1.psit_traj-r1.psic_t_traj)))
%cos(mean(adj_ang_a(r1.psit_traj-r1.psic_t_traj)))



S_t_mean_ana_sim = cos(-psic_1_mean)*1/N*sum(ka_kb(w(Nr_ind),K,rbar,omega_c));
C_t_mean_ana_sim = -sin(-psic_1_mean)*1/N*sum(ka_kb(w(Nr_ind),K,rbar,omega_c));


S_t_mean_ana_thmd = cos(-psic_1_sol)*1/N*sum(ka_kb(w(Nr_ind),K,r_sol,Om_sol));
C_t_mean_ana_thmd = -sin(-psic_1_sol)*1/N*sum(ka_kb(w(Nr_ind),K,r_sol,Om_sol));

disp([mean(r1.S_t_traj),mean(r1.C_t_traj)])

disp([S_t_mean_ana_sim,C_t_mean_ana_sim])

disp([S_t_mean_ana_thmd,C_t_mean_ana_thmd])

disp([S_t_mean_sol,C_t_mean_sol])

%%
% histogram of S_t with fitted Gaussian
x_data = r1.S_t_traj;
figure
histogram(r1.S_t_traj,'normalization','pdf','numbins',50)
%histogram_line(r1.S_t_traj)
hold on
zz = -5:0.01:5;
mu_x = mean(x_data);
var_x = var(x_data);
xx = mu_x +sqrt(var_x)*zz;
yy = 1/sqrt(2*pi*var_x)*exp(-1/2*zz.^2);
hold on
plot(xx,yy,'r','linewidth',1.9)
%[xx,yy] = fit_gaussian_dist(x_data,-5:0.01:5);
%plot(xx,yy)
%plot([S_t_mean_sol,S_t_mean_sol],[0,20],'--','color',[0.5,0.5,0.5],'linewidth',1.9)
%ax1 = plot(S_t_mean_sol,0,'o','markeredgecolor','#EDB120','markerfacecolor','#EDB120','MarkerSize',12)
%hold off
xlim([min(xx),max(xx)])
set(gca,'fontsize',15)
xlabel('$S$','interpreter','latex','fontsize',20)
ylabel('$\rho(S)$','interpreter','latex','fontsize',20)


%% histogram of C_t
figure
histogram(r1.C_t_traj,'normalization','pdf','numbins',50)
hold on
[xx,yy] = fit_gaussian_dist(r1.C_t_traj,-5:0.01:5);
plot(xx,yy,'r','linewidth',1.9)
%plot(C_t_mean_sol,0,'o','markerfacecolor','#EDB120','MarkerSize',12)
hold off
xlim([min(xx),max(xx)])
set(gca,'fontsize',15)
xlabel('$C$','interpreter','latex','fontsize',20)
ylabel('$\rho(C)$','interpreter','latex','fontsize',20)


%%
% selecting data to be plotted
%x_data = r1.rt_traj;
%x_data = r1.rc_t_traj;
%x_data = r1.S_t_traj;
x_data = r1.C_t_traj;
%%
% histogram of data with fitted Gaussian
figure
hist_with_gaussian_fit(x_data)
%xlim(1.15*[min(x_data),max(x_data)])

%
%figure

%figure

%%
% histogram of S_t as a line, with fitted Gaussian
figure
histogram_line(x_data)
hold on
gaussian_fit_line(x_data)
legend({'density from data','Gaussian fit'})

%% compute auto-covariance of S(t), C(t)

tauN0 = 0; dtauN = 1; tauNmax = 1000;
tauNt = tauN0:dtauN:tauNmax;
tauNl = length(tauNt);

dt = r1.dt;

tt_tauNt = tauNt*dt;

%tauN1 = 16e4;
tauN1 = 5e5;
t_tauN1 = tauN1*dt;
%%
% S_t_covar = x_covar(tauNt,r1.S_t);
% C_t_covar = x_covar(tauNt,r1.C_t);
% S_C_covar = x_y_covar(tauNt,r1.S_t,r1.C_t);
% C_S_covar = x_y_covar(tauNt,r1.C_t,r1.S_t);
% 
% figure
% plot(tauNt*r1.dt,S_t_covar,'.-')
% set(gca,'fontsize',15)
% xlabel('$\tau$','interpreter','latex','fontsize',20)
% ylabel('$\langle S(t+\tau)S(t)\rangle_t$','interpreter','latex','fontsize',20)
% %h = legend('')
% 
% figure
% plot(tauNt*r1.dt,C_t_covar,'.-');
% 
% figure
% plot(tauNt*r1.dt,S_C_covar,'.-')

%% calculate the auto-correlation function

%x_t = [r1.S_t(7e5:end),r1.C_t(7e5:end)];
x_t = sqrt(r1.N)*[r1.S_t_traj,r1.C_t_traj];
%x_t = [1/r1.N*r1.S_t,1/r1.N*r1.C_t];

x_t_vec_covar = x_vec_covar(tauNt,x_t(1:end,:),tauN1);

x_t_vec_covar_a = reshape(x_t_vec_covar,tauNl,4,1);

%% range of labels for plotting
label_range = {'\xi \xi ','\zeta \xi ','\xi \zeta ', '\zeta \zeta '};

%% plotting the calculated auto-correlaion function
ind_plot = 1;
figure
plot(tauNt*r1.dt,x_t_vec_covar_a(:,ind_plot),'.-','linewidth',2)
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel(['$R_{',label_range{ind_plot},'}(\tau)$'],'interpreter','latex','fontsize',20)

%legend({'$\tau_1 = 2000$','$\tau_1 = 4000$', '$\tau_1 = 8000$'},'interpreter','latex','location','northeast')

%% plotting the calculated auto-correlation function
figure
plot(tauNt*r1.dt,x_t_vec_covar_a(:,1),'.-')
hold on
%plot(tauNt*r1.dt,x_t_vec_covar(:,1,2),'.-')
%plot(tauNt*r1.dt,x_t_vec_covar(:,2,1),'.-')
%plot(tauNt*r1.dt,x_t_vec_covar_a(:,4),'.-')
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
%legend({'$R_{SS}(\tau)$','$R_{SC}(\tau)$','$R_{CS}(\tau)$','$R_{CC}(\tau)$'},'interpreter','latex')


% figure
% plot(tauNt*r1.dt,x_t_vec_covar(:,1,2),'.-')
% hold on
% plot(tauNt*r1.dt,x_t_vec_covar(:,2,1),'.-')
% plot(tauNt*r1.dt,1/2*(x_t_vec_covar(:,1,2)+x_t_vec_covar(:,2,1)),'.-')


%% auto-correlation for different length of data

% range of length of data
%tauN1_range = [2e4,4e4,8e4,16e4,32e4,64e4];
tauN1_range = [48e4,72e4,96e4];
tauN1_l = length(tauN1_range);
    
% setting up for storage
x_t_vec_covar_tauN1_range = zeros(tauNl,2,2,tauN1_l);

% calculate auto-correlation from different length of data
for tauN1_i = 1:tauN1_l
    tauN1 = tauN1_range(tauN1_i);
    
    x_t_vec_covar = x_vec_covar(tauNt,x_t,tauN1);
    
    x_t_vec_covar_tauN1_range(:,:,:,tauN1_i) = x_t_vec_covar;
end

%% plotting the auto-correlation calculated from different length of data

% range of tauN1 to be plotted
%tauN1_ind_plotting = [1:tauN1_l];
tauN1_ind_plotting = [1,2,3];

% plot auto-correlation calculated from different length of data
figure
for tauN1_i = tauN1_ind_plotting
    tauN1 = tauN1_range(tauN1_i);
    x_t_vec_covar = x_t_vec_covar_tauN1_range(:,:,:,tauN1_i);
    plot(tt_tauNt,x_t_vec_covar(:,1,1),'.-')
    hold on
end
hold off
%xlim([2.5,10])
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel('$C(\tau)$','interpreter','latex','fontsize',20)
%h = legend({'$t_1 = 4000$','$t_1 = 8000$', '$t_1 = 16000$'},'interpreter','latex','fontsize',15);


%%
% fitting the auto-correlaton function

tauNa_fit = 31;
tauNt_fit_range = [1,11:51];
tauNt_fit = tauNt(tauNt_fit_range);

tt_fit = tauNt_fit*dt;

x_t_vec_covar_a_fit = x_t_vec_covar_a(tauNt_fit_range,:);
x_t_vec_covar_t0 = x_t_vec_covar(1,:);


flg_var = 'a';

Lambda = 1000;
sgn_Lambda = '';
% Lambda = -1;
% sgn_Lambda = '_n';

% flg_sig = 'a3';
% f_fun = @(p,tt,flg_var) f_fun_OU_ac_1_a3(p,tt,flg_var);

flg_sig = 's3';
f_fun = @(p,tt,flg_var) f_fun_OU_ac_1_s3(p,tt,flg_var);
p_1 = [1,2,1,0,1];

% flg_sig = 's1';
% f_fun = @(p,tt,flg_var) f_fun_OU_ac_1_s1(p,tt,flg_var);
% p_1 = [1,2,0.3];


% tauNt_fit_range = 1:51;
% taut_data = tt(tauNt_fit_range);
% tautl_data = length(taut_data);
% x_t_autocov_data = N*x_t_autocov(tauNt_fit_range,:,:);
%     

[p_a,~,~,flg_out] = lsqnonlin(@(p) f_OU_curvefit(p,tt_fit,f_fun,flg_var,x_t_vec_covar_a_fit,Lambda,x_t_vec_covar_t0),p_1);

y_a = f_fun(p_a,tt_fit,'a');

ind_plot = 1;
figure
plot(tt_fit,x_t_vec_covar_a_fit(:,ind_plot),'b.-','linewidth',1.5,'markersize',9)
hold on
plot(tt_fit,y_a(:,ind_plot),'r-','linewidth',1.5)
set(gca,'fontsize',15)
%plot(tt_fit,y_a(:,1)/y_a(1,1)*x_t_vec_covar(1,1,1),'m-','linewidth',2)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel(['$R_{',label_range{ind_plot},'}(\tau)$'],'interpreter','latex','fontsize',20)
%title('$$','interpreter','latex')

figure
for ind_plot = 1:4
    subplot(2,2,ind_plot)
    plot(tt_fit,x_t_vec_covar_a_fit(:,ind_plot),'b.-','linewidth',1.5,'markersize',9)
    hold on
    plot(tt_fit,y_a(:,ind_plot),'r-','linewidth',1.5)
    set(gca,'fontsize',15)
    xlabel('$\tau$','interpreter','latex','fontsize',20)
    ylabel(['$R_{',label_range{ind_plot},'}(\tau)$'],'interpreter','latex','fontsize',20)
end
%plot(tt_fit,x_t_vec_covar_fit(:,1,1),'b.-','linewidth',1.5,'markersize',9)






%%


y_a_tt = f_fun(p_a,tt_tauNt,'a');

p_a_1 = p_a;
%p_a_1(1) = 0.8;
y_a_1_tt = f_fun(p_a_1,tt_tauNt,'a');

ind_plot = 1;
figure
plot(tt_tauNt,x_t_vec_covar_a(:,ind_plot))
hold on
plot(tt_tauNt,y_a_tt(:,ind_plot))
%plot(tt_tauNt,y_a_1_tt(:,1))

%xlim([0,10])

figure
for ind_plot = 1:4
    subplot(2,2,ind_plot)
    plot(tt_tauNt,x_t_vec_covar_a(:,ind_plot),'b-','linewidth',1.5,'markersize',9)
    hold on
    plot(tt_tauNt,y_a_tt(:,ind_plot),'r-','linewidth',1.5)
    plot(tt_tauNt,zeros(1,length(tt_tauNt)),'-','color',[0.5,0.5,0.5],'linewidth',1)
    xlim([0,15])
    set(gca,'fontsize',15)
    xlabel('$\tau$','interpreter','latex','fontsize',20)
    ylabel(['$R_{',label_range{ind_plot},'}(\tau)$'],'interpreter','latex','fontsize',20)
end

figure
ind_plot = 4;
plot(tt_tauNt,x_t_vec_covar_a(:,ind_plot))
hold on
plot(tt_tauNt,y_a_tt(:,ind_plot),'r-','linewidth',1.5)
plot(tt_tauNt,zeros(1,length(tt_tauNt)),'-','color',[0.5,0.5,0.5],'linewidth',1)
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)
ylabel(['$R_{',label_range{ind_plot},'}(\tau)$'],'interpreter','latex','fontsize',20)

%%
p_r_a = lsqnonlin(@(p)f_OU_curvefit(p,tt_fit,f_fun,'s',rt_vec_covar,Lambda),p_1);

figure
plot(tt_fit,rt_vec_covar)
hold on
rt_a = f_fun(p_r_a,tt_fit,flg_var);
plot(tt_fit,rt_a(:,1,1))

%% converting the fitting result to a standard form
switch flg_sig
    case 'a3'
        p_a_a3 = p_a;
    case 's3'
        Sig_a = [p_a(3),p_a(4);p_a(4),p_a(5)];
        SS_a = Sig_a*Sig_a';
        p_a_a3 = [p_a(1),p_a(2),SS_a(1,1),SS_a(1,2),SS_a(2,2)];
    case 's1'
        Sig_a = [p_a(3),0;0,p_a(3)];
        SS_a = Sig_a*Sig_a';
        p_a_a3 = [p_a(1),p_a(2),SS_a(1,1),SS_a(1,2),SS_a(2,2)];
end
%% save the result of fitting into a file
save(['fit_result_ks_sim_Delta_',setting_sgn,'_autocov_OU_',flg_sig,'_L',sgn_Lambda,'_tau1_0d5','_taumax_2d5','.mat'],'p_a','p_a_a3')

%%
xt = randn(1e3,100);

xt_vec_corr = x






%% MSQ for S(t) and C(t)

tauNt = 0:10:10000;

dt = r1.dt;

S_t_traj_msq = x_msq(tauNt,r1.S_t_traj);
C_t_traj_msq = x_msq(tauNt,r1.C_t_traj);

figure
plot(tauNt*dt,S_t_traj_msq,'linewidth',2.4)
hold on
plot(tauNt*dt,C_t_traj_msq,'linewidth',2.4)
set(gca,'fontsize',15)
xlim([0,20])
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$M(t)$','interpreter','latex','fontsize',20)
legend({'$S$','$C$'},'interpreter','latex','fontsize',15)


figure
plot(tauNt*dt,S_t_traj_msq,'linewidth',1)
set(gca,'fontsize',15)

plot_ind = 1:21;

figure
plot(tauNt(plot_ind)*dt,S_t_traj_msq(plot_ind,:),'linewidth',1)
set(gca,'fontsize',15)


%% MSQ for r(t) and r_c(t)
rt_traj_msq = x_msq(tauNt,r1.rt_traj);

rc_t_traj_msq = x_msq(tauNt,r1.rc_t_traj);

figure
plot(tauNt*dt,rt_traj_msq)
hold on
plot(tauNt*dt,rc_t_traj_msq)

%% MSQ for psi(t) and psi_c(t)

psit_traj_msq = x_msq(tauNt,adj_ang_a(r1.psit_traj-omega_c*r1.tt'));

psic_t_traj_msq = x_msq(tauNt,adj_ang_a(r1.psic_t_traj-omega_c*r1.tt'));

figure
plot(tauNt*dt,psit_traj_msq)

figure
plot(tauNt*dt,psic_t_traj_msq)
