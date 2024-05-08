%% Simulate an instance of the Kuramoto-Sakaguchi model



%%  
clearvars
%% Parameters

% N,K,A,l

%N = 39; % 39 79 159 319 639 1279 2559

%N=159; K=0.3; A=ones(N); l=pi/4;

N=159; K=3; A=ones(N); l=pi/4;
%N=39; K=2.05; A=ones(N); l=pi/8;
%N=159; K=1.85; A=ones(N); l=0;   % K = 1.85;

%N = 200; K = 1; A = ones(N); l = pi/4;  


% % generate intrinsic frequencies
N_1=N+1;
x_0=((1:(N_1-1))/N_1); % mid-point rule
%x_0 = (1/(2*N)+(0:(N-1))/N); % end-point rule
%x_0 = x_0 * (1+0.05*rand(1,N));

% random draws
% seed = 5;
% rng(seed);
% x_0 = sort(rand(1,N));

% Gaussian
mu=0; sig2=1;
cdf_normal=@(x)(1/2+1/2*erf((x-mu)/sqrt(2)/sig2));
w=fsolve(@(x)(cdf_normal(x)-x_0),zeros(1,N));

% uniform
% mu = 0; ga = 1;
% w = mu + ga*(2*x_0-1);

% Lorentzian
% ga = 0.5;
% w = ga*tan(pi*x_0-pi/2);


% storing the value of the parameters into a structure
param = struct('N',N,'K',K,'A',A,'l',l,'w',w);

%%
%setting_sgn = '10_c_l_0_8N';

setting_param_sgn = '1_a_K_1d25';

%setting_param_sgn = 'lorentzian_1_a';

setting_N_sgn = '_4N';

setting_rand_sgn = '';
%setting_rand_sgn = '_rand_5';

setting_sgn = [setting_param_sgn,setting_N_sgn,setting_rand_sgn];

%% saving the parameter settings into a file
save(['params_ks_sim_',setting_sgn,'.mat'],'-struct','param')

%% loading the parameter settings from a file
%setting_sgn = '1_a';
load(['params_ks_sim_',setting_sgn,'.mat'])

%% Defining the RHS of the differential equation
%f_rhs=@(t,x) kura_saka(t,N,K,A,l,w,x);
f_rhs = @(t,x) kura_saka_alltoall(t,N,K,l,w,x);
    

%% Numerically simulate the system 

% setting for the ode solver
options = odeset('reltol',1e-6,'abstol',1e-6);
%options = odeset('reltol',1e-12,'abstol',1e-12);

% Initial conditions
phi0=rand(1,N)*2*pi-pi;
%phi0 = zeros(1,N);

% Time-step size
dt=0.05;

% Transient
flg_transient = true; % flag of whether to first simulate the transient
t1=1e3;   % length of transient time

if (flg_transient == true)
    [tt1,phit_1]=rk_4_mod_2pi(f_rhs,[0:dt:t1],phi0,1:N,200);
    %[tt1,phit_1]=ode45(f_rhs,[0:dt:t1],phi0);
    phi0_1=phit_1(end,:);
else
    phi0_1=phi0;
end

%clear('phit_1')


% Simulate the system

% times
t0 = 0;  % initial time
%dt = 0.05; % time-step size
tmax = 1e3;  % max simulation time

% defining the vector of the time-steps
tt = t0:dt:tmax;

% no. of total time-steps
tl = length(tt);

% simulate the system
tic
[~,phit,dphit] = rk_4_mod_2pi(f_rhs,[tt],phi0_1);
%[~,phit] = ode45(f_rhs,[tt],phi0_1);
t = toc


%% Computing statistics of the dynamics

% Effective frequencies
% omegat = zeros(length(tt),N);
% for ti = 1:tl
%     omegat(ti,:) = f_rhs(tt(ti),xt(ti,:))';
% end
omegat = dphit;


wbart = mean(omegat,1);

% figure
% plot(tt,mean(omegat(:,Nl_ind),2))

figure
plot(wbart,'.')

% Location of synchronized cluster
Nl_ind = find(abs(diff(wbart))<1e-3);
if ~isempty(Nl_ind)
    Nl_ind = [Nl_ind, Nl_ind(end)+1];
end

%Nl_ind = unique([Nl_ind,Nl_ind+1]);

% ind_sync = zeros(1,N);
% ind_sync(Nl_ind) = 1;
% if ~isempty(Nl_ind)
%     ind_sync(Nl_ind+1) = 1;
% end
% %ind_sync([Nl_ind,Nl_ind+1]) = 1;
% Nl_ind = find(ind_sync == 1);

if isempty(Nl_ind)
    Nl_ind = [1];
end

Nr_ind=setdiff(1:N,Nl_ind);

Nc_l = length(Nl_ind);
Nr_l = length(Nr_ind);

if ~isempty(Nl_ind)
    disp(['Na = ',num2str(min(Nl_ind)),'      ','Nb = ',num2str(max(Nl_ind))])
else
    disp('incoherent')
end


% mean effective frequency of synchrpnized oscillators
if ~isempty(Nl_ind)
    omega_c = mean(wbart(Nl_ind));
else
    omega_c = 0;
end


%% compare effective frequencies with actual intrinsic frequencies
figure
plot(wbart,'bo','markersize',7.5)
hold on
plot(w,'rx','markersize',7.5)
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
legend({'$\hat{\omega}_i$','$\omega_i$'},'interpreter','latex','fontsize',15,'location','nw')

figure
plot(wbart-w,'o','markersize',7.5)
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
ylabel('$\hat{\omega}_i-\omega_i$','interpreter','latex','fontsize',20)


%% plotting the effective frequencies of the oscillators
figure
%plot(wbart,'x','markersize',7)
plot(Nl_ind,wbart(Nl_ind),'o','markersize',7)
hold on
plot(Nr_ind,wbart(Nr_ind),'x','markersize',7)
hold on
plot([Nl_ind(end)+0.5,Nl_ind(end)+0.5],[-2,2],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
text(mean(Nl_ind([1,end])),-1.5,'$i\in\mathcal{C}$','interpreter','latex','horizontalalignment','center','fontsize',15)
text(mean(Nr_ind([1,end])),0.6,'$i\in\mathcal{R}$','interpreter','latex','horizontalalignment','center','fontsize',15)
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
ylabel('$\Omega_i$','interpreter','latex','fontsize',20)

% figure
% plot(w,wbart,'x')

%% snapshot of phases of oscillators
t_ind = 1;
figure
%plot(phit(t_ind,:),'o')
plot(Nl_ind,phit(t_ind,Nl_ind),'o','markersize',7)
hold on
plot(Nr_ind,phit(t_ind,Nr_ind),'x','markersize',7)
hold on
plot([Nl_ind(end)+0.5,Nl_ind(end)+0.5],[-pi,pi],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
ylim([-pi,pi])
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
set(gca,'fontsize',15)
xlabel('$i$','interpreter','latex','fontsize',20)
ylabel('$\phi_i$','interpreter','latex','fontsize',20)


%% sample trajectories of oscillators

% in the rest frame
figure
plot(tt,adj_ang_a(phit))
xlim([0,30])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$\phi_i(t)$','interpreter','latex','fontsize',20)

% in the co-rotating frame
figure
plot(tt,adj_ang_a(phit(:,[117,118,119,129,139,149,159])-omega_c*tt'),'linewidth',2.4)
xlim([0,50])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$\theta_i(t)$','interpreter','latex','fontsize',20)
legend({'$i=117$','$i=118$','$i=119$','$i=129$','$i=139$','$i=149$','$i=159$'},...
    'interpreter','latex','fontsize',15,'location','nw')

% in the co-rotating frame, only first few
figure
plot(tt,adj_ang_a(phit(:,[117,118,119])-omega_c*tt'),'linewidth',2.4)
xlim([0,30])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$\theta_i(t)$','interpreter','latex','fontsize',20)
legend({'$i=117$','$i=118$','$i=119$'},...
    'interpreter','latex','fontsize',15,'location','nw')


%%
% Mean-field variables

% reipsit
reipsit = mean(exp(1i*phit),2);
rt=abs(reipsit);
psit=angle(reipsit);

% figure
% plot(tt,rt)

% figure
% plot(tt,mod(psit-omega_c*tt',2*pi))

% mean-field variables restricted to the synchronized cluster
if ~isempty(Nl_ind)
    reipsic_t = mean(exp(1i*phit(:,Nl_ind)),2);
    rc_t = abs(reipsic_t);
    psic_t = angle(reipsic_t);
else
    rc_t = zeros(tl,1);
    psic_t = zeros(tl,1);
end

% figure
% plot(tt,rc_t)

% figure
% plot(tt,mod(psic_t-omega_c*tt',2*pi))

%figure
%histogram(mod(psic_t-omega_c*tt',2*pi),'normalization','pdf')

% mean-field variables of the rogue oscillators
reipsir_t = mean(exp(1i*phit(:,Nr_ind)),2);
rr_t = abs(reipsir_t);
psir_t = angle(reipsir_t);

% reipsir_1_t = reipsir_t.*exp(-1i*omega_c*tt');
% rr_1_t = abs(reipsir_1_t);
% psir_1_t = angle(reipsir_1_t);

% Na_ind = 125:159;
% za = mean(exp(1i*phit(:,Na_ind)),2);
% 
% animate_reipsit(tt,abs(za),angle(za)-omega_c*tt',1,phit(:,Na_ind)-omega_c*tt')

% figure
% plot(tt,rr_t)

%%
figure
plot(tt,adj_ang_a(psit-omega_c*tt'))
hold on
plot(tt,adj_ang_a(psic_t-omega_c*tt'))
%xlim([0,10000])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
legend({'$\psi(t)$','$\psi_c(t)$'},'interpreter','latex','fontsize',15);

figure
%plot(tt,adj_ang_a(psit-omega_c*tt'))
%hold on
plot(tt,adj_ang_a(psic_t-omega_c*tt'))
%xlim([0,10000])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$\psi_c(t)$','interpreter','latex','fontsize',20)

figure
plot(tt,adj_ang_a(psit-omega_c*tt'))
%hold on
%plot(tt,adj_ang_a(psic_t-omega_c*tt'))
%xlim([0,10000])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$\psi(t)$','interpreter','latex','fontsize',20)

% figure
% plot(tt,adj_ang_a(psit-omega_c*tt'))
% hold on
% plot(tt,adj_ang_a(psic_t-omega_c*tt')-pi)
% %xlim([0,10000])
% set(gca,'fontsize',15)
% xlabel('$t$','interpreter','latex','fontsize',20)
% legend({'$\psi(t)$','$\psi_c(t)$'},'interpreter','latex','fontsize',15);

%figure
%plot(tt,adj_ang_a())

figure
plot(tt,adj_ang_a(psit-omega_c*tt')-adj_ang_a(psic_t-omega_c*tt'))
%xlim([0,10000])
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
legend({'$\psi(t)-\psi_c(t)$'},'interpreter','latex','fontsize',15);

%%
tauNt = 0:10:1000;

psit_msq = x_msq(tauNt,adj_ang_a(psit-omega_c*tt'));
psic_t_msq = x_msq(tauNt,adj_ang_a(psic_t-omega_c*tt'));

figure
plot(tauNt*dt,psit_msq)
set(gca,'fontsize',15)
xlabel('$\tau$','interpreter','latex','fontsize',20)    
ylabel('$\textrm{MSQ}(\psi(t))$','interpreter','latex','fontsize',20)


figure
plot(tauNt*dt,psic_t_msq)
set(gca,'fontsize',15)
xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$\textrm{MSQ}(\psi_c(t))$','interpreter','latex','fontsize',20)


%%
% sum_j sin(theta_j-psic-l) for j over rogues,  cos
S_t = sum(sin(phit(:,Nr_ind)-psic_t+l),2);
C_t = sum(cos(phit(:,Nr_ind)-psic_t+l),2);

figure
plot(tt,S_t)
hold on
plot(tt,C_t)

figure
plot(tt,S_t-mean(S_t))
hold on
plot(tt,C_t-mean(C_t))

figure
plot(S_t,C_t)

figure
plot(log(abs(fftshift(fft(S_t-mean(S_t))))))

Delta = 10000;
[tt_Delta,eta_S_t_Delta] = intervalintegrate(tt,S_t,Delta);

figure
histogram(eta_S_t_Delta/sqrt(Delta*dt),'normalization','pdf')

% mean of S_t and C_t
mu_S_t = mean(S_t,1);
mu_C_t = mean(C_t,1);

% moving average 
figure
plot(tt,cumsum(S_t*dt,1)./tt')
xlim([400,1000])


x_t = S_t;
mu_x_t = data_1.mu_S_t;

y_t = C_t;
mu_y_t = data_1.mu_C_t;

tauN0 = 0; dtauN = 1; tauNmax = 100;
tauNt = tauN0:dtauN:tauNmax;
tauNl = length(tauNt);

x_corr = zeros(tauNl,1);

t2_ind = 1;
for tauNi = 1:tauNl

    tauN = tauNt(tauNi);
    x_corr(tauNi,:) = mean((x_t(t2_ind,:,:)-mu_x_t).*(x_t(t2_ind+tauN,:,:)-mu_x_t),3);
end

figure
plot(tauNt*Delta*dt,x_corr(:,1))


% influence of a rogue oscillator on a synchronized oscillator
% index of rogue oscillator
Nj1_ind = 35;

% index of synchronized oscillator
Ni1_ind = 10;

f_r_j = sin(phit(:,Nj1_ind)-phit(:,Ni1_ind)-l);

figure
histogram(f_r_j)

Delta = 1000;
[tt_Delta,eta_f_r_j_Delta] = intervalintegrate(tt,f_r_j,Delta);

figure
histogram(eta_f_r_j_Delta/(Delta*dt),'normalization','pdf')

% sum of influence of rogues on each synchronized oscillator
f_r_i_m_t = zeros(tl,length(Nl_ind));
for Nl_i = 1:length(Nl_ind)
    f_r_i_m_t(:,Nl_i) = sum(sin(phit(:,Nr_ind)-phit(:,Nl_ind(Nl_i))-l),2);
end

% mean of f_r_i_m_t
mu_f_r_i_m_t = mean(f_r_i_m_t,1);

figure
plot(tt,f_r_i_m_t(:,Ni1_ind))

figure
histogram(f_r_i_m_t(:,Ni1_ind))


%%
rc_t_mean = mean(rc_t);

za_t = 1/N*(C_t+1i*S_t)*exp(-1i*l);

z1_t = Nc_l/N*rc_t_mean + za_t;

r1_t = abs(z1_t);

figure
histogram(rt,'normalization','pdf','numbins',50)
hold on
histogram(r1_t,'normalization','pdf','numbins',50)

figure
histogram_line(rt)
hold on
histogram_line(r1_t)

figure
qqplot(rt,r1_t)

%%
mu_s_c = mean([S_t,C_t]);
cov_s_c = cov([S_t,C_t]);

x_a = mvnrnd(mu_s_c,cov_s_c,tl);

zxa_t = 1/N*(x_a(:,2)+1i*x_a(:,1))*exp(-1i*l);

zx1_t = Nc_l/N*rc_t_mean + zxa_t;

rx1_t = abs(zx1_t);

figure
histogram_line(rt)
hold on
histogram_line(rx1_t)

figure
qqplot(rt,rx1_t)

%%
mu_s_c = mean([S_t,C_t]);
var_s_c = var([S_t,C_t]);

x_c = mvnrnd(mu_s_c,var_s_c,tl);

zxc_t = 1/N*(x_c(:,2)+1i*x_c(:,1))*exp(-1i*l);

zx3_t = Nc_l/N*rc_t_mean + zxc_t;

rx3_t = abs(zx3_t);

figure
histogram_line(rt)
hold on
histogram_line(rx3_t)

figure
qqplot(rt,rx3_t)

figure
histogram_line(rx1_t)
hold on
histogram_line(rx3_t)

figure
qqplot(rx1_t,rx3_t)


%%
% histogram of (instantaneous) positions of a rogue oscillator over time
% index of the rogue oscillator
Nj1_ind = 39;

figure
histogram(mod(phit(:,Nj1_ind)-tt'*omega_c-psit(1,:)+pi,2*pi)-pi,'binlimits',[-pi,pi],'normalization','pdf')
hold on
xx = -pi:0.01*pi:pi;
plot(xx,sqrt((w(Nj1_ind)-omega_c)^2-K^2*rbar^2)/(2*pi)*1./(w(Nj1_ind)-omega_c-K*rbar*sin(xx+l)),'r-')

%%

Ni1_ind = 25;

phi_i_t_a_diff = adjust_angles(adj_ang_a(phit(:,Nl_ind)) - adj_ang_a(psic_t));
%phi_i_t_a_diff_2 = phi_i_t_a_diff.^2;
% 
% figure
% plot(tt,phi_i_t_a_diff(:,Ni1_ind))
% 
% figure
% plot(tt,phi_i_t_a_diff_2(:,Ni1_ind).^2)

%phi_i_t_a_avg = mean(adj_ang_a(phit(:,Nl_ind)) - adj_ang_a(psic_t),1);
phi_i_t_a_diff_avg = mean(phi_i_t_a_diff,1);
phi_i_t_a_diff_2_avg = mean(phi_i_t_a_diff.^2,1);

phi_i_t_a_diff_var = phi_i_t_a_diff_2_avg - phi_i_t_a_diff_avg.^2;


figure
plot(phi_i_t_a_diff_avg,'+')
hold on
phi_c_approx = asin((w(Nl_ind)-omega_c)/(K*rbar));
psi_c_approx = angle(mean(exp(1i*phi_c_approx),2));
plot(phi_c_approx - psi_c_approx,'o')

figure
plot(phi_i_t_a_diff_2_avg,'.')

figure
plot(phi_i_t_a_diff_var,'.')

figure
histogram(phi_i_t_a_diff(:,Ni1_ind),'normalization','pdf')

%%

rbar = mean(rt);

figure
plot(Nr_ind,wbart(Nr_ind),'.')
hold on
plot(Nr_ind,omega_c+sqrt((w(Nr_ind)-omega_c).^2-K^2*rbar^2),'r.')
%plot(Nr_ind,omega_c+sqrt((w(Nr_ind)-omega_c).^2-K^2*max(rt)^2),'.')
%plot(Nr_ind,omega_c+sqrt((w(Nr_ind)-omega_c).^2-K^2*min(rt)^2),'.')

figure
plot(Nr_ind,omega_c+sqrt((w(Nr_ind)-omega_c).^2-K^2*rbar^2)-wbart(Nr_ind),'.')

a = w(Nj1_ind)-omega_c;
b = K*rbar;

phi_Nj1_t_a = adj_ang_a(phit(:,Nj1_ind))-omega_c*tt'-psit(1);

figure
plot(tt,adj_ang_a(phit(:,Nj1_ind))-omega_c*tt')
hold on

c_1 = @(a,b,x0) (2/sqrt(a^2-b^2)*atan(1/sqrt(a^2-b^2)*(a*tan(x0/2)-b)));
xt_analy_approx = @(t,a,b,c1) (2*atan(b/a+sqrt(a^2-b^2)/a*tan(sqrt(a^2-b^2)/2*(t+c1))));

c1 = c_1(a,b,phit(1,Nj1_ind)+l-psit(1));
phit_Nj1_ana_approx = xt_analy_approx(tt,w(Nj1_ind)-omega_c,K*rbar,c1)-l+psit(1);

plot(tt,adj_ang_a(phit_Nj1_ana_approx),'r')

figure
plot(tt,adj_ang_a(phit(:,Nj1_ind))-omega_c*tt'-adj_ang_a(phit_Nj1_ana_approx)')

figure
plot(tt,adjust_angles(adj_ang_a(phit(:,Nj1_ind))-omega_c*tt'-adj_ang_a(phit_Nj1_ana_approx)'))

figure
plot(tt,cos(adj_ang_a(phit(:,Nj1_ind))-omega_c*tt'-adj_ang_a(phit_Nj1_ana_approx)'))


figure
plot(tt,adjust_angles(phit(:,Nj1_ind)-omega_c*tt'))
hold on
plot(tt,adjust_angles(phit_Nj1_ana_approx),'r')

figure
histogram(adjust_angles(phit(:,Nj1_ind)-omega_c*tt'),'binlimits',[-pi,pi],'normalization','pdf')
hold on
histogram(adjust_angles(phit_Nj1_ana_approx),'binlimits',[-pi,pi],'normalization','pdf')

%%
tauNt = 1:20:1000;

psic_t_1 = adj_ang_a(psic_t-omega_c*tt');

psit_1_msq = x_msq(tauNt,adj_ang_a(psit-omega_c*tt'));
psic_t_1_msq = x_msq(tauNt,adj_ang_a(psic_t-omega_c*tt'));

phi_i_t_a_msq = x_msq(tauNt,adj_ang_a(phit(:,10)-psic_t));

%%
figure
plot(tauNt*dt,psit_1_msq,'b-','linewidth',2);

figure
plot(tauNt*dt,psic_t_1_msq,'r-','linewidth',2)
xlabel('$t$','interpreter','latex','fontsize',15)
ylabel('MSD of $\psi_c$ at interval length $t$','interpreter','latex','fontsize',15)
title('$\lambda=\pi/4$, one long time series','interpreter','latex','fontsize',20)
xlim([0,10])

figure
plot(tauNt*dt,phi_i_t_a_msq)
xlabel('$t$','interpreter','latex','fontsize',15)
ylabel('MSD of $\phi_{10}-\psi_c$','interpreter','latex','fontsize',15)
title('$\lambda = \pi/4$, long time series','interpreter','latex','fontsize',20)
xlim([0,30])

%%
%t_ind_range = 1:1e4;
%t_ind_range = 1e4+1:2e4;
%t_ind_range = 1:tl;
%t_ind_range = 1:floor(tl/2);
t_ind_range = floor(tl/2)+1:tl;

phit_c_a_mean = mean(adjust_angles(adj_ang_a(phit(t_ind_range,Nl_ind)-psic_t(t_ind_range))),1);

figure
plot(phit_c_a_mean,'o')

phit_c_a_var = var(adj_ang_a(phit(t_ind_range,Nl_ind)-psic_t(t_ind_range)));

phit_c_1_var = var(adj_ang_a(phit(t_ind_range,Nl_ind)-psic_t(1)-(omega_c+1e-2)*tt(t_ind_range)'));

figure
plot(phit_c_a_var,'o')

figure
plot(phit_c_1_var,'o')

%% save the time-series from a run within a time interval
save(['ks_sim_',setting_sgn,'_time_series','_1','.mat'],'setting_sgn',...
    'N','K','A','l','w','t0','dt','tmax','tt','flg_transient','t1',...
    'phi0','phit','wbart','Nl_ind','Nr_ind','omega_c',...
    'rt','psit','rc_t','psic_t','rr_t','psir_t','S_t','C_t')

%% Save a run with initial transient
save(['ks_sim_','with_transient_',setting_sgn,'_1','.mat'],'setting_sgn','N','K','A','l','w',...
    't0','dt','tmax','tt','phi0','rt','psit')

%% Save the data from the simulation
%save(['ks_sim_',setting_sgn,'.mat'])
save(['ks_sim_',setting_sgn,'_10.mat'],'N','K','A','l','w','t0','dt','tmax','t1','tt',...
    'wbart','Nl_ind','Nr_ind','omega_c','rt','psit','rc_t','psic_t','rr_t','psir_t','S_t','C_t')

%% Save some statistics from the simulation
save(['data_ks_sim_',setting_sgn,'_a','.mat'],'t0','dt','tmax','t1','Nl_ind','Nr_ind','omega_c','mu_S_t','mu_C_t','mu_f_r_i_m_t');
%save(['data_ks_sim_',setting_sgn,'_30','.mat'],'Nl_ind','Nr_ind','omega_c','mu_S_t','mu_C_t','mu_f_r_i_m_t','phi_i_t_a_diff_avg','phi_i_t_a_diff_2_avg');   