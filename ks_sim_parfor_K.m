%%


%% Clear variables
clearvars

%% Setting up
%krange = 1.51:0.01:2.5;  % kc
%krange = 7.91:0.01:8.90;  % kg
%range = 1:0.01:2;
krange = 0:0.05:5;
kl = length(krange);

% parameters
setting_sgn = '1_a_4N';
%setting_sgn = '10_c_l_0';
params_1 = load(['params_ks_sim_',setting_sgn,'.mat']);
N = params_1.N;
A = params_1.A;
w = params_1.w;
l = params_1.l;

% N = 200;
% A = ones(N);
% l = 0;

% x_1 = ((1:N)-1/2)/N;
% x_1 = (1:N)/(N+1);

% generate w_i according to a uniform g(w) with mean mu and half-width ga
% mu = 0; ga = 1;
% w = 2*ga*x_1 - ga;

% generate w_i according to a Lorentzian distribution
%mu = 0; ga = 0.5;
%w = mu + ga*tan(pi*(x_1-1/2));


% time
t0 = 0; dt = 0.05; tmax = 5e4;
tt = t0:dt:tmax;
tl = length(tt);

flg_transient = 'sim';
t1 = 5e4;

% storage for simulation

wbar_krange = zeros(kl,N);
rbar_krange = zeros(kl,1);

N12_krange = zeros(kl,2);

omega_c_krange = zeros(kl,1);

% rt_krange = zeros(tl,1,kl);
% psit_krange = zeros(tl,1,kl);
% 
% rc_t_krange = zeros(tl,1,kl);
% psic_t_krange = zeros(tl,1,kl);
% 
% rr_t_krange = zeros(tl,1,kl);
% psir_t_krange = zeros(tl,1,kl);

rt_mean_krange = zeros(kl,1);
r2t_mean_krange = zeros(kl,1);
rt_var_krange = zeros(kl,1);

rc_t_mean_krange = zeros(kl,1);
rc_2_t_mean_krange = zeros(kl,1);
rc_t_var_krange = zeros(kl,1);

rr_t_mean_krange = zeros(kl,1);
rr_2_t_mean_krange = zeros(kl,1);
rr_t_var_krange = zeros(kl,1);

S_t_mean_krange = zeros(kl,1);
S_2_t_mean_krange = zeros(kl,1);
S_t_var_krange = zeros(kl,1);

C_t_mean_krange = zeros(kl,1);
C_2_t_mean_krange = zeros(kl,1);
C_t_var_krange = zeros(kl,1);

%% Simulate the system at different K


tic
parfor ki = 1:kl
    
    % load the value of K
    K = krange(ki);
    
    disp(['K=',num2str(K)])
    
    % defining the rhs of the system
    %f_rhs = @(t,x) (kura_saka(t,N,K,A,l,w,x));
    f_rhs = @(t,x) (kura_saka_alltoall(t,N,K,l,w,x));
    
    % simulate the system
    
    % initial condition
    phi0 = rand(1,N)*2*pi-pi;
    %phi0 = zeros(1,N);
    
    % transient
    if (strcmp(flg_transient,"sim"))
        [~,phi_t1] = rk_4_mod_2pi(f_rhs,[0:dt:t1],phi0);
        phi0_1 = phi_t1(end,:);
    else
        phi0_1 = phi0;
    end
                
    % simulate the system
    [~,phit,dphit] = rk_4_mod_2pi(f_rhs,tt,phi0_1);
    
    % calculate statistics from the time series
    
    % effective frequencies
    omegat = dphit;
    wbar = mean(omegat,1);
    
    % location of synchronized cluster
    Nl_ind_a = find(abs(diff(wbar))<1e-3);
    if ~isempty(Nl_ind_a)
        Nl_ind = [Nl_ind_a,Nl_ind_a(end)+1];
        Nr_ind = setdiff(1:N,Nl_ind);
    else
        Nl_ind = [];
        Nr_ind = 1:N;
    end
    
    % mean frequency within the synchronized cluster
    if ~isempty(Nl_ind_a)
        omega_c = mean(wbar(Nl_ind));
    else
        omega_c = 0;
    end
    
    % mean field variables
    reipsit = mean(exp(1i*phit),2);
    rt = abs(reipsit);
    psit = angle(reipsit);
    
    if ~isempty(Nl_ind_a)
        reipsic_t = mean(exp(1i*phit(:,Nl_ind)),2);
        rc_t = abs(reipsic_t);
        psic_t = angle(reipsic_t);
    else
        rc_t = zeros(tl,1);
        psic_t = zeros(tl,1);
    end
    
    if ~isempty(Nr_ind)
        reipsir_t = mean(exp(1i*phit(:,Nr_ind)),2);
        rr_t = abs(reipsir_t);
        psir_t = angle(reipsir_t);
    else
        rr_t = zeros(tl,1);
        psir_t = zeros(tl,1);
    end
    
    if ~isempty(Nr_ind) && ~isempty(Nl_ind_a)
        S_t = 1/N*sum(sin(phit(:,Nr_ind)-psic_t+l),2);
        C_t = 1/N*sum(cos(phit(:,Nr_ind)-psic_t+l),2);
    else
        S_t = zeros(tl,1);
        C_t = zeros(tl,1);
    end
    
    % storing
    if ~isempty(Nl_ind_a)
        N12_krange(ki,:) = [Nl_ind_a(1),Nl_ind_a(end)+1];
    else
        N12_krange(ki,:) = [-1,-1];
    end
    
    omega_c_krange(ki,:) = omega_c;
    
    wbar_krange(ki,:) = wbar;
    
    rt_mean_krange(ki,1) = mean(rt);
    r2t_mean_krange(ki,1) = mean(rt.^2);
    rt_var_krange(ki,1) = var(rt);
    
    rc_t_mean_krange(ki,1) = mean(rc_t);
    rc_2_t_mean_krange(ki,1) = mean(rc_t.^2);
    rc_t_var_krange(ki,1) = var(rc_t);
       
    rr_t_mean_krange(ki,1) = mean(rr_t);
    rr_2_t_mean_krange(ki,1) = mean(rr_t.^2);
    rr_t_var_krange(ki,1) = var(rr_t);
       
    S_t_mean_krange(ki,1) = mean(S_t);
    S_2_t_mean_krange(ki,1) = mean(S_t.^2);
    S_t_var_krange(ki,1) = var(S_t);
        
    C_t_mean_krange(ki,1) = mean(C_t);
    C_2_t_mean_krange(ki,1) = mean(C_t.^2);
    C_t_var_krange(ki,1) = var(C_t);
          
    %rt_krange(:,1,ki) = rt;
    
    %rc_t_krange(:,1,ki) = rc_t;
    
    % fitting the auto-covariance of S(t) and C(t) against analytical
    % auto-covariance function for 2d OU process
    
    
    
    
end
t = toc

%% Load results
%load('')

%% Plotting results
figure
plot(krange,rt_mean_krange,'.')

%rt_var_krange = r2t_mean_krange-rt_mean_krange.^2;
rt_var_a_krange = r2t_mean_krange-rt_mean_krange.^2;

figure
plot(krange,r2t_mean_krange-rt_mean_krange.^2,'.') 

figure
plot(krange,rt_var_krange,'.')  

figure
plot(krange,log(rt_var_krange),'.')

figure
plot(log(krange),log(rt_var_krange),'.')

figure
semilogy(krange,rt_var_krange,'.')

figure
loglog(krange,rt_var_krange,'.')

figure
plot(krange,rt_mean_krange,'.')
hold on
plot(krange,rt_mean_krange-sqrt(rt_var_krange))
plot(krange,rt_mean_krange+sqrt(rt_var_krange))


ind_K_sync = find(N12_krange(:,2)>0);

figure
plot(krange(ind_K_sync),N12_krange(ind_K_sync,:),'.')
xlim([min(krange),max(krange)])
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$i_a,i_b$','interpreter','latex','fontsize',20)


figure
plot(krange(ind_K_sync),w(N12_krange(ind_K_sync,:)),'.')
xlim([min(krange),max(krange)])
ylim([-3,3])
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$\omega_a,\omega_b$','interpreter','latex','fontsize',20)

figure
plot(krange(ind_K_sync),N12_krange(ind_K_sync,1)/N,'bv','linewidth',0.6)
hold on
plot(krange(ind_K_sync),N12_krange(ind_K_sync,2)/N,'r^','linewidth',0.6)
xlim([min(krange),max(krange)])
ylim([-0.1,1.1])
h = legend({'$i_a/N$','$i_b/N$'},'interpreter','latex','fontsize',15,'location','east');
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
%ylabel('$i_a,i_b$','interpreter','latex','fontsize',20)


%
figure
plot(krange(ind_K_sync),omega_c_krange(ind_K_sync),'.')
hold on
plot(krange(ind_K_sync),omega_c_krange(ind_K_sync)-krange(ind_K_sync)'.*rt_mean_krange(ind_K_sync),'.')
plot(krange(ind_K_sync),omega_c_krange(ind_K_sync)+krange(ind_K_sync)'.*rt_mean_krange(ind_K_sync),'.')
plot(krange(ind_K_sync),ones(1,length(ind_K_sync))*(-ga))
plot(krange(ind_K_sync),ones(1,length(ind_K_sync))*ga)

%
kl = length(krange);

Nl_l_krange = zeros(kl,1);
Nr_l_krange = zeros(kl,1);

for ki = 1:kl
    if N12_krange(ki,2)> 0
        Nl_l_krange(ki,1) = N12_krange(ki,2)-N12_krange(ki,1)+1;
        Nr_l_krange(ki,1) = N - Nl_l_krange(ki,1);
    else
        Nl_l_krange(ki,1) = 0;
        Nr_l_krange(ki,1) = N;
    end
end

figure
plot(krange,Nl_l_krange/N,'bo-','linewidth',1)
hold on
plot(krange,Nr_l_krange/N,'rx-','linewidth',1)
ylim([-0.1,1.1])
legend({'$N_c/N$','$N_r/N$'},'interpreter','latex','fontsize',15,'location','east')
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
%ylabel('$N_c/N, N_r/N$','interpreter','latex','fontsize',20)



figure
plot(krange,omega_c_krange,'.')



figure
plot(krange,rt_mean_krange,'bo-','linewidth',0.6)
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$\bar{r}$','interpreter','latex','fontsize',20)


figure
plot(krange,rt_var_krange,'o')
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$\textrm{Var}(r(t))$','interpreter','latex','fontsize',20)


%% Estimating location of synchronized cluster from data
N12_a_krange = zeros(kl,2);
indic_sync_a_krange = false(kl,N);

for ki = 1:kl
    wbar = wbar_krange(ki,:);
    Nl_ind_1 = find(abs(diff(wbar))<1e-5);
    indic_sync = zeros(1,N);
    indic_sync(Nl_ind_1) = 1;
    if ~isempty(Nl_ind_1)
        indic_sync(Nl_ind_1+1) = 1;
    end
    Nl_ind = find(indic_sync>0.5);
    indic_sync_a_krange(ki,:) = indic_sync;
    
    if ~isempty(Nl_ind_1)
        N12_a_krange(ki,1) = min(Nl_ind);
        N12_a_krange(ki,2) = max(Nl_ind);
    else
        N12_a_krange(ki,:) = [-1,-1];
    end
end

omega_c_a_krange = zeros(kl,1);

for ki = 1:kl
    if N12_a_krange(ki,2)>0.5
        Nl_ind = N12_a_krange(ki,1):1:N12_a_krange(ki,2);
        omega_c_a = mean(wbar_krange(ki,Nl_ind),2);
        omega_c_a_krange(ki) = omega_c_a;
    end
end

%% plotting estimated cluster location
ind_K_sync_a = find(N12_a_krange(:,2)>0.5);

figure
plot(krange(ind_K_sync_a),N12_a_krange(ind_K_sync_a,:),'.')

figure
plot(krange,omega_c_a_krange,'.')


%% saving data
save(['ks_sim_parfor_1_uniform_N_200_l_0','.mat'],'N','A','l','w','krange','t0','dt','t1','flg_transient','tmax',...
   'wbar_krange','N12_krange','omega_c_krange','rt_mean_krange','r2t_mean_krange','rt_var_krange','rc_t_mean_krange','rc_2_t_mean_krange','rc_t_var_krange',...
   'rr_t_mean_krange','rr_2_t_mean_krange','rr_t_var_krange','S_t_mean_krange','S_2_t_mean_krange','S_t_var_krange','C_t_mean_krange','C_2_t_mean_krange','C_t_var_krange')