%%


%%
clearvars
%% Parameter setting
%setting_param_sgn = '1_a';
%setting_param_sgn = '1_a_K_7';
%setting_param_sgn = '1_a_K_1d75';
setting_param_sgn = '1_a_K_1d25';
%setting_param_sgn = '1_a_K_0d1';
%setting_param_sgn = '10_c_l_0';
%setting_param_sgn = 'lorentzian_1_a';

%setting_N_sgn = ''; % '_2N','_4N','_8N','_16N' 
setting_N_sgn = '_4N';

setting_rand_sgn = '';
%setting_rand_sgn = '_rand_5';

setting_sgn = [setting_param_sgn, setting_N_sgn, setting_rand_sgn];

%% loading the parameter settings from a file
load(['params_ks_sim_',setting_sgn,'.mat']);
%load(['params_ks_sim_1_a.mat']);
%load('params_ks_sim_10_a_l_0.mat')
%load('params_ks_sim_10_c_l_0.mat')

%% loading statistics from previous simulation of the system
%load(['data_ks_sim_',setting_sgn,'.mat']);
load(['data_ks_sim_',setting_sgn,'_a','.mat']);

Nc_l = length(Nl_ind);
Nr_l = length(Nr_ind);

%% Defining the RHS of the differential equation
%f_rhs=@(t,x) kura_saka(t,N,K,A,l,w,x);
f_rhs = @(t,x) kura_saka_alltoall(t,N,K,l,w,x);

%% Setting up for numerically simulating the system

% time
t0 = 0;
dt = 0.05;
tmax = 5e4;

tt = t0:dt:tmax;
tl = floor((tmax-t0+0.5*dt)/dt)+1;

% transient
flg_transient = true;
t1 = 5e4;

% size of interval
Delta = 2000;
tl_Delta = floor((tl-1+0.5)/Delta);

% value of time-steps at interval points
tt_Delta = (1:tl_Delta)*Delta*dt;

% setting up storage for storing data from simulation
phit_Delta = zeros(tl_Delta,N);

omegat_Delta_mean = zeros(tl_Delta,N);

rt_Delta_mean = zeros(tl_Delta,1);
rt_2_Delta_mean = zeros(tl_Delta,1);

rc_t_Delta_mean = zeros(tl_Delta,1);
rc_t_2_Delta_mean = zeros(tl_Delta,1);

rr_t_Delta_mean = zeros(tl_Delta,1);
rr_t_2_Delta_mean = zeros(tl_Delta,1);


S_t_Delta_mean = zeros(tl_Delta,1);
S_t_2_Delta_mean = zeros(tl_Delta,1);
C_t_Delta_mean = zeros(tl_Delta,1);
C_t_2_Delta_mean = zeros(tl_Delta,1);


S_t_Delta_int = zeros(tl_Delta,1);
C_t_Delta_int = zeros(tl_Delta,1);

S_t_traj = zeros(tl,1);
C_t_traj = zeros(tl,1);

rt_traj = zeros(tl,1);
psit_traj = zeros(tl,1);

rc_t_traj = zeros(tl,1);
psic_t_traj = zeros(tl,1);

rr_t_traj = zeros(tl,1);
psir_t_traj = zeros(tl,1);

%ind_phi = [1:5,10,50,100,112:116];
ind_phi = [1];
ind_phi_l = length(ind_phi);

phit_c_a_traj = zeros(tl,ind_phi_l);

Ni1_ind = 10;

f_r_j_m_t_Delta_mean = zeros(tl_Delta,length(Nr_ind),2);

f_r_j_m_t_Delta_int = zeros(tl_Delta,length(Nr_ind));

phit_c_a_Delta_mean = zeros(tl_Delta,Nc_l);
phit_c_a_2_Delta_mean = zeros(tl_Delta,Nc_l);

phit_c_a_1_Delta_mean = zeros(tl_Delta,Nc_l-1);
phit_c_a_1_2_Delta_mean = zeros(tl_Delta,Nc_l-1);

%% Numerically simulate the dynamics of the system and calculate statistics


% initial condition
phi0 = rand(1,N)*2*pi-pi;

% transient
if flg_transient == true
    [tt1,phit1] = rk_4_mod_2pi(f_rhs,[0:dt:t1],phi0,1:N,Delta);
    phii = phit1(end,:);
else
    phii = phi0;
end

% numerically simulate the dynamics of the system
tic

for Di = 1:tl_Delta
    
    tt_Di_ind = (Di-1)*Delta+1:Di*Delta+1;
    
    tt_Di = ((Di-1)*Delta:Di*Delta)*dt;
    
    [~,phit,dphit] = rk_4_mod_2pi(f_rhs,tt_Di,phii);
    
    reipsit = mean(exp(1i*phit),2);
    rt = abs(reipsit);
    psit = angle(reipsit);
    
    reipsic_t = mean(exp(1i*phit(:,Nl_ind)),2);
    rc_t = abs(reipsic_t);
    psic_t = angle(reipsic_t);
    
    reipsir_t = mean(exp(1i*phit(:,Nr_ind)),2);
    rr_t = abs(reipsir_t);
    psir_t = angle(reipsir_t);
    
    S_t = 1/N*sum(sin(phit(:,Nr_ind)-psic_t+l),2);
    C_t = 1/N*sum(cos(phit(:,Nr_ind)-psic_t+l),2);
    
    f_r_j_m_t = zeros(length(tt_Di),length(Nr_ind));
    for Ni = 1:length(Nr_ind)
        f_r_j_m_t(:,Ni) = sin(phit(:,Nr_ind(Ni))-phit(:,Ni1_ind)-l);
    end
    
    omegat = dphit;
    omegat_Delta_mean(Di,:) = mean(omegat(1:end-1,:),1);
    
    %phit_c_a_Delta_mean(Di,:) = adjust_angles(mean(adj_ang_a(phit(:,Nl_ind)-psic_t)));
    %phit_t_a_2_Delta_mean(Di,:) = mean((adj_ang_a(phit(:,Nl_ind)-psic_t) - ((phit(1,Nl_ind)-psic_t(1)) - adjust_angles(phit(1,Nl_ind)-psic_t(1)))).^2)
    
    %ang_adj_diff = (phit(1,Nl_ind)-psic_t(1)) - adjust_angles(phit(1,Nl_ind)-psic_t(1));
    %phit_c_a_Delta_mean(Di,:) = mean(adj_ang_a(phit(1:end-1,Nl_ind)-psic_t(1:end-1))-ang_adj_diff);
    %phit_c_a_2_Delta_mean(Di,:) = mean((adj_ang_a(phit(1:end-1,Nl_ind)-psic_t(1:end-1))-ang_adj_diff).^2);
    
    phit_c_a_Delta_mean(Di,:) = mean(adjust_angles(phit(1:end-1,Nl_ind) - psic_t(1:end-1)),1);
    phit_c_a_2_Delta_mean(Di,:) = mean(adjust_angles(phit(1:end-1,Nl_ind) - psic_t(1:end-1)).^2,1);
    
    % psic_t_1
    % i = 115, 116 : one oscillator that is synchronized in the full system
    % is not synchronized in the reduced system
    
    %psic_1_t = angle(mean(exp(1i*phit(:,1:115)),2));
    %psic_1_t = angle(mean(exp(1i*phit(:,Nl_ind)),2));
    
    %phit_c_a_1_Delta_mean(Di,:) = mean(adjust_angles(phit(1:end-1,1:115) - psic_1_t(1:end-1)),1);
    %phit_c_a_1_2_Delta_mean(Di,:) = mean(adjust_angles(phit(1:end-1,1:115) - psic_1_t(1:end-1)).^2,1);
    %phit_c_a_1_Delta_mean(Di,:) = mean(adjust_angles(phit(1:end-1,Nl_ind) - psic_1_t(1:end-1)),1);
    %phit_c_a_1_2_Delta_mean(Di,:) = mean(adjust_angles(phit(1:end-1,Nl_ind) - psic_1_t(1:end-1)).^2,1);
    
    S_t_traj(tt_Di_ind,:) = S_t;
    C_t_traj(tt_Di_ind,:) = C_t;
    
    rt_traj(tt_Di_ind,:) = rt;
    psit_traj(tt_Di_ind,:) = psit;
    
    rc_t_traj(tt_Di_ind,:) = rc_t;
    psic_t_traj(tt_Di_ind,:) = psic_t;
    
    rr_t_traj(tt_Di_ind,:) = rr_t;
    psir_t_traj(tt_Di_ind,:) = psir_t;
    
    phit_c_a_traj(tt_Di_ind,:) = adjust_angles(phit(:,ind_phi)-psic_t);
    
    rt_Delta_mean(Di) = mean(rt(1:end-1));
    rt_2_Delta_mean(Di) = mean(rt(1:end-1).^2);
    rc_t_Delta_mean(Di) = mean(rc_t(1:end-1));
    rc_t_2_Delta_mean(Di) = mean(rc_t(1:end-1).^2);
    rr_t_Delta_mean(Di) = mean(rr_t(1:end-1));
    rr_t_2_Delta_mean(Di) = mean(rr_t(1:end-1).^2);

    S_t_Delta_mean(Di) = mean(S_t(1:end-1));
    S_t_2_Delta_mean(Di) = mean(S_t(1:end-1).^2); 
    C_t_Delta_mean(Di,:) = mean(C_t(1:end-1));
    C_t_2_Delta_mean(Di) = mean(C_t(1:end-1).^2);
    
    f_r_j_m_t_Delta_mean(Di,:,:) = [mean(f_r_j_m_t)', mean(f_r_j_m_t.^2)'];
    
    S_t_Delta_int(Di,:) = trapz(tt_Di,S_t);
    C_t_Delta_int(Di,:) = trapz(tt_Di,C_t);
    
    f_r_j_m_t_Delta_int(Di,:) = trapz(tt_Di,f_r_j_m_t,1);
    
    phii = phit(end,:);

end

t = toc



%%
%%
wbart = mean(omegat_Delta_mean(1:end-1,:),1);
omega_c = mean(wbart(Nl_ind),2);


%%

%

%
xt_data_Delta = S_t_Delta_int;
%xt_data_Delta = C_t_Delta_int;
%xt_data_Delta = f_r_j_m_t_Delta_int(:,7);

%disp([length(xt_data_Delta)])

eta_xt_data_Delta = intervalintegrate(tt_Delta,xt_data_Delta,100);

figure
histogram(eta_xt_data_Delta(10:end),'normalization','pdf')

figure
qqplot(eta_xt_data_Delta)

[h,p] = kstest(normalize(eta_xt_data_Delta))

%
eta_f_r_j_m_t = intervalintegrate(tt_Delta,f_r_j_m_t_Delta_int,200);

figure
histogram(eta_f_r_j_m_t(:,7))

figure
qqplot(eta_f_r_j_m_t(:,7))




Delta_range = [1,10,20,50,100];
D_l = length(Delta_range);

eta_xt_data_Delta_range = cell(D_l);

for Di = 1:D_l
    Delta = Delta_range(Di);
    eta_xt_data_Delta = intervalintegrate(tt_Delta,xt_data_Delta,Delta);
    
    eta_xt_data_Delta_range{Di} = eta_xt_data_Delta;
end

figure
hold on
for Di = 1:D_l
    histogram(eta_xt_data_Delta_range{Di},'normalization','pdf')
    pause(1)
end


%%
figure
plot(S_t_Delta_mean,'.')

%%
Di_range = 1:100;
%Di_range = 101:200;

%%
disp([mean(S_t_Delta_mean(Di_range)),mean(S_t_2_Delta_mean(Di_range))])
disp([mean(C_t_Delta_mean(Di_range)),mean(C_t_2_Delta_mean(Di_range))])

%%
disp([mean(rt_Delta_mean(Di_range)),mean(rc_t_Delta_mean(Di_range)),mean(rr_t_Delta_mean(Di_range))])
disp([mean(rt_2_Delta_mean(Di_range)),mean(rc_t_2_Delta_mean(Di_range)),mean(rr_t_2_Delta_mean(Di_range))])

%%
disp([mean(phit_c_a_Delta_mean(Di_range,end))])
disp([mean(phit_c_a_2_Delta_mean(Di_range,end))])

%%
S_t_Delta_mu = mean(S_t_Delta_mean);
S_t_Delta_var = 1/(tl-2)*(tl-1)*(mean(S_t_2_Delta_mean)-mean(S_t_Delta_mean)^2);

C_t_Delta_mu = mean(C_t_Delta_mean);
C_t_Delta_var = 1/(tl-2)*(tl-1)*(mean(C_t_2_Delta_mean)-mean(C_t_Delta_mean)^2);

rt_Delta_mu = mean(rt_Delta_mean);
rt_Delta_var = 1/(tl-2)*(tl-1)*(mean(rt_2_Delta_mean) - (mean(rt_Delta_mean))^2);

rc_t_Delta_mu = mean(rc_t_Delta_mean);
rc_t_Delta_var = 1/(tl-2)*(tl-1)*(mean(rc_t_2_Delta_mean) - mean(rc_t_Delta_mean)^2);

rr_t_Delta_mu = mean(rr_t_Delta_mean);
rr_t_Delta_var = 1/(tl-2)*(tl-1)*(mean(rr_t_2_Delta_mean) - mean(rr_t_Delta_mean)^2);

phit_c_a_Delta_mu = mean(phit_c_a_Delta_mean);
phit_c_a_Delta_var = 1/(tl-2)*(tl-1)*(mean(phit_c_a_2_Delta_mean) - mean(phit_c_a_Delta_mean).^2);

phit_c_a_1_Delta_mu = mean(phit_c_a_1_Delta_mean);
phit_c_a_1_Delta_var = 1/(tl-2)*(tl-1)*(mean(phit_c_a_1_2_Delta_mean) - mean(phit_c_a_1_Delta_mean).^2);

%%
save(['data_ks_sim_Delta_',setting_sgn,'_15','.mat'],'setting_sgn','N','K','A','l','w',...
    't0','dt','tmax','flg_transient','t1','Delta','wbart','omega_c',...
    'rt_Delta_mu','rt_Delta_var','rc_t_Delta_mu','rc_t_Delta_var','rr_t_Delta_mu','rr_t_Delta_var',...
    'S_t_Delta_mu','S_t_Delta_var','C_t_Delta_mu','C_t_Delta_var',...
    'phit_c_a_Delta_mu','phit_c_a_Delta_var','phit_c_a_1_Delta_mu','phit_c_a_1_Delta_var')

%%
% save(['data_ks_sim_Delta_',setting_sgn,'_1','.mat'],'setting_sgn','N','K','A','l','w',...
%     't0','dt','tmax','t1','omegat_Delta_mean','rt_Delta_mean','rt_2_Delta_mean',...
%     'rc_t_Delta_mean','rc_t_2_Delta_mean','rr_t_Delta_mean','rr_t_2_Delta_mean',...
%     'S_t_Delta_mean','S_t_2_Delta_mean','C_t_Delta_mean','C_t_2_Delta_mean',...
%     'phit_c_a_Delta_mean','phit_c_a_2_Delta_mean')

%%
%filename = ['ks_sim_Delta_',setting_sgn,'_10','_tmax_',int2str(tmax),'.mat'];
filename = ['ks_sim_Delta_',setting_sgn,'_1','.mat'];
save(filename,'setting_sgn','N','K','A','l','w',...
    't0','dt','tmax','tt','t1','phi0','rt_traj','psit_traj','rc_t_traj','psic_t_traj',...
    'rr_t_traj','psir_t_traj','S_t_traj','C_t_traj','ind_phi','phit_c_a_traj','wbart','omega_c')