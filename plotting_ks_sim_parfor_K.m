%% Ploting results from simulations of the K-S Model at different K
%  asdasdasd


%% Clear variables
clearvars

%% Parameter Setting
setting_param_sgn = '1_a';

setting_N_sgn = '_4N';

setting_sgn = [setting_param_sgn,setting_N_sgn];

%% Load results from simulation
load(['ks_sim_parfor_1_',setting_sgn,'_long_transient','_phi0_rand','.mat'])

kl = length(krange);

%% load value of critical coupling strengths from simulation data
r1_kc_kg = load(['kc_kg_ks_sim_parfor_1_',setting_sgn,'_long_transient','_phi0_rand','.mat']);
Kc = r1_kc_kg.Kc;
Kg = r1_kc_kg.Kg;
%% Plotting

%% mean of r
figure
plot(krange(1:2:end),rt_mean_krange(1:2:end),'o','linewidth',1,'markersize',7)
hold on
plot([Kc,Kc],[-1,2],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
plot([Kg,Kg],[-1,2],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
plot([0,10],[0,0],'-','color',[0.5,0.5,0.5],'linewidth',1)
%plot([0,10],[1,1],'-','color',[0.5,0.5,0.5],'linewidth',1.5)
hold off
ylim([0,1.1])
set(gca,'fontsize',15)
text(Kc+0.15,0.05,'$K_c$','interpreter','latex','fontsize',17.5)
text(Kg+0.15,0.05,'$K_g$','interpreter','latex','fontsize',17.5)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$\bar{r}$','interpreter','latex','fontsize',20)

%% variance of r
rt_var_a_krange = r2t_mean_krange-rt_mean_krange.^2;

% normal scale
figure
%plot(krange,r2t_mean_krange-rt_mean_krange.^2,'o','linewidth',1)
plot(krange(1:2:end),rt_var_a_krange(1:2:end),'o','markersize',7,'linewidth',1)
hold on
plot([Kc,Kc],[-0.5,1],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
plot([Kg,Kg],[-0.5,1],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
plot(krange,zeros(1,length(krange)),'-','color',[0.5,0.5,0.5],'linewidth',1)
hold off
ylim([-0.2e-3,1.1*max(rt_var_a_krange)])
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$\textrm{Var}(r)$','interpreter','latex','fontsize',20)
ylabel('$\textrm{Var}(r)$','interpreter','latex','fontsize',20)

%figure
%plot(krange,rt_var_krange,'.')  

% log y
% figure
% plot(krange,log10(rt_var_krange),'o')


%% range of K correponding to synchronization
ind_K_sync = find(N12_krange(:,2)>0);

%% location of synchronized cluster (in i)
figure
plot(krange(ind_K_sync),N12_krange(ind_K_sync,:),'.')
xlim([min(krange),max(krange)])
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$i_a,i_b$','interpreter','latex','fontsize',20)

%% location of synchronized cluster (in w_i)
figure
plot(krange(ind_K_sync),w(N12_krange(ind_K_sync,:)),'.')
xlim([min(krange),max(krange)])
ylim([-3,3])
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$\omega_a,\omega_b$','interpreter','latex','fontsize',20)

%% relative position of synchronized cluster (in i/N)
figure
plot(krange(ind_K_sync),N12_krange(ind_K_sync,1)/N,'b^','linewidth',0.6)
hold on
plot(krange(ind_K_sync),N12_krange(ind_K_sync,2)/N,'rv','linewidth',0.6)
plot(krange,ones(1,length(krange))*0,'--','color',[0.5,0.5,0.5])
plot(krange,ones(1,length(krange))*1,'--','color',[0.5,0.5,0.5])
plot([Kc,Kc],[-1,2],'--','color',[0.5,0.5,0.5])
plot([Kg,Kg],[-1,2],'--','color',[0.5,0.5,0.5])
hold off
xlim([min(krange),max(krange)])
ylim([-0.1,1.1])
h = legend({'$i_a/N$','$i_b/N$'},'interpreter','latex','fontsize',15,'location','east');
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
ylabel('$i_a/N,i_b/N$','interpreter','latex','fontsize',20)


%% relative size of synchronized cluster

% calculate size of synchronized cluster

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

% plotting the relative size of the synchronized cluster
figure
plot(krange(1:2:end),Nl_l_krange(1:2:end)/N,'o-','markersize',7,'linewidth',1)
hold on
plot(krange(1:2:end),Nr_l_krange(1:2:end)/N,'x-','markersize',7,'linewidth',1)
%plot(krange,ones(1,length(krange))*0,'--','color',[0.5,0.5,0.5],'linewidth',1.5)
%plot(krange,ones(1,length(krange))*1,'--','color',[0.5,0.5,0.5],'linewidth',1.5)
plot([Kc,Kc],[-1,2],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
plot([Kg,Kg],[-1,2],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
hold off
ylim([-0.1,1.1])
legend({'$N_c/N$','$N_r/N$'},'interpreter','latex','fontsize',15,'location','east')
set(gca,'fontsize',15)
xlabel('$K$','interpreter','latex','fontsize',20)
%ylabel('$n_c$','interpreter','latex','fontsize',20)


%% mean and var of r_c
figure
plot(krange,rc_t_mean_krange,'o')

figure
plot(krange,rc_t_var_krange,'o')


%% mean and var of r_r
figure
plot(krange,rr_t_mean_krange,'o')

figure
plot(krange,rr_t_var_krange,'o')


%% Mean of S and C
figure
plot(krange,S_t_mean_krange,'o-')
hold on
plot(krange,C_t_mean_krange,'o-')

%% Variance of S and C

figure
plot(krange,S_t_var_krange,'o')
hold on
plot(krange,C_t_var_krange,'o')