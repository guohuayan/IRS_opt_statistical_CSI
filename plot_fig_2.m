close all
clear all

line=1.1;
load('down_random_P2.mat','snr_w','rate_w','M','N');
figure
plot(snr_w,rate_w,'ks-','LineWidth',line);
grid on
hold on
load('down_Stochastic_SCA_PA_v2.mat','snr_w','rate_w','M','N','rate_wv');
plot(snr_w,rate_w,'bv-','LineWidth',line);
rate_u=rate_w;
load('minibatch_adam_PA_v2.mat','snr_w','rate_w','M','N');
plot(snr_w,rate_w,'m^--','LineWidth',line);
load('perfect_SCA_P2.mat','snr_w','rate_w','M','N');
plot(snr_w,rate_w,'ro-','LineWidth',line);
%%
xlim([-10,0]);
%%
ylim([3,9]);
ylabel('R'); xlabel('P_T (dBm)');
it=10;
ht=plot(0:it,-1.*ones(1,it+1),'ro-',0:it,-1.*ones(1,it+1),'m^--',...
    0:it,-1.*ones(1,it+1),'bv-',0:it,-1.*ones(1,it+1),'ks-',...
    'linewidth',1);
g=legend(ht, 'Perfect CSI','Statistical CSI (SGD)','Statistical CSI (Online Learning)',...
    'Random Phase');
set(g,'Location','southeast','FontSize',10) ;