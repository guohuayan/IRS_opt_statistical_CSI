close all
clear all

line=1.1;
load('converge_SGD.mat','xl','rate_L');
rate_SGD=rate_L;
load('converge_SGD_bs_10.mat','xl','rate_L');
rate_SGD_10=rate_L;
load('converge_SGD_bs_1.mat','xl','rate_L');
rate_SGD_1=rate_L;
load('converge_SGD_bs_1e3.mat','xl','rate_L');
rate_SGD_1e3=rate_L;
load('converge_SCA.mat','xl','rate_L');
figure
ht=plot(xl,rate_SGD_1e3,'ro--',xl,rate_SGD,'k+--',xl,rate_SGD_10,'ks--',xl,rate_L,'bv-',xl,rate_SGD_1,'kx--','LineWidth',line);
grid on
ylabel('R'); xlabel('Iterations');
g=legend(ht, 'Statistical CSI (SGD I=1000)','Statistical CSI (SGD I=100)',...
'Statistical CSI (SGD I=10)',...    
'Statistical CSI (Online Learning)','Statistical CSI (SGD I=1)');
set(g,'Location','southeast','FontSize',10) ;