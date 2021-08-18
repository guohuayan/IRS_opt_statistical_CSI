close all
clear all

rng(123);
%%
Ld=path_LOS(10);
%%
noise=-170+10*log10(200*1e3);
path_d=10.^((-noise-2*Ld)/10);
%%
ite=1e2;
%%
N=20;
M=4;
K=2;
L=5;
%%
Phi_AoAw=zeros(1,K,ite);
Phi_AoDw=zeros(1,K,ite);
Phi_L=zeros(L,1,K,ite);
gain_w=zeros(L,1,K,ite);
theta_init=exp(1j.*rand(N*K,1,ite).*2.*pi);
pd=sqrt(path_d);
for j0=1:ite 
    for k0=1:K
        gain=sqrt(1/2).*(randn(L,1)+1j.*randn(L,1));
        gain=abs(gain).^2;
        gain=gain./sum(gain);
        gain_w(:,:,k0,j0)=gain;
    end
    %%
    Phi_AoAw(:,:,j0)=rand(1,K).*2.*pi;
    Phi_AoDw(:,:,j0)=rand(1,K).*2.*pi;
    Phi_L(:,:,:,j0)=rand(L,1,K).*2.*pi;
end
%%
save('user_channel.mat','Phi_AoAw','Phi_AoDw','Phi_L','gain_w','theta_init','pd','N','M','K','L','ite');









