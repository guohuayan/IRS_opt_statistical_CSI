close all
clear all

rng(123);
%%
load('user_channel.mat','Phi_AoAw','Phi_AoDw','Phi_L','gain_w','theta_init','pd','N','M','K','L','ite');
con=1e2;
rho=10;
A_test=zeros(N*K,N*K,con,ite);
for it0=1:ite
    Phi_AoA=Phi_AoAw(:,:,it0);
    Phi_AoD=Phi_AoDw(:,:,it0);
    Phi=Phi_L(:,:,:,it0);
    gaink=gain_w(:,:,:,it0);
    Gb=zeros(N,M,K);
    for k0=1:K
        Gb(:,:,k0)=ULA_fun( Phi_AoA(1,k0),N)*ULA_fun( Phi_AoD(1,k0),M)';
    end
    for c0=1:con
        gainkc=gaink;
        G=zeros(N,M,K);
        for k0=1:K
            G(:,:,k0)=sqrt(rho/(rho+1)).*Gb(:,:,k0)+sqrt(1/2/(rho+1)).*(randn(N,M)+1j.*randn(N,M));
            G(:,:,k0)=pd.*G(:,:,k0);
        end
        h=zeros(N,1,K);
        for k0=1:K
            ht=zeros(N,1);
            for l0=1:L
                ht=ht+sqrt(gainkc(l0,:,k0)).*(randn(1,1)+1j.*randn(1,1)).*ULA_fun( Phi(l0,1,k0),N);
            end
            h(:,:,k0)=ht;
        end
        H=zeros(N*K,M);
        for k0=1:K
            H(N*(k0-1)+1:N*k0,:)=diag(h(:,:,k0)')*G(:,:,k0);
        end
        A=H*H';
        A_test(:,:,c0,ite)=A;
    end
end
save('channel_test_rho_10.mat','A_test','con');