close all
clear all

load('user_channel.mat','Phi_AoAw','Phi_AoDw','Phi_L','gain_w','theta_init','pd','N','M','K','L','ite');
snr_w=linspace(-10,0,5);
rate_w=zeros(1,length(snr_w));
rho=10;
for s0=1:length(snr_w)
    %%
    fprintf('snr=%d\n',snr_w(s0));
    Pt=10.^(snr_w(s0)/10);  % 最大发射功率
    rate=zeros(1,ite);
    parfor it0=1:ite
        %%
        Phi_AoA=Phi_AoAw(:,:,it0);
        Phi_AoD=Phi_AoDw(:,:,it0);
        Phi=Phi_L(:,:,:,it0);
        theta_ini=theta_init(:,:,it0);
        gaink=gain_w(:,:,:,it0);
        %%
        Gb=zeros(N,M,K);
        for k0=1:K
            Gb(:,:,k0)=ULA_fun( Phi_AoA(1,k0),N)*ULA_fun( Phi_AoD(1,k0),M)';
        end
        %%
        con=1e2;
        ratec=zeros(1,con);
        for c0=1:con
            G=zeros(N,M,K);
            for k0=1:K
                G(:,:,k0)=sqrt(rho/(rho+1)).*Gb(:,:,k0)+sqrt(1/2/(rho+1)).*(randn(N,M)+1j.*randn(N,M));
                G(:,:,k0)=pd.*G(:,:,k0);
            end
            h=zeros(N,1,K);
            for k0=1:K
                ht=zeros(N,1);
                for l0=1:L
                    ht=ht+sqrt(gaink(l0,:,k0)).*(randn(1,1)+1j.*randn(1,1)).*ULA_fun( Phi(l0,1,k0),N);
                end
                h(:,:,k0)=ht;
            end
            H=zeros(N*K,M);
            for k0=1:K
                H(N*(k0-1)+1:N*k0,:)=diag(h(:,:,k0)')*G(:,:,k0);
            end
            %%
            A=H*H';
            f0=real(theta_ini'*A*theta_ini);
            %%
            ratec(c0)=log(1+Pt*f0);
        end
        rate(it0)=mean(ratec);
        %%
    end
    rate_w(s0)=mean(rate);
end
%%
figure
plot(snr_w,rate_w,'ro-');
save('down_random_P2.mat','snr_w','rate_w','M','N');