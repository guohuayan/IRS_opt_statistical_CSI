close all
clear all

load('user_channel.mat','Phi_AoAw','Phi_AoDw','Phi_L','gain_w','theta_init','pd','N','M','K','L','ite');
snr_w=linspace(-10,0,5);
rate_w=zeros(1,length(snr_w));
rate_wv=zeros(1,length(snr_w));
rho=10;
for s0=1:length(snr_w)
    %%
    fprintf('snr=%d\n',snr_w(s0));
    Pt=10.^(snr_w(s0)/10);  % 最大发射功率
    rate=zeros(1,ite);
    ratev=zeros(1,ite);
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
        t=0;
        thetaw=theta_ini;
        f0=0;
        f1=f0;
%         A=zeros(N*K,N*K);
        d=zeros(size(thetaw));
        while(1)
            t=t+1;
%             step=t^(-0.5);
            step=t^(-1);
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
            At=H*H';
            dt=At*thetaw;
            d=(1-step).*d+step.*dt;
            %%
            f_tmp=real(thetaw'*At*thetaw);
            x=angle(d);
            thetaw=exp(1j.*x);
            f1=real(thetaw'*At*thetaw);
            %%
            if f1<f_tmp
                fprintf('error sca\n');
            end
            f1=(1-step).*f0+step.*f1;
            if abs(f1-f0)<1e-2
                break
            end
            f0=f1;
        end
        ratev(it0)=log(1+Pt*f1);
        %%
        rate(it0)=test_rate(thetaw,Pt,pd,rho,Gb,gaink,Phi,N,M,K,L);
        %%
    end
    rate_w(s0)=mean(rate);
    rate_wv(s0)=mean(ratev);
end
%%
figure
plot(snr_w,rate_w,'ro-',snr_w,rate_wv,'bo-');
save('down_Stochastic_SCA_PA_v2.mat','snr_w','rate_w','M','N','rate_wv');