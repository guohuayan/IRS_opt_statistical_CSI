close all
clear all

load('user_channel.mat','Phi_AoAw','Phi_AoDw','Phi_L','gain_w','theta_init','pd','N','M','K','L','ite');
snr_w=linspace(-10,0,5);
rate_w=zeros(1,length(snr_w));
rate_wv=zeros(1,length(snr_w));
rho=10;
for s0=1:length(snr_w)
%     s0=1;
    %%
    fprintf('snr=%d\n',snr_w(s0));
    Pt=10.^(snr_w(s0)/10);  % 最大发射功率
    rate=zeros(1,ite);
    ratev=zeros(1,ite);
    parfor it0=1:ite
%         it0=1;
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
        phiw=angle(thetaw);
        beta1=0.9;beta2=0.999;
        mt=zeros(N*K,1);
        vt=mt;
%         figure
%         line_loss = animatedline;
        while(1)
            t=t+1;
            gainkc=gaink;
            batch_size=1e2;
            grad_inner=zeros(N*K,1,batch_size);
            for b0=1:batch_size
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
                At=H*H'*Pt;
                grad_inner(:,:,b0)=2*At*thetaw./(1+thetaw'*At*thetaw).*(-1j.*conj(thetaw));
            end
            grad=real(mean(grad_inner,3));
            %% adam
            mt=beta1.*mt+(1-beta1).*grad;
            vt=beta2.*vt+(1-beta2).*grad.^2;
            step=0.1;
            upd=step.*mt./(1-beta1^t)./(sqrt(vt./(1-beta2^t))+1e-8);
            %%
            phiw=phiw+upd;
            thetaw=exp(1j.*phiw);
            if norm(upd)<1e-1 && t>1e2
                break
            end
%             addpoints(line_loss,t,norm(step.*grad,2));
%             drawnow
        end       
        rate(it0)=test_rate(thetaw,Pt,pd,rho,Gb,gaink,Phi,N,M,K,L);
        %%
    end
    rate_w(s0)=mean(rate);
    rate_wv(s0)=mean(ratev);
end
%
figure
plot(snr_w,rate_w,'ro-');
save('minibatch_adam_PA_v2.mat','snr_w','rate_w','M','N');