function f = test_rate(thetaw,Pt,pd,rho,Gb,gaink,Phi,N,M,K,L)
    rng(123);
    con=1e2;
    ratec=zeros(1,con);
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
        %%
        A=H*H';
        f0=real(thetaw'*A*thetaw);
        ratec(c0)=log(1+Pt*f0);
    end
    f=mean(ratec);
end

