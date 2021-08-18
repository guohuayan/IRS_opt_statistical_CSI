function [ Ltheta ] = SCA_phi_log_step( A,N )
    %%
    A = (A+A')/2;
    %%
    Ut=abs(A);
    tmp1=zeros(1,N);
    for n0=1:N
        tmp=0;
        for i0=1:N
            tmp=tmp+Ut(n0,i0);
        end
        tmp1(n0)=tmp;
    end
    tmp2=zeros(1,N);
    for n0=1:N
        tmp=0;
        for i0=1:N
            if i0~=n0
                tmp=tmp+Ut(n0,i0);
            end
        end
        tmp2(n0)=tmp;
    end
    %%
    L=0;
    for n0=1:N
        for i0=1:N
            if i0==n0
                tmp=4*tmp1(n0)^2+2*tmp2(n0);
                L=L+tmp^2;
            else
                tmp=4*tmp1(n0)*tmp1(i0)+2*Ut(n0,i0);
                L=L+tmp^2;
            end
%             tmp
        end
    end
    %%
    Ltheta=sqrt(L);
end

