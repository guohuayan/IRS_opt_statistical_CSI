function [ f1 ] = QCQP_SCA_square( H,Pt,theta_ini,N,K ) 
    A=H*H'*Pt;  
    phiw=angle(theta_ini);
    x_theta=theta_ini; 
    f0=log(1+real(x_theta'*A*x_theta));
    f1=f0; 
    flag=0;
    [ Ltheta ] = SCA_phi_log_step( A,N*K );
    while(1)
        %%
        at=real(x_theta'*A*x_theta);          
        %%
        grad=-real((2*A*x_theta).*(-1j.*conj(x_theta)))/(1+at);
        phiw=phiw-grad/Ltheta*(1+at);
        x_theta=exp(1j.*phiw);
        f1=log(1+real(x_theta'*A*x_theta));
        %%
%         [ f1,phiw,x_theta ] = armijo_qcqp( grad,Ltheta,phiw,A);
        %%
        if f1<f0
            fprintf('error sca\n');
        elseif f1-f0<1e-6/max(1,Pt)
            flag=flag+1;
        else
            flag=0;
        end
        if flag>=3
            break
        end
        f0=f1;
    end
end

