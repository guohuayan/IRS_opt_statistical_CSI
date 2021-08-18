function [ f1,phi,x_theta ] = armijo_batch( grad,Ltheta,phiw,A_bw,bsize,f0)
    m=0;
    rhom=0.8;
    rho0=1/Ltheta*100;    
    len=norm(grad,2)^2;
    sig=0.4;
    x_theta=exp(1j.*phiw);
%     f0=fp_fun( Ut,vt,ct,x_theta );
%     f0=log(1+real(x_theta'*A*x_theta));
    while(1)
        rho=rho0*rhom^m;
        phi=phiw-rho*grad;
        x_theta=exp(1j.*phi);
%         f1=log(1+real(x'*At*x));
%         f1=fp_fun( Ut,vt,ct,x_theta );
        f1=batch_fun( A_bw,x_theta,bsize);
        if (f1-f0)>=sig*rho*len
            break
        end
        if (rho)<1/Ltheta/1.1
            break
        end
        m=m+1;
    end
end

