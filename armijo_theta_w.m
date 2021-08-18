function [ f1,phi,theta ] = armijo_theta_w( grad,L_old,phiw,At,f0)
    m=0;
    rhom=0.8;
    rho0=1/L_old*100;    
    len=norm(grad,2)^2;
    sig=0.4;
    while(1)
        rho=rho0*rhom^m;
        phi=phiw-rho*grad;
        x=exp(1j.*phi);
%         f1=log(1+real(x'*At*x));
        f1=real(x'*At*x);
        if (f1-f0)>=sig*rho*len
            break
        end
        if (rho)<1/L_old/10
            break
        end
        m=m+1;
    end
    theta=x;
end

