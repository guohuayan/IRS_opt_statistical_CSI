function [theta_new]=armijo_RMT(theta,G,D)
    Theta=diag(theta);
    dk=diag(2*D*Theta*G*G');
    beta=0.5; sigma=0.2;
    m=0; mmax=1e3;
    alpha=0;
    while(m<=mmax)
        if(real(-fun_RMT(theta+beta^m.*dk,G,D))<=real(-fun_RMT(theta,G,D)-sigma*beta^m.*dk'*dk))
            alpha=beta^m;
            break;
        end
        m=m+1;
    end    
    theta_new=theta+alpha.*dk;
end

