function [ f1,theta ] = QCQP_SCA( A,theta )    
    f0=real(theta'*A*theta);
    f1=f0;
    while(1)
        x=angle(A*theta);
        theta=exp(1j.*x);
        f1=real(theta'*A*theta);
        if f1<f0
            fprintf('error sca\n');
        elseif f1-f0<1e-4
            break
        end
        f0=f1;
    end
end

