function [ f0 ] = fun_B( thetaw,HH,K )
    f0=0;
    for k0=1:K
        for j0=1:K
            f0=f0+real(thetaw(:,:,k0)'*HH(k0,j0)*thetaw(:,:,j0));
        end
    end
end

