function [obj] = batch_fun_v2( A_bw,theta,bsize,btmp)
    obj=0;
    for b0=1:btmp
        A=A_bw(:,:,b0);
        obj=obj+log(1+real(theta'*A*theta));
    end
    obj=obj/btmp;
end

