function [obj] = batch_fun( A_bw,theta,bsize)
    obj=0;
    for b0=1:bsize
        A=A_bw(:,:,b0);
        obj=obj+log(1+real(theta'*A*theta));
    end
    obj=obj/bsize;
end

