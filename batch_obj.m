function [ f] = batch_obj( Ht_w,theta,bsize,Pt)
    f=0;
    for b0=1:bsize
        Ht=Ht_w(:,:,b0);
        H=theta'*Ht;
        w=sqrt(Pt).*H./norm(H,2);
        w=w';
        f=f+log(1+abs(theta'*Ht*w)^2);
    end
    f=f/bsize;
end

