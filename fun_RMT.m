function f=fun_RMT(theta,G,D)
    Theta=diag(theta);
    f=real(trace(G'*Theta'*D*Theta*G));
end

