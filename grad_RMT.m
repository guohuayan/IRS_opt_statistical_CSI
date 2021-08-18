function g=grad_RMT(theta,G,D)
    Theta=diag(theta);
    g=diag(2*D*Theta*G*G');
end