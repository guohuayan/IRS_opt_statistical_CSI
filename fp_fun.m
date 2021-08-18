function [ fout ] = fp_fun( U,v,C,x_theta )
    fout=real(2*v'*x_theta-x_theta'*U*x_theta+C);
end

