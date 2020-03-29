function v = shrinking_op(c_k, eta_k, lambda, rho, ind_r, ind_i)

    r = c_k + eta_k/rho;
%     r = c_k + eta_k;
    
    norm_r = sqrt([r(1,:).^2;
        r(ind_r,:).^2 + r(ind_i,:).^2;
        r(end,:).^2]);
    
    magnitude = norm_r - lambda./rho;
    magnitude( magnitude <= 0 ) = 0;
    
    aux_norm = [norm_r(1,:); 
        kron(norm_r(2:(end-1),:), [1;1]); 
        norm_r(end,:)];
    aux_mag = [magnitude(1,:); 
        kron(magnitude(2:(end-1),:), [1;1]); 
        magnitude(end,:)];
    
    v = aux_mag.*(r./aux_norm);

end