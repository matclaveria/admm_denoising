function [A, B, Bt, X, invX, invXB, BtinvX, atoms, gabor_mask, t] = ...
    compute_matrices(window_size, window_func, rho)

    t = (-(window_size/2-1) : window_size/2) - 0.5;

    if window_func == 1
        gabor_mask = ...
            sin ( pi/2 * sin(pi/(window_size)*(t+window_size/2)).^2);
    else
        gabor_mask = 0.5 + 0.5*cos(2*pi*t/window_size);
    end


    freqs = 2*pi * (0 : (window_size/2-1) ) / window_size;
    real_atoms = (cos(freqs'*t)).*gabor_mask;
    
    freqs = 2*pi * (1 : window_size/2 ) / window_size;
    imag_atoms = sin(freqs'*t).*gabor_mask;
    all_atoms = zeros(window_size, window_size);
    
    r_ind = [1, 2:2:(window_size-1)];
    all_atoms(r_ind,:) = real_atoms;
    i_ind = [3:2:window_size, window_size];
    all_atoms(i_ind,:) = imag_atoms;

    lower_half = 1:(window_size/2);
    upper_half = (window_size/2+1):window_size;

    % block diagonal
    B_0 = all_atoms*all_atoms';
    A = B_0 + eye(size(B_0)) * rho;

    % upper block off-diagonal
    B = all_atoms(:,upper_half) * all_atoms(:,lower_half)';
    % Bt_ = all_atoms(:,lower_half) * all_atoms(:,upper_half)';

    % lower block off-diagonal (B transpose)
    Bt = B';
    
    X = meini_algorithm(A, B);
    invX = inv(X);
    A_sol = X + Bt*invX*B;
    fprintf(['Diff between original (A) and reconstructed ', ...
        '(X + Bt*X^-1*B) matrices: %e\n'], max(max(abs(A_sol - A))))
    invXB = invX * B;
    BtinvX = Bt*invX;
    
    atoms.real_atoms = real_atoms;
    atoms.imag_atoms = imag_atoms;
    atoms.all_atoms = all_atoms;
    
end