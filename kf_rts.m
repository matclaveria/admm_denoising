function x_rts = kf_rts(yf, zz, K_kf, G_rts, G1, Gp, GP)
    
    % 1. dimensions
    L = size(G1,1);
    n_atoms_frame = size(G1,2)/3;
    P = size(yf,2);
    n_frames = length(zz)/L;
    
    % 2. Z in state-space form
    if min(size(zz)) == 1
        Z_hat = zeros(2*n_atoms_frame, P);
        zf = reshape(zz, [n_atoms_frame, n_frames]);
        for p = 1 : P
            if p == P
                Z_aux = [zf(:,end); ...
                         zeros(n_atoms_frame,1)];
            else
                ind_ = 2*p-1;
                Z_aux = [zf(:, ind_  );
                         zf(:, ind_+1)];
            end
            Z_hat(:,p) = Z_aux;
        end
    else
        Z_hat = [zz, zeros(n_atoms_frame,1)];
        Z_hat = reshape(Z_hat, [2*n_atoms_frame, P]);
    end

    % 4. calculating X through KF+RTS
    X_hat = zeros(3*n_atoms_frame, P);
    X_minus_all = zeros([size(X_hat,1), P]);
    
    % kf
    for p = 1 : P

       % 1 predictive mean
        if p == 1
            X_aux = [zeros(n_atoms_frame,1); Z_hat(:,1)];
        else
            X_aux = ...
                [X_hat((2*n_atoms_frame+1) : (3*n_atoms_frame), p-1); 
                 Z_hat(:,p)];
        end
        X_minus_all(:,p) = X_aux;

        if p == 1
            G = G1;
        elseif p == P
            G = GP;
        else
            G = Gp;
        end
        
        % Kalman gain
        % K_y = K_kf(:,:,p);

        % 5 updated mean
        X_aux = X_aux + K_kf(:,:,p)*(yf(:,p) - G * X_aux );
        X_hat(:,p) = X_aux;
        
     end

    %% RTS
    X_RTS_all = 0*X_hat;
    X_RTS_all(:,P) = X_hat(:,P);
    
    for p = (P-1) : -1: 1

        X_RTS_plus = X_RTS_all(:,p+1);
        X_p = X_hat(:, p);
        X_p_minus_plus = X_minus_all(:,p+1);

        % RTS gain matrix
        % RTS_G = ;
     
        % RTS estimate
        X_RTS_all(:,p) = X_p + G_rts(:,:,p)*(X_RTS_plus - X_p_minus_plus);
        
    end
    
    %% arrange solution into a single vector
    x_rts = zz*0;
    retrieve_ind = (n_atoms_frame+1):(n_atoms_frame*3);
    ind_p = 1:(n_atoms_frame*2);
    for p = 1 : P
        if p == P
            ind_p = ind_p(1):(ind_p(1)-1+n_atoms_frame);
            retrieve_ind = (n_atoms_frame+1) : (n_atoms_frame*2);
        end
        x_rts(ind_p) = X_RTS_all(retrieve_ind, p);
        ind_p = ind_p + (n_atoms_frame*2);
    end

end