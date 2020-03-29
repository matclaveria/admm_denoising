function [K_kf, G_rts, yf, G1, Gp, GP] = kf_rts_setup2(y_, zz, all_atoms, rho)
    
    % 1. Matrix "H"
    L = size(all_atoms,1);
    n_atoms_frame = size(all_atoms,2);
    if size(y_,2) == 1
        P = length(y_)/L;
        yf = reshape(y_, [L, P]);
    else
        P = size(y_,2);
        yf = y_;
    end
    n_frames = length(zz)/L;
    
    all_atoms_lower = [zeros(L/2,n_atoms_frame); ...
        all_atoms(1:(L/2),:)];

    all_atoms_upper = [all_atoms((L/2+1):L,:); ...
        zeros(L/2,n_atoms_frame)];
    
    G1 = [0*all_atoms, all_atoms, all_atoms_lower];
    Gp = [all_atoms_upper, all_atoms, all_atoms_lower];
    GP = [all_atoms_upper, all_atoms, 0*all_atoms];

    % 2. Kalman matrices
    eye_ = eye(n_atoms_frame);

    A = [0*eye_, 0*eye_,  1*eye_; ...
         0*eye_, 0*eye_,  0*eye_; ...
         0*eye_, 0*eye_,  0*eye_ ];

    Omega = ...
        [0*eye_, 1*eye_,  0*eye_; ...
         0*eye_, 0*eye_,  1*eye_ ];

    Q     = [0*eye_, 0*eye_        , 0*eye_         ; ...
             0*eye_, rho^-1 * eye_ , 0*eye_         ; ...
             0*eye_, 0*eye_        , rho^-1 * eye_ ];

    SigRho = [rho^-1 * eye_, 0*eye_        ; ...
              0*eye_       , rho^-1 * eye_];

    R = eye_;
    
    % 3. arrange pseudo-observations in state-space form
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
        end
    else
        Z_hat = [zz, zeros(n_atoms_frame,1)];
        Z_hat = reshape(Z_hat, [2*n_atoms_frame, P]);
    
    end
    
    % 4. calculating X through KF+RTS
    X_hat = zeros(3*n_atoms_frame, P);  
    
    P_minus_all = zeros([size(Q), P]);
    P_p_all = zeros([size(Q), P]);
    K_y_all = zeros([size(Gp'), P]);
%     X_minus_all = zeros([size(X_hat,1), P]);
    
    % kf
    for p = 1 : P

       % 1 predictive mean
%         X_aux = zeros(3*n_atoms_frame,1);
%         if p == 1
%             X_aux = [zeros(n_atoms_frame,1); Z_hat(:,1)];
%         else
%             X_aux = ...
%                 [X_hat((2*n_atoms_frame+1) : (3*n_atoms_frame), p-1); 
%                  Z_hat(:,p)];
%         end
%         X_minus_all(:,p) = X_aux;

        % 2 covariance matrix predictive mean (P)
        if p == 1
            P_minus = Q;
        else
            if p <= 20
            P_minus = A * P_p_all(:,:,p-1) * A' + Q;
            end
        end
        P_minus_all(:,:,p) = P_minus;

        % 3 S_y
        if p == 1
            G = G1;
        elseif p == P
            G = GP;
        else
            G = Gp;
        end
        if p <= 20
        S_y = G * P_minus * G' + R;
        end
        % 4 Kalman gain
        if p <= 20
        K_y = P_minus * G' * inv(S_y);
        end
        K_y_all(:,:,p) = K_y;

        % 5 updated mean
%         X_aux = X_aux + K_y*(yf(:,p) - G * X_aux );

        % 6 updated covariance matrix P
        if p <= 15
        P_y = P_minus - K_y*S_y*K_y';
        end
        % 7. save relevant values 
%         X_hat(:,p) = X_aux;
        P_p_all(:,:,p) = P_y;
        
    end
    fprintf('K OK\n')
    %% RTS
    RTS_G_all = 0*P_minus_all;
%     X_RTS_all = 0*X_hat;
    P_RTS_all = 0*P_minus_all;

%     X_RTS_all(:,P) = X_hat(:,P);
    P_RTS_all(:,:,P) = P_p_all(:,:,P);

    for p = (P-1) : -1: 1

        P_RTS_plus = P_RTS_all(:,:,p+1);
        P_p = P_p_all(:,:,p);
        P_minus_plus = P_minus_all(:,:,p+1);

%         X_RTS_plus = X_RTS_all(:,p+1);
%         X_p = X_hat(:, p);
%         X_p_minus_plus = X_minus_all(:,p+1);

        %1 RTS gain
        
        if p <= 15 || (p >= (P-8))
        RTS_G = P_p * A' *inv(P_minus_plus);
        end
        RTS_G_all(:,:,p) = RTS_G;

        %2 RTS estimate
%         X_RTS_all(:,p) = X_p + RTS_G*(X_RTS_plus - X_p_minus_plus);

        %3 RTS cov matrix
        %P_RTS_all(:,:,p) = P_p + RTS_G*(P_RTS_plus - P_minus_plus)*RTS_G';
        if p <= 15 || (p >= (P-8))
        P_RTS_all(:,:,p) = P_p + RTS_G*(P_RTS_plus - P_minus_plus)*RTS_G';
        end
    end
    fprintf('G OK\n')   
    
    %% rearrange solution in just one vector
%     x_rts = zz*0;
%     retrieve_ind = (n_atoms_frame+1):(n_atoms_frame*3);
%     ind_p = 1:(n_atoms_frame*2);
%     for p = 1 : P
%         if p == P
%             ind_p = ind_p(1):(ind_p(1)-1+n_atoms_frame);
%             retrieve_ind = (n_atoms_frame+1) : (n_atoms_frame*2);
%         end
%         x_rts(ind_p) = X_RTS_all(retrieve_ind, p);
%         ind_p = ind_p + (n_atoms_frame*2);
%     end

    K_kf = K_y_all;
    G_rts = RTS_G_all(:,1:n_atoms_frame,:);
    
end