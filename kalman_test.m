%% (load signal, calculate matrices)
admm_audio
tests = 0;
lazy_admm_loop;
%%
n_atoms_frame = size(all_atoms,2);

all_atoms_lower = [zeros(window_size/2,n_atoms_frame); ...
    all_atoms(1:(window_size/2),:)];

all_atoms_upper = [all_atoms((window_size/2+1):window_size,:); ...
    zeros(window_size/2,n_atoms_frame)];

G1 = [0*all_atoms, all_atoms, all_atoms_lower];
Gp = [all_atoms_upper, all_atoms, all_atoms_lower];
GP = [all_atoms_upper, all_atoms, 0*all_atoms];

if tests
load('results_11s_harpsic.mat')
%% test 1
c_sample = c_k(:,100:102);
y_p = Gp*c_sample(:);
close
figure
plot(y_p)
hold on
frame_number = 51;
init_ = (frame_number-1)*window_size+1;
last_ = (frame_number)*window_size;
plot(y_synth( init_: last_ ));


%% test 2
c_sample = c_k(:,1:2);
c_sample = [rand(window_size,1); c_sample(:)];
close all
% plot(c_sample);
frame_number = 1;
init_ = (frame_number-1)*window_size+1;
last_ = (frame_number)*window_size;
plot(y_synth( init_: last_ ))
hold on
plot(G1*c_sample(:),'.')

end

%% test 3:
c_k = malyshev_alg(B, BtinvX, invX, M, mu_tilde, n_frames);
ck_ = c_k(:);
P = length(y)/window_size;
X_hat = zeros(3*n_atoms_frame, P);
Z_hat = zeros(2*n_atoms_frame, P);
% building pseudo-obs
Z = v_k-eta_k/rho;

for p = 1 : P
    if p == P
        Z_aux = [Z(:,end); ...
                 zeros(n_atoms_frame,1)];
    else
        ind_ = 2*p-1;
        Z_aux = [Z(:, ind_  );
                 Z(:, ind_+1)];
    end
    Z_hat(:,p) = Z_aux;
end


%% test 4
close all
plot(Z(:))
hold on
plot(Z_hat((window_size+1):end,1),'.')
frame_number = n_frames;
init_ = (frame_number-2)*n_atoms_frame+1; 
fin_  = (frame_number)*n_atoms_frame;
plot(init_:fin_, ...
    Z_hat(1:(2*window_size),end),'.')

chunk_number = 20;
frame_number = chunk_number*2 - 1;
init_ = (frame_number-2)*n_atoms_frame+1; 
fin_  = (frame_number+1)*n_atoms_frame;
plot(init_:fin_, ...
    Z_hat(:, chunk_number),'.')
return


%% Kalman
eye_ = eye(n_atoms_frame);

A = [0*eye_, 0*eye_,  1*eye_; ...
     0*eye_, 0*eye_,  0*eye_; ...
     0*eye_, 0*eye_,  0*eye_ ];
Omega = ...
    [0*eye_, 1*eye_,  0*eye_; ...
     0*eye_, 0*eye_,  1*eye_ ];

Q     = [0*eye_, 0*eye_       , 0*eye_         ; ...
        0*eye_, rho^-1 * eye_ , 0*eye_         ; ...
        0*eye_, 0*eye_        , rho^-1 * eye_ ];

SigRho = [rho^-1 * eye_, 0*eye_        ; ...
          0*eye_       , rho^-1 * eye_];
      

%%
P_minus_all = zeros([size(Q), P]);
% P_y_all = zeros([size(Q), P]);
P_p_all = zeros([size(Q), P]);
K_y_all = zeros([size(Gp'), P]);
K_z_all = zeros([size(Omega'), P]);
X_minus_all = zeros([size(X_hat,1), P]);

%%
for p = 1 : P
    
   % 1 predictive mean
    X_aux = zeros(3*n_atoms_frame,1);
    if p == 1
        X_aux = Z_hat(:,1);
    else
        X_aux(1:n_atoms_frame,1) = X_hat((2*n_atoms_frame+1) : ...
            (3*n_atoms_frame),p-1);
    end
    X_minus_all(:,p) = X_aux;
    
    % 2 covariance matrix predictive mean (P)
    if p == 1
        P_minus = Q;
    else
        P_minus = A * P_prev * A' + Q;
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
    S_y = G * P_minus * G' + eye_;
    
    % 4 Kalman gain
    K_y = P_minus * G' * inv(S_y);
    K_y_all(:,:,p) = K_y;
    
    % 5 updated mean (without pseudo-measurement)
    X_aux = X_aux + K_y*(y_folded(:,p) - G * X_aux );
    
    % 6 updated covariance matrix P
    P_y = P_minus - K_y*S_y*K_y';
    
    % 7 S_z 
    S_z = Omega * P_y * Omega' + SigRho;
    
    % 8 Kalman gain (updated with pseudo-measurement)
    K_z = P_y * Omega' * inv(S_z);
    K_z_all(:,:,p) = K_z;
    
    % 9 mean value (final)
    X_aux = X_aux + K_y*(y_folded(:,p) - Omega * X_aux );
    X_hat(:,p) = X_aux;
    
    % 10 cov matrix (final)
    P_p = P_y - K_z * S_z * K_z';
    P_p_all(:,:,p) = P_p;
 end
 
% %% RTS
RTS_G_all = 0*P_minus_all;
X_RTS_all = 0*X_hat;
P_RTS_all = 0*P_minus_all;

X_RTS_all(:,P) = X_hat(:,P);
P_RTS_all(:,:,P) = P_p_all(:,:,P);

for p = (P-1) : -1: 1
    
    P_RTS_plus = P_RTS_all(:,:,p+1) 
    P_p = P_p_all(:,:,p);
    P_minus_plus = P_minus_all(:,:,p+1);
    
    X_RTS_plus = X_RTS_all(:,p+1);
    X_p = X_hat(:, p);
    X_p_minus_plus = X_minus_all(:,p+1);
    
    %1 RTS gain
    RTS_G = P_p * A \P_p_plus;
    RTS_G_all(:,:,p) = RTS_G;
    
    %2 RTS estimate
    X_RTS_all(:,p) = X_p + RTS_G*(X_RTS_plus - X_p_minus_plus);
    
    %3 RTS cov matrix
    P_RTS_all(:,:,p) = P_p + RTS_G*(P_RTS_plus - P_minus_plus)*RTS_G';

end