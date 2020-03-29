%%
admm_audio
opt_n_chunks = 6;
opt_n_frames = 2*opt_n_chunks-1;

%% Build linsys matrix
[T_, N_] = optional_matrices(A, B, Bt, X, window_size, opt_n_frames);

%% Build vector
f = 1:(opt_n_chunks*window_size);
dummyfreq = pi/180;
dummyfreq2 = pi/1700;
y_ = sin(f'*dummyfreq2).*(cos(f'*dummyfreq));
y_fol = reshape(y_, [window_size, length(y_)/window_size]);
GtY = g_t_product(all_atoms_T, y_fol);
Gty = GtY(:);
rng(10)
zz = 0.1*randn(opt_n_frames*window_size,1);
rng('shuffle')

f = Gty + rho*zz;

%% Solution using built-in linear solver
x_solver = T_\f;
return

%% 
[K_kf, G_rts, yf,  G1, Gp, GP] = kf_rts_setup(y_, zz, all_atoms, rho);
x_kalman = kf_rts(yf, zz, K_kf, G_rts, G1, Gp, GP);

%%
close all
plot(x_solver)
hold on
plot(x_kalman,'.')

%%
zz_ = 0.3*(cos(zz) - 0.995);
% plot(zz_)

%%
y_fol_ = randn(size(y_fol));
Gty_ = g_t_product(all_atoms_T, y_fol_);
Gty_ = Gty_(:);
f_ = Gty_ + rho*zz_;
x_solver2 = T_\f_;
x_kalman2 = kf_rts(y_fol_, zz_, K_kf, G_rts, G1, Gp, GP);

%%
close all
plot(x_solver2)
hold on
plot(x_kalman2,'.')

%%
close all
plot(x_solver2 - x_kalman2)