%% set basic ADMM params
% clear; 
close all
rho = 5;

%% load data

% selected:
% 1 harpsichord
% 2 glockenspiel
% 3 jazzTrio
selected = 2;
offset   = 0;  % seconds
duration = 2.65;  % seconds
% duration = 2.0;  % seconds
snrlvl   = 5; % SNR of the artificial noise
fs = 22050;    % sampling freq
hear_input =    false;
[raw_y, clean_y, sigma2] = ...
    load_signal(selected, offset, duration, snrlvl, fs);

%% set restoration parameters, pad original signal with zeros

window_size = 512;
[y, n_frames, n_samples] = pad_signal(raw_y, window_size);
[padded_clean_y, ~, ~] = pad_signal(clean_y, window_size);


%% relevant matrices
window_func = 1; % shape of the frame
[A, B, Bt, X, invX, invXB, BtinvX, atoms, gabor_mask, t] = ...
    compute_matrices(window_size, window_func, rho);

all_atoms_T = atoms.all_atoms;
all_atoms=atoms.all_atoms';

return
%%
admm_loop_kalman