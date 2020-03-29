%% set basic ADMM params
clear; close all
rho = 0.2;

%% load data

% selected:
% 1 harpsichord
% 2 glockenspiel
% 3 jazzTrio
selected = 1;
offset   = 0;  % seconds
% duration = 2.614;  % seconds
duration = 1.0;  % seconds
snrlvl   = 11; % SNR of the artificial noise
fs = 22050;    % sampling freq
hear_input = false;
[raw_y, clean_y, sigma2] = ...
    load_signal(selected, offset, duration, snrlvl, fs);

%% set restoration parameters, pad original signal with zeros

window_size = 512;
[y, n_frames, n_samples] = pad_signal(raw_y, window_size);
[padded_clean_y, ~, ~] = pad_signal(clean_y, window_size);

if hear_input
    %%
    sound(y,fs)
    %%
    sound(padded_clean_y,fs)
end

%% compute relevant matrices
window_func = 1; % shape of the frame
[A, B, Bt, X, invX, invXB, BtinvX, atoms, gabor_mask, t] = ...
    compute_matrices(window_size, window_func, rho);

all_atoms_T = atoms.all_atoms;
all_atoms=atoms.all_atoms';

Q = calcQ(X, invXB, BtinvX, n_frames);

M = BtinvX*(Q\(invXB));

return

%%
admm_loop