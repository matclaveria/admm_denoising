%% ADMM
last_sample = 0;
use_threshold = 1;
max_iter = 1500;
lambda0 = 1.2;

freqs = 2*pi * (0: (window_size/2))/window_size;
freqs = freqs';
alpha_f = 2 * pi * 0.5; %cut-off freq: 10,000 Hz
f_modulation = @(f) 1 ./ sqrt(1 + alpha_f^-1.5 * f.^2);
f_val = f_modulation(freqs);
lambda = lambda0*repmat(f_val, [1 n_frames]);

tic
if ~last_sample
c_k    = zeros(window_size, n_frames);
v_k    = zeros(window_size, n_frames);
eta_k  = zeros(window_size, n_frames);
end
y_folded = reshape(y, [window_size, length(y)/window_size]);
ind_r = 2:2:(window_size-1);
ind_i = 3:2:(window_size);

GtY = g_t_product(all_atoms_T, y_folded);

fprintf('-----------\n')
fprintf('Running SGS\n')
fprintf('\n          ')
for k = 1 : max_iter
    
    %fprintf('\b\b\b\b\b\b\b\b\b\b %09d', k);
    if mod(k,20) == 0 || k==1
%         fprintf('%09d, %d\n', k, ...
%             sum(sum(eta_k.*(c_k-v_k))));
        fprintf('%09d, %d\n', k, ...
            max(max(abs(c_k - v_k))));
    end
    
    % step 1: coefficient vector c:
%     mu_tilde = GtY + rho*(v_k-eta_k); 
    % eta/rho
    mu_tilde = GtY + rho*(v_k-eta_k/rho); % eta/rho
    if k == 100
        break
    end
    c_k = malyshev_alg(B, BtinvX, invX, M, mu_tilde, n_frames);
    
    % step 2: replica vector v:
    v_k = shrinking_op(c_k, eta_k, lambda, rho, ind_r, ind_i);
    
    % step 3: dual variable eta:
%     eta_k = eta_k + (c_k - v_k);
    eta_k = eta_k + rho*(c_k - v_k);
    if max(max(abs(c_k - v_k))) < 0.5*10^-7 && use_threshold
        break
    end
end

fprintf('\n\n')
toc
y_synth = reshape(sig_synthesis2(all_atoms, c_k), size(y));

%%
sound(y_synth,fs)
close all

plot(y); hold on; plot(y_synth)
return
%%
sound(y, fs)

%%
y_synth = reshape(sig_synthesis2(all_atoms, c_k), size(y));
plot(y); hold on; plot(y_synth)


%%
func1 = 0.5*norm(y_folded-sig_synthesis2(all_atoms, ck),'fro').^2 + rho/2*norm(ck(:)./(sqrt(V(:))+1e-6))^2;
func2 = 0.5*norm(y_folded-sig_synthesis2(all_atoms, z),'fro').^2 + rho/2*norm(ck(:)./(sqrt(V(:))+1e-6))^2;