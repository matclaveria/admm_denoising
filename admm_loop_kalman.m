%% ADMM preparation
last_sample = 0;
use_threshold = 1;
max_iter = 500;
lambda0 = 3.5;

freqs = 2*pi * (0: (window_size/2))/window_size;
freqs = freqs';
alpha_f = 2 * pi * 0.5; %cut-off freq: 10,000 Hz
f_modulation = @(f) 1 ./ sqrt(1 + alpha_f^-1.5 * f.^2);
f_val = f_modulation(freqs);
lambda = lambda0*repmat(f_val, [1 n_frames]);

if ~last_sample
c_k    = zeros(window_size, n_frames);
v_k    = zeros(window_size, n_frames);
eta_k  = zeros(window_size, n_frames);
end
y_folded = reshape(y, [window_size, length(y)/window_size]);
ind_r = 2:2:(window_size-1);
ind_i = 3:2:(window_size);

GtY = g_t_product(all_atoms_T, y_folded);

if ~exist('K_kf', 'var') && ~exist('Kstruct', 'var')
    fprintf('Calculating matrices...\n')
    [K_kf, G_rts, ~,  G1, Gp, GP] = ...
        kf_rts_setup(y_folded, randn(size(c_k)), all_atoms, rho);
end

%%
if ~exist('Kstruct', 'var')
    %%
    fprintf('Building K and G structs...\n')
    for p = 1 : size(K_kf,3)
        eval(strcat(['Kstruct.K',num2str(p) ,'= K_kf(:,:,p);']))
        eval(strcat(['Gstruct.G',num2str(p) ,'= G_rts(:,:,p);']))
    end
    clear K_kf;
    clear G_rts;
end
% return
%%
tic
fprintf('-----------\n')
fprintf('Running SGS\n')
fprintf('\n          ')
for k = 1 : max_iter
    
    %fprintf('\b\b\b\b\b\b\b\b\b\b %09d', k);
    if mod(k,10) == 0 || k==1
%         fprintf('%09d, %d\n', k, ...
%             sum(sum(eta_k.*(c_k-v_k))));
        fprintf('%09d, %d\n', k, ...
            max(max(abs(c_k - v_k))));
    end
    
    % step 1: coefficient vector c:
    z_k = (v_k-eta_k/rho);
    % c_k = kf_rts(y_folded, z_k, K_kf, G_rts, G1, Gp, GP);
    c_k = kf_rts_struct(y_folded, z_k, Kstruct, Gstruct, G1, Gp, GP);
        
    % step 2: replica vector v:
    v_k = shrinking_op(c_k, eta_k, lambda, rho, ind_r, ind_i);
    
    % step 3: dual variable eta:
%     eta_k = eta_k + (c_k - v_k);
    eta_k = eta_k + rho*(c_k - v_k);
    if max(max(abs(c_k - v_k))) < 2*10^-8 && use_threshold
        break
    end
end

fprintf('\n\n')

T_final = toc/60;
y_synth = reshape(sig_synthesis2(all_atoms, c_k), size(y));

%%
sound(y_synth,fs)
close all
timesteps = 1/fs * (1: (window_size*(n_frames+1)/2));
plot(timesteps,y); hold on; plot(timesteps,y_synth)
xlabel('Time [s]')
ylabel('y[t]')
axis([-0.05, timesteps(end)+0.05, -0.41 0.45])
legend('Input noisy signal', 'Reconstructed signal', ...
    'Location', 'Northwest')        
return
%%
sound(y, fs)

%%
y_synth = reshape(sig_synthesis2(all_atoms, c_k), size(y));
plot(y); hold on; plot(y_synth)

%% SNR
indie = 6000:length(y)-1000;
snr_in = 10*log10(sum(padded_clean_y(indie).^2)/sum((y(indie) - padded_clean_y(indie)).^2))
snr_out = 10*log10(sum(padded_clean_y(indie).^2)/sum((y_synth(indie) - padded_clean_y(indie)).^2))
%%
func1 = 0.5*norm(y_folded-sig_synthesis2(all_atoms, ck),'fro').^2 + rho/2*norm(ck(:)./(sqrt(V(:))+1e-6))^2;
func2 = 0.5*norm(y_folded-sig_synthesis2(all_atoms, z),'fro').^2 + rho/2*norm(ck(:)./(sqrt(V(:))+1e-6))^2;

%% make pdf
set(gcf, 'Color', 'w');
cd('export_fig')
str_ = strcat(['export_fig ../audio_denoising' ...
    '.pdf', ' -opengl']);
eval(str_)
cd('..')