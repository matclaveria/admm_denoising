% close all
coeffs_M = c_k;
n_windows = n_frames;
n_freqs = size(c_k,1)/2+1;
c_k_square = zeros(n_freqs,n_windows);
c_k_square(1,:) = c_k(1,:).^2;
c_k_square(end,:) = c_k(end,:).^2;
c_k_square(2:(end-1),:) = c_k(2:2:(end-1),:).^2 + c_k(3:2:end,:).^2;

    
figure;
subplot(1,2,1)
imagesc(c_k_square.^0.2)
set(gca,'YDir','normal')

c_sig = g_t_product(all_atoms', y_folded);
c_k_square = zeros(n_freqs,n_windows);
c_k_square(1,:) = c_sig(1,:).^2;
c_k_square(end,:) = c_sig(end,:).^2;
c_k_square(2:(end-1),:) = c_sig(2:2:(end-1),:).^2 + c_sig(3:2:end,:).^2;
subplot(1,2,2)
imagesc(c_k_square.^0.2)
set(gca,'YDir','normal')