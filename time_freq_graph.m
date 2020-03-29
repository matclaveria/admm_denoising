%coeffs_M = reshape(c_k, [2*n_freqs, n_windows])';
% figure
% coeffs_M = reshape(mean_coeffs2, [2*n_freqs, n_windows])';
close all
coeffs_M = c_k;
n_windows = n_frames;
n_freqs = size(c_k,1)/2+1;
c_k_square = zeros(n_freqs,n_windows);
c_k_square(1,:) = c_k(1,:).^2;
c_k_square(end,:) = c_k(end,:).^2;
c_k_square(2:(end-1),:) = c_k(2:2:(end-1),:).^2 + c_k(3:2:end,:).^2;

    
figure;
imagesc(c_k_square.^0.2)
set(gca,'YDir','normal')
xlabel('Window ID')
%return
% spaced_freqs = [1 100 200 300 400 500 640 800];
% spaced_freqs = round([1 100 200 300 400 500]/500*length(freqs));
%freq_printable = num2str(freqs(spaced_freqs)/2/pi); % ylabels
% yticks(spaced_freqs)
% yticklabels(split(num2str(round(freqs(spaced_freqs)/2/pi))));
%set(gca, 'YTick', 1:1:500, 'YTickLabel', freq_printable);
% ylabel('Frequency [Hz]')
return

T=linspace(1,40,200); 
datefmt = 'HH:MM';  %adjust as appropriate
xticks = 1:20:200;  %adjust as appropriate, positive integers only
xlabels = cellstr( datestr( T(xticks), datefmt ) );  %time labels
% xlabels = num2str( ( T(xticks)) ); % you don't need this 



%%
c_sig = g_t_product(all_atoms', y_folded);
c_k_square = zeros(n_freqs,n_windows);
c_k_square(1,:) = c_sig(1,:).^2;
c_k_square(end,:) = c_sig(end,:).^2;
c_k_square(2:(end-1),:) = c_sig(2:2:(end-1),:).^2 + c_sig(3:2:end,:).^2;
figure;
imagesc(c_k_square.^0.2)
set(gca,'YDir','normal')