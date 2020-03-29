KKK = 20;
T_all = zeros(KKK,1);
snr_out_all = zeros(KKK,1);
for kkk = 1 : KKK 
    admm_audio_kalman 
    admm_loop_kalman2
    T_all(kkk) = T_final;
    snr_out_all(kkk) = snr_out;
end

%%
fprintf('Average time [s]: %0.2f\n',mean(T_all)*60)
fprintf('Average SNR (out): %0.2f\n',mean(snr_out_all))