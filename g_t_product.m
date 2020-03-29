function c_replica2 = g_t_product(g_t, z1T_folded)
    % 
    window_size = size(g_t, 1);
    n_chunks = size(z1T_folded, 2);
    n_frames = 2*n_chunks-1;
    c_replica2 = zeros(window_size, n_frames);
    
    z1T_folded_shifted = reshape(z1T_folded, [window_size*n_chunks, 1]);
    z1T_folded_shifted = reshape(...
        z1T_folded_shifted(window_size/2 + ...
        (1:((n_chunks-1)*window_size))), ...
        [window_size, n_chunks-1]);
    
    c_replica2(:, 1:2:end) = g_t*z1T_folded;
    c_replica2(:, 2:2:end) = g_t*z1T_folded_shifted;
end


% nn = 6;
% g_t = magic(nn);
% n_chunks = 10;
% n_frames = 2*n_chunks-1;
% 
% G = zeros(nn*n_frames, nn*n_chunks);
% ind1 = 1:nn;
% ind2 = 1:nn;
% for i = 1 : n_frames
%     G(ind1,ind2) = g_t;
%     ind1 = ind1 + nn; 
%     ind2 = ind2 + nn/2;
% end
% 
% zet = 1:(nn*n_chunks);
% 
% dummyvec = sin(zet*0.2)';
% dummyvec_rs = reshape(dummyvec, [nn, n_chunks]);
% 
% 
% c = G*dummyvec;
% plot(c)
% c_replica2 = g_t_product2(g_t, dummyvec_rs);
% figure
% c_replica2_rs = reshape(c_replica2, [nn*n_frames,1]);
% plot(c_replica2_rs)