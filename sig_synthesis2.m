function y1 = sig_synthesis2(all_atoms, coeffs)
    % 
    n_frames = size(coeffs, 2);
    window_size = size(all_atoms, 1);
    n_y_chunks = (n_frames + 1) / 2 ;
    
    y1 = all_atoms*coeffs(:, 1:2:end);
    y2 = all_atoms*coeffs(:, 2:2:end);
    y2 = reshape(y2, [window_size*(n_y_chunks-1),1]);
    y2 = [y2; zeros(window_size,1)];
    y2 = circshift(y2, window_size/2);
    y2 = reshape(y2, size(y1));
    
    y1 = y1 + y2;    

end