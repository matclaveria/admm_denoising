function [y, n_frames, n_samples] = pad_signal(raw_y, window_size)
    
    sig_len = length(raw_y);

    n_frames = ceil(sig_len/window_size);

    samp_diff = n_frames*window_size - sig_len;

    if samp_diff > 0
        y = [zeros(ceil(samp_diff/2),1);
            raw_y;
            zeros(floor(samp_diff/2),1)];
    end
    n_samples = n_frames*window_size;
    n_frames = 2 * n_frames - 1;
    fprintf([inputname(1), ' padded with %d nought values.\n'],samp_diff)
end