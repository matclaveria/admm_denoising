function [y, clean_y, sigma2] = ...
    load_signal(selected, offset, duration, snrlvl, fs)
    
    % 1. load signal
    switch selected
        % case 0 
            % [xx,Fs_] = audioread('Piano1.wav');
        case 1
            filename = 'harpsi-cs.wav';
            fprintf('Harpsichord selected.\n')
        case 2
            filename = 'ieeesal_glocken_s.wav';
            fprintf('Glockenspiel selected.\n')
        case 3
            filename = 'JazzTrio.wav';
            fprintf('JazzTrio selected.\n')
        otherwise
            fprintf('Invalid value. No input was selected.\n')
            return
    end
    [input_y, input_fs] = audioread(filename);
    
    % 2. set fs to 22050 and the signal to mono
    
    aux_y = resample(input_y, fs, input_fs);
    
    aux_y = aux_y(:,1);
    
    % 3. set duration
    n_samples = round(duration * fs);
    n_offset  = round(offset * fs);
    if (n_offset + n_samples) > length(aux_y)
        fprintf(['Set valid time intervals for the input signal.'])
        return
    end
    
    y = aux_y(n_offset + 1:(n_samples));
    clean_y = y;
    % 4. add noise
    % typical snrlvl value: 10
    sigma2 = var(y)*10^(-snrlvl/10); % add noise
    y = y + sqrt(sigma2)*randn(size(y));
    
end