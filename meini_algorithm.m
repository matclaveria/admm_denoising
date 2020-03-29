function X = meini_algorithm(diag_M, off_diag_M)

    X = diag_M;
    A0 = diag_M;
    B0 = off_diag_M;

    nmax = 2000;
    epsi = 10^-22;
    
    for k = 1 : nmax
        Xprev = X;
        Bprev = B0;
        Aprev = A0;
        % B0 = Bprev*inv(Aprev)*Bprev;
        % A0 = Aprev - Bprev * inv(Aprev)* Bprev' - Bprev' * inv(Aprev)* Bprev;
        % X0 = Xprev - Bprev' * inv(Aprev)* Bprev;
        B0 = Bprev*(Aprev\Bprev);
        A0 = Aprev - Bprev * (Aprev\Bprev') - Bprev' * (Aprev \ Bprev);
        X = Xprev - Bprev' * (Aprev \ Bprev);
    
        if abs(abs(abs(X-Xprev))) < epsi
            break;
        end
        % fprintf('%f\n',norm_max)
        % disp(norm_max)
    end
    
end