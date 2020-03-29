function NE = calcQ(X, Q, P, n_blocks)

    % 1
    invX = inv(X);
    NE = zeros(size(X));
    Q_ = eye(size(X));
    P_ = eye(size(X));
    for s = 1 : (n_blocks+1)
        thisTerm = Q_ * (invX * P_);
        NE = NE + thisTerm;
        Q_ = Q_*Q;
        P_ = P_*P;
        if ~mod(s,10) || s == n_blocks
            fprintf('iter:%d; last Q val:%e\n',s, max(max(abs(thisTerm))));
        end
        if max(max(abs(thisTerm))) < 10^-250
           break;
        end
    end

end