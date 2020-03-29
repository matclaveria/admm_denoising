function f_ = fast_matrix_prod(A, B, Bt, x_)
    % T * x_ = f_
    x_folded = reshape(x_, [size(A,1), length(x_)/size(A,1)]);
    i0 = zeros(size(A,1),1);
    f1 = A*x_folded;
    f2 = [B*x_folded(:, 2:end), i0];
    f3 = [i0, Bt*x_folded(:, 1:(end-1))];
    f_ = reshape(f1+f2+f3, [size(A,1) * length(x_)/size(A,1),1]);
end