function [T, N] = optional_matrices(A, B, Bt, X, n, n_blocks)

    T = zeros(size(A)*n_blocks);

    ind_k = 1 : n;

    for k = 1 : n_blocks

        T(ind_k,ind_k) = A;

        if k > 1
            T(ind_k,ind_k-n) = Bt;
        end

        if k < n_blocks
            T(ind_k,ind_k+n) = B;
        end

        ind_k = ind_k + n;

    end

    ind_k = 1 : n;
    N(ind_k,ind_k) = X;
end