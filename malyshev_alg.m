function x = malyshev_alg(B, BtinvX, invX, M, f_reshaped, n_blocks)

    y = bt_solver1(BtinvX, invX, B, f_reshaped, n_blocks);
    % y_ast_1 = M*y(:,1);
    %z = bt_solver2(BtinvX, invX, B, y_ast_1, n_blocks);
    f_prime = f_reshaped*0;
    f_prime(:,1) = M*y(:,1);
    z = bt_solver1(BtinvX, invX, B, f_prime, n_blocks);
    x = y - z;
    
end

function z = bt_solver1(BtinvX, invX, B, f, n_frames)
    
    z = 0*f;
    
    % part 1: solve system L * w = f
    z(:, 1) = f(:, 1);

    for k = 2 : n_frames  
        z(:, k) = f(:, k) - BtinvX * z(:, k-1);
    end

    % part 2: solve system U * z = w
    % z is compute on matrix w for efficiency    
    z(:, end) = invX*z(:, end);

    for k = (n_frames-1):-1: 1 
        z(:, k) = invX*(z(:, k) - B*z(:, k+1));
    end
    
end

function z = bt_solver2(BtinvX, invX, B, y_ast_1, n_frames)
    
    z = zeros(size(B,1), n_frames);
    
    % part 1: solve system L * w = f
    z(:, 1) = y_ast_1;
    for k = 2 : n_frames  
        z(:, k) = - BtinvX * z(:, k-1);
    end

    % part 1: solve system U * z = w
    % z is compute on matrix w for efficiency
    z(:, end) = invX*z(:, end);
    for k = (n_frames-1):-1: 1 
        z(:, k) = invX*(z(:, k) - B*z(:, k+1));
    end
    
end