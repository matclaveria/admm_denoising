clear
window_size = 512;
nframes = 70;
c = rand(window_size,nframes);
eta = 0.4*randn(window_size,nframes);
lambda = 4;
rho = 4;
ind_r = 2:2:(window_size-1);
ind_i = 3:2:(window_size);

v = shrinking_op(c, eta, lambda, rho, ind_r, ind_i);

%%
clc
m = 256;
n = nframes;
if m == 1
    indc = 1;
elseif m == 257
    indc = 512;
else
    indc = [2*m-2,2*m-1];
end
val1 = shrink_just_one(c, eta, lambda, rho, ind_r, ind_i, m, n)
v(indc, n)

%%

function val = shrink_just_one(c, eta, lambda, rho, ind_r, ind_i, m, n)
    % extract c
    if m == 1
        indc = 1;
    elseif m == 257
        indc = 512;
    else
        indc = [2*m-2,2*m-1];
    end
    r = c(indc, n) + eta(indc, n)/rho; 
    norm_r = sqrt(sum(r.^2));
    
    if (norm_r - lambda/rho) > 0
        val = (norm_r - lambda/rho) .* r./norm_r; 
    else 
        val = indc'*0;
    end
end
