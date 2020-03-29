%%
% admm_audio
%% test 1: malyshev matrix inversion
% test1_steps = 1;
% 
% tic
% y_folded = reshape(y, [window_size, length(y)/window_size]);
% GtY = g_t_product(all_atoms_T, y_folded);
% for i = 1 : test1_steps
%     incog = GtY;
%     ck = malyshev_alg2(B, BtinvX, invX, M, ...
%         incog, n_frames); 
% end
% toc
% 
% return
%%
sol = ck;
fista_max = 1000;
lipschitz = 445;
alpha = randn(size(ck))*0.05 + ck*0;
z = alpha;
% V = exprnd(0.5, size(ck));
V = rho + alpha*0;

% close all
% plot(ck(:),'-')
% hold on
% plot(z(:),'.')

% %%

tic
for it = 1 : fista_max
    alphaold = alpha;
    % GtY = g_t_product(all_atoms_T, y_folded);
    alpha = z + g_t_product(all_atoms_T, ...
        y_folded - sig_synthesis2(all_atoms, alpha))/lipschitz;
    alpha = (V)./(V+rho/lipschitz) .* alpha;
    z = alpha + (it-1)/(it+3) * (alpha-alphaold);
    
    % assessing convergence criteria?
    % grad = -op.analysis(y - op.synthesis(alpha)) + lambda * alpha./(V+1e-6);
    % norm_grad =  norm(grad,'fro').^2/(M*N);
    % func = 0.5*norm(y-op.synthesis(alpha)).^2 + lambda/2*norm(alpha(:)./(sqrt(V(:))+1e-6))^2;
    err_rel = (norm(alpha-alphaold,'fro')/norm(alphaold,'fro'))^2;
    % fprintf('     in fista, it = %d -- func = %f -- ||Grad|| = %f -- err_rel = %f\n',it_ista,func,norm_grad,err_rel3);
    
    if err_rel < 0.5e-6
        it
        break;
    end
end
toc

close all
plot(z(:),'-')
hold on
plot(ck(:),'.')
figure
plot(z(:)-ck(:),'.')

% for it = 1 : max_it
%     alphaold = alpha;
%     alpha = z +   op.analysis(y - op.synthesis(z))/options.Lipschitz;
%     alpha = (V)./(V+lambda/options.Lipschitz) .* alpha;
%     z = alpha + (it_ista-1)/(it_ista+3) * (alpha-alphaold);
%     grad = -op.analysis(y - op.synthesis(alpha)) + lambda * alpha./(V+1e-6);
%     norm_grad =  norm(grad,'fro').^2/(M*N);
%     func = 0.5*norm(y-op.synthesis(alpha)).^2 + lambda/2*norm(alpha(:)./(sqrt(V(:))+1e-6))^2;
%     err_rel = (norm(alpha-alphaold,'fro')/norm(alphaold,'fro'))^2;
%     fprintf('     in fista, it = %d -- func = %f -- ||Grad|| = %f -- err_rel = %f\n',it_ista,func,norm_grad,err_rel);
%     if norm_grad < 5e-5 || err_rel < 5e-5
%         break;
%     end
% end


%%
y1 = sig_synthesis2(all_atoms,ck);
y2 = sig_synthesis2(all_atoms, z);
close all
plot(y1(:),'.'); hold on
plot(y2(:),'-')

%%
sound(y1(:),fs)

%%
sound(y2(:),fs)

%%
close all
plot(y1(:)-y2(:))