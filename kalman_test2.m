%%
admm_audio
tests = 0;
lazy_admm_loop;

%%

z_k = (v_k-eta_k/rho);  
[K_kf, G_rts, yf,  G1, Gp, GP] = ...
    kf_rts_setup(y_folded, z_k, all_atoms, rho);

c_kalman = kf_rts(y_folded, z_k, K_kf, G_rts, G1, Gp, GP);



