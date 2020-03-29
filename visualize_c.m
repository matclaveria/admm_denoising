norm_c = sqrt([c_k(1,:).^2;
    c_k(ind_r,:).^2 + c_k(ind_i,:).^2;
    c_k(end,:).^2]);

close
imagesc(norm_c)
set(gca,'YDir','normal')