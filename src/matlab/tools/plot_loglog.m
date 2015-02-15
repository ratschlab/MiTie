


plot(log(true_vals), log(expr_vals), '.'); hold on; 
plot(-ones(1, sum(true_vals==0)), log(expr_vals(true_vals==0)), '.'); 
plot(log(true_vals(expr_vals==0), -5*ones(1, sum(expr_vals==0)), '.')
