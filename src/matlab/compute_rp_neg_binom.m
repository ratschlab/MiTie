function [r, p] = compute_rp_neg_binom(mu,var)

p = 1-(mu/var);
%p = -(mu-2*var)/(2*var) - sqrt(((mu-2*var)/(2*var))^2 +mu/var -1);
r = mu*(1-p)/p;
return
