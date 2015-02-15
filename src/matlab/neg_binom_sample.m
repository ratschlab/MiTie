function [ret E] = neg_binom_sample(r, p)

assert(p<1 & p>=0)

val = 0;
cnt = 0;
while val<r
	cnt = cnt+1;
	if rand>p
		val = val+1;
	end
end

ret = cnt-val;

% expected value
E = p*r/(1-p);
