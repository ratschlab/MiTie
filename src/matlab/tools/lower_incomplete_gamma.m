function ret = lower_incomplete_gamma(s, x)


% sum_=0;
% using log(a+b) = log(a) + log(1+exp(log(b)-log(a)))
k=0;
sum_ = log(x)*k-gammaln(s+k+1);
for k = 1:100
	val = log(x)*k-gammaln(s+k+1);
	sum_ = sum_+log(1+exp(val-sum_));
	%sum_ = sum_+exp(log(x)*k-gammaln(s+k+1));
end

%ret = exp(s*log(x)+gammaln(s)-x)*sum_;
ret = exp(s*log(x)+gammaln(s)-x+sum_);
