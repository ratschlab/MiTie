function ret = gamma_sample(k, theta)

% variance is k*theta^2
% mean is k*theta

ret = 0;
for j = 1:k
	ret = ret - log(rand);
end
ret = ret*theta;
