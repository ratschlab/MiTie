fprintf('usage of segments in transcripts\n');
fprintf('initial tr:\t');
for k = 1:s
	fprintf('%3i', sum(admat(k, 1:k-1)>=-1)==0)
end
fprintf('\n')
fprintf('terminal tr:\t');
for k = 1:s
	fprintf('%3i', sum(admat(k, k+1:end)>=-1)==0)
end
fprintf('\n')

for j = 1:t
	fprintf('tr%2i(w_t=%.2f):\t', j, abs(res_W(j)));
	for k = 1:s
		fprintf('%3.0f', res_U((k-1)*t+j))
	end
	fprintf('\n')
end
fprintf('expected cov:\n');
for k = 1:t
	fprintf('  trans: %i:\t', k)
	for j = 1:s
		fprintf('%3.0f', result(E_idx((i-1)*s*t+(j-1)*t+k))*100);
	end
	fprintf('\n')
end
fprintf('segment cov:\t');
for k = 1:s
	fprintf('%3.0f', observed_cov(i,k)*100);
end
fprintf('\n')
fprintf('coverage dev:\t');
for j = 1:s
	fprintf('%3.0f', result(L_idx((i-1)*s+j))*100);
end
fprintf('\n')
fprintf('\n')

fprintf('Usage of introns (TODO: fix for multiple samples!!!): \n')
expected_cov=zeros(1,size(intron_conf,1)) ;
for x = 1:t
	fprintf('tr%2i(w_t=%.2f):\t', x, abs(res_W(x)));
	for j = 1:size(intron_conf, 1)
		fprintf('%3.0f', res_C((x-1)*c+j)*100)
        expected_cov(j) = expected_cov(j) + res_C((x-1)*c+j) ;
	end	
	fprintf('\n')
end
fprintf('exp intron cov:\t')
for j = 1:size(intron_conf, 1)
	fprintf('%3.0f', expected_cov(j)*100)
end	
fprintf('\n')
fprintf('obs intron cov:\t')
for j = 1:size(intron_conf, 1)
	fprintf('%3.0f', intron_conf(j,3)*100)
end	
fprintf('\n')
fprintf('Deviation:\t')
for j = 1:size(intron_conf, 1)
	fprintf('%3.f', res_D(j)*100)
end	
fprintf('\n')

if 0%use_pair
	fprintf('\n')
	for x = 1:t
		fprintf('PU: dist<is+tol (t=%i):\t', x)
		for j = 1:size(pair_observation{i}, 1)
			fprintf('%.2f ', abs(res_PU((j-1)*t+x)))
		end
		fprintf('\n')
	end
	for x = 1:t
		fprintf('PL: dist>is-tol (t=%i):\t', x)
		for j = 1:size(pair_observation{i}, 1)
			fprintf('%.2f ', abs(res_PL((j-1)*t+x)))
		end
		fprintf('\n')
	end
	for x=1:t
		fprintf('pair dist (t=%i):\t', x)
		for y=1:size(pair_observation{i}, 1)
			j = pair_observation{i}(y, 1);
			k = pair_observation{i}(y, 2);
			idx = (y-1)*t+x;
			tlen = 0.5*len(j)*res_U((j-1)*t+x);
			tlen = tlen+0.5*len(k)*res_U((k-1)*t+x);
			for l=j+1:k-1
				tlen = tlen+0.5*len(l)*res_U((l-1)*t+x);
			end
			fprintf('%4.0f ', tlen)
		end
		fprintf('\n')
	end
	for x = 1:t
		fprintf('usage of j and k (t=%i):\t', x)
		for y = 1:size(pair_observation{i}, 1)
			j = pair_observation{i}(y, 1);
			k = pair_observation{i}(y, 2);
			fprintf('%2i%2i ', round(res_U((j-1)*t+x)), round(res_U((k-1)*t+x)))
		end
		fprintf('\n')
	end

	fprintf('\n')
	for x = 1:t
		fprintf('expected pair cov (t=%i):\t', x)
		for y = 1:size(pair_observation{i}, 1)
			fprintf('%.2f ', abs(res_X((y-1)*t+x)))
		end
		fprintf('\n')
	end
	fprintf('\n')
	fprintf('observed pair coverage:\t\t')
	for j = 1:size(pair_observation{i}, 1)
		fprintf('%.2f ', pair_observation{i}(j,3))	
	end
	fprintf('\n')
	fprintf('deviation of pair cov:\t\t')
	for j = 1:size(pair_observation{i}, 1)
		fprintf('%.2f ', abs(res_DP(j)))	
	end
	fprintf('\n')
	fprintf('\n')
end

fprintf('Objective (sample%i): %.3f\n', i, result'*Q*result+f'*result);
idx = L_idx((i-1)*s+1:(i-1)*s+s);
fprintf('\t-> exon coverage: %.3f\n', res_L*Q(idx,idx)*res_L'+f(idx)'*res_L');
fprintf('\t-> transcript sum: %.3f\n', res_I*Q(I_idx,I_idx)*res_I'+f(I_idx)'*res_I');
idx = D_idx((i-1)*c+1:(i-1)*c+c);
fprintf('\t-> intron coverage: %.3f\n', res_D*Q(idx,idx)*res_D'+f(idx)'*res_D');
if 0%use_pair
	idx = DP_idx((i-1)*p+1:(i-1)*p+p);
	fprintf('\t-> pairs: %.3f\n', res_DP*Q(DP_idx,DP_idx)*res_DP'+f(DP_idx)'*res_DP');
end
if use_pair
	fprintf('\t-> pairs: %.3f\n', f(SP_idx)'*result(SP_idx));
end

