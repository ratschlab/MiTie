function intron_conf = compute_intron_list(admat)

	s = size(admat, 1);
	intron_conf = zeros(s*(s-1)/2, 3);
	cnt = 0;
	for j = 1:s
		for k = j+1:s
			if admat(j, k)>=0
				% valid intron not just a neighboring segment
				cnt = cnt+1;
				intron_conf(cnt,:) = [j, k, admat(j, k)];
			end	
		end
	end
	intron_conf(cnt+1:end, : ) = [];
return
