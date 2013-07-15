

l = load_pred('/cbio/grlab/nobackup/projects/mip/human_sim/conf_bhatt_mip_quant_exon_eta1_1.20_eta2_0.00_lambda_0_mm0/')

cnt = 1;
p_vals = -ones(1, length(l)*100);
for j = 1:length(l)
	for k = 2:length(l(j).conf_res)
		t = 2*l(j).conf_res(k)-2*l(j).conf_res(1);
	
		if t>100
			p_val = 0;
		else
			p_val = 1 - my_chi2cdf(t, 1);
		end
		p_vals(cnt) = p_val;
		cnt = cnt+1;
	end
end
p_vals(cnt:end) = [];
