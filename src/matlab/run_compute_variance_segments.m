fn_save = '~/tmp/mean_var_segments.mat';
fn_save = '~/tmp/mean_var_introns.mat';

fn_gtf = '/cbio/grlab/nobackup2/projects/mip/human_annotation/Homo_sapiens.GRCh37.68.chr.gtf';

out_dir='/cbio/grlab/nobackup/projects/mip/human_sim/anno_graph/';
fn_graph=sprintf('%sHomo_sapiens.GRCh37.68.chr.graphs.h5', out_dir);

nice_mkdir(out_dir);


if ~fexist(fn_graph)
	unix(sprintf('./generate_segment_graph %s --gtf %s', fn_graph, fn_gtf));
end

fn_bam = '/cbio/grlab/nobackup2/sequencing_runs/human/ENCODE/wgEncodeCshlLongRnaSeqK562CellPapAlnRep1.bam';
%fn_bam = '/cbio/grlab/home/jonas/tmp/TCGA-DK-A3IM-01A-11R-A20F-07.d0cbd85a-e244-4f99-9adf-71f7afa5f6ec.v3.star.sorted.bam';
%fn_bam = '/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired.bam';

if fexist(fn_save)
	load(fn_save, 'mean_', 'var_', 'num')
	cnt = find(mean_==-1, 1)
else
	cnt = 1;
	num = 1;
	mean_ = -ones(1, 1e5);
	var_ = -ones(1, 1e5);
end

for idx = num:100:40000

	num = idx;
	save(fn_save, 'mean_', 'var_', 'num')
	genes = load_graph_bin(fn_graph, fn_bam, idx:idx+99);

	for k = 1:length(genes)
		if size(genes(k).transcripts, 1)==1
			%cov_scale = [];
			%for r = 1:size(genes(k).seg_admat, 3)%loop over samples
			%	cov_scale(r) = max(genes(k).coverage(r, :));
			%	cov_scale(r) = max([cov_scale(r), max(genes(k).seg_admat(:, :, r))]);
			%	if cov_scale(r)>0
			%		sa = genes(k).seg_admat(:,:,r);
			%		sa(sa>0) = sa(sa>0)/cov_scale(r);
			%		genes(k).seg_admat(:,:,r) = sa;
			%		genes(k).coverage(r, :) = genes(k).coverage(r, :)/cov_scale(r);
			%	end
			%end
			len = genes(k).segments(2, :)-genes(k).segments(1, :); 

			readlen = 100;
			coverage = genes(k).coverage/readlen;% readlen

			intron_conf = compute_intron_list(genes(k).seg_admat);

			mean_(cnt) = mean(intron_conf(:, 3));
			var_(cnt) = var(intron_conf(:, 3));
			%mean_(cnt) = mean(coverage);
			%var_(cnt) = var(coverage);
			cnt = cnt+1;
		end
	end
end
mean_(cnt:end) = [];
var_(cnt:end) = [];
save(fn_save, 'mean_', 'var_', 'num')

return

fn_save = '~/tmp/mean_var_segments.mat';
%fn_save = '~/tmp/mean_var_introns.mat';
load(fn_save, 'mean_', 'var_')

m = mean_;
v = var_;
idx = find(m>0 & ~isnan(var_));
m = m(idx);
v = v(idx);
v = (sqrt(v)-m).^2; % 
%v = sqrt(v);

bins = my_prctile(m, 0:5:100);
vmedian = [];
for j = 1:length(bins)-1
	idx = find(m>=bins(j)&m<=bins(j+1)); 
	vmedian(j) = median(v(idx));
	mmedian(j) = median(m(idx));
end

%bins = (bins(1:end-1)+bins(2:end))/2;
bins = mmedian;

nmean = 0;
nstd = 2500;
bestfit = inf;
for eta1 = 0:0.01:3
	for eta2 = 0:0.01:3
		if 0
			y = (1+eta1)*m + eta2*m.^2;
			%fit = sum((y-v).^2.*normpdf(m, nmean, nstd));
			fit = sum((y-v).^2.*(1/(2*pi*nstd^2)*exp(-(m-nmean).^2./(2*nstd^2))));
		else
			y = (1+eta1)*bins + eta2*bins.^2;
			fit = sum((y-vmedian).^2.*(1/(2*pi*nstd^2)*exp(-(bins-nmean).^2./(2*nstd^2))));
		end
		if fit<bestfit;
			bestfit=fit;
			best_eta1 = eta1;
			best_eta2 = eta2;
		end
	end
end
best_eta1, best_eta2



close all; loglog(m,v,'.',m,(1+best_eta1)*m+best_eta2*m.^2,'r.', bins,vmedian,'g.') 
%print -depsc ~/tmp/plot_mean_var_mod.eps
print -depsc ~/tmp/plot_mean_var_segments.eps
%print -depsc ~/tmp/plot_mean_var_introns.eps

%close all; loglog(m,v./m.^2,'.',m,((1+best_eta1)*m+best_eta2*m.^2)./m.^2,'r.') 
%print -depsc ~/tmp/plot_mean_var_mod.eps
%print -depsc ~/tmp/plot_mean_scv_introns.eps

