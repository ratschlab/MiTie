

%load /cbio/grlab/nobackup2/projects/mip/human_annotation/hg19_annotations_merged_splice_graph_chr1.mat
load /cbio/grlab/nobackup2/projects/mip/human_annotation/hg19_annotations_merged_splice_graph.mat

if 0
	fn_bam = '/cbio/grlab/nobackup2/sequencing_runs/human/ENCODE/wgEncodeCshlLongRnaSeqK562CellPapAlnRep1.bam';
elseif 0
	%unix('ln -s /cbio/grlab/nobackup3/TCGA/PanCancer/bam_valid/BLCA/TCGA-DK-A3IM-01A-11R-A20F-07.d0cbd85a-e244-4f99-9adf-71f7afa5f6ec.v3.star.bam ~/tmp/')
	%unix('samtools sort ~/tmp/TCGA-DK-A3IM-01A-11R-A20F-07.d0cbd85a-e244-4f99-9adf-71f7afa5f6ec.v3.star.bam ~/tmp/TCGA-DK-A3IM-01A-11R-A20F-07.d0cbd85a-e244-4f99-9adf-71f7afa5f6ec.v3.star.sorted')
	%unix('samtools index ~/tmp/TCGA-DK-A3IM-01A-11R-A20F-07.d0cbd85a-e244-4f99-9adf-71f7afa5f6ec.v3.star.sorted.bam')
	fn_bam = '/cbio/grlab/home/jonas/tmp/TCGA-DK-A3IM-01A-11R-A20F-07.d0cbd85a-e244-4f99-9adf-71f7afa5f6ec.v3.star.sorted.bam'
else
	fn_bam = '/cbio/grlab/nobackup2/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired.bam'
end

for j = 1:length(genes), len(j) = length(genes(j).exons); end

genes = genes(len==1);

min_len = 30;

fn_save = '~/tmp/mean_var_data_TCGA.mat';
fn_save = '~/tmp/mean_var_data_sim.mat';
if fexist(fn_save)
	load(fn_save, 'mv', 'num')
else
	mv = zeros(length(genes), 2);
	num = 0;
end

readlen = 70; 

while num<length(genes)
	num = num+1;
	if mod(num, 100)==0
		fprintf('%i (%i) genes done\n', num, length(genes));
		save(fn_save, 'mv', 'num')
	end
	exons = genes(num).exons{1};
	num_exons = size(exons, 1);
	if num_exons==1
		%continue
	end
	c = 0;
	regs = zeros(0, 2);
	while c<20
		for k = 1:size(exons, 1)
			len = exons(k, 2)-exons(k, 1);
			if len<min_len+readlen
				continue
			end
			s = exons(k, 1) + floor(rand*(len-min_len-readlen));
			reg = [s s+min_len];

			idx = find(regs(:, 1)-readlen<=reg(2) & regs(:, 2)+readlen<=reg(1));
			if isempty(idx)
				regs = [regs; reg];
			end
		end
		c = c+1;
	end
	if size(regs, 1)<2
		continue;
	end

	counts = zeros(1, size(regs, 1));
	for j = 1:size(regs, 1)
		%reg1_str = sprintf('chr%s:%i-%i', genes(num).chr, regs(j,1), regs(j,2));
		reg1_str = sprintf('%s:%i-%i', genes(num).chr, regs(j,1), regs(j,2));
		[tmp ret1]= unix(sprintf('samtools view %s %s | cut -f 4', fn_bam, reg1_str));
			
		starts1 = str2num(ret1);
		counts(j) = length(find(starts1>=regs(j,1) & starts1<regs(j,2)));
	end
	fprintf('m: %.2f, var:%.2f num:%i, len:%i\n', mean(counts), var(counts), num, length(counts));
	mv(num, 1) = mean(counts);
	mv(num, 2) = var(counts);
end
mv(num+1:end, :) = [];

save(fn_save, 'mv', 'num')

figure
load(fn_save, 'mv')
m = mv(:, 1);
v = mv(:, 2);
idx = find(m>0);
m = m(idx);
v = v(idx);
v = (sqrt(v)-m).^2; % 

pct = prctile(m, 0:1:100);

nmean = 0;
nstd = 5000;
bestfit = 1e20;
for eta1 = 0%:0.1:2
	for eta2 = 0:0.01:3
		y = (1+eta1)*m + eta2*m.^2;
		fit = sum((y-v).^2.*normpdf(m, nmean, nstd));
		if fit<bestfit;
			bestfit=fit
			best_eta1 = eta1
			best_eta2 = eta2
		end
	end
end

loglog(m,v,'.',m,(1+best_eta1)*m+best_eta2*m.^2,'r.') 
%print -depsc ~/tmp/plot_mean_var_mod.eps
print -depsc ~/tmp/plot_mean_var_sim.eps

