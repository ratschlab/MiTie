function transcript_predictions_MAGIC(chr, strand)

%fn_graph = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/seg_graph20'
%fn_graph = sprintf('/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/seg_graph30_filtered.%s.bin', chr)
fn_graph = sprintf('/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/MiTie/graphs/seg_graph30_filtered_%s%s', chr, strand);
out_dir = sprintf('/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/MiTie/pred_%s%s', chr, strand); 

nice_mkdir(out_dir);

%% parse clusters from file
%fn_cluster = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/classify_bam_out.list';
%fn_mapping = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/library_ids'
%
%mapping = importdata(fn_mapping);
%WTCHG_hash = {};
%for j = 1:length(mapping)
%	c = separate(mapping{j});
%	num = str2num(c{2}(13:end));
%	%assert(length(WTCHG_hash)<num || isempty(WTCHG_hash{num}));
%	if length(WTCHG_hash)>=num && ~isempty(WTCHG_hash{num})
%		WTCHG_hash{num} = [WTCHG_hash{num} ',' c{1}];
%	else
%		WTCHG_hash{num} = c{1};
%	end
%end
%
%clust_str = importdata(fn_cluster);
%for j = length(clust_str):-1:1
%	c = separate(clust_str{j});
%	gene_names{j} = c{1};
%	for k = 2:length(c)
%		clusters{j, k-1} = c{k};
%	end
%end
%
%% load genome annotation because the graphs have no ATG names associated
%fn_annotation = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/consolidated_annotation.allmerged.Col_0.mat'; 
%l = load(fn_annotation, 'genes');

% change this flag to run transcript predictions 
% for each gene distributed on a cluster
parallel = 1

%fn_bam_tmp = dir('/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/alignments_magic/bam_mmr/*.bam');
%for j=1:length(fn_bam_tmp)
%	fn_bam{j} = ['/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/alignments_magic/bam_mmr/', fn_bam_tmp(j).name]; 
%end

% model selection parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
	C.num_transcripts = 500;
	C.introns = 100; 
	C.exon_cov = 100; 
end

C.num_transcripts_predef = 1.6;
C.pairs = 100; 
C.num_clusters=10;

param.insert_size_tol = 0.5; 

param.loss = 'nb';
param.num_clusters=1;

param.slope_threash = 3.1;

% nois level: this is the mean and variance of the poisson accounting for the noise
param.lambda = 1; 

% no model selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.number_of_additional_transcripts=1;
param.use_predef_trans=1;
param.len_cutoff=120*2;
param.time_limit=3600*5;

rand('seed', ceil(cputime*10));
jobid = ceil(rand*100);
%genes = genes([genes.id] == 29899);
jobinfo = rproc_empty(0);
cnt = 0;
gene_cnt = 0;

skip = -1;
for j = 0:100000
	fn_res = sprintf('%s/gene%i.mat', out_dir, j);
	if fexist(fn_res)
		cnt = 0;
		skip=0;
	else
		cnt = cnt+1;
		if skip==-1
			skip = j;
		end
	end
	if cnt>100
		break;
	end
end
	
% close file descriptor, because this is a static variable
% might be open from a previous call to this function
tmp = load_graph_bin_MAGIC(fn_graph, -1, 0);

cnt = 0;
% open binary file
while 1

	if skip>0
		gene_cnt = skip;
	end
	num = 100;
	genes = load_graph_bin_MAGIC(fn_graph, num, skip);
	skip = 0;

	for j = 1:length(genes)
		gene_cnt = gene_cnt+1;
		genes(j).id = gene_cnt;
	end

	if isempty(genes)
		break;
	end
	%% run prediction for all genes
	for k = 1:length(genes)
		if nargin==5 && ~ismember(genes(k).id, idx)
			continue;
		end
		if 1
			fn_res = sprintf('%s/gene%i.mat', out_dir, genes(k).id);
			if fexist(fn_res)
				x = load(fn_res, 'how');
				if isempty(strfind(x.how, 'Time limit exceeded'))
					continue;
				end
			end
		end
		% get list of bam files for samples
		

		opts = rproc_default_opts();
		opts.identifier = sprintf('MIP_%i_denovo', jobid);
		mem_req = 5000;
		time_req = 1e6;
		fprintf('\rsubmit gene:%i, mult:%i  ', genes(k).id, j)
		opts.addpaths = {fileparts(which('mip_paths')), fileparts(which('create_mip_simple_pair'))};
		opts.resubmit = 2;
		opts.priority = 17;
		opts.mem_req_resubmit = [15000 30000 45000];
		opts.time_req_resubmit = [1e6 1e6 1e6];
		opts.maxjobs = 100;
		
		cov_scale = [];
		for r = 1:size(genes(k).seg_admat, 3)%loop over samples
			cov_scale(r) = max(genes(k).coverage(r, :));
			cov_scale(r) = max([cov_scale(r), max(genes(k).seg_admat(:, :, r))]);
			if cov_scale(r)>0
				sa = genes(k).seg_admat(:,:,r);
				sa(sa>0) = sa(sa>0)/cov_scale(r);
				genes(k).seg_admat(:,:,r) = sa;
				genes(k).coverage(r, :) = genes(k).coverage(r, :)/cov_scale(r);
			end
		end

		if all(cov_scale==0)
			fprintf('no coverage found for gene: %i, %s%s:%i-%i\n', genes(k).id, genes(k).chr, genes(k).strand, genes(k).start, genes(k).stop)
			continue;
		end
		
		PAR.gene = genes(k);
		PAR.cov_scale = cov_scale;
		PAR.coverage = genes(k).coverage;
		PAR.len = genes(k).segments(2, :)-genes(k).segments(1,:)+1;
		PAR.seg_admat = genes(k).seg_admat;
		PAR.predef_trans = [];
		for j=1:size(genes(k).transcripts, 1)
			PAR.predef_trans{j} = genes(k).transcripts(j, :);
		end
		PAR.initial = genes(k).initial;
		PAR.terminal = genes(k).terminal;
		PAR.pair_list = [];
		PAR.insert_size = 0;
		PAR.param = param;
		PAR.strand = genes(k).strand;
		PAR.fn_res = fn_res;
		PAR.C = C;
		if isfield(genes, 'pair_list')
			PAR.pair_list = genes(k).pair_list;
		end

		if 0 % did not change much; slightly worse
			PAR.seg_admat = prune_graph(PAR.seg_admat, PAR.predef_trans, 1e4);
		end

		num_paths = count_paths(double(PAR.seg_admat>-2));
		if num_paths <= 5
			PAR.C.num_transcripts = C.num_transcripts_predef;
		elseif num_paths > 1e4
			PAR.C.num_transcripts = C.num_transcripts*10;
		end

		if nargin<4
			num_samples = size(PAR.seg_admat, 3);
			gpp = load('param/GP_param_samples.mat', 'res_orig', 'names');
			idx = min(size(gpp.res_orig, 1), num_samples);
			for n = 1:length(gpp.names)
				try
					eval(sprintf('PAR.%s;', gpp.names{n}))
					eval(sprintf('PAR.%s = gpp.res_orig(idx, n);', gpp.names{n}))
				catch
					%fprintf('not a valid field: %s\n', gpp.names{n})
				end
			end
			PAR.param.lambda = round(PAR.param.lambda);
		end

		if ~parallel
			create_mip_simple_pair(PAR);
		else
			cnt = cnt+1;
			jobinfo(cnt) = rproc('create_mip_simple_pair', PAR, mem_req, opts, time_req);
		end
	end
end
%if cnt>0
%	[jobinfo,num_crashed] = rproc_wait(jobinfo, 10, 1,-1, 1);
%end
return



