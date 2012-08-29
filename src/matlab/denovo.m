function denovo(fn_graph, fn_bam, out_dir)

if nargin==0
	fn_bam = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.paired.sorted.bam';
end
% model selection parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP opt num_transcripts:994.4662   pairs:99.8423   intron:99.6441   exon:99.7242    lambda:1.0000
C.num_transcripts_predef = 1.6;
C.num_transcripts = 500;
C.introns = 100; 
C.pairs = 100; 
C.exon_cov = 100; 
C.num_clusters=10;

% the tol is this value times the insert size 
param.insert_size_tol = 0.5; % optimized together with C.pairs using 2dms

param.loss = 'nb';
param.num_clusters=1;

param.slope_threash = 3.1;

% nois level: this is mean and variance of the poisson that explains the noise
param.lambda = 1; %2;

% no model selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%param.prune_graph = 1;
%param.subsample = 1000;
%param.use_pair=1;
param.number_of_additional_transcripts=5;
param.use_predef_trans=1;
%param.denovo = 0;
param.len_cutoff=120;
param.time_limit=3600;
%param.force=1;
%param.infer_segments = 0;
%param.fn_gp_train = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/nb_fit/GP_train_dat.mat';
%if param.use_pair>0
%	param.max_segment_length = 100;
%else
%	param.max_segment_length = 1e6;
%end

rand('seed', ceil(cputime*10));
jobid = ceil(rand*100);
%genes = genes([genes.id] == 29899);
jobinfo = rproc_empty(0);
cnt = 0;
gene_cnt = 0;

% open binary file
while 1

	num = 100;
	if ~isempty(fn_bam)
		genes = load_graph_bin(fn_graph, num, fn_bam);
	else
		genes = load_graph_bin(fn_graph, num);
	end
	for j = 1:length(genes)
		gene_cnt = gene_cnt+1;
		genes(j).id = gene_cnt;
	end

	if isempty(genes)
		break;
	end
	%% run prediction for all genes
	for k = 1:length(genes)
		if 1
			fn_res = sprintf('%s/gene%i.mat', out_dir, genes(k).id);
			if fexist(fn_res)
				x = load(fn_res, 'how');
				if isempty(strfind(x.how, 'Time limit exceeded'))
					continue;
				end
			end
		end

		opts = rproc_default_opts();
		opts.identifier = sprintf('MIP_%i_denovo', jobid);
		mem_req = 5000;
		time_req = 1e6;
		fprintf('\rsubmit gene:%i, mult:%i  ', genes(k).id, j)
		opts.addpaths = {fileparts(which('mip_paths')), fileparts(which('create_mip_simple_pair'))};
		opts.resubmit = 3;
		opts.priority = 17;
		opts.mem_req_resubmit = [15000 30000 45000];
		opts.time_req_resubmit = [1e6 1e6 1e6];
		
		cov_scale = max(genes(k).coverage);
		cov_scale = max([cov_scale, genes(k).seg_admat(:)']);
		if cov_scale==0
			fprintf('no coverage found for gene: %i\n', k)
			continue;
		end
		genes(k).seg_admat(genes(k).seg_admat>0) = genes(k).seg_admat(genes(k).seg_admat>0)/cov_scale;
		PAR.gene = genes(k);
		PAR.cov_scale = cov_scale;
		PAR.coverage = genes(k).coverage/cov_scale;
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

		if 0%count_paths(double(PAR.seg_admat>-2)) <= 5
			PAR.C.num_transcripts = 0;
		end

		%genes(k).chr
		%genes(k).strand
		%genes(k).start
		%genes(k).stop

		%create_mip_simple_pair(PAR);
	
		cnt = cnt+1;
		jobinfo(cnt) = rproc('create_mip_simple_pair', PAR, mem_req, opts, time_req);
	end
end
[jobinfo,num_crashed] = rproc_wait(jobinfo, 10, 1,-1, 1);
return



