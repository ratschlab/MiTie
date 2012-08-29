function genes = collect_results(res_dir, fn_save, add_weights, mmr)

if nargin<3
	add_weights=0;
end
if nargin<4
	mmr = 1;
end

if isstruct(res_dir)
	l = res_dir;
else
	l = load_pred(res_dir);
end
if nargin<2 && ~isstruct(res_dir)
	fn_save = [res_dir 'res_genes.mat'];
end

%t = cputime;
l = remove_dublicated_transcripts(l);
%rmtrtime = cputime-t

if mmr
	for j = 1:length(l)
		l(j).param.use_predef_trans = 1;
		l(j).gene.predef_trans = l(j).gene.transcripts;
		l(j).gene = rmfield(l(j).gene, 'transcripts');
	end
end

genes = [l.gene];
assert(length(genes)==length(l))
rmfield(l, 'gene');

for j = 1:length(l)

	num_predef = 0;
	if l(j).param.use_predef_trans
		num_predef = size(l(j).weights, 2)-l(j).param.number_of_additional_transcripts;
	end
	weights = sum(l(j).weights, 1)/size(l(j).weights, 1);
	idx = find(cellfun('isempty', l(j).transcripts)==0 & (weights>1e-2 | 1:length(weights)<=num_predef));

	l(j).transcripts = l(j).transcripts(idx);
	l(j).weights = l(j).weights(:, idx);
	if isfield(genes, 'exons')
		genes(j).exons_orig = genes(j).exons;
	else
		genes(j).exons_orig = [];
	end
	genes(j).exons = cell(1, length(l(j).transcripts));
	genes(j).weights = l(j).weights;
	genes(j).cov_scale = l(j).cov_scale;
	for k = 1:length(l(j).transcripts)
		exons = segments2exons(l(j).segments, l(j).transcripts{k});

		%if mmr && genes(j).strand == '+'
		%	exons = exons-1;
		%elseif mmr
		%	exons = exons+1;
		%end
		if k<=length(genes(j).exons_orig)
			%assert(isequal(genes(j).exons{k}, genes(j).exons_orig{k}));
			genes(j).exons{k} = genes(j).exons_orig{k};
			assert(length(genes(j).transcripts)>=k)
		else
			genes(j).exons{k} = exons;
			genes(j).transcripts{k} = sprintf('new_trans%i', k);
		end
	end
end

for j = 1:length(genes)
	genes(j).cds_exons = repmat({[]}, 1, length(genes(j).exons));
end

if ~mmr
	genes = half_open_to_closed(genes);
else
	%genes = genes([genes.strand] == '-');
end
names = {'exons_orig', 'introns', 'splicegraph', 'tss_info', 'tss_conf', 'cleave_info', 'cleave_conf', 'polya_info', 'cleave_conf', 'is_alt', 'is_alt_spliced', 'is_valid', 'transcript_complete', 'is_complete', 'is_correctly_gff3_referenced', 'confgenes_id', 'anno_id', 'cdsStop_conf', 'cdsStop_info', 'polya_conf', 'tis_conf', 'tis_info', 'exons_confirmed', 'cds_exons', 'utr5_exons', 'utr3_exons', 'transcript_info', 'transcript_status', 'transcript_valid', 'transcript_coding'};
for j = 1:length(names)
	if isfield(genes, names{j})
		genes = rmfield(genes, names{j});
	end
end

if ~isstruct(res_dir)
	save(fn_save, 'genes');
end
