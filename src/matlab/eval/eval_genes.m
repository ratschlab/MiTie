function [SN, SP, F, matched, correct, rank, ab_rank] = eval_genes(agenes, pgenes)

if ~isfield(pgenes, 'chr_num') || ~isfield(agenes, 'chr_num')
	chr = unique({agenes.chr});
	chr2 = unique({pgenes.chr});
	if isempty(intersect(chr, chr2))
		keyboard
	end
	chr = unique([chr, chr2]);
	agenes = set_chr_num(agenes, chr);
	pgenes = set_chr_num(pgenes, chr);
end


% split predicted genes into transcripts

num_trans = 0;
for j = 1:length(pgenes)
	if ~isfield(pgenes, 'weights') || isempty(pgenes(j).weights)
		pgenes(j).weights = zeros(1, length(pgenes(j).exons));
	end
	num_trans = num_trans+length(pgenes(j).exons);
end

% allocate memory
ptrans(num_trans) = pgenes(1);
cnt = 0;
for j = 1:length(pgenes)
	for k = 1:length(pgenes(j).exons)
		cnt = cnt+1;
		ptrans(cnt) = pgenes(j);
		ptrans(cnt).start = pgenes(j).exons{k}(1,1);
		ptrans(cnt).stop = pgenes(j).exons{k}(end,2);
		ptrans(cnt).exons = pgenes(j).exons(k);
		ptrans(cnt).transcripts = pgenes(j).transcripts(k);
		ptrans(cnt).weights = pgenes(j).weights(k);
	end
end

consider_strands = 1;
idx = find_overlapping_regions(agenes, ptrans, consider_strands);

% 
num_correct = 0;
num_annotated = 0;
num_predicted = 0;
matched = -ones(1, 1e4);
correct = -ones(1, 1e4);
rank    = -ones(1, 1e4);
ab_rank = -ones(1, 1e4);
cnt1 = 0;
cnt2 = 0;
for j = 1:length(agenes)
	aidx = find(idx(:, 1)==j);	
	pidx = idx(aidx, 2);

	if isempty(pidx)
		num_annotated = num_annotated+length(agenes(j).exons);
		continue
	end
	pgene = ptrans(pidx(1));
	tcnt = 0;
	for k = 1:length(pidx)
		idxm = find(idx(:, 2)==pidx(k));
		if 0%length(idxm)>1
			% this is a gene merge
			warning('found a gene merge: ignoring transcript');
			keyboard
			continue

		end
		tcnt = tcnt+1;
		pgene.exons(tcnt) = ptrans(pidx(k)).exons;
		pgene.transcripts(tcnt) = ptrans(pidx(k)).transcripts;
		pgene.weights(tcnt) = ptrans(pidx(k)).weights;
	end

	[num, m, c] = exon_eval(agenes(j), pgene);

	if cnt1+length(m)>length(matched)
		matched = [matched -ones(size(matched))];
		rank = [rank -ones(size(rank))];
	end
	if cnt2+length(c)>length(correct)
		correct = [correct -ones(size(correct))];
		ab_rank = [ab_rank -ones(size(ab_rank))];
	end

	assert(length(m)==length(agenes(j).exons));
	num_correct = num_correct+num;
	num_annotated = num_annotated+length(agenes(j).exons);
	num_predicted = num_predicted+length(pgene.exons);
	matched(cnt1+1:cnt1+length(m)) = m;
	if isfield(agenes, 'rank')
		rank(cnt1+1:cnt1+length(m)) = agenes(j).rank;
	elseif isfield(agenes, 'expr_orig')
		for k=1:length(agenes(j).exons)
			rank(cnt1+k) = sum(agenes(j).expr_orig>agenes(j).expr_orig(k))+1; 
		end
		assert(k==length(m))
	else
		rank(cnt1+1:cnt1+length(m)) = 0;
	end
	cnt1 = cnt1+length(m);
	correct(cnt2+1:cnt2+length(c)) = c;
	for k=1:length(pgene.exons)
		ab_rank(cnt2+k) = sum(pgene.weights>pgene.weights(k))+1;
	end
	assert(k==length(c))
	cnt2 = cnt2+length(c);

end
assert(sum(rank==-1)==sum(matched==-1))
assert(sum(ab_rank==-1)==sum(correct==-1))

matched(cnt1+1:end) = [];
rank(cnt1+1:end) = [];
correct(cnt2+1:end) = [];
ab_rank(cnt2+1:end) = [];

SN = num_correct/num_annotated;
SP = num_correct/num_predicted;
F = 2*SN*SP/(SN+SP);

return

