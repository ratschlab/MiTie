function l = remove_dublicated_transcripts(l)

	for j = 1:length(l)
		for k = 1:length(l(j).transcripts)
			if isempty(l(j).transcripts{k})
				continue
			end
			for m = k+1:length(l(j).transcripts)
				eq = isequal(l(j).transcripts{k}, l(j).transcripts{m});
				eq = eq || compare_trans_introns(l(j).segments, l(j).transcripts{k}, l(j).transcripts{m});
				if eq
					l(j).transcripts{m} = [];
					l(j).weights(:, k) = l(j).weights(:, k) + l(j).weights(:, m);
					l(j).weights(:, m) = 0;
				end
			end
		end
	end

return

function [eq tr] = compare_trans_introns(segments, tr1, tr2)

tr = [];

if isempty(tr1) || isempty(tr2)
	eq = 0;
	return
end

exons1 = segments2exons(segments, tr1);
exons2 = segments2exons(segments, tr2);
introns1 = [exons1(1:end-1, 2), exons1(2:end, 1)];
introns2 = [exons2(1:end-1, 2), exons2(2:end, 1)];
eq = isequal(introns1, introns2);

if isempty(introns1) && isempty(introns2)
	eq = ~isempty(intersect(tr1, tr2));
	%tr = union(tr1, tr2);
end

return
