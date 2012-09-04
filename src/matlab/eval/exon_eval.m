function [num, matched, correct] = exon_eval(agene, pgene)

matched = zeros(1, length(agene.exons));
correct = zeros(1, length(pgene.exons));
pairings = zeros(length(agene.exons), length(pgene.exons));

for j = 1:length(agene.exons)
	aexons = agene.exons{j};
	aintrons = [aexons(1:end-1, 2) aexons(2:end, 1)];
	for k = 1:length(pgene.exons)
		pexons = pgene.exons{k};
		pintrons = [pexons(1:end-1, 2) pexons(2:end, 1)];
		if size(aexons, 1)==1 && size(pexons, 1)==1
			matched(j) = matched(j) || overlap(aexons, pexons); 
			correct(k) = correct(k) || overlap(aexons, pexons);
			pairings(j,k) = overlap(aexons, pexons);
		elseif isequal(aintrons, pintrons)
			matched(j) = 1;
			correct(k) = 1;
			pairings(j,k) = 1;
		elseif isequal(size(aintrons), size(pintrons))
			% check for minor offsets
			m1 = aintrons(:)==pintrons(:);
			m2 = aintrons(:)==(pintrons(:)+1);
			m3 = aintrons(:)==(pintrons(:)-1);
			if all(m1 | m2 | m3)
				warning('match with inaccurate exon boundaries');
				matched(j) = 1;
				correct(k) = 1;
				pairings(j,k) = 1;
			end
		end
	end
end

% take care of multiple pairings
[assignment,cost] = munkres((1-pairings));
num = 0; 
for j = 1:length(agene.exons)
	for k = 1:length(pgene.exons)
		num = num + (assignment(j,k)==1 && pairings(j,k)==1);
	end
end

return
function ret=overlap(int1, int2)

	assert(size(int1, 1)==1)
	assert(size(int1, 2)==2)
	assert(size(int2, 1)==1)
	assert(size(int2, 2)==2)

	ret = int1(1)<=int2(1) && int1(2)>=int2(1);
	ret = ret || (int1(2)>=int2(2) && int1(1)<=int2(2));
	ret = ret || (int1(1)>=int2(1) && int1(2)<=int2(2));

return
