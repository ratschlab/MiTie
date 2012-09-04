function seg_admat = prune_graph(seg_admat, predef_trans, cutoff);

predef_admat = zeros(size(seg_admat));
for j = 1:length(predef_trans)
	seg = find(predef_trans{j});
	for k = 1:length(seg)-1
		predef_admat(seg(k), seg(k+1)) = 1;
	end
end

num_paths = count_paths(double(sum(seg_admat, 3)>-2));
fprintf('prune graph from %i paths ', num_paths)
cnt = 0;
while num_paths>cutoff
	num_paths = count_paths(double(sum(seg_admat, 3)>-2));
	tmp = sum(seg_admat, 3);
	bins = my_prctile(tmp(tmp>=0), 0:1:100);
	[row, col] = find(tmp<=bins(2)&tmp>=0&predef_admat<1);
	r = ceil(rand*length(row));

	seg_admat(row(r),col(r),:) = -2;
	seg_admat(col(r),row(r),:) = -2;
	cnt = cnt+1;
end
fprintf('to %i paths (%i edges removed)\n', num_paths, cnt);

