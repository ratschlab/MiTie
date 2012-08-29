function l = load_pred(res_dir, idx)

files = dir([res_dir '/gene*.mat']);

if nargin>1 && max(idx)<=length(files)
	files = files(idx);
end

% load predictions
fprintf('\n');
rm_idx = zeros(1, length(files));
l = [];
cnt = 0;
for j = 1:length(files)
    fprintf('\r load file %i(%i)', j, length(files));
	fname = [res_dir files(j).name];
	if 0
		%for s = 1:4
		%	if ~isempty(strfind(fname, sprintf('multiple%i', s)))
		%		fname2 = strrep(fname, sprintf('multiple%i', s), 'multiple4');
		%	end
		%end
		fname2 = sprintf('/fml/ag-raetsch/nobackup/projects/mip/human_sim//res_multiple4_alt10_C50_sigmoid_pred/%s', files(j).name);
		if fexist(fname2)
			x = load(fname2, 'how');
			if ~isempty(strfind(x.how, 'Time limit exceeded')) 
				continue
			end
		else
			continue
		end
	end
	cnt = cnt+1;
    if isempty(l)
		try
        	l = load(fname);
		catch
			fprintf('failed to load file: %s\n', fname);
			%unix(sprintf('rm %s', fname));

		end
    else
		try
        	x = load(fname);
		catch
			fprintf('failed to load file: %s\n', fname);
			%unix(sprintf('rm %s', fname));
		end
		if ~isempty(setdiff(fieldnames(l), fieldnames(x)))
			diff = setdiff(fieldnames(l), fieldnames(x));
			for f = 1:length(diff)
				x.(diff{f}) = [];
			end
		end
		if ~isempty(setdiff(fieldnames(x), fieldnames(l)))
			diff = setdiff(fieldnames(x), fieldnames(l));
			for f = 1:length(diff)
				l(1).(diff{f}) = [];
			end
		end
		try
			l(cnt) = x;
		catch
			x = orderfields(x, fieldnames(l));
			l(cnt) = x;
		end
    end
end

if ~isfield(l, 'segments')
	for j = 1:length(l)
		l(j).segments = l(j).PAR.gene.segments';
		l(j).gene = l(j).PAR.gene;
		l(j).param = l(j).PAR.param;
		l(j).seg_admat = l(j).PAR.seg_admat;
		l(j).cov_scale = l(j).PAR.cov_scale; 
	end
end

fprintf('\n');

