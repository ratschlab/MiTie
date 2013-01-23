

dir_ = '/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/MiTie/pred_Chr1+/'
files = dir([dir_ 'gene*.mat']);

starts = zeros(1, 100000);
stops = zeros(1, 100000);
for j = 1:length(files)
	l = load([dir_ files(j).name]);
	idx = find(starts==l.PAR.gene.start & stops == l.PAR.gene.stop);
	if ~isempty(idx)
		unix(sprintf('echo %s', [dir_ files(j).name]));
		unix(sprintf('rm %s', [dir_ files(j).name]));
	else
		starts(j) = l.PAR.gene.start;
		stops(j) = l.PAR.gene.stop;
	end
end
