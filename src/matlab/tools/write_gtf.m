function write_gtf(genes, gtf_fname, source)
%	write_gtf(genes, gtf_fname, source)
%
% -- input --
% genes: struct defining genes with start, stops (0-based), exons etc.
% gtf_fname: name of GTF file
% source: source of annotation


fprintf('creating file %s...\n', gtf_fname);
[fd msg] = fopen(gtf_fname, 'w+');
disp(msg);

for g = 1:length(genes),
	fprintf('writing gene %i...\r', g);	
	%if genes(g).is_valid==0, continue; end
	gene = genes(g);
	type	= 'gene';
	score = '.';
	phase = '.';
	for t = 1:length(gene.transcripts),
		%if gene.transcript_valid(t)==0, continue; end
		type = 'transcript';
		score = '.';
		phase = '.';
		start = min(gene.exons{t}(:,1));
		stop = max(gene.exons{t}(:,2));
		% transcript
		attr_str = sprintf('gene_id "%s"; transcript_id "%s";', gene.name, gene.transcripts{t});
		fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
		% exons
		attr_str = sprintf('gene_id "%s"; transcript_id "%s";', gene.name, gene.transcripts{t});
		exon_type = {'exons'};
		gff_types	= {'exon'};
		for tt = 1:length(exon_type),
			exons = gene.(exon_type{tt}){t};
			for e = 1:size(exons,1),
				type = gff_types{tt};
				score = '.';
				phase = '.';
				start = exons(e,1);
				stop = exons(e,2);
				fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
			end
		end
	end
end

fclose(fd);
