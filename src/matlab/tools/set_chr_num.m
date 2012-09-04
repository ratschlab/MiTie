function genes = set_chr_num(genes, chr)

for j = 1:length(genes)
	chr_num = find(strcmp(genes(j).chr, chr));
	assert(length(chr_num)==1)
	genes(j).chr_num = chr_num;
end

