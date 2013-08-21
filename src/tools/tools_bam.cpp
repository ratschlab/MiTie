
#include "tools_bam.h"
#include "bam_region.h"
#include "bam.h"

vector<Bam_Region*> parse_bam_regions(char* fn_regions)
{
	FILE* fd = fopen(fn_regions, "r");
	if (!fd)
	{
		fprintf(stderr, "tools: Could not open file: %s for reading\n", fn_regions);
	
	}
	int num_bytes = 1000;

	vector<Bam_Region*> regions;
	while (!feof(fd))
	{
		char line[num_bytes]; 
		if (!fgets(line, num_bytes, fd))
		{
			break;
		}
		if (line[0]=='%' || line[0]=='#')
			continue;
		char* chr = new char[100];
		char* strand = new char[100];
		int start = 0;
		int stop = 0;
		int num_read = sscanf(line, "%s\t%s\t%i\t%i", chr, strand, &start, &stop);
		if (num_read!=4)
		{
			fprintf(stderr, "tools: Error parsing line: %s\n", line);
			exit(-2);
		}
		Bam_Region* reg = new Bam_Region(start, stop, chr, strand[0]);
		regions.push_back(reg);
		delete[] strand; 
		delete[] chr;
	}

	fclose(fd);

	return regions;
}

void set_chr_num(Region* reg, bam_header_t* header)
{
	for (int j=0; j<header->n_targets; j++)
	{
		if (strcmp(reg->chr, header->target_name[j])==0)
		{
			reg->chr_num = j;
			return;
		}
	}
	fprintf(stderr, "tools: Did not find chr name in header: %s\n", reg->chr);
	exit(-2);
}

vector<vector<int> > region_overlap(vector<Bam_Region*> regions1, vector<Bam_Region*> regions2)
{
	// compute overlap
	vector<int> starts1;
	vector<int> stops1;
	vector<int> starts2;
	vector<int> stops2;

	for (uint i=0; i<regions1.size(); i++)
	{
		starts1.push_back(regions1[i]->start);
		stops1.push_back(regions1[i]->stop);
	}
	for (uint i=0; i<regions2.size(); i++)
	{
		starts2.push_back(regions2[i]->start);
		stops2.push_back(regions2[i]->stop);
	}

	vector<int> ov = interval_overlap(starts1, stops1, starts2, stops2);

	vector<vector<int> > ov_list(regions1.size());
	for (uint i=0; i<ov.size(); i+=2)
	{
		ov_list[ov[i]].push_back(ov[i+1]);
	}
	return ov_list;
}

