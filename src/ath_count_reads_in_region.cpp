#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <algorithm>
#include <assert.h>
#include "bam.h"
#include "region.h"
#include "get_reads_direct.h"
#include <fstream>
#include "tools.h"

//includes for samtools
////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "razf.h"
#include "bgzf.h"
#include "razf.h"
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	if (argc<4)
	{
			fprintf(stderr, "Usage: %s <fn_regions> <fn_out> <fn_bam1> <fn_bam2> ...\n", argv[0]);
			return -2;
	}
	char* fn_regions = argv[1];	
	char* fn_out = argv[2];	
	vector<char*> bam_files; 
    for (int i = 3; i < argc; i++)  
    {
		bam_files.push_back(argv[i]);
	}
	FILE* fd_out = fopen(fn_out, "w+");
	if (!fd_out)
	{
		fprintf(stderr, "[%s] Could not open file: %s for writing\n", argv[0], fn_out);
		return -2;
	}

	bamFile fd = bam_open(bam_files[0], "r");

	bam_header_t* header = bam_header_read(fd);
	if (header == 0)
	{
		fprintf(stderr, "[%s] Invalid BAM header.\n", argv[0]);
		return -2;
	}

	FILE* fd_null = fopen("/dev/null", "w");
	vector<Region*> regions = parse_regions(fn_regions);
	printf("read %i regions from file: %s\n", (int)regions.size(), fn_regions);

	int** read_counts = new int*[regions.size()];
	for (int i=0; i<regions.size(); i++)
		read_counts[i] = new int[bam_files.size()];

	for (int b=0; b<bam_files.size(); b++)
	{
		printf("starting with file %s\n", bam_files[b]);
		for (int i=0; i<regions.size(); i++)
		{
			printf("\r%i (%i) regions done", i, (int)regions.size());
			// get reads for the large region
			int num_bam = 1;
			int intron_len_filter = 200000;
			int filter_mismatch = 10;
			int exon_len_filter = 1;
			bool mm_filter = true;
			
			regions[i]->fd_out = fd_null;
			regions[i]->strand = '0';
			regions[i]->get_reads(&bam_files[b], num_bam, intron_len_filter, filter_mismatch, exon_len_filter, mm_filter);
			regions[i]->strand = '+';

			int num_reads=0;
			for (int r=0; r<regions[i]->reads.size(); r++)
				if (regions[i]->reads[r]->mismatches<1)
					num_reads++;
			
			read_counts[i][b] = num_reads;
			// cleanup
			regions[i]->clear_reads();
		}
	}
	printf("read %i regions from file\n", (int)regions.size());

	// output
	for (int i=0; i<regions.size(); i++)
	{
		fprintf(fd_out, "%s\t%c\t%i\t%i", regions[i]->chr, regions[i]->strand, regions[i]->start, regions[i]->stop);	
		for (int b=0; b<bam_files.size(); b++)
		{
			fprintf(fd_out, "\t%i", read_counts[i][b]);
		}
		fprintf(fd_out, "\n");
	}
	// cleanup
	for (int i=0; i<regions.size(); i++)
		delete regions[i];
	
	for (int i=0; i<regions.size(); i++)
		delete[] read_counts[i];
	delete[] read_counts;

	fclose(fd_out);
	bam_close(fd);
	bam_header_destroy(header);
}
