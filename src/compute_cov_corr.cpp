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
#include "tools.h"
//#include "genome.h"
//#include "gene.h"
//#include "gene_tools.h"
//#include "infer_genes.h"
//#include "Config.h"

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


struct Config
{
	vector<char*> bam_files;
	bool strand_specific;
};

void parse_args(int argc, char** argv,  Config* c)
{
	if (argc<2)
	{
		fprintf(stderr, "Usage: %s [-o <fn_out>] fn_bam1 [fn_bam2] ...", argv[0]);
		exit(1);
	}
	c->strand_specific = true;
    for (int i = 1; i < argc; i++)  
    {
        if (strcmp(argv[i], "--no-strand") == 0)
		{
			c->strand_specific = false;
		}
		else
		{
			c->bam_files.push_back(argv[i]);
		}	
	}
	assert(c->bam_files.size()==2);
}

int main(int argc, char* argv[])
{
	Config c;
	parse_args(argc, argv, &c);

	bamFile fd1 = bam_open(c.bam_files[0], "r");
	bamFile fd2 = bam_open(c.bam_files[1], "r");
	if (fd1==0)
	{
		fprintf(stderr, "[%s] Could not open bam file: %s", argv[0], c.bam_files[0]);
		exit(-1);
	}
	if (fd2==0)
	{
		fprintf(stderr, "[%s] Could not open bam file: %s", argv[0], c.bam_files[1]);
		exit(-1);
	}

	bam_header_t* header1 = bam_header_read(fd1);
	bam_header_t* header2 = bam_header_read(fd2);
	if (header1 == 0 || header2 == 0)
	{
		fprintf(stderr, "[%s] Invalid BAM header.", argv[0]);
		exit(-1);
	}


	int res=10;
	printf("number of chromosomes in header:%i\n", header1->n_targets);
	int num_chr = header1->n_targets;
	unsigned long int len[num_chr+1];
	len[0] = 0;
	// initialize maps
	for (int i=0; i<num_chr; i++)
	{
		len[i+1] = len[i]+header1->target_len[i];
		//printf("i: %i, len[i]: %lu, len[i+1]: %lu, chr_len: %i\n",i, len[i], len[i+1], header1->target_len[i]);
	}
	printf("total: %lu, res: %i\n", len[num_chr], res);
	int* map1 = new int[len[num_chr]/res];
	int* map2 = new int[len[num_chr]/res];
	memset(map1, 0, len[num_chr]/res*sizeof(int));
	memset(map2, 0, len[num_chr]/res*sizeof(int));

	bam1_t* b = (bam1_t*) calloc(1, sizeof(bam1_t));
	bam1_core_t* core = &b->core;


	for (int i=0; i<2; i++)
	{
		bamFile fd = fd1;
		if (i>0) fd = fd2;

		int cnt=0;
		while (bam_read1(fd, b)>=0)
		{
			cnt++;
			CRead* r = new CRead();
			parse_cigar(b, r, header1);
			int from = (r->start_pos+len[core->tid])/res;
			int to = (r->get_last_position()+len[core->tid])/res;
			//if (r->strand[0]=='-' && c.strand_specific)
			//	map_idx++;

			assert(from>=0);
			assert(to>from);
			assert(to<len[num_chr]/res);

			if (i==0)
			{
				for (int j=from; j<to; j++)
					map1[j]++;
			}
			else
			{
				for (int j=from; j<to; j++)
					map2[j]++;
			}
			delete r;
		}
		printf("num reads: %i \n", cnt);
	}

	float corr = pearson_corr(map1, map2, len[num_chr]/res);

	printf("pearson correlation: %.5f\n", corr);

	// cleanup 
	bam_destroy1(b);
	bam_header_destroy(header1);
	bam_header_destroy(header2);
	bam_close(fd1);
	bam_close(fd2);
}
