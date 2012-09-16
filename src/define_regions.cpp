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

void get_regions(vector<Region*>* regions, int* map, int len, char* chr, char strand, int map_resolution, bool shrink, bool cut)
{
	bool in_reg = false;
	int start = -1;
	int stop = -1;
	for (int i=0; i<len; i++)
	{
		if (!in_reg && map[i]>0)
		{
			start = i*map_resolution;
			stop = -1;
			in_reg = true;
		}
		else if (in_reg && map[i]==0)
		{
			stop = i*map_resolution;
			if (shrink)
			{
				//fprintf(stdout, "start:%i, stop:%i", start, stop);
				// move start
				for (int st=start/map_resolution; st<len && map[st]<2; st++)
					start = st*map_resolution;
				// move stop
				for (int st=i; st>start/map_resolution && map[st]<2; st--)
					stop = st*map_resolution;
				//fprintf(stdout, "  shrink-> start:%i, stop:%i\n", start, stop);
			}
			if (cut)
			{
				bool print = false;
				if (start==16190560 && stop==16245960)
				{
					print = true;
				}
				int low = 0;
				int new_stop = 0;
				int thresh = 2;
				for (int st=start/map_resolution+20; st<stop/map_resolution; st++)
				{
					if (st*map_resolution-start>200000&&stop-st*map_resolution>200000 && thresh == 2)
					{
						thresh = 4;
						printf("long region: %i, %i increase cut threshold to %i\n", start, stop, thresh);
					}
					else 
						thresh = 2;

					if (low==0 && map[st]<thresh)
					{
						new_stop = st*map_resolution;
						low++;
						if (print)
							printf("1: start: %i, stop: %i, new_stop: %i\n", start, stop, new_stop);
					}
					else if (map[st]>=thresh && low>50)
					{
						if (print)
							printf("2: start: %i, stop: %i, new_stop: %i, low:%i\n", start, stop, new_stop, low);

						// this is the end of a large lowly covered region
						assert(new_stop>0);
						// add region
						char* pchr = new char[strlen(chr)+1];
						strcpy(pchr, chr);
						Region* reg = new Region(start, new_stop, pchr, strand);
						regions->push_back(reg);

						// find next region
						start = st*map_resolution;
						new_stop = 0;
						low = 0;
					}
					else if (map[st]>=thresh)
					{
						// this is the end of a short lowly covered region
						new_stop = 0;
						low = 0;
					}
					else 
					{
						// within a lowly covered region
						low++;
					}
				}
			}

			char* pchr = new char[strlen(chr)+1];
			strcpy(pchr, chr);
			Region* reg = new Region(start, stop, pchr, strand);
			regions->push_back(reg);
			start = -1;
			stop = -1;
			in_reg = false;
		}
	}
}

struct Config
{
	vector<char*> bam_files;
	bool strand_specific;
	bool shrink;
	bool cut_regions; 
	FILE* fd_out;
	int max_len;
};

void parse_args(int argc, char** argv,  Config* c)
{
	if (argc<2)
	{
		fprintf(stderr, "Usage: %s [-o <fn_out>] fn_bam1 [fn_bam2] ...", argv[0]);
		exit(1);
	}
	c->strand_specific = true;
	c->shrink = false;
	c->cut_regions = false;
	c->max_len = 200000;

	c->fd_out = stdout;
    for (int i = 1; i < argc; i++)  
    {
        if (strcmp(argv[i], "-o") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option -o\n") ;
                exit(1);
            }
            i++;
			c->fd_out = fopen(argv[i], "w+");
			if (!(c->fd_out))
			{
				fprintf(stderr, " Could not open file: %s for writing.\n", argv[2]);
                exit(1);
			}
        }
	    else if (strcmp(argv[i], "--max-len") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --max-len\n") ;
                exit(1);
            }
            i++;
			c->max_len = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "--no-strand") == 0)
		{
			c->strand_specific = false;
		}
        else if (strcmp(argv[i], "--shrink") == 0)
		{
			c->shrink = true;
		}
        else if (strcmp(argv[i], "--cut-regions") == 0)
		{
			c->cut_regions = true;
		}
		else
		{
			c->bam_files.push_back(argv[i]);
		}	
	}
	assert(c->bam_files.size()>0);
}

int main(int argc, char* argv[])
{
	Config c;
	parse_args(argc, argv, &c);

	bamFile fd1 = bam_open(c.bam_files[0], "r");
	if (fd1==0)
	{
		fprintf(stderr, "[%s] Could not open bam file: %s", argv[0], c.bam_files[0]);
		exit(-1);
	}
	bam_header_t* header = bam_header_read(fd1);
	if (header == 0)
	{
		fprintf(stderr, "[%s] Invalid BAM header.", argv[0]);
		exit(-1);
	}


	printf("number of chromosomes in header:%i\n", header->n_targets);
	// initialize maps
	int res = 10;
	int** maps = new int*[header->n_targets*2];
	for (int i=0; i<header->n_targets; i++)
	{
		int len = header->target_len[i];
		maps[2*i] = new int[len/res];
		maps[2*i+1] = new int[len/res];
		memset(maps[2*i], 0, len/res*sizeof(int));
		memset(maps[2*i+1], 0, len/res*sizeof(int));
	}

	bam1_t* b = (bam1_t*) calloc(1, sizeof(bam1_t));
	bam1_core_t* core = &b->core;

	for (int i=0; i<c.bam_files.size(); i++)
	{

		printf("Running on bam file: %s\n", c.bam_files[i]);
		bamFile fd = bam_open(c.bam_files[i], "r");
		if (fd==0)
		{
			fprintf(stderr, "[%s] Could not open bam file: %s", argv[0], c.bam_files[i]);
			exit(-1);
		}
		header = bam_header_read(fd);
		if (header == 0)
		{
			fprintf(stderr, "[%s] Invalid BAM header.", argv[0]);
			exit(-1);
		}

		int last_tid=-1;
		int last_coor=0;
		// maps for exon, intron and pair coverage for each strand

		int cnt=0;
		int discarded = 0;
		int len;
		char* chr; 
		int num_read=1;
		while ((num_read = bam_read1(fd, b))>=0)
		{
			cnt++;
			if (last_tid < core->tid || (last_tid >= 0 && core->tid < 0))
			{
				// change of chromosomes
				last_tid = core->tid;
				len = header->target_len[core->tid];
				chr = header->target_name[core->tid];

				//if (strcmp(chr, "3")==0)
				//{
				//	int start = 16190560;
				//	int stop = 16190560+37000;
				//	int map_resolution = 10;
				//	printf("write map flat files: %i\n", last_tid);
				//	FILE* fd = fopen("/tmp/low_region_w.flat", "w");
				//	for (int d=start/map_resolution; d<stop/map_resolution; d++)
				//		fprintf(fd, "%i\n", maps[2][d]);
				//	fclose(fd);
				//	fd = fopen("/tmp/low_region_c.flat", "w");
				//	for (int d=start/map_resolution; d<stop/map_resolution; d++)
				//		fprintf(fd, "%i\n", maps[3][d]);
				//	fclose(fd);

				//}
				printf("running on region: %s (%i)\n", chr, len);
			} 
			else if ((uint32_t)last_tid > (uint32_t)core->tid)
			{
				fprintf(stderr, "the alignment is not sorted (%s): %d-th chr > %d-th chr\n",
						bam1_qname(b), last_tid+1, core->tid+1);
				return NULL;
			}
			else if ((int32_t)core->tid >= 0 && last_coor > core->pos) 
			{
				fprintf(stderr, "the alignment is not sorted (%s): %u > %u in %d-th chr\n",
						bam1_qname(b), last_coor, core->pos, core->tid+1);
				return NULL;
			}
			CRead* r = new CRead();
			parse_cigar(b, r, header);
			int from = r->start_pos;
			int to = r->get_last_position();
			int map_idx = 2*core->tid;
			if (r->strand[0]=='-' && c.strand_specific)
				map_idx++;

			assert(from>=0);
			if (to>=len)
			{
				fprintf(stderr, "Warning: read alignment exceeds contig length: %i>=%i\n", to, len);
				discarded++;
				delete r;
				last_coor = b->core.pos;
				continue;
			}
			assert(to>from);

			if (to-from>c.max_len)
			{
				discarded++;
				delete r;
				last_coor = b->core.pos;
				continue;
			}

			for (int j=from/res; j<to/res; j++)
				maps[map_idx][j]++;

			// add mate pair information to the map
			if (core->mtid==core->tid && core->mpos-core->pos>0 && core->mpos-core->pos<200000)
			{
				int mate_start = core->mpos;
				for (int j=from/res; j<mate_start/res; j++)
					maps[map_idx][j]++;
			
			}
			delete r;
			last_coor = b->core.pos;
		}
		bam_close(fd);
		printf("num reads: %i (%i discarded (length>%i))\n", cnt, discarded, c.max_len);
	}

	// get regions from maps
	vector<Region*> regions;
	for (int i=0; i<header->n_targets; i++)
	{

		int len = header->target_len[i]/res;
		char* chr = header->target_name[i];

		get_regions(&regions, maps[i*2]  , len, chr, '+', res, c.shrink, c.cut_regions);
		if (c.strand_specific)
			get_regions(&regions, maps[i*2+1], len, chr, '-', res, c.shrink, c.cut_regions);
	}

	// filter regions
	vector<Region*> regions_filtered;
	for (int i=0; i<regions.size(); i++)
	{
		if (regions[i]->stop-regions[i]->start>200)
			regions_filtered.push_back(regions[i]);
	}

	write_regions(regions_filtered, c.fd_out);
	printf("\tnumber regions: %u\n", (int)regions.size());
	printf("\tnumber regions filtered: %u\n", (int)regions_filtered.size());

	// cleanup 
	for (int i=0; i<regions.size(); i++)
		delete regions[i];
	bam_destroy1(b);
	for (int i=0; i<2*header->n_targets; i++)
		delete[] maps[i];
	delete[] maps;
	bam_header_destroy(header);
	bam_close(fd1);
	fclose(c.fd_out);
}
