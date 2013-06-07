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
#include "gtf_tools.h"
#include "tools.h"
#include "file_stats.h"

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
	char* fn_gtf;
	char* fn_regions;
	char* fn_out;
	float seg_filter;
	float tss_pval;
	bool split_chr;
	bool strand_specific;
	bool reads_by_chr;
	bool mm_filter;
	int intron_len_filter;
	int filter_mismatch;
	int exon_len_filter;
	int region_filter;
	int gtf_offset;
	int max_junk;
};

int parse_args(int argc, char** argv,  Config* c)
{
	if (argc<4)
	{
			fprintf(stdout, "Usage: %s <fn_out> [options] <fn_bam1> <fn_bam2> ...\n", argv[0]);
			fprintf(stdout, "\n");
			fprintf(stdout, "options:\n");
			fprintf(stdout, "\t--regions\t\t(file name) specify regions flat file (e.g. output of define_regions)\n");
			fprintf(stdout, "\t--few-regions \t\t(flag) load data not for whole chromosome, but for each region separate \n");
			fprintf(stdout, "\t--max-junk \t\t(default 2e7) if flag --few-regions not set, load reads for junk at once\n");
			fprintf(stdout, "\t\t\t\t -> more efficient if there are only a few region to predict on\n");
			fprintf(stdout, "\t--region-filter\t\t(default 1000) discard regions with less reads \n");
			fprintf(stdout, "\t--gtf\t\t\t(file name) specify gtf file\n");
			fprintf(stdout, "\t--gtf-offset\t\t(default 10000) region arround annotated gene to be included in search for alternative transcripts\n");
			fprintf(stdout, "\t--seg-filter\t\t(default 0.05) segments are filtered out if fraction of covered nucleotides is below\n");
			fprintf(stdout, "\t--tss-tts-pval\t\t(default 0.0001) p-value cutoff for transcription start and transcription termination site discovery\n");
			fprintf(stdout, "\t--split-chr\t\t(flag) generate one output file for each chromosom\n");
			fprintf(stdout, "\t\n");
			fprintf(stdout, "read related options:\n");
			fprintf(stdout, "\t--min-exonic-len\t(default 0) minimal number of aligned positions on each side of an intron\n");
			fprintf(stdout, "\t--mismatches\t\t(default 10) maximal number of mismatches\n");
			fprintf(stdout, "\t--best\t\t\t(flag) filter for best alignment\n");
			fprintf(stdout, "\t--no-strand \t\t(flag) data are not strand specific \n");
			fprintf(stdout, "\t\t\n");
			return -1;
	}

	// defaults
	c->fn_gtf = NULL;
	c->fn_regions = NULL;	
	c->fn_out = argv[1];	
	c->gtf_offset = 10000;
	c->strand_specific = true;
	c->reads_by_chr = true;
	c->intron_len_filter = 200000;
	c->filter_mismatch = 10;
	c->exon_len_filter = 0;
	c->mm_filter = false;
	c->region_filter = 1000;
	c->seg_filter = 0.05;
	c->tss_pval = 1e-4;
	c->split_chr=false;
	c->max_junk = 2e7;

    for (int i = 2; i < argc; i++)  
    {
    	if (strcmp(argv[i], "--no-strand") == 0)
		{
			c->strand_specific = false;
		}
		else if (strcmp(argv[i], "--split-chr") == 0)
        {
			c->split_chr = true;
        }
	    else if (strcmp(argv[i], "--best") == 0)
        {
			c->mm_filter = true;
        }
	    else if (strcmp(argv[i], "--max-junk") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --max-junk\n") ;
                return -1;
            }
            i++;
			c->max_junk = atoi(argv[i]);

        }
	    else if (strcmp(argv[i], "--gtf-offset") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --gtf-offset\n") ;
                return -1;
            }
            i++;
			c->gtf_offset = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--min-exonic-len") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --min-exonic-len\n") ;
                return -1;
            }
            i++;
			c->exon_len_filter = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--mismatches") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --mismatches\n") ;
                return -1;
            }
            i++;
			c->filter_mismatch = atoi(argv[i]);
        }
    	else if (strcmp(argv[i], "--few-regions") == 0)
		{
			c->reads_by_chr = false;
		}
	    else if (strcmp(argv[i], "--region-filter") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --region-filter\n") ;
                return -1;
            }
            i++;
			c->region_filter = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--regions") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --regions\n") ;
                return -1;
            }
            i++;
			c->fn_regions = argv[i];
        }
	    else if (strcmp(argv[i], "--gtf") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --gtf\n") ;
                return -1;
            }
            i++;
			c->fn_gtf = argv[i];
        }
	    else if (strcmp(argv[i], "--seg-filter") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --seg-filter\n") ;
                return -1;
            }
            i++;
			c->seg_filter = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--tss-tts-pval") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --tss_tts_pval\n") ;
                return -1;
            }
            i++;
			c->tss_pval = atof(argv[i]);
        }
		else
		{
			c->bam_files.push_back(argv[i]);
		}	
	}
	assert(c->bam_files.size()>0 || c->fn_gtf);
	if (!c->fn_regions && !c->fn_gtf)
	{
		printf("either --gtf, or --fn-regions, or both have to be specified\n");
	}
	return 0;
}

void process_gtf_region(Region* region)
{
	if (region->transcripts.size()==0)
	{
		printf("no transcript found\n");
		return;
	}
	int num_pos = region->stop-region->start+1;
	//printf("allocate mem array of length %i\n", num_pos);
	int* map = new int[num_pos];
	if (!map)
	{
		fprintf(stderr, "out of mem\n");
		exit(2);
	}

	memset(map, 0, num_pos*sizeof(int));

	// compute map
	for (int i=0; i<region->reads.size(); i++)
	{
		CRead* r = region->reads[i];
		int from = r->start_pos - region->start;
		int to = r->get_last_position() - region->start;


		// discard reads that span out of the region
		if (from<0)
		{
			//fprintf(stdout, "from: %i, to: %i\n", r->start_pos, r->get_last_position());
			continue;
		}

		if (to>num_pos)
		{
			//fprintf(stdout, "from: %i, to: %i\n", r->start_pos, r->get_last_position());
			continue;
		}

		assert(from<num_pos);
		if (to<0)
		{
			fprintf(stderr, "read last position (%i) smaller than region start (%i)\n", r->get_last_position(), region->start );
			exit(-1); 
			//continue;
		}

		for (int j=from; j<to; j++)
			map[j]++;
	}
	if (false)
	{
		FILE* fd = fopen("/fml/ag-raetsch/home/jonas/tmp/map", "w");
		for (int i=0; i<num_pos; i++)
			fprintf(fd, "%i\n", map[i]);
		fclose(fd);
	}

	int start = 1e9; 
	int stop = 0;
	for (int i=0; i<region->transcripts.size(); i++)
	{
		start = std::min(start, region->transcripts[i].front().first);
		stop = std::max(stop, region->transcripts[i].back().second);
	}
	if (!(start>=region->start))
	{
		//fprintf(stderr, "Start of transcript (%i) does not fit region start (%i)\n", start, region->start);
		exit(0);
	}
	assert(start<region->stop);
	int gap = 0;
	int new_start = region->start;
	for (int i=start; i>region->start; i--)
	{
		if (map[i-region->start]<1)
		{
			if (gap==0)
				new_start = i;
			gap++;
		}
		else
		{
			gap = 0;
		}
		
		if (gap>100)
		{
			//printf("move region start from %i to %i (anno: %i)\n", region->start, new_start, start);
			break;
		}
	}

	assert(stop>=region->start);
	if (stop>=region->stop)
	{
		fprintf(stderr, "Transcript end (%i) does not match region stop (%i)\n", stop, region->stop);
	}
	int new_stop = 0;
	gap = 0;
	for (int i=stop; i<region->stop; i++)
	{
		if (map[i-region->start]<2)
		{
			if (gap==0)
				new_stop = i;
			gap++;
		}
		else
		{
			gap = 0;
		}
		
		if (gap>100)
		{
			//printf("move region start from %i to %i (anno: %i) gap: %i\n", region->stop, new_stop, stop, gap);
			region->stop = new_stop;
			break;
		}

	}
	region->start = new_start;

	delete[] map;
}
bool compare_chr_and_strand(const Region* reg1, const Region* reg2)
{
	if (reg1->chr_num<reg2->chr_num)
		return true; 
	if (reg1->chr_num>reg2->chr_num)
		return false;

	if (reg1->strand=='+' && reg2->strand=='-')
		return true;
	if (reg1->strand=='-' && reg2->strand=='+')
		return false;
	
	if (reg1->start<reg2->start)
		return true;
	else 
		return false;
}

int main(int argc, char* argv[])
{
	Config c;
	int ret = parse_args(argc, argv, &c);
	if (ret!=0)
		return ret;

	std::ofstream* ofs = new std::ofstream();
	ofs->open(c.fn_out, std::ios::binary);

	if (!ofs->is_open())
	{
		fprintf(stderr, "[%s] Could not open file: %s for writing\n", argv[0], c.fn_out);
		return -2;
	}
	
	FILE* fd_null = fopen("/dev/null", "w");

	vector<Region*> gtf_regions;
	if (c.fn_gtf && c.bam_files.size()==0 )
	{
		// this is used e.g. when we only want to quantify known transcripts
		// compute graph for annotation only
		const char* format = determine_format(c.fn_gtf);
		printf("loading regions form %s file: %s\n", format, c.fn_gtf);
		if (strcmp(format, "gtf")==0)
			gtf_regions = parse_gtf(c.fn_gtf);
		else if (strcmp(format, "gff3")==0)
			gtf_regions = parse_gff(c.fn_gtf);
		else
		{
			printf("could not determine format of annotation file: %s\n", c.fn_gtf);
			exit(-1);
		}
		printf("number of regions from gtf file: %i\n", (int) gtf_regions.size());

		for (uint i=0;i<gtf_regions.size(); i++)
		{
			//gtf_regions[i]->fd_out = stdout;
			gtf_regions[i]->fd_out = fd_null;
			gtf_regions[i]->generate_segment_graph(c.seg_filter, c.tss_pval);

			// write region in binary file
			gtf_regions[i]->write_binary(ofs);
		}

		// cleanup
		ofs->close();
		delete ofs;
		return 0;
	}

	bamFile fd = bam_open(c.bam_files[0], "r");
	if (fd==0)
	{
		printf("failed to open bam file: %s\n", c.bam_files[0]);
		return -2;
	}
	bam_header_t* header = bam_header_read(fd);
	if (header == 0)
	{
		fprintf(stderr, "[%s] Invalid BAM header.\n", argv[0]);
		return NULL;
	}

	vector<Region*> regions;
	if (c.fn_regions)
	{
		printf("loading regions from flat file: %s\n", c.fn_regions);
		regions = parse_regions(c.fn_regions);
		printf("number of regions from flat file: %i\n", (int) regions.size());
	}
	int num_reg = regions.size();


	char strand_prev;
	vector<CRead*>::iterator curr; 
	Region* reg = NULL;
	int last_stop=0;
	char* chr_prev = (char*) "xxx";


	if (c.fn_gtf)
	{
		const char* format = determine_format(c.fn_gtf);
		printf("loading regions form %s file: %s\n", format, c.fn_gtf);
		if (strcmp(format, "gtf")==0)
			gtf_regions = parse_gtf(c.fn_gtf);
		else if (strcmp(format, "gff3")==0)
			gtf_regions = parse_gff(c.fn_gtf);
		else
		{
			printf("could not determine format of annotation file: %s\n", c.fn_gtf);
			exit(-1);
		}
		printf("number of regions from gtf file: %i\n", (int) gtf_regions.size());

		if (c.reads_by_chr)
		{
			for (int i=0; i<gtf_regions.size(); i++)
			{
				set_chr_num(gtf_regions[i], header);
			}
			sort(gtf_regions.begin(), gtf_regions.end(), compare_chr_and_strand);
		}

		printf("process gtf regions ... \n");
		for (int i=0; i<gtf_regions.size() && c.gtf_offset>0; i++)
		{
			printf("\rprocess gtf region %i (%i)", i, (int) gtf_regions.size());
			set_chr_num(gtf_regions[i], header);

			int chr_len = header->target_len[gtf_regions[i]->chr_num];
			gtf_regions[i]->start = std::max(0, gtf_regions[i]->start-c.gtf_offset);
			gtf_regions[i]->stop = std::min(chr_len-1, gtf_regions[i]->stop+c.gtf_offset);

			// shrink start and stop according to read coverage	
			int num_bam = c.bam_files.size();
			if (!c.strand_specific)
				gtf_regions[i]->strand = '0';
			
			gtf_regions[i]->fd_out = fd_null;
			if (c.reads_by_chr)// more fast if many regions are considered
			{
				bool get_reads = false;
				if (strcmp(gtf_regions[i]->chr, chr_prev)!=0 || (gtf_regions[i]->strand != strand_prev && c.strand_specific) || reg->stop<gtf_regions[i]->stop)
				{
					get_reads = true;
					chr_prev = gtf_regions[i]->chr;
					strand_prev = gtf_regions[i]->strand;
					printf("starting with chr: %s%c\n", chr_prev, strand_prev);
					delete reg;
					reg = new Region(gtf_regions[i]->start, gtf_regions[i]->start, chr_prev, strand_prev);
					set_chr_num(reg, header);
					reg->stop = std::min(gtf_regions[i]->start+c.max_junk, (int) header->target_len[reg->chr_num]);
				}

				if (get_reads)
				{
					printf("get reads from %i bam files for region %s:%i->%i\n", (int)c.bam_files.size(), reg->chr, reg->start, reg->stop);

					// get reads for the large region
					int num_bam = c.bam_files.size();
					if (!c.strand_specific)
						reg->strand = '0';
			
					reg->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);

					printf("sort reads by start position ... \n");
					sort(reg->reads.begin(), reg->reads.end(), CRead::compare_by_start_pos);
					printf("done\n");

					curr = reg->reads.begin();
					last_stop=0;
				}

				// get reads from chromosom region
				// if regions overlapp, decrement the read pointer accordingly
				if (gtf_regions[i]->start<last_stop)
				{
					if (curr == reg->reads.end() && curr != reg->reads.begin())
						curr--;
					while (curr != reg->reads.begin() && (*curr)->start_pos>=gtf_regions[i]->start)
						curr--;
				}
				last_stop = gtf_regions[i]->stop;

				//printf("add reads to gtf_region(%s%c:%i-%i)\n", gtf_regions[i]->chr, gtf_regions[i]->strand, gtf_regions[i]->start, gtf_regions[i]->stop);
				int num_reads = 0;
				while (curr != reg->reads.end())
				{
					if ((*curr)->start_pos>=gtf_regions[i]->stop)
						break;

					if ((*curr)->start_pos>=gtf_regions[i]->start && (*curr)->get_last_position()<gtf_regions[i]->stop)
					{
						gtf_regions[i]->reads.push_back(*curr);
						num_reads++;
					}
					curr++;
				}
				//printf("added %i reads\n", num_reads);
			}
			else
			{
				gtf_regions[i]->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);
			}

			process_gtf_region(gtf_regions[i]);
			if (c.reads_by_chr)
				gtf_regions[i]->reads.clear();
			else
				gtf_regions[i]->clear_reads();
			delete[] gtf_regions[i]->coverage;

		}
		printf("done\n");
		delete reg;
	}
	int num_gtf = gtf_regions.size();
	// add gtf regions
	for (int i=0; i<gtf_regions.size(); i++)
	{
		regions.push_back(gtf_regions[i]);
	}

	if (true)
	{
		printf("merging %i flat file regions with %i gtf regions (discard overlapping regions)\n", num_reg, num_gtf);
		bool change = true;
		int iter = 0;
		while (change)
		{
			printf("merge iteration: %i size:%i\n", iter++, (int) regions.size());
			if (iter >10)
				break;
			change = false;	
			vector<vector<int> > ov_list = region_overlap(regions, regions);
			for (int i=0; i<regions.size(); i++)
			{
				for (int j=0; j<ov_list[i].size(); j++)
				{
					// self overlap
					if (ov_list[i][j]==i)
						continue;

					if (strcmp(regions[i]->chr, regions[ov_list[i][j]]->chr)!=0 || regions[i]->strand != regions[ov_list[i][j]]->strand)
						continue;

					if (regions[ov_list[i][j]]->start==-1 || regions[i]->start==-1)
					{
						continue;
					}
					if (regions[ov_list[i][j]]->transcripts.size()>0 && regions[i]->transcripts.size()>0)
					{
						continue;
					}

					change = true;

					if (regions[ov_list[i][j]]->transcripts.size()==0 && regions[i]->transcripts.size()>0)
					{
						regions[ov_list[i][j]]->start = -1;
						continue;
					}
					if (regions[ov_list[i][j]]->transcripts.size()>0 && regions[i]->transcripts.size()==0)
					{
						regions[i]->start = -1;
						continue;
					}
					//printf("reg: %s%c:%i->%i\n", regions[i]->chr, regions[i]->strand, regions[i]->start, regions[i]->stop);
					//printf("ov: %s%c:%i->%i\n", regions[ov_list[i][j]]->chr, regions[ov_list[i][j]]->strand, regions[ov_list[i][j]]->start, regions[ov_list[i][j]]->stop);
					// merge regions
					regions[i]->start = std::min(regions[i]->start, regions[ov_list[i][j]]->start);
					regions[i]->stop = std::max(regions[i]->stop, regions[ov_list[i][j]]->stop);
					for (int k=0; k<regions[ov_list[i][j]]->transcripts.size(); k++)
						regions[i]->transcripts.push_back(regions[ov_list[i][j]]->transcripts[k]);

					regions[ov_list[i][j]]->start = -1;

				}
			}
			// remove merged regions
			vector<Region*> tmp;
			for (int i=0; i<regions.size(); i++)
				if (regions[i]->start>-1)
					tmp.push_back(regions[i]);
				else
					delete regions[i];

			regions = tmp;
		}
	}
	else if (num_reg>0 && num_gtf>0)
	{
		printf("merging %i flat file regions with %i gtf regions\n", num_reg, num_gtf);
		bool change = true;
		int iter = 0;
		while (change)
		{
			printf("merge iteration: %i size:%i\n", iter++, (int) regions.size());
			if (iter >10)
				break;
			change = false;	
			vector<vector<int> > ov_list = region_overlap(regions, regions);
			for (int i=0; i<regions.size(); i++)
			{
				for (int j=0; j<ov_list[i].size(); j++)
				{
					// self overlap
					if (ov_list[i][j]==i)
						continue;

					if (strcmp(regions[i]->chr, regions[ov_list[i][j]]->chr)!=0 || regions[i]->strand != regions[ov_list[i][j]]->strand)
						continue;

					if (regions[ov_list[i][j]]->start==-1 || regions[i]->start==-1)
					{
						continue;
					}
					change = true;

					//printf("reg: %s%c:%i->%i\n", regions[i]->chr, regions[i]->strand, regions[i]->start, regions[i]->stop);
					//printf("ov: %s%c:%i->%i\n", regions[ov_list[i][j]]->chr, regions[ov_list[i][j]]->strand, regions[ov_list[i][j]]->start, regions[ov_list[i][j]]->stop);
					// merge regions
					regions[i]->start = std::min(regions[i]->start, regions[ov_list[i][j]]->start);
					regions[i]->stop = std::max(regions[i]->stop, regions[ov_list[i][j]]->stop);
					for (int k=0; k<regions[ov_list[i][j]]->transcripts.size(); k++)
						regions[i]->transcripts.push_back(regions[ov_list[i][j]]->transcripts[k]);

					regions[ov_list[i][j]]->start = -1;

				}
			}
			// remove merged regions
			vector<Region*> tmp;
			for (int i=0; i<regions.size(); i++)
				if (regions[i]->start>-1)
					tmp.push_back(regions[i]);

			regions = tmp;
		}
	}
	if (c.reads_by_chr)
	{
		for (int i=0; i<regions.size(); i++)
		{
			set_chr_num(regions[i], header);
		}
		sort(regions.begin(), regions.end(), compare_chr_and_strand);
	}


	vector<int> bias_vector(100, 0);
	//filter regions
	int cnt = 0;
	//int last_stop=0;
	//char* chr_prev = (char*) "xxx";
	last_stop=0;
	chr_prev = (char*) "xxx";
	reg = NULL;

	for (int i=0; i<regions.size(); i++)
	{
		if (c.reads_by_chr)// more fast if many regions are considered
		{
			bool get_reads = false;
			if (strcmp(regions[i]->chr, chr_prev)!=0 || (regions[i]->strand != strand_prev && c.strand_specific))
			{
				if (c.split_chr)
				{
					char fn_out[1000]; 
					sprintf(fn_out, "%s_%s%c", c.fn_out, regions[i]->chr, regions[i]->strand);
					if (fexist(fn_out))
						continue;
					ofs->close();
					ofs->open(fn_out, std::ios::binary);
				}
				get_reads = true;
				chr_prev = regions[i]->chr;
				strand_prev = regions[i]->strand;
				printf("starting with chr: %s%c\n", chr_prev, strand_prev);
				delete reg;
				//reg = new Region(regions[i]);
				//reg->stop = std::min(regions[i]->start+c.max_junk, (int) header->target_len[reg->chr_num]);
				//reg->start = regions[i]->start;
				reg = new Region(regions[i]->start, regions[i]->start, chr_prev, strand_prev);
				set_chr_num(reg, header);
				reg->stop = std::min(regions[i]->start+c.max_junk, (int) header->target_len[reg->chr_num]);
			}

			if (reg->stop<regions[i]->stop)
			{
				int prev_st = reg->stop;
				delete reg;
				//reg = new Region(regions[i]);
				//reg->start = regions[i]->start;
				reg = new Region(regions[i]->start, regions[i]->stop, chr_prev, strand_prev);
				set_chr_num(reg, header);
				reg->stop = std::max(prev_st+c.max_junk, regions[i]->stop);
				reg->stop = std::min(reg->stop, (int) header->target_len[reg->chr_num]);

				get_reads = true;
			}

			if (get_reads)
			{
				printf("get reads from %i bam files for region %s:%i->%i\n", (int)c.bam_files.size(), reg->chr, reg->start, reg->stop);

				for (int i=0; i<(int)c.bam_files.size(); i++)
					printf("bf: %s", c.bam_files[i]);
				// get reads for the large region
				int num_bam = c.bam_files.size();
				//int num_bam = 1;
				if (!c.strand_specific)
					reg->strand = '0';
			
				reg->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);

				// sort reads
				//vector<CRead*>::iterator it;
				//for (int j=0; j<reg->reads.size(); j++)
				//int count = 0;
				//for (it = reg->reads.begin(); it != reg->reads.end(); it++)
				//{
				//	//printf("reg->reads[%i].start_pos: %i (%i)\n", j, reg->reads[j]->start_pos, (int) reg->reads.size());
				//	printf("[%i] it->start_pos: %i (%i)\n", count++, (**it).start_pos, (int) reg->reads.size());
				//}

				printf("sort reads by start position ... \n");
				sort(reg->reads.begin(), reg->reads.end(), CRead::compare_by_start_pos);
				printf("done\n");

				curr = reg->reads.begin();
				last_stop=0;
			}

			// get reads from chromosom region
			// if regions overlapp, decrement the read pointer accordingly
			if (regions[i]->start<last_stop)
			{
				if (curr == reg->reads.end() && curr != reg->reads.begin())
					curr--;
				while (curr != reg->reads.begin() && (*curr)->start_pos>=regions[i]->start)
					curr--;
			}
			last_stop = regions[i]->stop;

			printf("add reads to region(%s%c:%i-%i)\n", regions[i]->chr, regions[i]->strand, regions[i]->start, regions[i]->stop);
			int num_reads = 0;
			while (curr != reg->reads.end())
			{
				if ((*curr)->start_pos>=regions[i]->stop)
					break;

				if ((*curr)->start_pos>=regions[i]->start && (*curr)->get_last_position()<regions[i]->stop)
				{
					regions[i]->reads.push_back(*curr);
					num_reads++;
				}
				curr++;
			}
			printf("added %i reads\n", num_reads);
		}
		else
		{
			int num_bam = c.bam_files.size();
			if (!c.strand_specific)
				regions[i]->strand = '0';
			
			regions[i]->fd_out = fd_null;
			regions[i]->get_reads(&c.bam_files[0], num_bam, c.intron_len_filter, c.filter_mismatch, c.exon_len_filter, c.mm_filter);
		}
		//printf("%i reads in region: %s:%i-%i\n", (int)regions[i]->reads.size(), regions[i]->chr, regions[i]->start, regions[i]->stop);
		printf("\r region %i(%i)", i, (int)regions.size());

		int num_reads=0;
		for (int r=0; r<regions[i]->reads.size(); r++)
			if (regions[i]->reads[r]->mismatches<1)
				num_reads++;
		
		bool is_annotated = regions[i]->transcripts.size()>0;
		if (num_reads<=c.region_filter && !is_annotated)
		{
			if (c.reads_by_chr)
			{
				regions[i]->reads.clear();
			}
			else
				regions[i]->clear_reads();
			continue;
		}
		cnt++;	
		//regions[i]->fd_out = stdout;
		regions[i]->fd_out = fd_null;
		regions[i]->generate_segment_graph(c.seg_filter, c.tss_pval);
		//regions[i]->add_bias_counts(&bias_vector);

		// write region in binary file
		regions[i]->write_binary(ofs);

		// cleanup
		if (c.reads_by_chr)
			regions[i]->reads.clear();// dont delete the reads
		else
			regions[i]->clear_reads();

		//if (cnt>200)
		//	break;
	}
	printf("read %i regions from file\n", (int)regions.size());
	printf("%i regions after filtering\n", cnt);

	// cleanup
	for (int i=0; i<regions.size(); i++)
		delete regions[i];
	
	fclose(fd_null);
	delete reg;
	//fclose(fd_out);
	ofs->close();
	delete ofs;
	bam_close(fd);
	bam_header_destroy(header);

	return 0;
}
