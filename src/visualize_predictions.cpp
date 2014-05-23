#include "bam_region.h"
#include "gtf_tools.h"
#include <numeric>
#include <math.h>

int* get_bam_cnt(vector<char*>* bams)
{
	int* nums = new int[bams->size()];
	for (int j=0; j<bams->size(); j++)
	{
		char fname[1000];
		sprintf(fname, "%s.cnt", bams->at(j));
		FILE* fd = fopen(fname, "r"); 
		if (fd)
		{
			char line[1000]; 
			if (fgets(line, 1000, fd)==NULL) 
			{
				printf("failed to read file: %s\n", fname); 
				nums[j] = 1e6; 
				fclose(fd);
				continue; 
			}
			nums[j] = atof(line); 
			//printf("cnt[%i]: %i\n", j, nums[j]);
			fclose(fd);
		}
		else
		{
			printf("did not find file: %s\n", fname); 
			nums[j] = 1e6;
		}
	}
	return nums; 
}

int get_coverage(uint32_t**& cov, Region* reg, vector<char*>* bams)
{
	int num_pos = reg->stop-reg->start+1;
	cov = new uint32_t*[bams->size()];

	for (int k=0; k<bams->size(); k++)
	{
		Bam_Region* region = new Bam_Region(reg);
		region->get_reads(&bams->at(k), 1, 200000, 3, 8);
		region->compute_coverage();

		cov[k] = new uint32_t[num_pos];
		memcpy(cov[k], region->coverage, num_pos*sizeof(uint32_t)); 

		delete region;
	}
	return (int) bams->size(); 
}

int get_num_reads(Region* reg, char* fn_bam)
{
	int num_reads = -1; 
	Bam_Region* region = new Bam_Region(reg);
	region->get_reads(&fn_bam, 1, 200000, 3, 8);
	//num_reads = region->all_reads.size();  
	num_reads = region->reads.size(); // filtered
	delete region;
	return num_reads; 
}

void create_blocks(Region* reg, vector<pair<int, int> >* blocks)
{

	int num_pos = reg->stop-reg->start+1;
	bool* map = new bool[num_pos];
	assert(map);
	for (int i=0; i<num_pos; i++)
		map[i]=false;

	for (int i=0; i<reg->transcripts.size(); i++)
	{
		for (int j=0; j<reg->transcripts[i].size(); j++)
		{
			for (int k=reg->transcripts[i][j].first-100; k<=reg->transcripts[i][j].second+100; k++)
			{
				if (k-reg->start>=0 && k-reg->start<num_pos)
				map[k-reg->start] = true;
			}
		}
	}

	int start = 0;
	int stop = 0;
	for (int i=0; i<num_pos-1; i++)
	{
		if (map[i] && !map[i+1])
		{
			stop=i;
			pair<int, int> p(start, stop); 
			blocks->push_back(p);
		}
		else if (!map[i] && map[i+1])
		{
			start = i+1;
		}
	}
	if (map[num_pos-1])
	{
		stop=num_pos-1;
		pair<int, int> p(start, stop); 
		blocks->push_back(p);
	}

	//for (int i=0; i<blocks->size(); i++)
	//{
	//	printf("%i->%i\n", blocks->at(i).first, blocks->at(i).second);
	//}
}

int main(int argc, char** args)
{
	if (argc<4)
	{
		printf("Usage: %s <path_to_gnuplot> <out_dir> <fn_gtf_pred> <fn_gtf_anno> <fn_bam1> [fn_bam2 [...]]\n", args[0]);
		exit(-1);
	}

	char* path_to_gnuplot = args[1];
	char* out_dir = args[2]; 
	char* fn_gtf_pred = args[3]; 
	char* fn_gtf_anno = args[4]; 


	char fn_tmp[1000]; 
	char fn_gnu[1000]; 
	char fn_ps[1000]; 
	sprintf(fn_tmp, "%s/data.tmp", out_dir); 
	sprintf(fn_gnu, "%s/gnu.tmp", out_dir); 
	sprintf(fn_ps, "%s/ps.tmp", out_dir); 

	vector<char*> bam_files; 
	for (int i=5; i<argc; i++)
	{
		bam_files.push_back(args[i]); 
	}

	// parse gtfs 
	////////////////////////////////////////////////////////////////////////////////
	vector<Region*> pred_genes = parse_gtf(fn_gtf_pred);
	printf("number of predicted genes: %lu\n", pred_genes.size()); 
	vector<Region*> anno_genes = parse_gtf(fn_gtf_anno);
	printf("number of annotated genes: %lu\n", anno_genes.size()); 

	for (int i=0; i<pred_genes.size(); i++)
	{
		for (int j=0; j<pred_genes[i]->transcript_names.size(); j++)
		{
			pred_genes[i]->transcript_names[j] = string("pred_") + pred_genes[i]->transcript_names[j]; 
		}
	}
	for (int i=0; i<anno_genes.size(); i++)
	{
		for (int j=0; j<anno_genes[i]->transcript_names.size(); j++)
		{
			anno_genes[i]->transcript_names[j] = string("anno_") + anno_genes[i]->transcript_names[j]; 
		}
	}

	pred_genes.insert(pred_genes.end(), anno_genes.begin(), anno_genes.end()); 

	vector<Region*> regions = merge_overlapping_regions(pred_genes);

	

	for (int j=0; j<regions.size(); j++)
	{
		char fn_pdf[1000]; 
		sprintf(fn_pdf, "%s/region%i.pdf", out_dir, j); 

		bool skip = true; 
		for (int i=0; i<regions[j]->transcript_names.size(); i++)
		{
			if (strstr(regions[j]->transcript_names[i].c_str(), "pred"))
				skip = false; 
		}
		if (skip)
			continue; 

		if (regions[j]->start>1000)
			regions[j]->start -= 1000; 
		else
			regions[j]->start = 1;

		regions[j]->stop += 1000;

		uint32_t** all_cov;
		int num_bams = get_coverage(all_cov, regions[j], &bam_files);

		printf("got reads from %i bam files\n", num_bams); 

		// create blocks around exons to plot
		vector<pair<int, int> > blocks; 
		create_blocks(regions[j], &blocks);


		vector<int> cum_len; 
		int len = 0;
		for (int k=0; k<blocks.size(); k++)
		{
			len += blocks[k].second-blocks[k].first+1;
			cum_len.push_back(len);
		}

		int num_pos = regions[j]->stop-regions[j]->start+1;
		float max_val = 0;
		float max_exome = 0;
		int b = 0;
		FILE* fd_out = fopen(fn_tmp, "w");
		assert(fd_out); 
		for (int k=0; k<num_pos; k++)
		{
			while (k>blocks[b].second && b<blocks.size()-1)
				b++;

			if (k<blocks[b].first || k>blocks[b].second)
				continue; 
			
			int block_pos = cum_len[b]-(blocks[b].second-k);
			fprintf(fd_out, "%i", block_pos);
			for (int l=0; l<num_bams; l++)
			{
				if (log(all_cov[l][k])>max_val)
					max_val = log(all_cov[l][k]); 

				if (all_cov[l][k]>0)
					fprintf(fd_out, "\t%.3f", log(all_cov[l][k]));
				else
					fprintf(fd_out, "\t%.3f", 0.0);
			}
			fprintf(fd_out, "\n");
		}
		fclose(fd_out);


		// write gnu-plot script
		FILE* fd_plot = fopen(fn_gnu, "w");
		
		//fprintf(fd_plot, "set multiplot layout 2,1 rowsfirst\n");
		fprintf(fd_plot, "set terminal postscript color\n"); 
		fprintf(fd_plot, "set output \"%s\"\n", fn_ps);
		//fprintf(fd_plot, "set multiplot\n");
		fprintf(fd_plot, "set title \'%s:%i-%i (%c)\'\n", regions[j]->chr, regions[j]->start, regions[j]->stop, regions[j]->strand);
		fprintf(fd_plot, "unset key\n");// no legend
		//fprintf(fd_plot, "set size 1,0.5\n");
		//fprintf(fd_plot, "set origin 0, 0.5\n");
		//fprintf(fd_plot, "set yrange [] writeback\n");// use the same range for both coverage plots
		//fprintf(fd_plot, "set xrange [] writeback\n");// use the same range for both coverage plots
		//fprintf(fd_plot, "plot ");
		//for (int k=0; k<num_bams; k++)
		//{
		//	fprintf(fd_plot, "\"%s\"  using 1:%i with lines t \"col %i\", ", fn_tmp, k+2, k+2);
		//}
		//fprintf(fd_plot, "\n");

		// plot annotation
		//fprintf(fd_plot, "set size 1,0.99\n");
		//fprintf(fd_plot, "set origin 0, 0\n");
		//fprintf(fd_plot, "set xrange restore\n");
		fprintf(fd_plot, "set yrange [-%lu:%.2f]\n", regions[j]->transcripts.size()+2, max_val+0.1);

		int cnt=0; 
		for (int k=0; k<regions[j]->transcripts.size(); k++)
		{
			cnt++;
			int y = -k-1;
			fprintf(fd_plot, "set label %i \"%s\" at %i,%i\n", cnt, regions[j]->transcript_names[k].c_str(), 0, y);
			int b=0;
			for (int e=0; e<regions[j]->transcripts[k].size(); e++)
			{
				int exon_start = regions[j]->transcripts[k][e].first-regions[j]->start;
				int exon_stop = regions[j]->transcripts[k][e].second-regions[j]->start;
				int flag = regions[j]->transcripts[k][e].flag;
				while (exon_start>blocks[b].second &&  b<blocks.size()-1)
					b++;

				assert(exon_start>=blocks[b].first && exon_start<=blocks[b].second);

				int bstart = cum_len[b]-(blocks[b].second-exon_start);
				int bstop = cum_len[b]-(blocks[b].second-exon_stop);
				if (flag==4)
					fprintf(fd_plot, "set arrow from %i,%i to %i,%i nohead lc rgb \'red\'\n", bstart, y, bstop, y);
				else
					fprintf(fd_plot, "set arrow from %i,%i to %i,%i nohead lc rgb \'blue\'\n", bstart, y, bstop, y);
			}
		}
		fprintf(fd_plot, "plot ");
		for (int k=0; k<num_bams; k++)
		{
			fprintf(fd_plot, "\"%s\"  using 1:%i with lines t \"col %i\", ", fn_tmp, k+2, k+2);
		}
		fprintf(fd_plot, "\n");

		//fprintf(fd_plot, "plot NaN\n");
		//fprintf(fd_plot, "unset multiplot\n");

		fclose(fd_plot);


		char command[1000]; 
		sprintf(command, "cat %s | %s", fn_gnu, path_to_gnuplot); 
		//printf("run: %s\n", command); 
		system(command); 

		sprintf(command, "ps2pdf %s %s", fn_ps, fn_pdf); 
		//printf("run: %s\n", command); 
		system(command); 
		sprintf(command, "rm %s", fn_ps); 
		//printf("run: %s\n", command); 
		system(command); 
		printf("file %s written\n", fn_pdf); 

	}



	return 0;
}
