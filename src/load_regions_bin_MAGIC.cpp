#include <stdio.h>
#include <mex.h>
#include <assert.h>
#include <vector>
  using std::vector;
#include "region.h"
#include "get_var.h"
#include <algorithm>
	using std::max;
#include <fstream>
#include "tools.h"
#include <string>
#include <map>
	using std::map;
#include "gtf_tools.h"
#include "tools.h"


// parse the file mapping strain names to lane ids
map<string, vector<string> > parse_mapping()
{
	const char fn[] = "/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/library_ids";

	FILE* fd = fopen(fn, "r");

	map<string, vector<string> > mapping;

	while (~feof(fd))
	{
		char line[1000];
		char* ret = fgets(line, 1000, fd);
		if (!ret)
			break;
		vector<char*> fields = separate(line, '\t');

		string lane_id = fields[0];
		string strain_id = fields[1]+7; // skip the "Sample:" in "Sample:MAGIC176"

		//printf("%s, %s\n", lane_id.c_str(), strain_id.c_str());

		map<string, vector<string> >::iterator it; 
		it = mapping.find(strain_id);
		if (it !=mapping.end())
		{
			mapping[strain_id].push_back(lane_id);
		}
		else
		{
			vector<string> tmp;
			tmp.push_back(lane_id);
			mapping[strain_id] = tmp;
		}

	}
	fclose(fd);
	return mapping;
}

map<string, vector<vector<string> > > parse_clusters()
{
	map<string, vector<string> > mapping = parse_mapping();

	//map<string, vector<string> >::iterator it;
	//for (it=mapping.begin(); it!=mapping.end(); it++)
	//{
	//	printf("%s\t", it->first.c_str());
	//	for (uint i=0; i<it->second.size(); i++)
	//	{
	//		printf("%s\t", it->second[i].c_str());
	//	}
	//	printf("\n");
	//}

	const char bam_dir[] = "/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/alignments_magic/bam_mmr/"; 
	const char fn[] = "/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/classify_bam_out.list";

	FILE* fd = fopen(fn, "r");
	map<string, vector<vector<string> > > clusters;

	while (~feof(fd))
	{
		char line[10000];
		char* ret = fgets(line, 10000, fd);
		if (!ret)
			break;

		vector<char*> fields = separate(line, '\t');

		string gene_name = fields[0];

		//printf("%s, %s\n", lane_id.c_str(), strain_id.c_str());

		map<string, vector<vector<string> > >::iterator it; 
		it = clusters.find(gene_name);
		if (it !=clusters.end())
		{
			printf("Error: gene names should be unique: %s\n", gene_name.c_str());
			printf("Error: line %s\n", line);
			exit(-1);
		}
		else
		{
			vector<vector<string> > samples;
			for (uint i=1; i<fields.size(); i++)
			{
				//printf("sample: %u\t", i);
				vector<char*> strains = separate(fields[i], ',');
				vector<string> tmp;
				for (uint j=0; j<strains.size(); j++)
				{
					vector<char*> ff = separate(strains[j], '.');//get rid of the dot in "MAGIC.146"
					string strain;
					if (ff.size()!=2)
						strain = string(strains[j]);
					else
						strain = string(ff[0])+string(ff[1]);

					vector<string> lanes = mapping[strain];
					//assert(lanes.size()>0);
					//tmp.push_back(strain);

					//printf("%s (%lu):", strain.c_str(), lanes.size());
					for (int i=0; i<lanes.size(); i++)
					{
						char fn_bam[1000];
						sprintf(fn_bam, "%sall_merged.%s.sorted.mmr.bam", bam_dir, lanes[i].c_str());
						tmp.push_back(fn_bam);
					}
					//printf("\t");
				}
				//printf("\n");
				samples.push_back(tmp);
			}
			clusters[gene_name] = samples;
		}
	}
	fclose(fd);
	return clusters;
}

char* get_gene_name(vector<Region*> gtf_regions, Region* reg)
{
	vector<Region*> vreg;
	vreg.push_back(reg);

	vector<vector<int> > ov_list = region_overlap(vreg, gtf_regions);

	assert(ov_list.size() == 1);

	char* name = NULL; 
	for (uint j=0; j<ov_list[0].size(); j++)
	{
		if (strcmp(reg->chr, gtf_regions[ov_list[0][j]]->chr)!=0)
			continue;
		if (reg->strand != gtf_regions[ov_list[0][j]]->strand)
			continue;
		//if (name)
		//{
		//	mexErrMsgTxt("Found multiple overlapping genes");
		//}
		vector<char*> ff = separate((char*) gtf_regions[ov_list[0][j]]->transcript_names[0].c_str(), '.');
		name = ff[0];
		printf("%i, found gene name: %s\n", j, name);
	}
	assert(name);
	printf("found gene name: %s\n", name);
	return name;
}

char* get_gene_name_from_transcript_name(char* trans_name)
{
	vector<char*> ff = separate(trans_name, '.');
	char* gene_name = strstr(ff[0], "AT");
	printf("gene name: %s\n", gene_name);
	return gene_name;
}

/*
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs<3)
	{
		mexErrMsgTxt("Usage: load_regions_bin(fn_graph, gene_idx1, gene_idx2)");	
	}
	
	char* fname = get_string(prhs[0]);
	int num_genes = get_int(prhs[1]);
	int skip = get_int(prhs[2]);
	printf("load %i graphs (skip first %i)\n", num_genes, skip);

	// parse clusters file: obtain the list of bam files grouped into 
	// 20 clusters. The first 19 clusters correspond to those strains that 
	// share the same genotype as one of the 19 founder strains
	// the last cluster contains bam files from strains that had a unresolved genotype
	map<string, vector<vector<string> > > gene_samples = parse_clusters();

	// parse the original genome annotation to assign gene names to graphs
	//char fn_gtf[] = "/fml/ag-raetsch/nobackup2/projects/sequencing_runs/A_thaliana_magic/results/annotation/consolidated_annotation.allmerged.Col_0.gtf"; 

	//vector<Region*> gtf_genes = parse_gtf(fn_gtf);

	static std::ifstream ifs;
	if (num_genes<0 && ifs && ifs.is_open())
	{
		printf("load_regions_bin: closing binary file\n"); 
		ifs.close();
		return;
	}

	if (!ifs.is_open())
	{
		printf("load_regions_bin: opening binary file\n"); 
		ifs.open(fname, std::ios::binary);
	}
	FILE* fd_null = fopen("/dev/null", "w");

	vector<Region*> regions; 
	int cnt = 0;
	while (ifs.good() && cnt<num_genes+skip)
	{
		if ((cnt++)%100==0)
			printf("\r %i", cnt);
		Region* r = new Region();
		int ret = r->read_binary(&ifs);

		if (cnt<=skip)
		{
			delete r;
			continue;
		}
		if (!ifs.eof() && ret==0)
			regions.push_back(r);
		else
		{
			printf("read from file %s, %i, read %lu regions\n", fname, ret, regions.size());
			delete r;
			break;
		}
	}


	//ifs.close();
	//r->write_segment_graph(stdout);

	const char *field_names[] = {"chr", "strand", "start", "stop", "segments", "coverage", "seg_admat", "transcripts", "transcript_names", "pair_mat"};
	int nof_fields = 10;

	int num = regions.size();
	printf("\nparsed %i regions from file\n", num);

	mwSize dims[2] = {1, num};
	plhs[0] = mxCreateStructArray(2, dims, nof_fields, field_names);

	// get index of fields
	int chr_field = mxGetFieldNumber(plhs[0], "chr");
	int strand_field = mxGetFieldNumber(plhs[0], "strand");
	int start_field = mxGetFieldNumber(plhs[0], "start");
	int stop_field = mxGetFieldNumber(plhs[0], "stop");
	int segments_field = mxGetFieldNumber(plhs[0], "segments");
	int coverage_field = mxGetFieldNumber(plhs[0], "coverage");
	int seg_admat_field = mxGetFieldNumber(plhs[0], "seg_admat");
	int transcripts_field = mxGetFieldNumber(plhs[0], "transcripts");
	int transcript_names_field = mxGetFieldNumber(plhs[0], "transcript_names");
	int pair_field = mxGetFieldNumber(plhs[0], "pair_mat");


	cnt = 0;
	for (int j=0; j<num; j++)
	{

		// obtain the gene name
		//char* gene_name = get_gene_name(gtf_genes, regions[j]);
		char* gene_name = get_gene_name_from_transcript_name((char*) regions[j]->transcript_names[0].c_str());

		//assign the correct sample cluster
		//vector<vector<string> > samples = gene_samples.begin()->second;
		vector<vector<string> > samples = gene_samples[string(gene_name)];

		printf("samples: %lu\n", samples.size());
		//for (uint i=0; i<samples.size(); i++)
		//{
		//	printf("\tsample %i:\n", i);
		//	for (int k=0; k<samples[i].size(); k++)
		//	{
		//		printf("\t\t%s:\n", samples[i][k].c_str());
		//	}
		//	printf("\n");
		//}

		if ((cnt++)%100==0)
			printf("\r %i\t(%i)", cnt, num);

		// chromosome
		mxSetFieldByNumber(plhs[0], j, chr_field, mxCreateString(regions[j]->chr));

		// strand
		char strand[2];
		strand[0] = regions[j]->strand;
		strand[1] = '\0';
		mxSetFieldByNumber(plhs[0], j, strand_field, mxCreateString(strand));

		// start 
		mxArray* start = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(start) = regions[j]->start;
		mxSetFieldByNumber(plhs[0], j, start_field, start);

		// stop
		mxArray* stop = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(stop) = regions[j]->stop;
		mxSetFieldByNumber(plhs[0], j, stop_field, stop);

		// segments
		int num_seg = regions[j]->segments.size();
		mxArray* segments = mxCreateDoubleMatrix(2, num_seg, mxREAL);
		double* seg = mxGetPr(segments);
		for (int i=0; i<num_seg; i++)
		{
			seg[2*i] = (double) regions[j]->segments[i].first;
			seg[2*i+1] = (double) regions[j]->segments[i].second;
		}
		mxSetFieldByNumber(plhs[0], j, segments_field, segments);

		// coverage
		assert(regions[j]->seg_cov.size()==num_seg);
		mxArray* coverage = mxCreateDoubleMatrix(samples.size(), num_seg, mxREAL);
		double* cov = mxGetPr(coverage);
		mxSetFieldByNumber(plhs[0], j, coverage_field, coverage);

		// segment adjacency matrix
		int seg_admat_size = regions[j]->admat.size();
		//printf("admat.size(): %i, %i\n", len, seg_admat_field);
		//mxArray* seg_admat = mxCreateDoubleMatrix(len, len, mxREAL);
		int dims[3];
		dims[0] = seg_admat_size;
		dims[1] = seg_admat_size;
		dims[2] = samples.size();
		mxArray* seg_admat = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		double* seg_mat = mxGetPr(seg_admat);
		mxSetFieldByNumber(plhs[0], j, seg_admat_field, seg_admat);

		//transcripts
		int num_trans = regions[j]->transcript_paths.size();
		mxArray* trans_mat = mxCreateDoubleMatrix(num_trans, num_seg, mxREAL);
		double* trm = mxGetPr(trans_mat);
		for (int i=0; i<num_trans; i++)
		{
			for (int k=0; k<regions[j]->transcript_paths[i].size(); k++)
			{
				int p = regions[j]->transcript_paths[i][k];
				if (p>=num_seg || p<0)
				{
					printf("Error: trans: %i, segment: %i num_seg: %i\n", i, p, num_seg);
					continue;
				}
				trm[p*num_trans+i] = 1;
			}
		}
		mxSetFieldByNumber(plhs[0], j, transcripts_field, trans_mat);
		
		//transcript names
		mxArray* names_cell = mxCreateCellArray(1, &num_trans);
		for (int i=0; i<regions[j]->transcript_names.size(); i++)
		{
			const char* name = regions[j]->transcript_names[i].c_str();
			mxArray* sname = mxCreateString(name);
			mxSetCell(names_cell, i, sname);
		}
		mxSetFieldByNumber(plhs[0], j, transcript_names_field, names_cell);

		// segment adjacency matrix
		int len = regions[j]->pair_mat.size();
		dims[0] = len;
		dims[1] = len;
		dims[2] = samples.size();
		mxArray* pair_mat = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		//printf("admat.size(): %i, %i\n", len, seg_admat_field);
		//mxArray* pair_mat = mxCreateDoubleMatrix(len, len, mxREAL);
		mxSetFieldByNumber(plhs[0], j, pair_field, pair_mat);


		for (int k = 0; k<samples.size(); k++)
		{
			// each string in bam_files may be a comma separated list of files
			vector<char*> bams;
			for (int i=0; i<samples[k].size(); i++)
				bams.push_back((char*) samples[k][i].c_str());

			int intron_len_filter = 200000;
			int filter_mismatch = 10;
			int exon_len_filter = 0;
			bool mm_filter = false;
			//printf("chr: %s\n", regions[j]->chr);
			regions[j]->fd_out = fd_null;
			regions[j]->clear_reads();
			regions[j]->get_reads(&bams[0], bams.size(), intron_len_filter, filter_mismatch, exon_len_filter, mm_filter);
			//printf("admat.size():%i  segments.size():%i\n", (int) regions[j]->segments.size(), (int) regions[j]->admat.size());
			regions[j]->update_coverage_information();
			//printf("num_reads: %i\n", (int) regions[j]->reads.size());
			regions[j]->compute_coverage();
			regions[j]->compute_seg_cov();
			regions[j]->compute_pair_mat();
			//printf("admat.size():%i  segments.size():%i\n", (int) regions[j]->segments.size(), (int) regions[j]->admat.size());

			for (int i=0; i<num_seg; i++)
				cov[i*samples.size()+k] = regions[j]->seg_cov[i];

			for (int i=0; i<seg_admat_size; i++)
			{
				for (int l=0; l<seg_admat_size; l++)
				{
					seg_mat[k*seg_admat_size*seg_admat_size+seg_admat_size*i+l] = regions[j]->admat[i][l];
				}
			}

			double* pmat = mxGetPr(pair_mat);
			for (int i=0; i<len; i++)
			{
				for (int l=0; l<len; l++)
				{
					pmat[k*len*len+len*i+l] = regions[j]->pair_mat[i][l];
				}
			}

		}
	}

	// cleanup
	for (int j=0; j<num; j++)
		delete regions[j];
}
