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
#include "graph_tools.h"

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
	FILE* fd_out;
	char* fn_graph;
	bool mm_filter;
	int intron_len_filter;
	int filter_mismatch;
	int exon_len_filter;
	int max_num_trans;
	float eta1;
	float eta2;
	int lambda;
	int min_trans_len;
	bool use_pair;
};

template <typename T>
class semi_sparse_3d_matrix{

	private:
		map<pair<int, int>, vector<T> > mat;

	public:

		void add(int i, int j, T val)
		{
			pair<int, int> p(i, j);
			mat[p].push_back(val);
		}

		T get(int i, int j, int k)
		{
			pair<int, int> p(i, j);
			map<char,int>::iterator it2;
			typename map<pair<int, int>, vector<T> >::iterator it; 
			it = mat.find(p);
			if (it==mat.end())
			{
				return 0.0;
			}
			vector<T>* vec = &it->second;
			assert(vec->size()>k);

			return vec->at(k);
		}

		size_t size()
		{
			return mat.size();
		}
};


int parse_args(int argc, char** argv,  Config* c)
{
	if (argc<4)
	{
			fprintf(stdout, "Usage: %s <fn_graph> [options] <fn_bam1,fn_bam2> <fn_bam3> ...\n", argv[0]);
			fprintf(stdout, "(this corresponds to two samples where sample 1 has two bam files and sample two has one)\n");
			fprintf(stdout, "options:\n");
			fprintf(stdout, "\t--num-trans\t(default 5) maximal number of transcripts predicted in addition to the annotated transcripts");
			fprintf(stdout, "\t--param-eta1\t(default 0.0)");
			fprintf(stdout, "\t--param-eta2\t(default 0.42)");
			fprintf(stdout, "\t--param-lambda\t(default 3)");

			fprintf(stdout, "read related options:\n");
			fprintf(stdout, "\t--min-exonic-len\t(default 0) minimal number of aligned positions on each side of an intron\n");
			fprintf(stdout, "\t--mismatches\t\t(default 10) maximal number of mismatches\n");
			fprintf(stdout, "\t--best\t\t\t(flag) filter for best alignment\n");
			fprintf(stdout, "\t\t\n");
			return -1;
	}

	// defaults
	c->fn_graph = argv[1];	
	c->fd_out = stdout;
	c->intron_len_filter = 200000;
	c->filter_mismatch = 10;
	c->exon_len_filter = 0;
	c->mm_filter = false;
	c->max_num_trans = 5;
	c->eta1 = 0.0;
	c->eta2 = 0.42;
	c->lambda = 3;
	c->min_trans_len = 100;
	c-> use_pair = true;

    for (int i = 2; i < argc; i++)  
    {
	   	if (strcmp(argv[i], "--min-exonic-len") == 0)
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
	    else if (strcmp(argv[i], "--max-num-trans") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --max-num-trans\n") ;
                return -1;
            }
            i++;
			c->max_num_trans = atoi(argv[i]);
        }
	    else if (strcmp(argv[i], "--param-eta1") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --param-eta1\n") ;
                return -1;
            }
            i++;
			c->eta1 = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--param-eta2") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --param-eta2\n") ;
                return -1;
            }
            i++;
			c->eta2 = atof(argv[i]);
        }
	    else if (strcmp(argv[i], "--param-lambda") == 0)
        {
            if (i + 1 > argc - 1)
            {
                fprintf(stderr, "ERROR: Argument missing for option --param-lambda\n") ;
                return -1;
            }
            i++;
			c->lambda = atoi(argv[i]);
        }
		else
		{
			c->bam_files.push_back(argv[i]);
		}	
	}
	return 0;
}
int connect_neighbors(vector<vector<vector<float> > >* all_admat, vector<segment>* segments)
{
	int num_samples = all_admat->size();
	for (int j = 1; j<segments->size()-1; j++)
	{
		if (segments->at(j).second+1 == segments->at(j+1).first)
		{
			assert(all_admat->at(0).size()>j+2);
			//printf("segment[%i], %i,%i %.2f\n", j, segments->at(j).first, segments->at(j).second, all_admat->at(0)[j+1][j+2]);
			//printf("segment[%i], %i,%i\n", j, segments->at(j).first, segments->at(j).second);
			//printf("\n");

			for (int s=0; s<num_samples; s++)
			{
				assert(all_admat->at(s).size()>j+2);
				all_admat->at(s)[j+1][j+2] = -1;
			}
		}
	}

}

int union_connections(vector<vector<vector<float> > >* all_admat)
{
	// assuming that the first and the last segment are artificial source/sink nodes
	assert(all_admat->at(0).size()>0);
	assert(all_admat->at(0).size()==all_admat->at(0)[0].size());
	int num_nodes = all_admat->at(0).size();
	int num_samples = all_admat->size();

	for (int j=1; j<num_nodes-1; j++)
	{
		for (int k=1; k<num_nodes-1; k++)
		{
			float max = -2;
			for (int s=0; s<num_samples; s++)
			{
				if (all_admat->at(s)[j][k]>max)
					max = all_admat->at(s)[j][k];
			}
			if (max<0)
				continue;

			for (int s=0; s<num_samples; s++)
			{
				if (all_admat->at(s)[j][k]<0)
					all_admat->at(s)[j][k] = 0;
			}
		}
	}
	return 0;
}

semi_sparse_3d_matrix<float> compute_intron_list(vector<vector<vector<float> > >* all_admat)
{
	semi_sparse_3d_matrix<float> list; 
	int r = all_admat->size();
    int s = all_admat->at(0).size();;
    for (int j=1; j<s-1; j++)
	{
        for (int k=j+1; k<s-1; k++)
		{
            if (all_admat->at(0)[j][k]>=0)
            {
				// valid intron not just a neighboring segment
				for (int sample=0; sample<r; sample++)
				{
					//seg.features.push_back(all_admat->at(sample)[j][k]);
					list.add(j, k, all_admat->at(sample)[j][k]);
				}
			}
		}
	}
	return list;
}

semi_sparse_3d_matrix<float> compute_pair_list(vector<vector<vector<int> > >* all_pair_mat)
{
	semi_sparse_3d_matrix<float> list; 
	int r = all_pair_mat->size();
    int s = all_pair_mat->at(0).size();;
    for (int j=1; j<s-1; j++)
	{
        for (int k=j+1; k<s-1; k++)
		{
            if (all_pair_mat->at(0)[j][k]>=0)
            {
				for (int sample=0; sample<r; sample++)
				{
					list.add(j, k, (float) all_pair_mat->at(sample)[j][k]);
				}
			}
		}
	}
	return list;
}

vector<int> range(int lb, int ub)
{
	vector<int> ret;
	for (int i=lb; i<=ub; i++)
		ret.push_back(i);
	
	return ret;
}

int make_qp(Region* graph, vector<vector<vector<float> > >* all_admat, vector<vector<vector<int> > >* all_pair_mat, vector<vector<float> >* all_seg_cov, const Config* config)
{

	int len = 0;
	for (int i=0; i<graph->segments.size(); i++)
		len += graph->segments[i].second-graph->segments[i].first; 
	
	printf("length of all segments: %i\n", len);
	int min_trans_len = min(len, config->min_trans_len);


	semi_sparse_3d_matrix<float> mat;
	mat.add(35, 57, 5.0);
	printf("val: %.3f\n", mat.get(35, 57, 0));
	printf("val: %.3f\n", mat.get(35, 21, 0));


	semi_sparse_3d_matrix<float> intron_list = compute_intron_list(all_admat);
	semi_sparse_3d_matrix<float> pair_list = compute_pair_list(all_pair_mat);

	printf("intron_list.size(): %lu\n", intron_list.size());
	printf("pair_list.size(): %lu\n", pair_list.size());

	int num_introns = intron_list.size();
	int num_samples = all_admat->size();
	vector<int> binary_var_index_all;
	int n_of_equalities = 0;

	int max_num_trans = config->max_num_trans + graph->transcript_paths.size();

	// compute the number of variables
	//
	// U_st : usage of segment s in transcript t                                     (binary)
	// I_t  : indicator if transcript t has weight>=0 in any sample                  (binary)
	// E_str : expected coverage of segment s in transcript t in sample r                 
	// W_tr  : weight of transcript t in sample r
	// L_sr  : loss of segment s (deviation from the observed coverage) in sample r
	// C_ctr : Spliced read expected coverage (c is a subset of sxs) in sample r
	// D_cr  : deviation of expected intron coverage from observed intron cov in sample r
	//
	// var = [ U_st I_t E_str W_tr L_sr C_ctr D_cr]
	int c = num_introns;
	int t = max_num_trans;
	int s = graph->segments.size();
	int r = num_samples;
	//int nc = num_clusters;
	printf("#introns: %i, max #transcripts: %i(%lu+%i), #segments: %i, #samples: %i\n", c, t, graph->transcript_paths.size(), t-graph->transcript_paths.size(),  s, r);
	

	// compute the number of variables of the optimization problem
	int num_var = 0;
	// indices of variables in the parameter vector
	vector<int> U_idx = range(0, s*t-1); 					num_var += U_idx.size();
	vector<int> I_idx = range(num_var, num_var+t-1);        num_var += I_idx.size();
	vector<int> E_idx = range(num_var, num_var+s*t*r-1);    num_var += E_idx.size();
	vector<int> W_idx = range(num_var, num_var+t*r-1);      num_var += W_idx.size();
	vector<int> L_idx = range(num_var, num_var+s*r*2-1);    num_var += L_idx.size();
	vector<int> C_idx = range(num_var, num_var+c*t*r-1);    num_var += C_idx.size();
	vector<int> D_idx = range(num_var, num_var+c*r*2-1);    num_var += D_idx.size();
	
	//if use_cluster
	//    M_idx = num_var+1:num_var+r*nc; num_var = num_var+length(M_idx); % cluster assignment
	//    m_idx = num_var+1:num_var+t*nc; num_var = num_var+length(m_idx); % cluster centroids
	//    K_idx = num_var+1:num_var+nc; num_var = num_var+length(K_idx);   % cluster usage
	//end
	if (config->use_pair)
	{
	    int np = pair_list.size();
	    vector<int>  P_idx = range(num_var, num_var+t*np-1); num_var += P_idx.size();
	    vector<int> SP_idx = range(num_var, num_var+np-1);   num_var += SP_idx.size(); // slacks for pair observation
	}


	printf("QP has %i variables\n", num_var);

	return 0;
}


int main(int argc, char* argv[])
{
	Config c;
	int ret = parse_args(argc, argv, &c);
	if (ret!=0)
		return ret;

	printf("loading graphs from file: %s\n", c.fn_graph);
	
	static std::ifstream ifs;
	ifs.open(c.fn_graph, std::ios::binary);
	FILE* fd_null = fopen("/dev/null", "w");


	printf("found %i samples\n", (int)c.bam_files.size());
	vector<vector<char*> > samples;
	for (int i=0; i<c.bam_files.size(); i++)
	{
		printf("Sample %i:\n", i+1);
		vector<char*> bams = separate(c.bam_files[i], ',');
		for (int j = 0; j<bams.size(); j++)
			printf("\t%s\n", bams[j]);

		samples.push_back(bams);
	}

	// load graphs from binary file
	vector<Region*> graphs; 
	int cnt = 0;
	int num_return = 3;
	while (ifs.good()&& cnt<num_return)
	{
		if ((cnt++)%100==0)
			printf("\r %i", cnt);
		Region* r = new Region();
		int ret = r->read_binary(&ifs);
		if (!ifs.eof() && ret==0)
			graphs.push_back(r);
		else
		{
			printf("read from file %s, %i, read %lu graphs\n", c.fn_graph, ret, graphs.size());
			delete r;
			break;
		}
	}
	printf("\n");

	printf("obtained %i graphs\n", graphs.size());

	if (false)
	{
		for (uint j=0; j<graphs.size(); j++)
		{
			vector<int> initial;
			vector<int> terminal;
			initial.push_back(1);
			terminal.push_back(graphs[j]->admat.size()-2);
			for (uint k=2; k<graphs[j]->admat.size()-2; k++)
			{
				if (graphs[j]->is_initial(k+1))
					initial.push_back(k);
				if (graphs[j]->is_terminal(k+1))
					terminal.push_back(k);
			}
			int num_paths = 0;
			for (uint i=0; i<initial.size(); i++)
				for (uint k=0; k<terminal.size(); k++)
					num_paths += count_num_paths(graphs[j]->admat, initial[i]+1, terminal[k]+1);
			if (false && num_paths==0)
			{
				for (uint k=0; k<graphs[j]->admat.size(); k++)
				{
					for (uint l=0; l<graphs[j]->admat.size(); l++)
					{
						printf("%.3f, ", graphs[j]->admat[k][l]);
					}
					printf("\n");
				}
			}

			//printf("graph %u: initial:%lu, terminal:%lu, admat.size:%lu,\t%i\tpaths\n", j, initial.size(), terminal.size(), graphs[j]->admat.size(), num_paths);
			printf("xx\t%i\n", num_paths);
		}
	}

	for (int i=0; i<graphs.size(); i++)
		graphs[i]->write_HDF5();
	return 0;

	for (uint j=0; j<graphs.size(); j++)
	{
		// store coverage informaton for all samples
		vector<vector<vector<float> > > all_admat;
		vector<vector<vector<int> > > all_pair_mat;
		vector<vector<float> > all_seg_cov;

		for (uint k = 0; k<samples.size(); k++)
		{
			// each string in bam_files may be a comma separated list of files
			vector<char*> bams = samples[k];

			int intron_len_filter = 200000;
			int filter_mismatch = 10;
			int exon_len_filter = 0;
			bool mm_filter = false;//multi mappers
			//printf("chr: %s\n", graphs[j]->chr);
			graphs[j]->fd_out = fd_null;
			graphs[j]->clear_reads();
			graphs[j]->get_reads(&bams[0], bams.size(), intron_len_filter, filter_mismatch, exon_len_filter, mm_filter);
			printf("admat.size():%i  segments.size():%i\n", (int) graphs[j]->segments.size(), (int) graphs[j]->admat.size());
			graphs[j]->update_coverage_information();
			printf("num_reads: %i\n", (int) graphs[j]->reads.size());
			graphs[j]->compute_coverage();
			graphs[j]->compute_seg_cov();
			graphs[j]->compute_pair_mat();
			printf("admat.size():%i  segments.size():%i\n", (int) graphs[j]->segments.size(), (int) graphs[j]->admat.size());
			all_admat.push_back(graphs[j]->admat);
			all_pair_mat.push_back(graphs[j]->pair_mat);
			all_seg_cov.push_back(graphs[j]->seg_cov);
		}

		// process the connectivity matrix of the graph
		// add implicit connections between neighboring segments
		connect_neighbors(&all_admat, &graphs[j]->segments);

		// make sure connections are valid in each sample if 
		// they have evidence in one sample
		union_connections(&all_admat);

		make_qp(graphs[j], &all_admat, &all_pair_mat, &all_seg_cov, &c); 
	}
}


