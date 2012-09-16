
#ifndef _REGION_H__
#define _REGION_H__

//#define READ_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "genome.h"
#include "read.h"
#include <utility>
	using std::pair;
#include <vector>
	using std::vector;
#include <fstream>

typedef pair<int, int> segment;
#define NO_CONNECTION -2
#define NEIGHBOR -1
#define CONNECTION 0

//#define INNER=0
//#define INITIAL=1
//#define TERMINAL=2

class Region
{
	public:
		int start; 
		int stop;
		int chr_num;
		char* chr;
		char strand; 
		vector<CRead> all_reads; 
		vector<CRead*> reads; 
		bool reads_sorted;
		uint32_t* coverage;
		uint32_t* intron_coverage;
		vector<segment> intron_list;
		vector<segment> unique_introns;
		vector<int> intron_counts;
		vector<segment> segments;
		vector<float> seg_cov;
		vector<vector<float> > admat;
		vector<vector<int> > pair_mat;
		vector<vector<segment> > transcripts;
		vector<vector<int> > transcript_paths;
		FILE* fd_out;

		Genome* gio; 
    	
			
		/** default constructor*/	
		Region();

		/** constructor*/
		Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname);
		Region(int pstart, int pstop, int pchr_num, char pstrand);
		Region(int pstart, int pstop, char* chr, char pstrand);
		Region(Region* reg);

		/** destructor*/
		~Region();

		void clear_reads();

		char* get_sequence()
		{
			if (!seq)
				load_genomic_sequence();
			return seq;
		};

		bool check_region(){return (start>=0 && stop<gio->contig_len(chr_num));};

		float get_coverage_global(int pstart, int pstop);
		int get_intron_conf(int intron_start, int intron_stop);

		void set_gio(Genome* pgio);
	
		void load_genomic_sequence(); 			

		void get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter, bool mm_filter);

		void compute_coverage();

		void compute_intron_coverage();

		void compute_intron_list();

		char* get_region_str();

		void generate_segment_graph(float seg_filter, float tss_pval);

		void compute_admat(vector<int> starts, vector<int> stops);

		int compute_pair_mat();

		vector<int> get_parents(int node);

		vector<int> get_children(int node);

		void print_segments_and_coverage(_IO_FILE*& fd);

		virtual void print(_IO_FILE*& fd);

		void write_segment_graph(_IO_FILE*& fd);

		int write_binary(std::ofstream* ofs);
	
		int read_binary(std::ifstream* ifs);

		void add_bias_counts(vector<int>* vec);

		bool is_annotated(int i);

		bool is_acceptor_ss(int i);

		bool is_donor_ss(int i);

		bool is_initial(int i);

		bool is_terminal(int i);

		float get_coverage_seg(int i);

		void update_coverage_information();

		vector<int> find_max_path();
		
		void compute_seg_cov();
	
		int get_read_starts(int from, int to);

		int get_read_ends(int from, int to);

		void find_tss_and_cleave(vector<int>* pos, vector<int>* starts, vector<int>* stops, float pval);

	private: 
		char* seq;

};

#endif
