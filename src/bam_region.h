#ifndef _BAM_REGION_H__
#define _BAM_REGION_H__

#include "region.h"
#include "read.h"

class Bam_Region: public Region
{
	public:
		vector<CRead> all_reads; 
		vector<CRead*> reads; 
		bool reads_sorted;
		uint32_t* coverage;
		uint32_t* intron_coverage;
		vector<int> intron_counts;
		vector<float> seg_cov;
		vector<vector<float> > admat;
		vector<vector<int> > pair_mat;


		Bam_Region();
		Bam_Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname);
		Bam_Region(int pstart, int pstop, int pchr_num, char pstrand);
		Bam_Region(int pstart, int pstop, char* chr, char pstrand);
		Bam_Region(Bam_Region* reg);
		Bam_Region(Region* reg);
		
		~Bam_Region();

		void init();

		void clear_reads();
		float get_coverage_global(int pstart, int pstop);
		int get_intron_conf(int intron_start, int intron_stop);

		void get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter);
		void get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter, bool mm_filter);

		void compute_coverage();

		void compute_intron_coverage();

		void compute_intron_list();

		int write_HDF5(char* filename);

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
		void write_segment_graph(_IO_FILE*& fd);

		void generate_segment_graph(float seg_filter, float tss_pval);

		void compute_admat(vector<int> starts, vector<int> stops);

		int compute_pair_mat();

		vector<int> get_parents(int node, bool no_neighbors);

		vector<int> get_children(int node, bool no_neighbors);

		void print_segments_and_coverage(_IO_FILE*& fd);

};
#endif
