#ifndef _TRANSCRIPT_PREDICTION_H__
#define _TRANSCRIPT_PREDICTION_H__

struct Config
{
	vector<char*> bam_files;
	char* fn_quant;
	char* fn_graph;
	char* fn_gtf;
	bool mm_filter;
	int intron_len_filter;
	int filter_mismatch;
	int exon_len_filter;
	int max_num_trans;
	int min_trans_len;
	bool use_pair;

	// loss
	float eta1;
	float eta2;
	int lambda;

	// regularizer
	float C_exon; 
	float C_intron;
	float C_pair;
	float C_num_trans;
	float C_num_trans_predef;

};

class Tr_Pred
{
	public:
		QP* qp;
		Bam_Region* graph;
		vector<vector<vector<float> > >* all_admat;
		vector<vector<vector<int> > >* all_pair_mat;
		vector<vector<float>* >* all_seg_cov;
		const Config* config;
		vector<vector<double> >* loss_param;


		void make_qp();
		void get_loss_coef(float (&ret)[4], float obs);
};

#endif
