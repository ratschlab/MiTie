#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "region.h"
#include "tools.h"
#include "gtf_tools.h"
#include "eval.h"

struct Config
{
	char* fn_pred;
	char* fn_anno;
	FILE* fd_out;
};

void parse_args(int argc, char** argv,  Config* c)
{
	if (argc<3)
	{
		fprintf(stderr, "Usage: %s fn_pred fn_anno ...", argv[0]);
		exit(1);
	}
	c->fn_pred = argv[1];
	c->fn_anno = argv[2];
	c->fd_out = stdout;
    for (int i = 2; i < argc; i++)  
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
	}
}

int main(int argc, char* argv[])
{
	Config c;
	parse_args(argc, argv, &c);

	// parse prediction
	//vector<Region*> pred_genes; 
	//parse_gff_file(c.fn_pred, &pred_genes);
	vector<Region*> pred_genes = parse_gtf(c.fn_pred);
	printf("Read %i genes from prediction file\n", (int) pred_genes.size());

	vector<Region*> anno_genes = parse_gtf(c.fn_anno);	
	printf("Read %i genes from annotation file\n", (int) anno_genes.size());


	vector<vector<int> > ov_list = region_overlap(pred_genes, anno_genes);

	int* correct = new int[pred_genes.size()];
	int* num_pred = new int[pred_genes.size()];
	int* matched = new int[anno_genes.size()];
	int* num_anno = new int[anno_genes.size()];

	memset(correct, 0, pred_genes.size()*sizeof(int));
	memset(num_pred, 0, pred_genes.size()*sizeof(int));
	memset(matched, 0, anno_genes.size()*sizeof(int));
	memset(num_anno, 0, anno_genes.size()*sizeof(int));

	printf("found %i overlapping genes\n", ov_list.size()); 

	int cnt=0;
	for (int i=0; i<ov_list.size(); i++)
	{
		num_pred[i] = pred_genes[i]->transcripts.size();
		for (int j=0; j<ov_list[i].size(); j++)
		{
			//if (pred_genes[i]->strand!=anno_genes[ov_list[i][j]]->strand)
			//	continue;
			if (strcmp(pred_genes[i]->chr, anno_genes[ov_list[i][j]]->chr)!=0)
				continue;
			int m = 0; 
			int c = 0; 
			int compatible = 0;
			//eval_genes(pred_genes[i], anno_genes[ov_list[i][j]], &m, &c, &compatible);
			eval_all_intron_level(pred_genes[i], anno_genes[ov_list[i][j]], &m, &c);
			correct[i] += c;
			matched[ov_list[i][j]] += m;
			num_anno[ov_list[i][j]] = anno_genes[ov_list[i][j]]->transcripts.size();

			if (false && c==0)
			{
				printf("Prediction: %i\n", c);
				print_exons(pred_genes[i]);
				printf("Annotation: %i\n", m);
				print_exons(anno_genes[ov_list[i][j]]);
				printf("\n");
			}
		}
	}
	int num_correct = 0;
	int num_trans_pred = 0;
	for (int i=0; i<pred_genes.size(); i++)
	{
		num_trans_pred += num_pred[i];
		num_correct+=correct[i];
	}
	int num_matched = 0;
	int num_trans_anno = 0;
	for (int i=0; i<anno_genes.size(); i++)
	{
		//num_trans_anno += num_anno[i];
		num_trans_anno += anno_genes[i]->transcripts.size();
		num_matched+=matched[i];
	}
	printf("%i/%i (%.4f) correct, %i/%i (%.4f) matched\n", num_correct, num_trans_pred, ((float) num_correct)/num_trans_pred, num_matched, num_trans_anno, ((float) num_matched)/num_trans_anno);
}



