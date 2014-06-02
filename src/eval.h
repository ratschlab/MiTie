
#ifndef __EVAL_GTF_H__
#define __EVAL_GTF_H__

void parse_gff_file(char* fn_gff, vector<Region*>* genes);
bool compare_intron(segment intr1, segment intr2, int tol); 
bool all_intron_compare(vector<segment> trans1, vector<segment> trans2);
void eval_all_intron_level(Region* pred, Region* anno, int* match, int* correct); 
void eval_genes(Region* pred, Region* frag, int* match, int* correct, int* compatible); 
void print_exons(Region* reg);
#endif
