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

void parse_gff_file(char* fn_gff, vector<Region*>* genes)
{
	FILE* fd = fopen(fn_gff, "r");

	vector<segment> transcript; 

	int num_exons_prev = 0;
	char prev_chr[100];
	char prev_strand;
	while (!feof(fd))
	{
		char line[1000];
		char* ret = fgets(line, 1000, fd);
		if (!ret)
			break;

		vector<char*> fields = separate(line, '\t');	

		if (fields.size()!=13)
		{
			printf("Warning: found invalid gff line: %i fields\n", (int) fields.size());
			for (int i=0; i<fields.size(); i++)
				printf("\t%s", fields[i]);
			printf("\n");
			break;
		}

		int num_exons = atoi(fields[10]);
		int exon_id = atoi(fields[9]);
		int start = atoi(fields[3]);
		int stop = atoi(fields[4]);

		if (exon_id==0 && transcript.size()>0)
		{
			//printf("Found transcript with %i exons\n", (int) transcript.size());
			assert(transcript.size()==num_exons_prev);

			int start = transcript.front().first;
			int stop = transcript.back().second;
			char *chr = new char[100];
			strcpy(chr, prev_chr);
			Region* reg = new Region(start, stop, chr, prev_strand);
			reg->transcripts.push_back(transcript);
			genes->push_back(reg);
			transcript.clear();
		}	
		segment* seg = new segment(start, stop);
		assert(start>0 && stop>0);
		transcript.push_back(*seg);
		num_exons_prev = num_exons;
		prev_strand = fields[5][0];
		strcpy(prev_chr, fields[0]);
	}

	if (transcript.size()>0)
	{
		int start = transcript.front().first;
		int stop = transcript.back().second;
		char *chr = new char[100];
		strcpy(chr, prev_chr);
		Region* reg = new Region(start, stop, chr, prev_strand);
		reg->transcripts.push_back(transcript);
		genes->push_back(reg);
	}
	fclose(fd);
}
bool compare_intron(segment intr1, segment intr2, int tol)
{
	bool acc = false;
	bool don = false;
	for (int i=-tol; i<tol; i++)
	{
		don = don || intr1.first+i==intr2.first;
		acc = acc || intr1.second+i==intr2.second;
	}
	return don && acc;
}
bool all_intron_compare(vector<segment> trans1, vector<segment> trans2)
{
		
	vector<segment> introns1; 
	vector<segment> introns2; 
	
	for (uint i=0; i<trans1.size()-1; i++)
	{
		if (trans1[i].second==trans1[i+1].first || trans1[i].second==trans1[i+1].first-1)
			continue; 
		introns1.push_back(segment(trans1[i].second, trans1[i+1].first)); 
	}
	for (uint i=0; i<trans2.size()-1; i++)
	{
		if (trans2[i].second==trans2[i+1].first || trans2[i].second==trans2[i+1].first-1)
			continue; 
		introns2.push_back(segment(trans2[i].second, trans2[i+1].first)); 
	}

	if (introns1.size() != introns2.size())
		return false;

	for (uint i=0; i<introns1.size(); i++)
	{
		if (introns1[i].first!=introns2[i].first)
			return false;
		if (introns1[i].second!=introns2[i].second)
			return false;
	}
	return true;
}

void eval_all_intron_level(Region* pred, Region* anno, int* match, int* correct)
{
	bool all_correct[pred->transcripts.size()]; 
	bool all_matched[anno->transcripts.size()]; 

	for (uint i=0; i<pred->transcripts.size(); i++)
		all_correct[i] = false;
	for (uint j=0; j<anno->transcripts.size(); j++)
		all_matched[j] = false;
	
	for (uint i=0; i<pred->transcripts.size(); i++)
	{
		for (uint j=0; j<anno->transcripts.size(); j++)
		{
			if (all_intron_compare(pred->transcripts[i], anno->transcripts[j]))
			{
				all_correct[i] = true;
				all_matched[j] = true;
			}
		}
	}

	for (uint i=0; i<pred->transcripts.size(); i++)
		*correct += all_correct[i];
	for (uint j=0; j<anno->transcripts.size(); j++)
		*match += all_matched[j];
}

void eval_genes(Region* pred, Region* frag, int* match, int* correct, int* compatible)
{
	assert(pred->transcripts.size()==1);
	assert(frag->transcripts.size()==1);

	// compute overlap fraction
	int ov_start = std::max(pred->start, frag->start);
	int ov_stop = std::min(pred->stop, frag->stop);

	int ov_len = 0;
	int pred_len = 0;
	for (int i=0; i<pred->transcripts[0].size(); i++)
	{
		pred_len += pred->transcripts[0][i].second-pred->transcripts[0][i].first;

		if (pred->transcripts[0][i].first<=ov_start && pred->transcripts[0][i].second>=ov_stop)
			ov_len += ov_stop-ov_start;// ov completely within one exon
		else if (pred->transcripts[0][i].first<=ov_start && pred->transcripts[0][i].second>=ov_start)
			ov_len += pred->transcripts[0][i].second-ov_start; // first overlapping exon
		else if ((pred->transcripts[0][i].first<=ov_stop && pred->transcripts[0][i].second>=ov_stop))
			ov_len += ov_stop-pred->transcripts[0][i].first; // last overlapping exon
		else if (pred->transcripts[0][i].first>=ov_start && pred->transcripts[0][i].second<=ov_stop)
			ov_len += pred->transcripts[0][i].second-pred->transcripts[0][i].first;
	}

	float frag_len = 0;
	for (int i=0; i<frag->transcripts[0].size(); i++)
	{
		frag_len += frag->transcripts[0][i].second-frag->transcripts[0][i].first;
	}

	float ov1 = ((float) ov_len)/pred_len;
	float ov2 = ((float) ov_len)/frag_len;

	
	// single exon genes
	if (pred->transcripts[0].size()==1 && frag->transcripts[0].size())
	{

		if (ov1>0.8)
			*correct = 1;

		if (ov2>0.8)
			*match = 1;

		return;
	}

	vector<segment> pintrons;
	for (int i=1; i<pred->transcripts[0].size(); i++)
	{
		segment intron(pred->transcripts[0][i-1].second, pred->transcripts[0][i].first);
		pintrons.push_back(intron);
	}
	vector<segment> fintrons;
	for (int i=1; i<frag->transcripts[0].size(); i++)
	{
		segment intron(frag->transcripts[0][i-1].second, frag->transcripts[0][i].first);
		fintrons.push_back(intron);
	}

	int i = 0;
	int j = 0;
	while (i<pintrons.size() && j<fintrons.size())
	{
		if (pintrons[i].second<fintrons[j].first)
		{
			i++;
			continue;
		}
		if (fintrons[j].second<pintrons[i].first)
		{
			j++;
			continue;
		}
		break;
	}
	int intron_tol = 3;
	int ov_cnt=0;
	while (i<pintrons.size() && j<fintrons.size())
	{
		//if (pintrons[i].first==fintrons[j].first && pintrons[i].second==fintrons[j].second)
		if (compare_intron(pintrons[i], fintrons[j], intron_tol))
		{
			ov_cnt++;
			i++;
			j++;
			continue;
		}
		else
		{
			// incompatible
			return;
		}
		break;
	}

	*compatible = 1;

	assert(ov_len<=ov_stop-ov_start);
	
	if (ov_cnt==pintrons.size() && ov_cnt==fintrons.size())
	{
		// perfect match
		*match = 1;
		*correct = 1;
		return;
	}

	// do not count as correct if fragment intron is outside of match
	for (int i=0; fintrons.size(); i++)
	{
		if (fintrons[i].first < ov_start || fintrons[i].second>ov_stop)
		{
			return;
		}
	}

	if (ov1>0.8)
		*correct = 1;

	if (ov2>0.8)
		*match = 1;

	return;
}

void print_exons(Region* reg)
{
	for (int j=0; j<reg->transcripts.size(); j++)
	{
		printf("\ttranscript %i\n", j); 
		for (int i=0; i<reg->transcripts[j].size(); i++)
		{
			printf("\t\t%i->%i\n", reg->transcripts[j][i].first, reg->transcripts[j][i].second);
		}
	}
}


