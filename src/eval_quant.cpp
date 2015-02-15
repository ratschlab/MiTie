
#include <stdio.h>
#include <stdlib.h>
#include <string>
	using std::string;
#include <assert.h>
#include <vector>
	using std::vector;
#include <map>
	using std::map;
	using std::pair;
#include "tools/basic_tools.h"

int main(int argc, char** args)
{
	if (argc<3)
	{
		printf("Usage: %s <fn_quant_label> <fn_quant_pred>\n", args[0]);
		exit(0);
	}

	char* fn_quant_label = args[1];
	char* fn_quant_pred = args[2];


	FILE* fd = fopen(fn_quant_label, "r");
	if (!fd)
	{
		printf("could not open file: %s\n", fn_quant_label); 
		exit(-1);
	}

	map<string, float> label;
	while (~feof(fd)) 
	{
		char line[1000];
		if (fgets(line, 1000, fd)==NULL) break;

		vector<char*> fields = get_fields(line);

		assert(fields.size()==2);

		pair<string, float> tmp(string(fields[0]), atof(fields[1]));

		label.insert(tmp);
	}
	fclose(fd);

	printf("num label trans: %lu\n", label.size());


	// parse prediction
	FILE* fd2 = fopen(fn_quant_pred, "r");
	if (!fd2)
	{
		printf("could not open file: %s\n", fn_quant_pred); 
		exit(-1);
	}

	//vector<pair<float, float> > trans_pairs;
	int correct = 0;
	int num_label = 0;
	int num_pred = 0;
	while (~feof(fd2)) 
	{
		char line[1000];
		if (fgets(line, 1000, fd2)==NULL) break;

		vector<char*> fields = get_fields(line);

		assert(fields.size()==2);

		map<string, float>::iterator it;
		it = label.find(string(fields[0]));
		float val = atof(fields[1]);

		if (val>1e-3)
			num_pred++;

		if (it!=label.end())
		{
			//printf("found %.3f, %.3f\n", it->second, atof(fields[1]));
			//trans_pairs.push_back(pair<float, float>(it->second, atof(fields[1])));
			//
			if (it->second>1e-3)
				num_label++;

			if (it->second>1e-3 && val>1e-3)
				correct++;
		}
		else
		{
			printf("not found\n");
		}
	}
	float SN = ((float) correct)/num_label;
	float SP = ((float) correct)/num_pred;

	printf("SN,%.3f, SP%.3f\n", SN, SP);
	fclose(fd2);
}
