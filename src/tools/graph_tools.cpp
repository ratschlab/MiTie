#include "graph_tools.h"
#include <stdio.h>
#include <stdlib.h>

int count_num_paths(vector<vector<float> > admat, int node1, int node2)
{
	if (node2<node1)
		return 0;
	if (node2==node1)
		return 1;
	if (node1>=(int) admat.size())
	{
		fprintf(stderr, "node1 out of bounds");
		exit(1);
	}
	if (node2>=(int) admat.size())
	{
		fprintf(stderr, "node2 out of bounds");
		exit(1);
	}
	int num_paths[node2-node1+1];
	num_paths[0] = 1;
	for (int i=node1+1; i<=node2; i++)
	{
		num_paths[i-node1] = 0;
		// sum over parents
		for (int j=node1; j<i; j++)
		{
			assert(num_paths[j-node1]>=0); 
			if (admat[j][i]>=-1)
				num_paths[i-node1] += num_paths[j-node1];
		}
	}
	return num_paths[node2-node1];
}

