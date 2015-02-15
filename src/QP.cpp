
#include <stdio.h>
#include <assert.h>
#include "QP.h"


QP::QP(int num)
{
	num_var = num;
}
double QP::compute_obj()
{
	double res = 0;
	Q.reset_it();
	while (true)
	{
		int i=0;
		int j=0;
		float val = Q.next(&i, &j);
		//printf("%i %i %.2f\n", i, j, val);
		if (i==-1)
			break;
		assert(result.size()>i && result.size()>j);
		res += result[i]*result[j]*val;
	}

	assert(F.size()==0 || result.size()==F.size());
	for (int i=0; i<F.size(); i++)
		res+= result[i]*F[i];

	return res;
}

int QP::line_search()
{
	double eps = 1e-4;

	result = vector<double>(lb.size(), 1);

	printf("obj:%.5f\n", compute_obj());
	for (int i=0; i<num_var; i++)
	{
		
		if (lb[i]+eps>ub[i])
			continue;

		double best = compute_obj();
		double best_val = lb[i];
		double from = lb[i];
		double to = ub[i];
		if (from<-10)
			from = -10;

		if (to>10)
			to = 10;

		double step = (to-from)/10;
		for (double val=from; val<to; val+=step)
		{
			result[i] = val;
			double res = compute_obj();
			if (res<best)
			{
				best_val = val;
				best = res;
			}
		}
		result[i] = best_val;
		//printf("obj:%.5f\n", compute_obj());
	}
}
