
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

