#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include <assert.h>
#include "get_var.h"


double gammaln(double x)
{
	double y,tmp,ser;
	static double cof[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5 };
 
 	y=x;
  	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
 	ser = 1.000000000190015;
	for (int j=0; j<=5; j++)
		ser+= cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}

/*
 *
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs!=1)
		mexErrMsgTxt("Expected 1argument\n");

	double n = get_double(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* ret = mxGetPr(plhs[0]);

	if (n<0)
		mexErrMsgTxt("gammaln function only defined for positive values");
	else
		*ret = gammaln(n); 
}



