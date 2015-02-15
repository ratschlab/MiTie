
#ifndef __LOSS_TANGENT_H__
#define __LOSS_TANGENT_H__

#include "math_tools.h"
void compute_mus(vector<float>* mus, float eta1, float eta2, float obs)
{
		float std_= sqrt(eta1*(obs+1) + pow(eta2*obs, 2));
		for (int j=-5; j<5; j++)
		{
			mus->push_back(j*std_);
		}
}

double compute_loss(float eta1, float eta2, int lambda, float obs, float mu)
{

	float var_ = eta1*mu + pow(eta2*mu, 2);
	
	// compute parameters of negative binomial
	double p = 1-(mu/var_);
	double r = mu*(1-p)/p;
	
	double logrp1 = factln(r-1); 
	double logrp2 = r*log(1-p); 
	double res; 
	double t1, t2, t3, t4, t5 = 0.0; 
	// A ~ Pois(lambda)
	// B ~ NB(r,p)
	// obs = A+B 
	//
	for (int i=0; i<100 && i<obs; i++)
	{
		double pois;
		if (lambda==0 && i==0)
		{
			pois = 0.0;
		}
		else if (lambda==0)
		{
			break;
		}
		else
		{
			pois = i*log(lambda) - lambda - factln(i); 
		}
	
		double nb;
		if (eta1==1 && eta2==0)
		{
			// special case: in the limit of var->mu we get the poisson distrib
			nb = (obs-i)*log(mu) - mu - factln(obs-i);
		}
		else
		{
			nb = factln(obs-i+r-1) - factln(obs-i) - logrp1 + logrp2 + (obs-i)*log(p); 
			//printf("%.3f %.3f %.3f %.3f %.3f\n", factln(obs-i+r-1),  -factln(obs-i), -factln(r-1), r*log(1-p), (obs-i)*log(p)); 
			//t1 += factln(obs-i+r-1); 
			//t2 += - factln(obs-i); 
			//t3 += -logrp1;
			//t4 += logrp2; 
			//t5 += (obs-i)*log(p); 

		}
	
		// naive: Y(k) = Y(k) + exp(pois+nb); 
		// using log(a+b) = log(a) + log(1+exp(log(b)-log(a)))
		if (i==0)
		{
			res = -(nb+pois);
			//if (obs==mu)
				//printf("nb:%.3f, pois:%.3f\n", nb, pois); 
		}
		else
		{
			res = res - (log(1+exp(nb+pois-res)));
		}
	}

	//printf("%.3f %.3f %.3f %.3f %.3f\n", t1, t2, t3, t4, t5); 
	return res; 
}


double compute_loss_deriv(float eta1, float eta2, int lambda, float obs, float mu, float h)
{

	double val1 = compute_loss(eta1, eta2, lambda, obs, mu-h); 
	double val2 = compute_loss(eta1, eta2, lambda, obs, mu+h); 

	//printf("compute_loss_deriv: val2:%.3f val1:%.3f, ret:%.3f\n", val2, val1, (val2-val1)/(2*h)); 

	return (val2-val1)/(2*h);

}

#endif
