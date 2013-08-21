#ifndef _VECTOR_OP_H__
#define _VECTOR_OP_H__

#include <assert.h>

template<typename T>
T max(vector<T>* vec)
{
	assert(vec->size()>=1);
	T ret = vec->at(0);
	for (unsigned int i=0; i<vec->size(); i++)
	{
		if (vec->at(i)>ret)
			ret = vec->at(i);
	}
	return ret;
}

template<typename T>
void mult(vector<T>* vec, T val)
{
	for (int i=0; i<vec->size(); i++)
	{
		(vec->at(i))*=val;
	}
}

vector<int> range(int lb, int ub)
{
	vector<int> ret;
	for (int i=lb; i<=ub; i++)
		ret.push_back(i);
	
	return ret;
}
#endif
