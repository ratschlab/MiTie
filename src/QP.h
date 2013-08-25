#ifndef _QP_H__
#define _QP_H__

#include <vector>
	using std::vector;
#include <map>
	using std::map;
	using std::pair;
#include <assert.h>

template <typename T>
class semi_sparse_matrix{

	private:
		vector<vector<int> > idx;
		vector<vector<T> > coef;
		vector<vector<int> >::iterator it; 
		vector<vector<int> >::iterator it2; 

	public:
		semi_sparse_matrix()
		{
			reset_it();
		}
		void reset_it()
		{
			it = idx.begin(); 
			it2 = idx.begin(); 
		}
		int next(vector<int>* p_idx, vector<T>* p_coef)
		{
			if (it==idx.end())
			{
				return -1;
			}

			p_idx = &(*it);
			p_coef = &(*it2);
			it++;
			it2++;
			return 1;
		}
		int get(vector<int>** p_idx, vector<T>** p_coef, int j)
		{
			if (j>=idx.size())
			{
				return -1;
			}
			*p_idx = &(idx[j]);
			*p_coef = &(coef[j]);
			return 1;
		}


		void set(int i, int j, T val)
		{
			assert(i>=0);
			assert(j>=0);
			if (i>=idx.size())
			{
				assert(i==idx.size());
				idx.push_back(vector<int>());
				coef.push_back(vector<T>());
			}
			idx[i].push_back(j);
			coef[i].push_back(val);
			assert(idx[i].size()==coef[i].size());
		}

		size_t size()
		{
			return idx.size();
		}
};

template <typename T>
class sparse_matrix{

	private:
		map<pair<int, int>, T> mat;
		typename map<pair<int, int>, T>::iterator it; 

	public:
		sparse_matrix()
		{
			reset_it();
		}
		void reset_it()
		{
			it = mat.begin(); 
		}
		T next(int* i, int* j)
		{
			if (it==mat.end())
			{
				*i = -1;
				*j = -1; 
				return -1.0;
			}

			*i = it->first.first;
			*j = it->first.second;
			T ret = it->second;
			it++;
			return ret;
		}

		void set(int i, int j, T val)
		{
			pair<int, int> p(i, j);
			mat[p] = val;
		}

		T get(int i, int j)
		{
			pair<int, int> p(i, j);
			typename map<pair<int, int>, T>::iterator it; 
			it = mat.find(p);
			if (it==mat.end())
			{
				return 0.0;
			}
			return it->second;
		}

		size_t size()
		{
			return mat.size();
		}
};
template <typename T>
class semi_sparse_3d_matrix{

	private:
		map<pair<int, int>, vector<T> > mat;
		typename map<pair<int, int>, vector<T> >::iterator it; 

	public:
		semi_sparse_3d_matrix()
		{
			reset_it();
		}

		void add(int i, int j, T val)
		{
			pair<int, int> p(i, j);
			mat[p].push_back(val);
		}

		T get(int i, int j, int k, bool& exist)
		{
			exist = true;
			pair<int, int> p(i, j);
			typename map<pair<int, int>, vector<T> >::iterator it; 
			it = mat.find(p);
			if (it==mat.end())
			{
				exist = false;
				return 0.0;
			}
			vector<T>* vec = &it->second;
			assert(vec->size()>k);

			return vec->at(k);
		}
		void reset_it()
		{
			it = mat.begin(); 
		}
		vector<T>* next(int* i, int* j)
		{
			if (it==mat.end())
			{
				*i = -1;
				*j = -1;
				return NULL;
			}

			*i = it->first.first;
			*j = it->first.second;
			vector<T>* ret = &it->second;
			it++;
			return ret;
		}

		size_t size()
		{
			return mat.size();
		}
};

class QP
{
	// notation: min_x x'Qx + F'x
	//   s.t.    A_{eq}x = b_{eq}
	//   		 A_{ieq}x <= b_{ieq}
	//			 lb<=x<=ub
	//  the partitioning of A and b into 
	//  equality and inequality constraints 
	//  is implemented using an index eq_idx 
	public:
		QP(int num);

		int num_var;
		sparse_matrix<float> Q;
		vector<float> F;
		semi_sparse_matrix<float> A;
		vector<float> b;
		vector<float> lb;
		vector<float> ub;
		vector<int> eq_idx;
		vector<int> binary_idx;

		vector<double> result;

		double compute_obj();

		int line_search();
};
#endif
