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

		sparse_matrix(sparse_matrix<T>* mat)
		{
			mat->reset_it(); 
			int i=0; 
			int j=0; 
			T val; 
			for (;;)
			{
				val = mat->next(&i,&j);
				if (i<0)
					break; 
				set(i, j, val); 
			}
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

		void mult(vector<T>* vec, vector<T>* result)
		{
			mult(vec, result, 1.0); 
		}

		void mult(vector<T>* vec, vector<T>* result, T weight)
		{
			reset_it(); 
			int i=0; 
			int j=0; 
			T val; 
			for (;;)
			{
				val = next(&i,&j);
				if (i<0)
					break; 

				assert(i<result->size()); 
				assert(j<vec->size()); 

				result->at(i) += vec->at(j)*val*weight; 
			}
		}

		void mult(sparse_matrix<T>* mat, sparse_matrix<T>* result)
		{
			reset_it(); 
			int i1=0; 
			int j1=0; 
			int i2=0; 
			int j2=0; 

			for (;;)
			{
				T val1 = next(&i1,&j1); 
				if (i1<0)
					break; 
				
				//printf("%i, %i, %f\n", i1, j1, val1); 
					
				mat->reset_it(); 
				for (;;)
				{
					T val2 = mat->next(&i2,&j2);
					if (i2<0)
						break; 

					if (i2<j1)
						continue; 
					else if (i2>j1)
						break; 
					else if (i2==j1)
					{
						T val = result->get(i1, j2); 
						result->set(i1, j2, val+val1*val2); 
						//printf("%i, %i, %f, %i %i %f\n", i1, j1, val1, i2, j2, val2); 
					}
				}
			}
		}

		void print()
		{
			reset_it(); 
		
			int row=0; 
			int col=0; 
			int i=0; 
			int j=0; 
			T val; 
			for (;;)
			{
				val = next(&i,&j);
				if (i<0)
				{
					printf("\n"); 
					break; 
				}
				if (row<i)
				{
					col=0; 
					row++; 
					printf("\n"); 
				}

				while (row<i)
				{
					printf("%i:\t0.00, ...\n", row); 
					row++; 
				}
				//printf("%i:\t", row); 
				while (col++<j)
					printf("0.00\t"); 
				printf("%.2f\t", val); 

			}
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
