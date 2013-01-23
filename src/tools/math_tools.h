#ifndef _MATH_TOOLS_H__
#define _MATH_TOOLS_H__
#include <vector>


double gammaln(double x);

double factln(int n);

double bino_pdf(int k, int n, float p);

float bino_cfd(int cnt, int num, float p);

template <typename T>
std::vector<T> prctile(std::vector<T> vec, std::vector<float> pos);

#endif
