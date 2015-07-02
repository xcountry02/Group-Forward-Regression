#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED
#endif

#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif

#include <iostream>
#include <armadillo>

using namespace arma;

template<class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}

template void print_matrix<arma::fmat>(arma::fmat matrix);
template void print_matrix<arma::cx_fmat>(arma::cx_fmat matrix);
// call print_matrix<arma::Mat<double> >(matrix) 

float rnorm(float, float);
float rgamma(float, float);
fmat rmvn_Choleski(int n, const fvec &mu, const fmat &sigma, int dim);
double RSS(const fvec &y, const fmat &X, int dim);
int which(const vec &q, int dim, double num);
void gFR(int numc, const fmat &X, const fvec &y, int n, int ncov, int size, int BICc, float** par, int* lenp, float** grp, int* leng);
float NORMSINV(float p, float mu, float sigma);
