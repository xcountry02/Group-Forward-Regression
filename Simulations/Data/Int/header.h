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

int nCr(int n, int r);
fvec rnorm(float, float, int);
int nCr(int n, int r);
fvec coef(const fvec &y, const fmat &X, int dimX);
float MSE(const fvec &est, const fvec &tru, int dim);
// float rgamma(float, float);
// fmat rmvn_Choleski(int n, const fvec &mu, const fmat &sigma, int dim);
fmat rmvn(int n, int p, int rho);
float NORMSINV(float p, float mu, float sigma);
void print(int p0, float n_grp, int n, float SNR, int name, int var_ind, int oracle);
