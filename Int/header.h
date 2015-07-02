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
void quicksort(int *x, int first, int last);
void del(fmat X, int dimX, int n, int pos, int size, double *mFR, int **delpIND, int *delpLEN, int **delgIND, int *delgLEN);
fvec coefm(const fvec &y, const fmat &X, int dimX, int n, int nump, int *mFR);
void coefi(const fvec &y, const fmat &X, int n, int ncov, int nump, int numg, int numip, int *size, int *cNum, int *MmFR, int *ImFR, double *beta);
double RSS(const fvec &y, const fmat &X, int dim);
int which(double *q, int dim, double num);
double min(double *q, int dim);
