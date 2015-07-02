#include <cstdlib>
#include <cmath>
#include <time.h>
#include "header.h"

using namespace std;
using namespace arma;


        /*-------------------------------------------------- nCR --------------
         |  Function nCr
         |
         |  Purpose:  This is the combination function.  It calculates how many
         |            subgroups of size r can be formed from the total number 
         |            of groups (n).
         |
         |  Parameters: n (IN) -- total number in group
         |              r (IN) -- size of subgroup to be chosen
         |
         |  Returns: number is calculated
         *-------------------------------------------------------------------*/




int nCr(int n, int r) {
   if(r>n) {
      printf("FATAL ERROR"); return 0;
     }
   if(n==0 || r==0 || n==r) {
      return 1;
   } else {
      return (int)lround( ((double)n/(double)(n-r)/(double)r) * exp(lgamma(n) - lgamma(n-r) - lgamma(r)));
   }
}






        /*-------------------------------------------------- MSE --------------
         |  Function MSE
         |
         |  Purpose:  This function calculates the Mean Square Error 
         |            between the extimated vector and true vector
         |
         |  Parameters: est (IN) -- estimate vector
         |              tru (IN) -- true vector
         |              dim (IN) -- length of the vectors
         |
         |  Returns: number is calculated
         *-------------------------------------------------------------------*/


float MSE(const fvec &est, const fvec &tru, int dim) {

  float sum=0;
  for (int i = 0; i < dim; i++) {
    float tmp = est(i) - tru(i);
    sum = pow(tmp, 2) + sum;
  }

  return sum;
}

        /*-------------------------------------------------- coef -------------
         |  Function coef
         |
         |  Purpose:  This function calculates the beta coefficients, using
         |            RSS method.
         |
         |  Parameters: y (IN) -- response vector
         |              X (IN) -- Design matrix
         |              dimX (IN) -- Dimension of beta vector
         |
         |  Returns: vector is returned
         *-------------------------------------------------------------------*/

fvec coef(const fvec &y, const fmat &X, int dimX) {
  // The number of parameters is the number of groups * group size

  fvec coef(dimX);
  coef = inv_sympd(X.t() * X) * X.t() * y; // Calculate coefficients

  return coef;
}

        /*-------------------------------------------------- rnorm -------------
         |  Function rnorm
         |
         |  Purpose:  This function is a random number generator from a 
         |            normal distribution.
         |
         |  Parameters: mu (IN) -- mean
         |              sigma (IN) -- standard deviation
         |              n (IN) -- sample size being generated
         |
         |  Returns: vector of numbers
         *-------------------------------------------------------------------*/

fvec rnorm(float mu, float sigma, int n){

  fvec ret(n);
  float tmp1;

  for (int i=0; i<n; i++) {
    tmp1 = rand() / (double)RAND_MAX;
    tmp1 = (pow(tmp1, 0.135) - pow(1 - tmp1, 0.135)) / 0.1975;
    ret(i)=mu + sigma * tmp1;
  }

  return ret;
}

        /*-------------------------------------------------- rgamma ------------
         |  Function rgamma
         |
         |  Purpose:  This function is a random number generator from a 
         |            gamma distribution.
         |
         |  Parameters: alpha (IN) -- 
         |              beta (IN) -- 
         |
         |  Returns: returns one number
         *-------------------------------------------------------------------*/

/*
float rgamma(float alpha=1, float beta=1){

  int flag = 0;
  if (alpha < 1) {
    alpha = alpha + 1;
    flag = 1;
  }

  float d = alpha - 1./3;
  float c = 1 / sqrt(9 * d);
  float Z = rnorm(0, 1);
  float V = pow(1 + c*Z, 3);


// float U = rand() / (double)RAND_MAX;
// float u1 = log(U);
// float v1 = log(V);
// float z1 = 0.5*pow(Z,2) + d - d*V + d * v1;
// float c1 = -1/c;
//
// while(!((Z > c1) && (u1 < z1))) {
//   U = rand() / (double)RAND_MAX;
//   Z = rnorm();
//   V = pow(1 + c*Z, 3);
// }


  float X = d*V;
  X = X/beta;
  if (flag == 0)
    return X;

  float U = rand() / (double)RAND_MAX;
  float X_p = X * pow(U, 1/(alpha-1));
  return X_p;
}
*/

/*
fmat rmvn_Choleski(int n, const fvec &mu, const fmat &sigma, int dim) {

  clock_t t=clock();
  srand((unsigned int) t);

  fmat Q(dim, dim);

  Q = chol(sigma);
  fmat A(n, dim), Z(n, dim);
  
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < n; j++) {
      Z(j, i) = rnorm(0, 1);
    }
  }
 
  for (int j = 0; j < n; j++) {
    A.row(j) = mu.t();
  }
  
  A = Z * Q + A;

  return (A);
}
*/

// Can I change mean and sd?
        /*-------------------------------------------------- rnorm -------------
         |  Function rnorm
         |
         |  Purpose:  This creates Multivariate Normal with mean 0, sd 1.  It is
         |            not as much of a memory hog, and is a lot faster.  It is
         |            stuck with a mean of 0 and a "nice: IID or AR(1) 
         |            covariance structure.
         |
         |  Parameters: n (IN) -- sample size
         |              p (IN) -- dimension of multivariate normal
         |              rho (IN) -- float b/t (0,1)
         |                     -- 0   -- IID variance structure
         |                     -- 0.5 -- AR(1) variance structure
         |
         |  Returns: matrix of numbers
         *-------------------------------------------------------------------*/

fmat rmvn(int n, int p, int rho) {

  fmat A(p, n);
  
  for (int i = 0; i < n; i++) {
    A.col(i) = rnorm(0, 1, p);
  }

  return A.t();
}

