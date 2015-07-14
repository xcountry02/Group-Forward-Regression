#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

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

       /*-------------------------------------------------- quicksort --------
         |  Function quicksort
         |
         |  Purpose:  The function takes 'x', a vector of doubles and sorts
         |            the vector from least to most using the algorithm
         |            quicksort.  The vector between 'first' and 'last' is
         |            sorted.
         |
         |  Parameters: x (IN/OUT) -- a pointer to a vector of doubles
         |              first (IN) -- an integer indicating the first element
         |                            of the vector one wants sorted
         |              last (IN) -- an integer indicating the last element
         |                           of the vector to be sorted
         |
         |  Returns:  x -- as a sorted vector.
         *-------------------------------------------------------------------*/

void quicksort(int *x,int first,int last){

    int pivot,j,temp,i;

    if(first<last){
        pivot=first;
        i=first;
        j=last;

        while(i<j){
            while(x[i]<=x[pivot]&&i<last)
                i++;
            while(x[j]>x[pivot])
                j--;
            if(i<j){
                temp=x[i];
                x[i]=x[j];
                x[j]=temp;
            }
        }
  
        temp=x[pivot];
        x[pivot]=x[j];
        x[j]=temp; 
        quicksort(x,first,j-1);
        quicksort(x,j+1,last);
    } 
} 

        /*-------------------------------------------------- del --------------
         |  Function del
         |
         |  Purpose:  Given a design matrix (X), significant parameters (mFR) and 
         |            a response the function makes sure there are no duplicate 
         |            columns.  Those duplicate columns are removed and their 
         |            position reported.
         |
         |  Parameters: X (IN) -- The design matrix -- exluding covariates
         |              dimX (IN) -- original number of parameters in X
         |              n (IN) -- sample size
         |              pos (IN) -- The number of groups chosen from GFR
         |              size (IN) -- The group size
         |              mFR (IN) -- The parameters found significant from the
         |                          group forward regression.  The parameter columns
         |                          will be taken from X.
         |              delpIND (OUT) -- The parameters deleted from X because
         |                               of duplicate columns.
         |              delpLEN (OUT) -- The number of parameters deleted from
         |                               X because of deplicate columns
         |              delgIND (OUT) -- The groups deleted from X because
         |                               duplicate columns
         |              delgLEN (OUT) -- The number of groups delete from X 
         |                               because duplicate columns
         |
         |  Returns:  Nothing
         *-------------------------------------------------------------------*/

void del(fmat X, int dimX, int n, int pos, int size, double *mFR, int **delpIND, int *delpLEN, int **delgIND, int *delgLEN) {
  // The number of parameters is the number of groups * group size
  int nump = pos*size;

  fmat tmpX(n, nump);

  int *delp, *delg;
  delp = (int *) malloc(sizeof(int) * nump); // Keeps track of deleted parameters
  delg = (int *) malloc(sizeof(int) * pos); // Keeps track of deleted groups

  int cntDELp = 0; // Counts the number of deleted parameters
  int cntDELg = 0; // Counts the number of deleted groups

  int cnt = 0;
  // Go through every column of design matrix X
  for (int i = 0; i < dimX; i++) {
    // Check the column number (i) with mFR (j), the columns to be kept
    for (int j = 0; j < nump; j++) {
      // If the column is to be saved, then put it in tmpX
      if (mFR[j] == i) {
        // Check to see if column hasn't already been added
        int flag = 0;
        if (cnt > 0) { // Not the first column to be added
          for (int k = 0; k < cnt; k++) {
            int sum = 0;
            for (int a = 0; a < n; a++) {
              if (tmpX(a, k) == X(a, i)) {
                sum++;
              }
            }
            // If the column has been added, add this to the delp vector
            if (sum == n) {
              flag = 1;
              // Add the group to the delg vector
              int g = (j - (j % size)) / size;
              delg[cntDELg++] = g;

              // Must delete all columns in the same group as 'j'
              delp[cntDELp++] = j;
              int tmpJ = j;
              int tmpFLAG = (tmpJ % size) + 1;
              while (tmpFLAG != 1) {
                tmpJ--;
                delp[cntDELp++] = tmpJ;
                tmpFLAG = (tmpJ % size) + 1;
                cnt--; // Already added previous column in group, must be deleted
                tmpX.shed_col(cnt); // It is the last column added to tmpX
              }

              tmpFLAG = (j % size) + 1;
              while (tmpFLAG != 1) {
                j++;
                i++; // If more in group skip the rest of the group
                delp[cntDELp++] = j;
                tmpFLAG = (j % size) + 1;
              }
              break;
            }
          }
        }
        if (flag == 0) { // If no identical column then add it
          tmpX.col(cnt++) = X.col(i);
        }
      }
    }
  }
  delp[cntDELp] = 0;
  delg[cntDELg] = 0;
  *delpIND = delp;
  *delpLEN = cntDELp;
  *delgIND = delg;
  *delgLEN = cntDELg;

  X = tmpX;
}

        /*-------------------------------------------------- coefm ------------
         |  Function coefm
         |
         |  Purpose:  Given a design matrix (X), significant parameters (mFR) and 
         |            a response the beta coefficients are calculated. 
         |
         |  Parameters: y (IN) -- The response vector
         |              X (IN) -- The design matrix -- exluding covariates
         |              dimX (IN) -- original number of parameters in X
         |              n (IN) -- sample size
         |              pos (IN) -- The number of groups chosen from GFR
         |              size (IN) -- The group size
         |              mFR (IN) -- The parameters found significant from the
         |                          group forward regression.  The parameter columns
         |                          will be taken from X.
         |
         |  Returns:  The beta coefficients that are calculated.
         *-------------------------------------------------------------------*/

fvec coefm(const fvec &y, const fmat &X, int dimX, int n, int nump, int *mFR) {
  // The number of parameters is the number of groups * group size
  quicksort(mFR, 0, nump - 1);
  fmat tmpX(n, nump);

  int colNUM=0;
  for (int i=0; i < dimX; i++) {
    for (int j=0; j < nump; j++) {
      if (i == mFR[j]) {
        tmpX.col(colNUM) = X.col(i);
        colNUM++;
        break;
      }
      if (i < j){
        break;
      }
    }
  }

  fvec coef(nump);
  coef = inv_sympd(tmpX.t() * tmpX) * tmpX.t() * y; // Calculate coefficients

  fvec ans(dimX); // Most of the beta's are zero

  for (int i = 0; i < dimX; i++) {
    ans(i) = 0;
    for (int j = 0; j < nump; j++) { // For those in mFR, beta_j's have been calculated
      if (mFR[j] == i) {
        ans(i) = coef(j);
      }
    }
  }

  return ans;
}




        /*-------------------------------------------------- coefi ------------
         |  Function coefi
         |
         |  Purpose:  Given a design matrix (X), significant parameters (mFR) 
         |            and a response the beta coefficients are calculated.  
         |            This specific coefficient function is when interactions
         |            are in the model.
         |
         |  Parameters: y (IN) -- The response vector
         |              X (IN) -- The design matrix -- exluding covariates
         |              ncov (IN) -- original number of parameters in X
         |              n (IN) -- sample size
         |              nump (IN) -- The number of main parameters chosen from gFR
         |              numg (IN) -- The number of main groups chosen from gFR
         |              numip (IN) -- The number of interaction parameters chosen 
         |                            from gFR
         |              size (IN) -- A vector containing the size of each group,
         |                           found significant
         |              cNum (IN) -- A vector containing the first subscript of 
         |                           each group found significant
         |              MmFR (IN) -- The main parameters found significant 
         |                        from the group forward regression.  
         |              ImFR (IN) -- The interaction parameters found significant 
         |                        from the group forward regression.  
         |              beta (OUT) -- The beta coefficients are returned 
         |
         |  Returns:  The beta coefficients that are calculated.
         *-------------------------------------------------------------------*/
void coefi(const fvec &y, const fmat &X, int n, int ncov, int nump, int numg, int numip, int *size, int *cNum, int *MmFR, int *ImFR, double *beta) {

  // The number of parameters is the number of groups * group size
//  int nump = numg*size;
//  int isize = pow(size, 2);
//  int numi = numig*isize;
  int flag = 0;
//printf("dimX=%d\n", ncov);
//printf("n=%d\n", n);
//printf("numg=%d\n", numg);
//printf("size=%d\n", size);
//printf("numip=%d\n", numip);
//printf("nump=%d\n", nump);
//printf("icov=%d\n", icov);
//printf("isize=%d\n", isize);
//printf("X:\n");
//print_matrix(X);
  fvec ans(nump + numip);
//  fvec coef(nump + numi);

  if (numip > 0) {

    while(flag == 0) {
      fmat tmpX(n, nump + numip);

      for (int i = 0; i < nump; i++) {
        tmpX.col(i) = X.col(i);
      }
// Note: ImFR and MmFR are sorted.
      int cnt=0;
      int cnt2=0;
      for (int i = 0; i < numg; i++) {
        for (int j = i+1; j < numg; j++) {
          for (int a = cNum[i]; a < cNum[i] + size[i]; a++) {
            for (int b = cNum[j]; b < cNum[j] + size[j]; b++) {
              if (ImFR[cnt] == cnt2) {
                tmpX.col(cnt+nump) = X.col(i) % X.col(j);
                cnt++;
              }
              if (cnt == numip) {
                break;
              }
              cnt2++;
            }
          }
        }
      }

//printf("X:\n");
//print_matrix(tmpX);
      ans = inv_sympd(tmpX.t() * tmpX) * tmpX.t() * y; // Calculate coefficients
//printf("ans:\n");
//print_matrix(ans);
      fvec check = abs(ans);

      for (int i = 0; i < nump + numip; i++) {
        if ( check(i) > 25 ) {
          flag = 0;
        } else {
          flag = 1;
        }
      }

      tmpX.reset();
//
//      for (int i = 0; i < nump; i++) {
//        ans(i) = 0;
//        for (int j = 0; j < nump; j++) { // For those in mFR, beta_j's have been calculated
//          if (MmFR[j] == i) {
//            ans(i) = coef(j);
//          }
//        }
//      }
//
//      for (int i = 0; i < numi; i++) {
//        ans(nump + i) = coef(i+nump);
//      }
    }

  } else {

    ans = inv_sympd(X.t() * X) * X.t() * y; // Calculate coefficients
  }
//printf("ans:\n");
//print_matrix(ans);
//for (int j=0; j< nump; j++) {
//  printf("MmFR[%d]=%d\n", j, MmFR[j]);
//}

  for (int i=0; i<ncov+numip; i++) {
    beta[i]=0;
    for (int j=0; j < nump; j++) {
//    printf("i=%d\tj=%d\n", i, j);
      if (MmFR[j] == i) {
        beta[i] = ans(j);
        break;
//  printf("MmFR[%d]=%d\tbeta[%d]=%lf\n", j, MmFR[j], i, beta[i]);
      }
    }
    if (i>=ncov) {
      beta[i] = ans(i);
    }
  }
}

//fvec coefi(const fvec &y, const fmat &X, int dimX, int n, int numg, int size, double *mMFR, double *mFR, int numig ) {
//  // The number of parameters is the number of groups * group size
//  int nump = numg*size;
//  int isize = pow(size, 2);
//  int numi = numig*isize;
//  int flag = 0;
////printf("dimX=%d\n", dimX);
////printf("n=%d\n", n);
////printf("numg=%d\n", numg);
////printf("size=%d\n", size);
////printf("numig=%d\n", numig);
////printf("nump=%d\n", nump);
////printf("numi=%d\n", numi);
////printf("isize=%d\n", isize);
////printf("X1:\n");
////print_matrix(X);
//  fvec ans(dimX + numig*isize);
//  fvec coef(nump + numi);
//
//  if (numig > 0) {
//
//    while(flag == 0) {
//      fmat tmpX(n, nump + numi);
//
//      for (int i = 0; i < nump; i++) {
//        tmpX.col(i) = X.col(i);
//      }
//
//      int cnt=0;
//      int cnt2=0;
//      for (int i = 0; i < nump; i+=size) {
//        for (int j = i; j < nump; j+=size) {
//          if (i != j) {
//            for (int a=0; a< size; a++) {
//              for (int b=0; b<size; b++) {
//                if (mFR[cnt] == cnt2) {
//                  tmpX.col(cnt+nump) = X.col(i + a) % X.col(j + b);
//                  cnt++;
//                }
//                cnt2++;
//              }
//            }
//          }
//        }
//      }
//
////printf("X:\n");
////print_matrix(tmpX);
//      coef = inv_sympd(tmpX.t() * tmpX) * tmpX.t() * y; // Calculate coefficients
////printf("coef:\n");
////print_matrix(coef);
//      fvec check = abs(coef);
//
//      for (int i = 0; i < nump + numi; i++) {
//        if ( check(i) > 25 ) {
//          flag = 0;
//        } else {
//          flag = 1;
//        }
//      }
//
//      tmpX.reset();
//
//      for (int i = 0; i < dimX; i++) {
//        ans(i) = 0;
//        for (int j = 0; j < nump; j++) { // For those in mFR, beta_j's have been calculated
//          if (mMFR[j] == i) {
//            ans(i) = coef(j);
//          }
//        }
//      }
//
//      for (int i = 0; i < numi; i++) {
//        ans(dimX + i) = coef(i+nump);
//      }
//    }
//
//  } else {
//
//    coef = inv_sympd(X.t() * X) * X.t() * y; // Calculate coefficients
//
//    for (int i = 0; i < dimX; i++) {
//      ans(i) = 0;
//      for (int j = 0; j < nump; j++) { // For those in mFR, beta_j's have been calculated
//        if (mMFR[j] == i) {
//          ans(i) = coef(j);
//        }
//      }
//    }
//  }
//  return ans;
//}


        /*-------------------------------------------------- RSS --------------
         |  Function RSS
         |
         |  Purpose:  The Residual Sum of Squares.  Essentially this is the 
         |            error term.  A model with a smallest RSS is the better
         |            fitting model.
         |
         |  Parameters: y (IN) -- Response vector
         |              X (IN) -- Explanatory variable (design matrix)
         |              dim (IN) -- number of parameters
         |
         |  Returns:  a value (double) called RSS
         *-------------------------------------------------------------------*/

double RSS(const fvec &y, const fmat &X, int dim){

  fmat I;
  I.eye(dim, dim);
  fmat H(dim, dim);
  double err;
 
  //  First try the regular method of calculating RSS, it will work
  //    99% of the time.
 
  try
  {
    H = X * inv_sympd(X.t() * X) * X.t();

    err = as_scalar(y.t() * (I - H) * y);
  }
  catch (const std::exception& e)
  {
    //  Every once in awhile (< 1%) there is an X matrix that has a row
    //    of zeros, if this happens must calculate the RSS in a slightly
    //    convoluted way.  Sometimes (< .001%) the inverse cannot be
    //    calculated, if that happens the RSS cannot be calculated, return
    //    a large number.
    try
    {
      fmat tmp(X.n_cols, X.n_cols);
      fmat inver(X.n_cols, X.n_cols);
      fmat tmp2(X.n_cols, dim);

      tmp = X.t() * X;
      inver = inv_sympd(tmp);
      tmp.reset();
      tmp2 = inver * X.t();
      inver.reset();
      H = X * tmp2;
      tmp2.reset();

      err = as_scalar(y.t() * (I - H) * y);
    }
    catch (const std::exception& e)
    {
      err = std::numeric_limits<double>::max();
    }
  }

  return err;
}


        /*-------------------------------------------------- which ------------
         |  Function which
         |
         |  Purpose:  Finds the  index at which a value occurs
         |
         |  Parameters: q (IN) -- the vector of doubles
         |              dim (IN) -- length of q vector
         |              num (IN) -- value to be found in vector
         |
         |  Returns:  the index value, if not found -1
         *-------------------------------------------------------------------*/


int which(double *q, int dim, double num) {
  num = roundf(100000000000*num)/100000000000;
  for (int i = 0; i < dim; i++) {
    if (num == roundf(100000000000 * q[i])/100000000000)
      return i;
  }
  return -1;
}


        /*-------------------------------------------------- min --------------
         |  Function min
         |
         |  Purpose:  Finds the minimum value in a vector
         |
         |  Parameters: q (IN) -- the vector of doubles
         |              dim (IN) -- length of q vector
         |
         |  Returns:  returns the minimum value
         *-------------------------------------------------------------------*/

double min(double *q, int dim) {
        double tm;

        for (int i = 0; i < dim; i++) {
                if (i == 0) {
                        tm = q[i];
                } else {
                        tm = min(tm, q[i]);
                }
        }
        return tm;
}

