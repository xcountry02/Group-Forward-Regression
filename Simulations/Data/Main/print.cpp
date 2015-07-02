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

        /*-------------------------------------------------- print ------------
         |  Function print
         |
         |  Purpose:  Given the correct inputs this function 'prints' the response
         |            vector and design matrix to specific simulation settings.  
         |            This function works with only a fixed group size.  
         |
         |  Parameters: p0      (IN) -- Tne total number of groups
         |              n_grp   (IN) -- Group size (3-5) (fixed)
         |              n       (IN) -- Sample size
         |              sigma   (IN) -- Variance on sample (>0)
         |              name    (IN) -- Subscript name
         |              var_ind (IN) -- 0 -- AR1 variance structure
         |                           -- 1 -- IID variance structue
         |              oracle  (IN) -- 0 -- Ouput Data Set
         |                           -- 1 -- Output Result Cards
         |
         |  Returns:  If oracle==0, then the output is two files, a response (Y)
         |                and the explanatory variables (X), which are all factors.
         |                If oracle==1, then th output is still two files, but now
         |                the files are result cards, which have a very specific
         |                format.
         *-------------------------------------------------------------------*/


void print(int p0, float n_grp, int n, float sigma, int name, int var_ind, int oracle)
{

  int i, j, k;
  int size = n_grp - 1; // Group size
  int rho = 0; // The variance structure is IID be default
  int iter = 100; // Number of iteractions for oracle results
  int tNgrp = 0; // True number of groups (changes depending on group size)

  fmat samp(n, p0); // The sample data set created
  fvec b(size); // True beta
  fvec epsilon(n); // Error
  fvec y(n); // Response
  fvec hbeta(12); // For Oracle, the estimated beta, size 12 b/c always 12 parameters
                  // 1*12, 2*6, 3*4, 4*3,
  fvec calc(iter); // For oracle, MSE results stored here

  /*
    Setup the distributions.

    'var' is 0 then AR1 structure, if 1 independent
  */
 
  if (var_ind == 0) {
    rho = 0.5;
  } else {
    rho = 0;
  }

  //  Separate into groups based on normal distributions

  for (float a = 0.0; a < size; a+=1.0) {
    float tmp = (float) (a + 1) / n_grp;
    b(a) = NORMSINV(tmp, 0, 1);
  }

  // The key to a random seed is to do it once per program.  Don't do any
  //   of the randomness in stat.cpp

  struct timeval time;
  gettimeofday(&time,NULL);

  // microsecond has 1 000 000
  // Assuming you did not need quite that accuracy
  // Also do not assume the system clock has that accuracy.
  srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

  // The trouble here is that the seed will repeat every
  // 24 days or so.

  // If you use 100 (rather than 1000) the seed repeats every 248 days.

  // Do not make the MISTAKE of using just the tv_usec
  // This will mean your seed repeats every second.


  /*
      This sets the true beta, must be hand entered.
  */
  
  fvec beta(p0*size);

  for (i = 0; i < p0*size; i++) {
    beta(i) = 0;
  }

  if (size == 1) {
    beta(0) = 3;
    beta(1) = 3;
    beta(2) = 3;
    beta(3) = 3;
    beta(4) = 3;
    beta(5) = 3;
    beta(6) = 2;
    beta(7) = 2;
    beta(8) = 2;
    beta(9) = 2;
    beta(10) = 2;
    beta(11) = 2;
    tNgrp = 12;
  } else if(size == 2) {
    beta(0) = 3;
    beta(1) = 3;

    beta(2) = 3;
    beta(3) = 3;

    beta(4) = 3;
    beta(5) = 3;

    beta(6) = 3;
    beta(7) = 3;

    beta(8) = 3;
    beta(9) = 3;

    beta(10) = 3;
    beta(11) = 3;

    tNgrp = 6;
  } else if(size == 3) {
    beta(0) = 3;
    beta(1) = 3;
    beta(2) = 3;

    beta(3) = 3;
    beta(4) = 3;
    beta(5) = 3;

    beta(6) = 3;
    beta(7) = 3;
    beta(8) = 3;

    beta(9) = 3;
    beta(10) = 3;
    beta(11) = 3;

    tNgrp = 4;
  } else if(size == 4) {
    beta(0) = 3;
    beta(1) = 3;
    beta(2) = 3;
    beta(3) = 3;

    beta(4) = 3;
    beta(5) = 3;
    beta(6) = 3;
    beta(7) = 3;

    beta(8) = 3;
    beta(9) = 3;
    beta(10) = 3;
    beta(11) = 3;

    tNgrp = 3;
  }

  if (oracle == 0) {


// This section is to create the Sigma matrix when using the multivariate
//   normal distribution.  It is not necessary because I am using a much
//   faster way to sample from the rmvnorm which is not as memory intensive

/*
  fmat g_sigma(p0, p0);
  fvec mu(p0);
  if (var_ind == 0) {
    for (i = 0; i< p0; i++) {
      mu(i) = 0;
      for (j = 0; j < p0; j++) {
        g_sigma(i, j) = pow(0.5, abs(i - j));
      }
    }
  } else {
    for (i = 0; i< p0; i++) {
      mu(i) = 0;
      for (j = 0; j < p0; j++) {
        if (i == j) {
          g_sigma(i, j) = 1;
        } else {
          g_sigma(i, j) = 0;
        }
      }
    }
  }
*/

//printf("b(0) = %lf\n", b(0));

// This part is used if we have SNR -- Signal to Noise Ratio.  We sample
//   from this model built and try to estimate the true variance.
//
// This part will not be used because instead of using SNR, I am changing
//   the standard deviation of the epsilon's directly
//
//  Only run this the first time (Finding true_var, very slow)
/*  

    int sdn = 20000;
    fmat sdd(sdn, p0);

//    sdd = rmvn_Choleski(sdn, mu, g_sigma, p0);
    sdd = rmvn(sdn, p0, rho);

    for ( k = 0; k < sdn; k++) {
      for ( i = 0; i < p0; i++) {
        for ( j = 0; j < size; j++) {
          if (sdd(k, i) <= b(j)) {
            sdd(k, i) = j;
            break;
          } else if (j == size - 1) {
            sdd(k, i) = j + 1;
            break;
          }
        }
      }
    }

    fmat sdX(sdn, p0*(size));

    for ( i = 0; i < sdn; i++) {
     int cnt = 0;
      for ( j = 0; j < p0*(size); j+=size){
        for ( k = 0; k < size; k++) {
          if (sdd(i, cnt) == k){
            sdX(i, j + k) = 1;
          } else {
            sdX(i, j + k) = 0;
          }
        }
        cnt++;
      }
    }
 
    fvec sdY(sdn);

    sdY = sdX * beta;
    float samp_var = var(sdY);
    true_var = samp_var / SNR;

    printf("True Variance: %f\n", true_var);

    if (true_var == -1) {
      if (size == 1) {
        true_var = 0.3663755;
      } else if(size == 2) {
        true_var = 1.061375;
      } else if(size == 3) {
        true_var = 1.108877;
      } else if(size == 4) {
        true_var = 1.059319;
      }
    }
*/

// Create random sample from multivariate normal
    samp = rmvn(n, p0, rho);

// Based upon previus limits already found, this section dichotomizes
// (or trichotomizes) the current continuous variable into and integer
// thus creating the group structure, by factorizing the samples.
    for ( k = 0; k < n; k++) {
      for ( i = 0; i < p0; i++) {
        for ( j = 0; j < size; j++) {
          if (samp(k, i) <= b(j)) {
            samp(k, i) = j;
            break;
          } else if (j == size - 1) {
            samp(k, i) = j + 1;
            break;
          }
        }
      }
    }

// In order to simulate the response the design matrix must be created.
// The difference between samp and X is X contains dummy variables and
// all 0's and 1's
    fmat X(n, p0*(size));

    if ( size == 1 ){
      X = samp;
    } else {
      for ( i = 0; i < n; i++) {
        int cnt = 0;
        for ( j = 0; j < p0*(size); j+=size){
          for ( k = 0; k < size; k++) {
            if (samp(i, cnt) == k){
              X(i, j + k) = 1;
            } else {
              X(i, j + k) = 0;
            }
          }
          cnt++;
        }
      }
    }

// Epsilon is the variance noise.  If sigma is smaller than that noise is
// smaller.
    epsilon = rnorm(0, sigma, n);

// y is created, using the linear model.
    y = X * beta + epsilon;

    epsilon.reset();
    beta.reset();

    ofstream ofp;
    stringstream NameX, Namey;
  
    Namey << "Y";
    Namey << name;
    Namey << ".";
    Namey << p0;
    Namey << ".";
    Namey << n_grp;
    Namey << ".";
    Namey << n;
    Namey << ".";
    Namey << sigma;
    Namey << "_";
    Namey << var_ind;
  
    ofp.open(Namey.str().c_str(), ios::out);
    ofp << y; // Return y
    ofp.close();
  
    NameX << "X";
    NameX << name;
    NameX << ".";
    NameX << p0;
    NameX << ".";
    NameX << n_grp;
    NameX << ".";
    NameX << n;
    NameX << ".";
    NameX << sigma;
    NameX << "_";
    NameX << var_ind;
  
    ofp.open(NameX.str().c_str(), ios::out);
    ofp << samp; // Return sample matrix
    ofp.close();

  } else {

    // The process done here is exactly the same as above, only instead of
    // returning data sets the information is stored and calculated
    for (int c = 0; c < iter; c++) {

      samp = rmvn(n, p0, rho);
    
      for ( k = 0; k < n; k++) {
        for ( i = 0; i < p0; i++) {
          for ( j = 0; j < size; j++) {
            if (samp(k, i) <= b(j)) {
              samp(k, i) = j;
              break;
            } else if (j == size - 1) {
              samp(k, i) = j + 1;
              break;
            }
          }
        }
      }

      fmat X(n, p0*size);
      if ( size == 1 ){
        X = samp;
      } else {
        for ( i = 0; i < n; i++) {
          int cnt = 0;
          for ( j = 0; j < p0*(size); j+=size){
            for ( k = 0; k < size; k++) {
              if (samp(i, cnt) == k){
                X(i, j + k) = 1;
              } else {
                X(i, j + k) = 0;
              }
            }
            cnt++;
          }
        }
      }

      samp.reset();
    
      epsilon = rnorm(0, sigma, n);
    
      y = X * beta + epsilon;

      epsilon.reset();

      fmat oX(n, 12); 

      for (i = 0; i < 12; i++) {
        oX.col(i) = X.col(i);
      }
      
      X.reset();

      hbeta = coef(y, oX, 12);

      oX.reset();

      calc(c) = MSE(hbeta, beta, 12);

      y.reset();
      hbeta.reset();
    }

   float m=mean(calc);
   float s=stddev(calc);

   if (var_ind == 0) {
     var_ind = 3;
   } else {
     var_ind = 2;
   }

   // Create Result Cards
   for (int mbic = 0; mbic <= 1; mbic++) {
//     for (int ibic = 0; ibic <= 1; ibic++) {

        ofstream ofp;
        stringstream nam;

        nam << "summary_result.";
        nam << p0;
        nam << ".";
        nam << n;
        nam << ".";
        nam << n_grp;
        nam << ".";
        nam << sigma;
        nam << ".";
        nam << mbic;
        nam << ".";
        nam << var_ind;
        nam << ".out";

        ofp.open(nam.str().c_str(), ios::out);

        if (mbic == 1){
          ofp << "Main BIC = Small BIC calculation" << endl;
        } else {
          ofp << "Main BIC = Large BIC calculation" << endl;
        }
        ofp << "Total number of Groups = " << p0 << endl;
        ofp << "Total number of Parameters = " << p0*size << endl;
        ofp << "Group Size = " << n_grp << endl;
        ofp << "Sample Size = " << n << endl;
        ofp << "Iterations = " << iter << endl;
        if (var_ind == 3){
          ofp << "Variance Structure = AOR" << endl;
        } else {
          ofp << "Variance Structure = IOR" << endl;
        }
        ofp << "Sigma = " << sigma << endl;
        ofp << "Coverage probability = 1.00 (0.00)" << endl;
        ofp << "Percentage of correct zeros = 1.00 (0.00)" << endl;
        ofp << "Percentage of incorrect zeros = 0.00 (0.00)" << endl;
        ofp << "Exact select probability = 1.00 (0.00)" << endl;
        ofp << "Model Size = " << tNgrp << " (0.00)" << endl;
        ofp << "MSE = " << m << " (" << s << ")" << endl;
        ofp << "Total time = 0.00 (0.00)" << endl;
           
        ofp.close();

//      }
   }

  }
}

