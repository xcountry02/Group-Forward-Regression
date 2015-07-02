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
         |              SNR     (IN) -- Signal-to-noise-ratio on sample (>0)
         |              name    (IN) -- Subscript name
         |              var_int (IN) -- 0 -- AR1 variance structure
         |                           -- 1 -- IID variance structue
         |              oracle  (IN) -- 0 -- Ouput Data Set
         |                           -- 1 -- Output Result Cards
         |
         |  Returns:  If oracle==0, then the output is two files, a response (Y)
         |            and the explanatory variables (X), which are all factors.
         |            If oracle==1, then the output is four files, and
         |            the files are result cards, which have a very specific
         |            format.
         *-------------------------------------------------------------------*/

// There are not a lot of comments on this file.  To see more look at main.


void print(int p0, float n_grp, int n, float SNR, int name, int var_ind, int oracle)
{

  int i, j, k;
  int size = n_grp - 1; // Main group Size
  int rho = 0; // The variance structure IID is default
  int iter = 25; // Number of iterations for oracle results
  int tNgrp = 4; // True number of main groups (constant 4)
  int iSize = pow(size, 2); // Interaction group size
  int sub[48];
  int sdn = 50000;
  int sdp0 = 100;
  float true_var;
  string cSNR; // SNR
  int iSNR; // sigma

  // What this is a conversion between sigma (iSNR) and SNR (cSNR)

  float acceptedDiff = 0.0000001;
  if (fabsf(SNR-3) < acceptedDiff) {
    iSNR = 1;
    cSNR = "3";
  } else if (fabsf(SNR-1) < acceptedDiff) {
    iSNR = 2;
    cSNR = "1";
  } else {
    iSNR = 3;
    cSNR = "0.33";
  }

  fvec tbeta(size*tNgrp + iSize*3);
  fmat samp(n, p0);
  fvec b(size);
  fvec epsilon(n);
  fvec y(n);
  fvec hbeta(size*tNgrp + iSize*3);
  fvec calc(iter);

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
      This sets the true beta, must be hand entered
  */
  
  fvec beta(p0*size);
  // Initialize beta to be 0
  for (i = 0; i < p0*size; i++) {
    beta(i) = 0;
  }

  fvec sdbeta(sdp0*size);

  for (i = 0; i < sdp0*size; i++) {
    sdbeta(i) = 0;
  }
  // Enter everything else
  if (size == 1) {
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

  } else if(size == 2) {
    beta(0) = 2;
    beta(1) = 2;

    beta(2) = 2;
    beta(3) = 2;

    beta(4) = 2;
    beta(5) = 2;

    beta(6) = 2;
    beta(7) = 2;

    sdbeta(0) = 2;
    sdbeta(1) = 2;

    sdbeta(2) = 2;
    sdbeta(3) = 2;

    sdbeta(4) = 2;
    sdbeta(5) = 2;

    sdbeta(6) = 2;
    sdbeta(7) = 2;
//    beta(8) = 3;
//    beta(9) = 3;
//
//    beta(10) = 3;
//    beta(11) = 3;

    tbeta(0) = 2;
    tbeta(1) = 2;

    tbeta(2) = 2;
    tbeta(3) = 2;

    tbeta(4) = 2;
    tbeta(5) = 2;

    tbeta(6) = 2;
    tbeta(7) = 2;

//    tbeta(8) = 3;
//    tbeta(9) = 3;
//
//    tbeta(10) = 3;
//    tbeta(11) = 3;

  } else if(size == 3) {
    beta(0) = 2;
    beta(1) = 2;
    beta(2) = 2;

    beta(3) = 2;
    beta(4) = 2;
    beta(5) = 2;

    beta(6) = 2;
    beta(7) = 2;
    beta(8) = 2;

    beta(9) = 2;
    beta(10) = 2;
    beta(11) = 2;

    sdbeta(0) = 2;
    sdbeta(1) = 2;
    sdbeta(2) = 2;

    sdbeta(3) = 2;
    sdbeta(4) = 2;
    sdbeta(5) = 2;

    sdbeta(6) = 2;
    sdbeta(7) = 2;
    sdbeta(8) = 2;

    sdbeta(9) = 2;
    sdbeta(10) = 2;
    sdbeta(11) = 2;
//    beta(12) = 3;
//    beta(13) = 3;
//    beta(14) = 3;
//
//    beta(15) = 3;
//    beta(16) = 3;
//    beta(17) = 3;

    tbeta(0) = 2;
    tbeta(1) = 2;
    tbeta(2) = 2;

    tbeta(3) = 2;
    tbeta(4) = 2;
    tbeta(5) = 2;

    tbeta(6) = 2;
    tbeta(7) = 2;
    tbeta(8) = 2;

    tbeta(9) = 2;
    tbeta(10) = 2;
    tbeta(11) = 2;

//    tbeta(12) = 3;
//    tbeta(13) = 3;
//    tbeta(14) = 3;
//
//    tbeta(15) = 3;
//    tbeta(16) = 3;
//    tbeta(17) = 3;

  } else if(size == 4) {
    beta(0) = 2;
    beta(1) = 2;
    beta(2) = 2;
    beta(3) = 2;

    beta(4) = 2;
    beta(5) = 2;
    beta(6) = 2;
    beta(7) = 2;

    beta(8) = 2;
    beta(9) = 2;
    beta(10) = 2;
    beta(11) = 2;

    beta(12) = 2;
    beta(13) = 2;
    beta(14) = 2;
    beta(15) = 2;

    sdbeta(0) = 2;
    sdbeta(1) = 2;
    sdbeta(2) = 2;
    sdbeta(3) = 2;

    sdbeta(4) = 2;
    sdbeta(5) = 2;
    sdbeta(6) = 2;
    sdbeta(7) = 2;

    sdbeta(8) = 2;
    sdbeta(9) = 2;
    sdbeta(10) = 2;
    sdbeta(11) = 2;

    sdbeta(12) = 2;
    sdbeta(13) = 2;
    sdbeta(14) = 2;
    sdbeta(15) = 2;
//    beta(16) = 3;
//    beta(17) = 3;
//    beta(18) = 3;
//    beta(19) = 3;
//
//    beta(20) = 3;
//    beta(21) = 3;
//    beta(22) = 3;
//    beta(23) = 3;

    tbeta(0) = 2;
    tbeta(1) = 2;
    tbeta(2) = 2;
    tbeta(3) = 2;

    tbeta(4) = 2;
    tbeta(5) = 2;
    tbeta(6) = 2;
    tbeta(7) = 2;

    tbeta(8) = 2;
    tbeta(9) = 2;
    tbeta(10) = 2;
    tbeta(11) = 2;

    tbeta(12) = 2;
    tbeta(13) = 2;
    tbeta(14) = 2;
    tbeta(15) = 2;

//    tbeta(16) = 3;
//    tbeta(17) = 3;
//    tbeta(18) = 3;
//    tbeta(19) = 3;
//
//    tbeta(20) = 3;
//    tbeta(21) = 3;
//    tbeta(22) = 3;
//    tbeta(23) = 3;

  }

  int iNum = nCr(tNgrp, 2);

  fvec ibeta(iNum*iSize);

  for (i = 0; i < iNum*iSize; i++) {
    ibeta(i) = 0;
  }

  if (iSize == 1) {
    ibeta(0) = 2;
    ibeta(1) = 2;
  } else if(iSize == 4) {
    sub[0] = 0;
    sub[1] = 1;
    sub[2] = 2;
    sub[3] = 3;
    
    sub[4] = 12;
    sub[5] = 13;
    sub[6] = 14;
    sub[7] = 15;
    
    sub[8] = 20;
    sub[9] = 21;
    sub[10] = 22;
    sub[11] = 23;

    ibeta(sub[0]) = 3;
    ibeta(sub[1]) = 3;
    ibeta(sub[2]) = 3;
    ibeta(sub[3]) = 3;

    ibeta(sub[4]) = 3;
    ibeta(sub[5]) = 3;
    ibeta(sub[6]) = 3;
    ibeta(sub[7]) = 3;

    ibeta(sub[8]) = 3;
    ibeta(sub[9]) = 3;
    ibeta(sub[10]) = 3;
    ibeta(sub[11]) = 3;

    tbeta(8) = 3;
    tbeta(9) = 3;
    tbeta(10) = 3;
    tbeta(11) = 3;

    tbeta(12) = 3;
    tbeta(13) = 3;
    tbeta(14) = 3;
    tbeta(15) = 3;

    tbeta(16) = 3;
    tbeta(17) = 3;
    tbeta(18) = 3;
    tbeta(19) = 3;
  } else if(iSize == 9) {
    
    sub[0] = 0;
    sub[1] = 1;
    sub[2] = 2;
    sub[3] = 3;
    sub[4] = 4;
    sub[5] = 5;
    sub[6] = 6;
    sub[7] = 7;
    sub[8] = 8;
    
    sub[9] = 27;
    sub[10] = 28;
    sub[11] = 29;
    sub[12] = 30;
    sub[13] = 31;
    sub[14] = 32;
    sub[15] = 33;
    sub[16] = 34;
    sub[17] = 35;
    
    sub[18] = 45;
    sub[19] = 46;
    sub[20] = 47;
    sub[21] = 48;
    sub[22] = 49;
    sub[23] = 50;
    sub[24] = 51;
    sub[25] = 52;
    sub[26] = 53;

    ibeta(sub[0]) = 3;
    ibeta(sub[1]) = 3;
    ibeta(sub[2]) = 3;
    ibeta(sub[3]) = 3;
    ibeta(sub[4]) = 3;
    ibeta(sub[5]) = 3;
    ibeta(sub[6]) = 3;
    ibeta(sub[7]) = 3;
    ibeta(sub[8]) = 3;

    ibeta(sub[9]) = 3;
    ibeta(sub[10]) = 3;
    ibeta(sub[11]) = 3;
    ibeta(sub[12]) = 3;
    ibeta(sub[13]) = 3;
    ibeta(sub[14]) = 3;
    ibeta(sub[15]) = 3;
    ibeta(sub[16]) = 3;
    ibeta(sub[17]) = 3;

    ibeta(sub[18]) = 3;
    ibeta(sub[19]) = 3;
    ibeta(sub[20]) = 3;
    ibeta(sub[21]) = 3;
    ibeta(sub[22]) = 3;
    ibeta(sub[23]) = 3;
    ibeta(sub[24]) = 3;
    ibeta(sub[25]) = 3;
    ibeta(sub[26]) = 3;

    tbeta(12) = 3;
    tbeta(13) = 3;
    tbeta(14) = 3;
    tbeta(15) = 3;
    tbeta(16) = 3;
    tbeta(17) = 3;
    tbeta(18) = 3;
    tbeta(19) = 3;
    tbeta(20) = 3;

    tbeta(21) = 3;
    tbeta(22) = 3;
    tbeta(23) = 3;
    tbeta(24) = 3;
    tbeta(25) = 3;
    tbeta(26) = 3;
    tbeta(27) = 3;
    tbeta(28) = 3;
    tbeta(29) = 3;

    tbeta(30) = 3;
    tbeta(31) = 3;
    tbeta(32) = 3;
    tbeta(33) = 3;
    tbeta(34) = 3;
    tbeta(35) = 3;
    tbeta(36) = 3;
    tbeta(37) = 3;
    tbeta(38) = 3;
  } else if(iSize == 16) {
    sub[0] = 0;
    sub[1] = 1;
    sub[2] = 2;
    sub[3] = 3;
    sub[4] = 4;
    sub[5] = 5;
    sub[6] = 6;
    sub[7] = 7;
    sub[8] = 8;
    sub[9] = 9;
    sub[10] = 10;
    sub[11] = 11;
    sub[12] = 12;
    sub[13] = 13;
    sub[14] = 14;
    sub[15] = 15;
    
    sub[16] = 48;
    sub[17] = 49;
    sub[18] = 50;
    sub[19] = 51;
    sub[20] = 52;
    sub[21] = 53;
    sub[22] = 54;
    sub[23] = 55;
    sub[24] = 56;
    sub[25] = 57;
    sub[26] = 58;
    sub[27] = 59;
    sub[28] = 60;
    sub[29] = 61;
    sub[30] = 62;
    sub[31] = 63;
    
    sub[32] = 80;
    sub[33] = 81;
    sub[34] = 82;
    sub[35] = 83;
    sub[36] = 84;
    sub[37] = 85;
    sub[38] = 86;
    sub[39] = 87;
    sub[40] = 88;
    sub[41] = 89;
    sub[42] = 90;
    sub[43] = 91;
    sub[44] = 92;
    sub[45] = 93;
    sub[46] = 94;
    sub[47] = 95;

    ibeta(sub[0]) = 3;
    ibeta(sub[1]) = 3;
    ibeta(sub[2]) = 3;
    ibeta(sub[3]) = 3;
    ibeta(sub[4]) = 3;
    ibeta(sub[5]) = 3;
    ibeta(sub[6]) = 3;
    ibeta(sub[7]) = 3;
    ibeta(sub[8]) = 3;
    ibeta(sub[9]) = 3;
    ibeta(sub[10]) = 3;
    ibeta(sub[11]) = 3;
    ibeta(sub[12]) = 3;
    ibeta(sub[13]) = 3;
    ibeta(sub[14]) = 3;
    ibeta(sub[15]) = 3;

    ibeta(sub[16]) = 3;
    ibeta(sub[17]) = 3;
    ibeta(sub[18]) = 3;
    ibeta(sub[19]) = 3;
    ibeta(sub[20]) = 3;
    ibeta(sub[21]) = 3;
    ibeta(sub[22]) = 3;
    ibeta(sub[23]) = 3;
    ibeta(sub[24]) = 3;
    ibeta(sub[25]) = 3;
    ibeta(sub[26]) = 3;
    ibeta(sub[27]) = 3;
    ibeta(sub[28]) = 3;
    ibeta(sub[29]) = 3;
    ibeta(sub[30]) = 3;
    ibeta(sub[31]) = 3;

    ibeta(sub[32]) = 3;
    ibeta(sub[33]) = 3;
    ibeta(sub[34]) = 3;
    ibeta(sub[35]) = 3;
    ibeta(sub[36]) = 3;
    ibeta(sub[37]) = 3;
    ibeta(sub[38]) = 3;
    ibeta(sub[39]) = 3;
    ibeta(sub[40]) = 3;
    ibeta(sub[41]) = 3;
    ibeta(sub[42]) = 3;
    ibeta(sub[43]) = 3;
    ibeta(sub[44]) = 3;
    ibeta(sub[45]) = 3;
    ibeta(sub[46]) = 3;
    ibeta(sub[47]) = 3;

    tbeta(16) = 3;
    tbeta(17) = 3;
    tbeta(18) = 3;
    tbeta(19) = 3;
    tbeta(20) = 3;
    tbeta(21) = 3;
    tbeta(22) = 3;
    tbeta(23) = 3;
    tbeta(24) = 3;
    tbeta(25) = 3;
    tbeta(26) = 3;
    tbeta(27) = 3;
    tbeta(28) = 3;
    tbeta(29) = 3;
    tbeta(30) = 3;
    tbeta(31) = 3;

    tbeta(32) = 3;
    tbeta(33) = 3;
    tbeta(34) = 3;
    tbeta(35) = 3;
    tbeta(36) = 3;
    tbeta(37) = 3;
    tbeta(38) = 3;
    tbeta(39) = 3;
    tbeta(40) = 3;
    tbeta(41) = 3;
    tbeta(42) = 3;
    tbeta(43) = 3;
    tbeta(44) = 3;
    tbeta(45) = 3;
    tbeta(46) = 3;
    tbeta(47) = 3;

    tbeta(48) = 3;
    tbeta(49) = 3;
    tbeta(50) = 3;
    tbeta(51) = 3;
    tbeta(52) = 3;
    tbeta(53) = 3;
    tbeta(54) = 3;
    tbeta(55) = 3;
    tbeta(56) = 3;
    tbeta(57) = 3;
    tbeta(58) = 3;
    tbeta(59) = 3;
    tbeta(60) = 3;
    tbeta(61) = 3;
    tbeta(62) = 3;
    tbeta(63) = 3;
  }

  if (oracle == 0) {
//print_matrix(beta); 


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
//  Only run this the first time (Finding true_var)

    fmat sdd(sdn, sdp0);

//  sdd = rmvn_Choleski(sdn, mu, g_sigma, p0);
    sdd = rmvn(sdn, sdp0, rho);

    for ( k = 0; k < sdn; k++) {
      for ( i = 0; i < sdp0; i++) {
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

    fmat sdX(sdn, sdp0*(size));

    for ( i = 0; i < sdn; i++) {
     int cnt = 0;
      for ( j = 0; j < sdp0*(size); j+=size){
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
 
//  fmat sdiX(sdn, iNum*iSize);
//
//  int cnt=0;
//  for (i = 0; i < tNgrp*size; i+=size) {
//    for (j = i; j < tNgrp*size; j+=size) {
//      if (i != j) {
//        for (int a = 0; a < size; a++) {
//          for (int b = 0; b < size; b++) {
//            sdiX.col(cnt) = sdX.col(i+a) % sdX.col(j+b);
//            cnt++;
//          }
//        }
//      }
//    }
//  }

    fvec sdY(sdn);

//  sdY = sdX * sdbeta + sdiX * ibeta;
    sdY = sdX * sdbeta;

    float samp_var = var(sdY);

    true_var = samp_var / SNR;

    float true_sigma = sqrt(true_var);

//  printf("Sample Variance (n = %d): %lf\n", sdn, samp_var);
//  printf("SNR: %f\n", SNR);
//  printf("True Variance: %lf\tTrue Sigma: %lf\n", true_var, true_sigma);
/*
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
 
    fmat sdiX(sdn, iNum*iSize);

    int cnt=0;
    for (i = 0; i < tNgrp*size; i+=size) {
      for (j = i; j < tNgrp*size; j+=size) {
        if (i != j) {
          for (int a = 0; a < size; a++) {
            for (int b = 0; b < size; b++) {
              sdiX.col(cnt) = sdX.col(i+a) % sdX.col(j+b);
              cnt++;
            }
          }
        }
      }
    }

    fvec LsdY(sdn);

    LsdY = sdX * beta + sdiX * ibeta;

    float Lsamp_var = var(LsdY);

    true_var = Lsamp_var / SNR;

    printf("Sample Variance (po = %d): %f\n", p0, Lsamp_var);
*/
/*
    if (samp_var == -1) {
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

//    samp = rmvn_Choleski(n, mu, g_sigma, p0);
    samp = rmvn(n, p0, rho);
//printf("samp (before):\n");
//print_matrix(samp);

//    printf("Here 1\n");
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
//printf("samp (after):\n");
//print_matrix(samp);

//    printf("Here 2\n");
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

    samp.reset();
    fmat iX(n, iNum*iSize);

    int cnt=0;
    for (i = 0; i < tNgrp*size; i+=size) {
      for (j = i; j < tNgrp*size; j+=size) {
        if (i != j) {
          for (int a = 0; a < size; a++) {
            for (int b = 0; b < size; b++) {
              iX.col(cnt) = X.col(i+a) % X.col(j+b);
              cnt++;
            }
          }
        }
      }
    }

    epsilon = rnorm(0, true_sigma, n);

//    printf("Here 4\n");

    y = X * beta + iX * ibeta + epsilon;

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
    Namey << SNR;
    Namey << "_";
    Namey << var_ind;
  
    ofp.open(Namey.str().c_str(), ios::out);
    ofp << y;
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
    NameX << SNR;
    NameX << "_";
    NameX << var_ind;
  
    ofp.open(NameX.str().c_str(), ios::out);
    ofp << X;
    ofp.close();

  } else {

    for (int c = 0; c < iter; c++) {

      fmat sdd(sdn, sdp0);
  
      sdd = rmvn(sdn, sdp0, rho);
  
      for ( k = 0; k < sdn; k++) {
        for ( i = 0; i < sdp0; i++) {
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
  
      fmat sdX(sdn, sdp0*(size));
  
      for ( i = 0; i < sdn; i++) {
       int cnt = 0;
        for ( j = 0; j < sdp0*(size); j+=size){
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
   
//    fmat sdiX(sdn, iNum*iSize);
//
//    int cnt=0;
//    for (i = 0; i < tNgrp*size; i+=size) {
//      for (j = i; j < tNgrp*size; j+=size) {
//        if (i != j) {
//          for (int a = 0; a < size; a++) {
//            for (int b = 0; b < size; b++) {
//              sdiX.col(cnt) = sdX.col(i+a) % sdX.col(j+b);
//              cnt++;
//            }
//          }
//        }
//      }
//    }
  
      fvec sdY(sdn);
  
//    sdY = sdX * sdbeta + sdiX * ibeta;
      sdY = sdX * sdbeta;
  
      float samp_var = var(sdY);
  
      true_var = samp_var / SNR;
  
      float true_sigma = sqrt(true_var);

      sdY.reset();
      sdd.reset();
      sdX.reset();
//    sdiX.reset(); 

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

      fmat iX(n, iNum*iSize);

      int cnt=0;
      for (i = 0; i < tNgrp*size; i+=size) {
        for (j = i; j < tNgrp*size; j+=size) {
          if (i != j) {
            for (int a = 0; a < size; a++) {
              for (int b = 0; b < size; b++) {
                iX.col(cnt) = X.col(i+a) % X.col(j+b);
                cnt++;
              }
            }
          }
        }
      }

      epsilon = rnorm(0, true_sigma, n);
    
      y = X * beta + iX * ibeta + epsilon;
//    y = X * beta + epsilon;

      epsilon.reset();

      fmat oX(n, size*tNgrp + 3*iSize); 
//    fmat oX(n, size*tNgrp); 

      for (i = 0; i < size*tNgrp; i++) {
        oX.col(i) = X.col(i);
      }

      for (i = 0; i < 3*iSize; i++) {
        oX.col(size*tNgrp + i) = iX.col(sub[i]);
      }
      
      X.reset();
      iX.reset();

      hbeta = coef(y, oX, tNgrp*size + 3*iSize);
//printf("hbeta: \n");
//print_matrix(hbeta);
//      hbeta = coef(y, oX, tNgrp*size);

      if (sum(hbeta) == 0) {
        c = c - 1;
      } else {
        calc(c) = MSE(hbeta, tbeta, tNgrp*size + 3*iSize);
//        calc(c) = MSE(hbeta, tbeta, tNgrp*size);
      }
//printf("calc(%d) = %lf\n", c, calc(c));
      oX.reset();
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

   for (int mbic = 0; mbic <= 1; mbic++) {
     for (int ibic = 0; ibic <= 1; ibic++) {

        ofstream ofp;
        stringstream nam;

        nam << "summary_result.";
        nam << p0;
        nam << ".";
        nam << n;
        nam << ".";
        nam << n_grp;
        nam << ".";
        nam << iSNR;
        nam << ".";
        nam << mbic;
        nam << ".";
        nam << ibic;
        nam << ".";
        nam << var_ind;
        nam << ".out";

        ofp.open(nam.str().c_str(), ios::out);

        if (mbic == 1){
          ofp << "Main BIC = Small BIC calculation" << endl;
        } else {
          ofp << "Main BIC = Large BIC calculation" << endl;
        }
        if (ibic == 1){
          ofp << "Interaction BIC = Small BIC calculation" << endl;
        } else {
          ofp << "Interaction BIC = Large BIC calculation" << endl;
        }
        ofp << "Total number of main Groups = " << p0 << endl;
        ofp << "Total number of main Parameters = " << p0*size << endl;
        ofp << "Main Group Size = " << n_grp << endl;
        ofp << "Sample Size = " << n << endl;
        ofp << "Iterations = " << iter << endl;
        if (var_ind == 3){
          ofp << "Variance Structure = AOR" << endl;
        } else {
          ofp << "Variance Structure = IOR" << endl;
        }
        ofp << "SNR = " << cSNR << endl;
        ofp << "Coverage probability = 1.00 (0.00)" << endl;
        ofp << "Percentage of correct zeros = 1.00 (0.00)" << endl;
        ofp << "Percentage of incorrect zeros = 0.00 (0.00)" << endl;
        ofp << "Exact select probability = 1.00 (0.00)" << endl;
        ofp << "Model Size = " << tNgrp << " (0.00)" << endl;
        ofp << "Interaction Coverage probability = 1.00 (0.00)" << endl;
        ofp << "Interaction Percentage of correct zeros = 1.00 (0.00)" << endl;
        ofp << "Interaction Percentage of incorrect zeros = 0.00 (0.00)" << endl;
        ofp << "Interaction Exact select probability = 1.00 (0.00)" << endl;
        ofp << "Interaction Model Size = 3.00 (0.00)" << endl;
        ofp << "MSE = " << m << " (" << s << ")" << endl;
        ofp << "Total time = 0.00 (0.00)" << endl;

        ofp.close();

     }
   }

//   printf("p0: %d\tn: %d\tngrp: %d\tsigma: %d\tvar: %d\tmean MSE: %lf\tsd MSE: %lf\n", p0, n, (int) n_grp, sigma, m, s);
  }
}
