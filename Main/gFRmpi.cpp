/*=============================================================================
 |      Program:  gFRmpi.cpp
 |
 |       Author:  K. Michels
 |     Language:  C++ (mpic++ compiler)
 |   To Compile:  mpic++ -g -Wall -larmadillo gFRmpi.cpp -o gFRmpi (t is assumed 
 |                that the stack header file, header.h, is in the same directory.)      
 |
 +-----------------------------------------------------------------------------
 |
 |  Description:  Given the correct input files, this program does my algorithm
 |                of group forward regression to find the significant groups in
 |                in the input file.
 |
 |        Input:  There are two input options.  The first is PLINK format, 
 |                (assume no missing data).  PLINK has two data
 |                files as input.  First the *.ped, explained at this website:
 |                http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
 |                The second being *.map, explained at this website:
 |                http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
 |                The following inputs are required for PLINK format:
 |
 |                argv[1] -- PM indicates PLINK format TRUE
 |                argv[2] -- .ped file name
 |                argv[3] -- .map file name
 |                Need tags, either 0 or 1; FALSE or TRUE
 |                argv[4] -- Family ID included
 |                argv[5] -- Parents included
 |                argv[6] -- Sex included
 |
 |                The second input option is simply termed standard format.  
 |                This format includes two data files, the response (referred
 |                to as y) and the explanatory (referred to as X).  The X, is a 
 |                design matrix where the continuous covariates are in the first
 |                'numc' columns.  Then each column after that is non-continuous,
 |                or factors, each column is a group.  The following inputs are 
 |                required for standard input:
 |
 |                argv[1] -- STD indicatins STANDARD format TRUE
 |                argv[2] -- response (y) file name
 |                argv[3] -- explanatory (X) file name.
 |                Tags, indicating sizes
 |                argv[4] -- 'numc', number of covariates included
 |                argv[5] -- '.csv output', (either yes (1) or no)
 |                argv[6] -- 'var', variance structure of data // if .csv true
 |                argv[7] -- 'sigma', error of epsilon         //   then include
 |                                                             //   both of these
 |                                                             //   parameters
 |
 |       Output:  The output produced includes 3 lists of significant groups.  The 
 |                first is ALL groups found after the first step.  The second is
 |                significant groups after large BIC on second step.  And the third
 |                is significant groups after small BIC on second step.  If NOT 
 |                PLINK format then the beta.hat estimates are also calculated and 
 |                reported.
 |
 |    Algorithm:  The grouping forward regression algorithm is used.  It is 
 |                detailed out in the included paper.
 |
 |   Known Bugs:  None; all operations work as they should.
 |
 *===========================================================================*/

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "header.h"

using namespace std;
using namespace arma;

int main( int argc, char *argv[] )
{

  int i, j, k; // iterators
  wall_clock timer; // times program

  double t;
  // The send and receive buffers for mpi
  double *sendb = (double *)malloc(4*sizeof(double));
  double *recvb = (double *)malloc(4*sizeof(double));
  int np; // Total number of nodes
  int me; // Current node

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  // Every node reads in the data in parallel //
  string type = argv[1];
  int n = 0; // Sample size
  int p0 = 0; // Number of parameters in X batrix
  int numc = 0; // Number of continuous covariates
  int csv = 0; // Two csv choices for output
  int tgrp = 0; // Total number of groups

  vector<string> FiD; // Family ID
  int cntFiD = 0; // Number of members in family

  vector<string> PiD; // Paternal ID
  int cntPiD = 0; // Number of fathers

  vector<string> MiD; // Maternal ID
  int cntMiD = 0; // Number of mothers

  vector<int> SiD; // Sex
  int cntSiD = 0;

  int nuc = 0; // If 1 then A, G, C, T; if 0 then 1, 2, 3, 4
  vector<string> strNUC;
  vector<int> zeroNUC;
  vector<int> iNUC;



  if (type.compare("PM") == 0){

    int fam = atoi(argv[4]);
    int par = atoi(argv[5]);
    int sex = atoi(argv[6]);

    // Get the sample size
    char *ped = argv[2];
    ifstream fp;
    fp.open(ped, ios::in);
    string l;

    getline(fp, l);
    stringstream pheadline; pheadline<<l;
    vector<string> Plist;

    int cmpFiD = 0;
    int cmpPiD = 0;
    int cmpMiD = 0;

    // Look through first 7 columns of .ped file
    for(int i=0; i<7; i++)
    {
      pheadline>>l;
      Plist.push_back(l);
      if ((i == 0) && (fam == 1)) { // The first column is Family ID
        FiD.push_back(Plist.at(i));
      } else if ((i == 2) && (par == 1)) { // The third column is Paternal ID
        PiD.push_back(Plist.at(i));
      } else if ((i == 3) && (par == 1)) { // The fourth column is Maternal ID
        MiD.push_back(Plist.at(i));
      } else if (i == 6) { // The 7th column is the SNP, this determines what type SNP is, this isn't failproof though
        while (l.compare("0") == 0) {
          pheadline>>l;
        }
        if (l.compare("A") || l.compare("G") || l.compare("C") || l.compare("T")) {
          nuc = 1;
        }
      }
    }

    // Because the ID's are factors, first determines all possible values for the
    //  various ID's.
    while(!fp.eof()) {
      getline(fp, l);
      stringstream line; line<<l;
      vector<string> Llist;

      for(int i=0; i<6; i++)
      {
        line>>l;
        Llist.push_back(l);
        if ((i == 0) && (fam == 1)) {
          if (std::find(FiD.begin(), FiD.end(), Llist.at(i)) != FiD.end())
          {
            cmpFiD = 1;
          }
          if (cmpFiD == 0) {
            FiD.push_back(Llist.at(i));
            cntFiD++;
          }
          cmpFiD = 0;
        } else if ((i == 2) && (par == 1)) {
          if (std::find(PiD.begin(), PiD.end(), Llist.at(i)) != PiD.end())
          {
            cmpPiD = 1;
          }
          if (cmpPiD == 0) {
            PiD.push_back(Llist.at(i));
            cntPiD++;
          }
          cmpPiD = 0;
        } else if ((i == 3) && (par == 1)) {
          if (std::find(MiD.begin(), MiD.end(), Llist.at(i)) != MiD.end())
          {
            cmpMiD = 1;
          }
          if (cmpMiD == 0) {
            MiD.push_back(Llist.at(i));
            cntMiD++;
          }
          cmpFiD = 0;
        }
      }
      n++;
    }

    if (fam == 1) {
      cntFiD--;
    } else if (par == 1) {
      cntPiD--;
      cntMiD--;
    } else if (sex == 1) {
      cntSiD++;
    }

    fp.close();

    // Get the number of SNPs
    char *map = argv[3];

    fp.open(map, ios::in);
    tgrp=-1; // Number of SNPs

    while(!fp.eof()) {
      getline(fp, l);
      tgrp++;
    }

    fp.close();

    numc = cntFiD + cntPiD + cntMiD + cntSiD; // For PLINK, number of covariates
    p0 = numc + p0;
  } else if (type.compare("STD") == 0){

    // Get the sample size
    char *response = argv[2];
    ifstream fp;
    fp.open(response, ios::in);
    string l;

    while(!fp.eof()){
      getline(fp, l);
      n++;
    }
    n--;
    fp.close();

    // Get the 'X' data set configured
    char *fixed = argv[3]; // First get number of fixed parameters
    fp.open(fixed, ios::in);

    getline(fp, l);
    stringstream headerline; headerline<<l;
    vector<string> Xlist;

    while(!headerline.eof())
    {
        headerline>>l;
        Xlist.push_back(l);
    }

    tgrp = Xlist.size(); // Number of fixed parameters
    numc = atoi(argv[4]);
    tgrp = tgrp - numc;
//printf("tgrp=%d\n", tgrp);
    csv = atoi(argv[5]);
  }

  // Stores group members for each group (or SNP)
  vector< vector<string> > groups(tgrp, vector<string>(n));
  fvec nINgrp(tgrp); // Number of members in group - 1
  fvec nINgrp_org(tgrp); // Same as above only isn't altered
  fvec cNum(tgrp); // Specifies matrix column in which group starts
  fvec cNum_org(tgrp); // Same as above only isn't altered

  if (type.compare("PM") == 0) {
    char *ped = argv[2];
    ifstream fp;
    fp.open(ped, ios::in);
    string l;

    for (int j=0; j<n; j++) { // Go through by individual
      getline(fp, l); // Get full line
      stringstream line; line << l;
      int iter=0;
      for(int i=0; i<2*tgrp+6; i++) // Remember must be every other column
      {
        line>>l; // Pick of member
        int cnt=0;
        int tmp=0;
        if ((i > 5) && (i%2 == 0)) {
          if (j == 0) { // If first individual, then first time.  Put member in group
            groups[iter][cnt] = l;
            cnt++; // Increment count
            nINgrp(iter)=cnt; // Store how many members
            nINgrp_org(iter)=cnt; // Store how many members
          } else {
            cnt=nINgrp(iter); // Initialize count to current groups count
            cnt=nINgrp_org(iter); // Initialize count to current groups count
            if (std::find(groups[iter].begin(), groups[iter].end(), l) != groups[iter].end())
            { // Check to see if 'l' is a new member, if 'l' is in group then tmp=1
              tmp = 1;
            }
            if (tmp == 0) { // If 'l' is not in group
              groups[iter][cnt]=l; // Add 'l' to group
              cnt++; // Increment count
              nINgrp(iter)=cnt; // Update group size
              nINgrp_org(iter)=cnt; // Update group size
            }
          }
          iter++; // Iterate
        }
      }
    }
    fp.close();
    int tp0 = 0; // Need to initialize a few vectors 
    for (int i=0; i<tgrp; i++) {
      nINgrp(i) = nINgrp(i) - 1; // Group size is size - 1
      nINgrp_org(i) = nINgrp_org(i) - 1; // Group size is size - 1
      cNum(i) = tp0;
      cNum_org(i) = tp0;
      tp0 = nINgrp(i) + tp0;
    }

//if (me == 0) {
//  for (int i = 0; i < p0; i++) {
//    printf("nINgrp[%d] = %d\n", i, (int) nINgrp(i));
//  }
//  for (int i = 0; i < p0; i++) {
//    printf("cNum[%d] = %d\n", i, (int) cNum(i));
//  }
//  printf("tp0 = %d\n", tp0);
//}

    p0 = tp0 + numc;

  } else if (type.compare("STD") == 0){
    char *fixed = argv[3];
    ifstream fp;
    fp.open(fixed, ios::in);
    string l;
    // Exact same idea is with PLINK format, I'm not going to add comments, look in
    //   PLINK section for comments.
    for (int j=0; j<n; j++) {
      getline(fp, l);
      stringstream line; line << l;
      for (i=0; i<numc; i++) {
        line>>l;
      }
      for(int i=numc; i<tgrp+numc; i++)
      {
        line>>l;
        int cnt=0;
        int tmp=0;
        if (j == 0) {
          groups[i-numc][cnt] = l;
          cnt++;
          nINgrp(i-numc)=cnt;
          nINgrp_org(i-numc)=cnt;
        } else {
          cnt=nINgrp(i-numc);
          if (std::find(groups[i-numc].begin(), groups[i-numc].end(), l) != groups[i-numc].end())
          {
            tmp = 1;
          }
          if (tmp == 0) {
            groups[i-numc][cnt]=l;
            cnt++;
            nINgrp(i-numc)=cnt;
            nINgrp_org(i-numc)=cnt;
          }
        }
      }
    }

    fp.close();
    int tp0 = 0;
    for (int i=0; i<tgrp; i++) {
      nINgrp(i) = nINgrp(i) - 1;
      nINgrp_org(i) = nINgrp_org(i) - 1;
      cNum(i) = tp0;
      cNum_org(i) = tp0;
      tp0 = nINgrp(i) + tp0;
    }

    p0 = tp0 + numc;
  }

//if (me == 0) {
//  for (i=0; i<tgrp; i++) {
//    for (j=0; j<nINgrp(i); j++) {
//      cout << "groups[" << i << "][" << j << "]: " << groups[i][j] << endl;
//    }
//  }
//}

  //  Create the X 'matrix'
  fmat X(n, p0);

  //  Create the y 'matrix'
  fvec y(n);

  if (type.compare("PM") == 0){

    int fam = atoi(argv[4]);
    int par = atoi(argv[5]);
    int sex = atoi(argv[6]);

    char *ped = argv[2];
    ifstream fp;
    fp.open(ped, ios::in);
    string l;

    // Go through .ped, add to X file
    for (int j=0; j < n; j++) {

      getline(fp, l);
      stringstream line; line<<l;
      string s;
      double d;
      int a;
      int cnt = 0;

      line>>s; // Family ID   
      if (fam == 1) {
        for (int i = 0; i < cntFiD; i++) {
          if (FiD[i].compare(s) == 0) {
            X(j, cnt++) = 1;
          } else {
            X(j, cnt++) = 0;
          }
        }
      }

      line>>s; // Individual ID

      line>>s; // Paternal ID
      if (par == 1) {
        for (int i = 0; i < cntPiD; i++) {
          if (PiD[i].compare(s) == 0) {
            X(j, cnt++) = 1;
          } else {
            X(j, cnt++) = 0;
          }
        }
      }

      line>>s; // Maternal ID 
      if (par == 1) {
        for (int i = 0; i < cntMiD; i++) {
          if (MiD[i].compare(s) == 0) {
            X(j, cnt++) = 1;
          } else {
            X(j, cnt++) = 0;
          }
        }
      }

      line>>a; // Sex
      if (sex == 1) {
        X(j, cnt++) = a;
      }

      line>>d; // Phenotype
      y(j) = d;

      // This is a very simplistic setup, and it doesn't account for missing data
      //   if there is missing data it is simple to deal with that contingency 
      //   using other means
      int iter=0;
      // Go through all SNPs for an individual
      for ( i = 6; i < 2*tgrp + 6; i++){
        line>>s;
        // Since PLINK format has two columns for every SNP
        if (i % 2 == 0) {
          // Go through members of SNP
          for ( k = 0; k < nINgrp(iter); k++) {
            // If the comparison is equal 1 o/w 0
            if (groups[iter][k].compare(s) == 0) {
              X(j, cnt++) = 1;
            } else {
              X(j, cnt++) = 0;
            }
          }
          iter++; // Iterate
        }
      }

//  Below is for if there is missing data in the PED file, the following can
//    be used.  It was used for the main effects, but it has not been tested
//    sufficiently

//        if (nuc == 1) {
//          line>>s;
//        } else {
//          line>>a;
//        }
//        if (i % 2 == 0) {
//          if (nuc == 1) { // This means A, G, C and T
//            // If '0' is first then make a note of it.  '0' indicates a missing
//            //   value.  So '0' (clearly) is not A, G, C or T
//            if ((j == 0) && (s.compare("0") == 0)) {
//              zeroNUC.insert(zeroNUC.begin() + cntNUC, 1);
//            } else if (j == 0) {
//              zeroNUC.insert(zeroNUC.begin() + cntNUC, 0);
//            }
//            // This value is the comparison value.  The column is biallelic,
//            //   meaning it is this value OR a different one.  If this value 
//            //   is '0', this means it is missing.  So it must be replaced.
//            if (j == 0) {
//              strNUC.push_back(s);
//            }
//            // Here if '0' is first, then count number of zeroes in the column,
//            //   store this number in zeroNUC
//            if ((zeroNUC[cntNUC] >= 1) && (s.compare("0") == 0)) {
//              zeroNUC[cntNUC]++;
//            }
//            // Here the first value was '0', now we've come to the first value
//            //   in the column which is not '0', strNUC needs to be replaced
//            //   by that value.
//            if ((zeroNUC[cntNUC] > 0) && (s.compare("0") != 0)) {
//              strNUC.insert(strNUC.begin() + cntNUC, s);
//              // What I can do here is randomly choose 1 or 0, here I just
//              //   insert '1' for every value that was missing before.
//              int tmp = zeroNUC[cntNUC] - 1;
//              for (int a = 0; a < tmp; a++) {
//                X(a, cnt) = 0;
//                zeroNUC[cntNUC]--;           
//              }
//              zeroNUC[cntNUC]--;           
//              // Now the value is no longer '0', needs to be replaced.
//            }
//
//            // The current value either is the same as the first or different.
//            if (strNUC[cntNUC++].compare(s) == 0) {
//              X(j, cnt++) = 1;
//            } else {
//              X(j, cnt++) = 0;
//            }
//
//          } else { // This means 1, 2, 3 and 4
//            // If '0' is first then make a note of it.  '0' indicates a missing
//            //   value.  So '0' (clearly) is not 1, 2, 3 or 4
//            if ((j == 0) && (a == 0)) {
//              zeroNUC.insert(zeroNUC.begin() + cntNUC, 1);
//            } else if (j == 0) {
//              zeroNUC.insert(zeroNUC.begin() + cntNUC, 0);
//            }
//            // This value is the comparison value.  The column is biallelic,
//            //   meaning it is this value OR a different one.  If this value 
//            //   is '0', this means it is missing.  So it must be replaced.
//            if (j == 0) {
//              iNUC.push_back(a);
//            }
//            // Here if '0' is first, then count number of zeroes in the column,
//            //   store this number in zeroNUC
//            if ((zeroNUC[cntNUC] >= 1) && (a == 0)) {
//              zeroNUC[cntNUC]++;
//            }
//            // Here the first value was '0', now we've come to the first value
//            //   in the column which is not '0', strNUC needs to be replaced
//            //   by that value.
//            if ((zeroNUC[cntNUC] > 0) && (a != 0)) {
//              iNUC.insert(iNUC.begin() + cntNUC, a);
//              // What I can do here is randomly choose 1 or 0, here I just
//              //   insert '1' for every value that was missing before.
//              int tmp = zeroNUC[cntNUC] - 1;
//              for (int b = 0; b < tmp; b++) {
//                X(b, cnt) = 0;
//                zeroNUC[cntNUC]--;           
//              }
//              zeroNUC[cntNUC]--;           
//              // Now the value is no longer '0', needs to be replaced.i
//            }
//            if (iNUC[cntNUC++] == a) {
//              X(j, cnt++) = 1;
//            } else {
//              X(j, cnt++) = 0;
//            }
//          }
//        }
//      }
    }
    fp.close();
  // If standard form (STD)
  } else if (type.compare("STD") == 0){

    // Read in the design matrix (X)
    char *fixed = argv[3];
    ifstream fp;
    fp.open(fixed, ios::in);
    string l;

    // Go through .ped, add to X file
    for (j=0; j < n; j++) {

      getline(fp, l);
      stringstream line; line<<l;
      string s;
      int cnt = 0;

      // This is exactly the same as PLINK format comments, if questions look
      //   at those.
      for (int i = 0; i < numc; i++){
        line>>s;
//printf("Cts. Element: %d\t%lf\n", i, atof(s.c_str()));
        X(j, cnt++) = atof(s.c_str());
      }

      int iter=0;
      // This is exactly the same as PLINK format comments, if questions look
      //   at those.
      for (int i = numc; i < tgrp + numc; i++){
        line>>s;
        for (int k = 0; k < nINgrp(iter); k++) {
          if (groups[iter][k].compare(s) == 0) {
            X(j, cnt++) = 1;
          } else {
            X(j, cnt++) = 0;
          }
        }
        iter++;
      }
    }
    fp.close();

    // Read in response vector (y)
    char *response = argv[2];
    fp.open(response, ios::in);

    // Put data into vector form
    for (j=0; j < n; j++) {
      getline(fp, l);
      stringstream line; line<<l;
      double d;
      line>>d;
      y(j) = d;
    }
    fp.close();

  }


//if (me == 0) {
//  printf("X:\n");
//  print_matrix(X);
//  printf("n=%d\tp0=%d\tnumc=%d\n", n, p0, numc);
//}

//  These are the first things added to the model.  It makes sense to me, because the data
//    obviously has this structure.  The problem is 17/36 are chosen, with some being very
//    strong effects, while when the covariates are not done SIS finds 23/36, but does not
//    find them to be really strong.
// 
  // This first adds the covariates to the model
  if (numc > 0) {
    fmat covX(n, numc);
    for ( i = 0; i < numc; i++) {
      covX.col(i) = X.col(i) / norm(X.col(i), 2);
    }
    for ( i = 0; i < numc; i++){
      X.shed_col(0);
    }

//    if (type.compare("PM") != 0) {
      y = y - covX * (covX.t() * y);
      X = X - covX * (covX.t() * X);
//    }
  }

  // The parameters left is the total minus the number of covariates
  //  When I don't need to calculate a coeffecient I can get rid of X0 and y0
  int ncov = p0 - numc;
  // This part is memory intensive for large data files.  This is used to 
  //   calculate the beta coefficients only.  Not necessary
  fmat X0(n, ncov);
  X0 = X;
  fvec y0(n);
  y0 = y;
 
  // This is for Node 0 only.  Node 0 is the main Node that controls everything
  //   it controls all sending and receiving of messages.
  if (me == 0){

    // Start the group FR //
    timer.tic();

    // 'nStep' at the end will be total number of groups selected in Step 1
    int nStep = numc; // Covariates added to model

    vec L2(n + 1); // Likelihood vector
    vec size(n + 1); // Group sizes for each group
    fvec mFR(n + 1); // SNP vector (keeps parameters)
    fvec gFR(n + 1); // Significant groups
    fvec ID(ncov); // Stores all SNP ids

    for (i=0; i<ncov; i++)
        ID(i) = i;

    int len2 = tgrp; // Number of groups

    // Collect the RSS values
    double *r = (double *)malloc(len2*sizeof(double));

    // Collect the SNP associated with the RSS
    double *subs = (double *)malloc(len2*sizeof(double));

    int cnt=0;

    double minRSS;
    int minID;

    // The following loop is the Forward Regression loop.  Each step of the 
    //   loop creates the RSS vector (r) using MPI.  Each step finds the 
    //   most significant group
    int k=0;
    // The loop goes until no more parameters can be put in (nStep < n),
    //   or when all parameters have been added (len2 > 0)
    while ((nStep < n) && ( len2 > 0 )) {

      // After the first step, because on the first step we know the Nodes's
      //   are ready, Node 0 must wait for the other Nodes to confirm
      //   they are ready before proceeding any further.
      if (k != 0) {
        for ( i = 1; i < np; i++ ){
          MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }

      // Calculates the last Likelihood.  Proven in the paper.
      L2(k) = norm(y, 2);

//double tmpL2 = L2(k);
//printf("L[%d] = %.10lf\n", k, tmpL2);

      cnt = 0;
      // If the number of groups left in X is less the the number of MPI
      //   processors
      if ( len2 < np) {
        cnt++;
      
        for ( i = 0; i < len2; i++){ // Only allocate 'len2' MPI processors
          // Prepare send buffer
          sendb[0] = (double) nINgrp(cnt - 1); // Group size
          sendb[1] = (double) 0; // Code that it's for computing RSS
          sendb[2] = (double) cNum(cnt - 1); // Column number 
          sendb[3] = (double) cnt - 1; // Group number
          MPI_Send(sendb, 4, MPI_DOUBLE, cnt++, 1, MPI_COMM_WORLD);

//if (k == 0) {
//  printf("Group Size: %d\tColumn Number: %d\tGroup Number: %d\n", (int) sendb[0], (int) sendb[2], (int) sendb[3]);
//}
        }
        // Wait to receive the RSS from the 'len2' processes
        for ( i = 0; i < len2; i++ ){
          MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          r[i] = recvb[0]; // Put the RSS in the RSS vector
          subs[i] = recvb[1]; // The corresponding group (to RSS) is recorded
        }

      // If the number of groups left is equal to or more than the number
      //   of MPI processes
      } else {
        cnt++;
        // allocate to all MPI processes
        for ( i = 1; i < np; i++){
          // Prepare send buffer
          sendb[0] = (double) nINgrp(cnt - 1);
          sendb[1] = (double) 0; // Code that it's for computing RSS
          sendb[2] = (double) cNum(cnt - 1); // column number
          sendb[3] = (double) cnt - 1; // Group number
          MPI_Send(sendb, 4, MPI_DOUBLE, cnt++, 1, MPI_COMM_WORLD);
        }

        // Wait until one of the processes is done.  When it is done save the 
        //   receiver buffer.  To save time send information from Node 0 to 
        //   calculate the next RSS.  Keep going in this manner until
        //   all RSS values are calculated.
        while (cnt <= len2) {
          // Receive message, put in buffer
          MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          // Prepare send buffer
          sendb[0] = nINgrp(cnt - 1);
          sendb[1] = (double) 0; // Code for computing RSS
          sendb[2] = (double) cNum(cnt - 1); // column number
          sendb[3] = (double) cnt - 1; // Group number

          // Send to the same MPI process that was received from (recvb[2])
          MPI_Send(sendb, 4, MPI_DOUBLE, (int) recvb[2], 1, MPI_COMM_WORLD);

          r[cnt - np] = recvb[0]; // put RSS in the RSS vector
          subs[cnt++ - np] = recvb[1]; // The corresponding group (to RSS)
        }

        // Wait to receive RSS from the last if the MPI processes
        for ( i = 1; i < np; i++ ){
          MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          r[cnt-np] = recvb[0]; // Put RSS in RSS vector
          subs[cnt++-np] = recvb[1]; // Column number
        }
      }

//if (k == 0) {
//  for (i=0; i<len2; i++) {
//    printf("r[%d] = %lf\n", i, r[i]);
//  }
//}

      // Find the minimum RSS value
      minRSS = min(r, len2);

      // And the associated subscript (corresponding group to the minimum
      //   RSS value.
      int wID = which(r, len2, minRSS);
      minID = subs[wID];
      size(k) = nINgrp(minID);
      nStep = nStep + size(k);

      // If this is the last step of the Forward Regression loop, then the
      //   MPI processes are done so send the kill signal to all nodes.
      if ((nStep >= n) || (len2 == 1)) {
        for ( i = 1; i < np; i++ ){
          sendb[1] = -1.0; // Kill code
          sendb[0] = (double) minID;
          sendb[2] = (double) k;
          MPI_Send(sendb, 4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
      // If not the last step then send the minimum RSS value AND the associated
      //   subscript (group that is to be added to the bag) to all nodes.  The 
      //   reason for this is every node has a separate X (design matrix).  
      //   Every node must delete the column (or group of columns) associated
      //   with the minimum RSS value.
      } else {
        for ( i = 1; i < np; i++ ){
          sendb[0] = (double) size(k); // Group size
          sendb[1] = 1; // Code to delete group from design matrix X
          sendb[2] = (double) cNum(minID); // Column Number
          sendb[3] = (double) minID; // Group identifier
          MPI_Send(sendb, 4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
      }

      if (nStep < n) {

        int ak = ID(minID); // Significant group number
        gFR(k) = ak; // Record significant group number
        fmat tmpX(n, size(k));
  
        for (i = 0; i < size(k); i++) {
          mFR(nStep - size(k) + i) = cNum_org(ak) + i; // Add the group of subscripts to mFR
          
          // Create matrix with columns of minimum RSS
          tmpX.col(i) = X.col(cNum(minID) + i) / norm(X.col(cNum(minID) + i), 2); 
        }
  
        // Must remove those columns from the design matrix
        for (i = 0; i < size(k); i++) {
          X.shed_col(cNum(minID));
        }
  
        ID.shed_row(minID); // Remove group of subscripts from subscript vector
        nINgrp.shed_row(minID);
        cNum.shed_row(minID);
  
        for (int h=minID; h<len2 - 1; h++) {
          cNum(h) = cNum(h) - size(k);
        }
  
        // Recalculate both y and X
        y = y - tmpX * (tmpX.t() * y);
        X = X - tmpX * (tmpX.t() * X);
        tmpX.reset();
  
        len2 = len2-1; // the number of parameters taken out of X is size
      }
      k++;
    }

    int len=k-1;

    if (nStep < n) {
      L2(k) = norm(y, 2);
      len=k;
//double tmpL2 = L2(k);
//printf("L[%d] = %.10lf\n", k, tmpL2);
    }

    // All significant groups stored.
    int *all;
    all = (int *) malloc(sizeof(int) * len);

    for (i=0; i<len; i++) {
      all[i] = (int) gFR(i);
    }

    // Calculate BIC
    double *BICl = (double *)malloc((len + 1)*sizeof(double));
    double *BICs = (double *)malloc((len + 1)*sizeof(double));

    // This vectore keeps track of model size for each model.
    fvec newvec(len+1);
    newvec(0)=numc; // Covariates added to model
    for (i=1; i<=len; i++){
      newvec(i) = size(i-1) + newvec(i-1);
    }

    // Small BIC calculation
    for ( i = 0; i <= len; i++) {
      BICs[i] = 2 * log(L2(i)) + ((newvec(i) + 1) * log(n)) / n;
    }
    // Large BIC calculation
    for ( i = 0; i <= len; i++) {
      BICl[i] = 2 * log(L2(i)) + ((newvec(i) + 1) * (log(n) + 2 * log(ncov))) / n;
    }

    // Find the minimum BIC value and the associated subscript.  This subscript
    //   is the `best' model according to BIC

    // Find the large BIC model
    double minb = min(BICl, len + 1);
    int numgl = which(BICl, len + 1, minb); // Number of groups when large BIC

    int numpl=0; // Number of parameters when large BIC
    int cNumL[numgl];
    for (i=0; i<numgl; i++){
      numpl = size(i) + numpl;
      cNumL[i] = numpl;
    }
 
    t = timer.toc();

    // Significant (large BIC) parameters
    int *anspl;
    anspl = (int *) malloc(sizeof(int) * numpl);

    // Significant (large BIC) groups
    int *ansgl;
    ansgl = (int *) malloc(sizeof(int) * numgl);

    if (numgl > 0) {
      // Takes the first `numgl' from mFR to ansp
      int tmp=0;
      for (i = 0; i < numpl; i++) {
        anspl[i] = (int) mFR(i);
        if ((i + 1) == cNumL[tmp]){
          ansgl[tmp] = (int) gFR(tmp);
          tmp++;
        }
      }
    } else {
      anspl[0] = 0;
      ansgl[0] = 0;
    }

    // Find the small BIC model
    minb = min(BICs, len + 1);
    int numgs = which(BICs, len + 1, minb); // Number of groups when small BIC

    int numps=0; // Number of parameters when small BIC
    int cNumS[numgs];
    for (i=0; i<numgs; i++){
      numps = size(i) + numps;
      cNumS[i] = numps;
    }

    // Significant (small BIC) parameters
    int *ansps;
    ansps = (int *) malloc(sizeof(int) * numps);

    // Significant (small BIC) groups
    int *ansgs;
    ansgs = (int *) malloc(sizeof(int) * numgs);

    if (numgs > 0) {
      // Takes the first `numgs' from mFR to ansp
      int tmp=0;
      for (i = 0; i < numps; i++) {
        ansps[i] = (int) mFR(i);
        if ((i + 1) == cNumS[tmp]){
          ansgs[tmp] = (int) gFR(tmp);
          tmp++;
        }
      }
    } else {
      ansps[0] = 0;
      ansgs[0] = 0;
    }

    // If PEDMAP or PLINK format.  In 'ansp' are parameter numbers.  With PLINK
    //   format the .map file contains the SNP name.  We return the SNP names
    //   with PLINK format instead of just parameter numbers.
    if (type.compare("PM") == 0){
// Do this for alll SNPs found
      // The parameters must be sorted
      quicksort(all, 0, len - 1);

      vector<string> ChromoName; // Chromosome Name
      int cntChromoName = 0;

      int total = 0; // Iteration of while loop
   
      // Open .map file
      char *map = argv[3];
      ifstream fp;
      fp.open(map, ios::in);
      string l;

      // Go through every line of the .map file.  When the iteration of the
      //   while loop `total' is equal to the parameter number, add that
      //   name to ChromoName.  Since ansp is sorted, only need to go 
      //   through .map file once. 
      while(!fp.eof()) {

        getline(fp, l);
        stringstream line; line<<l;
        string s1, s2;

        line>>s1; // Chromosome
        line>>s2; // String (2nd column, SNP identifier

        // If ChromoName isn't full
        if (cntChromoName < len) {
          // Check parameters
          if (total == (int) all[cntChromoName]) {
            ChromoName.push_back(s2);
            cntChromoName++;
          }
        // If ChromoName is full, break
        } else {
          break;
        }
        total++;
      }
      fp.close();
// Do this for large BIC
      // The parameters must be sorted
      quicksort(ansgl, 0, numgl - 1);

      vector<string> ChromoLName; // Chromosome Name
      cntChromoName = 0;

      total = 0; // Iteration of while loop
   
      // Open .map file
      fp.open(map, ios::in);

      // Go through every line of the .map file.  When the iteration of the
      //   while loop `total' is equal to the parameter number, add that
      //   name to ChromoName.  Since ansp is sorted, only need to go 
      //   through .map file once. 
      while(!fp.eof()) {

        getline(fp, l);
        stringstream line; line<<l;
        string s1, s2;

        line>>s1; // Chromosome
        line>>s2; // String (2nd column, SNP identifier

        // If ChromoName isn't full
        if (cntChromoName < numgl) {
          // Check parameters
          if (total == (int) ansgl[cntChromoName]) {
            ChromoLName.push_back(s2);
            cntChromoName++;
          }
        // If ChromoName is full, break
        } else {
          break;
        }
        total++;
      }
      fp.close();
// Do this for small BIC
      // The parameters must be sorted
      quicksort(ansgs, 0, numgs - 1);

      vector<string> ChromoSName; // Chromosome Name
      cntChromoName = 0;

      total = 0; // Iteration of while loop
   
      // Open .map file
      fp.open(map, ios::in);
      // Go through every line of the .map file.  When the iteration of the
      //   while loop `total' is equal to the parameter number, add that
      //   name to ChromoName.  Since ansp is sorted, only need to go 
      //   through .map file once. 
      while(!fp.eof()) {

        getline(fp, l);
        stringstream line; line<<l;
        string s1, s2;

        line>>s1; // Chromosome
        line>>s2; // String

        // If ChromoName isn't full
        if (cntChromoName < numgs) {
          // Check parameters
          if (total == (int) ansgs[cntChromoName]) {
            ChromoSName.push_back(s2);
            cntChromoName++;
          }
        // If ChromoName is full, break
        } else {
          break;
        }
        total++;
      }
      fp.close();

      printf("SNPs chosen (before BIC)\n");
      printf("Number of SNPs chosen: %d\n", len);
      printf("\nTotal Number of SNPs: %d\tSample size: %d\n", tgrp, n);
      printf("SNPs:\t Name\t\tNumber\n");
      for ( i = 0; i < len; i++) {
        printf("\tSNP %d: ", i);
        cout << ChromoName.at(i) << "\t" << all[i] << endl;
      }
      free(all);

      printf("\nResults for Large BIC\n");
      printf("Number of SNPs chosen: %d\n", numgl);
      printf("\nTotal Number of SNPs: %d\tSample size: %d\tLarge BIC\n", tgrp, n);
      printf("SNPs:\t Name\t\tNumber\n");
      for ( i = 0; i < numgl; i++) {
        printf("\tSNP %d: ", i);
        cout << ChromoLName.at(i) << "\t" << ansgl[i] << endl;
      }
      free(anspl);
      free(ansgl);

      printf("\nResults for Small BIC\n");
      printf("Number of SNPs chosen: %d\n", numgs);
      printf("\nTotal Number of SNPs: %d\tSample size: %d\tSmall BIC\n", tgrp, n);
      printf("SNPs:\t Name\t\tNumber\n");
      for ( i = 0; i < numgs; i++) {
        printf("\tSNP %d: ", i);
        cout << ChromoSName.at(i) << "\t" << ansgs[i] << endl;
      }
      free(ansps);
      free(ansgs);

      cout << "\nTime: " << t << " seconds" << endl;

    } else {

//      int *delp, *delg;
//      int delpLEN, delgLEN;
//
//      del(X0, ncov, n, pos, size, ansp, &delp, &delpLEN, &delg, &delgLEN);

//    Below is the correct form to calculate beta.  Go to line 958 and
//      uncomment X0 and y0.  This part is here just to keep the program
//      from throwing an error
      fvec betal = coefm(y0, X0, ncov, n, numpl, anspl);
//      fvec betal = coefm(y, X, ncov, n, numpl, anspl);

      fvec betas = coefm(y0, X0, ncov, n, numps, ansps);
//      fvec betas = coefm(y, X, ncov, n, numps, ansps);

//    printf("Delete Numbers\n");
//    for ( i = 0; i < delpLEN; i++) {
//      printf("delp[%d] = %d\tans[delp[%d]] = %.10lf\n", i, delp[i], i, ansp[delp[i]]);
//    }


//      cnt = 0;
//      int flag = 0;
//      for (i = 0; i < pos*size; i++ ){
//        for (j = 0; j < delpLEN; j++) {
//          if ((i + cnt) == delp[j]) {
//            cnt++;
//            flag = 1;
//            break;
//          }
//        }
//        if (i + cnt >= pos*size) {
//          ansp[i] = 0;
//        } else {
//          if (cnt > 0) {
//            if (flag == 1) {
//              int flag2 = 0;
//              while(flag2 ==0) {
//                int flag3 = 0;
//                for (j = 0; j < delpLEN; j++) {
//                  if ((i + cnt) == delp[j]) {
//                    flag3 = 1;
//                    break;
//                  }
//                }
//                if (flag3 == 0) {
//                  flag2 = 1;
//                } else {4
//                  cnt++;
//                }
//              }
//            }
//            ansp[i] = ansp[i + cnt];
//            flag = 0;
//          }
//        }
//      }
//
//      nump = pos*size - cnt;
//    printf("\n(Before Sort) Number of Parameters chosen: %d\n", nump);
//    printf("Parameters:\t Name\t\tNumber\n");
//    for ( i = 0; i < nump; i++) {
//      printf("\tparam %d: %.10lf\n", i, ansp[i]);
//    }

      quicksort(anspl, 0, numpl - 1);

      quicksort(ansps, 0, numps - 1);

//      cnt = 0;
//      for (i = 0; i < pos; i++ ){
//        for (j = 0; j < delgLEN; j++) {
//          if ((i + cnt) == delg[j]) {
//            cnt++;
//            break;
//          }
//        }
//        if (i + cnt >= pos) {
//          ansg[i] = 0;
//        } else {
//          if (cnt > 0) {
//            int flag = 0;
//            for (j = 0; j < delgLEN; j++) {
//              if ((i + cnt) == delg[j]) {
//                flag = 1;
//                break;
//              }
//            }
//            if (flag == 0) {
//              ansg[i] = ansg[i + cnt];
//            }
//          }
//        }
//      }
//
//      numg = pos -cnt;

      quicksort(all, 0, len - 1);
      quicksort(ansgl, 0, numgl - 1);
      quicksort(ansgs, 0, numgs - 1);

      if (csv == 0 ) {
      printf("Total number of groups: %d\tSample size: %d\n", tgrp, n);
      printf("\nNumber of groups chosen (before BIC): %d\n", len);
      printf("Groups:\n");
      for ( i = 0; i < len; i++) {
//        printf("\tSNP %d: ", i);
//        cout << ChromoName.at(i) << "\t" << all[i] << endl;
          printf("\tgroup %d: %d\n", i, all[i]);
      }
      free(all);
//        printf("Results for Large BIC\n");
//        printf("Number of Parameters chosen: %d\n", numpl);
//        printf("\nTotal Number of Parameters: %d\tSample size: %d\tLarge BIC\n", ncov, n);
//        printf("Parameters:\n");
//        for ( i = 0; i < numpl; i++) {
//          printf("\tparam %d: %lf\n", i, anspl[i]);
//        }
//
//      free(delp);
        free(anspl);

        printf("\nNumber of groups chosen (large BIC): %d\n", numgl);
        printf("Groups:\n");
        for ( i = 0; i < numgl; i++) {
          printf("\tgroup %d: %d\n", i, ansgl[i]);
        }

//      free(delg);
        free(ansgl);

//        printf("\nBeta Hat estimates:\n");
//        print_matrix(beta);

//        printf("\nResults for Small BIC\n");
//        printf("Number of Parameters chosen: %d\n", numps);
//        printf("\nTotal Number of Parameters: %d\tSample size: %d\tSmall BIC\n", ncov, n);
//        printf("Parameters:\n");
//        for ( i = 0; i < numps; i++) {
//          printf("\tparam %d: %lf\n", i, ansps[i]);
//        }
//
//      free(delp);
        free(ansps);

        printf("\nNumber of groups chosen (small BIC): %d\n", numgs);
        printf("Groups:\n");
        for ( i = 0; i < numgs; i++) {
          printf("\tgroup %d: %d\n", i, ansgs[i]);
        }

//      free(delg);
        free(ansgs);

//        printf("\nBeta Hat estimates:\n");
//        print_matrix(beta);

        cout << "\nTime: " << t << " seconds" << endl;
      } else {

        int vartype = atoi(argv[6]);
        int sigma = atoi(argv[7]);
        int n_grp=ncov/tgrp;
        // Large BIC results
        printf("%d,1,%d,%d,%d,%d,%d,%d,", n, ncov, tgrp, sigma, vartype, n_grp, numgl);

        for ( i = 0; i < numgl; i++) {
          if (i != (numgl -1)) {
            printf("%d,", (int) ansgl[i]);
          } else {
            printf("%d\n", (int) ansgl[i]);
          }
        }

//      free(delg);
        free(ansgl);
//      free(delp);
        free(anspl);

        for ( i = 0; i < ncov; i++) {
          if (i != (ncov - 1)) {
            printf("%lf,", betal(i));
          } else {
            printf("%lf\n", betal(i));
          }
        }
        printf("%lf\n", t);
        // Small BIC results
        printf("%d,0,%d,%d,%d,%d,%d,%d,", n, ncov, tgrp, sigma, vartype, n_grp, numgs);

        for ( i = 0; i < numgs; i++) {
          if (i != (numgs -1)) {
            printf("%d,", (int) ansgs[i]);
          } else {
            printf("%d\n", (int) ansgs[i]);
          }
        }

//      free(delg);
        free(ansgs);
//      free(delp);
        free(ansps);
        free(all);
        for ( i = 0; i < ncov; i++) {
          if (i != (ncov - 1)) {
            printf("%lf,", betas(i));
          } else {
            printf("%lf\n", betas(i));
          }
        }
        printf("%lf\n", t);
      }
    }

    free(BICl);
    free(BICs);
    free(r);
    free(subs);

  } else {  
    int id, job, size, grp;
    int flag = 1;

    // Every node does the samee thing.  Every node has a design matrix (X),
    //   and for Foward Regression to work correctly X must be identical.
    //   After finding the signficant factor (each step) the columns of those
    //   significant factors are to be removed and the X is recalculated.
    while (flag == 1) { // The node is waiting for a message
      MPI_Recv(recvb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      grp = (int) recvb[3]; // Unique group identifier
      id = (int) recvb[2]; // Column number where group starts
      job = (int) recvb[1]; // What is being calculated?
      size = (int) recvb[0]; // Group size
//printf("n = %d\tsize = %d\tid = %d\tjob = %d\tgrp = %d\tnode = %d\n", n, size, id, job, grp, me);
      fmat tmpX(n, size);

      // If the code received is 0 then the RSS needs to be calculated
      if ( job == 0) {

        // Using the id that was received get the correct matrix
        for (j = 0; j < size; j++) {
          tmpX.col(j) = X.col(id + j);
        }

        // Create send buffer
        sendb[0] = RSS(y, tmpX, n); // Calculate RSS based on that matrix
        sendb[1] = (double) grp; // Send back the group whose RSS was calc.
        sendb[2] = (double) size; // Send back the size of the group
        sendb[3] = (double) me; // And the node number
//printf("n = %d\tsize = %d\tid = %d\tjob = %d\tgrp = %d\tnode = %d\tRSS = %lf\n", n, size, id, job, grp, me, sendb[0]);

        // Send buffer to Node 0
        MPI_Send(sendb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

      // If code received is 1 then X needs to be edited
      } else if (job == 1) {

        // Create matrix with id that was received
        for ( i = 0; i < size; i++) {
          tmpX.col(i) = X.col(id + i) / norm(X.col(id + i), 2);
        }

        // Remove those columns from X
        for ( i = 0; i < size; i++){
          X.shed_col(id);
        }

        // Recalculate X and y
        y = y - tmpX * (tmpX.t() * y);
        X = X - tmpX * (tmpX.t() * X);

        // Send back confirmation
        MPI_Send(sendb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

      // If code is something else, kill node
      } else {
        flag = -1;
      }

      tmpX.reset();
    }
  }

  free(sendb);
  free(recvb);
  MPI_Finalize();
  return 0;

}
