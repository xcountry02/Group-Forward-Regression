/*=============================================================================
 |      Program:  gFRmpiINT.cpp
 |
 |       Author:  K. Michels
 |     Language:  C++ (mpic++ compiler)
 |   To Compile:  mpic++ -g -Wall -larmadillo gFRmpiINT.cpp -o gFRmpiINT (it is 
 |                assumed that the stack header file, header.h, is in the same 
 |                directory.)      
 |
 +-----------------------------------------------------------------------------
 |
 |  Description:  Given the correct input files, this program does my algorithm
 |                of group forward regression to find the significant groups
 |                in the input file and significant interactions.
 |
 |        Input:  There are two input options.  The first is PLINK format,
 |                which has two data files as input.  First the *.ped, 
 |                explained at this website:
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
 |                design matrix where the continous covariates are in the first 
 |                'numc' columns.  Then each column after that is non-continuou,
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
 |       Output:  The output produced includes a list of significant groups for
 |                all combinations of main BIC and interaction BIC (4).  The 
 |                beta.hat estimates from those parameters are also calculated
 |                and reported.
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

  int i, j, k;
  wall_clock timer; // Keeps track of computation time
  int tgrp; // Total number of main groups
  int ig; // Total number of interaction groups
  int prnt=2; // Means print nothing

  double t;
  // The send and receive buffers
  double *sendb = (double *)malloc(4*sizeof(double));
  double *recvb = (double *)malloc(4*sizeof(double));
  // The integer send receive buffer for interaction
  int *sendib = (int *)malloc(2500*sizeof(int));
  int *recvib = (int *)malloc(2500*sizeof(int));
  int np; // Total number of nodes
  int me; // Current node

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  // Every node reads in the data in parallel //
  string type = argv[1]; // Either PM or STD
  int n = 0; // Sample size
  int p0 = 0; // Number of main parameters
  int numc = 0; // Number of continuous covariates
  int csv = 0; // Two csv choices for output

  vector<string> FiD; // Family ID
  int cntFiD = 0;

  vector<string> PiD; // Paternal ID
  int cntPiD = 0;

  vector<string> MiD; // Maternal ID
  int cntMiD = 0;

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
//  printf("\nX:\tnumc=%d\n", numc);
//  print_matrix(X);
//}
// I need to calculate those every time.

  // Number of dummy variables
//  int size = n_grp - 1;
//  int isize = pow(size, 2);
//
//  fmat tmpX(n, size);
//  fmat itmpX(n, isize);


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
  int icov=0;
  fmat X0(n, ncov);
  X0 = X;
  fvec y0(n);
  y0 = y;

  // This is for Node 0 only.  Node 0 is the main Node that controls everything
  //   it controls all sending and receiving of messages.
  if (me == 0){

    // Start the group FR //
    timer.tic();

    // 'nStep' is the number of groups that will be added to the model.
//    int nStep = min(ncov/size, (n - 1)/size);
    int nStep = numc; // Covariates added to model

//printf("nStep: %d\n", nStep);
    vec L2(n + 1); // Likelihood vector
    fvec MmFR(n + 1); // SNP vector (keeps parameters)
    fvec MgFR(n + 1); // Significant groups
    fvec ID(ncov); // Stores all SNP ids

    for (i=0; i<ncov; i++)
        ID(i) = i;

//    int len = ncov; // Number of parameters
    int len2 = tgrp; // Number of groups

    vector< vector<string> > name(pow(n,2), vector<string>(2));
    int **ansp, **ansg, ***ansip, ***ansig;
    double ***beta;
    ansp = (int **)malloc(2*sizeof(int*));
    ansp[0] = (int *)malloc(n*sizeof(int));
    ansp[1] = (int *)malloc(n*sizeof(int));
    ansg = (int **)malloc(2*sizeof(int*));
    ansg[0] = (int *)malloc(n*sizeof(int));
    ansg[1] = (int *)malloc(n*sizeof(int));
    ansip = (int ***)malloc(2*sizeof(int**));
    ansip[0] = (int **)malloc(sizeof(int*));
    ansip[1] = (int **)malloc(sizeof(int*));
    ansip[0][0] = (int *)malloc(n*sizeof(int));
    ansip[0][1] = (int *)malloc(n*sizeof(int));
    ansip[1][0] = (int *)malloc(n*sizeof(int));
    ansip[1][1] = (int *)malloc(n*sizeof(int));
    ansig = (int ***)malloc(2*sizeof(int**));
    ansig[0] = (int **)malloc(sizeof(int*));
    ansig[1] = (int **)malloc(sizeof(int*));
    ansig[0][0] = (int *)malloc(n*sizeof(int));
    ansig[0][1] = (int *)malloc(n*sizeof(int));
    ansig[1][0] = (int *)malloc(n*sizeof(int));
    ansig[1][1] = (int *)malloc(n*sizeof(int));
    beta = (double ***)malloc(2*sizeof(double**));
    beta[0] = (double **)malloc(sizeof(double*));
    beta[1] = (double **)malloc(sizeof(double*));
    beta[0][0] = (double *)malloc((ncov+n)*sizeof(double));
    beta[0][1] = (double *)malloc((ncov+n)*sizeof(double));
    beta[1][0] = (double *)malloc((ncov+n)*sizeof(double));
    beta[1][1] = (double *)malloc((ncov+n)*sizeof(double));
    
    int nump[2], numg[2], numip[2][2], numig[2][2];
    int *size = (int *)malloc(n*sizeof(int));
    int *newCnum = (int *)malloc(n*sizeof(int));
    // Collect the RSS values
    double *r = (double *)malloc(len2*sizeof(double));


    double *subs = (double *)malloc(len2*sizeof(double));

    int cnt=0;

    double minRSS;
    int minID;

    // The following loop is the Forward Regression loop.  Each step of the 
    //   loop creates the RSS vector (r) using MPI.  Each step finds the 
    //   most significant group
    int k=0;

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
//printf("len2 = %d\tnp = %d\n", len2, np);
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
//cout << "here1" << endl;

        // Wait to receive the RSS from the 'len2' processes
        for ( i = 0; i < len2; i++ ){
//printf("i=%d\n",i);
          MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          r[i] = recvb[0]; // Put the RSS in the RSS vector
          subs[i] = recvb[1]; // The corresponding group (to RSS) is recorded
        }
//cout << "here2" << endl;

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
//printf("RSS = %lf\tgroup num = %d\n", recvb[0], (int) recvb[1]);
        }

        // Wait to receive RSS from the last if the MPI processes
        for ( i = 1; i < np; i++ ){
          MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          r[cnt-np] = recvb[0]; // Put RSS in RSS vector
          subs[cnt++-np] = recvb[1]; // Column number
//printf("RSS = %lf\tgroup num = %d\n", recvb[0], (int) recvb[1]);
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
      size[k] = nINgrp(minID);
      nStep = nStep + size[k];
//printf("nStep = %d\n", nStep);
      // If this is the last step of the Forward Regression loop, then the
      //   MPI processes are done so send the kill signal to all nodes.
      if ((nStep >= n) || (len2 == 1)) {
        for ( i = 1; i < np; i++ ){
          sendb[1] = 2; // Kill code
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
          sendb[0] = (double) size[k]; // Group size
          sendb[1] = 1; // Code to delete group from design matrix X
          sendb[2] = (double) cNum(minID); // Column Number
          sendb[3] = (double) minID; // Group identifier
          MPI_Send(sendb, 4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
      }


//printf("k = %d\tnStep = %d\tn = %d\n", k, nStep, n);

      if (nStep < n) {

        int ak = ID(minID); // Significant group number
        MgFR(k) = ak; // Record significant group number
        fmat tmpX(n, size[k]);

        for (i = 0; i < size[k]; i++) {
          MmFR(nStep - size[k] + i) = cNum_org(ak) + i; // Add the group of subscripts to mFR

          // Create matrix with columns of minimum RSS
          tmpX.col(i) = X.col(cNum(minID) + i) / norm(X.col(cNum(minID) + i), 2);
        }

        // Must remove those columns from the design matrix
        for (i = 0; i < size[k]; i++) {
          X.shed_col(cNum(minID));
        }

        ID.shed_row(minID); // Remove group of subscripts from subscript vector
        nINgrp.shed_row(minID);
        cNum.shed_row(minID);

        for (int h=minID; h<len2 - 1; h++) {
          cNum(h) = cNum(h) - size[k];
        }    
  
        // Recalculate both y and X
        y = y - tmpX * (tmpX.t() * y);
        X = X - tmpX * (tmpX.t() * X);
        tmpX.reset();
  
        len2 = len2-1; // the number of parameters taken out of X is size
//if (k==-1) {
//printf("ID = %d\tsize = %d\t", minID, (int) size(k));
//for (int h=0; h<len2; h++){
//  printf("cNum[%d] = %d\n", h, (int) cNum(h));
//}
//}
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

    // Calculate model size at each step
    fvec newvec(len + 1);
    newvec(0)=numc; // Covariates added to model
    for (i=1; i<=len; i++){
      newvec(i) = size[i-1] + newvec(i-1);
    }

    // Calculate both small BIC and large BIC.  All possible values are returned
    // Calculate BIC (using small)
    double *mgBICs = (double *)malloc((len + 1)*sizeof(double));

    // Calculate BIC (using small)
    double *mgBICl = (double *)malloc((len + 1)*sizeof(double));

    // Small BIC calculation
    for ( i = 0; i <= len; i++) {
      mgBICs[i] = 2 * log(L2(i)) + ((newvec(i) + 1) * log(n)) / n;
    }

    // Large BIC calculation
    for ( i = 0; i <= len; i++) {
      mgBICl[i] = 2 * log(L2(i)) + ((newvec(i) + 1) * (log(n) + 2 * log(ncov))) / n;
    }
  
    L2.reset();
    ID.reset();
    newvec.reset();

    // Find the minimum BIC value and the associated subscript.  This subscript
    //   is the `best' model according to BIC
    double minbs = min(mgBICs, len + 1);
    double minbl = min(mgBICl, len + 1);


    // The values recored pertain to both small BIC and large BIC.  They must
    //   be entered into the beginning of the next step.
    numg[0] = which(mgBICs, len + 1, minbs); // Small BIC
    numg[1] = which(mgBICl, len + 1, minbl); // Small BIC

    free(mgBICs);
    free(mgBICl);

    nump[0]=0;
    nump[1]=0;

    // Calculate cNum for both.
    int cNumB[2][n];
    for (int a=0; a<2; a++){
      for (i=0; i<numg[a]; i++){
        nump[a] = size[i] + nump[a];
        cNumB[a][i] = nump[a];
      }
    }

    // This section is adding the interaction terms to the model.  I need to
    //   have two iterators because the main terms were selected by both
    //   small BIC (a=0) and large bic (a=1).  Everything created in the loop
    //   is locally only.
    for (int a=0; a<2; a++){
      // Parameter vector
      fvec ImFR(n+1);
      // Output the main results into ansp[a] and ansg[a] where a=0,1
      if (numg[a] > 0) {
        // Takes the first `Tlen' from mFR to ansp[a]
        int tmp=0;
        for (i = 0; i < nump[a]; i++) {
          ansp[a][i] = MmFR(i);

          if ((i+1) == cNumB[a][tmp]){
            ansg[a][tmp] = MgFR(tmp);
            tmp++;
          }
        }
      } else {
        ansp[a][0] = 0;
        ansg[a][0] = 0;
      }
      
      // Quicksort before gathering interacts  
      quicksort( ansp[a], 0, nump[a] - 1);

      quicksort( ansg[a], 0, numg[a] - 1);
 
      // Finding the starting position of each group is crucial 
      newCnum[0]=0;
      size[0]=nINgrp_org(ansg[a][0]);

      int tmp=0;
      for ( i = 1; i < numg[a]; i++) {
        size[i] = nINgrp_org(ansg[a][i]);
        tmp=size[i-1];
        newCnum[i]=newCnum[i-1]+tmp;
      }
  
      fmat mMAT(n, nump[a]); // Create main matrix

      for ( i=0; i < nump[a]; i++) {  
        mMAT.col(i) = X0.col(ansp[a][i]);
      }

      // Check the interactions and get rid of any groups that produce
      //   a column of zeroes  
      if (numg[a] > 1) { // More than one group (can be interactions)
        for ( k = 0; k < numg[a]; k++) { // Iterate over main groups
          // Interaction only with different column
          for ( int l = k+1; l < numg[a]; l++) { 
            int tmp=size[k]*size[l]; // size of single interaction group
            fmat test(n, tmp); // Create temporary test interaction matrix
    
            int intcnt=0; // Iterator
            int iflag=0; // Flag, check if any vector of zeroes
            // Iterate over two main groups
            for ( i = (int) newCnum[k]; i < (int) newCnum[k] + (int) size[k]; i++) {
              for ( j = (int) newCnum[l]; j < (int) newCnum[l] + (int) size[l]; j++) {
                test.col(intcnt) = mMAT.col(i) % mMAT.col(j); // Multiply vectors
                if (sum(test.col(intcnt))==0) { // Check if zero
                  iflag=1; // If so do not include group
                  break;
                }
                intcnt++;
              }
            }
            if (iflag == 0) { // If all nonzero, count interaction term
              ig++;
              icov=icov+tmp; // Find total number of parmeters
            }
            test.reset();
          }
        }
      }
//  printf("Number of interaction groups: %d\tNumber of interaction parameters: %d\n", ig, icov);
  
  
  //    X0.reset();
  //    X.reset();
      //  If the total number of groups is 1 or less no interactions can be
      //    found.  Kill all processes
      if ( numg[a] <= 1) {
        if (a==0) {
          break;
        }
        t = timer.toc(); // End time
        for ( i = 1; i < np; i++ ){ // Send kill signal
          sendib[0] = -1;
          sendib[1] = nump[a]; // Number of main parameters
          sendib[2] = numg[a]; // Number of main groups
          sendib[3] = icov; // Number of POSSIBLLE interaction parameters
          MPI_Send(sendib, 4, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
        //  No interactions.
        for (int b=0; b<2; b++) {
          numig[a][b]=0;
          numip[a][b]=0;
          ansip[a][b][0] = 0;
          ansig[a][b][0] = 0;
coefi(y0, mMAT, (int) n, (int) ncov, (int) nump[a], (int) numg[a], (int) numip[a][b], (int*) size, (int*) newCnum, (int *) ansp[a], (int *) ansip[a][b], beta[a][b]);
        }
      } else {
  
        //  If the number of interactions (when counting) is zero, then all
        //    processes must be killed.
        if (ig == 0) {
          if (a==0) { // Still need to do other
            break;
          }
          t = timer.toc(); // End time
          for ( i = 1; i < np; i++ ){ // Send kill signal
            sendib[0] = -1;
            sendib[1] = nump[a]; // Number of main parameters
            sendib[2] = numg[a]; // Number of main groups
            sendib[3] = icov; // Number of POSSIBLLE interaction parameters
            MPI_Send(sendib, 4, MPI_INT, i, 1, MPI_COMM_WORLD);
          }
          //  No interactions.
          for (int b=0; b<2; b++) {
            numig[a][b]=0;
            numip[a][b]=0;
            ansip[a][b][0] = 0;
            ansig[a][b][0] = 0;
  // I really need to fix the coefficient vectors  
coefi(y0, mMAT, (int) n, (int) ncov, (int) nump[a], (int) numg[a], (int) numip[a][b], (int*) size, (int*) newCnum, (int *) ansp[a], (int *) ansip[a][b], beta[a][b]);
          }
      
        } else {
          // Create interaction matrix  
          fmat iMAT(n, nump[a] + icov);
    
          // Enter parameters for matrix
          for ( i = 0; i < nump[a]; i++) {
            iMAT.col(i) = mMAT.col(i);
          }
  
          int cntp=0; // Count interaction parameters
          int cntg=0; // Count interaction groups
    
          int inum=0;
          fvec inINgrp(ig); // Interaction group size
          fvec icNum(ig); // Subscript interaction group
          fvec icNum_org(ig);
          for ( k = 0; k < numg[a]; k++) {
            for ( int l = k+1; l < numg[a]; l++) {
              int tmp=size[k]*size[l];
              fmat test(n, tmp);
    
              int intcnt=0;
              int iflag=0;
              for ( i = (int) newCnum[k]; i < (int) newCnum[k] + (int) size[k]; i++) {
                for ( j = (int) newCnum[l]; j < (int) newCnum[l] + (int) size[l]; j++) {
                  test.col(intcnt) = mMAT.col(i) % mMAT.col(j);
                  if (sum(test.col(intcnt))==0) {
                    iflag=1;
                    break;
                  }
                  intcnt++;
                }
              }
              if (iflag == 0) {
                for (int b = 0; b < tmp; b++) {
                  iMAT.col(cntp) = test.col(b);
                  cntp++;
                }
                icNum(cntg)=inum;
                icNum_org(cntg)=inum;
                inINgrp(cntg) = tmp;
                inum = inINgrp(cntg) + inum;
                stringstream s;
                s << k;
                s << ",";
                s << l;
                name[cntg][a]=s.str().c_str(); // Record interaction name
                cntg++;
              }
              test.reset();
            }
          }

          y = y0;

if ((me == 0) && (prnt==0)) {
  for (i = 0; i < numg[a]; i++) {
    printf("Node: %d\tGroup: %d\tsize[%d]=%d\tcNum[%d]=%d\n", me, (int) ansg[a][i], i, (int) size[i], i, (int) newCnum[i]);
  }
  printf("\n(1) iMAT:\tNumber of INT groups: %d\tNumber of INT parameters: %d\n", ig, icov);
  print_matrix(iMAT);
  printf("\n(1) y\n");
  print_matrix(y);
}

          // Add main parameters to model  
          for (j = 0; j < numg[a]; j++) {
            fmat tmpX(n, size[j]);
            for (i = 0; i < size[j]; i++) {
              tmpX.col(i) = iMAT.col(i) / norm(iMAT.col(i), 2);
            }
    
//          for (i = 0; i < size(j); i++) {
//            iMAT.shed_col(0);
//  
//            y = y - tmpX.col(i) * (tmpX.col(i).t() * y);
//            iMAT = iMAT - tmpX.col(i) * (tmpX.col(i).t() * iMAT);
//          }
  
            // Must remove those columns from the design matrix
            for (i = 0; i < size[j]; i++) {
              iMAT.shed_col(0);
            }
    
            // Recalculate both y and X
            y = y - tmpX * (tmpX.t() * y);
            iMAT = iMAT - tmpX * (tmpX.t() * iMAT);
            tmpX.reset();
          }
// cout << "here4" << endl; 
  
if ((me == 0) && (prnt ==0)) {
  printf("\n(2) iMAT:\n");
  print_matrix(iMAT);
  printf("\n(2) y:\n");
  print_matrix(y);
}
  
          for ( i = 1; i < np; i++ ){
            sendib[0] = 1;
            sendib[1] = nump[a]; // Number of main parameters
            sendib[2] = numg[a]; // Number of group parameters
            sendib[3] = icov; // Number of POSSIBLE interaction parameters
            for (j = 0; j < numg[a]; j++) {
              sendib[j+4] = ansg[a][j];
            }
            for (j = 0; j < nump[a]; j++) {
              sendib[j+numg[a]+4] = ansp[a][j];
            }
            MPI_Send(sendib, numg[a] + nump[a] + 4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
          }
    
 




//newCnum.reset();

 
  
  
  
  
// cout << "here5" << endl; 
 
          // 'iStep' is the number of groups that will be added to the model.
      //    int iStep = min(ncov/size, (n - 1)/size);
      //printf("nump[a] = %d\n", nump[a]);
          int iStep = nump[a] + numc; // Main groups in model and covariates
      
      //printf("iStep: %d\n", iStep);
          vec L2(n + 1); // Likelihood vector
          vec isize(n + 1); // Group sizes for each group
          fvec IgFR(n + 1); // Significant groups
          fvec ID(icov); // Stores all SNP ids
      
          for (i=0; i<icov; i++)
              ID(i) = i;
      
      //    int len = ncov; // Number of parameters
          int len2 = ig; // Number of groups
      
      
          // Collect the RSS values
          double *ri = (double *)malloc(len2*sizeof(double));
      
          // Collect the SNP associated with the RSS
          double *subsi = (double *)malloc(len2*sizeof(double));
      
          int cnt=0;
      
          double minRSS;
          int minID;
      
          // The following loop is the Forward Regression loop.  Each step of the 
          //   loop creates the RSS vector (r) using MPI.  Each step finds the 
          //   most significant group
          int k=0;
      
          while ((iStep < n) && ( len2 > 0 )) {
      
            // After the first step, because on the first step we know the Nodes's
            //   are ready, Node 0 must wait for the other Nodes to confirm
            //   they are ready before proceeding any further.
// cout << "here6" << endl; 
      
            for ( i = 1; i < np; i++ ){
              MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
// cout << "here7" << endl; 
      
            // Calculates the last Likelihood.  Proven in the paper.
            L2(k) = norm(y, 2);
      
//double tmpL2 = L2(k);
//printf("L[%d] = %.10lf\n", k, tmpL2);
      
            cnt = 0;
//      printf("len2 = %d\tnp = %d\n", len2, np);
            // If the number of groups left in X is less the the number of MPI
            //   processors
            if ( len2 < np) {
              cnt++;
      
              for ( i = 0; i < len2; i++){ // Only allocate 'len2' MPI processors
                // Prepare send buffer
                sendb[0] = (double) inINgrp(cnt - 1); // Group size
                sendb[1] = (double) 0; // Code that it's for computing RSS
                sendb[2] = (double) icNum(cnt - 1); // Column number 
                sendb[3] = (double) cnt - 1; // Group number
                MPI_Send(sendb, 4, MPI_DOUBLE, cnt++, 1, MPI_COMM_WORLD);
      
//      if (k == 0) {
//        printf("Group Size: %d\tColumn Number: %d\tGroup Number: %d\n", (int) sendb[0], (int) sendb[2], (int) sendb[3]);
//      }
            }
      //cout << "here1" << endl;
              // Wait to receive the RSS from the 'len2' processes
      //printf("RSS = %lf\tgroup num = %d\n", recvb[0], (int) recvb[1]);
              for ( i = 0; i < len2; i++ ){
      //printf("i=%d\n",i);
                MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                ri[i] = recvb[0]; // Put the RSS in the RSS vector
                subsi[i] = recvb[1]; // The corresponding group (to RSS) is recorded
      //printf("RSS = %lf\tgroup num = %d\n", recvb[0], (int) recvb[1]);
              }
      //cout << "here2" << endl;
      
      
            // If the number of groups left is equal to or more than the number
            //   of MPI processes
            } else {
              cnt++;
              // allocate to all MPI processes
              for ( i = 1; i < np; i++){
                // Prepare send buffer
                sendb[0] = (double) inINgrp(cnt - 1);
                sendb[1] = (double) 0; // Code that it's for computing RSS
                sendb[2] = (double) icNum(cnt - 1); // column number
                sendb[3] = (double) cnt - 1; // Group number
                MPI_Send(sendb, 4, MPI_DOUBLE, cnt++, 1, MPI_COMM_WORLD);
//      if (k == 0) {
//        printf("Group Size: %d\tColumn Number: %d\tGroup Number: %d\n", (int) sendb[0], (int) sendb[2], (int) sendb[3]);
//      }
              }
      
              // Wait until one of the processes is done.  When it is done save the 
              //   receiver buffer.  To save time send information from Node 0 to 
              //   calculate the next RSS.  Keep going in this manner until
              //   all RSS values are calculated.
              while (cnt <= len2) {
                // Receive message, put in buffer
                MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     cout << "hello" << endl;
                // Prepare send buffer
                sendb[0] = inINgrp(cnt - 1);
                sendb[1] = (double) 0; // Code for computing RSS
                sendb[2] = (double) icNum(cnt - 1); // column number
                sendb[3] = (double) cnt - 1; // Group number
      
                // Send to the same MPI process that was received from (recvb[2])
                MPI_Send(sendb, 4, MPI_DOUBLE, (int) recvb[2], 1, MPI_COMM_WORLD);
      
                ri[cnt - np] = recvb[0]; // put RSS in the RSS vector
                subsi[cnt++ - np] = recvb[1]; // The corresponding group (to RSS)
//      printf("RSS = %lf\tgroup num = %d\n", recvb[0], (int) recvb[1]);
              }
      
              // Wait to receive RSS from the last if the MPI processes
              for ( i = 1; i < np; i++ ){
                MPI_Recv(recvb, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                ri[cnt-np] = recvb[0]; // Put RSS in RSS vector
                subsi[cnt++-np] = recvb[1]; // Column number
      //printf("RSS = %lf\tgroup num = %d\n", recvb[0], (int) recvb[1]);
              }
            }
      
//      if (k == 0) {
//        for (i=0; i<len2; i++) {
//          printf("ri[%d] = %lf\n", i, ri[i]);
//        }
//      }
      
      
            // Find the minimum RSS value
            minRSS = min(ri, len2);
      
            // And the associated subscript (corresponding group to the minimum
            //   RSS value.
            int wID = which(ri, len2, minRSS);
            minID = subsi[wID];
            isize(k) = inINgrp(minID);
            iStep = iStep + isize(k);
      //printf("iStep = %d\n", iStep);
            // If this is the last step of the Forward Regression loop, then the
            //   MPI processes are done so send the kill signal to all nodes.

//printf("iStep=%d\tlen2=%d\ta=%d\n", iStep, len2, a);
//if (a==1) {
//  printf("a=%d\tn=%d\tlen2=%d\tiStep=%d\n", a, n, len2, iStep);
//}
            if ((iStep >= n) || (len2 == 1)) {
              for ( i = 1; i < np; i++ ){
                if (a==0) {
                  sendb[1] = 2.0; // Kill code
                } else {
//printf("Kill Code Sent\ta=%d\n", a);
                  sendb[1] = -1.0; // Kill code
                }
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
                sendb[0] = (double) isize(k); // Group size
                sendb[1] = 1; // Code to delete group from design matrix X
                sendb[2] = (double) icNum(minID); // Column Number
                sendb[3] = (double) minID; // Group identifier
                MPI_Send(sendb, 4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
              }
            }
      
      
      
      //printf("k = %d\tiStep = %d\tn = %d\n", k, iStep, n);
      
            if (iStep < n) {
      
              int ak = ID(minID); // Significant group number
              IgFR(k) = ak; // Record significant group number
              fmat tmpX(n, isize(k));
      
              for (i = 0; i < isize(k); i++) {
                ImFR(iStep - isize(k) + i) = icNum_org(ak) + i; // Add the group of subscripts to mFR
                // Create matrix with columns of minimum RSS
                tmpX.col(i) = iMAT.col(icNum(minID) + i) / norm(iMAT.col(icNum(minID) + i), 2);
              }
      
              // Must remove those columns from the design matrix
              for (i = 0; i < isize(k); i++) {
                iMAT.shed_col(icNum(minID));
              }
      
              ID.shed_row(minID); // Remove group of subscripts from subscript vector
              inINgrp.shed_row(minID);
              icNum.shed_row(minID);
      
              for (int h=minID; h<len2 - 1; h++) {
                icNum(h) = icNum(h) - isize(k);
              }
      
              // Recalculate both y and X
              y = y - tmpX * (tmpX.t() * y);
              iMAT = iMAT - tmpX * (tmpX.t() * iMAT);
              tmpX.reset();
      //cout << " Did I make it here?" << endl;
      
              len2 = len2-1; // the number of parameters taken out of X is size
      //if (k==-1) {
      //printf("ID = %d\tsize = %d\t", minID, (int) isize(k));
      //for (int h=0; h<len2; h++){
      //  printf("cNum[%d] = %d\n", h, (int) icNum(h));
      //}
      //}
            }
            k++;
      //cout << " Did I make it here?" << endl;
          }
      //cout << " Did I make it here?" << endl;
          int len=k-1;
          if (iStep < n) {
            L2(k) = norm(y, 2);
            len=k;
//double tmpL2 = L2(k);
//printf("L[%d] = %.10lf\n", k, tmpL2);
          }

//printf("L2:\n");
//print_matrix(L2);
          //free(ri);
          //free(subsi);
//cout << "here1" << endl;
          fvec newvec(len + 1);
          newvec(0)=nump[a]+numc;
          for (i=1; i<=len; i++){
            newvec(i) = isize(i-1) + newvec(i-1);
          }
//cout << "here2" << endl;
      
          // Calculate BIC (using small)
          double *igBICs = (double *)malloc((len + 1)*sizeof(double));
          // Calculate BIC (using large)
          double *igBICl = (double *)malloc((len + 1)*sizeof(double));
      
//cout << "here3" << endl;
          // Small BIC calculation
          for ( i = 0; i <= len; i++) {
            igBICs[i] = 2 * log(L2(i)) + ((newvec(i) + 1) * log(n)) / n;
          }
          // Large BIC calculation
          for ( i = 0; i <= len; i++) {
            igBICl[i] = 2 * log(L2(i)) + ((newvec(i) + 1) * (log(n) + 2*log(ncov+icov) / n;
          }

          // Find the minimum BIC value and the associated subscript.  This subscript
          //   is the `best' model according to BIC
          double minbs = min(igBICs, len + 1);
          double minbl = min(igBICl, len + 1);
          numig[a][0] = which(igBICs, len + 1, minbs);
          numig[a][1] = which(igBICl, len + 1, minbl);
          free(igBICs);
          free(igBICl);
      
          numip[a][0]=0;
          numip[a][0]=0;
          int icNumB[2][n]; // cNumB instead of cNumL
          for (int b=0; b<2; b++) {
            for (i=0; i<numig[a][b]; i++){
              numip[a][b] = isize(i) + numip[a][b];
              icNumB[b][i] = numip[a][b];
            }
          }
        
          for (int b=0; b<2; b++) {
            if (numig[a][b] > 0) {
              // Takes the first `Tlen' from mFR to ansip
              int tmp=0;
              for (i = 0; i < numip[a][b]; i++) {
                ansip[a][b][i] = ImFR(i);
                if ((i+1) == icNumB[b][tmp]){
                  ansig[a][b][tmp] = IgFR(tmp);
                  tmp++;
                }
              }
            } else {
              ansip[a][b][0] = 0;
              ansig[a][b][0] = 0;
            }
//cout << "here8" << endl;
//for (i=0; i< nump[a]; i++) {
//  printf("ansp[%d]=%d\n", i, (int) ansp[a][i]);
//}
            quicksort( ansip[a][b], 0, numip[a][b] - 1);
            quicksort( ansig[a][b], 0, numig[a][b] - 1);
coefi(y0, mMAT, (int) n, (int) ncov, (int) nump[a], (int) numg[a], (int) numip[a][b], (int*) size, (int*) newCnum, (int *) ansp[a], (int *) ansip[a][b], beta[a][b]);
//cout << "here9" << endl;
          }
free(ri);
free(subsi);
ID.reset();
          ImFR.reset();
          IgFR.reset();
isize.reset();
      L2.reset();
      ID.reset();
      newvec.reset();
      mMAT.reset();
      iMAT.reset();
      inINgrp.reset();
//      icNUM.reset();
//      icNUM_org.reset();
        }
//cout << "here10" << endl;
      }
//cout << "here11" << endl;
    }
    t = timer.toc();


//for (int a=0; a<2; a++) {
//  printf("numg[%d]=%d\tnump[%d]=%d\n", a, (int) numg[a], a, (int) nump[a]);
////  for (int b=0; b<numg[a]; b++) {
////    printf("ansg[%d][%d]=%d\n", a, b, (int) ansg[a][b]);
////  }
//  for (int c=0; c<2; c++) {
//    printf("numig[%d][%d]=%d\tnumip[%d][%d]=%d\n", a, c, (int) numig[a][c], a, c, (int) numip[a][c]);
//  }
//}

    // If PEDMAP or PLINK format.  In 'ansp' are parameter numbers.  With PLINK
    //   format the .map file contains the SNP name.  We return the SNP names
    //   with PLINK format instead of just parameter numbers.
    if (type.compare("PM") == 0){
   
      // Open .map file
      char *map = argv[3];
      ifstream fp;
      string l;

      printf("Total Number of SNPs: %d\tSample size: %d\n", tgrp, n);

      for (int a=0; a<2; a++) {
        fp.open(map, ios::in);

        vector<string> ChromoName; // Chromosome Name
        int cntChromoName = 0;
  
        int total = 0; // Iteration of while loop
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
          if (cntChromoName < numg[a]) {
            // Check parameters
            if (total == (int) ansg[a][cntChromoName]) {
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

        for (int b=0; b<2; b++) {

          if ((a==0) && (b==0)) {
            printf("\nMain Effects (small BIC):\n");
            printf("Number of main SNPs chosen: %d\n", numg[a]);
            printf("\tSNPs:\tName\t\t\tNumber\n");
            for (i = 0; i < numg[a]; i++) {
              for (j=newCnum[i]; j<newCnum[i]+size[i]; j++) {
                if ( j == newCnum[i]) {
                  printf("\tSNP %d:\t", i);
                  cout << ChromoName.at(i) << "\t\t\t" << ansg[a][i];
                  printf("\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansg[a][j]]);
                } else {
                  printf("\t\t\t\t\t\t\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansp[a][j]]); 
                }
              }
            }
            printf("\nInteraction Effects (small BIC):\n");
            if (numig[a][b] == 0) {
              printf("None\n");
            } else {
              printf("Number of interation SNPs chosen: %d\n", numig[a][b]);
              printf("\tINTs:\tName\t\t\tNumber\n");
              for (i = 0; i < numig[a][b]; i++) {
                std::string input = name[i][a];
                std::istringstream ss(input);
                std::string token;
                  
                int flag=0;
                while(std::getline(ss, token, ',')) {
                  int tmp=atoi(token.c_str());
                  if (flag == 0) {
                    printf("\tINT %d:\t%s--", i, ChromoName.at(tmp).c_str());
                    flag++;
                  } else {
                    printf("%s\t%d\n", ChromoName.at(tmp).c_str(), (int) ansig[a][b][i]);
                    flag=0;
                  }
                }
              }
            }
          } else if ((a==0) && (b==1)){
            printf("\nMain Effects (small BIC):\n");
            printf("Number of main SNPs chosen: %d\n", numg[a]);
            printf("\tSNPs:\tName\t\t\tNumber\n");
            for (i = 0; i < numg[a]; i++) {
              for (j=newCnum[i]; j<newCnum[i]+size[i]; j++) {
                if ( j == newCnum[i]) {
                  printf("\tSNP %d:\t", i);
                  cout << ChromoName.at(i) << "\t\t\t" << ansg[a][i];
                  printf("\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansg[a][j]]);
                } else {
                  printf("\t\t\t\t\t\t\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansp[a][j]]); 
                }
              }
            }
            printf("\nInteraction Effects (large BIC):\n");
            if (numig[a][b] == 0) {
              printf("None\n");
            } else {
              printf("Number of interation SNPs chosen: %d\n", numig[a][b]);
              printf("\tINTs:\tName\t\t\tNumber\n");
              for (i = 0; i < numig[a][b]; i++) {
                std::string input = name[i][a];
                std::istringstream ss(input);
                std::string token;
                  
                int flag=0;
                while(std::getline(ss, token, ',')) {
                  int tmp=atoi(token.c_str());
                  if (flag == 0) {
                    printf("\tINT %d:\t%s--", i, ChromoName.at(tmp).c_str());
                    flag++;
                  } else {
                    printf("%s\t%d\n", ChromoName.at(tmp).c_str(), (int) ansig[a][b][i]);
                    flag=0;
                  }
                }
              }
            }
          } else if ((a==1) && (b==0)){
            printf("\nMain Effects (large BIC):\n");
            printf("Number of main SNPs chosen: %d\n", numg[a]);
            printf("\tSNPs:\tName\t\t\tNumber\n");
            for (i = 0; i < numg[a]; i++) {
              for (j=newCnum[i]; j<newCnum[i]+size[i]; j++) {
                if ( j == newCnum[i]) {
                  printf("\tSNP %d:\t", i);
                  cout << ChromoName.at(i) << "\t\t\t" << ansg[a][i];
                  printf("\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansg[a][j]]);
                } else {
                  printf("\t\t\t\t\t\t\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansp[a][j]]); 
                }
              }
            }
            printf("\nInteraction Effects (small BIC):\n");
            if (numig[a][b] == 0) {
              printf("None\n");
            } else {
              printf("Number of interation SNPs chosen: %d\n", numig[a][b]);
              printf("\tINTs:\tName\t\t\tNumber\n");
              for (i = 0; i < numig[a][b]; i++) {
                std::string input = name[i][a];
                std::istringstream ss(input);
                std::string token;
                  
                int flag=0;
                while(std::getline(ss, token, ',')) {
                  int tmp=atoi(token.c_str());
                  if (flag == 0) {
                    printf("\tINT %d:\t%s--", i, ChromoName.at(tmp).c_str());
                    flag++;
                  } else {
                    printf("%s\t%d\n", ChromoName.at(tmp).c_str(), (int) ansig[a][b][i]);
                    flag=0;
                  }
                }
              }
            }
          } else if ((a==1) && (b==1)){
            printf("\nMain Effects (large BIC):\n");
            printf("Number of main SNPs chosen: %d\n", numg[a]);
            printf("\tSNPs:\tName\t\t\tNumber\tCoefficient\n");
            for (i = 0; i < numg[a]; i++) {
              for (j=newCnum[i]; j<newCnum[i]+size[i]; j++) {
                if ( j == newCnum[i]) {
                  printf("\tSNP %d:\t", i);
                  cout << ChromoName.at(i) << "\t\t\t" << ansg[a][i];
                  printf("\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansg[a][j]]);
                } else {
                  printf("\t\t\t\t\t\t\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansp[a][j]]); 
                }
              }
            }
            printf("\nInteraction Effects (large BIC):\n");
            if (numig[a][b] == 0) {
              printf("None\n");
            } else {
              printf("Number of interation SNPs chosen: %d\n", numig[a][b]);
              printf("\tINTs:\tName\t\t\tNumber\n");
              for (i = 0; i < numig[a][b]; i++) {
                std::string input = name[i][a];
                std::istringstream ss(input);
                std::string token;
                  
                int flag=0;
                while(std::getline(ss, token, ',')) {
                  int tmp=atoi(token.c_str());
                  if (flag == 0) {
                    printf("\tINT %d:\t%s--", i, ChromoName.at(tmp).c_str());
                    flag++;
                  } else {
                    printf("%s\t%d\n", ChromoName.at(tmp).c_str(), (int) ansig[a][b][i]);
                    flag=0;
                  }
                }
              }
            }
          }
        }
      }
  cout << "\nTime: " << t << " seconds" << endl;
        //free(ansp[a]);
        //free(ansg[a]);
//      // The parameters must be sorted
//      quicksort(ansp, 0, nump - 1);
//
//      vector<string> ChromoName; // Chromosome Name
//      int cntChromoName = 0;
//
//      int total = 0; // Iteration of while loop
//   
//      // Open .map file
//      char *map = argv[3];
//      ifstream fp;
//      fp.open(map, ios::in);
//      string l;
//
//      // Go through every line of the .map file.  When the iteration of the
//      //   while loop `total' is equal to the parameter number, add that
//      //   name to ChromoName.  Since ansp is sorted, only need to go 
//      //   through .map file once. 
//      while(!fp.eof()) {
//
//        getline(fp, l);
//        stringstream line; line<<l;
//        int i;
//        string s;
//
//        line>>i; // Chromosome
//        line>>s; // String
//
//        // If ChromoName isn't full
//        if (cntChromoName < nump) {
//          // Check parameters
//          if (total == (int) ansp[cntChromoName]) {
//            ChromoName.push_back(s);
//            cntChromoName++;
//          }
//        // If ChromoName is full, break
//        } else {
//          break;
//        }
//        total++;
//      }
//      fp.close();
//
//      printf("\nNumber of Parameters chosen: %d\n", nump);
//      printf("\nTotal Number of Parameters: %d\tSample size: %d\tBIC: %d\n", ncov, n, mBIC);
//      printf("Parameters:\t Name\t\tNumber\n");
//      for ( i = 0; i < nump; i++) {
//        printf("\tparam %d: ", i);
//        cout << ChromoName.at(i) << "\t" << ansp[i] << endl;
//      }
//      //////free(ansp);
//      //////free(ansg);
//
//      cout << "\nTime: " << t << " seconds" << endl;

    } else {

//      int *delp, *delg;
//      int delpLEN, delgLEN;

//      del(X0, ncov, n, pos, size, ansp, &delp, &delpLEN, &delg, &delgLEN);    
//printf("iBIC:\n");

//      fvec beta = coefi(y0, X0, ncov, n, pos, size, ansp);
//      fvec beta = coefi(y0, mMAT, ncov, n, numg, size(0), ansp, ansip, numig);
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
//                } else {
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

//      quicksort(ansp, 0, nump - 1);


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

//      quicksort(ansg, 0, numg - 1);

      if (csv == 0 ) {
        for ( int a=0; a<2; a++) {
          for ( int b=0; b<2; b++) {
            printf("\nMain Effects");
            if (a == 0) {
              printf(" (small BIC):\n");
            } else {
              printf(" (large BIC):\n");
            }
            for (i = 0; i < numg[a]; i++) {
              for (j=newCnum[i]; j<newCnum[i]+size[i]; j++) {
                if ( j == newCnum[i]) {
                  printf("Grp: %d\tBeta[%d]=%lf\n", ansg[a][i], ansp[a][j], beta[a][b][ansp[a][j]]); 
                } else {
                  printf("\t\tBeta[%d]=%lf\n", ansp[a][j], beta[a][b][ansp[a][j]]); 
                }
              }
            }
            printf("Interactions");
            if (b == 0) {
              printf(" (small BIC):\n");
            } else {
              printf(" (large BIC):\n");
            }
            if (numig[a][b] == 0) {
              printf("None\n");
            } else {
              for (i = 0; i < numig[a][b]; i++) {
                std::string input = name[i][a];
                std::istringstream ss(input);
                std::string token;

                int flag=0;
                while(std::getline(ss, token, ',')) {
                  int tmp=atoi(token.c_str());
                  if (flag == 0) {
                    printf("\tINT %d:\tG%d--", i, ansg[a][tmp]);
                    flag++;
                  } else {
                    printf("G%d\t%d\n", ansg[a][tmp], (int) ansig[a][b][i]);
                    flag=0;
                  }
                }
              }
            }
          }
        }
        cout << "\nTime: " << t << " seconds" << endl;
//        printf("\nNumber of Parameters chosen: %d\n", nump);
//        printf("\nTotal Number of Parameters: %d\tSample size: %d\tBIC: %d\n", ncov, n, mBIC);
//        printf("Parameters:\n");
//        for ( i = 0; i < nump; i++) {
//          printf("\tparam %d: %lf\n", i, ansp[i]);
//        }
//
////      //////free(delp);
//        //////free(ansp);
//
//        printf("\nNumber of Groups chosen: %d\n", numg);
//        printf("Groups:\n");
//        for ( i = 0; i < numg; i++) {
//          printf("\tgroup %d: %lf\n", i, ansg[i]);
//        }
//
////      //////free(delg);
//        //////free(ansg);
//
////        printf("\nBeta Hat estimates:\n");
////        print_matrix(beta);
//
      } else {

        int vartype = atoi(argv[6]);
        int sigma = atoi(argv[7]);
        int n_grp = ncov/tgrp;
        int in_grp = icov/ig;
        
        for (int a=0; a<2; a++) { // mBIC
          for (int b=0; b<2; b++) { // iBIC
            printf("-1,%d,%d,%d,", b, in_grp, numig[a][b]);
            if (numig[a][b] != 0) {
              for (i = 0; i < numig[a][b]; i++) {
                cout << name[ansig[a][b][i]][a] << ",";
              }
            }
  
            printf("%d,%d,%d,%d,%d,%d,%d,%d,", n, a, ncov, tgrp, sigma, vartype, (int) n_grp, numg[a]);
  
            for ( i = 0; i < numg[a]; i++) {
              if (i != (numg[a] -1)) {
                printf("%d,", (int) ansg[a][i]);
              } else {
                printf("%d\n", (int) ansg[a][i]);
              }
            }
  
  
            for ( i = 0; i < ncov + numip[a][b]; i++) {
              if (i != (ncov + numip[a][b] - 1)) {
                printf("%lf,", beta[a][b][i]);
              } else {
                printf("%lf\n", beta[a][b][i]);
              }
            }
            printf("%lf\n", t);
          }
        }
      }
    }

//cout << "here" << endl;
    ////free(ansip);
    ////free(ansg);
    ////free(ansig);
    ////free(ansp);

  } else {  
    int id, job, size, grp;
    int flag = 1;
    int flag2 = 0;
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
      fmat tmpX(n, size);

      // If the code received is 0 then the RSS needs to be calculated
      if ( job == 0) {

        // Using the id that was received get the correct matrix
        for (j = 0; j < size; j++) {
          tmpX.col(j) = X.col(id + j);
        }

        // Create send buffer
        sendb[0] = RSS(y, tmpX, n); // Calculate RSS based on that matrix
        sendb[1] = (double) grp; // Send back the associated subscript
        sendb[2] = (double) size; // Send back the associated subscript
        sendb[3] = (double) me; // And the node number

        // Send buffer to Node 0
        MPI_Send(sendb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

      // If code received is 1 then X needs to be edited
      } else if (job == 1) {

        // Create matrix with id that was received
        for ( i = 0; i < size; i++) {
          tmpX.col(i) = X.col(id + i) / norm(X.col(id + i), 2);
        }

//      for (i = 0; i < size; i++) {
//        X.shed_col(id);
//
//        y = y - tmpX.col(i) * (tmpX.col(i).t() * y);
//        X = X - tmpX.col(i) * (tmpX.col(i).t() * X);
//      }

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
      } else if (job == 2) { // This kills current jobs nodes are doing
        flag2 = 2;
        flag = -1;
      } else { // This kills ALL nodes immediately
        flag2 = -1;
        flag = -1;
      }
      tmpX.reset();
    }

    X.reset();
//printf("flag2=%d\n", flag2);
    while (flag2 == 2) { // The node is waiting for a message
  
      MPI_Recv(recvib, 2500, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
      flag = recvib[0];
      int nump = recvib[1]; // Number of main parameters
      int numg = recvib[2]; // Number of main groups
      int icov = recvib[3]; // Number of POSSIBLE interaction parameters
  
      fmat mMAT(n, nump);
  
      fmat iMAT(n, nump + icov);
      
      if (flag != -1){
        ivec size(numg);
        fvec newCnum(numg);
        newCnum(0)=0;
        size(0)=nINgrp_org(recvib[4]);
        int tmp=0;
        for ( i = 1; i < numg; i++) {
          size(i) = nINgrp_org(recvib[i+4]);
          tmp=size(i-1);
          newCnum(i)=newCnum(i-1)+tmp;
        }
  //      for ( i = 0; i < numg; i++) {
  //printf("size[%d] = %d\n", i, (int) size[i]);
  //printf("newCnum[%d] = %d\n", i, (int) newCnum[i]);
  //      }
  
        for ( i = 0; i < nump; i++) {
          mMAT.col(i) = X0.col(recvib[i+numg+4]);
        }
  
        for ( i = 0; i < nump; i++) {
          iMAT.col(i) = mMAT.col(i);
        }
  
  //      X0.reset();
  
        int cntp=0;
        int ig=0;
  
        for ( k = 0; k < numg; k++) {
          for ( int l = k+1; l < numg; l++) {
            int tmp=size(k)*size(l);
            fmat test(n, tmp);
  
            int intcnt=0;
            int iflag=0;
  
            for ( i = (int) newCnum(k); i < (int) newCnum(k) + (int) size(k); i++) {
              for ( j = (int) newCnum(l); j < (int) newCnum(l) + (int) size(l); j++) {
                test.col(intcnt) = mMAT.col(i) % mMAT.col(j);
                if (sum(test.col(intcnt))==0) {
                  iflag=1;
                  break;
                }
                intcnt++;
              }
            }
            if (iflag == 0) {
              for (int a = 0; a < tmp; a++) {
                iMAT.col(cntp+nump) = test.col(a);
                cntp++;
              }
              ig++;
            }
            test.reset();
          }
        }
  
  
  //      for ( i = 0; i < nump; i+=size) {
  //        for (j = i; j < nump; j+=size) {
  //          if (i != j) {
  //            int tcnt=0;
  //            int tflag=0;
  //            for (int a=0; a< size; a++) {
  //              for (int b=0; b<size; b++) {
  //                test.col(tcnt) = mMAT.col(i + a) % mMAT.col(j + b);
  //                if (0 == sum(test.col(tcnt))) {
  //                  tflag=1;
  //                }
  //                tcnt++;
  //              }
  //            }
  //            if ( 0 == tflag) {
  //              for (int a = 0; a < isize; a++) {
  //                iMAT.col(cntp+nump) = test.col(a);
  //                cntp++;
  //              }
  //            }
  //          }
  //        }
  //      }
  //      test.reset();
  
        y = y0;
  //      y0.reset();
  
  if ((me == 1) && (prnt==1)) {
    for (i = 0; i < numg; i++) {
      printf("Node: %d\tGroup: %d\tsize[%d]=%d\tcNum[%d]=%d\n", me, (int) recvib[i+4], i, (int) size(i), i, (int) newCnum(i));
    }
    printf("\n(1) iMAT:\tNumber of INT groups: %d\tNumber of INT parameters: %d\n", ig, icov);
    print_matrix(iMAT);
    printf("\n(1) y\n");
    print_matrix(y);
  }
  
        for (j = 0; j < numg; j++) {
          fmat tmpX(n, size(j));
          for (i = 0; i < size(j); i++) {
            tmpX.col(i) = iMAT.col(i) / norm(iMAT.col(i), 2);
          }
  
  //        for (i = 0; i < size(j); i++) {
  //          iMAT.shed_col(0);
  //
  //          y = y - tmpX.col(i) * (tmpX.col(i).t() * y);
  //          iMAT = iMAT - tmpX.col(i) * (tmpX.col(i).t() * iMAT);
  //        }
          // Must remove those columns from the design matrix
          for (i = 0; i < size(j); i++) {
            iMAT.shed_col(0);
          }
  
          // Recalculate both y and X
          y = y - tmpX * (tmpX.t() * y);
          iMAT = iMAT - tmpX * (tmpX.t() * iMAT);
          tmpX.reset();
        }
  
  if ((me == 1) && (prnt == 1)) {
    printf("\n(2) iMAT:\n");
    print_matrix(iMAT);
    printf("\n(2) y:\n");
    print_matrix(y);
  }
  
        size.reset();
        newCnum.reset();
        MPI_Send(sendb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      } else {
        flag2=-1;
      }
  
      // Every node does the same thing.  Every node has an interaction matrix (iMAT),
      //   and for Foward Regression to work correctly iMAT must be identical on all
      //   nodes. After finding the signficant factor (each step) the columns of those
      //   significant factors are to be removed and the iMAT is recalculated.
  
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
            tmpX.col(j) = iMAT.col(id + j);
          }
  
          // Create send buffer
          sendb[0] = RSS(y, tmpX, n); // Calculate RSS based on that matrix
          sendb[1] = (double) grp; // Send back the associated subscript
          sendb[2] = (double) size; // Send back the associated subscript
          sendb[3] = (double) me; // And the node number
  //printf("n = %d\tsize = %d\tid = %d\tjob = %d\tgrp = %d\tnode = %d\tRSS = %lf\n", n, size, id, job, grp, me, sendb[0]);
  
          // Send buffer to Node 0
          MPI_Send(sendb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  
        // If code received is 1 then X needs to be edited
        } else if (job == 1) {
  
          // Create matrix with id that was received
          for ( i = 0; i < size; i++) {
            tmpX.col(i) = iMAT.col(id + i) / norm(iMAT.col(id + i), 2);
          }
  
  //        for (i = 0; i < size; i++) {
  //          iMAT.shed_col(id);
  //
  //          y = y - tmpX.col(i) * (tmpX.col(i).t() * y);
  //          iMAT = iMAT - tmpX.col(i) * (tmpX.col(i).t() * iMAT);
  //        }
  
          // Must remove those columns from the design matrix
          for (i = 0; i < size; i++) {
            iMAT.shed_col(id);
          }
  
          // Recalculate both y and X
          y = y - tmpX * (tmpX.t() * y);
          iMAT = iMAT - tmpX * (tmpX.t() * iMAT);
  
          // Send back confirmation
          MPI_Send(sendb, 4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  
        // If code is something else, kill node
        } else if (job == 2) { // This kills current jobs nodes are doing
          flag2 = 2;
          flag = -1;
        } else { // This kills ALL nodes immediately
          flag2 = -1;
          flag = -1;
        }
        tmpX.reset();
      }
    }
  }
//printf("Node: %d made it here.\n", me);
  //free(sendib);
  //free(recvib);
  //free(sendb);
  //free(recvb);
  MPI_Finalize();
  return 0;
}
