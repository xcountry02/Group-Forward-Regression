/*=============================================================================
 |      Program:  main.cpp
 |
 |       Author:  K. Michels
 |     Language:  C++ (g++ compiler)
 |   To Compile:  g++ -Wall -g  -I /home/u2/kamichels/gFRmpi/arma/usr/local/
 |                include/ -L /home/u2/kamichels/gFRmpi/lib/ -L /home/u2/kamichels/
 |                gFRmpi/arma/usr/local/lib64/ -larmadillo  main.C -o ans (it is 
 |                assumed that the stack header file, header.h, is in the same 
 |                directory.)
 |
 +-----------------------------------------------------------------------------
 |
 |  Description:  Taking the data simulations, this program takes those
 |                results and quantifies it in summary statistics.
 |
 |        Input:  argv[1] -- (int) -- The number of total groups
 |                argv[2] -- (int) -- Sample size
 |                argv[3] -- (int) -- The number of iterations on simulations
 |
 |       Output:  The output are many files that contain all of the necessary
 |                information for the simulation.  These files have a very
 |                specific format and they are referred to as result cards.
 |
 |   Known Bugs:  None; all operations work as they should.
 |
 *===========================================================================*/

// Note:  The 

#include <cstdlib>
#include <cmath>
#include <time.h>
#include "header.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
int p0 = atoi(argv[1]);
int n = atoi(argv[2]);
int iter = atoi(argv[3]);
int ngrp = 0;
int sigma=0;
int bic = 0;
int varind = 0;
int total = 0;
int size = 0;

int tLEN = 0;

float trueB[1000000];

int trueGRP[12];

for ( ngrp=3; ngrp<=5; ngrp++) { // ngrp
  for ( varind=0; varind<=1; varind++) { // variance structure
    for ( sigma=1; sigma<=3; sigma++) { // sigma
      for ( bic=0; bic<=1; bic++) { // bic

//   printf("p0 = %d\n", p0);
//   printf("n = %d\n", n);
//   printf("ngrp = %d\n", ngrp);
//   printf("var = %d\n", varind);
//   printf("sigma = %d\n", sigma);
//   printf("bic = %d\n", bic);
//   printf("iter = %d\n", iter);
  
        stringstream execNAME;
        execNAME << "./gather.sh ";
        execNAME << p0;
        execNAME << " ";
        execNAME << n;
        execNAME << " ";
        execNAME << ngrp;
        execNAME << " ";
        execNAME << varind;
        execNAME << " ";
        execNAME << sigma;
        execNAME << " ";
        execNAME << bic;
        execNAME << " ";
        execNAME << iter;
  
  // printf("execute gather.sh = %s\n", execNAME.str().c_str());
  
        system(execNAME.str().c_str());
  
        stringstream name;
        name << "results.";
        name << p0;
        name << ".";
        name << n;
        name << ".";
        name << ngrp;
        name << ".";
        name << sigma;
        name << ".";
        name << bic;
        name << ".";
        name << varind;
        name << ".out";
  
  // printf("result filename = %s\n", name.str().c_str());
  
        ifstream fp;
        fp.open(name.str().c_str(), ios::in); // Open file
        string l;
  
        vector<double> all;
  
        size = ngrp - 1;
        total = p0*size;
  
        for (int a=0; a<total; a++) {
          trueB[a] = 0;
        }
  
        if (ngrp == 2) {
          trueB[0] = 3;
          trueB[1] = 3;
          trueB[2] = 3;
          trueB[3] = 3;
          trueB[4] = 3;
          trueB[5] = 3;
          trueB[6] = 2;
          trueB[7] = 2;
          trueB[8] = 2;
          trueB[9] = 2;
          trueB[10] = 2;
          trueB[11] = 2;
  
          tLEN = 12;
          trueGRP[0] = 0;
          trueGRP[1] = 1;
          trueGRP[2] = 2;
          trueGRP[3] = 3;
          trueGRP[4] = 4;
          trueGRP[5] = 5;
          trueGRP[6] = 6;
          trueGRP[7] = 7;
          trueGRP[8] = 8;
          trueGRP[9] = 9;
          trueGRP[10] = 10;
          trueGRP[11] = 11;
        } else if (ngrp == 3) {
          trueB[0] = 3;
          trueB[1] = 3;
  
          trueB[2] = 3;
          trueB[3] = 3;
  
          trueB[4] = 3;
          trueB[5] = 3;
  
          trueB[6] = 3;
          trueB[7] = 3;
  
          trueB[8] = 3;
          trueB[9] = 3;
  
          trueB[10] = 3;
          trueB[11] = 3;
  
          tLEN = 6;
          trueGRP[0] = 0;
          trueGRP[1] = 1;
          trueGRP[2] = 2;
          trueGRP[3] = 3;
          trueGRP[4] = 4;
          trueGRP[5] = 5;
         } else if (ngrp == 4) {
          trueB[0] = 3;
          trueB[1] = 3;
          trueB[2] = 3;
  
          trueB[3] = 3;
          trueB[4] = 3;
          trueB[5] = 3;
  
          trueB[6] = 3;
          trueB[7] = 3;
          trueB[8] = 3;
  
          trueB[9] = 3;
          trueB[10] = 3;
          trueB[11] = 3;
  
          tLEN = 4;
          trueGRP[0] = 0;
          trueGRP[1] = 1;
          trueGRP[2] = 2;
          trueGRP[3] = 3;
        } else if (ngrp == 5) {
          trueB[0] = 3;
          trueB[1] = 3;
          trueB[2] = 3;
          trueB[3] = 3;
  
          trueB[4] = 3;
          trueB[5] = 3;
          trueB[6] = 3;
          trueB[7] = 3;
  
          trueB[8] = 3;
          trueB[9] = 3;
          trueB[10] = 3;
          trueB[11] = 3;
  
          tLEN = 3;
          trueGRP[0] = 0;
          trueGRP[1] = 1;
          trueGRP[2] = 2;
        }
  
        int cnt = 0;
        int eLEN = 0;
        int Size[100];
        float Cov[100];
        float Inc0[100];
        float Ext[100];
        float Cor0[100];
        float arrMSE[100];
        float T[100];
        
        int N=0;
        int estGRP[500];
        int numGRP = 0;
        int test=0;
  
        float avgSize = 0;
        float avgCov = 0;
        float avgInc0 = 0;
        float avgExt = 0;
        float avgCor0 = 0;
        float avgMSE = 0;
        float avgT = 0;
        float MSE = 0;
        //vector<string> Xlist;
        while(!fp.eof())
        {
          getline(fp, l); // First line is names, just read and count
          stringstream line; line<<l;
          istringstream ss(l);
          string token;
  
          if ((cnt % 3) == 0) {
            getline(ss, token, ',');
            test = atoi(token.c_str());
  
            if ((cnt > 0) && (test == 0)) {
              break;
            }
  
            n = test;
  
            getline(ss, token, ',');
            bic = atoi(token.c_str());
           
            getline(ss, token, ',');
            total = atoi(token.c_str());
        
            getline(ss, token, ',');
            p0 = atoi(token.c_str());
        
            getline(ss, token, ',');
            sigma = atoi(token.c_str());
        
            getline(ss, token, ',');
            varind = atoi(token.c_str());
        
            getline(ss, token, ',');
            ngrp = atoi(token.c_str());

//   printf("p0 = %d\n", p0);
//   printf("n = %d\n", n);
//   printf("ngrp = %d\n", ngrp);
//   printf("var = %d\n", varind);
//   printf("sigma = %d\n", sigma);
//   printf("bic = %d\n", bic);
//   printf("iter = %d\n", iter);
//   printf("cnt = %d\n", cnt);
        
            getline(ss, token, ',');
            eLEN = atoi(token.c_str());
            Size[N] = eLEN;
            for (int a=0; a < eLEN; a++) {
              getline(ss, token, ',');
              estGRP[a] = atoi(token.c_str());
            }
        
            for (int a=0; a<tLEN; a++) {
              for (int b=0; b<eLEN; b++) {
                if (trueGRP[a] == estGRP[b]) {
                  numGRP++;
                }
              }
            }

            if (numGRP == tLEN) {
              Cov[N] =1;
              Inc0[N]=0;
              if (eLEN == tLEN) {
                Ext[N] = 1;
                Cor0[N]=1;
              } else {
                Ext[N]=0;
                float tmp1 = p0 - tLEN;
                float tmp2 = eLEN - tLEN;
                float tmpTOP = tmp1 - tmp2;
                Cor0[N] = tmpTOP / tmp1;
              }
            } else {
              float top = tLEN - numGRP;
              Inc0[N]= top / tLEN;
              Cov[N] = 0;
              Ext[N]=0;
              float tmp1 = p0 - tLEN;
              float tmp2 = eLEN - numGRP;
              float tmpTOP = tmp1 - tmp2;
              Cor0[N] = tmpTOP / tmp1;
            }

            numGRP = 0;
            avgSize = avgSize + Size[N];
            avgCov = avgCov + Cov[N];
            avgCor0 = avgCor0 + Cor0[N];
            avgInc0 = avgInc0 + Inc0[N];
            avgExt = avgExt + Ext[N];

          } else if ((cnt % 3) == 1) {
  
            float tmp = 0;
            int TMPcnt = 0;
            while(getline(ss, token, ',')) {
               tmp = atof(token.c_str()) - trueB[TMPcnt++];
               tmp = pow(tmp, 2);
               MSE = MSE + tmp;
            }
            arrMSE[N] = MSE;
            MSE = 0;
            avgMSE = avgMSE + arrMSE[N];
          } else {
            getline(ss, token, ',');
            T[N] = atof(token.c_str());
            avgT = avgT + T[N];
            N++;
          }
          cnt++;
        }
        
        fp.close();
        
        avgSize = avgSize / N;
        avgCov = avgCov / N;
        avgCor0 = avgCor0 / N;
        avgInc0 = avgInc0 / N;
        avgExt = avgExt / N;
        avgMSE = avgMSE / N;
        avgT = avgT / N;

        float sdSize = 0;
        float sdCov = 0;
        float sdInc0 = 0;
        float sdExt = 0;
        float sdCor0 = 0;
        float sdMSE = 0;
        float sdT = 0;
        
        float tmp1=0;
        float tmp2=0;
        float tmp3=0;
        float tmp4=0;
        float tmp5=0;
        float tmp6=0;
        float tmp7=0;
        
        for (int a=0; a<N; a++) {
          tmp1=Cov[a] - avgCov;
          sdCov=sdCov + pow(tmp1, 2);
          tmp2=Cor0[a] - avgCor0;
          sdCor0=sdCor0 + pow(tmp2, 2);
          tmp3=Inc0[a] - avgInc0;
          sdInc0=sdInc0 + pow(tmp3, 2);
          tmp4=Ext[a] - avgExt;
          sdExt=sdExt + pow(tmp4,2);
          tmp5=T[a] - avgT;
          sdT=sdT + pow(tmp5,2);
          tmp6=arrMSE[a] - avgMSE;
          sdMSE=sdMSE + pow(tmp6,2);
          tmp7=Size[a] - avgSize;
          sdSize=sdSize + pow(tmp7,2);
          tmp1=0;
          tmp2=0;
          tmp3=0;
          tmp4=0;
          tmp5=0;
          tmp6=0;
          tmp7=0;
        }
        
        sdSize = sdSize / N;
        sdCov = sdCov / N;
        sdCor0 = sdCor0 / N;
        sdInc0 = sdInc0 / N;
        sdExt = sdExt / N;
        sdMSE = sdMSE / N;
        sdT = sdT / N;
        
        sdSize = sqrt(sdSize);
        sdCov = sqrt(sdCov);
        sdCor0 = sqrt(sdCor0);
        sdInc0 = sqrt(sdInc0);
        sdExt = sqrt(sdExt);
        sdMSE = sqrt(sdMSE);
        sdT = sqrt(sdT);
        
        ofstream ofp;
        stringstream nam;
        
        nam << "summary_result.";
        nam << p0;
        nam << ".";
        nam << n;
        nam << ".";
        nam << ngrp;
        nam << ".";
        nam << sigma;
        nam << ".";
        nam << bic;
        nam << ".";
        nam << varind;
        nam << ".out";
        
        ofp.open(nam.str().c_str(), ios::out);
        
        if (bic == 1){
          ofp << "BIC = Small BIC calculation" << endl;
        } else {
          ofp << "BIC = Large BIC calculation" << endl;
        }
        ofp << "Total number of Groups = " << p0 << endl;
        ofp << "Total number of Parameters = " << total << endl;
        ofp << "Group Size = " << ngrp << endl;
        ofp << "Sample Size = " << n << endl;
        ofp << "Iterations = " << iter << endl;
        if (varind == 0){
          ofp << "Variance Structure = AR1" << endl;
        } else {
          ofp << "Variance Structure = IID" << endl;
        }
        ofp << "Sigma = " << sigma << endl;      
        ofp << "Coverage probability = " << avgCov << " (" << sdCov << ")" << endl;
        ofp << "Percentage of correct zeros = " << avgCor0 << " (" << sdCor0 << ")" << endl;
        ofp << "Percentage of incorrect zeros = " << avgInc0 << " (" << sdInc0 << ")" << endl;
        ofp << "Exact select probability = " << avgExt << " (" << sdExt << ")" << endl;
        ofp << "Model Size = " << avgSize << " (" << sdSize << ")" << endl;
        ofp << "MSE = " << avgMSE << " (" << sdMSE << ")" << endl;
        ofp << "Total time = " << avgT << " (" << sdT << ")" << endl;
        
        ofp.close();
      }
    }
  }
}
return 0;
}
