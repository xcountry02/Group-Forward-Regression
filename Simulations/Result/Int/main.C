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
 |        Input:  argv[1] -- (int) -- The total number of main groups
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



#include <cstdlib>
#include <cmath>
#include <time.h>
#include "header.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
  int mg = atoi(argv[1]); // Total main groups
  int ig = nCr(mg, 2); // Total number of POSSIBLE interaction groups
  int n = atoi(argv[2]); // Sample size
  int iter = atoi(argv[3]); // Number or iteractions per simulation
  vector<string> cSNR; // SNR values
  cSNR.push_back("3");
  cSNR.push_back("1");
  cSNR.push_back("0.33");
  float fSNR[3]; // Number SNR values (in case calculations)
  fSNR[0] = 3;
  fSNR[1] = 1;
  fSNR[2] = .33;
  int mbic = 0; // Main BIC (0 or 1) small or large
  int ibic = 0; // Interaction BIC (0 or 1) small or large
  int varind = 0; // Variance Structure (0 or 1) AR(1) or IID
  int mcov = 0; // Number of main covariates
  int mgs = 0; // Main group size
  int igs = 0; // Interaction group size
  int is2nr = 0; // SNR iterator
  float fs2nr = 0;
  int tLEN = 0; // Main true model length
  int tiLEN=0; // Interaction true model length
  
  float trueB[1000000]; // True Beta (connected directly to simulation data set)
  float trueiB[3][1000000]; // True interaction beta 
  
  int trueGRP[25]; // True Group Numbers
  int trueiGRP[50][2]; // True interaction group numbers
  
  
  for ( mgs=2; mgs<=4; mgs++) { // mgs
    for ( varind=0; varind<=1; varind++) { // variance structure
      for ( is2nr=0; is2nr<=2; is2nr++) { // SNR
        for ( mbic=0; mbic<=1; mbic++) { // main bic
          for ( ibic=0; ibic<=1; ibic++) { // interaction bic
    
//  printf("\nBefore\nmg = %d\n", mg);
//  printf("n = %d\n", n);
//  printf("mgs = %d\n", mgs);
//  printf("var = %d\n", varind);
//  printf("sigma = %d\n", sigma);
//  printf("mbic = %d\n", mbic);
//  printf("ibic = %d\n", ibic);
//  printf("iter = %d\n", iter);
      
            stringstream execNAME;
            execNAME << "./gather.sh ";
            execNAME << mg;
            execNAME << " ";
            execNAME << n;
            execNAME << " ";
            execNAME << mgs+1;
            execNAME << " ";
            execNAME << varind;
            execNAME << " ";
            execNAME << cSNR.at(is2nr);
            execNAME << " ";
            execNAME << mbic;
            execNAME << " ";
            execNAME << ibic;
            execNAME << " ";
            execNAME << iter;
      
// printf("execute gather.sh = %s\n", execNAME.str().c_str());
      
            system(execNAME.str().c_str());
      
            stringstream name;
            name << "results.";
            name << mg;
            name << ".";
            name << n;
            name << ".";
            name << mgs+1;
            name << ".";
            name << cSNR.at(is2nr);
            name << ".";
            name << mbic;
            name << ".";
            name << ibic;
            name << ".";
            name << varind;
            name << ".out";
      
//  printf("result filename = %s\n", name.str().c_str());
      
            ifstream fp;
            fp.open(name.str().c_str(), ios::in); // Open file
            string l;
      
            vector<double> all;
      
            mcov = mg*mgs;
// Create true B, must be hand put in      
            for (int a=0; a<mcov; a++) {
              trueB[a] = 0;
            }
      
            if (mgs == 1) {
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
            } else if (mgs == 2) {
              trueB[0] = 2;
              trueB[1] = 2;
     
              trueB[2] = 2;
              trueB[3] = 2;
      
              trueB[4] = 2;
              trueB[5] = 2;
      
              trueB[6] = 2;
              trueB[7] = 2;
      
//              trueB[8] = 3;
//              trueB[9] = 3;
//      
//              trueB[10] = 3;
//              trueB[11] = 3;
      
              tLEN = 4;
              trueGRP[0] = 0;
              trueGRP[1] = 1;
              trueGRP[2] = 2;
              trueGRP[3] = 3;
//              trueGRP[4] = 4;
//              trueGRP[5] = 5;
    
              trueiB[0][0] = 3;
              trueiB[0][1] = 3;
              trueiB[0][2] = 3;
              trueiB[0][3] = 3;
      
              trueiB[1][0] = 3;
              trueiB[1][1] = 3;
              trueiB[1][2] = 3;
              trueiB[1][3] = 3;

              trueiB[2][0] = 3;
              trueiB[2][1] = 3;
              trueiB[2][2] = 3;
              trueiB[2][3] = 3;

              tiLEN = 3;
    
              trueiGRP[0][0] = 0;
              trueiGRP[0][1] = 1;

              trueiGRP[1][0] = 1;
              trueiGRP[1][1] = 2;

              trueiGRP[2][0] = 2;
              trueiGRP[2][1] = 3;
             } else if (mgs == 3) {
              trueB[0] = 2;
              trueB[1] = 2;
              trueB[2] = 2;
      
              trueB[3] = 2;
              trueB[4] = 2;
              trueB[5] = 2;
      
              trueB[6] = 2;
              trueB[7] = 2;
              trueB[8] = 2;
      
              trueB[9] = 2;
              trueB[10] = 2;
              trueB[11] = 2;
      
//              trueB[12] = 3;
//              trueB[13] = 3;
//              trueB[14] = 3;
//      
//              trueB[15] = 3;
//              trueB[16] = 3;
//              trueB[17] = 3;
      
              tLEN = 4;
              trueGRP[0] = 0;
              trueGRP[1] = 1;
              trueGRP[2] = 2;
              trueGRP[3] = 3;
//              trueGRP[4] = 4;
//              trueGRP[5] = 5;
    
              trueiB[0][0] = 3;
              trueiB[0][1] = 3;
              trueiB[0][2] = 3;
              trueiB[0][3] = 3;  
              trueiB[0][4] = 3;
              trueiB[0][5] = 3;
              trueiB[0][6] = 3;
              trueiB[0][7] = 3;
              trueiB[0][8] = 3;
    
              trueiB[1][0] = 3;
              trueiB[1][1] = 3;
              trueiB[1][2] = 3;
              trueiB[1][3] = 3;  
              trueiB[1][4] = 3;
              trueiB[1][5] = 3;
              trueiB[1][6] = 3;
              trueiB[1][7] = 3;
              trueiB[1][8] = 3;
    
              trueiB[2][0] = 3;
              trueiB[2][1] = 3;
              trueiB[2][2] = 3;
              trueiB[2][3] = 3;  
              trueiB[2][4] = 3;
              trueiB[2][5] = 3;
              trueiB[2][6] = 3;
              trueiB[2][7] = 3;
              trueiB[2][8] = 3;
      
              tiLEN = 3;
    
              trueiGRP[0][0] = 0;
              trueiGRP[0][1] = 1;

              trueiGRP[1][0] = 1;
              trueiGRP[1][1] = 2;

              trueiGRP[2][0] = 2;
              trueiGRP[2][1] = 3;
            } else if (mgs == 4) {
              trueB[0] = 2;
              trueB[1] = 2;
              trueB[2] = 2;
              trueB[3] = 2;
      
              trueB[4] = 2;
              trueB[5] = 2;
              trueB[6] = 2;
              trueB[7] = 2;
      
              trueB[8] = 2;
              trueB[9] = 2;
              trueB[10] = 2;
              trueB[11] = 2;
      
              trueB[12] = 2;
              trueB[13] = 2;
              trueB[14] = 2;
              trueB[15] = 2;
      
              tLEN = 4;
              trueGRP[0] = 0;
              trueGRP[1] = 1;
              trueGRP[2] = 2;
              trueGRP[3] = 3;
//              trueGRP[4] = 4;
//              trueGRP[5] = 5;
    
              trueiB[0][0] = 3;
              trueiB[0][1] = 3;
              trueiB[0][2] = 3;
              trueiB[0][3] = 3;  
              trueiB[0][4] = 3;
              trueiB[0][5] = 3;
              trueiB[0][6] = 3;
              trueiB[0][7] = 3;
              trueiB[0][8] = 3;
              trueiB[0][9] = 3;
              trueiB[0][10] = 3;
              trueiB[0][11] = 3;
              trueiB[0][12] = 3;  
              trueiB[0][13] = 3;
              trueiB[0][14] = 3;
              trueiB[0][15] = 3;
      
              trueiB[1][0] = 3;
              trueiB[1][1] = 3;
              trueiB[1][2] = 3;
              trueiB[1][3] = 3;  
              trueiB[1][4] = 3;
              trueiB[1][5] = 3;
              trueiB[1][6] = 3;
              trueiB[1][7] = 3;
              trueiB[1][8] = 3;
              trueiB[1][9] = 3;
              trueiB[1][10] = 3;
              trueiB[1][11] = 3;
              trueiB[1][12] = 3;  
              trueiB[1][13] = 3;
              trueiB[1][14] = 3;
              trueiB[1][15] = 3;
      
              trueiB[2][0] = 3;
              trueiB[2][1] = 3;
              trueiB[2][2] = 3;
              trueiB[2][3] = 3;  
              trueiB[2][4] = 3;
              trueiB[2][5] = 3;
              trueiB[2][6] = 3;
              trueiB[2][7] = 3;
              trueiB[2][8] = 3;
              trueiB[2][9] = 3;
              trueiB[2][10] = 3;
              trueiB[2][11] = 3;
              trueiB[2][12] = 3;  
              trueiB[2][13] = 3;
              trueiB[2][14] = 3;
              trueiB[2][15] = 3;
      
              tiLEN = 3;
    
              trueiGRP[0][0] = 0;
              trueiGRP[0][1] = 1;

              trueiGRP[1][0] = 1;
              trueiGRP[1][1] = 2;

              trueiGRP[2][0] = 2;
              trueiGRP[2][1] = 3;
            }
// These values are the ones kept and stored, they are the summary
// results, explained in paper
            int cnt = 0;
            int ilen = 0;
            int mlen = 0;
            int Size[100];
            float Cov[100];
            float Inc0[100];
            float Ext[100];
            float Cor0[100];
            int iSize[100];
            float iCov[100];
            float iInc0[100];
            float iExt[100];
            float iCor0[100];
            float arrMSE[100];
            float T[100];
            
            int N=0;
            int estGRP[500];
            int numGRP = 0;
            int estiGRP[500][2];
            int numiGRP = 0;
            int test=0;
      
            float avgSize = 0;
            float avgCov = 0;
            float avgInc0 = 0;
            float avgExt = 0;
            float avgCor0 = 0;
            float avgiSize = 0;
            float avgiCov = 0;
            float avgiInc0 = 0;
            float avgiExt = 0;
            float avgiCor0 = 0;
            float avgMSE = 0;
            float avgT = 0;
            float MSE = 0;
// Open the results.* file.  Pick off line by line, summarize results
            while(!fp.eof())
            {
              getline(fp, l); // First line is names, just read and count
              stringstream line; line<<l;
              istringstream ss(l);
              string token;
// The results from the analysis were set in a predictable manner,
// that way each line can be read, and the results recorded.      
              if ((cnt % 3) == 0) {
                getline(ss, token, ',');
                test = atoi(token.c_str());
      
                if ((cnt > 0) && (test == 0)) {
                  break;
                }
     
                if (test == -1) {
    
                  getline(ss, token, ','); // Interaction BIC
                  ibic = atoi(token.c_str());
              
                  getline(ss, token, ','); // Interaction group size
                  igs = atoi(token.c_str());

                  getline(ss, token, ','); // Number of significant interactions
                  ilen = atoi(token.c_str());
      
                  if (ilen != 0) {
                    iSize[N] = ilen;
                    for (int a=0; a < ilen; a++) {
                      for (int b=0; b < 2; b++) {
                        getline(ss, token, ',');
                        estiGRP[a][b] = atoi(token.c_str());
                      }
                    }
// Compare the results from the file to the truth.  Then make appropriate
// comparisons based on that.
                    for (int a=0; a<tiLEN; a++) {
                      for (int b=0; b<ilen; b++) {
                        if (((trueiGRP[a][0] == estiGRP[b][0]) && 
                             (trueiGRP[a][1] == estiGRP[b][1])) || 
                            ((trueiGRP[a][1] == estiGRP[b][0]) && 
                             (trueiGRP[a][0] == estiGRP[b][1]))){
                          numiGRP++;
                        }
                      }
                    }
        
                    if (numiGRP == tiLEN) {
                      iCov[N] = 1;
                      iInc0[N] = 0;
                      if (ilen == tiLEN) {
                        iExt[N] = 1;
                        iCor0[N] = 1;
                      } else {
                        iExt[N] = 0;
                        float tmp1 = ig - tiLEN;
                        float tmp2 = ilen - tiLEN;
                        float tmpTOP = tmp1 - tmp2;
                        iCor0[N] = tmpTOP / tmp1;
                      }
                    } else {
                      float top = tiLEN - numiGRP;
                      iInc0[N]= top / tiLEN;
                      iCov[N] = 0;
                      iExt[N]=0;
                      float tmp1 = ig - tiLEN;
                      float tmp2 = ilen - numiGRP;
                      float tmpTOP = tmp1 - tmp2;
                      iCor0[N] = tmpTOP / tmp1;
                    }
        
                    numiGRP = 0;
                    avgiSize = avgiSize + iSize[N];
                    avgiCov = avgiCov + iCov[N];
                    avgiCor0 = avgiCor0 + iCor0[N];
                    avgiInc0 = avgiInc0 + iInc0[N];
                    avgiExt = avgiExt + iExt[N];
                  } else {
                    numiGRP = 0;
                    iSize[N] = 0;
                    iCov[N] = 0;
                    iCor0[N] = 1;
                    iInc0[N] = 1;
                    iExt[N] = 0;
                    avgiSize = avgiSize + iSize[N];
                    avgiCov = avgiCov + iCov[N];
                    avgiCor0 = avgiCor0 + iCor0[N];
                    avgiInc0 = avgiInc0 + iInc0[N];
                    avgiExt = avgiExt + iExt[N];
                  }
              
                  getline(ss, token, ','); // Sample mgs
                  n = atoi(token.c_str());
      
                } else {
                  n = test;
                }
    
                getline(ss, token, ','); // Main BIC
                mbic = atoi(token.c_str());
               
                getline(ss, token, ','); // Total number of main covariates
                mcov = atoi(token.c_str());
            
                getline(ss, token, ','); // Total number of groups
                mg = atoi(token.c_str());
            
                getline(ss, token, ','); // Regression error
                fs2nr = atof(token.c_str());
            
                getline(ss, token, ','); // Variance type (IID or AR1)
                varind = atoi(token.c_str());
            
                getline(ss, token, ','); // Main group mgs
                mgs = atoi(token.c_str());
                mgs = mgs - 1;    

//  printf("\nmg = %d\n", mg);
//  printf("n = %d\n", n);
//  printf("mgs = %d\n", mgs);
//  printf("var = %d\n", varind);
//  printf("sigma = %d\n", sigma);
//  printf("mbic = %d\n", mbic);
//  printf("ibic = %d\n", ibic);
//  printf("iter = %d\n", iter);
//  printf("cnt = %d\n", cnt);
//  printf("ilen = %d\n", ilen);               
//  printf("igs = %d\n", igs);                            
//  printf("ig = %d\n", ig);               
//  printf("ibic = %d\n", ibic);               
            
                getline(ss, token, ','); // Number of main sig. parameters
                mlen = atoi(token.c_str());
                Size[N] = mlen;
                for (int a=0; a < mlen; a++) {
                  getline(ss, token, ',');
                  estGRP[a] = atoi(token.c_str());
                }
            
                for (int a=0; a<tLEN; a++) {
                  for (int b=0; b<mlen; b++) {
                    if (trueGRP[a] == estGRP[b]) {
                      numGRP++;
                    }
                  }
                }
    
                if (numGRP == tLEN) {
                  Cov[N] =1;
                  Inc0[N]=0;
                  if (mlen == tLEN) {
                    Ext[N] = 1;
                    Cor0[N]=1;
                  } else {
                    Ext[N]=0;
                    float tmp1 = mg - tLEN;
                    float tmp2 = mlen - tLEN;
                    float tmpTOP = tmp1 - tmp2;
                    Cor0[N] = tmpTOP / tmp1;
                  }
                } else {
                  float top = tLEN - numGRP;
                  Inc0[N]= top / tLEN;
                  Cov[N] = 0;
                  Ext[N]=0;
                  float tmp1 = mg - tLEN;
                  float tmp2 = mlen - numGRP;
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
                int Icnt=0;
                int finish=0;
                while(getline(ss, token, ',')) {
                  finish=0;
                  if (TMPcnt < mcov) {
                    tmp = atof(token.c_str()) - trueB[TMPcnt++];
                    tmp = pow(tmp, 2);
                    MSE = MSE + tmp;
                  } else {
                    for (int c=0; c<tiLEN; c++) {
                      int flag=0;
                      float diff=0;
                      if (finish == 0){
                        for (int b=0; b<2; b++) {
                          for (int d=0; d<2; d++) {
                            if (estiGRP[Icnt][b] == trueiGRP[c][d]) {
                              flag=flag + 1;
                              if (flag == 2) {
                                for (int a=0; a<igs; a++) {
                                  diff = atof(token.c_str()) - trueiB[c][a];
                                  diff = pow(diff,2);
                                  MSE = MSE + diff;
                                  if (a != (igs - 1)) {
                                    getline(ss, token, ',');
                                  }
                                }
                                finish=1;
                              } 
                            }
                          }
                        }
                            if ((c == (tiLEN - 1)) && (finish == 0)){
                              for (int a=0; a<igs; a++) {
                                diff = atof(token.c_str());
                                diff = pow(diff,2);
                                MSE = MSE + diff;
                                if (a != (igs - 1)) {
                                  getline(ss, token, ',');
                                }
                              }
                              finish=1;
                            }
                           
                      } else {
                        break;
                      }
                    }
                    Icnt++;
                  }
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
//  printf("\navgiSize = %lf\n", avgiSize);
//  printf("avgiCov = %lf\n", avgiCov);
//  printf("avgiCor0 = %lf\n", avgiCor0);
//  printf("avgiInc0 = %lf\n", avgiInc0);
//  printf("avgiExt = %lf\n", avgiExt);
//  printf("avgSize = %lf\n", avgSize);
//  printf("avgCov = %lf\n", avgCov);
//  printf("avgCor0 = %lf\n", avgCor0);
//  printf("avgInc0 = %lf\n", avgInc0);
//  printf("avgExt = %lf\n", avgExt);               
//  printf("avgMSE = %lf\n", avgMSE);                            
//  printf("avgT = %lf\n", avgT);               
              cnt++;
            }
            
            fp.close();
            avgiSize = avgiSize / N;
            avgiCov = avgiCov / N;
            avgiCor0 = avgiCor0 / N;
            avgiInc0 = avgiInc0 / N;
            avgiExt = avgiExt / N;
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
            float sdiSize = 0;
            float sdiCov = 0;
            float sdiInc0 = 0;
            float sdiExt = 0;
            float sdiCor0 = 0;
            
            float tmp1=0;
            float tmp2=0;
            float tmp3=0;
            float tmp4=0;
            float tmp5=0;
            float tmp6=0;
            float tmp7=0;
            float tmp8=0;
            float tmp9=0;
            float tmp10=0;
            float tmp11=0;
            float tmp12=0;
            
            for (int a=0; a<N; a++) {
              tmp1=Cov[a] - avgCov;
              sdCov=sdCov + pow(tmp1, 2);
              tmp2=Cor0[a] - avgCor0;
              sdCor0=sdCor0 + pow(tmp2, 2);
              tmp3=Inc0[a] - avgInc0;
              sdInc0=sdInc0 + pow(tmp3, 2);
              tmp4=Ext[a] - avgExt;
              sdExt=sdExt + pow(tmp4,2);
              tmp5=Size[a] - avgSize;
              sdSize=sdSize + pow(tmp5,2);
              tmp6=iCov[a] - avgiCov;
              sdiCov=sdiCov + pow(tmp6, 2);
              tmp7=iCor0[a] - avgiCor0;
              sdiCor0=sdiCor0 + pow(tmp7, 2);
              tmp8=iInc0[a] - avgiInc0;
              sdiInc0=sdiInc0 + pow(tmp8, 2);
              tmp9=iExt[a] - avgiExt;
              sdiExt=sdiExt + pow(tmp9,2);
              tmp10=iSize[a] - avgiSize;
              sdiSize=sdiSize + pow(tmp10,2);
              tmp11=T[a] - avgT;
              sdT=sdT + pow(tmp11,2);
              tmp12=arrMSE[a] - avgMSE;
              sdMSE=sdMSE + pow(tmp12,2);
              tmp1=0;
              tmp2=0;
              tmp3=0;
              tmp4=0;
              tmp5=0;
              tmp6=0;
              tmp7=0;
              tmp8=0;
              tmp9=0;
              tmp10=0;
              tmp11=0;
              tmp12=0;
            }
            
            sdSize = sdSize / N;
            sdCov = sdCov / N;
            sdCor0 = sdCor0 / N;
            sdInc0 = sdInc0 / N;
            sdExt = sdExt / N;
            sdiSize = sdiSize / N;
            sdiCov = sdiCov / N;
            sdiCor0 = sdiCor0 / N;
            sdiInc0 = sdiInc0 / N;
            sdiExt = sdiExt / N;
            sdiSize = sqrt(sdiSize);
            sdiCov = sqrt(sdiCov);
            sdiCor0 = sqrt(sdiCor0);
            sdiInc0 = sqrt(sdiInc0);
            sdiExt = sqrt(sdiExt);
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
           
// Create result card based on that information gathered.
 
            nam << "summary_result.";
            nam << mg;
            nam << ".";
            nam << n;
            nam << ".";
            nam << mgs + 1;
            nam << ".";
            nam << is2nr + 1;
            nam << ".";
            nam << mbic;
            nam << ".";
            nam << ibic;
            nam << ".";
            nam << varind;
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
            ofp << "Total number of main Groups = " << mg << endl;
            ofp << "Total number of main Parameters = " << mcov << endl;
            ofp << "Main Group Size = " << mgs + 1 << endl;
            ofp << "Sample Size = " << n << endl;
            ofp << "Iterations = " << iter << endl;
            if (varind == 0){
              ofp << "Variance Structure = AR1" << endl;
            } else {
              ofp << "Variance Structure = IID" << endl;
            }
            ofp << "SNR = " << cSNR[is2nr] << endl;      
            ofp << "Coverage probability = " << avgCov << " (" << sdCov << ")" << endl;
            ofp << "Percentage of correct zeros = " << avgCor0 << " (" << sdCor0 << ")" << endl;
            ofp << "Percentage of incorrect zeros = " << avgInc0 << " (" << sdInc0 << ")" << endl;
            ofp << "Exact select probability = " << avgExt << " (" << sdExt << ")" << endl;
            ofp << "Model Size = " << avgSize << " (" << sdSize << ")" << endl;
            ofp << "Interaction Coverage probability = " << avgiCov << " (" << sdiCov << ")" << endl;
            ofp << "Interaction Percentage of correct zeros = " << avgiCor0 << " (" << sdiCor0 << ")" << endl;
            ofp << "Interaction Percentage of incorrect zeros = " << avgiInc0 << " (" << sdiInc0 << ")" << endl;
            ofp << "Interaction Exact select probability = " << avgiExt << " (" << sdiExt << ")" << endl;
            ofp << "Interaction Model Size = " << avgiSize << " (" << sdiSize << ")" << endl;
            ofp << "MSE = " << avgMSE << " (" << sdMSE << ")" << endl;
            ofp << "Total time = " << avgT << " (" << sdT << ")" << endl;
            
            ofp.close();
          }
        }
      }
    }
  }
  return 0;
}
