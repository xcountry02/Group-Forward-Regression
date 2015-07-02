/*=============================================================================
 |      Program:  main.cpp
 |
 |       Author:  K. Michels
 |     Language:  C++ (g++ compiler)
 |   To Compile:  g++ -c -Wall -g  -I /home/u2/kamichels/gFRmpi/arma/usr/local/
 |                include/ -L /home/u2/kamichels/gFRmpi/lib/ -L /home/u2/kamichels/
 |                gFRmpi/arma/usr/local/lib64/ -larmadillo main.cpp (it is assumed 
 |                that the stack header file, header.h, is in the same directory.)
 |
 +-----------------------------------------------------------------------------
 |
 |  Description:  Given the correct inputs this function 'prints' the response
 |                vector and design matrix to specific simulation settings.  
 |                This function works with only a fixed group size.
 |
 |        Input:  argv[1] -- (int) -- The number of total groups
 |                argv[2] -- (int) -- Sample size
 |                argv[3] -- (int) -- Group size (3-5)
 |                argv[4] -- (float) -- Variance on samples
                          -- standard values (1, 2, 3) -- could be anything >0
 |                argv[5] -- (int) -- subscript name (integer)
 |                argv[6] -- (bool) -- 0 -- AR1 variance structure
 |                                  -- 1 -- IID variance structure
 |                argv[7] -- (bool) -- 0 -- Output Data Set
 |                                  -- 1 -- Create Oracle data
 |
 |       Output:  If argv[7]==0, then the output is two files, a response (Y)
 |                and the explanatory variables (X), which are all factors.
 |                If args[7]==1, then th output is still two files, but now
 |                the files are result cards, which have a very specific format.
 |
 |   Known Bugs:  None; all operations work as they should.
 |
 *===========================================================================*/


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

#include "header.h"

using namespace std;
using namespace arma;

int main( int argc, char *argv[] )
{

  int p0 = atoi(argv[1]);
  int n = atoi(argv[2]);
  float n_grp = atof(argv[3]);
  float sigma = atof(argv[4]);
  int name = atoi(argv[5]);
  int var_ind = atoi(argv[6]);
  int oracle = atoi(argv[7]);

  wall_clock timer;
  double t;

  timer.tic();
  print(p0, n_grp, n, sigma, name, var_ind, oracle);
  t = timer.toc();

  cout << "Time: " << t << " seconds" << endl;

  return 0;
}
