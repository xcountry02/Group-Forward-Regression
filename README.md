# Group-Forward-Regression
Group variable screening for ultra-high dimensional data

*INTRODUTION

This program contains multiple execution scripts.  This program is the dissertation of Kurt Michels.  The main programs are in the folders Main/ and Int/.  The Main/ folder contains the group forward regression model using MPI, to find the most significant groups in for ultra-high-dimensional data.  The Int/ folder contains the group forward regression interaction model using MPI.  This finds the most significant groups AND interactions in the ultra-high dimensional data.

The Simulation/ folder contains Data/ which will create simulated data sets based upon 6 parameters.  It also conatians Result/ which gathers and reports results from the analysis done on the simulated data.  The folder PBS/ contains all PBS scripts to make life easier.

The folder also contains an armadillo tarball.  The tar ball won't work right out of the box, that is why the folder armadillo/ because from that folder one can set the correct requirements (LAPACK) to install Armadillo.  Then Armadillo can be installed to a folder, using a non-root way to install it.

*REQUIREMENTS

In order to compile the model one needs access to MPI, and also the linear library Armadillo.
