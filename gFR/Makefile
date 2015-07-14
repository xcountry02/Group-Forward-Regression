CC=mpic++
INC_DIR = /home/u2/kamichels/Program/arma/usr/local/include/
LIB_DIR = /home/u2/kamichels/Program/lib/
LIB_DIR1 = /home/u2/kamichels/Program/arma/usr/local/lib64/
DEBUG = -g #-ggdb3 
CFLAG = -c
WFLAG = -Wall
LIB = -larmadillo
CCFLAGS = $(CFLAG) $(WFLAG) $(DEBUG) -I $(INC_DIR) -L $(LIB_DIR) -L $(LIB_DIR1) $(LIB) #$(OFLAG)
CCFLAGS1 = $(WFLAG) $(DEBUG) -I $(INC_DIR) -L $(LIB_DIR) -L $(LIB_DIR1) $(LIB) #$(OFLAG)

all: gFRmpi

gFRmpi: misc.o gFRmpi.o
	$(CC) $(CCFLAGS1) misc.o gFRmpi.o -o gFRmpi

gFRmpi.o: gFRmpi.cpp
	$(CC) $(CCFLAGS) gFRmpi.cpp

misc.o: misc.cpp
	$(CC) $(CCFLAGS) misc.cpp

clean:
	rm -fr *.o
	rm -f gFRmpi
