CC=g++
INC_DIR = /home/u2/kamichels/Program/arma/usr/local/include/
LIB_DIR = /home/u2/kamichels/Program/lib/
LIB_DIR1 = /home/u2/kamichels/Program/arma/usr/local/lib64/
DEBUG = -g #-ggdb3 
# OFLAG = -O2
CFLAG = -c
WFLAG = -Wall
LIB = -larmadillo
# CPP11 = -std=c99 #-std=c++11
CCFLAGS = $(CFLAG) $(WFLAG) $(DEBUG) -I $(INC_DIR) -L $(LIB_DIR) -L $(LIB_DIR1) $(LIB) #$(OFLAG)
CCFLAGS1 = $(WFLAG) $(DEBUG) -I $(INC_DIR) -L $(LIB_DIR) -L $(LIB_DIR1) $(LIB) #$(OFLAG)

all: ans 

ans: main.o stat.o
	$(CC) $(CCFLAGS1) stat.o main.o -o ans

main.o: main.C
	$(CC) $(CCFLAGS) main.C

stat.o: stat.cpp
	$(CC) $(CCFLAGS) stat.cpp

clean:
	rm -fr *.o
	rm -f ans
