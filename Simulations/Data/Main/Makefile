CC=g++
INC_DIR = /home/u2/kamichels/Program/arma/usr/local/include/
LIB_DIR = /home/u2/kamichels/Program/lib/
LIB_DIR1 = /home/u2/kamichels/Program/arma/usr/local/lib64/
DEBUG = -g #-ggdb3 
CFLAG = -c
WFLAG = -Wall
LIB = -larmadillo
CCFLAGS = $(CFLAG) $(WFLAG) $(DEBUG) -I $(INC_DIR) -L $(LIB_DIR) -L $(LIB_DIR1) $(LIB) #$(OFLAG)
CCFLAGS1 = $(WFLAG) $(DEBUG) -I $(INC_DIR) -L $(LIB_DIR) -L $(LIB_DIR1) $(LIB) #$(OFLAG)

all: print

print: print.o stat.o main.o qnorm.o 
	$(CC) $(CCFLAGS1) print.o stat.o qnorm.o main.o -o print

main.o: main.cpp
	$(CC) $(CCFLAGS) main.cpp

print.o: print.cpp
	$(CC) $(CCFLAGS) print.cpp

stat.o: stat.cpp
	$(CC) $(CCFLAGS) stat.cpp

qnorm.o: qnorm.cpp
	$(CC) $(CCFLAGS) qnorm.cpp

clean:
	rm -fr *.o
	rm -f print
