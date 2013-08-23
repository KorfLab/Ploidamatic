#######################
# Makefile for Ploidy #
#######################

LIB = -lm -lstochhmm
INC = -IStochHMM/src
LNC = -LStochHMM/src
CC  = g++
CFLAGS = -O2 -Wall


APP = ploidy
SRC = ploidy.cpp
OBJ = ploidy.o

DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########


$(APP): $(OBJ)
	$(CC) -o $(APP) $(CFLAGS) $(LNC) $(OBJ) $(LIB)

clean:
	rm -f *.o $(APP1)


###################
# Inference Rules #
###################

%.o: %.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<


