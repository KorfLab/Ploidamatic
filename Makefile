#######################
# Makefile for Ploidy #
#######################

APP1 = ploidy
SRC1 = ploidy.cpp

APP2 = cpw
SRC2 = cpw.c

###########
# Targets #
###########

default:
	make $(APP1)
	make $(APP2)


$(APP1): $(SRC1)
	g++ -o $(APP1) $(SRC1) -IStochHMM/src -LStochHMM/src -lm -lstochhmm

$(APP2): $(SRC2)
	gcc -o $(APP2) $(SRC2) -O2 -Wall -Werror

clean:
	rm -f *.o $(APP1) $(APP2)


