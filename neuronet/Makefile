CC = gcc
CCPARAM = -Wall
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
UTILINC = ../utility
NNETINITSRC = nn_neuronet.c nnetinit.c
NNETINITBIN = nnetinit
NNETTEACHSRC = nn_neuronet.c nnetteach.c
NNETTEACHBIN = nnetteach
NNETTESTSRC = nn_neuronet.c readsym.c
NNETTESTBIN = readsym
all:
	$(CC) $(NNETINITSRC) $(UTILSRC) -I$(UTILINC) -o $(NNETINITBIN) $(CCPARAM) $(LIB)
	$(CC) $(NNETTEACHSRC) $(UTILSRC) -I$(UTILINC) -o $(NNETTEACHBIN) $(CCPARAM) $(LIB)
	$(CC) $(NNETTESTSRC) $(UTILSRC) -I$(UTILINC) -o $(NNETTESTBIN) $(CCPARAM) $(LIB)
