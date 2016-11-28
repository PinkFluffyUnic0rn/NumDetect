CC = gcc
CCPARAM = -Wall
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
UTILINC = ../utility
NNETTEACHSRC = neuronet.c nnetteach.c
NNETTEACHBIN = nneteach
NNETTESTSRC = neuronet.c nnettest.c
NNETTESTBIN = readsym
all:
	$(CC) $(NNETTEACHSRC) $(UTILSRC) -I$(UTILINC) -o $(NNETTEACHBIN) $(CCPARAM) $(LIB)
	$(CC) $(NNETTESTSRC) $(UTILSRC) -I$(UTILINC) -o $(NNETTESTBIN) $(CCPARAM) $(LIB)
