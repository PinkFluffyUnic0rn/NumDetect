CC = gcc
CCPARAM = -Wall -O3
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
NORMPLATELIB = `pkg-config --cflags --libs gtk+-3.0` -lm
UTILINC = ../utility
NORMPLATESRC = ed_fft.c ed_edgedetect.c ed_findborder.c normplate.c
NORMPLATEBIN = normplate

all:
	$(CC) $(NORMPLATESRC) $(UTILSRC) -I $(UTILINC) -o $(NORMPLATEBIN) $(CCPARAM) $(NORMPLATELIB)
