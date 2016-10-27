CC = gcc
CCPARAM = -Wall -O3
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
NORMPLATELIB = `pkg-config --cflags --libs gtk+-3.0` -lm
UTILINC = ../utility
NORMPLATESRC = np_edgedetect.c normplate.c fft.c
NORMPLATEBIN = normplate

all:
	$(CC) $(NORMPLATESRC) $(UTILSRC) -I $(UTILINC) -o $(NORMPLATEBIN) $(CCPARAM) $(NORMPLATELIB)
