CC = gcc
CCPARAM = -Wall -O3
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
NORMPLATELIB = `pkg-config --cflags --libs gtk+-3.0` -lm
CVPERSPLIB = `pkg-config --cflags --libs gtk+-3.0 opencv` -lm
SRCINC = ../utility
NORMPLATESRC = np_edgedetect.c normplate.c fft.c
NORMPLATEBIN = normplate
CVPERSPSRC = cv_persp.c
CVPERSPBIN = cvpersp

all:
	$(CC) $(NORMPLATESRC) $(UTILSRC) -I $(SRCINC) -o $(NORMPLATEBIN) $(CCPARAM) $(NORMPLATELIB)
	$(CC) $(CVPERSPSRC) $(UTILSRC) -I $(SRCINC) -o $(CVPERSPBIN) $(CCPARAM) $(CVPERSPLIB)
