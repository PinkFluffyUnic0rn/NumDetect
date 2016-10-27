CC=gcc
SRC=segplate.c
EDSRC=../edgedetect/np_edgedetect.c ../edgedetect/fft.c
UTILSRC=../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
OPT=-Wall -O3
LIB=`pkg-config --cflags --libs gtk+-3.0` -lm
BIN=segplate
UTILINC=../utility
EDINC=../edgedetect

all:
	$(CC) $(SRC) $(UTILSRC) $(EDSRC) -I$(UTILINC) -I$(EDINC) \
	-o $(BIN) $(CCPARAM) $(LIB)
