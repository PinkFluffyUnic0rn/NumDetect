CC = gcc
CCPARAM = -Wall -O3
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
INC = ../utility
SRC = segplate.c np_edgedetect.c fft.c
BIN = segplate

all:
	$(CC) $(SRC) $(UTILSRC) -I $(INC) -o $(BIN) $(CCPARAM) $(LIB)
