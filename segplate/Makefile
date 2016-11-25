CC=gcc
SRC=numseg.c
HCSRC=../haarcascade/hc_hcascade.c ../haarcascade/hc_rect.c\
../haarcascade/hc_scanimgpyr.c
EDSRC=../edgedetect/ed_fft.c ../edgedetect/ed_edgedetect.c \
../edgedetect/ed_findborder.c
UTILSRC=../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c\
../utility/nd_procsync.c
OPT=-Wall
LIB=-lavcodec -lavutil -lavformat -lswscale -pthread -lm \
`pkg-config --cflags --libs gtk+-3.0`
BIN=numseg
UTILINC=../utility
HCINC=../haarcascade
EDINC=../edgedetect

all:
	$(CC) $(SRC) $(HCSRC) $(EDSRC) $(UTILSRC) \
	 -I$(HCINC) -I$(EDINC) -I$(UTILINC) -o $(BIN) $(OPT) $(LIB) 
