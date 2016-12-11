CC=gcc
SRC=numseg.c
NNSRC=../neuronet/nn_neuronet.c
HCSRC=../haarcascade/hc_hcascade.c ../haarcascade/hc_rect.c
UTILSRC=../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c\
../utility/nd_procsync.c
OPT=-Wall
LIB=-lavcodec -lavutil -lavformat -lswscale -pthread -lm \
`pkg-config --cflags --libs gtk+-3.0`
BIN=numseg
NNINC=../neuronet
HCINC=../haarcascade
UTILINC=../utility

all:
	$(CC) $(SRC) $(NNSRC) $(HCSRC) $(UTILSRC) \
	 -I$(NNINC) -I$(HCINC) -I$(UTILINC) -o $(BIN) $(OPT) $(LIB) 
