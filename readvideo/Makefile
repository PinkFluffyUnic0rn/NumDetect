CC=gcc
SRC=readvideo.c
HCSRC=../haarcascade/hc_hcascade.c ../haarcascade/hc_rect.c\
../haarcascade/hc_scanimgpyr.c
UTILSRC=../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c\
../utility/nd_procsync.c
OPT=-Wall
LIB=-lavcodec -lavutil -lavformat -lswscale -pthread -lm \
`pkg-config --cflags --libs gtk+-3.0`
BIN=readvideo
UTILINC=../utility
HCINC=../haarcascade

all:
	$(CC) $(SRC) $(HCSRC) $(UTILSRC) -I$(UTILINC) -I$(HCINC) \
	-o $(BIN) $(OPT) $(LIB) 
