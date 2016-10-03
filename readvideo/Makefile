CC=gcc
SRC=readvideo.c hc_hcascade.c hc_rect.c hc_scanimgpyr.c
OPT=-Wall
LIB=-lavcodec -lavutil -lavformat -lswscale -pthread -lm \
`pkg-config --cflags --libs gtk+-3.0`
BIN=readvideo
SRCINC = ../utility
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c

all:
	$(CC) $(SRC) $(UTILSRC) -I $(SRCINC) -o $(BIN) $(OPT) $(LIB) 
