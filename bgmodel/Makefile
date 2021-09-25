CC = gcc

UTILINC = ./

UTILOBJ = $(UTILINC)/nd_error.o $(UTILINC)/nd_image.o $(UTILINC)/nd_vecmat.o\
$(UTILINC)/nd_procsync.o
READVIDEOOBJ= $(UTILOBJ) $(HCOBJ) $(EDOBJ) readvideo.o gui.o em3d.o bgmodel.o

READVIDEOBIN = ../bin/bgmodel

INC = -I$(UTILINC) -I$(HCINC) -I$(EDINC)
LIB=-lavcodec -lavdevice -lavutil -lavformat -lswscale -pthread -lm  \
-lxcb -lxcb-keysyms `pkg-config --cflags --libs gtk+-3.0`
CFLAGS = $(INC) -Wall $(LIB)

all : $(READVIDEOBIN)

$(READVIDEOBIN): $(READVIDEOOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

%.o : %.c $(INC)/*.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : clean

clean :
	rm -f $(READVIDEOOBJ) 
