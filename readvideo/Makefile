CC = gcc

UTILINC = ../utility
BGMINC=../bgmodel
HCINC=../haarcascade
EDINC=../edgedetect

UTILOBJ = $(UTILINC)/nd_error.o $(UTILINC)/nd_image.o $(UTILINC)/nd_vecmat.o\
$(UTILINC)/nd_procsync.o
BGMOBJ=$(BGMINC)/bgmodel.o $(BGMINC)/em3d.o
HCOBJ= $(HCINC)/hc_hcascade.o $(HCINC)/hc_rect.o $(HCINC)/hc_scanimgpyr.o
EDOBJ = $(EDINC)/ed_edgedetect.o $(EDINC)/ed_findborder.o $(EDINC)/ed_fft.o
READVIDEOOBJ= $(UTILOBJ) $(BGMOBJ) $(HCOBJ) $(EDOBJ) readvideo.o gui.o

READVIDEOBIN = ../bin/readvideo

INC = -I$(UTILINC) -I$(BGMINC) -I$(HCINC) -I$(EDINC)
LIB=-lavcodec -lavutil -lavformat -lswscale -pthread -lm  \
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
