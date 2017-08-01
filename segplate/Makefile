CC = gcc

UTILINC = ../utility
HCINC =../haarcascade
NNINC =../neuronet

UTILOBJ = $(UTILINC)/nd_error.o $(UTILINC)/nd_image.o $(UTILINC)/nd_vecmat.o
HCOBJ= $(HCINC)/hc_hcascade.o $(HCINC)/hc_rect.o
NNOBJ= $(NNINC)/nn_neuronet.o
NUMSEGOBJ= $(UTILOBJ) $(HCOBJ) $(NNOBJ) numseg.o

NUMSEGBIN = ../bin/numseg

LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
INC = -I$(UTILINC) -I$(HCINC) -I$(NNINC)
CFLAGS = $(INC) -Wall $(LIB)

all : $(NUMSEGBIN)

$(NUMSEGBIN) : $(NUMSEGOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

%.o : %.c $(INC)/*.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : clean

clean :
	rm -f $(NUMSEGOBJ) 
