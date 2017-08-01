CC = gcc

UTILINC = ../utility
UTILOBJ = $(UTILINC)/nd_error.o $(UTILINC)/nd_image.o $(UTILINC)/nd_vecmat.o
NORMPLATEOBJ = $(UTILOBJ) ed_edgedetect.o ed_findborder.o normplate.o
NORMPLATEBIN = ../bin/normplate

LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
INC = -I$(UTILINC)
CFLAGS = $(INC) -Wall $(LIB)

all : $(NORMPLATEBIN)

$(NORMPLATEBIN) : $(NORMPLATEOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

%.o : %.c $(INC)/*.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : clean

clean :
	rm -f $(NORMPLATEOBJ)
