CC = gcc

UTILINC = ../utility
UTILOBJ = $(UTILINC)/nd_error.o $(UTILINC)/nd_image.o $(UTILINC)/nd_vecmat.o
NNETINITOBJ = $(UTILOBJ) nn_neuronet.o nnetinit.o
NNETTEACHOBJ = $(UTILOBJ) nn_neuronet.o nnetteach.o
READSYMOBJ = $(UTILOBJ) nn_neuronet.o readsym.o
NNETINITBIN = ../bin/nnetinit
NNETTEACHBIN = ../bin/nnetteach
READSYMBIN = ../bin/readsym

LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
INC = -I$(UTILINC)
CFLAGS = $(INC) -Wall $(LIB)

all : $(NNETINITBIN) $(NNETTEACHBIN) $(READSYMBIN)

$(NNETINITBIN) : $(NNETINITOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(NNETTEACHBIN) : $(NNETTEACHOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(READSYMBIN) : $(READSYMOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

%.o : %.c $(INC)/*.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : clean

clean :
	rm -f $(NNETINITOBJ) $(NNETTEACHOBJ) $(READSYMOBJ)
