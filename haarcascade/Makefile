CC = gcc
UTILINC = ../utility
UTILOBJ = $(UTILINC)/nd_error.o $(UTILINC)/nd_image.o $(UTILINC)/nd_vecmat.o
INITHCOBJ = $(UTILOBJ) initcascade.o hc_hcascade.o
FINDWCOBJ = $(UTILOBJ) findwcs.o hc_hcascade.o
BUILDHCOBJ = $(UTILOBJ) buildcascade.o hc_hcascade.o
SCANIMGOBJ = $(UTILOBJ) scanimage.o hc_rect.o hc_hcascade.o hc_scanimgpyr.o
INITHCBIN = ../bin/initcascade
FINDWCBIN = ../bin/findwcs
BUILDHCBIN = ../bin/buildcascade
SCANIMGBIN = ../bin/scanimage

LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
INC = -I$(UTILINC)
CFLAGS = $(INC) -Wall -O3 $(LIB)

all : $(INITHCBIN) $(FINDWCBIN) $(BUILDHCBIN) $(SCANIMGBIN)

$(INITHCBIN) : $(INITHCOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(FINDWCBIN) : $(FINDWCOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(BUILDHCBIN) : $(BUILDHCOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(SCANIMGBIN) : $(SCANIMGOBJ)
	$(CC) -o $@ $^ $(CFLAGS)

%.o : %.c $(INC)/*.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : clean

clean :
	rm -f $(INITHCOBJ) $(FINDWCOBJ) $(BUILDHCOBJ) $(SCANINGOBJ)
