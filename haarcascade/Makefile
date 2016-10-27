CC = gcc
CCPARAM = -Wall -O3
UTILSRC = ../utility/nd_error.c ../utility/nd_image.c ../utility/nd_vecmat.c
LIB = `pkg-config --cflags --libs gtk+-3.0` -lm
UTILINC = ../utility
INITHCSRC = initcascade.c hc_hcascade.c
INITHCBIN = initcascade
FINDWCSRC = findwcs.c hc_hcascade.c
FINDWCBIN = findwcs
BUILDHCSRC = buildcascade.c hc_hcascade.c
BUILDHCBIN = buildcascade
SCANIMGSRC = scanimage.c hc_rect.c hc_hcascade.c hc_scanimgpyr.c
SCANIMGBIN = scanimage
all:
	$(CC) $(INITHCSRC) $(UTILSRC) -I$(UTILINC) -o $(INITHCBIN) $(CCPARAM) $(LIB)
	$(CC) $(FINDWCSRC) $(UTILSRC) -I$(UTILINC) -o $(FINDWCBIN) $(CCPARAM) $(LIB)
	$(CC) $(BUILDHCSRC) $(UTILSRC) -I$(UTILINC) -o $(BUILDHCBIN) $(CCPARAM) $(LIB)
	$(CC) $(SCANIMGSRC) $(UTILSRC) -I$(UTILINC) -o $(SCANIMGBIN) $(CCPARAM) $(LIB)
