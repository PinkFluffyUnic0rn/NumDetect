CC = gcc
CCPARAM = -Wall -O3
UTILSRC = utility/nd_error.c utility/nd_image.c
SRCINC = utility
SCANIMGLIB = `pkg-config --cflags --libs gtk+-3.0` -lm
SCANIMGSRC = haarcascade/scanimage.c haarcascade/hc_rect.c \
	haarcascade/hc_hcascade.c haarcascade/hc_scanimgpyr.c
SCANIMGBIN = bin/scanimage
NORMPLATELIB = `pkg-config --cflags --libs gtk+-3.0` -lm
NORMPLATESRC = normplate/np_edgedetect.c normplate/normplate.c normplate/fft.c
NORMPLATEBIN = bin/normplate
CVPERSPLIB = `pkg-config --cflags --libs opencv gtk+-3.0` -lm
CVPERSPSRC = normplate/cv_persp.c
CVPERSPBIN = bin/cvpersp

all:
	$(CC) $(SCANIMGSRC) $(UTILSRC) -I $(SRCINC) -o $(SCANIMGBIN) $(CCPARAM) $(SCANIMGLIB)
	$(CC) $(NORMPLATESRC) $(UTILSRC) -I $(SRCINC) -o $(NORMPLATEBIN) $(CCPARAM) $(NORMPLATELIB)
	$(CC) $(CVPERSPSRC) $(UTILSRC) -I $(SRCINC) -o $(CVPERSPBIN) $(CCPARAM) $(CVPERSPLIB)
