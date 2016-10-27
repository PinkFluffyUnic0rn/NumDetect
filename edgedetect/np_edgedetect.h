#ifndef EDGEDETECT_H
#define EDGEDETECT_H

#include <stdlib.h>
#include <stdio.h>

#include "nd_image.h"

struct lineseg {
	double x0;
	double y0;
	double x1;
	double y1;
	int pointc;
};

int gaussblur(struct nd_image *img, double sigma);

int np_removelowfreq(struct nd_image *img, double wthres, double hthres);

int np_canny(struct nd_image *img, int *res, double thres1, double thres2);

int np_hough(int *img, int imgw, int imgh, double dang,
	double *lines, int linemaxc);

int np_houghseg(int *img, int imgw, int imgh, double *lines, int linec,
	struct lineseg *lineseg, int linesegc, double maxgap);

#endif
