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

int ed_gaussblur(struct nd_image *img, double sigma);

int ed_removelowfreq(struct nd_image *img, double wthres, double hthres);

int ed_removehorfreq(struct nd_image *img, struct nd_image *res);

int ed_canny(struct nd_image *img, int *res, double thres1, double thres2,
	int onlyver);

int ed_hough(int *img, int imgw, int imgh, double dang,
	double *lines, int linemaxc);

int ed_houghseg(int *img, int imgw, int imgh, double *lines, int linec,
	struct lineseg *lineseg, int linesegc, double maxgap);

#endif
