#ifndef ND_IMAGE_H
#define ND_IMAGE_H

#include "nd_vecmat.h"

struct nd_image {
	int w;
	int h;
	int chans;
	double *data;
};

int nd_imgisvalid(const struct nd_image *img);

int nd_imgcreate(struct nd_image *img, int w, int h, int chans);

int nd_imgdestroy(struct nd_image *img);

int nd_imgread(const char *imgpath, struct nd_image *img);

int nd_imgwrite(const char *imgpath, const struct nd_image *img);

int nd_imghsvval(struct nd_image *img);

int nd_imggrayscale(struct nd_image *img);

int nd_imgnormalize(struct nd_image *img, int avr, int dev);

int nd_imgcrop(const struct nd_image *imgin, int x0, int y0, int w, int h,
	struct nd_image *imgout);

int nd_imgscale(const struct nd_image *inimg, double wrel, double hrel,
	struct nd_image *outimg);

int nd_imgscalebicubic(const struct nd_image *inimg, double wrel, double hrel,
	struct nd_image *outimg);

int nd_imgscalebilinear(const struct nd_image *inimg, double wrel, double hrel,
	struct nd_image *outimg);

int nd_getpersptransform(double *inpoints, double *outpoints,
	struct nd_matrix3 *mr);

int nd_imgapplytransform(struct nd_image *img, const struct nd_matrix3 *m);

#endif
