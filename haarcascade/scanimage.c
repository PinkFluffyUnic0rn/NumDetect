#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <pthread.h>

#include <cairo.h>

#include "nd_image.h"
#include "nd_error.h"
#include "hc_hcascade.h"
#include "hc_rect.h"

#include "hc_scanimgpyr.h"

int imgcroptofile(struct nd_image *img, struct hc_rect *r, const char *path)
{
	struct nd_image imginwin;

	if (nd_imgcreate(&imginwin, abs(r->x1 - r->x0),
		abs(r->y1 - r->y0), img->format) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}
	
	if (nd_imgcrop(img, r->x0, r->y0,
		abs(r->x1 - r->x0),
		abs(r->y1 - r->y0), &imginwin) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}
		
	if (nd_imgwrite(&imginwin, path) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	return 0;
}

int main(int argc, char **argv)
{	
	struct hc_hcascade hc;
	struct nd_image img;
	struct hc_scanconfig scanconf;

	struct hc_rect *newr;
	int newrc;

	int rn;
	
	if (hc_hcascaderead(&hc, argv[1]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
	
	if (nd_imgread(argv[2], &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imggrayscale(&img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	hc_confbuild(&scanconf, hc.ww, hc.wh, img.w, img.h, 0.9, 1, 1, 8);

	struct timespec ts, te;
	clock_gettime(CLOCK_REALTIME, &ts);

	hc_imgpyramidscan(&hc, &img, &newr, &newrc, &scanconf);
	
	clock_gettime(CLOCK_REALTIME, &te);
	printf("%lf\n", te.tv_sec + te.tv_nsec * 1e-9
		- (ts.tv_sec + ts.tv_nsec * 1e-9));

	for (rn = 0; rn < newrc; ++rn) {
		char a[255];

		sprintf(a, "%s/%u.png", argv[3], rn);
		
		double rh = abs(newr[rn].y0 - newr[rn].y1);

		newr[rn].y0 = (int) (newr[rn].y0) - (rh * 0.15);
		newr[rn].y1 = newr[rn].y1 + (rh * 0.15);

		if (newr[rn].y0 >= 0 && newr[rn].y1 < img.h)
			if (imgcroptofile(&img, newr + rn, a) < 0) {
				return 1;
			}
	}

	return 0;
}
