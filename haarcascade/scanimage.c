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

int main(int argc, char **argv)
{	
	struct hc_hcascade hc;
	struct nd_image img;
	struct hc_scanconfig scanconf;

	struct hc_rect *newr;
	int newrc;

	int rn;
	
	time_t t;

	if (hc_hcascaderead(&hc, argv[1]) < 0) {
		fprintf(stderr, "hc_hcascaderead: %s.\n",
			argv[1]);
		return 1;
	}
	
	if (nd_imgread(argv[2], &img) < 0) {
		fprintf(stderr, "nd_imgread: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	if (nd_imggrayscale(&img) < 0) {
		fprintf(stderr, "nd_imggrayscale: %s.\n",
			nd_strerror(nd_error));
		return 1;
	}

	scanconf.scalestep = 0.9;
	scanconf.winhstep = 1;
	scanconf.winwstep = 1;

	t = time(0);

	imgpyramidscan(&hc, &img, &newr, &newrc, &scanconf);

	printf("%lu\n", time(0) - t);
	
	for (rn = 0; rn < newrc; ++rn) {
		char a[255];

		sprintf(a, "%s/%u.png", argv[3], rn);
		
		double rh = abs(newr[rn].y0 - newr[rn].y1);

		newr[rn].y0 = (int) (newr[rn].y0) - (rh * 0.15);
		newr[rn].y1 = newr[rn].y1 + (rh * 0.15);

		if (newr[rn].y0 >= 0 && newr[rn].y1 < img.h) {
			struct nd_image imginwin;

			nd_imgcreate(&imginwin, abs(newr[rn].x1 - newr[rn].x0),
				abs(newr[rn].y1 - newr[rn].y0), img.chans);
			
			nd_imgcrop(&img,
				newr[rn].x0, newr[rn].y0,
				abs(newr[rn].x1 - newr[rn].x0),
				abs(newr[rn].y1 - newr[rn].y0), &imginwin);
			
			nd_imgwrite(a, &imginwin);
		}
	}

	return 0;
}
