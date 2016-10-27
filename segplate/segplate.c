#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nd_image.h"
#include "nd_error.h"
#include "ed_edgedetect.h"

#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

int main(int argc, char **argv)
{
	struct nd_image img;
	struct nd_image canny;
	struct nd_image hist;
	int x, y;
	int *mask;
	
	if (argc < 2) {
		fprintf(stderr, "Not enough arguments.\n");
		return 1;
	}
	
	if (nd_imgread(argv[1], &img) < 0) {
		fprintf(stderr, "nd_imgread: %s\n",
			nd_strerror(nd_error));
		return 1;
	}

	if (nd_imghsvval(&img) < 0) {
		fprintf(stderr, "nd_imghsvval: %s\n",
			nd_strerror(nd_error));
		return 1;
	}
	
	if (nd_imgnormalize(&img, 1, 1) < 0)
		return (-1);

	if ((mask = malloc(sizeof(int) * img.w * img.h)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if (ed_canny(&img, mask, -1.0, -1.0) < 0)
		return (-1);

	canny.data = malloc(sizeof(double) * img.w * img.h);
	canny.w = img.w;
	canny.h = img.h;
	canny.format = ND_PF_GRAYSCALE;

	for (y = 0; y < canny.h; ++y)
		for (x = 0; x < canny.w; ++x)
			canny.data[y * canny.w + x]
				= mask[y * canny.w + x] ? 1.0 : 0.0;
	
	nd_imgwrite(&canny, argv[2]);

	hist.w = img.w;
	hist.h = 100;
	hist.format = ND_PF_GRAYSCALE;
	hist.data = calloc(hist.w * hist.h, sizeof(double));

	for (x = 0; x < img.w; ++x) {
		int s;

		s = 0;
		
		for (y = 0; y < img.h; ++y)
			s += 2 * img.data[y * img.w + x];
	
		for (y = 0; y < MIN(s, img.h); ++y)
			hist.data[y * hist.w + x] = 1;
	}
	
	nd_imgwrite(&hist, "hist.png");

	return 0;
}
