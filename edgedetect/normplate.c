#include <stdio.h>

#include "ed_findborder.h"
#include "nd_image.h"
#include "nd_error.h"

int main(int argc, char **argv)
{
	struct nd_image img;
	struct nd_matrix3 m;
	double inpoints[8];
	double outpoints[8];

// loading image
	if (argc < 2) {
		fprintf(stderr, "Not enough arguments.\n");
		return 1;
	}
	
	if (nd_imgread(argv[1], &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imghsvval(&img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
/*
	ed_removelowfreq(&img, 0.05, 0.0);
*/

	if (ed_findborder(&img, inpoints) < 0) {
		fprintf(stderr, "Border not found.\n");
		return 1;
	}
	
	outpoints[0] = 0.0; outpoints[1] = 0.0;
	outpoints[2] = img.w; outpoints[3] = 0.0;
	outpoints[4] = img.w; outpoints[5] = img.h;
	outpoints[6] = 0.0; outpoints[7] = img.h;

	printf("%f %f\n%f %f\n%f %f\n%f %f\n\n", 
		inpoints[0], inpoints[1],
		inpoints[2], inpoints[3],
		inpoints[4], inpoints[5],
		inpoints[6], inpoints[7]);

	if (nd_getpersptransform(inpoints, outpoints, &m) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imgapplytransform(&img, &m, &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	nd_imgwrite(&img, argv[2]);

	return 0;
}
