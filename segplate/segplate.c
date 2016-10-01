#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nd_image.h"
#include "nd_error.h"
#include "np_edgedetect.h"

int main(int argc, char **argv)
{
	struct nd_image img;
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

	if (np_canny(&img, mask, -1.0, -1.0) < 0)
		return (-1);


///////////////////////////////////////////////////////////////////////////////
	struct nd_image test;
	int x, y;

	test.data = malloc(sizeof(double) * img.w * img.h);
	test.w = img.w;
	test.h = img.h;
	test.chans = 1;

	for (y = 0; y < test.h; ++y)
		for (x = 0; x < test.w; ++x)
			test.data[y * test.w + x]
				= mask[y * test.w + x] ? 1.0 : 0.0;
	
	nd_imgwrite(argv[2], &test);

///////////////////////////////////////////////////////////////////////////////


	return 0;
}
