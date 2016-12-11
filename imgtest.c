#include <stdio.h>

#include "utility/nd_image.h"
#include "utility/nd_error.h"

int main(int argc, char **argv)
{	
	struct nd_image img;

	if (nd_imgread(argv[1], &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imgscalebicubic(&img, 0.5, 0.5, &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imgwrite(&img, argv[2]) < 0)
		return (-1);


	return 0;
}
