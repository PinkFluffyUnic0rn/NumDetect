#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nd_image.h"
#include "nd_error.h"
#include "np_edgedetect.h"

#include <gtk/gtk.h>
#include <cairo.h>

int main(int argc, char **argv)
{
	struct nd_image img;
	struct nd_image *imgparts;
	double *linestop;
	double *linesbot;
	double *linesleft;
	double *linesright;
	int linescount;
	struct lineseg ltop;
	struct lineseg lbot;
	struct lineseg lleft;
	struct lineseg lright;
	int li1, li2;
	double rho, theta;	
	double interx, intery;

// loading image
	if (argc < 3) {
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

	if (nd_imgwrite(argv[2], &img) < 0) {
		fprintf(stderr, "nd_imgwrite: %s\n",
			nd_strerror(nd_error));
	}

	return 0;
}
