#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cairo.h>

#include "nd_error.h"
#include "hc_hcascade.h"

int main(int argc, const char **argv)
{
	struct hc_hcascade hc;
	struct hc_trainingset ts;

	if (argc < 4) {
		fprintf(stderr, "%s", "Too few arguments.");
		return 1;
	}

	if (hc_readtrset(&ts, argv[1]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (hc_hcascaderead(&hc, argv[2]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
	
	if (hc_buildcascade(&hc, &ts, 0.5, 0.97) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (hc_hcascadewrite(&hc, argv[3]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	return 0;
}
