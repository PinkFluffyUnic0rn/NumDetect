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
		printf("%s", "Too few arguments.");
		return 1;
	}

	if (hc_readtrset(&ts, argv[1]) < 0) {
		printf("hc_readtrset: %s\n", nd_strerror(nd_error));
		return 1;
	}

	if (hc_hcascaderead(&hc, argv[2]) < 0) {
		printf("hc_hcascaderead: %s\n", nd_strerror(nd_error));
		return 1;
	}

	if (hc_buildcascade(&hc, &ts, 0.5, 1.0) < 0) {
		printf("hc_buildcascade: %s\n", nd_strerror(nd_error));
		return 1;
	}

	if (hc_hcascadewrite(&hc, argv[3]) < 0) {
		printf("hc_hcascadewrite: %s\n", nd_strerror(nd_error));
		return 1;
	}

	return 0;
}
