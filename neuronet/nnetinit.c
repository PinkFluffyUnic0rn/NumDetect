#include <stdio.h>
#include <stdlib.h>

#include "nd_error.h"
#include "nn_neuronet.h"

int main(int argc, const char **argv)
{	
	struct nn_neuronet nnet;
	int imgw, imgh;
	int *neuronc;
	int levelc;
	char *next;
	int i;

	imgw = strtol(*(++argv), &next, 0);
	if (*argv == next) { 
		fprintf(stderr, "Wrong format of width.\n");
		return 1;
	}

	imgh = strtol(*(++argv), &next, 0);
	if (*argv == next) { 
		fprintf(stderr, "Wrong format of height.\n");
		return 1;
	}

	levelc = strtol(*(++argv), &next, 0);
	if (*argv == next) { 
		fprintf(stderr, "Wrong format of level count.\n");
		return 1;
	}

	if ((neuronc = malloc(sizeof(int) * levelc)) == NULL) {
		fprintf(stderr, "cannot allocate memory");
		return 1;
	}

	for (i = 0; i < levelc; ++i) {
		neuronc[i] = strtol(*(++argv), &next, 0);
		if (*argv == next) { 
			fprintf(stderr, "Wrong format of neuron count.\n");
			return 1;
		}
	}

	if (nn_createneuronet(&nnet, imgw * imgh, levelc, neuronc) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nn_neuronettofile(&nnet, *(++argv))) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	return 0;
}
