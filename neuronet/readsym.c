#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nd_image.h"
#include "nd_error.h"
#include "nn_neuronet.h"

#define WIN_W 9
#define WIN_H 20

int loadimgtoinput(const char *imgpath, double *input)
{
	struct nd_image img;
	int pixn;
	
	if (nd_imgread(imgpath, &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}
	
	if (nd_imggrayscale(&img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	for (pixn = 0; pixn < img.w * img.h; ++pixn) {
		input[pixn] = 2.0 * img.data[pixn] - 1.0;
	}

	free(img.data);

	return 0;
}

int findpattern(struct nn_neuronet *nnet, double *input)
{
	double *out;
	uint minoutid;
	int i;
	uint outputc;
	
	outputc = nnet->nodec[nnet->levelc - 1];
	
	if ((out = (double *) malloc(sizeof(double) * outputc)) == NULL)
		return (-1);
	
	nn_neuroneteval(nnet, input, out);

	minoutid = 0;
	for (i = 1; i < outputc; ++i)
		if (out[i] > out[minoutid])
			minoutid = i;

	// should do something with threshold	
	if (out[minoutid] > 0.5) {
		free(out);
		return minoutid;
	} else {
		free(out);
		return -1;
	}
}

int main(int argc, const char **argv)
{
	double *imgdata;
	struct nn_neuronet nnet;

	if (argc < 3) {
		printf("Not enough arguments.\n");
		exit(1);
	}

	imgdata = (double *) malloc(sizeof(double *) * WIN_W * WIN_H);
	
	if (loadimgtoinput(argv[2], imgdata) < 0)
		return 1;

	if (nn_neuronetfromfile(&nnet, argv[1]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	int res;
	char a[23] = "0123456789abcehkmoptxy";
	
	res = findpattern(&nnet, imgdata);

	if (res != -1)
		printf("%c", a[res] );
	
	fflush(stdout);

	return 0;
}
