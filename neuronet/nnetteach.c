#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cairo.h>

#include "nd_image.h"
#include "nd_error.h"
#include "nn_neuronet.h"

#define TEXTBUFSZ 1024

struct rgb {
	unsigned char b, g, r, a;
};

struct exampleslist {
	int imgw;
	int imgh;
	int examplec;
	int iterc;
	int classc;
	char **imgpath;
	int *imgclass;
};

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

static void shufflearr(int *idx, int len)
{
	int p1;

	for (p1 = 0; p1 < len; ++p1) {
		int p2;
		int tmp;
	
		p2 = rand() % len;
		tmp = idx[p1];
		
		idx[p1] = idx[p2];
		idx[p2] = tmp;	
	}
}

int teachneuronet(struct nn_neuronet *nnet, struct nn_trainingset *ts)
{
	int *imgidx;
	double *target;
	int i;
	
	if ((target = (double *) malloc(sizeof(double) * ts->classcount))
		== NULL)
		return (-1);

	if ((imgidx = malloc(sizeof(int) * ts->imgc)) == NULL)
		return (-1);

	for (i = 0; i < ts->imgc; ++i)
		imgidx[i] = i;

	for (i = 0; i < ts->iterc; ++i) {
		int imgn;

		shufflearr(imgidx, ts->imgc);

		for (imgn = 0; imgn < ts->imgc; ++imgn) {	
			int cn;
			int idx;
				
			idx = imgidx[imgn];
				
			for (cn = 0; cn < ts->classcount; ++cn)
				target[cn] = (cn == ts->imgclass[idx])
					? 1.0 : 0.0;

			if (nn_teachneuronet(nnet, ts->img[idx].data,
				target, 1.0) < 0) {
				fprintf(stderr, nd_geterrormessage());
				return (-1);
			}
		}

		printf("%1.3f%% done.\n", 100.0 * (i + 1) / ts->iterc);
	}

	free(imgidx);

	return 0;
}

int main(int argc, const char **argv)
{	
	struct nn_neuronet nnet;
	struct nn_trainingset ts;
	
	if (argc < 4) {
		fprintf(stderr, "%s", "Too few agruments.");
		return 1;
	}

	if (nn_readtrset(&ts, argv[1]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
	
	if (nn_neuronetfromfile(&nnet, argv[2]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (teachneuronet(&nnet, &ts) < 0)
		return 1;

	if (nn_neuronettofile(&nnet, argv[3])) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	return 0;
}
