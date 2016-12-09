#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cairo.h>

#include "nd_image.h"
#include "neuronet.h"

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

void loadimgtoinput(const char *imgpath, double *input)
{
	struct nd_image img;
	int pixn;
	
	nd_imgread(imgpath, &img); 
	nd_imggrayscale(&img);

	for (pixn = 0; pixn < img.w * img.h; ++pixn) {
		input[pixn] = 2.0 * img.data[pixn] - 1.0;
	}

	free(img.data);
}

void loadexlist(const char *elpath, struct nn_neuronet *nnet,
	struct exampleslist *el)
{
	char buf[TEXTBUFSZ];
	int imgw, imgh;
	int levelc;
	int *neuronc;
	int leveln, en, cn;

	FILE *elfile;
	
	elfile = fopen(elpath, "r");

	fgets(buf, TEXTBUFSZ, elfile);
	buf[strlen(buf)-1] = buf[strlen(buf)-1] == '\n'
		? '\0' : buf[strlen(buf)-1];

	imgw = atoi(strtok(buf, " \t"));
	imgh = atoi(strtok(NULL, " \t"));

	fgets(buf, TEXTBUFSZ, elfile);
	buf[strlen(buf)-1] = buf[strlen(buf)-1] == '\n'
		? '\0' : buf[strlen(buf)-1];

	levelc = atoi(strtok(buf, " \t"));
	
	neuronc = (int *) malloc(sizeof(int) * levelc);
	
	neuronc[0] = atoi(strtok(NULL, " \t"));

	for (leveln = 1; leveln < levelc; ++leveln)
		neuronc[leveln] = atoi(strtok(NULL, " \t"));
	
	nn_createneuronet(nnet, imgw * imgh, levelc, neuronc);

	el->imgw = imgw;
	el->imgh = imgh;

	fgets(buf, TEXTBUFSZ, elfile);

	el->examplec = atoi(strtok(buf, " \t"));
	el->classc = atoi(strtok(NULL, " \t"));
	el->iterc = atoi(strtok(NULL, " \t"));
	
	el->imgpath = (char **) malloc(sizeof(char *) * el->examplec);
	el->imgclass = (int *) malloc(sizeof(int) * el->examplec);

	en = 0;	
	cn = -1;
	int badn = 0;
	while (en < el->examplec) {
		size_t len;
		ssize_t readc;
		char *tmps;

		tmps = NULL;
		len = 0;
		readc = getline(&tmps, &len, elfile);			
		
		tmps[readc - 1] = tmps[readc - 1] == '\n'
			? '\0' : tmps[readc - 1];

		if (tmps[0] == '-') {
			if (strcmp("bad", tmps + 1) == 0)
				badn = 1;
			else {
				++cn;
				badn = 0;
			}
		}
		else {
			el->imgpath[en] = tmps;
			el->imgclass[en] = (badn ? -1 : cn);	
			++en;
		}
	}
	
	el->classc = cn + 1;
	
	fclose(elfile);
}

static void shufflearr(int *pathidx, int len)
{
	int p1;

	for (p1 = 0; p1 < len; ++p1) {
		int p2;
		int tmp;
	
		p2 = rand() % len;
		tmp = pathidx[p1];
		
		pathidx[p1] = pathidx[p2];
		pathidx[p2] = tmp;	
	}
}

void teachneuronet(struct nn_neuronet *nnet, struct exampleslist *el)
{
	int *pathidx;
	int i, en;
	int cn;

	double *input = (double *) malloc(sizeof(double) * el->imgw * el->imgh);
	double *target = (double *) malloc(sizeof(double) * el->classc);

	pathidx = malloc(sizeof(int) * el->examplec);
	for (i = 0; i < el->examplec; ++i)
		pathidx[i] = i;

	for (i = 0; i < el->iterc; ++i) {
		shufflearr(pathidx, el->examplec);

		for (en = 0; en < el->examplec; ++en) {	
			int idx;
			
			idx = pathidx[en];
				
			loadimgtoinput(el->imgpath[idx], input);

			for (cn = 0; cn < el->classc; ++cn)
				target[cn] = (cn == el->imgclass[idx])
					? 1.0 : 0.0;
			nn_teachneuronet(nnet, input, target, 1.0);
		}

		printf("%1.3f%% done.\n", 100.0 * (i + 1) / el->iterc);
	}

	free(pathidx);
}

int main(int argc, const char **argv)
{	
	struct nn_neuronet nnet;
	struct exampleslist el;

	loadexlist(argv[1], &nnet, &el);
	
	teachneuronet(&nnet, &el);

	nn_neuronettofile(&nnet, argv[2]);

	return 0;
}
