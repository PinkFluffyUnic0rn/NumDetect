#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cairo.h>

#include "nd_error.h"
#include "hc_hcascade.h"

int loadweights(double **weight, int *wcount)
{
	FILE *weightfile;
	char *cur;
	char *next;
	char *str;
	size_t strsz;
	int wn;

	if ((weightfile = fopen("weight", "r")) == NULL) {
		*weight = NULL;
		return 1;
	}

	strsz = 0;	
	if (getline(&str, &strsz, weightfile) <= 0)
		return (-1);

	*wcount = strtol(str, &next, 0);

	if (str == next)
		return (-1);

	free(str);
	str = NULL;
	
	if ((*weight = malloc(sizeof(double) * *wcount)) == NULL)
		return (-1);

	if (getline(&str, &strsz, weightfile) <= 0)
		return (-1);

	cur = str;

	for (wn = 0; wn < *wcount; ++wn) {
		(*weight)[wn] = strtod(cur, &next);

		if (cur == next)
			return (-1);

		cur = next;
	}
	
	fclose(weightfile);

	return 0;
}

int saveweights(double *weight, int wcount)
{
	FILE *weightfile;
	int wn;

	if ((weightfile = fopen("weight", "w")) == NULL)
		return 1;

	fprintf(weightfile, "%d\n", wcount);

	for (wn = 0; wn < wcount; ++wn)
		fprintf(weightfile, "%lf%c", weight[wn],
			wn != (wcount - 1) ? ' ' : '\n');
	
	fclose(weightfile);

	return 0;
}

int main(int argc, const char **argv)
{
	struct hc_hcascade hc;
	struct hc_trainingset ts;
	double *weight;
	int iters;
	char *next;
	int wcount;

	if (argc < 5) {
		printf("%s", "Too few agruments.");
		return 1;
	}

	if (hc_readtrset(&ts, argv[1]) < 0) {
		printf("hc_readtrset: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	if (hc_hcascaderead(&hc, argv[2]) < 0) {
		printf("hc_hcascaderead: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	iters = strtol(argv[4], &next, 0);

	if (argv[4] == next) {
		printf("%s", "Wrong format of iterations count.\n");
		return 1;
	}

	if (loadweights(&weight, &wcount) < 0
		|| (weight != NULL && wcount != ts.imgc)) {
		printf("%s", "Error while loading weights.\n");
		return 1;	
	}

	if (hc_findwc(&hc, &ts, &weight, iters) < 0) {
		printf("hc_findwcs: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	saveweights(weight, ts.imgc);

	if (hc_hcascadewrite(&hc, argv[3]) < 0) {
		printf("hc_hcascadewrite: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	return 0;
}
