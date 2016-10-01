#ifndef HCASCADE_H
#define HCASCADE_H

#include <stdlib.h>

#include "nd_image.h"

struct hc_wclassifier
{
	int fn;
	int h;
	int w;
	int x;
	int y;

	double ineqdir;
	double thres;
};

struct hc_stage {
	int wcc;
	double thres;
};

struct hc_hcascade {
	int wh;
	int ww;

	struct nd_image *feature;
	int featurec;

	struct hc_wclassifier *wc;
	double *wccoef;
	int wccount;

	struct hc_stage *stage;
	int stagecount;
};

struct hc_trainingset {
	struct nd_image *img;
	int *isgood;
	int imgc;
};

int hc_readtrset(struct hc_trainingset *ts, const char *hcpath);

int hc_create(struct hc_hcascade *hc, int wh, int ww,
	const struct nd_image *feature, int featurec);

int hc_hcascaderead(struct hc_hcascade *hc, const char *hcpath);

int hc_hcascadewrite(struct hc_hcascade *hc, const char *hcpath);

int hc_imgintegral(struct nd_image *img);

int hc_findwc(struct hc_hcascade *hc, struct hc_trainingset *ts,
	double **weights, int iterc);

int hc_buildcascade(struct hc_hcascade *hc, struct hc_trainingset *ts,
	double fprtarget, double fntarget);

int hc_imgclassify(const struct hc_hcascade *hc, const struct nd_image *img);

int hc_fastimgclassify(const struct hc_hcascade *hc,
	const struct nd_image *img, int imgstride, double sd);

#endif
