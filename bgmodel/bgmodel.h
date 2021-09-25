#ifndef BGMODEL_H
#define BGMODEL_H

#include "em3d.h"
#include "nd_image.h"

struct bgmodel {
	struct em_gaussmix *gm;
	int compcount;
	double thres;
	double lrate;
	int *isforeground;
	int w;
	int h;
};

struct bg_context
{
	struct bgmodel bgm;
	int compcount;
	double thres;
	double lrate;
	struct nd_image *initframes;
	int initframescount;
	int initframesstep;
	int framestoinit;
	int c;
};

struct bg_rect {
	int x0;
	int y0;
	int x1;
	int y1;
};

int bg_allocbgmodel(struct bgmodel *bgm, int w, int h);

int bg_initbgmodel(struct nd_image *img, int imgcnt, struct bgmodel *bgm,
	int compcount, double thres, double lrate);

int bg_updatebgmodel(const struct nd_image *img, struct bgmodel *bgm);

int bg_isforegroundtolowres(int *isfg, int w, int h,
	int *isfglowres, int lrw, int lrh);

int bg_isforegroundrectgroup(int *isforeground, int w, int h,
	struct bg_rect **r, int *rc);

int bg_initcontext(struct bg_context *ctx,
	int compcount, double thres, double lrate, int bgmw, int bgmh,
	int initframescount, int initframesstep);

int bg_putimgtocontext(struct bg_context *ctx, struct nd_image *img,
	struct bg_rect **r, int *rc);

int bg_destroycontext(struct bg_context *ctx);

#endif
