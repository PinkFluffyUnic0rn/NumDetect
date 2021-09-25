#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <pthread.h>

#include "nd_error.h"

#include "bgmodel.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define THREADCOUNT 8

int bg_allocbgmodel(struct bgmodel *bgm, int w, int h)
{
	if ((bgm->gm = malloc(sizeof(struct em_gaussmix) * w * h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto gmmallocerror;
	}

	if ((bgm->isforeground = malloc(sizeof(int) * w * h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto isforegroundmallocerror;
	}

	bgm->w = w;
	bgm->h = h;

	return 0;

isforegroundmallocerror:
	free(bgm->gm);
gmmallocerror:
	return (-1);
}

static int bg_getcellsamples(struct nd_vector3 *sample, struct nd_image *img,
	int imgcnt, int cellx, int celly, int bgwidth, int bgheight)
{
	int idx;
	int x0, x1, y0, y1;
	int i, x, y;

	idx = 0;

	for (i = 0; i < imgcnt; ++i) {
		struct nd_image *curimg;

		curimg = img + i;

		x0 = cellx * curimg->w / bgwidth;
		x1 = (cellx + 1) * curimg->w / bgwidth;
		y0 = celly * curimg->h / bgheight;
		y1 = (celly + 1) * curimg->h / bgheight;
	
		for (y = y0; y < y1; ++y) {
			for (x = x0; x < x1; ++x) {
				int pixidx;
				
				pixidx = y * curimg->w + x;
			
				sample[idx].x = curimg->data[pixidx * 3 + 0];
				sample[idx].y = curimg->data[pixidx * 3 + 1];
				sample[idx].z = curimg->data[pixidx * 3 + 2];
				
				++idx;
			}
		}
	}

	return idx;
}

static int bg_printexpectation(struct bgmodel *bgm, const char *path)
{
	struct nd_image outimg;
	int i, j;

	nd_imgcreate(&outimg, bgm->w, bgm->h, ND_PF_RGB);

	for (i = 0; i < bgm->h * bgm->w; ++i) {
		struct nd_vector3 expectation;
		struct em_gaussmix *gm;

		gm = bgm->gm + i;
		
		expectation.x = expectation.y = expectation.z = 0.0;
		for (j = 0; j < gm->pn; ++j) {
			expectation.x += gm->w[j] * gm->m[j].x;
			expectation.y += gm->w[j] * gm->m[j].y;
			expectation.z += gm->w[j] * gm->m[j].z;
		}

		outimg.data[i * 3 + 0] = expectation.x;
		outimg.data[i * 3 + 1] = expectation.y;
		outimg.data[i * 3 + 2] = expectation.z;
	}

	nd_imgwrite(&outimg, path);
	nd_imgdestroy(&outimg);

	return 0;
}

int bg_initbgmodel(struct nd_image *img, int imgcnt, struct bgmodel *bgm,
	int compcount, double thres, double lrate)
{
	int i;	
	struct nd_vector3 *sample;

	for (i = 0; i < imgcnt; ++i)
		assert(img[i].format == ND_PF_RGB);

	bgm->compcount = compcount;
	bgm->thres = thres;
	bgm->lrate = lrate;

	int maxsamples = 0;
	for (i = 0; i < imgcnt; ++i) {
		int cnt;

		cnt = (img[i].w / bgm->w + 1) * (img[i].h / bgm->h + 1);

		maxsamples = (maxsamples < cnt) ? cnt : maxsamples;
	}

	if ((sample = malloc(sizeof(struct nd_vector3) * maxsamples * imgcnt))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto mallocxerror;
	}
	
	for (i = 0; i < bgm->w * bgm->h; ++i) {
		int samplecnt;

		samplecnt = bg_getcellsamples(sample, img, imgcnt,
			i % bgm->w, i / bgm->w, bgm->w, bgm->h);

		em_allocgaussmix(bgm->gm + i, bgm->compcount);
		
		if (em_gauss3d(sample, samplecnt, bgm->gm + i, bgm->compcount,
			NULL, NULL, 0.1) < 0)
			goto emgauss3derror;

		printf("%1.3f%%\n", i / (double)(bgm->w * bgm->h) * 100.0);
	}

	free(sample);

	return 0;

emgauss3derror:
	free(sample);
mallocxerror:
	return (-1);
}

static int bg_matchtest(struct nd_vector3 *pix, struct em_gaussmix *gm)
{
	int i;

	for (i = 0; i < gm->pn; ++i) {
		struct nd_vector3 md;
		
		nd_v3lincomb(pix, -1.0, gm->m + i, &md);

		if (((pow(md.x, 2.0) + pow(md.y, 2.0) + pow(md.z, 2.0))
			< 6.25 * gm->c[i]))
			return i;
	}

	return (-1);
}

static int bg_leastprobcomp(struct em_gaussmix *gm)
{
	int lpcompidx;
	int i;

	lpcompidx = 0;

	for (i = 0; i < gm->pn; ++i)
		lpcompidx = (gm->w[i] < gm->w[lpcompidx]) ? i : lpcompidx;

	return lpcompidx;
}

static int bg_updategausscomp(double *w, struct nd_vector3 *m, double *c,
	struct nd_vector3 *pix, int matched, double lrate)
{
	struct nd_vector3 md;
	double p;	

	if (matched == 0) {
		*w = (1.0 - lrate) * *w;
		return 0;
	}

	*w = (1.0 - lrate) * *w + lrate;

	// muliple of pow(EM_MINVARIANCE, 1.5) is
	// needed for p to stay in  range of [0, 1]
	p = lrate * em_probability_gauss3d(m, *c, pow(*c, 3.0), pix)
		* pow(EM_MINVARIANCE, 1.5);

	m->x = (1.0 - p) * m->x + p * pix->x;
	m->y = (1.0 - p) * m->y + p * pix->y;
	m->z = (1.0 - p) * m->z + p * pix->z;

	nd_v3lincomb(pix, -1.0, m, &md);

	*c = (1.0 - p) * *c
		+ p * (pow(md.x, 2.0) + pow(md.y, 2.0) + pow(md.z, 2.0));
	*c = (*c < EM_MINVARIANCE) ? EM_MINVARIANCE : *c;

	return 0;
}

static int sortgm(struct em_gaussmix *gm)
{
	struct nd_vector3 tmpv;
	double tmps;
	int i, j;

	for (i = 0; i < gm->pn; ++i) {
		int maxp;
		double maxr;

		maxp = i;
		maxr = gm->w[i] / gm->c[i];
		for (j = i + 1; j < gm->pn; ++j) {
			double r;

			r = gm->w[j] / gm->c[j];
			if (r > maxr) {
				maxr = r;
				maxp = j;
			}
		}

		tmps = gm->w[i];
		gm->w[i] = gm->w[maxp];
		gm->w[maxp] = tmps;

		tmpv = gm->m[i];
		gm->m[i] = gm->m[maxp];
		gm->m[maxp] = tmpv;

		tmps = gm->c[i];
		gm->c[i] = gm->c[maxp];
		gm->c[maxp] = tmps;
	}

	return 0;
}

struct bg_updatethreadarg
{
	int p0;
	int p1;
	const struct nd_image *img;
	struct bgmodel *bgm;
};

static void *updatebgfunc(void *arg)
{
	struct bg_updatethreadarg *a;
	int i, j;

	a = arg;

	for (i = a->p0; i < a->p1; ++i) {
		struct nd_vector3 pix;
		struct em_gaussmix *gm;
		int cellx, celly;
		int gi;

		cellx = (i % a->img->w) * a->bgm->w / a->img->w;
		celly = (i / a->img->w) * a->bgm->h / a->img->h;

		pix.x = a->img->data[i * 3 + 0];
		pix.y = a->img->data[i * 3 + 1];
		pix.z = a->img->data[i * 3 + 2];

		gm = a->bgm->gm + (celly * a->bgm->w + cellx);

		if ((gi = bg_matchtest(&pix, gm)) < 0) {
			int lpcompidx;
		
			a->bgm->isforeground[celly * a->bgm->w + cellx] = 1;

			lpcompidx = bg_leastprobcomp(gm);

			gm->m[lpcompidx] = pix;
			gm->c[lpcompidx] = 1e9;
		}
		else {
			double wsum;
			int c;
	
			wsum = 0.0;
			c = 0;
			while (wsum < a->bgm->thres)
				wsum += gm->w[c++];

			a->bgm->isforeground[celly * a->bgm->w + cellx]
				= (gi < c) ? 0 : 1;

			wsum = 0.0;
			for (j = 0; j < gm->pn; ++j)
				wsum += gm->w[j];

			for (j = 0; j < gm->pn; ++j)
				gm->w[j] /= wsum;

			for (j = 0; j < gm->pn; ++j)
				bg_updategausscomp(gm->w + j,
					gm->m + j, gm->c + j, &pix,
					(j == gi) ? 1 : 0, a->bgm->lrate);
		}

		sortgm(gm);
	}

	return NULL;
}

static int bg_isfgfilter(int **isfg, int w, int h)
{
	int x, y;
	int *newisforeground;

	if ((newisforeground = malloc(sizeof(int) * w * h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (y = 0; y < h; ++y)
		for (x = 0; x < w; ++x) {
			int c;

			if (y == 0 || x == 0
				|| y == (h - 1) || x == (w - 1)) {
				newisforeground[y * w + x] = 0;
				continue;
			}

			c = 0;

			if ((*isfg)[(y - 1) * w + (x - 1)])
				++c;
			if ((*isfg)[(y - 1) * w + x])
				++c;
			if ((*isfg)[(y - 1) * w + (x + 1)])
				++c;
			if ((*isfg)[y * w + (x - 1)])
				++c;
			if ((*isfg)[y * w + x])
				++c;
			if ((*isfg)[y * w + (x + 1)])
				++c;
			if ((*isfg)[(y + 1) * w + (x - 1)])
				++c;
			if ((*isfg)[(y + 1) * w + x])
				++c;
			if ((*isfg)[(y + 1) * w + (x + 1)])
				++c;

			newisforeground[y * w + x] = (c > 4) ? 1 : 0;
		}

	free(*isfg);
	*isfg = newisforeground;

	return 0;
}

int bg_updatebgmodel(const struct nd_image *img, struct bgmodel *bgm)
{
	pthread_t *thread;
	struct bg_updatethreadarg *arg;
	int i;

	assert(img != NULL && bgm != NULL);
	assert(img->format == ND_PF_RGB);
	
	if ((thread = malloc(sizeof(pthread_t) * THREADCOUNT)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto threadmallocerror;
	}

	if ((arg = malloc(sizeof(struct bg_updatethreadarg) * THREADCOUNT))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto argmallocerror;
	}

	for (i = 0; i < THREADCOUNT; ++i) {
		arg[i].p0 = i * img->w * img->h / THREADCOUNT;
		arg[i].p1 = (i + 1) * img->w * img->h / THREADCOUNT;
		arg[i].img = img;
		arg[i].bgm = bgm;
	}

	for (i = 0; i < THREADCOUNT; ++i)
		pthread_create(thread + i, NULL, &updatebgfunc,
			(void *)(arg + i));
		
	for (i = 0; i < THREADCOUNT; ++i)
		pthread_join(thread[i], NULL);

	free(arg);

	if (bg_isfgfilter(&(bgm->isforeground), bgm->w, bgm->h))
		return (-1);

///////////////////////////////////////////////////////////////////////////////
/*
	static int cc = 0;
	char path[255];

	sprintf(path, "res/%d.png", cc++);

	bg_printexpectation(bgm, path);
*/
///////////////////////////////////////////////////////////////////////////////

	return 0;

argmallocerror:
	free(thread);
threadmallocerror:
	return (-1);
}

int bg_isforegroundtolowres(int *isfg, int w, int h,
	int *isfglowres, int lrw, int lrh)
{
	int x, y;
	
	memset(isfglowres, 0, sizeof(int) * lrw * lrh);

	for (y = 0; y < h; ++y)
		for (x = 0; x < w; ++x) {
			int lrx, lry;
			int *isfglr;

			lry = lrh * y / h;
			lrx = lrw * x / w;
			
			isfglr = isfglowres + lry * lrw + lrx;
		
			*isfglr = (isfg[y * w + x] == 1) ? 1 : *isfglr;
		}

	return 0;
}

int bg_isforegroundrectgroup(int *isforeground, int w, int h,
	struct bg_rect **r, int *rc)
{
	int maxrc;
	int i;
	int *isfg;
	int *stack;	
	int stacksize;

	if (((*r) = malloc(sizeof(struct bg_rect))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}
	(*rc) = 0;
	maxrc = 1;

	if ((stack = malloc(sizeof(int) * w * h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if ((isfg = malloc(sizeof(int) * w * h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	memcpy(isfg, isforeground, sizeof(int) * w * h);

	for (i = 0; i < h * w; ++i)
		if (isfg[i] == 1) {
			struct bg_rect *rect;
			
			stack[0] = i;
			stacksize = 1;

			rect = (*r) + (*rc);
			
			rect->x0 = w;
			rect->y0 = h;
			rect->x1 = 0;
			rect->y1 = 0;

			while (stacksize > 0) {
				int x, y;

				--stacksize;

				y = stack[stacksize] / w;
				x = stack[stacksize] % w;

				rect->x0 = (x < rect->x0) ? x : rect->x0;
				rect->y0 = (y < rect->y0) ? y : rect->y0;
				rect->x1 = (x > rect->x1) ? x : rect->x1;
				rect->y1 = (y > rect->y1) ? y : rect->y1;

				if (y > 0 && x > 0 && isfg[(y - 1) * w + (x - 1)] == 1) {
					stack[stacksize++] = (y - 1) * w + (x - 1);
					isfg[(y - 1) * w + (x - 1)] = 0;
				}

				if (y > 0 && isfg[(y - 1) * w + x] == 1) {
					stack[stacksize++] = (y - 1) * w + x;
					isfg[(y - 1) * w + x] = 0;
				}

				if (y > 0 && (x + 1) < w && isfg[(y - 1) * w + (x + 1)] == 1) {
					stack[stacksize++] = (y - 1) * w + (x + 1);
					isfg[(y - 1) * w + (x + 1)] = 0;
				}
			
				if (x > 0 && isfg[y * w + (x - 1)] == 1) {
					stack[stacksize++] = y * w + (x - 1);
					isfg[y * w + (x - 1)] = 0;
				}
				if ((x + 1) < w && isfg[y * w + (x + 1)] == 1) {
					stack[stacksize++] = y * w + (x + 1);
					isfg[y * w + (x + 1)] = 0;
				}

				if ((y + 1) < h && x > 0 && isfg[(y + 1) * w + (x - 1)] == 1) {
					stack[stacksize++] = (y + 1) * w + (x - 1);
					isfg[(y + 1) * w + (x - 1)] = 0;
				}

				if ((y + 1) < h && isfg[(y + 1) * w + x] == 1) {
					stack[stacksize++] = (y + 1) * w + x;
					isfg[(y + 1) * w + x] = 0;
				}

				if ((y + 1) < h && (x + 1) < w && isfg[(y + 1) * w + (x + 1)] == 1) {
					stack[stacksize++] = (y + 1) * w + (x + 1);
					isfg[(y + 1) * w + (x + 1)] = 0;
				}
			}

			++(*rc);
			if (*rc >= maxrc) {
				maxrc *= 2;
				if (((*r) = realloc(*r,
					sizeof(struct bg_rect) * maxrc))
					== NULL) {
					nd_seterrormessage(ND_MSGALLOCERROR,
						__func__);
					return (-1);
				}
			}
		}

	free(stack);
	free(isfg);

	return 0;
}

int bg_initcontext(struct bg_context *ctx,
	int compcount, double thres, double lrate, int bgmw, int bgmh,
	int initframescount, int initframesstep)
{
	if (bg_allocbgmodel(&(ctx->bgm), bgmw, bgmh) < 0)
		return (-1);

	ctx->compcount = compcount;
	ctx->thres = thres;
	ctx->lrate = lrate;

	ctx->initframescount = initframescount;
	ctx->initframesstep = initframesstep;

	ctx->framestoinit = initframescount;
	ctx->c = 0;

	if ((ctx->initframes
		= malloc(sizeof(struct nd_image) * initframescount)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	return 0;	
}

int bg_putimgtocontext(struct bg_context *ctx, struct nd_image *img,
	struct bg_rect **r, int *rc)
{
	assert(img->w == ctx->bgm.w && img->h == ctx->bgm.h
		&& img->format == ND_PF_RGB);

	*rc = 0;

	ctx->c = (ctx->c < ctx->initframesstep) ? ctx->c + 1 : 0;

	if (ctx->c == 0 && ctx->framestoinit > 0) {
		--(ctx->framestoinit);

		nd_imgcopy(img, ctx->initframes + ctx->framestoinit);
	}
	if (ctx->framestoinit == 0) {
		--(ctx->framestoinit);

		if (bg_initbgmodel(ctx->initframes, ctx->initframescount,
			&(ctx->bgm), ctx->compcount, ctx->thres, ctx->lrate)
			< 0)
			return (-1);

		free(ctx->initframes);
	}
	if (ctx->framestoinit < 0) {
		int *isfglowres;
		int lrw, lrh;
		int sqrsz;
		int i;

		if (bg_updatebgmodel(img, &(ctx->bgm)) < 0)
			return (-1);

		sqrsz = MAX(ctx->bgm.w, ctx->bgm.h) / 25;

		lrw = ctx->bgm.w / sqrsz;
		lrh = ctx->bgm.h / sqrsz;

		if ((isfglowres = malloc(sizeof(double) * lrw * lrh))
			== NULL) {
			nd_seterrormessage(ND_MSGALLOCERROR, __func__);
			return (-1);
		}

		if (bg_isforegroundtolowres(ctx->bgm.isforeground,
			ctx->bgm.w, ctx->bgm.h, isfglowres, lrw, lrh) < 0)
			return (-1);
		
		if (bg_isforegroundrectgroup(isfglowres, lrw, lrh, r, rc) < 0)
			return (-1);
	
		for (i = 0; i < *rc; ++i) {
			struct bg_rect *rect;

			rect = (*r) + i;
			
			rect->y0 = rect->y0 * ctx->bgm.h / lrh;
			rect->y1 = (rect->y1 + 1) * ctx->bgm.h / lrh;
			rect->x0 = rect->x0 * ctx->bgm.w / lrw;
			rect->x1 = (rect->x1 + 1) * ctx->bgm.w / lrw;
		}

		free(isfglowres);
	}

	return 0;
}

int bg_destroycontext(struct bg_context *ctx)
{

	return 0;
}
