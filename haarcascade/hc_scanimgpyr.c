#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <pthread.h>

#include <cairo.h>

#include "nd_image.h"
#include "nd_error.h"
#include "hc_hcascade.h"
#include "hc_rect.h"

#include "hc_scanimgpyr.h"

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))

int hc_confbuild(struct hc_scanconfig *conf, int ww, int wh,
	int w, int h, double d, int wstep, int hstep, int tc)
{
	int n;
	double a;
	int b, e;
	
	conf->scalestep = d;
	conf->winwstep = wstep;
	conf->winhstep = hstep;
	
	conf->tc = 0;
	
	n = floor(log(MAX((double) ww / w, (double) wh / h)) / log(d));

	d = pow(d, 2.0);
	
	a = (1 - pow(d, n)) / tc / (1 - d);

	conf->map = malloc(sizeof(int *) * (n + 1));

	conf->tc = 0;
	b = 0;
	e = n;

	while (b < e) {
		double s;
		int c;

		c = 1;
		conf->map[conf->tc] = malloc(sizeof(int) * c);
		conf->map[conf->tc][c - 1] = b;

		s = pow(d, b);
		
		while(s + pow(d, e) < a && b < e) {
			conf->map[conf->tc] = realloc(conf->map[conf->tc],
				sizeof(int) * (++c));
			conf->map[conf->tc][c - 1] = e;
			
			s += pow(d, e);
			--e;
		}

		conf->map[conf->tc] = realloc(conf->map[conf->tc],
			sizeof(int) * (++c));
		conf->map[conf->tc][c - 1] = -1;
	
		++(conf->tc);
		++b;
	}
	
	return 0;
};

static int hc_fastimgscan(struct nd_image *img, double *sd,
	struct hc_hcascade *hc, struct hc_rect **r, int *roffset,
	const struct hc_scanconfig *conf)
{
	struct nd_image imginwin;
	int x;
	int y;
	for (y = 1; y < img->h - hc->wh; y += conf->winhstep)
		for (x = 1; x < img->w - hc->ww; x += conf->winwstep) {
			imginwin.w = hc->ww + 1;
			imginwin.h = hc->wh + 1;
			imginwin.format = ND_PF_GRAYSCALE;
			imginwin.data = img->data + (y - 1) * img->w + (x - 1);

			if (hc_fastimgclassify(hc, &imginwin,
				img->w, sd[y * img->w + x])) {
				++(*roffset);
			
				if ((*r = realloc(*r, sizeof(struct hc_rect)
					* (*roffset))) == NULL)
					return (-1);
			
				(*r)[*roffset - 1].x0 = (double) x;
				(*r)[*roffset - 1].y0 = (double) y;
				(*r)[*roffset - 1].x1 = (double) (x + hc->ww);
				(*r)[*roffset - 1].y1 = (double) (y + hc->wh);
			}
		}
	
	return 0;
}

static int imgcheck(struct hc_hcascade *hc, struct nd_image *img,
	struct hc_rect **r, int *rc, const struct hc_scanconfig *conf)
{
	struct nd_image imgsq;
	double *sd;
	int x, y;
	
	if (nd_imgcreate(&imgsq, img->w, img->h,
		img->format) < 0) {
		return (-1);
	}

	for (y = 0; y < imgsq.h; ++y)
		for (x = 0; x < imgsq.w; ++x) {
			imgsq.data[y * imgsq.w + x] = img->data[y * img->w + x]
				* img->data[y * img->w + x];
	}

	if (hc_imgintegral(img) < 0) {
		nd_imgdestroy(&imgsq);
		return (-1);
	}

	if (hc_imgintegral(&imgsq) < 0) {
		nd_imgdestroy(&imgsq);
		return (-1);
	}

	if ((sd = malloc(sizeof(double) * img->h * img->w))
		== NULL) {
		nd_imgdestroy(&imgsq);		
		return (-1);
	}

	for (y = 0; y < img->h - hc->wh; y += 1)
		for (x = 0; x < img->w - hc->ww; x += 1) {
			double pixsum;
			double sqpixsum;
			int pixcount;
			int imgx0;
			int imgy0;
			int imgx1;
			int imgy1;

			double sx1y0;
			double sx0y1;
			double sx0y0;;
			double sx1y1;

			pixcount = hc->ww * hc->wh;

			pixsum = 0.0;
			sqpixsum = 0.0;

			imgx0 = x - 1;
			imgy0 = y - 1;	
			imgx1 = x + hc->ww - 1;
			imgy1 = y + hc->wh - 1;
			
			sx1y1 = img->data[imgy1
				* img->w + imgx1];
		
			if (imgx0 < 0 || imgy0 < 0) {
				sx1y0 = (imgy0 >= 0)
					? img->data[imgy0
					* img->w + imgx1] : 0.0;
				sx0y1 = (imgx0 >= 0)
					? img->data[imgy1
					* img->w + imgx0] : 0.0;
				sx0y0 = 0.0;
			}
			else {
				sx1y0 = img->data[imgy0
					* img->w + imgx1];
				sx0y1 = img->data[imgy1
					* img->w + imgx0];
				sx0y0 = img->data[imgy0
					* img->w + imgx0];
			}			
			
			pixsum = sx0y0 + sx1y1 - sx1y0 - sx0y1;		

			sx1y1 = imgsq.data[imgy1
				* imgsq.w + imgx1];
		
			if (imgx0 < 0 || imgy0 < 0) {
				sx1y0 = (imgy0 >= 0)
					? imgsq.data[imgy0
					* imgsq.w + imgx1] : 0.0;
				sx0y1 = (imgx0 >= 0)
					? imgsq.data[imgy1
					* imgsq.w + imgx0] : 0.0;
				sx0y0 = 0.0;
			}
			else {
				sx1y0 = imgsq.data[imgy0
					* imgsq.w + imgx1];
				sx0y1 = imgsq.data[imgy1
					* imgsq.w + imgx0];
				sx0y0 = imgsq.data[imgy0
					* imgsq.w + imgx0];
			}			

			sqpixsum = sx0y0 + sx1y1 - sx1y0 - sx0y1;

			sd[y * img->w + x] = sqrt(1.0
				/ (pixcount - 1.0) * (sqpixsum
				- pixsum * pixsum / pixcount));
		}

	*r = NULL;
	*rc = 0;

	if (hc_fastimgscan(img, sd, hc, r, rc, conf) < 0) {
		free(sd);
		return -1;
	}


	free(sd);
	nd_imgdestroy(&imgsq);
		
	return 0;
}

struct hc_pyrcheckarg {
	struct hc_hcascade *hc;
	const struct hc_scanconfig *conf;
	
	const struct nd_image *img;
	const int *map;
};

struct hc_pyrcheckret {
	struct hc_rect *r;
	int rc;
};

void *pyrcheckfunc(void *a)
{
	struct hc_pyrcheckarg *arg;
	struct hc_pyrcheckret *ret;
	const int *curidx;
	
	arg = a;

	ret = malloc(sizeof(struct hc_pyrcheckret));
	ret->r = NULL;
	ret->rc = 0;

	curidx = arg->map;
	while (*curidx >= 0) {
		struct hc_rect *tmpr;
		int tmprc;
		int rn;

		struct nd_image scaledimg;

		if (nd_imgscalebilinear(arg->img,
			pow(arg->conf->scalestep, *curidx),
			pow(arg->conf->scalestep, *curidx), &scaledimg) < 0)
			return NULL;

		if (imgcheck(arg->hc, &scaledimg, &tmpr, &tmprc, arg->conf)
			< 0)
			return NULL;

		
		for (rn = 0; rn < tmprc; ++rn) {
			tmpr[rn].x0 /= pow(arg->conf->scalestep, *curidx);
			tmpr[rn].y0 /= pow(arg->conf->scalestep, *curidx);
			tmpr[rn].x1 /= pow(arg->conf->scalestep, *curidx);
			tmpr[rn].y1 /= pow(arg->conf->scalestep, *curidx);
		}

		ret->r = realloc(ret->r, sizeof(struct hc_rect)
			* (ret->rc + tmprc));
		memcpy(ret->r + ret->rc, tmpr, sizeof(struct hc_rect) * tmprc);
		ret->rc += tmprc;

		nd_imgdestroy(&scaledimg);
		
		++curidx;
	}

	return ret;
}

int hc_imgpyramidscan(struct hc_hcascade *hc, const struct nd_image *img,
	struct hc_rect **newr, int *newrc, const struct hc_scanconfig *conf)
{
	struct hc_rect *r;
	int rc;
	
	r = NULL;
	rc = 0;

	pthread_t *thread;
	struct hc_pyrcheckarg *threadarg;
	struct hc_pyrcheckret **threadret;
	int i;

	thread = malloc(sizeof(pthread_t) * conf->tc);
	threadarg = malloc(sizeof(struct hc_pyrcheckarg) * conf->tc);
	threadret = malloc(sizeof(struct hc_rect *) * conf->tc);

	for (i = 0; i < conf->tc; ++i) {
		threadarg[i].hc = hc;
		threadarg[i].conf = conf;
		threadarg[i].img = img;
		threadarg[i].map = conf->map[i];

		pthread_create(thread + i, NULL, &pyrcheckfunc,
			(void *) (threadarg + i));
	}

	for (i = 0; i < conf->tc; ++i)
		pthread_join(thread[i], (void **) (threadret + i));

	for (i = 0; i < conf->tc; ++i) {
		r = realloc(r, sizeof(struct hc_rect)
			* (rc + threadret[i]->rc));
		memcpy(r + rc, threadret[i]->r, sizeof(struct hc_rect)
			* threadret[i]->rc);
		rc += threadret[i]->rc;
	}


	free(thread);
	free(threadarg);
	
	for (i = 0; i < conf->tc; ++i) {
		free(threadret[i]->r);
		free(threadret[i]);
	}
	free(threadret);

	*newr = NULL;
	*newrc = 0;
	if (rc > 0 && hc_conrect(r, rc, newr, newrc) < 0)
		return (-1);


	return 0;
}
