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

const int threadindex[4][6] = {{0, 4, 5},
				{1, 17, 16, 15, 14, 13},
				{2, 12, 11, 10, 9},
				{3, 8, 7, 6}};

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

int imgcheck(struct hc_hcascade *hc, const struct nd_image *img, double d,
	struct hc_rect **r, int *rc, const struct hc_scanconfig *conf)
{
	struct nd_image scaledimg;
	struct nd_image scaledimgsq;
	double *sd;
	int x, y;
	int rn;
	
	if (nd_imgscalebilinear(img, d, d, &scaledimg) < 0) {		
		return (-1);
	}

	if (nd_imgcreate(&scaledimgsq, scaledimg.w, scaledimg.h,
		img->format) < 0) {
		nd_imgdestroy(&scaledimg);
		return (-1);
	}

	for (y = 0; y < scaledimgsq.h; ++y)
		for (x = 0; x < scaledimgsq.w; ++x) {
			scaledimgsq.data[y * scaledimgsq.w + x]
				= scaledimg.data[y * scaledimg.w + x]
				* scaledimg.data[y * scaledimg.w + x];
	}

	if (hc_imgintegral(&scaledimg) < 0) {
		nd_imgdestroy(&scaledimg);
		nd_imgdestroy(&scaledimgsq);
	
		return (-1);
	}

	if (hc_imgintegral(&scaledimgsq) < 0) {
		nd_imgdestroy(&scaledimg);
		nd_imgdestroy(&scaledimgsq);
	
		return (-1);
	}

	if ((sd = malloc(sizeof(double) * scaledimg.h * scaledimg.w))
		== NULL) {
		nd_imgdestroy(&scaledimg);
		nd_imgdestroy(&scaledimgsq);		
	
		return (-1);
	}

	for (y = 0; y < scaledimg.h - hc->wh; y += 1)
		for (x = 0; x < scaledimg.w - hc->ww; x += 1) {
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
			
			sx1y1 = scaledimg.data[imgy1
				* scaledimg.w + imgx1];
		
			if (imgx0 < 0 || imgy0 < 0) {
				sx1y0 = (imgy0 >= 0)
					? scaledimg.data[imgy0
					* scaledimg.w + imgx1] : 0.0;
				sx0y1 = (imgx0 >= 0)
					? scaledimg.data[imgy1
					* scaledimg.w + imgx0] : 0.0;
				sx0y0 = 0.0;
			}
			else {
				sx1y0 = scaledimg.data[imgy0
					* scaledimg.w + imgx1];
				sx0y1 = scaledimg.data[imgy1
					* scaledimg.w + imgx0];
				sx0y0 = scaledimg.data[imgy0
					* scaledimg.w + imgx0];
			}			
			
			pixsum = sx0y0 + sx1y1 - sx1y0 - sx0y1;		

			sx1y1 = scaledimgsq.data[imgy1
				* scaledimgsq.w + imgx1];
		
			if (imgx0 < 0 || imgy0 < 0) {
				sx1y0 = (imgy0 >= 0)
					? scaledimgsq.data[imgy0
					* scaledimgsq.w + imgx1] : 0.0;
				sx0y1 = (imgx0 >= 0)
					? scaledimgsq.data[imgy1
					* scaledimgsq.w + imgx0] : 0.0;
				sx0y0 = 0.0;
			}
			else {
				sx1y0 = scaledimgsq.data[imgy0
					* scaledimgsq.w + imgx1];
				sx0y1 = scaledimgsq.data[imgy1
					* scaledimgsq.w + imgx0];
				sx0y0 = scaledimgsq.data[imgy0
					* scaledimgsq.w + imgx0];
			}			

			sqpixsum = sx0y0 + sx1y1 - sx1y0 - sx0y1;

			sd[y * scaledimg.w + x] = sqrt(1.0
				/ (pixcount - 1.0) * (sqpixsum
				- pixsum * pixsum / pixcount));
		}

	*r = NULL;
	*rc = 0;

	if (hc_fastimgscan(&scaledimg, sd, hc, r, rc, conf) < 0) {
		free(sd);
		nd_imgdestroy(&scaledimg);
	
		return -1;
	}

	for (rn = 0; rn < *rc; ++rn) {
		(*r)[rn].x0 /= d;
		(*r)[rn].y0 /= d;
		(*r)[rn].x1 /= d;
		(*r)[rn].y1 /= d;
	}
	
	free(sd);
	nd_imgdestroy(&scaledimg);
	nd_imgdestroy(&scaledimgsq);
		
	return 0;
}

int hc_imgpyramidscan(struct hc_hcascade *hc, const struct nd_image *img,
	struct hc_rect **newr, int *newrc, const struct hc_scanconfig *conf)
{
	double d;
	
	struct hc_rect *r;
	int rc;

	r = NULL;
	rc = 0;

	d = 1.0;
	do {
		struct hc_rect *tmpr;
		int tmprc;

		if ((int) ceil(img->w * d) < hc->ww
			|| (int) ceil(img->h * d) < hc->wh)
			break;
	
		if (imgcheck(hc, img, d, &tmpr, &tmprc, conf) < 0)
			return (-1);
		
		r = realloc(r, sizeof(struct hc_rect) * (rc + tmprc));
		memcpy(r + rc, tmpr, sizeof(struct hc_rect) * tmprc);
		rc += tmprc;

		d *= conf->scalestep;
	} while (1);

	*newr = NULL;
	*newrc = 0;

	if (rc > 0 && hc_conrect(r, rc, newr, newrc) < 0)
		return (-1);

	return 0;
}
