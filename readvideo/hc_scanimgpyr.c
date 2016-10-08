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

static int nd_fastimgscan(struct nd_image *img, double *sd,
	struct hc_hcascade *hc, struct hc_rect **r, int *roffset, int *rmax,
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
				if (*roffset >= *rmax) {
					*rmax = (*rmax > 0) ? (*rmax * 2) : 1;
	
					if ((*r = realloc(*r,
						sizeof(struct hc_rect)
						* *rmax)) == NULL)
						return (-1);
				}
			
				(*r)[*roffset].x0 = (double) x;
				(*r)[*roffset].y0 = (double) y;
				(*r)[*roffset].x1 = (double) (x + hc->ww);
				(*r)[*roffset].y1 = (double) (y + hc->wh);
	
				++(*roffset);
				
			}
		}

	return 0;
}

int nd_imgpyramidscan(struct hc_hcascade *hc, struct nd_image *img,
	struct hc_rect **newr, int *newrc, const struct hc_scanconfig *conf)
{
	double d;
	
	struct hc_rect *r;
	int rmax;
	int rc;
	int rn;

	r = NULL;
	rc = 0;
	rmax = 0;
	*newrc = 0;

	d = 1.0;
	do {
		struct nd_image scaledimg;
		struct nd_image scaledimgsq;
		double *sd;
		int prevrc;
		int x, y;

		if (nd_imgscalebilinear(img, d, d, &scaledimg) < 0) {
		
			return (-1);
		}

		if (scaledimg.w < hc->ww || scaledimg.h < hc->wh) {
			nd_imgdestroy(&scaledimg);
			break;
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

		prevrc = rc;
		if (nd_fastimgscan(&scaledimg, sd, hc, &r, &rc, &rmax, conf)
			< 0) {
			free(sd);
			nd_imgdestroy(&scaledimg);
			nd_imgdestroy(&scaledimgsq);

			return -1;
		}
		
		for (rn = prevrc; rn < rc; ++rn) {
			r[rn].x0 /= d;
			r[rn].y0 /= d;
			r[rn].x1 /= d;
			r[rn].y1 /= d;
		}

		d *= conf->scalestep;
		
		free(sd);
		nd_imgdestroy(&scaledimg);
		nd_imgdestroy(&scaledimgsq);
	} while (1);
	
	if (rc > 0 && hc_conrect(r, rc, newr, newrc) < 0)
		return (-1);

	return 0;
}
