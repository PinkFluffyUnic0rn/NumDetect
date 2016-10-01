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

int imgscan(struct nd_image *img, struct hc_hcascade *hc,
	struct hc_rect **r, int *roffset, int *rmax)
{
	struct nd_image imginwin;
	int x;
	int y;

	nd_imgcreate(&imginwin, hc->ww, hc->wh, img->chans);

	for (y = 0; y < img->h - hc->wh; y += 1)//hc->ww / 3)
		for (x = 0; x < img->w - hc->ww; x += 1) {//hc->wh / 3) {
			if (nd_imgcrop(img, x, y,
				hc->ww, hc->wh, &imginwin) < 0)
				return (-1);

			if (nd_imgnormalize(&imginwin, 0, 1) < 0)
				return (-1);
			
			if (hc_imgintegral(&imginwin) < 0)
				return (-1);


			if (hc_imgclassify(hc, &imginwin)) {	
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

	nd_imgdestroy(&imginwin);
	
	return 0;
}

int fastimgscan(struct nd_image *img, double *sd, struct hc_hcascade *hc,
	struct hc_rect **r, int *roffset, int *rmax)
{
	struct nd_image imginwin;
	int x;
	int y;

	nd_imgcreate(&imginwin, hc->ww, hc->wh, img->chans);

	for (y = 0; y < img->h - hc->wh; y += 1)//hc->ww / 3)
		for (x = 0; x < img->w - hc->ww; x += 1) {//hc->wh / 3) {
			int pixn;
			
			if (nd_imgcrop(img, x, y,
				hc->ww, hc->wh, &imginwin) < 0)
				return (-1);
			
			for (pixn = 0; pixn < imginwin.w * imginwin.h; ++pixn)
				imginwin.data[pixn] /= sd[y * img->w + x];

			if (hc_imgintegral(&imginwin) < 0)
				return (-1);


			if (hc_imgclassify(hc, &imginwin)) {	
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

	nd_imgdestroy(&imginwin);
	
	return 0;
}

int imgpyramidscan(struct hc_hcascade *hc, struct nd_image *img,
	struct hc_rect **newr, int *newrc)
{
	double d;
	
	struct hc_rect *r;
	int rmax;
	int rc;
	int rn;

	r = NULL;
	rc = 0;
	rmax = 0;

	d = 1.0;	
	do {
		struct nd_image scaledimg;
		int prevrc;

		if (nd_imgscale(img, d, d, &scaledimg))
			exit(2);
	
		printf("scale = %f\n", d);
		if (scaledimg.w < hc->ww || scaledimg.h < hc->wh)
			break;

		double *sd;
		int x, y;

		sd = malloc(sizeof(double) * scaledimg.h * scaledimg.w);

		for (y = 0; y < scaledimg.h - hc->wh; y += 1)
			for (x = 0; x < scaledimg.w - hc->ww; x += 1) {
				double pixsum;
				double sqpixsum;
				int pixcount;
				int yy, xx;
			
				pixcount = hc->ww * hc->wh;

				pixsum = 0.0;
				sqpixsum = 0.0;

				for (yy = y; yy < y + hc->wh; ++yy)
					for (xx = x; xx < x + hc->ww; ++xx) {
						pixsum += scaledimg.data[yy * scaledimg.w + xx];
						sqpixsum += scaledimg.data[yy * scaledimg.w + xx]
							* scaledimg.data[yy * scaledimg.w + xx];
					}

				sd[y * scaledimg.w + x] = sqrt(1.0 / (pixcount - 1.0)
					* (sqpixsum - pixsum * pixsum / pixcount));

//				printf("%f\n", sd[y * scaledimg.w + x]);
			}	
/*
		if (hc_imgintegral(&scaledimg) < 0)
				return (-1);
*/

		prevrc = rc;
		if (fastimgscan(&scaledimg, sd, hc, &r, &rc, &rmax) < 0) {
			fprintf(stderr, "nd_imgscan: %s.\n",
				nd_strerror(nd_error));
		
			return -1;
		}
	
		for (rn = prevrc; rn < rc; ++rn) {
			r[rn].x0 /= d;
			r[rn].y0 /= d;
			r[rn].x1 /= d;
			r[rn].y1 /= d;
		}

		free(sd);

/*
		prevrc = rc;
		if (imgscan(&scaledimg, hc, &r, &rc, &rmax) < 0) {
			fprintf(stderr, "nd_imgscan: %s.\n",
				nd_strerror(nd_error));
		
			return -1;
		}
	
		for (rn = prevrc; rn < rc; ++rn) {
			r[rn].x0 /= d;
			r[rn].y0 /= d;
			r[rn].x1 /= d;
			r[rn].y1 /= d;
		}
*/
		d *= 0.9;
	} while (1);
	
	if (rc <= 0)
		return (-1);

	if (hc_conrect(r, rc, newr, newrc) < 0)
		return (-1);

	return 0;
}

int main(int argc, char **argv)
{	
	struct hc_hcascade hc;
	struct nd_image img;

	struct hc_rect *newr;
	int newrc;

	int rn;
	
	time_t t;

	if (hc_hcascaderead(&hc, argv[1]) < 0) {
		fprintf(stderr, "hc_hcascaderead: %s.\n",
			argv[1]);
		return 1;
	}
	
	if (nd_imgread(argv[2], &img) < 0) {
		fprintf(stderr, "nd_imgread: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	if (nd_imggrayscale(&img) < 0) {
		fprintf(stderr, "nd_imggrayscale: %s.\n",
			nd_strerror(nd_error));
		return 1;
	}

	t = time(0);

	imgpyramidscan(&hc, &img, &newr, &newrc);

	printf("%lu\n", time(0) - t);
	
	for (rn = 0; rn < newrc; ++rn) {
		char a[255];

		sprintf(a, "%s/%u.png", argv[3], rn);
		
		double rh = abs(newr[rn].y0 - newr[rn].y1);

		newr[rn].y0 = (int) (newr[rn].y0) - (rh * 0.15);
		newr[rn].y1 = newr[rn].y1 + (rh * 0.15);

		if (newr[rn].y0 >= 0 && newr[rn].y1 < img.h) {
			struct nd_image imginwin;

			nd_imgcreate(&imginwin, abs(newr[rn].x1 - newr[rn].x0),
				abs(newr[rn].y1 - newr[rn].y0), img.chans);
			
			nd_imgcrop(&img,
				newr[rn].x0, newr[rn].y0,
				abs(newr[rn].x1 - newr[rn].x0),
				abs(newr[rn].y1 - newr[rn].y0), &imginwin);
			
			nd_imgwrite(a, &imginwin);
		}
	}

	return 0;
}
