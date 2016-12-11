#include <stdlib.h>
#include <stdio.h>

#include "nd_image.h"
#include "nd_error.h"
#include "hc_hcascade.h"
#include "nn_neuronet.h"

void imgtoinput(struct nd_image *img, double *input)
{
	int pixn;
	
	for (pixn = 0; pixn < img->w * img->h; ++pixn) {
		input[pixn] = 2.0 * img->data[pixn] - 1.0;
	}
}

int findpattern(struct nn_neuronet *nnet, double *input)
{
	double *out;
	uint minoutid;
	int i;
	uint outputc;
	
	outputc = nnet->nodec[nnet->levelc - 1];
	
	out = (double *) malloc(sizeof(double) * outputc);	
	
	nn_neuroneteval(nnet, input, out);

	minoutid = 0;
	for (i = 1; i < outputc; ++i)
		if (out[i] > out[minoutid])
			minoutid = i;

	// should do something with threshold	
	if (out[minoutid] > 0.5) {
		free(out);
		return minoutid;
	} else {
		free(out);
		return -1;
	}
}

int main(int argc, char **argv)
{
	struct hc_hcascade hc;
	struct nn_neuronet nnet;
	struct nd_image img;
	int *xfound;
	int xfoundcnt;
	int *xjoined;
	int xjoinedcnt;
	int x, i;

	if (hc_hcascaderead(&hc, argv[1]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nn_neuronetfromfile(&nnet, argv[2]) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imgread(argv[3], &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_imgscalebilinear(&img, (double) hc.wh / img.h,
		(double) hc.wh / img.h, &img) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
	
	if (nd_imggrayscale(&img)) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
	
	int rn = 0;

	xfoundcnt = 0;
	xfound = NULL;

	for (x = 0; x < img.w - hc.ww; ++x) {
		struct nd_image imginwin;

		nd_imgcreate(&imginwin, hc.ww, hc.wh, img.format);
		nd_imgcrop(&img, x, 0, hc.ww, hc.wh, &imginwin);
		nd_imgnormalize(&imginwin, 0, 1);
		hc_imgintegral(&imginwin);
		

		if (hc_imgclassify(&hc, &imginwin)) {
			if ((xfound = realloc(xfound,
				sizeof(int) * (++xfoundcnt))) == NULL)
				return 1;

			xfound[xfoundcnt - 1] = x;
		}	
		
		nd_imgdestroy(&imginwin);
	}

	i = 0;

	xjoinedcnt = 0;
	xjoined = NULL;

	do {
		int xsum;
		int j;
		
		j = i;
		xsum = xfound[j];
		while ((xfound[j + 1] - xfound[j] <= 2)
			&& (j + 1 < xfoundcnt)) {
			++j;
			xsum += xfound[j];
		}

		xsum = (double) xsum / (j - i + 1);

		if ((xjoined = realloc(xjoined,
			sizeof(int) * (++xjoinedcnt))) == NULL)
			return 1;

		xjoined[xjoinedcnt - 1] = xsum;

		i = j + 1;

	} while (i < xfoundcnt);

	for (i = 0; i < xjoinedcnt; ++i) {
		struct nd_image imginwin;
		char str[1024];
		double *nninput;
		char a[23] = "0123456789abcehkmoptxy";
		int res;
	
		sprintf(str, "%s/%d.png", argv[3], rn++);

		nd_imgcreate(&imginwin, hc.ww, hc.wh, img.format);
		nd_imgcrop(&img, xjoined[i], 0, hc.ww, hc.wh, &imginwin);
	
		nninput = malloc(sizeof(double) * imginwin.w * imginwin.h);
		
		imgtoinput(&imginwin, nninput);
		
		res = findpattern(&nnet, nninput);

		if (res != -1)
			printf("%c", a[res]);
		
		nd_imgwrite(&imginwin, str);
		nd_imgdestroy(&imginwin);
	}

	printf("\n");

	return 0;
}
