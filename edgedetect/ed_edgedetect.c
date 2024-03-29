#include <math.h>
#include <assert.h>

#include "ed_edgedetect.h"
#include "nd_error.h"
#include "ed_fft.h"

#define HISTSIZE 500

struct houghline {
	double ang;
	double d;
	
	int votes;
};

struct point {
	int x;
	int y;
};

static void ed_safefree(void **p)
{
	if (*p != NULL) {
		free(*p);
		*p = NULL;
	}
}

static int gblurx(double *data, double *tmp, int y, int x, int w, int h,
	double *ws, int sigma3)
{
	double s;	
	int i, j;

	s = 0.0;
	for (i = 0; i <= 2 * sigma3; i++)
		s += ws[i];

	double c = 0.0;
	
	for (j = -sigma3; j <= sigma3; j++) {
		int curx = x + j;
		curx = (curx < 0) ? 0 : curx;
		curx = (curx > (w - 1)) ? (w - 1) : curx;
	
		c += ws[j + sigma3] * data[w * y + curx];
	}
		
	tmp[w * y + x] = c / s;

	return 0;
}

static int gblury(double *data, double *tmp, int y, int x, int w, int h,
	double *ws, int sigma3)
{
	double s;	
	int i, j;

	s = 0.0;
	for (i = 0; i <= 2 * sigma3; i++)
		s += ws[i];

	double c = 0.0;
		
	for (j = -sigma3; j <= sigma3; j++) {
		int cury = y + j;
		cury = (cury < 0) ? 0 : cury;
		cury = (cury > (h - 1)) ? (h - 1) : cury;
	
		c += ws[j + sigma3] * data[w * cury + x];
	}
		
	tmp[w * y + x] = c / s;
	
	return 0;
}

static double gfunc(double x, double sigma)
{
	return exp((-1) * (x * x) / (2.0 * sigma * sigma));
}

int ed_gaussblur(struct nd_image *img, double sigma)
{
	int x, y, i, sigma3;
	int w, h;
	double *ws, *tmpx, *tmpy;
	double *data;
	
	assert(img != NULL && sigma > 0.0);
	assert(nd_imgisvalid(img) && img->format == ND_PF_GRAYSCALE);

	w = img->w;
	h = img->h;
	sigma3 = (int)(3.0 * sigma);
	
	if ((ws = malloc(sizeof(double) * (2 * sigma3 + 1))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	data = (double *)img->data;

	if ((tmpx = (double *)calloc(w * h, sizeof(double))) == NULL) {
		ed_safefree((void **)&tmpx);

		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if ((tmpy = (double *)calloc(w * h, sizeof(double))) == NULL) {
		ed_safefree((void **)&tmpx);
		ed_safefree((void **)&tmpy);
		
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (i = 0; i <= 2*sigma3; i++)
		ws[i] = gfunc(i - sigma3, sigma);

	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			gblurx(data, tmpx, y, x, w, h, ws, sigma3);

	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			gblury(tmpx, tmpy, y, x, w, h, ws, sigma3);

	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			img->data[y * w + x]  = (int)(tmpy[y * w + x]);
	
	free(ws);
	free(tmpx);
	free(tmpy);

	return 0;
}

int ed_removehorfreq(struct nd_image *img, struct nd_image *imgout)
{
	int oldw, oldh;

	oldw = img->w;
	oldh = img->h;

	nd_imgscalebicubic(img, 256 / (double) oldw,
		512 / (double) oldh, imgout);

	struct complexd *datac;

	datac = malloc(sizeof(struct complexd) * imgout->w * imgout->h);

	int pixn;
	for (pixn = 0; pixn < imgout->w * imgout->h; ++pixn) {
		datac[pixn].real = imgout->data[pixn]
			* pow(-1.0, pixn / imgout->w + pixn % imgout->w);
		datac[pixn].imag = 0.0;
	}

	fft2d(datac, imgout->w, imgout->h, 1);

	for (pixn = 0; pixn < imgout->w * imgout->h; ++pixn) {
		int x = pixn % imgout->w;
		if (x >= imgout->w / 2 - 1 && x <= imgout->w + 1)
			datac[pixn].real = datac[pixn].imag = 0.0;
	}

	fft2d(datac, imgout->w, imgout->h, -1);
	
	for (pixn = 0; pixn < imgout->w * imgout->h; ++pixn) {
		imgout->data[pixn] = datac[pixn].real
			* pow(-1.0, pixn / imgout->w + pixn % imgout->w);
	}

	nd_imgscalebicubic(imgout, oldw / (double) 256,
		oldh / (double) 512, imgout);

	return 0;
}

static double imgconvophelp(struct nd_image *img,
	double *op, int opw, int oph, int resposy, int resposx)
{
	int x, y;
	double res;

	res = 0.0;

	for (x = 0; x < oph; ++x)
		for (y = 0; y < opw; ++y) {
			if ( resposy - y  < 0
				|| resposx - x < 0
				|| resposy - y >= img->h
				|| resposx - x >= img->w )
				continue;
						
			res += img->data[(resposy - y)
				* img->w + (resposx - x)] * op[x * opw + y];
	}

	return res;
}

static int imgconvop(struct nd_image *img, double *op, int opw, int oph,
	double *res)
{
	int resposy, resposx;

	for (resposy = 0; resposy < img->h; ++resposy)
		for (resposx = 0; resposx < img->w; ++resposx)
			res[resposy * img->w + resposx]
				= imgconvophelp(img, op, opw, oph,
					resposy + 1, resposx + 1); 	

	return 0;
}

static int sobel(struct nd_image *img, double *gradval, double *graddir)
{
	double sobel3y[9] = { -1.0, -2.0, -1.0,
			0.0, 0.0, 0.0,
			1.0, 2.0, 1.0 };
	double sobel3x[9] = { -1.0, 0.0, 1.0,
			-2.0, 0.0, 2.0,
			-1.0, 0.0, 1.0 };

	double *gy, *gx;
	int imgx, imgy;

	if ((gy = malloc(sizeof(double) * img->w * img->h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if ((gx = malloc(sizeof(double) * img->w * img->h)) == NULL) {
		ed_safefree((void **)&gx);
	
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	imgconvop(img, sobel3y, 3, 3, gy);
	imgconvop(img, sobel3x, 3, 3, gx);

	for (imgy = 0; imgy < img->h; ++imgy)
		for (imgx = 0; imgx < img->w; ++imgx) {
			gradval[imgy * img->w + imgx]
				= sqrt(pow(gx[imgy * img->w + imgx], 2.0)
				+ pow(gy[imgy * img->w + imgx], 2.0));	

			graddir[imgy * img->w + imgx]
				= atan2(gy[imgy * img->w + imgx],
				 gx[imgy * img->w + imgx]);
		}

	free(gx);
	free(gy);
	
	return 0;
}

static double imgmaxval(double *pix, int imgsize)
{
	double maxval;
	int pixn;

	maxval = pix[0];
	for (pixn = 1; pixn < imgsize; ++pixn) {
		double val;

		val = pix[pixn];
		maxval = (val > maxval) ? val : maxval;
	}

	return maxval;
}

static int buildhist(int **hist, int histsize,
	double *pix, int imgsize, double maxval)
{
	double rel;
	int pixn;

	if (((*hist) = calloc(histsize, sizeof(int))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	rel = (double) (histsize - 1) / maxval;
	
	for (pixn = 1; pixn < imgsize; ++pixn) {
		int histel;

		histel = (int) ceil(pix[pixn] * rel);

// !!!
		if (histel < 0 || histel >= histsize)
			continue;
//
	
		++((*hist)[histel]);
	}

	return 0;
}

static int histwrite(int *hist, int histsize, double histscale,
	int imghisth, const char *path)
{
	int x, y;
	
	struct nd_image imghist;
	
	if (nd_imgcreate(&imghist, histsize, imghisth, ND_PF_GRAYSCALE) < 0)
		return (-1);

	for (y = 0; y < imghisth; ++y)
		for (x = 0; x < histsize; ++x) 
			imghist.data[y * histsize + x]
				= ((imghisth - y) < hist[x] * histscale)
				? 1.0 : 0.0;

	if (nd_imgwrite(&imghist, path) < 0)
		return (-1);

	return 0;
}

static int otsu(double *pix, int imgsize, int histsize, double *thres)
{
	double maxval;
	int *hist;

	int m, n;
	int a1, b1;	
	int t;

	double maxsigma;
	int thr;

	maxval = imgmaxval(pix, imgsize);

	if (buildhist(&hist, histsize, pix, imgsize, maxval) < 0)
		return (-1);

	m = n = 0;

	for (t = 0; t <= histsize; ++t) {
		m += t * hist[t];
	}
	
	n = imgsize;
	
	
	maxsigma = -1.0;
	thr = 0;
	a1 = b1 = 0;
	
	for (t = 0; t < histsize - 1; ++t) {
		double w1;
		double a;
		double sigma;
		
		a1 += t * hist[t];
		b1 += hist[t];

		w1 = (double) b1 / n;

		a = (double) a1 /  b1 - (m - a1) / (n - b1);
	
		sigma = w1 * (1 - w1) * a * a;

		if (sigma > maxsigma) {
			maxsigma = sigma;
			thr = t;
		}
	}

	*thres = maxval * thr / histsize;

	free(hist);

	return 0;
}

static int ed_edgenbours(int *res, int resw, int imgy, int imgx)
{
	return (res[(imgy - 1) * resw + imgx]
		&& res[(imgy + 1) * resw + imgx]
		&& res[(imgy - 1) * resw + imgx]
		&& res[(imgy - 1) * resw + (imgx - 1)]
		&& res[(imgy - 1) * resw + (imgx + 1)]
		&& res[(imgy + 1) * resw + (imgx + 1)]
		&& res[(imgy + 1) * resw + (imgx + 1)]
		&& res[imgy * resw + (imgx - 1)]
		&& res[imgy * resw + (imgx + 1)]);
}

static int removehorgrad(double *gradval, double *graddir, int pixcount)
{
	int pixn;

	for (pixn = 0; pixn < pixcount; ++pixn) {
		int sectind;
		
		sectind = 4.0 * (graddir[pixn]) / M_PI + 1.0 / 2.0;
		
		if (sectind == 4 || sectind == 0 )
			gradval[pixn] = 0.0;
	
		if (sectind == 2 || sectind == 6)
			gradval[pixn] *= 2.0;
	}

	return 0;
}

int ed_canny(struct nd_image *img, int *outmask, double thres1, double thres2,
	int onlyver)
{
	double *gradval, *graddir;
	int imgx, imgy;

	assert(img != NULL && outmask != NULL);
	assert(nd_imgisvalid(img) && img->format == ND_PF_GRAYSCALE);

	if ((gradval = malloc(sizeof(double) * img->w * img->h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if ((graddir = malloc(sizeof(double) * img->w * img->h)) == NULL) {
		ed_safefree((void **)&gradval);
		
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if (sobel(img, gradval, graddir) < 0) {
		ed_safefree((void **)&gradval);
		ed_safefree((void **)&graddir);
	
		return (-1);
	}

	if (onlyver)
		removehorgrad(gradval, graddir, img->w * img->h);	
	
	for (imgy = 0; imgy < img->h; ++imgy)
		for (imgx = 0; imgx < img->w; ++imgx)
			outmask[imgy * img->w + imgx] = 0;

	for (imgy = 1; imgy < img->h - 1; ++imgy)
		for (imgx = 1; imgx < img->w - 1; ++imgx) {
			int xleft, yleft;
			int xright, yright;

			if (graddir[imgy * img->w + imgx] < 0.0) 
				graddir[imgy * img->w + imgx] += 2.0 * M_PI;

			int sectind =
				4.0 * (graddir[imgy * img->w + imgx]) / M_PI
				+ 1.0 / 2.0;

			sectind %= 8;;

			xleft = imgx;
			yleft = imgy;
			xright = imgx;
			yright = imgy;

			switch (sectind) {
			case 0:
			case 4:
				yleft += 1;
				yright -= 1;
				break;
			case 2:
			case 6:
				xleft -= 1;
				xright += 1;
				break;
			case 1:
			case 5:
				xleft -= 1;
				yleft -= 1;
				xright += 1;
				yright += 1;
				break;	
			case 3:
			case 7:
				xleft -= 1;
				yleft += 1;
				xright += 1;
				yright -= 1;
				break;
			}

			if ((gradval[imgy * img->w + imgx]
				< gradval[yleft * img->w + xleft])
				|| (gradval[imgy * img->w + imgx]
				< gradval[yright * img->w + xright])) {

				outmask[imgy * img->w + imgx] = 0;
				gradval[imgy * img->w + imgx] = 0.0;
			}
			else
				outmask[imgy * img->w + imgx] = 2;	
		}

	if (thres2 < 0.0 || thres1 < 0.0) {
		if (otsu(gradval, img->w * img->h, HISTSIZE, &thres2) < 0)
			return (-1);
		
		thres1 = thres2 * 0.5;
	}
	
	for (imgy = 1; imgy < img->h - 1; ++imgy)
		for (imgx = 1; imgx < img->w - 1; ++imgx)
			if (gradval[imgy * img->w + imgx] < thres2) {
				if (gradval[imgy * img->w + imgx]
					> thres1)
					outmask[imgy * img->w + imgx] = 1;
				else
					outmask[imgy * img->w + imgx] = 0;
			}

// should create new buffer
	for (imgy = 1; imgy < img->h - 1; ++imgy)
		for (imgx = 1; imgx < img->w - 1; ++imgx)
			if (outmask[imgy * img->w + imgx] == 1) {
				if (!ed_edgenbours(outmask,
					img->w, imgy, imgx))
					outmask[imgy * img->w + imgx] = 0;
			}

	return 0;
}

static int compline(const void *hl0, const void *hl1)
{
	return (((struct houghline *) hl1)->votes
		- ((struct houghline *) hl0)->votes);
}

static int ed_houghismax(int *acc, int arange, int drange, int ang, int d)
{
	if (ang - 1 > 0 && acc[d * arange + ang]
		<= acc[d * arange + ang - 1])
		return 0;

	if (ang + 1 < arange && acc[d * arange + ang]
		< acc[d * arange + ang + 1])
		return 0;

	if (d - 1 > 0 && acc[d * arange + ang]
		<= acc[(d - 1) * arange + ang])
		return 0;

	if (d + 1 < drange && acc[d * arange + ang]
		< acc[(d + 1) * arange + ang])
		return 0;

	return 1;
}

int ed_hough(int *img, int imgw, int imgh, double dang,
	double *lines, int linemaxc)
{
	int arange, drange;
	int *acc;
	int imgx, imgy;
	int ang, d;

	struct houghline *hl;
	int linec;
	int linen;

	assert(img != NULL && imgw > 0 && imgh > 0
		&& dang > 0.0 && lines != NULL && linemaxc > 0);

	arange = M_PI / dang;
	drange = (imgw + imgh) * 2 + 1;

	if ((acc = calloc(arange * drange, sizeof(int))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (imgy = 0; imgy < imgh; ++imgy)
		for (imgx = 0; imgx < imgw; ++imgx)
			if (img[imgy * imgw + imgx])
				for (ang = 0; ang < arange; ++ang) {
					double rang = (double) ang * dang;

					d = round(imgx * cos(rang)
						+ imgy * sin(rang));
				
					d += drange / 2;

					++acc[d * arange + ang];
				}

	linec = 0;
	for (d = 0; d < drange; ++d)
		for (ang = 0; ang < arange; ++ang)
			if (acc[d * arange + ang]
				&& ed_houghismax(acc, arange, drange, ang, d))
				++linec;

	if ((hl = malloc(sizeof(struct houghline) * linec)) == NULL) {
		ed_safefree((void **)&acc);
	
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	linen = 0;
	for (d = 0; d < drange; ++d)
		for (ang = 0; ang < arange; ++ang)
			if (acc[d * arange + ang]
				&& ed_houghismax(acc, arange,
				drange, ang, d)) {
				hl[linen].ang = (double) ang * dang;
				hl[linen].d = (double) (d - drange / 2);
				hl[linen].votes = acc[d * arange + ang];

				++linen;
			}

	qsort(hl, linec, sizeof(struct houghline), compline);

	linen = 0;
	for (linen = 0; linen < linemaxc; ++linen) {
		lines[2 * linen + 0] = hl[linen].ang;
		lines[2 * linen + 1] = hl[linen].d;
	}

	ed_safefree((void **)&acc);
	ed_safefree((void **)&hl);

	return 0;
}
