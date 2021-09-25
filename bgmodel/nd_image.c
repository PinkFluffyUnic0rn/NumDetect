#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <cairo.h>

#include "nd_image.h"
#include "nd_error.h"
#include "nd_vecmat.h"

int nd_imgisvalid(const struct nd_image *img)
{
	if (img->data == NULL || img->w <= 0 || img->h <= 0 || img->format < 0)
		return 0;
	else
		return 1;
}

int nd_imgchanscount(enum ND_PIXELFORMAT format)
{
	switch (format) {
	case ND_PF_GRAYSCALE:
		return 1;
	case ND_PF_RGB:
		return 3;
	case ND_PF_ARGB:
		return 3;
	default:
		return (-1);
	}
}

int nd_imgcreate(struct nd_image *img, int w, int h,
	enum ND_PIXELFORMAT format)
{
	img->w = w;
	img->h = h;
	img->format = format;

	assert(w > 0 && h > 0 && format >= 0);
	
	if ((img->data = malloc(sizeof(double) * w * h
		* nd_imgchanscount(format))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	return 0;
}

int nd_imgcopy(const struct nd_image *imgsrc, struct nd_image *imgdest)
{
	imgdest->w = imgsrc->w;
	imgdest->h = imgsrc->h;
	imgdest->format = imgsrc->format;

	assert(nd_imgisvalid(imgsrc));

	if ((imgdest->data = malloc(sizeof(double) * imgdest->w * imgdest->h
		* nd_imgchanscount(imgdest->format))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	memcpy(imgdest->data, imgsrc->data, sizeof(double)
		* imgsrc->w * imgsrc->h * nd_imgchanscount(imgsrc->format));

	return 0;
}

int nd_imgdestroy(struct nd_image *img)
{
	img->w = img->h = img->format = 0;

	free(img->data);

	return 0;
}

int nd_imgread(const char *imgpath, struct nd_image *img)
{
	cairo_surface_t *sur;	
	unsigned char *data;
	int pixsz;
	int imgx, imgy, nchan;
	
	assert(img != NULL && imgpath != NULL);

	sur = cairo_image_surface_create_from_png(imgpath);
	
	if (cairo_surface_status(sur) != CAIRO_STATUS_SUCCESS) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		return (-1);
	}
	
	data = cairo_image_surface_get_data(sur);
	img->w = cairo_image_surface_get_width(sur);
	img->h = cairo_image_surface_get_height(sur);

	switch (cairo_image_surface_get_format(sur)) {
	case CAIRO_FORMAT_RGB24:
	case CAIRO_FORMAT_ARGB32:
		pixsz = 4;
		img->format = ND_PF_RGB;
		break;
	
	default:
		nd_seterrormessage("image has a wrong format", __func__);
		return (-1);
	}

	if ((img->data = malloc(sizeof(double *)
		* img->w * img->h * nd_imgchanscount(img->format)))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}
	
	for (imgy = 0; imgy < img->h; ++imgy)
		for (imgx = 0; imgx < img->w; ++imgx)
			for (nchan = 0;
				nchan < nd_imgchanscount(img->format);
				++nchan) {
				int surdatapos, imgdatapos;
				
				surdatapos = (imgy * img->w + imgx)
					* pixsz + nchan;
				
				imgdatapos = (imgy * img->w + imgx)
					* nd_imgchanscount(img->format)
					+ nchan;	

				img->data[imgdatapos]
					= data[surdatapos] / 255.0;
			}

	cairo_surface_destroy(sur);
	
	return 0;
}

int nd_imgwrite(const struct nd_image *img, const char *imgpath)
{
	cairo_surface_t *sur;	
	unsigned char *surdata;
	int imgy, imgx, nchan;

	assert(img != NULL && imgpath != NULL);
	assert(nd_imgisvalid(img));

	sur = cairo_image_surface_create(CAIRO_FORMAT_RGB24, img->w, img->h);	
	
	if (cairo_surface_status(sur) != CAIRO_STATUS_SUCCESS) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		return (-1);
	}
	
	surdata = cairo_image_surface_get_data(sur);
	
	switch(img->format) {
	case ND_PF_GRAYSCALE:
		for (imgy = 0; imgy < img->h; ++imgy)
			for (imgx = 0; imgx < img->w; ++imgx) {
				int val;
			
				val = ceil(img->data[(imgy * img->w + imgx)]
					* 255.0);
				val = (val <= 255)
					? (val >= 0 ? val : 0) : 255;
					
				surdata[(imgy * img->w + imgx) * 4 + 0] = val;
				surdata[(imgy * img->w + imgx) * 4 + 1] = val;
				surdata[(imgy * img->w + imgx) * 4 + 2] = val;
			}
		break;

	case ND_PF_RGB:
		for (imgy = 0; imgy < img->h; ++imgy)
			for (imgx = 0; imgx < img->w; ++imgx)
				for (nchan = 0;
					nchan < nd_imgchanscount(img->format);
					++nchan) {
					int val;
					size_t imgdatapos;	
		
					imgdatapos = (imgy * img->w + imgx)
						* nd_imgchanscount(img->format)
						+ nchan;

					val = ceil(img->data[imgdatapos]
						* 255.0);
					val = (val <= 255)
						? (val >= 0 ? val : 0) : 255;
					
					surdata[(imgy * img->w + imgx) * 4
						+ nchan] = val;
			}
		break;
	
	default:
		nd_seterrormessage("image has a wrong format", __func__);
		return (-1);
	}

	if (cairo_surface_write_to_png(sur, imgpath) != CAIRO_STATUS_SUCCESS) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		return (-1);
	}
	
	cairo_surface_destroy(sur);

	return 0;
}

int nd_imghsvval(struct nd_image *img)
{
	double *hsvval;
	int npix;
	
	assert(img != NULL);
	assert(nd_imgisvalid(img));

	if ((hsvval = malloc(sizeof(double) * img->w * img->h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (npix = 0; npix < img->w * img->h; ++npix) {
		int nchan;
		
		hsvval[npix]
			= img->data[npix * nd_imgchanscount(img->format) + 0];

		for (nchan = 1;
			nchan < nd_imgchanscount(img->format); ++nchan) {
			double chanval;
			
			chanval = img->data[npix
				* nd_imgchanscount(img->format) + nchan];
			hsvval[npix] = (chanval > hsvval[npix])
				? chanval : hsvval[npix];
		}
	}

	free(img->data);
	
	img->format = ND_PF_GRAYSCALE;
	img->data = hsvval;

	return 0;
}

int nd_imggrayscale(struct nd_image *img)
{
	double *grayscale;
	int npix;
	
	assert(img != NULL);
	assert(nd_imgisvalid(img));

	if ((grayscale = malloc(sizeof(double) * img->w * img->h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (npix = 0; npix < img->w * img->h; ++npix) {
		int nchan;

		grayscale[npix] = 0.0;
	
		for (nchan = 0; nchan < nd_imgchanscount(img->format); ++nchan)
			grayscale[npix]
				+= img->data[npix
					* nd_imgchanscount(img->format)
					+ nchan];

		grayscale[npix] /= nd_imgchanscount(img->format);
	}

	free(img->data);
	
	img->format = ND_PF_GRAYSCALE;
	img->data = grayscale;

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

/*
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
*/

int nd_histequalization(struct nd_image *img, int histsize)
{
	int *hist;
	double maxval;
	double *newval;
	int i;

	maxval = imgmaxval(img->data, img->w * img->h);

	if (buildhist(&hist, histsize, img->data, img->w * img->h, maxval) < 0)
		return (-1);

	if ((newval = malloc(sizeof(double) * histsize)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	newval[0] = hist[0] / (double) (img->w * img->h);
	for (i = 1; i < histsize; ++i)
		newval[i] = newval[i - 1] + hist[i]
		/ (double) (img->w * img->h);

	for (i = 0; i < img->h * img->w; ++i) {
		int histel;

		histel = (int) ceil(img->data[i] / maxval * (histsize - 1));

		if (histel < 0 || histel >= histsize)
			continue;

		img->data[i] = newval[histel];
	}

	return 0;
}




int nd_imgnormalize(struct nd_image *img, int normavr, int normdev)
{
	int npix;
	double pixsum;
	int pixcount;
	double sd;

	assert(img != NULL);
	assert(nd_imgisvalid(img) && img->format == ND_PF_GRAYSCALE);

	pixcount = img->w * img->h;
	sd = 1.0;

	pixsum = 0.0;
	for (npix = 0; npix < pixcount; ++npix)
		pixsum += img->data[npix];

	if (normdev) {
		double sqpixsum;
	
		sqpixsum = 0.0;
		for (npix = 0; npix < pixcount; ++npix)
			sqpixsum += img->data[npix] * img->data[npix];
	
		sd = sqrt(1.0 / (pixcount - 1.0)
			* (sqpixsum - pixsum * pixsum / pixcount));
	}

	if (normavr) {
		double avr;
		
		avr = pixsum / pixcount;
		
		for (npix = 0; npix < pixcount; ++npix)
			img->data[npix] -= avr;
	}

	if (normdev) {
		if (fabs(sd) < 0.00001) {
			nd_seterrormessage(
				"Standard deviation is close to zero",
				__func__);
			return (-1);
		}

		for (npix = 0; npix < pixcount; ++npix)
			img->data[npix] /= sd;
	}

	return 0;
}

int nd_imgcrop(const struct nd_image *imgin, int x0, int y0, int w, int h,
	struct nd_image *imgout)
{
	int y;

	assert(imgin != NULL
		&& x0 >= 0 && x0 < imgin->w 
		&& y0 >= 0 && y0 < imgin->h
		&& w > 0 && h > 0
		&& x0 + w <= imgin->w && y0 + h <= imgin->h
		&& imgout != NULL
		&& imgout->w >= w
		&& imgout->h >= h
		&& imgout->format == imgin->format);
	assert(nd_imgisvalid(imgin));
	assert(nd_imgisvalid(imgout));

	for (y = 0; y < h; ++y)
		memcpy(imgout->data
			+ y * imgout->w * nd_imgchanscount(imgout->format),
			imgin->data
			+ (y + y0) * imgin->w * nd_imgchanscount(imgin->format)
			+ x0 * nd_imgchanscount(imgin->format),
			sizeof(double) * w * nd_imgchanscount(imgin->format));

	return 0;
}

static double nd_cubicinterpx(const struct nd_image *inimg,
	double x, int y, int cn)
{
	double a, b, c, d;
	double p0, p1, p2, p3;
	double xr;

	p1 = inimg->data[(y * inimg->w + (int)floor(x))
		* nd_imgchanscount(inimg->format) + cn];
	p0 = (((int)(floor(x) - 1)) >= 0)
		? inimg->data[(y * inimg->w + (int)(floor(x) - 1))
		* nd_imgchanscount(inimg->format) + cn] : p1;
	p2 = (((int)ceil(x)) < (inimg->w))
		? inimg->data[(y * inimg->w + (int)ceil(x))
		* nd_imgchanscount(inimg->format) + cn] : p1;
	p3 = ((int)(ceil(x) + 1) < (inimg->w))
		? inimg->data[(y * inimg->w + (int)(ceil(x) + 1))
		* nd_imgchanscount(inimg->format) + cn] : p2;

	a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
	b = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
	c = -0.5 * p0 + 0.5 * p2;
	d = p1;

	xr = x - floor(x);

	return a * xr * xr * xr + b * xr * xr + c * xr + d;
}

static double nd_cubicinterp(const struct nd_image *inimg,
	double x, double y, int cn)
{
	double a, b, c, d;
	double p0, p1, p2, p3;
	double yr;
	
	p1 = nd_cubicinterpx(inimg, x, (int)floor(y), cn);
	p0 = (((int)(floor(y) - 1)) >= 0)
		? nd_cubicinterpx(inimg, x, (int)(floor(y) - 1), cn) : p1;
	p2 = ((int)ceil(y)) < (inimg->h)
		? nd_cubicinterpx(inimg, x, (int)ceil(y), cn) : p1;
	p3 = ((int)(ceil(y) + 1) < (inimg->h))
		? nd_cubicinterpx(inimg, x, (int)(ceil(y) + 1), cn) : p2;

	a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
	b = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
	c = -0.5 * p0 + 0.5 * p2;
	d = p1;

	yr = y - floor(y);
	return a * yr * yr * yr + b * yr * yr + c * yr + d;
}

int nd_imgscalebicubic(const struct nd_image *inimg, double wrel, double hrel,
	struct nd_image *outimg)
{
	int x, y;
	struct nd_image tmpimg;

	assert(inimg != NULL && hrel > 0.0 && wrel > 0.0 && outimg != NULL);
	assert(nd_imgisvalid(inimg));

	tmpimg.w = ceil(wrel * inimg->w);	
	tmpimg.h = ceil(hrel * inimg->h);
	tmpimg.format = inimg->format;

	if ((tmpimg.data = malloc(sizeof(double) * tmpimg.w * tmpimg.h
		* nd_imgchanscount(tmpimg.format))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (y = 0; y < tmpimg.h; ++y)
		for (x = 0; x < tmpimg.w; ++x) {
			double iny, inx;
			int c;

			iny = ((double) y) / hrel;
			inx = ((double) x) / wrel;
			
			
			for (c = 0; c < nd_imgchanscount(inimg->format); ++c) {
				tmpimg.data[(y * tmpimg.w + x)
					* nd_imgchanscount(tmpimg.format) + c]
					= nd_cubicinterp(inimg, inx, iny, c);
			}
		}

	if (inimg == outimg)
		free(outimg->data);

	outimg->w = tmpimg.w;
	outimg->h = tmpimg.h;
	outimg->format = tmpimg.format;
	outimg->data = tmpimg.data;

	return 0;
}

static double nd_linearinterpx(const struct nd_image *inimg,
	double x, int y, int c)
{
	double a, b;
	double p0, p1;
	double xr;

	p0 = inimg->data[(y * inimg->w + (int)floor(x))
		* nd_imgchanscount(inimg->format) + c];
	p1 = (((int)ceil(x)) < (inimg->w))
		? inimg->data[(y * inimg->w + (int)ceil(x))
		* nd_imgchanscount(inimg->format) + c] : p0;

	a = p1 - p0;
	b = p0;

	xr = x - floor(x);

	return a * xr + b;
}

static double nd_linearinterp(const struct nd_image *inimg,
	double x, double y, int c)
{
	double a, b;
	double p0, p1;
	double yr;
	
	p0 = nd_linearinterpx(inimg, x, (int)floor(y), c);
	p1 = ((int)ceil(y) < (inimg->h))
		? nd_linearinterpx(inimg, x, (int)ceil(y), c) : p0;

	a = p1 - p0;
	b = p0;

	yr = y - floor(y);
	return a * yr + b;
}

int nd_imgscalebilinear(const struct nd_image *inimg, double wrel, double hrel,
	struct nd_image *outimg)
{
	int x, y;
	struct nd_image tmpimg;

	assert(inimg != NULL && hrel > 0.0 && wrel > 0.0 && outimg != NULL);
	assert(nd_imgisvalid(inimg));

	tmpimg.w = ceil(wrel * inimg->w);	
	tmpimg.h = ceil(hrel * inimg->h);
	tmpimg.format = inimg->format;

	if ((tmpimg.data = malloc(sizeof(double) * tmpimg.w * tmpimg.h
		* nd_imgchanscount(tmpimg.format))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}
	
	for (y = 0; y < tmpimg.h; ++y)
		for (x = 0; x < tmpimg.w; ++x) {
			double iny, inx;
			int c;
	
			iny = ((double) y) / hrel;
			inx = ((double) x) / wrel;
				
			for (c = 0; c < nd_imgchanscount(inimg->format); ++c) {
				tmpimg.data[(y * tmpimg.w + x)
					* nd_imgchanscount(tmpimg.format) + c]
					= nd_linearinterp(inimg, inx, iny, c);
			}
		}

	if (inimg == outimg)
		free(outimg->data);

	outimg->w = tmpimg.w;
	outimg->h = tmpimg.h;
	outimg->format = tmpimg.format;
	outimg->data = tmpimg.data;

	return 0;
}

int nd_getpersptransform(double *inpoints, double *outpoints,
	struct nd_matrix3 *mr)
{
	struct nd_matrix3 a;
	struct nd_matrix3 b;
	struct nd_vector3 tmpv;
	double r[3];

	assert(inpoints != NULL && outpoints != NULL && mr != NULL);

	a._11 = inpoints[0 * 2 + 0];	a._12 = inpoints[1 * 2 + 0];
		a._13 = inpoints[2 * 2 + 0];
	a._21 = inpoints[0 * 2 + 1];	a._22 = inpoints[1 * 2 + 1];
		a._23 = inpoints[2 * 2 + 1];
	a._31 = 1.0;			a._32 = 1.0;
		a._33 = 1.0;

	tmpv.x = inpoints[3 * 2 + 0];
	tmpv.y = inpoints[3 * 2 + 1];
	tmpv.z = 1.0;

	if (nd_m3nonhomsolve(&a, &tmpv, r) < 0)
		return (-1);

	a._11 *= r[0]; a._12 *= r[1]; a._13 *= r[2]; 
	a._21 *= r[0]; a._22 *= r[1]; a._23 *= r[2]; 
	a._31 *= r[0]; a._32 *= r[1]; a._33 *= r[2]; 

	b._11 = outpoints[0 * 2 + 0];	b._12 = outpoints[1 * 2 + 0];
		b._13 = outpoints[2 * 2 + 0];
	b._21 = outpoints[0 * 2 + 1];	b._22 = outpoints[1 * 2 + 1];
		b._23 = outpoints[2 * 2 + 1];
	b._31 = 1.0;			b._32 = 1.0;
		b._33 = 1.0;

	tmpv.x = outpoints[3 * 2 + 0];
	tmpv.y = outpoints[3 * 2 + 1];
	tmpv.z = 1.0;

	if (nd_m3nonhomsolve(&b, &tmpv, r) < 0)
		return (-1);
	
	b._11 *= r[0]; b._12 *= r[1]; b._13 *= r[2]; 
	b._21 *= r[0]; b._22 *= r[1]; b._23 *= r[2]; 
	b._31 *= r[0]; b._32 *= r[1]; b._33 *= r[2]; 

	if (nd_m3inverse(&b, &b) < 0)
		return (-1);

	nd_m3mult(&a, &b, mr);

	return 0;
}

int nd_imgapplytransform(struct nd_image *imgin, const struct nd_matrix3 *m,
	struct nd_image *imgout)
{
	int x, y;
	double *newdata;

	assert(imgin != NULL && m != NULL && imgout);

	assert(nd_imgisvalid(imgin)
		&& (imgin->format == ND_PF_GRAYSCALE
		|| imgin->format == ND_PF_RGB));

	if ((newdata = malloc(sizeof(double) * imgin->w * imgin->h
		* nd_imgchanscount(imgin->format))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	for (y = 0; y < imgin->h; ++y)
		for (x = 0; x < imgin->w; ++x) {
			struct nd_vector3 v;
			double inx, iny;
			int c;

			v.x = (double) x;
			v.y = (double) y;
			v.z = 1.0;
			
			nd_v3m3mult(&v, m, &v);
		
			inx = v.x / v.z;
			iny = v.y / v.z;

			for (c = 0; c < nd_imgchanscount(imgin->format); ++c) {
				newdata[(y * imgin->w + x)
					* nd_imgchanscount(imgin->format) + c]
					= (iny >= 0 && iny < imgin->h
					&& inx >= 0 && inx < imgin->w)
					? nd_linearinterp(imgin, inx, iny, c)
					: 0.0;

			}
		}
	if (imgin != imgout) {
		memmove(imgout, imgin, sizeof(struct nd_image));
		imgout->data = newdata;
	}
	else {
		memcpy(imgout->data, newdata, sizeof(double) * imgin->w
			* imgin->h * nd_imgchanscount(imgin->format));
		
		free(newdata);
	}

	return 0;
}

static double nd_pixrgbtogray(double *pix)
{
	return 0.299 * pix[0] + 0.587 * pix[1] + 0.114 * pix[2];
}

static void nd_sortgraypix(double *pix, int *idx, int pixc)
{
	int i, j;

	for (i = 0; i < pixc; ++i) {
		int minp;
		double tmps;
		int tmpi;
		
		minp = i;
		for (j = i + 1; j < pixc; ++j)
			minp = (pix[j] < pix[minp]) ? j : minp;

		tmps = pix[i];
		pix[i] = pix[minp];
		pix[minp] = tmps;

		tmpi = idx[i];
		pix[i] = idx[minp];
		idx[minp] = tmpi;
	}
}

int nd_imgmedianfilter3(struct nd_image *img)
{
	int x, y;
	double *newdata;

	if ((newdata = malloc(sizeof(double) * img->w * img->h
		* nd_imgchanscount(img->format))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	assert(nd_imgisvalid(img) && img->format == ND_PF_RGB);

	for (y = 1; y < img->h - 1; ++y)
		for (x = 1; x < img->w - 1; ++x) {
			double grayscale[9];
			int idx[9];
			int i;

			idx[0] = (y - 1) * img->w + (x - 1);
			grayscale[0] = nd_pixrgbtogray(img->data + idx[0] * 3);
			idx[1] = (y - 1) * img->w + x;
			grayscale[1] = nd_pixrgbtogray(img->data + idx[1] * 3);
			idx[2] = (y - 1) * img->w + (x + 1);
			grayscale[2] = nd_pixrgbtogray(img->data + idx[2] * 3);
			
			idx[3] = y * img->w + (x - 1);
			grayscale[3] = nd_pixrgbtogray(img->data + idx[3] * 3);
			idx[4] = y * img->w + x;
			grayscale[4] = nd_pixrgbtogray(img->data + idx[4] * 3);
			idx[5] = y * img->w + (x + 1);
			grayscale[5] = nd_pixrgbtogray(img->data + idx[5] * 3);
				
			idx[6] = (y + 1) * img->w + (x - 1);
			grayscale[6] = nd_pixrgbtogray(img->data + idx[6] * 3);
			idx[7] = (y + 1) * img->w + x;
			grayscale[7] = nd_pixrgbtogray(img->data + idx[7] * 3);
			idx[8] = (y + 1) * img->w + (x + 1);
			grayscale[8] = nd_pixrgbtogray(img->data + idx[8] * 3);
		
			nd_sortgraypix(grayscale, idx, 9);

			for (i = 0; i < 9; ++i)
				if (i != 4)
					memcpy(newdata + idx[i] * 3,
						img->data + idx[4] * 3,
						sizeof(double) * 3);
		}

	memcpy(img->data, newdata, sizeof(double) * img->w * img->h
		* nd_imgchanscount(img->format));
	free(newdata);

	return 0;
}
