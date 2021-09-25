#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "nd_image.h"
#include "nd_error.h"
#include "ed_edgedetect.h"

#include <gtk/gtk.h>
#include <cairo.h>

#define MAXLINES 10
#define MAXANGDIF 2.5 * M_PI / 180.0
#define VBANDLEN 0.125
#define THETAMINVAL 0.001

enum LINEORIENTATION {
	LO_TOP,
	LO_BOTTOM,
	LO_LEFT,
	LO_RIGHT
};

static void safefree(void **p)
{
	if (*p != NULL) {
		free(*p);
		*p = NULL;
	}
}

static int imghsplit(struct nd_image *img,
	struct nd_image *imgparts, int partc)
{
	int partn;

	for (partn = 0; partn < partc; ++partn) {
		if (nd_imgcreate(imgparts + partn,
			img->w, img->h / partc, ND_PF_GRAYSCALE) < 0) {
			int partnn;

			for (partnn = 0; partnn < partn; ++partnn)
				nd_imgdestroy(imgparts + partn);
		
			return (-1);
		}
		
		memcpy(imgparts[partn].data,
			img->data + (img->w * (img->h * partn / partc)),
			sizeof(double) * img->w * (img->h / partc));
	}

	return 0;
}

static int imgvsplit(struct nd_image *img,
	struct nd_image *imgparts, int partc)
{
	int partw;
	int partn;

	partw = img->w / partc;

	for (partn = 0; partn < partc; ++partn) {
		int imgy;	
		
		if (nd_imgcreate(imgparts + partn,
			partn != partc - 1 ? partw : img->w - partw * partn,
			img->h, ND_PF_GRAYSCALE) < 0) {
			int partnn;

			for (partnn = 0; partnn < partn; ++ partnn)
				nd_imgdestroy(imgparts + partn);
		
			return (-1);
		}
		
		for (imgy = 0; imgy < img->h; ++imgy)
			memcpy(imgparts[partn].data + imgy * imgparts[partn].w,
				img->data + (imgy * img->w) + partw * partn,
				sizeof(double) * imgparts[partn].w);
	}

	return 0;
}

static double bandsum(struct nd_image *img,
	double rho0, double rho1, double theta)
{	
	int x, y;
	double s;

	s = 0.0;

	for (y = 0; y < img->h; ++y)
		for (x = 0; x < img->w; ++x) {
			double pabs;
			double parg;
			double prho;

			parg = atan2((double) y, (double) x);
			pabs = sqrt(pow(x, 2.0) + pow(y, 2.0));

			prho = pabs * cos(parg - theta);
		
			if (prho > rho0 && prho < rho1)
				s += img->data[y * img->w + x];
		}

	return s;
}

static int checkline(struct nd_image *img, double theta, double rho,
	enum LINEORIENTATION lo)
{	
	if ((lo == LO_TOP || lo == LO_BOTTOM)
		&& (theta > M_PI * 0.25 && theta < M_PI * 0.75))
		return 1;	

	if ((lo == LO_LEFT || lo == LO_RIGHT)
		&& !(theta > M_PI * 0.25 && theta < M_PI * 0.75)) {	
		double sumborder, suml, sumr;
		double ivbandlen;
		double valdif;

//		return 1;

		ivbandlen = img->w * VBANDLEN;
		
		suml = bandsum(img, rho - ivbandlen * 1.5 - 0.01,
			rho - ivbandlen * 0.5 + 0.01, theta);
		
		sumr = bandsum(img, rho + ivbandlen * 0.5 - 0.01,
			rho + ivbandlen * 1.5 + 0.01, theta);
		
		sumborder = bandsum(img, rho - ivbandlen * 0.5 - 0.01,
			rho + ivbandlen * 0.5 + 0.01, theta);

		valdif = ((lo == LO_LEFT) ? sumr : suml) - sumborder;

		if (theta > M_PI / 2.0)
			valdif *= -1.0;

		return valdif > 0.0;
	}

	return 0;
}

static int imgdrawline(struct nd_image *img, struct lineseg *l)
{
	int deltax, deltay;
	int error;
	int signy;
	int x0, y0, x1, y1;

	assert(img->format == ND_PF_RGB);

	x0 = l->x0;
	y0 = l->y0;
	x1 = l->x1;
	y1 = l->y1;

	if (x0 > x1) {
		double tmp;

		tmp = x0;
		x0 = x1;
		x1 = tmp;

		tmp = y0;
		y0 = y1;
		y1 = tmp;
	}

	deltax = abs(x1 - x0);
	deltay = abs(y1 - y0);
	signy = y0 < y1 ? 1 : -1;
	error = deltax - deltay;

	while (x0 != x1 || y0 != y1) {
		int error2;

		img->data[(y0 * img->w + x0) * 3 + 0] = 0.0;
		img->data[(y0 * img->w + x0) * 3 + 1] = 1.0;
		img->data[(y0 * img->w + x0) * 3 + 2] = 0.0;

		error2 = error * 2;

		if (error2 > -deltay) {
			error -= deltay;
			x0 += 1;
		}

		if (error2 < deltax) {
			error += deltax;
			y0 += signy;
		}
	}

	return 0;
}

static int printmask(int *mask, int w, int h, const char *path)
{
	struct nd_image img;
	int x, y;

	nd_imgcreate(&img, w, h, ND_PF_GRAYSCALE);

	for (y = 0; y < h; ++y)
		for (x = 0; x < w; ++x)
			img.data[y * w + x] = mask[y * w + x] ? 255 : 0;

	nd_imgwrite(&img, path);

	nd_imgdestroy(&img);

	return 0;
}

static int imgfindlines(struct nd_image *img, double *lines, int *linescount,
	enum LINEORIENTATION lineorient)
{
	int *mask;
	double *hlines;
	int maxlinescount;
	int linen;

	if (nd_imgnormalize(img, 1, 1) < 0)
		return (-1);

	if ((mask = malloc(sizeof(int) * img->w * img->h)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if (ed_canny(img, mask, -1.0, -1.0) < 0)
		return (-1);

////////////////////////////////////////////////

	const char *path;

	switch(lineorient) {
	case LO_TOP:
		path = "top.png";
		break;
	case LO_BOTTOM:
		path = "bot.png";
		break;
	case LO_LEFT:
		path = "left.png";
		break;
	case LO_RIGHT:
		path = "right.png";
		break;
	}

	printmask(mask, img->w, img->h, path);

////////////////////////////////////////////////


	if ((hlines = malloc(sizeof(double) * (*linescount) * 2)) == NULL) {
		safefree((void **)&mask);

		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if (ed_hough(mask, img->w, img->h, M_PI / 180.0, hlines, (*linescount))
		< 0) {
		safefree((void **)&mask);
		safefree((void **)&lines);
		
		return (-1);
	}

	maxlinescount = *linescount;
	
	*linescount = 0;
	for (linen = 0; linen < maxlinescount; ++linen) {
		double rho, theta;
	
		theta = hlines[linen * 2];
		rho = hlines[linen * 2 + 1];
		
		if (checkline(img, theta, rho, lineorient)) {
			lines[*linescount * 2] = theta;
			lines[*linescount * 2 + 1] = rho;
		
			++(*linescount);
		}
	}
	
	safefree((void **)&mask);
	safefree((void **)&hlines);

	return 0;
}

static int lintersect(struct lineseg *l1, struct lineseg *l2,
	double *interx, double *intery)
{
	double m1, c1, m2, c2;
	double dx, dy;
	int isvert1, isvert2;

	isvert1 = isvert2 = 0;

	dx = l1->x1 - l1->x0;

	if (dx < 0.001)
		isvert1 = 1;
	else {
		dy = l1->y1 - l1->y0;
		m1 = dy / dx;
		c1 = l1->y0 - m1 * l1->x0;
	}

	dx = l2->x1 - l2->x0;

	if (dx < 0.001)
		isvert2 = 1;
	else {
		dy = l2->y1 - l2->y0;
		m2 = dy / dx;
		c2 = l2->y0 - m2 * l2->x0;
	}

	if (isvert1 || isvert2) {
		if (isvert1 && isvert2)
			return 0;
		if (isvert1 && !isvert2) {	
			*interx = (double) (l1->x0);
			*intery = (double) (l2->y0) + m2 * (double) (l1->x0);
		
			return 1;
		}

		if (isvert2 && !isvert1) {
			*interx = (double) (l2->x0);
			*intery = (double) (l1->y0) + m1 * (double) (l2->x0);
			
			return 1;
		}
	}
	else {	
		if( (m1 - m2) == 0)
			return 0;
		else {
			*interx = (c2 - c1) / (m1 - m2);
			*intery = m1 * *interx + c1;
		
			return 1;
		}
	}

	return 0;
}

static int getparallel(double *lines1, int lines1count, double *lines2,
	int lines2count, int *li1, int *li2)
{
	int minisum;
	int i1, i2;
	int minidif;

	minisum = lines1count + lines2count - 2;
	minidif = (lines1count > lines2count) ? lines1count : lines2count;

	for (i1 = 0; i1 < lines1count; ++i1)
		for (i2 = 0; i2 < lines2count; ++i2)
			if (fabs(lines1[i1 * 2] - lines2[i2 * 2]) < MAXANGDIF
				&& (i1 + i2 < minisum
				|| (i1 + i2 == minisum
				&& abs(i2 - i1) < minidif))) {
				*li1 = i1;
				*li2 = i2;
				minisum = i1 + i2;
				minidif = abs(i2 - i1);
			}

	return 1;
}

int linetoseg(double theta, double rho, int imgw, int imgh,
	struct lineseg *lseg)
{
	if (theta < THETAMINVAL) {	
		lseg->x0 = rho;
		lseg->y0 = 0.0;
		lseg->x1 = rho;
		lseg->y1 = (double) imgh;
	}
	else {		
		lseg->x0 = 0.0;
		lseg->y0 = rho / sin(theta);
		lseg->x1 = (double) imgw;
		lseg->y1 = (rho - imgw * cos(theta)) / sin(theta);
	}

	return 0;
}

int ed_findborder(struct nd_image *img, double inpoints[8])
{
	struct nd_image *imgparts;
	double *lines0;
	double *lines1;
	int linescount;
	struct lineseg ltop, lbot;
	struct lineseg lleft, lright;
	int li1, li2;

	assert(img != NULL && inpoints != NULL);

// split image into two parts by a horizontal axis
	if ((imgparts = malloc(sizeof(struct nd_image) * 2)) == NULL)
		return (-1);

	if (imghsplit(img, imgparts, 2) < 0)
		return (-1);

	if ((lines0 = malloc(sizeof(double) * 2 * MAXLINES)) == NULL)
		return (-1);

	if ((lines1 = malloc(sizeof(double) * 2 * MAXLINES)) == NULL)
		return (-1);

	li1 = li2 = 0;
	
	linescount = MAXLINES;
	if (imgfindlines(imgparts + 0, lines0, &linescount, LO_TOP) < 0)
		return (-1);

	if (linescount == 0)
		return (-1);

	linescount = MAXLINES;
	if (imgfindlines(imgparts + 1, lines1, &linescount, LO_BOTTOM) < 0)
		return (-1);

	if (linescount == 0)
		return (-1);

	getparallel(lines0, linescount, lines1, linescount, &li1, &li2);

	linetoseg(lines0[li1 * 2], lines0[li1 * 2 + 1], img->w, img->h, &ltop);
	linetoseg(lines1[li2 * 2], lines1[li2 * 2 + 1], img->w, img->h, &lbot);

	nd_imgtorgb(imgparts);
	imgdrawline(imgparts, &ltop);

	nd_imgwrite(imgparts, "toplines.png");

	lbot.y0 += img->h / 2;
	lbot.y1 += img->h / 2;

// split image into four parts by a vertical axis
	if ((imgparts = malloc(sizeof(struct nd_image) * 5)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}

	if (imgvsplit(img, imgparts, 5) < 0)
		return (-1);

	linescount = MAXLINES;
	if (imgfindlines(imgparts + 0, lines0, &linescount, LO_LEFT) < 0)
		return (-1);

	if (linescount == 0)
		return (-1);

	linescount = MAXLINES;
	if (imgfindlines(imgparts + 4, lines1, &linescount, LO_RIGHT) < 0)
		return (-1);
	
	if (linescount == 0)
		return (-1);

	getparallel(lines0, linescount, lines1, linescount, &li1, &li2);

	linetoseg(lines0[li1 * 2], lines0[li1 * 2 + 1],
		img->w, img->h, &lleft);
	linetoseg(lines1[li2 * 2], lines1[li2 * 2 + 1],
		img->w, img->h, &lright);
	
	lright.x0 += 4 * img->w / 5;
	lright.x1 += 4 * img->w / 5;

	lintersect(&lleft, &ltop, inpoints + 0, inpoints + 1);
	lintersect(&lright, &ltop, inpoints + 2, inpoints + 3);
	lintersect(&lright, &lbot, inpoints + 4, inpoints + 5);
	lintersect(&lleft, &lbot, inpoints + 6, inpoints + 7);

	return 0;
}
