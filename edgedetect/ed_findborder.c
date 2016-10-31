#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nd_image.h"
#include "nd_error.h"
#include "ed_edgedetect.h"

#include <gtk/gtk.h>
#include <cairo.h>

#define HOUGHMAXLINES 100
#define MAXANGDIF 1.0 * M_PI / 180.0

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

static int checkline(struct nd_image *img, double rho, double theta,
	enum LINEORIENTATION lo)
{	
	if ((lo == LO_TOP || lo == LO_BOTTOM)
		&& (theta > M_PI * 0.25 && theta < M_PI * 0.75))
		return 1;	

	if ((lo == LO_LEFT || lo == LO_RIGHT)
		&& !(theta > M_PI * 0.25 && theta < M_PI * 0.75)) {	
/*
		double vbandlen;
		double sumborder;
		double suminner;

		vbandlen = img->w * 0.025 * ((lo == LO_LEFT) ? 1.0 : -1.0);

		sumborder = bandsum(img, rho - 0.01,
			rho + vbandlen * 1.0 + 0.01, theta);
			
		suminner = bandsum(img, rho + vbandlen * 1.0 - 0.01,
			rho + vbandlen * 2.0 + 0.01, theta);

		return (sumborder > suminner) ? 1 : 0;
*/
		double vbandlen = img->w * 0.125;
		double valdif;
		double sumborder;
		double suml;
		double sumr;

		suml = bandsum(img, rho - vbandlen * 1.5 - 0.01,
			rho - vbandlen * 0.5 + 0.01, theta);
		
		sumr = bandsum(img, rho + vbandlen * 0.5 - 0.01,
			rho + vbandlen * 1.5 + 0.01, theta);
		
		sumborder = bandsum(img, rho - vbandlen * 0.5 - 0.01,
			rho + vbandlen * 0.5 + 0.01, theta);

		valdif = ((lo == LO_LEFT) ? sumr : suml) - sumborder;

		if (theta > M_PI / 2.0)
			valdif *= -1.0;

		return valdif > 0.0;
	}

	return 0;
}

static int imgfindlines(struct nd_image *img, double *lines, int *linescount,
	enum LINEORIENTATION lineorient)
{
	int *mask;
	double *hlines;
	int maxlinescount;
	int linen;

	if (lineorient == LO_LEFT || lineorient == LO_RIGHT)
		if (ed_removelowfreq(img, 0.0, 0.0) < 0)
			return (-1);

	if (lineorient == LO_TOP || lineorient == LO_BOTTOM)
		if (ed_removelowfreq(img, 0.0, 0.0) < 0)
			return (-1);

	if (nd_imgnormalize(img, 1, 1) < 0)
		return (-1);
/*
	int i, j;
	for (i = 0; i < img->h; ++i)
		for (j = 0; j < img->w; ++j)
			if (!(img->data[i * img->w + j] > 0.0))
				printf("%lf\n", img->data[i * img->w + j]);

*/
	if ((mask = malloc(sizeof(int) * img->w * img->h)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if (ed_canny(img, mask, -1.0, -1.0) < 0)
		return (-1);

	if ((hlines = malloc(sizeof(double) * HOUGHMAXLINES * 2)) == NULL) {
		safefree((void **)&mask);

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	if (ed_hough(mask, img->w, img->h,
		M_PI / 180.0, hlines, HOUGHMAXLINES) < 0) {
		safefree((void **)&mask);
		safefree((void **)&lines);
		
		return (-1);
	}

	maxlinescount = *linescount;

	*linescount = 0;
	for (linen = 0; linen < HOUGHMAXLINES && *linescount < maxlinescount;
		++linen) {
		double rho, theta;
	
		theta = hlines[linen * 2];
		rho = hlines[linen * 2 + 1];
		
		if (checkline(img, rho, theta, lineorient)) {
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

int ed_findborder(struct nd_image *img, double inpoints[8])
{
	struct nd_image *imgparts;
	double *linestop;
	double *linesbot;
	double *linesleft;
	double *linesright;
	int linescount;
	struct lineseg ltop;
	struct lineseg lbot;
	struct lineseg lleft;
	struct lineseg lright;
	int li1, li2;
	double rho, theta;
/*
	ed_removelowfreq(&img, 0.05, 0.0);

	nd_imgwrite("test.png", &img);
	exit(1);
*/

// split image into two parts by a horizontal axis
	if ((imgparts = malloc(sizeof(struct nd_image) * 2)) == NULL)
		return (-1);

	if (imghsplit(img, imgparts, 2) < 0)
		return (-1);

	linescount = 100;
	linestop = malloc(sizeof(double) * 2 * linescount);
	linesbot = malloc(sizeof(double) * 2 * linescount);

	imgfindlines(imgparts + 0, linestop, &linescount, LO_TOP);
	imgfindlines(imgparts + 1, linesbot, &linescount, LO_BOTTOM);

	li1 = li2 = 0;
	getparallel(linestop, linescount, linesbot, linescount, &li1, &li2);

	theta = linestop[li1 * 2];
	rho = linestop[li1 * 2 + 1];
	
	if (theta < 0.00001) {	
		ltop.x0 = rho;
		ltop.y0 = 0.0;
		ltop.x1 = rho;
		ltop.y1 = (double) img->h;
	}
	else {		
		ltop.x0 = 0.0;
		ltop.y0 = rho / sin(theta);
		ltop.x1 = img->w;
		ltop.y1 = (rho - img->w * cos(theta)) / sin(theta);
	}

	theta = linesbot[li2 * 2];
	rho = linesbot[li2 * 2 + 1];

	if (theta < 0.00001) {	
		lbot.x0 = rho;
		lbot.y0 = 0.0;
		lbot.x1 = rho;
		lbot.y1 = img->h;
	}
	else {
		lbot.x0 = 0.0;
		lbot.y0 = rho / sin(theta);
		lbot.x1 = img->w;
		lbot.y1 = (rho - img->w * cos(theta)) / sin(theta);
	}

	lbot.y0 += img->h / 2;
	lbot.y1 += img->h / 2;

// split image into four parts by a vertical axis
	if ((imgparts = malloc(sizeof(struct nd_image) * 5)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if (imgvsplit(img, imgparts, 5) < 0)
		return (-1);

	linescount = 100;
	
	if ((linesleft = malloc(sizeof(double) * 2 * linescount)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((linesright = malloc(sizeof(double) * 2 * linescount)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	imgfindlines(imgparts + 0, linesleft, &linescount, LO_LEFT);
	imgfindlines(imgparts + 4, linesright, &linescount, LO_RIGHT);
	
	getparallel(linesleft, linescount, linesright, linescount, &li1, &li2);

	theta = linesleft[li1 * 2];
	rho = linesleft[li1 * 2 + 1];

	if (theta < 0.00001) {	
		lleft.x0 = rho;
		lleft.y0 = 0.0;
		lleft.x1 = rho;
		lleft.y1 = (double) img->h;
	}
	else {		
		lleft.x0 = 0.0;
		lleft.y0 = rho / sin(theta);
		lleft.x1 = img->w;
		lleft.y1 = (rho - img->w * cos(theta)) / sin(theta);
	}

	theta = linesright[li2 * 2];
	rho = linesright[li2 * 2 + 1];

	if (theta < 0.00001) {	
		lright.x0 = rho;
		lright.y0 = 0.0;
		lright.x1 = rho;
		lright.y1 = img->h;
	}
	else {
		lright.x0 = 0.0;
		lright.y0 = rho / sin(theta);
		lright.x1 = img->w;
		lright.y1 = (rho - img->w * cos(theta)) / sin(theta);
	}
	
	lright.x0 += 4 * img->w / 5;
	lright.x1 += 4 * img->w / 5;

	lintersect(&lleft, &ltop, inpoints + 0, inpoints + 1);
	lintersect(&lright, &ltop, inpoints + 2, inpoints + 3);
	lintersect(&lright, &lbot, inpoints + 4, inpoints + 5);
	lintersect(&lleft, &lbot, inpoints + 6, inpoints + 7);

	return 0;
}
