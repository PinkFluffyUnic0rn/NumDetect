#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cv.h>
#include <highgui.h>

#define ESC 27
#define SPACE 32
#define KEY_P 112

struct line {
	double x0;
	double y0;
	double x1;
	double y1;
};

struct perspquad {
	CvPoint2D32f points[4];
};

CvMat *persp(CvMat *in, struct perspquad *perspquad)
{
	CvPoint2D32f outp[4];
	CvMat *mpersp;
	CvMat *copy;
	
	copy = cvCreateMat(in->rows, in->cols, CV_8UC1);
	mpersp = cvCreateMat(3, 3, CV_32FC1);

	outp[0] = cvPoint2D32f(in->cols, in->rows);
	outp[2] = cvPoint2D32f(in->cols, 0);
	outp[1] = cvPoint2D32f(0, in->rows);
	outp[3] = cvPoint2D32f(0, 0);

	cvGetPerspectiveTransform(perspquad->points, outp, mpersp);

	cvWarpPerspective(in, copy, mpersp,
		CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll(0));
	
	cvReleaseMat(&mpersp);

	return copy;
}

void drawhoughline(CvMat *img, double rho, double theta)
{
	CvPoint pt1;
	CvPoint pt2;
	double a; 
	double b;
	
	a = cos(theta);
	b = sin(theta);

	pt1.x = cvRound(0.0);
	pt1.y = cvRound(rho / b);

	pt2.x = cvRound((double) (img->cols));
	pt2.y = cvRound((rho - (double) (img->cols) * a) / b);

	cvLine(img, pt1, pt2, CV_RGB(255, 255, 255), 1, 8, 0);
}

int main(int argc, char **argv)
{
	CvMat *img;

	img = cvLoadImageM(argv[1], CV_LOAD_IMAGE_GRAYSCALE);

	int x0, y0;
	int x1, y1;
	int x2, y2;
	int x3, y3;

	scanf("%d %d\n", &x0, &y0);
	scanf("%d %d\n", &x1, &y1);
	scanf("%d %d\n", &x2, &y2);
	scanf("%d %d\n", &x3, &y3);

	struct perspquad pq;

	pq.points[0].x = x0;
	pq.points[0].y = y0;

	pq.points[1].x = x1;
	pq.points[1].y = y1;
	
	pq.points[2].x = x2;
	pq.points[2].y = y2;
	
	pq.points[3].x = x3;
	pq.points[3].y = y3;
	
	img = persp(img, &pq);

	cvSaveImage(argv[2], img, 0);
	
	cvReleaseMat(&img);

	return 0;
}

