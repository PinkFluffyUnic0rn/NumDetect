#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "nd_error.h"
#include "nd_vecmat.h"

#include "em3d.h"

int em_allocgaussmix(struct em_gaussmix *gm, int pn)
{
	gm->pn = pn;

	if ((gm->m = malloc(sizeof(struct nd_vector3) * pn)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto mmamllocerror;
	}

	if ((gm->c = malloc(sizeof(double) * pn)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto cmallocerror;
	}

	if ((gm->w = malloc(sizeof(double) * pn)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto wmallocerror;
	}

	return 0;

wmallocerror:
	free(gm->w);
cmallocerror:
	free(gm->m);
mmamllocerror:
	return (-1);
}

double em_probability_gauss3d(const struct nd_vector3 *m,
	double c, double cidet,
	const struct nd_vector3 *x)
{
	struct nd_vector3 md;
	double s;

	nd_v3lincomb(x, -1.0, m, &md);

	s = (pow(md.x, 2.0) + pow(md.y, 2.0) + pow(md.z, 2.0)) / c;

	return exp(-0.5 * s - 0.5 * log(cidet) - 3.0 / 2.0 * log(2.0 * M_PI));
}

int em_gauss3dstep(const struct nd_vector3 *x, int xn,
	struct em_gaussmix *gm, double *gout, double *ll)
{
	double **g;
	double *cdet;
	int i, j;

	// precomputing determinant
	if ((cdet = malloc(sizeof(double) * gm->pn)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto cdetmallocerror;
	}

	for (i = 0; i < gm->pn; ++i)
		cdet[i] = pow(gm->c[i], 3.0); 

	// allocating g[][]
	if ((g = malloc(sizeof(double *) * xn)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto gmallocerror;
	}

	g[0] = gout;

	for (i = 1; i < xn; ++i)
		g[i] = g[i - 1] + gm->pn;
	
	// computing g[][] and logarithm of the current
	// likelihood function value
	*ll = 0.0;

	for (i = 0; i < xn; ++i) {
		double s;

		s = 0.0;
		for (j = 0; j < gm->pn; ++j) {
			s += gm->w[j] * em_probability_gauss3d(gm->m + j,
				gm->c[j], cdet[j], x + i);
		}

		for (j = 0; j < gm->pn; ++j) {
			g[i][j] = gm->w[j] * em_probability_gauss3d(gm->m + j,
					gm->c[j], cdet[j], x + i) / s;
		}

		*ll += log(s);
	}

	// computing w[]
	for (i = 0; i < gm->pn; ++i) {
		gm->w[i] = 0.0;

		for (j = 0; j < xn; ++j)
			gm->w[i] += g[j][i];

		gm->w[i] /= xn;
	}

	// computing m[]
	for (i = 0; i < gm->pn; ++i) {	
		gm->m[i].x = gm->m[i].y = gm->m[i].z = 0.0;

		for (j = 0; j < xn; ++j)
			nd_v3lincomb(gm->m + i, g[j][i], x + j, gm->m + i);

		nd_scalarv3mult(gm->m + i, 1.0 / (xn * gm->w[i]), gm->m + i);
	}

	// computing v[]
	for (i = 0; i < gm->pn; ++i) {
		gm->c[i] = 0.0;

		for (j = 0; j < xn; ++j) {
			struct nd_vector3 tmpv;

			nd_v3lincomb(x + j, -1.0, gm->m + i, &tmpv);
		
			gm->c[i] += g[j][i] * (pow(tmpv.x, 2.0)
				+ pow(tmpv.y, 2.0) + pow(tmpv.z, 2.0));
		}

		gm->c[i] /= (xn * gm->w[i]);

		gm->c[i] = (gm->c[i] < EM_MINVARIANCE)
			? EM_MINVARIANCE : gm->c[i];
	}
	
	free(cdet);
	free(g);

	return 0;

gmallocerror:
	free(cdet);
cdetmallocerror:
	return (-1);
}

int em_gauss3d(const struct nd_vector3 *x, int xn,
	struct em_gaussmix *gm, int pn, double *gout, double *llout,
	double eps)
{
	double prevll, ll;
	double *g;
	int i;

	gm->pn = pn;

	for (i = 0; i < pn; ++i) {
		gm->m[i] = x[rand() % xn];
			
		gm->c[i] = 1.0;
		
		gm->w[i] = 1.0 / gm->pn;
	}

	if (gout != NULL)
		g = gout;
	else if ((g = malloc(sizeof(double) * xn * gm->pn)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);
	}
	
	if (em_gauss3dstep(x, xn, gm, g, &ll) < 0)
		return (-1);
	
	do {
		prevll = ll;

		if (em_gauss3dstep(x, xn, gm, g, &ll) < 0)
			return (-1);

	} while (fabs(ll - prevll) >= eps);

	if (gout == NULL)
		free(g);

	if (llout != NULL)
		*llout = ll;
	
	return 0;
}

int em3dprintp(struct nd_vector3 *m, struct nd_matrix3 *c, double *w, int pn)
{
	int i;

	printf("pn: %d\n", pn);
	
	printf("m:\n");
	for (i = 0; i < pn; ++i)
		printf("%f %f %f\n", m[i].x, m[i].y, m[i].z);
	printf("\n");

	printf("c:\n");
	for (i = 0; i < pn; ++i)
		printf("%f %f %f\n%f %f %f\n%f %f %f\n\n",
			c[i]._11, c[i]._12, c[i]._13,
			c[i]._21, c[i]._22, c[i]._23,
			c[i]._31, c[i]._32, c[i]._33);
		
	printf("w: ");
	for (i = 0; i < pn; ++i)
		printf("%f%c", w[i], i != (pn - 1) ? ' ' : '\n');

	return 0;
}
