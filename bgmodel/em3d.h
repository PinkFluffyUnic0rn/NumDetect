#ifndef ND_EM3D
#define ND_EM3D

#include "nd_vecmat.h"

#define EM_MINVARIANCE 1.0e-2

struct em_gaussmix {
	struct nd_vector3 *m;
	double *c;
	double *w;
	int pn;
};

int em_allocgaussmix(struct em_gaussmix *gm, int pn);

double em_probability_gauss3d(const struct nd_vector3 *m,
	double c, double cidet,
	const struct nd_vector3 *x);

int em_gauss3d(const struct nd_vector3 *x, int xn,
	struct em_gaussmix *gm, int pn, double *gout, double *llout,
	double eps);

int em3dprintp(struct nd_vector3 *m, struct nd_matrix3 *c, double *w, int pn);

#endif
