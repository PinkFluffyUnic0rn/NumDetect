#include <math.h>
#include <string.h>

#include "nd_vecmat.h"
#include "nd_error.h"

int nd_m3scale(struct nd_matrix3 *r, double x, double y, double z)
{
	if (r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	r->_11 = x;	r->_12 = 0.0;	r->_13 = 0.0;
	r->_21 = 0.0;	r->_22 = y;	r->_23 = 0.0;
	r->_31 = 0.0;	r->_32 = 0.0;	r->_33 = z;

	return 0;
}

int nd_m3rotatex(struct nd_matrix3 *r, double ang)
{
	if (r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	r->_11 = 1.0;	r->_12 = 0.0;		r->_13 = 0.0;
	r->_21 = 0.0;	r->_22 = cos(ang);	r->_23 = sin(ang);
	r->_31 = 0.0;	r->_32 = -sin(ang);	r->_33 = cos(ang);

	return 0;
}

int nd_m3rotatey(struct nd_matrix3 *r, double ang)
{
	if (r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	r->_11 = cos(ang);	r->_12 = 0.0;		r->_13 = sin(ang);
	r->_21 = 0.0;		r->_22 = 1.0;		r->_23 = 0.0;
	r->_31 = -sin(ang);	r->_32 = 0.0;		r->_33 = cos(ang);

	return 0;
}

int nd_m3rotatez(struct nd_matrix3 *r, double ang)
{
	if (r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	r->_11 = cos(ang);	r->_12 = sin(ang);	r->_13 = 0.0;
	r->_21 = -sin(ang);	r->_22 = cos(ang);	r->_23 = 0.0;
	r->_31 = 0.0;		r->_32 = 0.0;		r->_33 = 1.0;

	return 0;
}

int nd_m3translate(struct nd_matrix3 *r, double x, double y, double z)
{
	if (r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	r->_11 = 1.0;	r->_12 = 0.0;	r->_13 = 0.0;
	r->_21 = 0.0;	r->_22 = 1.0;	r->_23 = 0.0;
	r->_31 = x;	r->_32 = y;	r->_33 = z;

	return 0;
}

int nd_v3m3mult(const struct nd_vector3 *v, const struct nd_matrix3 *m,
	struct nd_vector3 *r)
{
	struct nd_vector3 tmp;

	if (v == NULL || m == NULL || r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	tmp.x = (v->x * m->_11) + (v->y * m->_12) + (v->z * m->_13);
	tmp.y = (v->x * m->_21) + (v->y * m->_22) + (v->z * m->_23);
	tmp.z = (v->x * m->_31) + (v->y * m->_32) + (v->z * m->_33);

	memcpy(r, &tmp, sizeof(struct nd_vector3));

	return 0;
}

int nd_m3mult(const struct nd_matrix3 *m0, const struct nd_matrix3 *m1,
	struct nd_matrix3 *r)
{
	struct nd_matrix3 tmp;

	if (m0 == NULL || m1 == NULL || r == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	tmp._11 = (m1->_11 * m0->_11)
		+ (m1->_21 * m0->_12) + (m1->_31 * m0->_13);
	tmp._12 = (m1->_12 * m0->_11)
		+ (m1->_22 * m0->_12) + (m1->_32 * m0->_13);
	tmp._13 = (m1->_13 * m0->_11)
		+ (m1->_23 * m0->_12) + (m1->_33 * m0->_13);

	tmp._21 = (m1->_11 * m0->_21)
		+ (m1->_21 * m0->_22) + (m1->_31 * m0->_23);
	tmp._22 = (m1->_12 * m0->_21)
		+ (m1->_22 * m0->_22) + (m1->_32 * m0->_23);
	tmp._23 = (m1->_13 * m0->_21)
		+ (m1->_23 * m0->_22) + (m1->_33 * m0->_23);

	tmp._31 = (m1->_11 * m0->_31)
		+ (m1->_21 * m0->_32) + (m1->_31 * m0->_33);
	tmp._32 = (m1->_12 * m0->_31)
		+ (m1->_22 * m0->_32) + (m1->_32 * m0->_33);
	tmp._33 = (m1->_13 * m0->_31)
		+ (m1->_23 * m0->_32) + (m1->_33 * m0->_33);

	memcpy(r, &tmp, sizeof(struct nd_matrix3));

	return 0;
}

static int nd_chooseroworder3v(double m[3][3], double v[3])
{
	if (fabs(m[0][0]) >= EPSILON && fabs(m[1][1]) >= EPSILON
		&& fabs(m[2][2]) >= EPSILON)
		return 0;

	if (fabs(m[0][0]) >= EPSILON && fabs(m[2][1]) >= EPSILON
		&& fabs(m[1][2]) >= EPSILON) {
		double tmpr[3];
		double tmps;

		memcpy(tmpr, m[2], sizeof(double) * 3);
		memcpy(m[2], m[1], sizeof(double) * 3);
		memcpy(m[1], tmpr, sizeof(double) * 3);

		tmps = v[2];
		v[2] = v[1];
		v[1] = tmps;

		return 0;
	}
	
	if (fabs(m[1][0]) >= EPSILON && fabs(m[0][1]) >= EPSILON
		&& fabs(m[2][2]) >= EPSILON) {
		double tmpr[3];
		double tmps;

		memcpy(tmpr, m[1], sizeof(double) * 3);
		memcpy(m[1], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		tmps = v[1];
		v[1] = v[0];
		v[0] = tmps;

		return 0;
	}

	if (fabs(m[1][0]) >= EPSILON && fabs(m[2][1]) >= EPSILON
		&& fabs(m[1][2]) >= EPSILON) {
		double tmpr[3];
		double tmps;

		memcpy(tmpr, m[1], sizeof(double) * 3);
		memcpy(m[1], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, m[2], sizeof(double) * 3);
		memcpy(m[2], m[1], sizeof(double) * 3);
		memcpy(m[1], tmpr, sizeof(double) * 3);

		tmps = v[1];
		v[1] = v[0];
		v[0] = tmps;

		tmps = v[2];
		v[2] = v[1];
		v[1] = tmps;

		return 0;
	}

	if (fabs(m[2][0]) >= EPSILON && fabs(m[0][1]) >= EPSILON
		&& fabs(m[1][2]) >= EPSILON) {
		double tmpr[3];
		double tmps;

		memcpy(tmpr, m[1], sizeof(double) * 3);
		memcpy(m[1], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, m[2], sizeof(double) * 3);
		memcpy(m[2], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		tmps = v[1];
		v[1] = v[0];
		v[0] = tmps;

		tmps = v[2];
		v[2] = v[0];
		v[0] = tmps;

		return 0;
	}

	if (fabs(m[2][0]) >= EPSILON && fabs(m[1][1]) >= EPSILON
		&& fabs(m[0][2]) >= EPSILON) {
		double tmp[3];
		double tmps;

		memcpy(tmp, m[2], sizeof(double) * 3);
		memcpy(m[2], m[0], sizeof(double) * 3);
		memcpy(m[0], tmp, sizeof(double) * 3);

		tmps = v[2];
		v[2] = v[0];
		v[0] = tmps;

		return 0;
	}

	return (-1);
}

int nd_m3nonhomsolve(const struct nd_matrix3 *m, const struct nd_vector3 *v,
	double *r)
{
	double tmpM[3][3];
	double tmpV[3];
	double coef;
	int i, j, k;

	memcpy(tmpM, m, sizeof(struct nd_matrix3));
	memcpy(tmpV, v, sizeof(struct nd_vector3));

	if (nd_chooseroworder3v(tmpM, tmpV) < 0) {
		nd_seterror(ND_WRONGMATRIX);
		return (-1);
	}

	for (k = 0; k < 3; ++k) {
		coef = 1/tmpM[k][k];

		for (j = 0; j < 3; ++j)
			tmpM[k][j] *= coef;
		tmpV[k] *= coef;

		for (i = 0; i < 3; ++i)
			if (i != k) {
				coef = -(tmpM[i][k] / tmpM[k][k]);
				
				for (j = 0; j < 3; ++j)
					tmpM[i][j] += tmpM[k][j] * coef;
			
				tmpV[i] += tmpV[k] * coef;
			}
	}

	memcpy(r, tmpV, sizeof(double) * 3);
	
	return 0;
}


static int nd_chooseroworder3m(double m[3][3], double mi[3][3])
{
	if (fabs(m[0][0]) >= EPSILON && fabs(m[1][1]) >= EPSILON
		&& fabs(m[2][2]) >= EPSILON)
		return 0;

	if (fabs(m[0][0]) >= EPSILON && fabs(m[2][1]) >= EPSILON
		&& fabs(m[1][2]) >= EPSILON) {
		double tmpr[3];

		memcpy(tmpr, m[2], sizeof(double) * 3);
		memcpy(m[2], m[1], sizeof(double) * 3);
		memcpy(m[1], tmpr, sizeof(double) * 3);
		
		memcpy(tmpr, mi[2], sizeof(double) * 3);
		memcpy(mi[2], mi[1], sizeof(double) * 3);
		memcpy(mi[1], tmpr, sizeof(double) * 3);

		return 0;
	}
	
	if (fabs(m[1][0]) >= EPSILON && fabs(m[0][1]) >= EPSILON
		&& fabs(m[2][2]) >= EPSILON) {
		double tmpr[3];

		memcpy(tmpr, m[1], sizeof(double) * 3);
		memcpy(m[1], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, mi[1], sizeof(double) * 3);
		memcpy(mi[1], mi[0], sizeof(double) * 3);
		memcpy(mi[0], tmpr, sizeof(double) * 3);

		return 0;
	}

	if (fabs(m[1][0]) >= EPSILON && fabs(m[2][1]) >= EPSILON
		&& fabs(m[1][2]) >= EPSILON) {
		double tmpr[3];

		memcpy(tmpr, m[1], sizeof(double) * 3);
		memcpy(m[1], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, m[2], sizeof(double) * 3);
		memcpy(m[2], m[1], sizeof(double) * 3);
		memcpy(m[1], tmpr, sizeof(double) * 3);

		memcpy(tmpr, mi[1], sizeof(double) * 3);
		memcpy(mi[1], mi[0], sizeof(double) * 3);
		memcpy(mi[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, mi[2], sizeof(double) * 3);
		memcpy(mi[2], mi[1], sizeof(double) * 3);
		memcpy(mi[1], tmpr, sizeof(double) * 3);

		return 0;
	}

	if (fabs(m[2][0]) >= EPSILON && fabs(m[0][1]) >= EPSILON
		&& fabs(m[1][2]) >= EPSILON) {
		double tmpr[3];

		memcpy(tmpr, m[1], sizeof(double) * 3);
		memcpy(m[1], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, m[2], sizeof(double) * 3);
		memcpy(m[2], m[0], sizeof(double) * 3);
		memcpy(m[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, mi[1], sizeof(double) * 3);
		memcpy(mi[1], mi[0], sizeof(double) * 3);
		memcpy(mi[0], tmpr, sizeof(double) * 3);

		memcpy(tmpr, mi[2], sizeof(double) * 3);
		memcpy(mi[2], mi[0], sizeof(double) * 3);
		memcpy(mi[0], tmpr, sizeof(double) * 3);

		return 0;
	}

	if (fabs(m[2][0]) >= EPSILON && fabs(m[1][1]) >= EPSILON
		&& fabs(m[0][2]) >= EPSILON) {
		double tmp[3];

		memcpy(tmp, m[2], sizeof(double) * 3);
		memcpy(m[2], m[0], sizeof(double) * 3);
		memcpy(m[0], tmp, sizeof(double) * 3);

		memcpy(tmp, mi[2], sizeof(double) * 3);
		memcpy(mi[2], mi[0], sizeof(double) * 3);
		memcpy(mi[0], tmp, sizeof(double) * 3);

		return 0;
	}

	return (-1);
}

int nd_m3inverse(const struct nd_matrix3 *m, struct nd_matrix3 *r)
{
	double invM[3][3], tmpM[3][3];
	double coef;
	int i, j, k;

	memcpy(tmpM, m, sizeof(struct nd_matrix3));

	for (i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			invM[i][j] = (i == j) ? 1.0f : 0.0;

	if (nd_chooseroworder3m(tmpM, invM) < 0) {
		nd_seterror(ND_WRONGMATRIX);
		return (-1);
	}

	for (k = 0; k < 3; ++k) {
		coef = 1/tmpM[k][k];

		for ( j = 0; j < 3; ++j ) {
			tmpM[k][j] *= coef;
			invM[k][j] *= coef;
		}

		for (i = 0; i < 3; ++i)
			if (i != k) {
				coef = -(tmpM[i][k] / tmpM[k][k]);
				
				for (j = 0; j < 3; ++j) {
					tmpM[i][j] += tmpM[k][j]*coef;
					invM[i][j] += invM[k][j]*coef;
				}
			}
	}

	memcpy(r, invM, sizeof(struct nd_matrix3));
	
	return 0;
}
