#ifndef VECMAT_H
#define VECMAT_H

#define EPSILON 0.00001

struct nd_vector2 {
	double x, y;
};

struct nd_vector3 {
	double x, y, z;
};

struct nd_matrix3
{
	double _11, _12, _13,
		_21, _22, _23,
		_31, _32, _33;
};

int nd_m3scale(struct nd_matrix3 *r, double x, double y, double z);

int nd_m3rotatex(struct nd_matrix3 *r, double ang);

int nd_m3rotatey(struct nd_matrix3 *r, double ang);

int nd_m3rotatez(struct nd_matrix3 *r, double ang);

int nd_m3translate(struct nd_matrix3 *r, double x, double y, double z);

int nd_v3m3mult(const struct nd_vector3 *v, const struct nd_matrix3 *m,
	struct nd_vector3 *r);

int nd_m3mult(const struct nd_matrix3 *m0, const struct nd_matrix3 *m1,
	struct nd_matrix3 *r);

int nd_m3nonhomsolve(const struct nd_matrix3 *m, const struct nd_vector3 *v,
	double *r);

int nd_m3inverse(const struct nd_matrix3 *m, struct nd_matrix3 *r);

#endif
