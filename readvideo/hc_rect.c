#include "hc_rect.h"
#include "nd_error.h"

static void hc_safefree(void **p)
{
	if (*p != NULL) {
		free(*p);
		*p = NULL;
	}
}

static int hc_rectinter(const struct hc_rect *r0, const struct hc_rect *r1)
{
	int xninter;
	int yninter;

	xninter = (r0->x0 > r1->x1) || (r0->x1 < r1->x0);
	yninter = (r0->y0 > r1->y1) || (r0->y1 < r1->y0);

	return !(xninter || yninter);
}

static int hc_recttograph(const struct hc_rect *r, int recc,
	int ***ve, int **vec)
{
	int r0, r1;
	int *e;
	int ec;

	if ((*ve = malloc(sizeof(int *) * recc)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((*vec = malloc(sizeof(int) * recc)) == NULL) {
		hc_safefree((void **) ve);

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	ec = 0;
	for (r0 = 0; r0 < recc; ++r0)
		for (r1 = 0; r1 < recc; ++ r1)
			if (hc_rectinter(r + r0, r + r1))
				++ec;

	if ((e = (int *) malloc(sizeof(int) * ec)) == NULL) {
		hc_safefree((void **)ve);
		hc_safefree((void **)vec);
		
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	ec = 0;
	for (r0 = 0; r0 < recc; ++ r0) {
		(*ve)[r0] = NULL;
		(*vec)[r0] = 0;
	
		for (r1 = 0; r1 < recc; ++ r1)
			if (hc_rectinter(r + r0, r + r1) && r0 != r1) {
				if ((*vec)[r0] == 0)	
					(*ve)[r0] = e + ec;

				e[ec++] = r1;
				++((*vec)[r0]);
			}
		}

	return 0;
}

static int hc_getcomponent(int **ve, int *vec, int vn, int *isreached)
{
	uint en;

	isreached[vn] = 1;

	for (en = 0; en < vec[vn]; ++en) {
		if (!isreached[ve[vn][en]])
			hc_getcomponent(ve, vec, ve[vn][en], isreached);
	}

	return 0;
}

static int hc_splitgraph(int **ve, int *vec, int vc,
	int **v, int **cvc, int *cc)
{
	int *isreached;
	int *isreachedg;
	int rn;
	int curvc;
	int nextv;
	int maxc;

	maxc = 1;

	if ((*v = (int *) malloc(sizeof(int) * vc)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((*cvc = (int *) malloc(sizeof(int) * maxc)) == NULL) {
		free(*v);

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((isreached = (int *) malloc(sizeof(int) * vc)) == NULL) {
		hc_safefree((void **)v);
		hc_safefree((void **)cvc);
		
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	if ((isreachedg = (int *) malloc(sizeof(int) * vc)) == NULL) {
		hc_safefree((void **)v);
		hc_safefree((void **)cvc);
		hc_safefree((void **)&isreached);

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	for (rn = 0; rn < vc; ++rn)
		isreachedg[rn] = 0;	

	nextv = 0;
	curvc = 0;
	*cc = 0;
	while (1) {
		if (nextv < vc) {
			int prevc;

			for (rn = 0; rn < vc; ++rn)
				isreached[rn] = 0;
	
			hc_getcomponent(ve, vec, nextv, isreached);
		
			prevc = curvc;
			for (rn = 0; rn < vc; ++rn)
				if (isreached[rn])
					(*v)[curvc++] = rn;
			
			(*cvc)[*cc] = curvc - prevc;
			++(*cc);

			if (*cc >= maxc) {
				int *newcvc;

				maxc *= 2;

				if ((newcvc = realloc(*cvc,
					sizeof(int) * maxc)) == NULL) {
					hc_safefree((void **)v);
					hc_safefree((void **)cvc);
					hc_safefree((void **)&isreached);

					nd_seterror(ND_ALLOCFAULT);
					return (-1);
				}

				*cvc = newcvc;
			}

			for (rn = 0; rn < vc; ++rn)
				isreachedg[rn] = isreachedg[rn]
					|| isreached[rn];
	
		
			while (isreachedg[nextv] && nextv < vc)
				++nextv;	
		}
		else
			break;
	}

	return 0;
}

int hc_conrect(const struct hc_rect *r, int recc,
	struct hc_rect **newr, int *newrc)
{
	int **ve;
	int *vec;

	int *v;
	int *cvc;
	int cc;

	int cn;
	int vn;
	int *curv;

	int res;

	if (r == NULL || r <= 0 || newr == NULL || newrc == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	if ((res = hc_recttograph(r, recc, &ve, &vec)) < 0) {
		return (-1);
	}

	if ((res = hc_splitgraph(ve, vec, recc, &v, &cvc, &cc)) < 0) {
		return (-1);
	}

	if ((*newr = (struct hc_rect *) malloc(sizeof(struct hc_rect) * cc))
		== NULL) {
		nd_error = ND_ALLOCFAULT;
		return (-1);
	}

	*newrc = cc;
	curv = v;
	for (cn = 0; cn < *newrc; ++cn) {
		double sx0, sy0;
		double sx1, sy1;

		sx0 = sy0 = sx1 = sy1 = 0.0;

		for (vn = 0; vn < cvc[cn]; ++vn) {
			sx0 += r[*curv].x0;
			sy0 += r[*curv].y0;
			sx1 += r[*curv].x1;
			sy1 += r[*curv].y1;
	
			curv++;
		}

		(*newr)[cn].x0 = sx0 / (double) cvc[cn];
		(*newr)[cn].y0 = sy0 / (double) cvc[cn];
		(*newr)[cn].x1 = sx1 / (double) cvc[cn];
		(*newr)[cn].y1 = sy1 / (double) cvc[cn];
	}

	return 0;
}
