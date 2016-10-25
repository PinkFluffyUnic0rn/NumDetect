#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>

#include "hc_hcascade.h"
#include "nd_error.h"

struct hc_findwcthreadarg {
	struct hc_wclassifier *wc;
	int wcc;
	const struct nd_image *feature;
	const struct hc_trainingset *ts;
	double *w;
};	

struct hc_findwcthreadret {
	struct hc_wclassifier minerrwc;
	int *iserror;
	double minwcerr;
};	

struct hprimvalue {
	double val;
	int imgn;
};

static void hc_safefree(void **p)
{
	if (*p != NULL) {
		free(*p);
		*p = NULL;
	}
}

static int hc_strtoui(const char *cur, char **next)
{
	long int v;
	
	v = strtol(cur, next, 0);

	if (cur == *next)
		return (-1);

	if (v == LONG_MIN || v == LONG_MAX)
		return (-1);
	
	if (v > INT_MAX)
		return (-1);

	return v;
}

int hc_tsisvalid(struct hc_trainingset *ts)
{
	if (ts->img == NULL || ts->isgood == NULL || ts->imgc <= 0)
		return 0;

	return 1;
}

int hc_readtrset(struct hc_trainingset *ts, const char *hcpath)
{
	FILE *file;
	char *str;
	size_t strsz;
	char *cur;
	int goodc;
	int badc;
	int imgn;

	if (ts == NULL || hcpath == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	if ((file = fopen(hcpath, "r")) == NULL) {
		nd_seterror(ND_FOPENERROR);
		return (-1);
	}

	strsz = 0;
	if (getline(&str, &strsz, file) <= 0) {
		nd_seterror(ND_READFILEERROR);
		return (-1);
	}

	cur = str;
	
	if ((goodc = hc_strtoui(cur, &cur)) < 0) {
		hc_safefree((void **)&str);
		
		nd_seterror(ND_READFILEERROR);
		return (-1);
	}

	if ((badc = hc_strtoui(cur, &cur)) < 0) {
		hc_safefree((void **)&str);
		
		nd_seterror(ND_READFILEERROR);
		return (-1);
	}

	ts->imgc = goodc + badc;

	hc_safefree((void **)&str);
	
	if ((ts->img = malloc(sizeof(struct nd_image) * ts->imgc)) == NULL) {
		hc_safefree((void **)&str);
		
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	for (imgn = 0; imgn < ts->imgc; ++imgn) {
		if (getline(&str, &strsz, file) <= 0) {
			int imgnn;

			for (imgnn = 0; imgnn < imgn; ++imgnn)
				nd_imgdestroy(ts->img + imgnn);

			hc_safefree((void **)&(ts->img));
		
			nd_seterror(ND_READFILEERROR);
			return (-1);
		}
	
		str[strlen(str) - 1] = '\0';
			
		if (nd_imgread(str, ts->img + imgn) < 0) {
			int imgnn;
			
			for (imgnn = 0; imgnn < imgn; ++imgnn)
				nd_imgdestroy(ts->img + imgnn);
	
			hc_safefree((void **)&(ts->img));
			hc_safefree((void **)&str);
			
			return (-1);
		}

		if (nd_imggrayscale(ts->img + imgn) < 0) {
			int imgnn;
			
			for (imgnn = 0; imgnn < imgn; ++imgnn)
				nd_imgdestroy(ts->img + imgnn);

			hc_safefree((void **)&(ts->img));
			hc_safefree((void **)&str);

			return (-1);
		}

		hc_safefree((void **)&str);
	}

	if ((ts->isgood = malloc(sizeof(int) * ts->imgc)) == NULL) {
		int imgnn;
		
		for (imgnn = 0; imgnn < imgn; ++imgnn)
			nd_imgdestroy(ts->img + imgnn);

		hc_safefree((void **)&(ts->img));

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	for (imgn = 0; imgn < ts->imgc; ++imgn)
		ts->isgood[imgn] = imgn < goodc ? 1 : 0;

	return 0;
}

int hc_create(struct hc_hcascade *hc, int ww, int wh, 
	const struct nd_image *feature, int featurec)
{
	if (hc == NULL || ww <= 0 || wh <= 0
		|| feature == NULL || featurec <= 0) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	hc->wh = wh;
	hc->ww = ww;

	if ((hc->feature = malloc(sizeof(struct nd_image) * featurec))
		== NULL) {	
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	hc->featurec = featurec;
	memcpy(hc->feature, feature, sizeof(struct nd_image) * featurec);

	hc->wc = NULL;
	hc->wccoef = NULL;
	hc->wccount = 0;
	hc->stagecount = 0;

	return 0;
}

static int hc_hcisvalid(const struct hc_hcascade *hc)
{
	if (hc->wh <= 0 || hc->ww <= 0)
		return 0;

	if (hc->featurec < 0 || hc->wccount < 0 || hc->stagecount < 0)
		return 0;
	
	if ((hc->featurec > 0 && hc->feature == NULL)
		|| (hc->wccount > 0 && (hc->wc == NULL || hc->wccoef == NULL))
		|| (hc->stagecount > 0 && hc->stage == NULL))
		return 0;

	if ((hc->stagecount == 0 || hc->featurec == 0) && hc->wccount > 0)
		return 0;

	return 1;
}

static int hc_readfeatures(FILE *file, struct hc_hcascade *hc)
{
	int fn;

	if ((hc->feature = malloc(sizeof(struct nd_image) * hc->featurec))
		== NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	for (fn = 0; fn < hc->featurec; ++fn) {
		char *str;
		size_t strsz;
		char *cur;
		int w, h;
		int pixn;

		strsz = 0;
		if (getline(&str, &strsz, file) <= 0) {
			hc_safefree((void **)&(hc->feature));
			
			nd_seterror(ND_READFILEERROR);
			return (-1);
		}
	
		cur = str;

		if ((w = hc_strtoui(cur, &cur)) < 0) {
			free(hc->feature);
			hc->feature = NULL;
			hc_safefree((void **)&(hc->feature));
			
			nd_seterror(ND_READFILEERROR);
			return (-1);
		}

		if ((h = hc_strtoui(cur, &cur)) < 0) {
			hc_safefree((void **)&(hc->feature));

			nd_seterror(ND_READFILEERROR);
			return (-1);
		}
 
		if (nd_imgcreate(hc->feature + fn, w, h, 1) < 0) {
			hc_safefree((void **)&(hc->feature));

			return (-1);
		}
		
		for (pixn = 0; pixn < w * h; ++pixn) {
			char *prev;

			prev = cur;
			hc->feature[fn].data[pixn] = strtod(cur, &cur);
			
			if (prev == cur) {
				hc_safefree((void **)&(hc->feature));	
				nd_imgdestroy(hc->feature + fn);
				
				nd_seterror(ND_READFILEERROR);
				return (-1);
			}
		}

		free(str);
		str = NULL;
	}
	
	return 0;
}

static int hc_readwcs(FILE *file, struct hc_hcascade *hc)
{
	int wcn;

	if ((hc->wc = malloc(sizeof(struct hc_wclassifier) * hc->wccount))
		== NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	if ((hc->wc != 0 && hc->wccoef = malloc(sizeof(double) * hc->wccount))
		== NULL) {
		hc_safefree((void **)&(hc->wc));
		
		nd_seterror(ND_ALLOCFAULT);
		return (-1);	
	}

	for (wcn = 0; wcn < hc->wccount; ++wcn)
		if (fscanf(file, "%d %d %d %d %d %lf %lf\n", &(hc->wc[wcn].fn),
			&(hc->wc[wcn].x), &(hc->wc[wcn].y),
			&(hc->wc[wcn].w), &(hc->wc[wcn].h),
			&(hc->wc[wcn].ineqdir), &(hc->wc[wcn].thres)) == EOF) {
			hc_safefree((void **)&(hc->wc));
			hc_safefree((void **)&(hc->wccoef));
		
			nd_seterror(ND_ALLOCFAULT);
			return (-1);	
		}

	for (wcn = 0; wcn < hc->wccount; ++wcn)
		fscanf(file, "%lf", hc->wccoef + wcn);

	return 0;
}


static int hc_readstages(FILE *file, struct hc_hcascade *hc)
{
	int stn;

	if ((hc->stage = malloc(sizeof(struct hc_stage) * hc->stagecount))
		== NULL) {	
		nd_seterror(ND_ALLOCFAULT);
		return (-1);	
	}

	for (stn = 0; stn < hc->stagecount; ++stn)
		if (fscanf(file, "%d %lf", &(hc->stage[stn].wcc),
			&(hc->stage[stn].thres)) == EOF) {
				hc_safefree((void **)&(hc->stage));
				
				nd_seterror(ND_READFILEERROR);
				return (-1);
		}

	return 0;
}

int hc_hcascaderead(struct hc_hcascade *hc, const char *hcpath)
{
	FILE *file;

	if (hc == NULL || hcpath == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	if ((file = fopen(hcpath, "r")) == NULL) {
		nd_seterror(ND_FOPENERROR);
		return (-1);
	}

	fscanf(file, "%d %d %d %d %d\n", &(hc->ww), &(hc->wh), &(hc->featurec),
		&(hc->wccount), &(hc->stagecount));

	if (hc->featurec > 0) {
		if (hc_readfeatures(file, hc) < 0)
			return (-1);
	}
	else
		hc->feature = NULL;

	if (hc->wccount > 0) {
		if (hc_readwcs(file, hc) < 0) {
			hc_safefree((void **)&(hc->feature));

			return (-1);
		}
	}
	else {
		hc->wc = NULL;
		hc->wccoef = NULL;
	}
	
	if (hc->stagecount > 0) {
		if (hc_readstages(file, hc) < 0) {
			hc_safefree((void **)&(hc->feature));
			hc_safefree((void **)&(hc->wc));
	
			return (-1);
		}
	}
	else
		hc->stage = NULL;

	if (fclose(file) == EOF) {
		nd_seterror(ND_FCLOSEERROR);
		return (-1);
	}

	return 0;
}

int hc_hcascadewrite(struct hc_hcascade *hc, const char *hcpath)
{
	FILE *file;
	int fn;
	int wcn;
	int stn;
	
	if (hc == NULL || hcpath == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	if (!hc_hcisvalid(hc)) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	if ((file = fopen(hcpath, "w")) == NULL) {
		nd_seterror(ND_FOPENERROR);
		return (-1);
	}

	if (fprintf(file, "%d %d %d %d %d\n", hc->ww, hc->wh,
		hc->featurec, hc->wccount, hc->stagecount) < 0) {
		nd_seterror(ND_WRITETOFILEFERROR);
		return (-1);
	}

	for (fn = 0; fn < hc->featurec; ++fn) {
		int pixn;
		int pixc;

		if (fprintf(file, "%d %d ",
			hc->feature[fn].w, hc->feature[fn].h) < 0) {
			nd_seterror(ND_WRITETOFILEFERROR);
			return (-1);
		}

		pixc = hc->feature[fn].w * hc->feature[fn].h;

		for (pixn = 0; pixn < pixc; ++pixn)
			if (fprintf(file, "%lf%c", hc->feature[fn].data[pixn],
				(pixn != pixc - 1) ? ' ' : '\n') < 0) {
				nd_seterror(ND_WRITETOFILEFERROR);
				return (-1);
			}
	}

	for (wcn = 0; wcn < hc->wccount; ++wcn)
		if (fprintf(file, "%d %d %d %d %d %lf %lf\n", hc->wc[wcn].fn,
			hc->wc[wcn].x, hc->wc[wcn].y,
			hc->wc[wcn].w, hc->wc[wcn].h, 
			hc->wc[wcn].ineqdir, hc->wc[wcn].thres) < 0) {
				nd_seterror(ND_WRITETOFILEFERROR);
				return (-1);
			}

	for (wcn = 0; wcn < hc->wccount; ++wcn)
		if (fprintf(file, "%lf%c", hc->wccoef[wcn],
			(wcn != hc->wccount - 1) ? ' ' : '\n') < 0) {
				nd_seterror(ND_WRITETOFILEFERROR);
				return (-1);
			}

	for (stn = 0; stn < hc->stagecount; ++stn)
		if (fprintf(file, "%d %lf\n",
			hc->stage[stn].wcc, hc->stage[stn].thres) < 0) {
				nd_seterror(ND_WRITETOFILEFERROR);
				return (-1);
			}

	if (fclose(file) == EOF) {
		nd_seterror(ND_FCLOSEERROR);
		return (-1);
	}

	return 0;
}

static int hc_next_wclassifier(struct hc_wclassifier *wc,
	struct nd_image *feature, int featurec, int wh, int ww)
{
	if (wc->fn == -1) {
		wc->fn = 0;

		wc->x = 0;
		wc->y = 0;

		wc->h = feature[wc->fn].h;
		wc->w = feature[wc->fn].w;
	}
	else if (wc->x + wc->w  < ww)
		++wc->x;
	else if (wc->y + wc->h < wh) {
		wc->x = 0;
		++wc->y;
	}
	else if (wc->h < wh) {
		wc->x = 0;
		wc->y = 0;
		
		++wc->h;
	}
	else if (wc->w < ww) {
		wc->h = feature[wc->fn].h;
	
		wc->x = 0;
		wc->y = 0;
	
		++wc->w;
	}
	else if (wc->fn + 1 < featurec) {
		++wc->fn;

		wc->x = 0;
		wc->y = 0;
	
		wc->h = feature[wc->fn].h;
		wc->w = feature[wc->fn].w;
	}
	else
		return 0;

	return 1;
}

static inline double hc_gethprimval(const struct hc_wclassifier *wc,
	const struct nd_image *feat, const struct nd_image *img)
{
	int fx, fy;
	double s;
	double wrel;
	double hrel;
	int ywoffset;
	double relfxoffset;
	double relfyoffset;
	
	wrel = (double) wc->w / (double) feat->w;
	hrel = (double) wc->h / (double) feat->h;

	s = 0.0;

	ywoffset = 0;
	relfyoffset = 0.0;

	for (fy = 0; fy < feat->h; ++fy) {	
		relfxoffset = 0.0;
		for (fx = 0; fx < feat->w; ++fx) {
			int imgx0;
			int imgy0;
			int imgx1;
			int imgy1;

			double sx1y0;
			double sx0y1;
			double sx0y0;;
			double sx1y1;

			imgx0 = (int) (relfxoffset) + wc->x - 1;
			imgy0 = img->w
				* ((int) (relfyoffset) + wc->y - 1);	

			imgx1 = (int) (relfxoffset + wrel) + wc->x - 1;
			imgy1 = img->w
				* ((int) (relfyoffset + hrel) + wc->y - 1);

			
			sx1y1 = img->data[imgy1 + imgx1];

			if (imgx0 >= 0 && imgy0 >= 0) {
				sx1y0 = img->data[imgy0 + imgx1];
				sx0y1 = img->data[imgy1 + imgx0];
				sx0y0 = img->data[imgy0 + imgx0];
			}
			else {	
				sx1y0 = (imgy0 >= 0)
					? img->data[imgy0 + imgx1]
					: 0.0;

				sx0y1 = (imgx0 >= 0)
					? img->data[imgy1 + imgx0]
					: 0.0;
			
				sx0y0 = 0.0;
			}

			s += feat->data[ywoffset + fx]
				* (sx0y0 + sx1y1 - sx1y0 - sx0y1);
		
			relfxoffset += wrel;
		}

		ywoffset += feat->w;
		relfyoffset += hrel;
	}
	
	return s;
}

static int hc_wclassify(const struct hc_wclassifier *wc,
	const struct nd_image *feat, const struct nd_image *img)
{
	return ((wc->ineqdir * hc_gethprimval(wc, feat, img))
		< (wc->ineqdir * wc->thres));
}

int hc_imgintegral(struct nd_image *img)
{
	double *c;
	int x;
	int y;
	int cn;

	c = malloc(sizeof(double) * img->w * img->h);

	for (cn = 0; cn < img->h * img->w; ++cn)
		c[cn] = 0.0;

	for (y = 0; y < img->h; ++y) {
		for (x = 0; x < img->w; ++x) {
			double sum;
			double t, a;
			
			sum = 0.0;
			a = img->data[y * img->w + x];
			
			if (y > 0) {
				sum += img->data[(y - 1) * img->w + x];
				a -= c[(y - 1) * img->w + x];
			}

			if (x > 0) {
				sum += img->data[y * img->w + (x - 1)];
				a -= c[y * img->w + (x - 1)];
			}

			if (x > 0 && y > 0) {
				sum -= img->data[(y - 1) * img->w + (x - 1)];	
				a += c[(y - 1) * img->w + (x - 1)];
			}

			t = sum + a;
			c[y * img->w + x] = (t - sum) - a;
			sum = t;
			img->data[y * img->w + x] = sum;	
		}
	}

	free(c);

	return 0;
}

static int hc_initweights(double *w, struct hc_trainingset *ts)
{
	int imgn;
	int goodec;
	int badec;

	goodec = badec = 0;	
	for (imgn = 0; imgn < ts->imgc; ++imgn)
		if (ts->isgood[imgn])
			++goodec;
		else
			++badec;

	for (imgn = 0; imgn < ts->imgc; ++imgn)
		w[imgn] = 1.0
			/ (2.0 * (double) (ts->isgood[imgn] ? goodec : badec));
	
	return 0;
}

static int comphval(const void *hv0, const void *hv1)
{
	return ((((struct hprimvalue *) hv0)->val
		> ((struct hprimvalue *) hv1)->val) ? 1 : -1);

}

static double hc_trainwclassifier(struct hc_wclassifier *curwc,
	const struct nd_image *feature, const struct hc_trainingset *ts,
	const double *w)
{
	int imgn;
	struct hprimvalue *hpval;
	double minerr;
	double minerrthres;
	double minerrdir;
	double thres;

	int curimgn;
	double *goodwi;
	double *badwi;

	if ((hpval = malloc(sizeof(struct hprimvalue) * ts->imgc)) == NULL)
		return (-1.0);
		
	if ((goodwi = malloc(sizeof(double) * ts->imgc)) == NULL)
		return (-1.0);

	if ((badwi = malloc(sizeof(double) * ts->imgc)) == NULL)
		return (-1.0);

	for (imgn = 0; imgn < ts->imgc; ++imgn) {
		hpval[imgn].val = hc_gethprimval(curwc,
			feature + curwc->fn, ts->img + imgn);
		hpval[imgn].imgn = imgn;
	}

	qsort(hpval, ts->imgc, sizeof(struct hprimvalue), comphval);

	curimgn = hpval[0].imgn;
	
	goodwi[0] = (ts->isgood[curimgn]) ? w[curimgn] : 0.0;
	badwi[0] = (ts->isgood[curimgn]) ? 0.0 : w[curimgn];

	for (imgn = 1; imgn < ts->imgc; ++imgn) {
		curimgn = hpval[imgn].imgn;

		goodwi[imgn] = goodwi[imgn - 1]
			+ ((ts->isgood[curimgn]) ? w[curimgn] : 0.0);
		badwi[imgn] = badwi[imgn - 1]
			+ ((ts->isgood[curimgn]) ? 0.0 : w[curimgn]);
	}


	minerr = 1.0;
	minerrthres = 0.0;
	minerrdir = 1.0;

	for (imgn = 1; imgn < ts->imgc; ++imgn) {
		double err;

		thres = (hpval[imgn - 1].val + hpval[imgn].val) * 0.5;

		err = badwi[imgn - 1]
			+ goodwi[ts->imgc - 1] - goodwi[imgn - 1];

		if (err < minerr) {
			minerr = err;
			minerrthres = thres;
			minerrdir = 1.0;
		}

		err = goodwi[imgn - 1]
			+ badwi[ts->imgc - 1] - badwi[imgn - 1];

		if (err < minerr) {
			minerr = err;
			minerrthres = thres;
			minerrdir = -1.0;
		}
	}
	
	hc_safefree((void **) &goodwi);
	hc_safefree((void **) &badwi);
	hc_safefree((void **) &hpval);

	curwc->thres = minerrthres;
	curwc->ineqdir = minerrdir;

	return minerr;
}


static void *findwcthreadfunc(void *a)
{
	int wcn;
	double minwcerr;
	struct hc_wclassifier *minerrwc;
	struct hc_findwcthreadarg *arg;
	struct hc_findwcthreadret *ret;
	
	arg = (struct hc_findwcthreadarg *) a;

	minwcerr = 1.0;
	minerrwc = NULL;
	for (wcn = 0; wcn < arg->wcc; ++wcn) {
		double wcerror;

		if ((wcerror = hc_trainwclassifier(arg->wc + wcn, arg->feature,
			arg->ts, arg->w)) < 0.0)
			pthread_exit(NULL);
		
		if ((wcerror < minwcerr) || (wcn == 0)) {
			minwcerr = wcerror;
			minerrwc = arg->wc + wcn;
		}
	}

	if ((ret = malloc(sizeof(struct hc_findwcthreadret))) == NULL)
		pthread_exit(NULL);

	ret->minerrwc = *minerrwc;
	ret->minwcerr = minwcerr;

	pthread_exit(ret);
}

static int hc_findbestwc(const struct hc_hcascade *hc,
	const struct hc_trainingset *ts, double *weight,
	struct hc_wclassifier *minerrwc, double *minwcerr, int threadwcc,
	int tc)
{
	struct hc_wclassifier curwc;
	int islastpart;
	struct hc_wclassifier *wc;
	pthread_t *thrs;
	struct hc_findwcthreadarg *threadarg;
	struct hc_findwcthreadret **threadret;	

	if ((wc = malloc(sizeof(struct hc_wclassifier) * threadwcc * tc))
		== NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}
	
	if ((thrs = malloc(sizeof(pthread_t) * tc)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		hc_safefree((void **)&wc);

		return (-1);
	}
	
	if ((threadarg = (struct hc_findwcthreadarg *)
		malloc(sizeof(struct hc_findwcthreadarg) * tc)) == NULL) {	
		hc_safefree((void **)&wc);
		hc_safefree((void **)&thrs);
		
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((threadret = (struct hc_findwcthreadret **)
		malloc(sizeof(struct hc_findwcthreadret *) * tc)) == NULL) {
		hc_safefree((void **)&wc);
		hc_safefree((void **)&thrs);
		hc_safefree((void **)&threadarg);

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	curwc.fn = -1;
	curwc.w = 0;
	curwc.h = 0;
	curwc.x = 0;
	curwc.y = 0;
	curwc.ineqdir = 1.0;
	curwc.thres = 0.0;

	*minwcerr = 1.0;
	*minerrwc = curwc;

	islastpart = 0;

	while(!islastpart) {
		int tn;
		int wcn;
		int partwcc;
	
		partwcc = 0;
		for (wcn = 0; wcn < (threadwcc * tc); ++wcn) {
			int res;
			
			if ((res = hc_next_wclassifier(&curwc, hc->feature,
				hc->featurec, hc->wh, hc->ww)) == 0) {
				islastpart = 1;
				break;
			}

			++partwcc;
		
			wc[wcn] = curwc;
		}

		threadwcc = partwcc / tc;

		for (tn = 0; tn < tc; ++tn) {	
			threadarg[tn].wc = wc + tn * threadwcc;
			threadarg[tn].wcc = (tn != (tc - 1))
				? threadwcc : (partwcc - tn * threadwcc);
			threadarg[tn].feature = hc->feature;
			threadarg[tn].ts = ts;
			threadarg[tn].w = weight;
	
			pthread_create(thrs + tn, NULL, &findwcthreadfunc,
				(void *)(threadarg + tn));
		}
	
		for (tn = 0; tn < tc; ++tn)
			pthread_join(thrs[tn], (void **) (threadret + tn));

		for (tn = 0; tn < tc; ++tn)
			if (threadret[tn] == NULL) {
				int tnn;
				
				for (tnn = 0; tnn < tc; ++tnn)
					hc_safefree((void **)
						(threadret + tnn));
	
				hc_safefree((void **)&wc);
				hc_safefree((void **)&thrs);
				hc_safefree((void **)&threadarg);
				hc_safefree((void **)&threadret);			
		
				return (-1);
			}
			else if (threadret[tn]->minwcerr < *minwcerr) {
				*minwcerr = threadret[tn]->minwcerr;
				*minerrwc = threadret[tn]->minerrwc;
			}

		for (tn = 0; tn < tc; ++tn)
			hc_safefree((void **)(threadret + tn));
	}

	hc_safefree((void **)&wc);
	hc_safefree((void **)&thrs);
	hc_safefree((void **)&threadarg);
	hc_safefree((void **)&threadret);

	return 0;
}

int hc_findwc(struct hc_hcascade *hc, struct hc_trainingset *ts,
	double **weights, int iterc)
{
	double *w;
	int wn;
	int imgn;
	int i;
	int *iserror;
	int wcn;

	if (hc == NULL || ts == NULL || weights == NULL || iterc == 0) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	if (!hc_hcisvalid(hc)) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	if (!hc_tsisvalid(ts)) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	for (imgn = 0; imgn < ts->imgc; ++imgn) {
		if (nd_imgnormalize(ts->img + imgn, 0, 1) < 0)
			return (-1);

		if (hc_imgintegral(ts->img + imgn) < 0)
			return (-1);
	}

	if ((iserror = malloc(sizeof(int) * ts->imgc)) == NULL) {
		nd_error = ND_ALLOCFAULT;
		return (-1);
	}

	if (*weights == NULL) {
		if ((w = malloc(sizeof(double) * ts->imgc)) == NULL) {
			hc_safefree((void **)&iserror);

			nd_error = ND_ALLOCFAULT;
			return (-1);
		}
	
		hc_initweights(w, ts);
	}
	else
		w = *weights;
	
	for (i = 0; i < iterc; ++i) {
		double sumw;
		double minwcerr;
		struct hc_wclassifier minerrwc;

		struct hc_wclassifier *newwc;
		double *newwccoef;

		int j;

		sumw = 0;
		for (wn = 0; wn < ts->imgc; ++wn)
			sumw += w[wn];

		for (wn = 0; wn < ts->imgc; ++wn)
			w[wn] /= sumw;
		
		if (hc_findbestwc(hc, ts, w,
			&minerrwc, &minwcerr, 500, 8) < 0) {
			hc_safefree((void **)&iserror);
			
			return (-1);
		}

		for (imgn = 0; imgn < ts->imgc; ++imgn) {
			iserror[imgn] = abs(ts->isgood[imgn]
				- hc_wclassify(&minerrwc, 
					hc->feature + minerrwc.fn,
					ts->img + imgn));
		}

		double beta = minwcerr / (1.0 - minwcerr);

		for (j = 0; j < ts->imgc; ++j)
			w[j] *= pow(beta, 1 - iserror[j]);

	
		++hc->wccount;

		if ((newwc = realloc(hc->wc, sizeof(struct hc_wclassifier)
			* hc->wccount)) == NULL
			|| (newwccoef = realloc(hc->wccoef, sizeof(double)
			* hc->wccount)) == NULL) {
			hc_safefree((void **)&iserror);
			hc_safefree((void **)&(hc->wc));
			hc_safefree((void **)&(hc->wccoef));

			--hc->wccount;
		
			nd_error = ND_ALLOCFAULT;
			return (-1);
		}
		
		hc->wc = newwc;

		hc->wccoef = newwccoef;

		hc->wc[hc->wccount - 1] = minerrwc;
		hc->wccoef[hc->wccount - 1] = log(1.0 / beta);
	}
	
	*weights = w;

	if (hc->stage != NULL) {
		free(hc->stage);
		hc->stage = NULL;
	}

	hc->stagecount = 1;
	
	if ((hc->stage = malloc(sizeof(struct hc_stage))) == NULL) {
		hc_safefree((void **)&iserror);
		hc_safefree((void **)&(hc->wc));
		hc_safefree((void **)&(hc->wccoef));
	}

	hc->stage[0].wcc = hc->wccount;

	hc->stage[0].thres = 0.0;
	for (wcn = 0; wcn < hc->wccount; ++wcn)
		hc->stage[0].thres += hc->wccoef[wcn];

	return 0;
}

static int hc_stageclassify(const struct hc_wclassifier *wc, const double *wccoef,
	const struct hc_stage *stage, const struct nd_image *feature,
	const struct nd_image *img)
{
	int wcn;
	double res;

	res = 0.0;
	for (wcn = 0; wcn < stage->wcc; ++wcn)
		if (hc_wclassify(wc + wcn, feature + wc[wcn].fn, img))
			res += wccoef[wcn];
	
	return ((res >= stage->thres) ? 1 : 0);
}

static int hc_imgtest(const struct hc_wclassifier *wc,
	const double *wccoef, const struct hc_stage *stage, int stagec,
	const struct nd_image *feature, const struct nd_image *img)
{
	int sn;

	int res = 1;

	for (sn = 0; sn < stagec; ++sn) {
		if (hc_stageclassify(wc, wccoef, stage + sn,
			feature, img) == 0) {
			res = 0;
			break;
		}
	}

	return res;
}

const double hc_stagefprate(struct hc_wclassifier *wc, double *wccoef,
	struct hc_stage *stage, int curstage, struct hc_trainingset *ts,
	struct nd_image *feature, int *fpcount)
{
	int fncount;
	int imgn;
	int laststimgc;

	fncount = 0;
	(*fpcount) = 0;
	laststimgc = 0;

	for (imgn = 0; imgn < ts->imgc; ++imgn) {
		int res;
	
		res = hc_imgtest(wc, wccoef, stage, curstage,
			feature, ts->img + imgn);

		if (res) {
			res = hc_stageclassify(wc, wccoef,
				stage + curstage, feature, ts->img + imgn);

			if (ts->isgood[imgn] == 0)
				++laststimgc;

			res = ts->isgood[imgn] - res;
			
			if (res > 0)
				++fncount;
			if (res < 0)
				++(*fpcount);
		}
	}

	return (double) (*fpcount) / (double) laststimgc;
}

static double hc_getstagethres(double *hvsum, int goodimgc, double tprate)
{
	int thresid;

	thresid = (int) ((1.0 - tprate) * goodimgc);
	
	return hvsum[thresid] + 0.00001;
}

static int sorthvsum(const void *hvsum0, const void *hvsum1)
{
	return (*((double *)hvsum0) > *((double *)hvsum1)) ? 1 : -1;
}

int hc_buildcascade(struct hc_hcascade *hc, struct hc_trainingset *ts,
	double fprtarget, double fntarget)
{
	int imgn;
	int stagemaxc;
	int fpcount;
	int wcn;
	double *hvsum;
	double *goodhvsum;
	
	if (hc == NULL || ts == NULL) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	if (!hc_hcisvalid(hc)) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	if (!hc_tsisvalid(ts)) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}

	if (fprtarget < 0.0 || fprtarget > 1.0) {
		nd_seterror(ND_INVALIDARG);
		return (-1);
	}
	
	stagemaxc = 1;
	hc->stagecount = 0;

	if ((hc->stage = (struct hc_stage *) malloc(sizeof(struct hc_stage)
		* stagemaxc)) == NULL) {
		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((hvsum = (double *) malloc(sizeof(double) * ts->imgc)) == NULL) {
		hc_safefree((void **)&(hc->stage));

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}

	if ((goodhvsum = (double *) malloc(sizeof(double) * ts->imgc))
		== NULL) {
		hc_safefree((void **)&(hc->stage));
		hc_safefree((void **)&(hvsum));

		nd_seterror(ND_ALLOCFAULT);
		return (-1);
	}


	for (imgn = 0; imgn < ts->imgc; ++imgn)
		hvsum[imgn] = 0.0;

	for (imgn = 0; imgn < ts->imgc; ++imgn) {
		if (nd_imgnormalize(ts->img + imgn, 0, 1) < 0) {
			hc_safefree((void **)&(hc->stage));
			hc_safefree((void **)&hvsum);

			return (-1);
		}

		if (hc_imgintegral(ts->img + imgn) < 0) {
			hc_safefree((void **)&(hc->stage));
			hc_safefree((void **)&hvsum);

			return (-1);
		}
	}
	
	hc->stagecount = 0;
	hc->stage[hc->stagecount].wcc = 1;
	do {	
		double fprate;
		int wcc;
		int goodimgn;
	
		wcc = hc->stage[hc->stagecount].wcc;
		
		for (imgn = 0; imgn < ts->imgc; ++imgn)
			if (hc_wclassify(hc->wc + wcc,
				hc->feature + hc->wc[wcc].fn, ts->img + imgn))
				hvsum[imgn] += hc->wccoef[wcc];

		goodimgn = 0;
		for (imgn = 0; imgn < ts->imgc; ++imgn)
			if (ts->isgood[imgn])
				goodhvsum[goodimgn++] = hvsum[imgn];

		qsort(goodhvsum, goodimgn, sizeof(double), sorthvsum);

		hc->stage[hc->stagecount].thres = hc_getstagethres(goodhvsum,
			goodimgn, fntarget);

		fprate = hc_stagefprate(hc->wc, hc->wccoef, hc->stage,
			hc->stagecount, ts, hc->feature,
			&fpcount);

		if (fprate > fprtarget)
			++(hc->stage[hc->stagecount].wcc);	
		else {
			++hc->stagecount;

			if (hc->stagecount >= stagemaxc) {
				struct hc_stage *newstage;

				stagemaxc *= 2;
				
				if ((newstage = (struct hc_stage *)
					realloc(hc->stage,
					sizeof(struct hc_stage)
					* stagemaxc)) == NULL) {
					hc_safefree((void **)&(hc->stage));
					hc_safefree((void **)&hvsum);
					
					nd_seterror(ND_ALLOCFAULT);
					return (-1);
				}
			
				hc->stage = newstage;
			}

			hc->stage[hc->stagecount].wcc
				= hc->stage[hc->stagecount - 1].wcc + 1;
		}

		if (fpcount == 0)
			break;

		if (hc->stage[hc->stagecount].wcc >= hc->wccount) {
			++hc->stagecount;
			break;			
		}

	} while (1);

	if (hc->stage[hc->stagecount - 1].wcc < hc->wccount) {	
		++hc->stagecount;
		
		if (hc->stagecount >= stagemaxc) {
			struct hc_stage *newstage;
			
			stagemaxc *= 2;

			if ((newstage = (struct hc_stage *) realloc(hc->stage,
				sizeof(struct hc_stage) * stagemaxc))
				== NULL) {
				hc_safefree((void **)&(hc->stage));
				hc_safefree((void **)&hvsum);
				
				nd_seterror(ND_ALLOCFAULT);
				return (-1);
			}
		
			hc->stage = newstage;
		}

		hc->stage[hc->stagecount - 1].wcc = hc->wccount;	
	}
	
	hc->stage[hc->stagecount - 1].thres = 0.0;

	for (wcn = 0; wcn < hc->stage[hc->stagecount - 1].wcc; ++wcn)
		hc->stage[hc->stagecount - 1].thres += hc->wccoef[wcn];

	hc->stage[hc->stagecount - 1].thres *= 0.5;

	return 0;
}

int hc_imgclassify(const struct hc_hcascade *hc, const struct nd_image *img)
{
	int sn;
	int res;

	res = 1;
	for (sn = 0; sn < hc->stagecount; ++sn)
		if (hc_stageclassify(hc->wc, hc->wccoef, hc->stage + sn, 
			hc->feature, img) == 0) {
			res = 0;
			break;
		}

	return res;
}

static inline double hc_fastgethprimval(const struct hc_wclassifier *wc,
	const struct nd_image *feat, const struct nd_image *img,
	int imgstride, double sd)
{
	int fx, fy;
	double s;
	double wrel;
	double hrel;
	int ywoffset;
	double relfxoffset;
	double relfyoffset;
	
	wrel = (double) wc->w / (double) feat->w;
	hrel = (double) wc->h / (double) feat->h;

	s = 0.0;

	ywoffset = 0;
	relfyoffset = 0.0;

	for (fy = 0; fy < feat->h; ++fy) {	
		relfxoffset = 0.0;
		for (fx = 0; fx < feat->w; ++fx) {
			int imgx0;
			int imgy0;
			int imgx1;
			int imgy1;

			double sx1y0;
			double sx0y1;
			double sx0y0;;
			double sx1y1;

			imgx0 = (int) (relfxoffset) + wc->x;
			imgy0 = imgstride
				* ((int) (relfyoffset) + wc->y);	

			imgx1 = (int) (relfxoffset + wrel) + wc->x;
			imgy1 = imgstride
				* ((int) (relfyoffset + hrel) + wc->y);

			
			sx1y1 = img->data[imgy1 + imgx1];
			sx1y0 = img->data[imgy0 + imgx1];
			sx0y1 = img->data[imgy1 + imgx0];
			sx0y0 = img->data[imgy0 + imgx0];

			s += feat->data[ywoffset + fx]
				* (sx0y0 + sx1y1 - sx1y0 - sx0y1) / sd;
		
			relfxoffset += wrel;
		}

		ywoffset += feat->w;
		relfyoffset += hrel;
	}
	
	return s;
}

static int hc_fastwclassify(const struct hc_wclassifier *wc,
	const struct nd_image *feat, const struct nd_image *img,
	int imgstride, double sd)
{
	return ((wc->ineqdir * hc_fastgethprimval(wc, feat, img,
		imgstride, sd)) < (wc->ineqdir * wc->thres));
}

static int hc_faststageclassify(const struct hc_wclassifier *wc,
	const double *wccoef, const struct hc_stage *stage,
	const struct nd_image *feature, const struct nd_image *img,
	int imgstride, double sd)
{
	int wcn;
	double res;

	res = 0.0;
	for (wcn = 0; wcn < stage->wcc; ++wcn)
		if (hc_fastwclassify(wc + wcn, feature + wc[wcn].fn, img,
			imgstride, sd))
			res += wccoef[wcn];
	
	return ((res >= stage->thres) ? 1 : 0);
}

int hc_fastimgclassify(const struct hc_hcascade *hc,
	const struct nd_image *img, int imgstride, double sd)
{
	int sn;
	int res;

	res = 1;
	for (sn = 0; sn < hc->stagecount; ++sn)
		if (hc_faststageclassify(hc->wc, hc->wccoef, hc->stage + sn, 
			hc->feature, img, imgstride, sd) == 0) {
			res = 0;
			break;
		}

	return res;
}
