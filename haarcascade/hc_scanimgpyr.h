#include "nd_image.h"
#include "nd_error.h"
#include "hc_hcascade.h"
#include "hc_rect.h"

struct hc_scanconfig {
	double scalestep;
	int winwstep;
	int winhstep;
};

int hc_imgpyramidscan(struct hc_hcascade *hc, const struct nd_image *img,
	struct hc_rect **newr, int *newrc, const struct hc_scanconfig *conf);
