#ifndef CONRECT_H
#define CONRECT_H

#include <stdlib.h>
#include <stdio.h>

struct hc_rect {
	uint x0;
	uint y0;
	uint x1;
	uint y1;
};

int hc_conrect(const struct hc_rect *r, int recc,
	struct hc_rect **newr, int *newrc);

#endif
