#include <nd_image.h>

int ed_findborder(struct nd_image *img, struct nd_image *models,
	int modelcount, double inpoints[8], int *bestmodel);
