#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <dirent.h>

#include <cairo.h>

#include "nd_error.h"
#include "hc_hcascade.h"

int imgtohfeature(struct nd_image *img)
{
	int pixn;

	for (pixn = 0; pixn < img->w * img->h; ++pixn) {
		img->data[pixn] *= 2.0;
		img->data[pixn] -= 1.0;
	}

	return 0;
}

int loadfeatures(const char *dirpath, struct nd_image **f, int *fcount)
{
	DIR *d;
	struct dirent *de;
	int32_t maxfcount;

	maxfcount = 1;
	*fcount = 0;

	*f = malloc(sizeof(struct nd_image) * maxfcount);

	if ((d = opendir(dirpath)) == NULL) {
		closedir(d);
		return (-1);
	}
	
	while ((de = readdir(d)) != NULL) {
		char fullpath[255];
		struct nd_image img;
	
		sprintf(fullpath, "%s/%s", dirpath, de->d_name);
		
		if (nd_imgread(fullpath, &img) >= 0) {
			if (nd_imggrayscale(&img))
				return (-1);

			imgtohfeature(&img);

			(*f)[*fcount] = img;
			++(*fcount);
	
			if (*fcount >= maxfcount) {
				struct nd_image *newf;
	
				maxfcount *= 2;
				
				if ((newf = realloc(*f,
					sizeof(struct nd_image) * maxfcount))
					== NULL) {
					free(*f);
					return (-1);
				}
				
				*f = newf;
			}
		}
	}

	closedir(d);

	return 0;
}

int main(int argc, const char **argv)
{
	struct hc_hcascade hc;
	struct nd_image *f;
	int fcount;
	int w, h;
	char *next;


	if (argc != 5) {
		printf("%s\n", "Wrong argument count.");
		return 1;
	}

	w = strtol(argv[1], &next, 0);

	if (argv[1] == next) {
		printf("%s", "Wrong format of width.");
		return 1;
	}

	h = strtol(argv[2], &next, 0);

	if (argv[2] == next) {
		printf("%s", "Wrong format of height.");
		return 1;
	}

	loadfeatures(argv[3], &f, &fcount);

	if (hc_create(&hc, w, h, f, fcount) < 0) {
		printf("%s\n", "Cannot create cascade.");
		return 1;
	}

	if (hc_hcascadewrite(&hc, argv[4]) < 0) {
		printf("%s\n", "Cannot write cascade to file.");
		return 1;	
	}


	return 0;
}
