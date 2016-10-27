#ifndef FFT_H
#define FFT_H

#include <stdint.h>
#include <stdlib.h>

struct complexd {
	double real;
	double imag;
};

int log2uint(uint32_t a);

struct complexd *fft(struct complexd *data, int samples_count, int inv);

int fft2d(struct complexd *data, int w, int h, int inv);

#endif
