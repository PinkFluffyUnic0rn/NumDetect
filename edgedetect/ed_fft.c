#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ed_fft.h"

int log2uint(uint32_t a)
{
	int c = 0;

	while (a >>= 1) 
		++c;

	return c;
}

struct complexd *fft( struct complexd *data, int samples_count, int inv )
{
	struct complexd *out;
	int bits_count;
	int i, j, k;
		
	int lp;
	int lp2;

	out = (struct complexd *) malloc(sizeof( struct complexd)*samples_count);
	bits_count = log2( samples_count );

	for ( i = 0; i < samples_count; ++i )
	{
		int b = 0;
		int a = i;
			
		for ( j = 0; j < bits_count; ++j )
		{
			b <<= 1;
			b |= (a & 1);
			a >>= 1;
		}
	
		out[i] = data[b];
	}

	for ( k = 1; k <= bits_count; ++k )
	{
		double darg;
		double arg;

		lp = pow( 2, k );
		lp2 = lp / 2;
		darg = -inv * M_PI / lp2;
		arg = 0;

		for ( j = 0; j < lp2; ++j )
		{
			double c;
			double s;

			c = cos(arg);
			s = sin(arg);

			arg += darg;

			for ( i = j; i < samples_count; i += lp )
			{
				int iw;
				double wr, wi;

				iw = i + lp2;
		
				wr = out[iw].real * c - out[iw].imag * s;
				wi = out[iw].real * s + out[iw].imag * c;
				
				out[iw].real = out[i].real - wr;
				out[iw].imag = out[i].imag - wi;
				out[i].real = out[i].real + wr;
				out[i].imag = out[i].imag + wi;
			}
		}
	}

	if ( inv == 1 )
		for ( i = 0; i < samples_count; ++i )
		{
			out[i].real = out[i].real / samples_count;
			out[i].imag = out[i].imag / samples_count;
		}
		
	return out;
}

int fft2d(struct complexd *data, int w, int h, int inv)
{
	struct complexd *tmph
		= (struct complexd *) malloc(sizeof(struct complexd) * h );
	struct complexd *out;
	int i, j;

	for (i = 0; i < h; ++i) {	
		out = fft(data + i * w, w, inv);
		memmove(data + i*w, out, sizeof(struct complexd) * w);
		free(out);
	}
	
	for (i = 0; i < w; ++i) {
		for (j = 0; j < h; ++j)
			tmph[j] = data[j * w + i];
		
		out = fft(tmph, h, inv);

		for (j = 0; j < h; ++j)
			data[j * w + i] = out[j];

		free(out);
	}
	
	free(tmph);

	return 1;
}
