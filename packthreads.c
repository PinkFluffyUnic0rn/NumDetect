#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIN_W 12
#define WIN_H 56 

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))

static int fillpacket(int *input, int ic, int **curp, int s, int ib, int ie)
{
	int ps;
	
	*((*curp)++) = ib;
	ps = input[ib];

	while (ps < s && ib != ie) {
		*((*curp)++) = ie;
		ps += input[ie];

		--ie;
	}

	return ie;
}

static int fillpacketarray(int *input, int ic, int **curp, int pc)
{
	int ib, ie;
	int s;
	int pn;
	int res;

	ib = 0;
	ie = ic - 1;

	s = 0.0;

	do {
		int pn;
	
		s += input[ib];
	
		pn = 0;
		do {
			ie = fillpacket(input, ic, curp + pn, s, ib, ie);

			if (ib == ie)
				break;

			++pn;
			++ib;
		} while (pn < pc);
	} while (ib != ie);

	return 0;
}

int pack(int *input, int ic, int **p, int pc)
{
	int i;

	int **curp = malloc(sizeof(int *) * pc);

	for (i = 0; i < pc; ++i)
		curp[i] = p[i] = malloc(sizeof(int) * ic);
	
	fillpacketarray(input, ic, curp, pc);

	for (i = 0; i < pc; ++i)
		*(curp[i]) = -1;
	
	return 0;
}

int main(int argc, const char **argv)
{
	int w, h;
	double sc;
	int tc;
	int *input;
	int i;
	double b;

	w = atoi(argv[1]);
	h = atoi(argv[2]);
	sc = atof(argv[3]);
	tc = atoi(argv[4]);

	b = log(sc);

	input = malloc(sizeof(int)
		* (int) (floor(log(MAX((double) WIN_W / w, (double) WIN_H / h))
		/ b) + 1));

	i = 0;
	do {	
		if (w < WIN_W || h < WIN_H)
			break;

		input[i++] = w * h;
		
		w = (double) w * sc;
		h = (double) h * sc;
	} while (1);

	int **p;

	p = malloc(sizeof(int *) * tc);
	pack(input, i, p, tc);

	int pn;
	for (pn = 0; pn < tc; ++pn) {
		for (i = 0; p[pn][i] >= 0; ++i)
			printf("%d ", p[pn][i]);
		printf("\n");
	}

	return 0;
}
