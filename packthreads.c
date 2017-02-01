#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIN_W 56
#define WIN_H 12 

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))

int main(int argc, const char **argv)
{
	int w, h;
	double sc;
	double a;
	int tc;
	int n;

	w = atoi(argv[1]);
	h = atoi(argv[2]);
	sc = atof(argv[3]);
	tc = atoi(argv[4]);

	n = floor(log(MAX((double) WIN_W / w, (double) WIN_H / h)) / log(sc));
	
	a = (1 - pow(sc, n)) / tc / (1 - sc);

	int b,e;
	int tn;
	
	b = 0;
	e = n;
	
	while (b < e) {
		double s;

		s = pow(sc, b);

		printf("%1.3f ", pow(sc, b));
		while(s + pow(sc, e) < a && b < e) {
			s += pow(sc, e);
			printf("%1.3f ", pow(sc, e));
			--e;
		}

		printf("\n");
		
		++b;
	}
	
	return 0;
}
