#ifndef NEURONET_H
#define NEURONET_H

#include <stdlib.h>

struct nn_neuronet
{
	int levelc;
	int *nodec;
	int *weightc;
	int inputc;
	double *weight;
	double *out;
	double **vout;
};

int nn_createneuronet(struct nn_neuronet *nnet, int inputc,
	int levelc, const int *neuronc);

int nn_neuroneteval(const struct nn_neuronet *nnet,
	double *input, double *output);

int nn_teachneuronet(const struct nn_neuronet *nnet, double *input,
	double *target, double lrate);

int nn_printneuronet(const struct nn_neuronet *nnet);

int nn_neuronettofile(const struct nn_neuronet *nnet, const char *fname);

int nn_neuronetfromfile(struct nn_neuronet *nnet, const char *fname);

#endif
