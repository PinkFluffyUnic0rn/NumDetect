#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "nn_neuronet.h"
#include "nd_error.h"

int nn_createneuronet(struct nn_neuronet *nnet, int inputc,
	int levelc, const int *neuronc)
{
	int totaln, totalw;
	int i;

	assert(neuronc != NULL && nnet != NULL);

	nnet->levelc = levelc;
	
	if ((nnet->nodec = malloc(sizeof(int) * levelc)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto nodecmallocerror;
	}

	if ((nnet->weightc = malloc(sizeof(int) * levelc)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto weightcmallocerror;
	}

	nnet->inputc = inputc;

	for (i = 0; i < levelc; ++i) {
		nnet->nodec[i] = neuronc[i];
		nnet->weightc[i] = ((i == 0) ? inputc : neuronc[i - 1])
			* neuronc[i];
	}

	totaln = totalw = 0;
	for (i = 0; i < levelc; ++i) {
		totaln += nnet->nodec[i];
		totalw += nnet->weightc[i];
	}

	if ((nnet->weight = malloc(sizeof(double) * totalw)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto weightmallocerror;
	}
	
	if ((nnet->out = malloc(sizeof(double) * totaln)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto outmallocerror;
	}

	if ((nnet->vout = malloc(sizeof(double *) * (levelc + 1))) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto voutmallocerror;
	}

	for (i = 0; i < totalw; ++i)
		nnet->weight[i] = (double) (rand() % 10001) / 5000.0 - 1.0;

	for (i = 0; i < levelc; ++i)
		if ((nnet->vout[i + 1]
			= malloc(sizeof(double) * nnet->nodec[i])) == NULL) {
			nd_seterrormessage(ND_MSGALLOCERROR, __func__);
			goto voutelmallocerror;
		}

	return 0;

voutelmallocerror:
	for (; i >= 0; --i)
		free(nnet->vout[i]);
voutmallocerror:
	free(nnet->out);
outmallocerror:
	free(nnet->weight);
weightmallocerror:
	free(nnet->weightc);
weightcmallocerror:
	free(nnet->nodec);
nodecmallocerror:
	return (-1);
}

int nn_neuroneteval(const struct nn_neuronet *nnet,
	double *input, double *output)
{
	int nl, nn, nw;
	double *pw;

	assert(nnet != NULL && input != NULL);

	nnet->vout[0] = input;

	pw = nnet->weight;
	for (nl = 0; nl < nnet->levelc; ++nl) {
		double *pinput, *pout;
		int wcount;
		
		pout = nnet->vout[nl+1];
		wcount = (nl == 0) ? nnet->inputc : nnet->nodec[nl - 1];
		
		for (nn = 0; nn < nnet->nodec[nl]; ++nn) {
			double val;

			pinput = nnet->vout[nl];
	
			val = 0.0;
			for (nw = 0; nw < wcount; ++nw)
				val += *pw++ * *pinput++;
			
			*pout++ = 1.0 / (1.0 + exp(-val));
		}
	}

	if (output != NULL)
		memcpy(output, nnet->vout[nnet->levelc],
			sizeof(double) * nnet->nodec[nnet->levelc - 1]);

	return 0;
}

int nn_teachneuronet(const struct nn_neuronet *nnet,
	double *input, double *target, double lrate)
{
	int nl, nn, nw;
	double *pw;
	
	assert(nnet != NULL && input != NULL && target != NULL);
	
	if (nn_neuroneteval(nnet, input, NULL) < 0)
		return (-1);

	pw = nnet->weight;
	for (nl = 0; nl < nnet->levelc; ++nl)
		pw += nnet->weightc[nl];

	for (nl = nnet->levelc - 1; nl >= 0; --nl) {
		double *pin, *pout, *nextlevw;
		int inputc;

		pout = nnet->vout[nl+1] + nnet->nodec[nl] - 1;
		inputc = (nl == 0) ? nnet->inputc : nnet->nodec[nl - 1];
		
		for (nn = nnet->nodec[nl] - 1; nn >= 0; --nn) {
			double error;

			if (nl == nnet->levelc - 1)
				error = *pout * (1.0 - *pout) * (target[nn] - *pout);
			else {
				double *perr, *perrw;
				
				perrw = nextlevw + nn;
				perr = nnet->vout[nl+1+1];

				error = 0.0;
				for (nw = 0; nw < nnet->nodec[nl + 1]; ++nw) {
					error += (*perr++ * *perrw);
					perrw += nnet->nodec[nl];
				}
				
				error = *pout * (1.0 - *pout) * error;
			}
			
			pin = nnet->vout[nl] + inputc;
			
			for (nw = 0; nw < inputc; ++nw)
				*(--pw) += lrate * error * *(--pin);	
			
			// after output value (in pout) of neuron was used, 
			// it replaced with error value of this neuron 
			*(pout--) = error;
		}

		nextlevw = pw;
	}

	return 0;
}

int nn_printneuronet(const struct nn_neuronet *nnet)
{
	int nl, nn, nw;
	double *pw = nnet->weight;

	if (nnet == NULL)
		return 0x10;
	
	for (nl = 0; nl < nnet->levelc; ++nl) {
		printf("level %d:\n", nl);
		for (nn = 0; nn < nnet->nodec[nl]; ++nn) {
			int weightc = (nl == 0)
				? nnet->inputc : nnet->nodec[nl - 1];
			
			printf("\tneuron %d weights: ", nn);
			for (nw = 0; nw < weightc; ++nw)
				printf("%f ", *pw++);
			printf("\n");
		}
		
		printf("\n");
	}

	return 0;
}

int nn_neuronettofile(const struct nn_neuronet *nnet, const char *fname)
{
	int nl, nw;
	FILE *nnfile;
	
	assert(nnet != NULL && fname != NULL);

	nnfile = fopen(fname, "w");

	if (fprintf(nnfile, "%u ",nnet->levelc) < 0) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		return (-1);
	}
	
	if (fprintf(nnfile, "%u ",nnet->inputc) < 0) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		return (-1);
	}
	
	for (nl = 0; nl < nnet->levelc; ++nl)
		if (fprintf(nnfile, "%u ", nnet->nodec[nl]) < 0) {
			nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
			return (-1);
		}

	if (fprintf(nnfile, "\n") < 0) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		return (-1);
	}

	int totalw = 0;
	for (nl = 0; nl < nnet->levelc; ++nl)
		totalw += nnet->weightc[nl];

	double *pw = nnet->weight;
	for (nw	= 0; nw < totalw; ++nw)
		if (fprintf(nnfile, "%.10lf\n", *pw++) < 0) {
			nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
			return (-1);
		}

	return 0;
}

int nn_neuronetfromfile(struct nn_neuronet *nnet, const char *fname)
{
	int nl, nw;
	FILE *nnfile;
	int levelc, inputc;
	int *nodec;

	assert(nnet != NULL && fname != NULL);

	nnfile = fopen(fname, "rb");

	if (fscanf(nnfile, "%u", &levelc) == EOF) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		goto scanflevelcerror;
	}
	
	if (fscanf(nnfile, "%u", &inputc) == EOF) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		goto scanfinputcerror;
	}

	if ((nodec = (int *) malloc(sizeof(int) * levelc)) == NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto nodecmallocerror;
	}

	for (nl = 0; nl < levelc; ++nl) {
		if (fscanf(nnfile, "%u", nodec + nl) == EOF) {
			nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
			goto scanfnodecerror;
		}
	}

	if (nn_createneuronet(nnet, inputc, levelc, nodec) < 0)
		goto createneuroneterror;

	int totalw = 0;
	for (nl = 0; nl < levelc; ++nl)
		totalw += nnet->weightc[nl];

	double *pw = nnet->weight;
	for (nw	= 0; nw < totalw; ++nw)
		fscanf(nnfile, "%lf\n", pw++);

	free(nodec);

	return 0;

createneuroneterror:
scanfnodecerror:
	free(nodec);

nodecmallocerror:
scanfinputcerror:
scanflevelcerror:
	return (-1);
}
