#ifndef ND_ERROR_H
#define ND_ERROR_H

enum ND_ERROR
{
	ND_INVALIDARG,
	ND_ALLOCFAULT,
	ND_FOPENERROR,
	ND_FCLOSEERROR,
	ND_OPENPNGERROR,
	ND_WRONGFORMAT,
	ND_CAIROERROR,
	ND_INVALIDIMAGE,
	ND_INVALIDMATRIX,
	ND_READFILEERROR,
	ND_WRITETOFILEFERROR
};

extern enum ND_ERROR nd_error;

void nd_seterror(enum ND_ERROR err);
const char *nd_strerror(enum ND_ERROR err);

#endif
