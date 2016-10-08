#include "nd_error.h"

enum ND_ERROR nd_error;

void nd_seterror(enum ND_ERROR err)
{
	nd_error = err;	
}

const char *nd_strerror(enum ND_ERROR err)
{
	switch(err) {
	case ND_INVALIDARG:
		return "Wrong function argument";

	case ND_ALLOCFAULT:
		return "Cannot allocate memory";

	case ND_FOPENERROR:
		return "Cannot open file";

	case ND_FCLOSEERROR:
		return"Cannot close file";

	case ND_OPENPNGERROR:
		return "Cannot open .png file";

	case ND_WRONGFORMAT:
		return "Wrong file format";

	case ND_CAIROERROR:
		return "Error in cairo function";

	case ND_INVALIDIMAGE:
		return "Image argument has wrong parameters";

	case ND_INVALIDMATRIX:
		return "Matrix argument has wrong parameters";

	case ND_READFILEERROR:
		return "Cannot read file";

	default:
		return "Unknown error";
	}
}
