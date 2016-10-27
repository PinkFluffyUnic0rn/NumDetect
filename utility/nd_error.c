#include "nd_error.h"

enum ND_ERROR nd_error;

const char *terr[] = {
	"Wrong function argument",
	"Cannot allocate memory",
	"Cannot open file",
	"Cannot close file",
	"Cannot open .png file",
	"Wrong file format",
	"Error in cairo function",
	"Image argument has wrong parameters",
	"Matrix argument has wrong parameters",
	"Cannot read from file",
	"Cannot write to file",
	"Pipe call failed",
	"Read call failed",
	"Write call failed",
	"Select call failed",
	"Close call failed",
};

void nd_seterror(enum ND_ERROR err)
{
	nd_error = err;	
}

const char *nd_strerror(enum ND_ERROR err)
{
	if (err < 0 || err > ND_ERRORCOUNT)
		return "Unknown error";

	return terr[err];
}
