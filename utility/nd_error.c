#include <stdio.h>

#include "nd_error.h"

void nd_seterrormessage(const char *msg, const char *funcname)
{
	sprintf(errormessage, "%s: %s.\n", funcname, msg);
}

const char *nd_geterrormessage()
{
	return errormessage;
}
