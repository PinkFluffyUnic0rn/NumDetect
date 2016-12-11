#ifndef ND_ERROR_H
#define ND_ERROR_H

/*
static char const *ND_MSGALLOCERROR = "cannot alocate memory";
static char const *ND_MSGFILEIOERROR = "cannot read or write file";
*/

#define ND_MSGALLOCERROR "cannot alocate memory"
#define ND_MSGFILEIOERROR "cannot read or write file"

char errormessage[1024];

void nd_seterrormessage(const char *msg, const char *funcname);
const char *nd_geterrormessage();

/*
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
	ND_WRITEFILEERROR,
	ND_PIPEERROR,
	ND_READERROR,
	ND_WRITEERROR,
	ND_SELECTERROR,
	ND_CLOSEERROR,
	ND_ERRORCOUNT
};

extern enum ND_ERROR nd_error;

void nd_seterror(enum ND_ERROR err);
const char *nd_strerror(enum ND_ERROR err);
*/

#endif
