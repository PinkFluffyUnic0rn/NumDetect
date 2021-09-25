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

#endif
