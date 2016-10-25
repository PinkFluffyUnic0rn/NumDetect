#ifndef HC_SYNC_H
#define HC_SYNC_H

int nd_psinitprefork();
int nd_psinitposkfork(int pn);
int nd_pstrylock(int readfd);
int nd_pslock(int readfd);
int nd_psunlock(int writefd);
int nd_psclose();

#endif
