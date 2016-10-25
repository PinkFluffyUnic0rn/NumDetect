#include "hc_procsync.h"

#include <unistd.h>
#include <sys/time.h>

static int pfdc[2];
static int pfdp[2];

int nd_psinitprefork()
{	
	pipe(pfdc);	
	pipe(pfdp);	
	
	nd_psunlock(0);

	return 0;
}

int nd_psinitprostfork(int pn)
{
	if (pn == 0) {
		close(pfdc[0]);
		close(pfdp[1]);
	}
	else {
		close(pfdc[1]);
		close(pfdp[0]);
	}

	return 0;
}

int nd_pstrylock(int pn)
{
	fd_set readfs;
	struct timeval tv;
	int retval;
	int fd;

	fd = (pn == 0) ? pfdp[0] : pfdc[0];

	FD_ZERO(&readfs);
	FD_SET(fd, &readfs);

	tv.tv_sec = 0;
	tv.tv_usec = 0;

	retval = select(fd + 1, &readfs, NULL, NULL, &tv);

	if (retval > 0) {
		char buf[1024];
		if (read(fd, buf, 1024) < 0)
			return -1;

		return 1;
	}

	return retval;
}

int nd_pslock(int pn)
{
	int c;
	char buf[1024];
	int fd;

	fd = (pn == 0) ? pfdp[0] : pfdc[0];

	if ((c = read(fd, buf, 1024)) < 0)
		return -1;

	return 0;
}

int nd_psunlock(int pn)
{
	int c;
	int fd;

	fd = (pn == 0) ? pfdc[1] : pfdp[1];
	
	if ((c = write(fd, "u", 1)) < 0)
		return -1;

	return 0;
}

int nd_psclose()
{
	close(pfdp[0]);
	close(pfdp[1]);

	close(pfdc[0]);
	close(pfdc[1]);

	return 0;
}
