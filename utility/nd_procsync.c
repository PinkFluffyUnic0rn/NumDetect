#include <unistd.h>
#include <sys/time.h>
#include <assert.h>

#include "nd_error.h"
#include "nd_procsync.h"

static int pfdc[2];
static int pfdp[2];

int nd_psinitprefork()
{	
	if (pipe(pfdc) < 0) {
		nd_seterrormessage("pipe error", __func__);
		return (-1);
	}

	if (pipe(pfdp) < 0) {
		nd_seterrormessage("pipe error", __func__);
		return (-1);
	}
	
	if (nd_psunlock(0) < 0) {
		nd_seterrormessage("pipe error", __func__);
		return (-1);
	}

	return 0;
}

int nd_psinitpostfork(int pn)
{
	assert(pn <= 1 && pn >= 0);

	if (pn == 0) {
		if (close(pfdc[0]) < 0 || close(pfdp[1]) < 0) {
			nd_seterrormessage("pipe error", __func__);
			return (-1);
		}
	}
	else {
		if (close(pfdc[1]) < 0 || close(pfdp[0]) < 0) {
			nd_seterrormessage("pipe error", __func__);
			return (-1);
		}
	}

	return 0;
}

int nd_pstrylock(int pn)
{
	fd_set readfs;
	struct timeval tv;
	int retval;
	int fd;

	assert(pn <= 1 && pn >= 0);

	fd = (pn == 0) ? pfdp[0] : pfdc[0];

	FD_ZERO(&readfs);
	FD_SET(fd, &readfs);

	tv.tv_sec = 0;
	tv.tv_usec = 0;

	if ((retval = select(fd + 1, &readfs, NULL, NULL, &tv)) < 0) {
		nd_seterrormessage("pipe error", __func__);
		return (-1);
	}

	if (retval > 0) {
		if (nd_pslock(pn) < 0)
			return (-1);

		return 1;
	}
	
	return retval;
}

int nd_pslock(int pn)
{
	int c;
	char buf[1024];
	int fd;

	assert(pn <= 1 && pn >= 0);

	fd = (pn == 0) ? pfdp[0] : pfdc[0];

	if ((c = read(fd, buf, 1024)) < 0) {
		nd_seterrormessage("pipe error", __func__);
		return (-1);
	}

	return 0;
}

int nd_psunlock(int pn)
{
	int c;
	int fd;

	assert(pn <= 1 && pn >= 0);

	fd = (pn == 0) ? pfdc[1] : pfdp[1];
	
	if ((c = write(fd, "u", 1)) < 0) {
		nd_seterrormessage("pipe error", __func__);
		return (-1);
	}

	return 0;
}

int nd_psclose(int pn)
{
	assert(pn <= 1 && pn >= 0);
	
	if (pn == 0) {
		if (close(pfdp[0]) < 0 || close(pfdc[1]) < 0) {
			nd_seterrormessage("pipe error", __func__);
			return (-1);
		}
	}
	else {
		if (close(pfdp[1]) < 0 || close(pfdc[0]) < 0) {
			nd_seterrormessage("pipe error", __func__);
			return (-1);
		}
	}

	return 0;
}
