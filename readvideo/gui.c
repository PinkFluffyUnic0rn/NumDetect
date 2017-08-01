#include <stdlib.h>
#include <math.h>

#include <xcb/xcb.h>
#include <xcb/xcb_keysyms.h>
#include <cairo/cairo-xcb.h>

#include <sys/time.h>
#include <signal.h>

#include "gui.h"

void *timerud;

static xcb_visualtype_t *
find_visual(xcb_connection_t *c, xcb_visualid_t v)
{
	xcb_screen_iterator_t si;

	si = xcb_setup_roots_iterator(xcb_get_setup(c));

	for (; si.rem; xcb_screen_next(&si)) {
		xcb_depth_iterator_t di;

		di = xcb_screen_allowed_depths_iterator(si.data);

		for (; di.rem; xcb_depth_next(&di)) {
			xcb_visualtype_iterator_t vi;

			vi = xcb_depth_visuals_iterator(di.data);

			for (; vi.rem; xcb_visualtype_next(&vi))
				if (v == vi.data->visual_id)
					return vi.data;
 		}
	}

	return NULL;
}

static int
blocksigalrm(sigset_t *oldmask)
{
	sigset_t newmask;

	sigemptyset(&newmask);
	sigaddset(&newmask, SIGALRM);

	if (sigprocmask(SIG_BLOCK, &newmask, oldmask) < 0)
		return (-1);

	return 0;
}

static int
unblocksigalrm(sigset_t *oldmask)
{
	if (sigprocmask(SIG_SETMASK, oldmask, NULL) < 0)
		return (-1);

	return 0;
}

int
initgui(struct xdata *guidata, int w, int h)
{
	uint32_t mask;
	uint32_t values[2];

	guidata->connection = xcb_connect(NULL, NULL);

	guidata->screen = xcb_setup_roots_iterator(
		xcb_get_setup(guidata->connection)).data;

	guidata->win = xcb_generate_id(guidata->connection);

	mask = XCB_CW_BACK_PIXEL | XCB_CW_EVENT_MASK;
	values[0] = guidata->screen->white_pixel;
	values[1] = XCB_EVENT_MASK_EXPOSURE | XCB_EVENT_MASK_KEY_PRESS
		| XCB_EVENT_MASK_POINTER_MOTION | XCB_EVENT_MASK_BUTTON_PRESS
		| XCB_EVENT_MASK_BUTTON_RELEASE;

	xcb_create_window(guidata->connection, 24, guidata->win,
		guidata->screen->root, 0, 0, w, h, 0,
		XCB_WINDOW_CLASS_INPUT_OUTPUT, guidata->screen->root_visual,
		mask, values);

	xcb_map_window(guidata->connection, guidata->win);

	guidata->xcbsur = cairo_xcb_surface_create(guidata->connection,
		guidata->win,
		find_visual(guidata->connection, guidata->screen->root_visual),
		w, h);
	guidata->sur = cairo_image_surface_create(CAIRO_FORMAT_RGB24, w, h);

	guidata->cr = cairo_create(guidata->xcbsur);

	guidata->defaultcallback = NULL;
	guidata->drawcallback = NULL;
	guidata->keypresscallback = NULL;
	guidata->motioncallback = NULL;
	guidata->buttonpresscallback = NULL;
	guidata->buttonreleasecallback = NULL;

	xcb_flush(guidata->connection);

	return 0;
}

int
windowforceredraw(xcb_connection_t *c, xcb_window_t *win)
{
	xcb_expose_event_t *ee;

	ee = calloc(32, 1);

	ee->window = *win;
	ee->response_type = XCB_EXPOSE;

	ee->x = 0;
	ee->y = 0;
	ee->width = 0;
	ee->height = 0;
	ee->count = 1;

	xcb_send_event(c, 0, *win, XCB_EVENT_MASK_EXPOSURE, (char *)ee);

	xcb_flush(c);

	free(ee);

	return 0;
}

int
mainloop(struct xdata *guidata, void *userdata)
{
	while (1) {
		xcb_generic_event_t *event;
		sigset_t oldmask;
		xcb_key_symbols_t *ks;
		xcb_keysym_t keysym;

		event = xcb_poll_for_event(guidata->connection);

		if (event == NULL) {
			guidata->defaultcallback(userdata);

			windowforceredraw(guidata->connection,
				&(guidata->win));

			continue;
		}

		blocksigalrm(&oldmask);

		switch (event->response_type & ~0x80) {
		case XCB_EXPOSE:
			guidata->drawcallback(guidata->sur, userdata);

			cairo_set_source_surface(guidata->cr,
				guidata->sur, 0, 0);
			cairo_paint(guidata->cr);

			cairo_surface_flush(guidata->sur);

			break;

		case XCB_KEY_PRESS:
			ks = xcb_key_symbols_alloc(guidata->connection);

			keysym = xcb_key_symbols_get_keysym(ks,
				((xcb_key_press_event_t *) event)->detail, 0);

//			guidata->keypresscallback(keysym, userdata);

			xcb_key_symbols_free(ks);

			break;

		case XCB_MOTION_NOTIFY:
/*
			guidata->motioncallback(
				((xcb_motion_notify_event_t *) event)->event_x,
				((xcb_motion_notify_event_t *) event)->event_y,
				userdata);
*/
			break;

		case XCB_BUTTON_PRESS:
/*
			guidata->buttonpresscallback(
				((xcb_button_press_event_t *)event)->detail,
				userdata);
*/
			break;


		case XCB_BUTTON_RELEASE:
/*
				guidata->buttonreleasecallback(
				((xcb_button_press_event_t *)event)->detail,
				userdata);
*/
			break;

		default:
			break;
		}

		unblocksigalrm(&oldmask);

		free(event);
		xcb_flush(guidata->connection);
	}

	return 0;
}

int
settimer(void (*timerfunc)(int), void *userdata, double secs)
{
	struct itimerval tv;

	timerud = userdata;

	if (signal(SIGALRM, timerfunc) == SIG_ERR)
		return 1;

	tv.it_value.tv_sec = tv.it_interval.tv_sec = (int) floor(secs);
	tv.it_value.tv_usec = tv.it_interval.tv_usec
		= (int) (1e6 * (secs - floor(secs)));

	setitimer(ITIMER_REAL, &tv, NULL);

	return 0;
}
