#include <stdlib.h>
#include <libavcodec/avcodec.h>

#include <libavformat/avformat.h>
#include <libavformat/avio.h>

#include <libswscale/swscale.h>

#include "libavutil/channel_layout.h"
#include "libavutil/frame.h"

#include <gtk/gtk.h>
#include <cairo.h>

#include "nd_image.h"
#include "hc_hcascade.h"
#include "hc_scanimgpyr.h"

#define IMGWIDTH 640

struct gtkdata {
	GtkWidget *mainwindow;
	GtkWidget *drawarea;
};

struct avdata {
	AVFormatContext *s;
	AVCodecContext *vcodecc;
	AVFrame *frame;
	int vstreamid;
};

struct hcdata {
	struct hc_hcascade hc;
	struct nd_image img;
	struct hc_scanconfig scanconf;
	pthread_mutex_t framemutex;
	int continuescan;
};

struct alldata {
	struct gtkdata gui;
	struct avdata av;
	struct hcdata hc;
	char *outputdir;
};

int imgtocairosur(const struct nd_image *img, cairo_surface_t **sur)
{
	unsigned char *surdata;
	int imgy, imgx;

	*sur = cairo_image_surface_create(CAIRO_FORMAT_RGB24, img->w, img->h);	
	
	if (cairo_surface_status(*sur) != CAIRO_STATUS_SUCCESS)
		return (-1);
	
	if ((surdata = cairo_image_surface_get_data(*sur)) == NULL)
		return (-1);

	for (imgy = 0; imgy < img->h; ++imgy)
		for (imgx = 0; imgx < img->w; ++imgx) {
			int val;
		
			val = ceil(img->data[(imgy * img->w + imgx)]
				* 255.0);
			val = (val <= 255)
				? (val >= 0 ? val : 0) : 255;
				
			surdata[(imgy * img->w + imgx) * 4 + 0] = val;
			surdata[(imgy * img->w + imgx) * 4 + 1] = val;
			surdata[(imgy * img->w + imgx) * 4 + 2] = val;
		}

	return 0;
}

int openinput(AVFormatContext **s, const char *pathname, int *vstreamid)
{
	int ret;
	
	*s = NULL;

	if ((ret = avformat_open_input(s, pathname, NULL, NULL)) < 0) {
		char buf[255];
		av_strerror(ret, buf, 255);
		fprintf(stderr, "%s\n", buf);

		return (-1);
	}
	
	if ((ret = avformat_find_stream_info(*s, NULL)) < 0) {
		char buf[255];
		av_strerror(ret, buf, 255);
		fprintf(stderr, "%s\n", buf);

		return (-1);
	}

	int sn;

	for (sn = 0; sn < (*s)->nb_streams; ++sn) {
		if ((*s)->streams[sn]->codec->codec_type == AVMEDIA_TYPE_VIDEO)
			*vstreamid = sn;
	}
	
	return 0;
}

int openvideocodec(AVFormatContext *s, AVCodecContext **vcodecc, int vstreamid)
{
	int ret;
	AVCodec *vcodec;
	
	vcodec = NULL;

	if ((vcodec = avcodec_find_decoder(
		s->streams[vstreamid]->codec->codec_id)) == NULL) {
		
		return (-1);
	}

	*vcodecc = avcodec_alloc_context3(vcodec);

	if (vcodec->capabilities & CODEC_CAP_TRUNCATED)
		(*vcodecc)->flags |= CODEC_FLAG_TRUNCATED;

	avcodec_copy_context(*vcodecc, s->streams[vstreamid]->codec);

	if ((ret = avcodec_open2(*vcodecc, vcodec, NULL)) < 0) {
		char buf[255];
		av_strerror(ret, buf, 255);
		fprintf(stderr, "%s\n", buf);

		return (-1);
	}

	vcodec = NULL;

	return 0;
}

int readframe(AVFormatContext *s, AVCodecContext *vcodecc, int vstreamid,
	AVFrame *frame, int *isdecoded)
{
	AVPacket pack;
	int ret;
	
	av_init_packet(&pack);	

	memset(&pack, 0, sizeof(AVPacket));
	pack.data = NULL;
	pack.buf = NULL;
	pack.size = 0;

	if ((ret = av_read_frame(s, &pack)) < 0) {
		if (ret != AVERROR_EOF) {
			char buf[255];
			
			av_strerror(ret, buf, 255);
			fprintf(stderr, "%s\n", buf);
		}

		return (-1);
	}
		
	if (pack.stream_index == vstreamid) {
		int cret;

		if ((cret = avcodec_decode_video2(vcodecc, frame,
			isdecoded, &pack)) < 0) {
			char buf[255];
			
			av_strerror(cret, buf, 255);
			fprintf(stderr, "%s\n", buf);

			return (-1);
		}
	}

	av_packet_unref(&pack);

	return 0;
}

int draw(GtkWidget *wgt, cairo_t *cr, gpointer data)
{
	struct hcdata *hcd;
	cairo_surface_t *sur;
	
	hcd = data;
	
	if (pthread_mutex_lock(&(hcd->framemutex))) {
		perror("pthread_mutex_lock: ");
		return 1;
	}
			
	if (imgtocairosur(&(hcd->img), &sur))
		return 1;

	if (pthread_mutex_unlock(&(hcd->framemutex))) {
		cairo_surface_destroy(sur);

		perror("pthread_mutex_unlock: ");
		return 1;
	}
		
	cairo_set_source_surface(cr, sur, 0, 0);
	cairo_paint(cr);

	cairo_surface_destroy(sur);

	return 0;
}

int frametorgb(AVFrame *frame, int rgbw, int rgbh, uint8_t **rgbdata,
	int *rgblinesize)
{
	struct SwsContext *swsc;
	
	if ((swsc = sws_getContext(frame->width, frame->height, frame->format,
		rgbw, rgbh, AV_PIX_FMT_RGB24, SWS_BICUBIC, NULL, NULL, NULL))
		== NULL)
		return (-1);
	
	if ((*rgbdata = malloc(sizeof(uint8_t) * rgbw * rgbh * 3)) == NULL) {
		printf("Cannot allocate memory.\n");
		return (-1);
	}

	*rgblinesize = rgbw * 3;

	if (sws_scale(swsc, (const uint8_t **) frame->data, frame->linesize,
		0, frame->height, rgbdata, rgblinesize) <= 0)
		return (-1);

	sws_freeContext(swsc);
	
	return 0;
}


int detectedtofile(struct nd_image *img, struct hc_rect *r, int rc,
	const char *outputdir)
{
	char imgpath[255];
	static int foundn = 0;
	int rn;

	if (sprintf(imgpath, "%s/%d.png", outputdir, foundn) < 0) {
		fprintf(stderr, "sprintf: %s.\n", strerror(errno));
		return (-1);
	}
	
	if (nd_imgwrite(img, imgpath) < 0) {
		fprintf(stderr, "nd_imgwrite: %s.\n",
			nd_strerror(nd_error));
		return (-1);
	}
	
	for (rn = 0; rn < rc; ++rn) {
		if (sprintf(imgpath, "%s/%d_%d.png", outputdir, foundn, rn)
			< 0) {
			fprintf(stderr, "sprintf: %s.\n", strerror(errno));
			return (-1);
		}
		
		double rh = abs(r[rn].y0 - r[rn].y1);

		r[rn].y0 = (int) (r[rn].y0) - (rh * 0.15);
		r[rn].y1 = r[rn].y1 + (rh * 0.15);

		if (r[rn].y0 >= 0 && r[rn].y1 < img->h) {
			struct nd_image imginwin;

			if (nd_imgcreate(&imginwin, abs(r[rn].x1 - r[rn].x0),
				abs(r[rn].y1 - r[rn].y0), img->format) < 0) {
				fprintf(stderr, "nd_imgcreate: %s.\n",
					nd_strerror(nd_error));
				return (-1);
			}
			
			if (nd_imgcrop(img, r[rn].x0, r[rn].y0,
				abs(r[rn].x1 - r[rn].x0),
				abs(r[rn].y1 - r[rn].y0), &imginwin) < 0) {
				nd_imgdestroy(&imginwin);
				
				fprintf(stderr, "nd_imgcrop: %s.\n",
					nd_strerror(nd_error));
				return (-1);
			}
			
			if (nd_imgwrite(&imginwin, imgpath) < 0) {
				nd_imgdestroy(&imginwin);
				
				fprintf(stderr, "nd_imgwrite: %s.\n",
					nd_strerror(nd_error));
				return (-1);
			}

			nd_imgdestroy(&imginwin);
		}
	}

	++foundn;

	return 0;
}

void *scanframe(void *ud)
{
	struct alldata *data;
	struct hc_rect *r;
	int rc;

	data = ud;

	while (1) {
		struct nd_image img;
	
		if (pthread_mutex_lock(&(data->hc.framemutex))) {
			perror("pthread_mutex_lock: ");
			return NULL;
		}
	
		if (!(data->hc.continuescan))
			return NULL;

		if (nd_imgcreate(&img, data->hc.img.w, data->hc.img.h,
			data->hc.img.format) < 0) {
			fprintf(stderr, "nd_imgcreate: %s.\n",
				nd_strerror(nd_error));
		}

		memcpy(img.data, data->hc.img.data,
			sizeof(double) * data->hc.img.w * data->hc.img.h);
	
		if (pthread_mutex_unlock(&(data->hc.framemutex))) {
			nd_imgdestroy(&img);

			perror("pthread_mutex_unlock: ");
			return NULL;
		}
/*
		struct timespec ts, te;
		clock_gettime(CLOCK_REALTIME, &ts);
*/
		if (nd_imgpyramidscan(&(data->hc.hc), &img, &r, &rc,
			&(data->hc.scanconf)) < 0) {
			nd_imgdestroy(&img);
			
			fprintf(stderr, "nd_imgpyramidscan: %s.\n",
				nd_strerror(nd_error));
		}
/*		
		clock_gettime(CLOCK_REALTIME, &te);
		printf("%lf\n", (te.tv_sec - ts.tv_sec)
			+ (te.tv_nsec - ts.tv_nsec) * 1e-9);
*/
/*
		if (rc > 0)
			printf("found!\n");
*/
		if (rc) {
			detectedtofile(&img, r, rc, data->outputdir);
			free(r);
		}
		
		nd_imgdestroy(&img);
	}

	return NULL;
}

int persptransformframe(struct nd_image *frame)
{
	struct nd_matrix3 m;
	double inpoints[8];
	double outpoints[8];

	inpoints[0] = 305; inpoints[1] = 242;
	inpoints[2] = 367; inpoints[3] = 257;
	inpoints[4] = 367; inpoints[5] = 272;
	inpoints[6] = 305; inpoints[7] = 255;

	outpoints[0] = 100 + frame->w / 2 - 4 * 31 / 3;
	outpoints[1] = 100 + frame->h / 2 - 4 * 7 / 3;
	
	outpoints[2] = 100 + frame->w / 2 + 4 * 31 / 3;
	outpoints[3] = 100 + frame->h / 2 - 4 * 7 / 3;
	
	outpoints[4] = 100 + frame->w / 2 + 4 * 31 / 3;
	outpoints[5] = 100 + frame->h / 2 + 4 * 7 / 3;
	
	outpoints[6] = 100 + frame->w / 2 - 4 * 31 / 3;
	outpoints[7] = 100 + frame->h / 2 + 4 * 7 / 3;

	if (nd_getpersptransform(inpoints, outpoints, &m) < 0) {
		fprintf(stderr, "nd_getpersptransform: %s.\n",
			nd_strerror(nd_error));
		return (-1);
	}

	if (nd_imgapplytransform(frame, &m) < 0) {
		fprintf(stderr, "nd_imgapplytransform: %s.\n",
			nd_strerror(nd_error));
		return (-1);
	}

	return 0;
}

int waitfortimestamp(AVRational *tbase, struct timespec *tstart,
	struct timespec *tcur, int64_t *timestamp)
{
	double tcurd;
	double tframed;

	tcurd = (double) tbase->den / tbase->num * *timestamp;
	tframed = (tcur->tv_sec - tstart->tv_sec)
		+ (tcur->tv_nsec - tstart->tv_nsec) * 1.0e-9;

	if (tcurd > tframed) {
		struct timespec twait;
		
		twait.tv_sec = floor(tcurd) - floor(tframed);
		twait.tv_nsec = ((tcurd - floor(tcurd))
			- (tframed - floor(tframed))) * 1.0e9;

		nanosleep(&twait, NULL);
	}

	return 0;
}

void *decodeframe(void *ud)
{
	struct alldata *data;
	int isdecoded;
	uint8_t *rgbdata;
	int rgblinesize;
	int x, y;
/*
	struct timespec tstart;
	struct timespec tcur;
	int64_t timestamp;
*/
	data = ud;
/*
	clock_gettime(CLOCK_REALTIME, &tstart);
	timestamp = 0;
*/
	while (1) {
		do {
			if (readframe(data->av.s, data->av.vcodecc,
				data->av.vstreamid, data->av.frame,
				&(isdecoded)) < 0)
				return NULL;
		} while (!isdecoded);
		
		if (frametorgb(data->av.frame, data->hc.img.w, data->hc.img.h,
			&rgbdata, &rgblinesize) < 0)
			continue;

		if (pthread_mutex_lock(&(data->hc.framemutex))) {
			perror("pthread_mutex_lock: ");
			return NULL;
		}
		
		for (y = 0; y < data->hc.img.h; ++y)
			for(x = 0; x < data->hc.img.w; ++x) {
				data->hc.img.data[y * data->hc.img.w + x]
					= (rgbdata[y * rgblinesize + x * 3 + 0]
					+ rgbdata[y * rgblinesize + x * 3 + 1]
					+ rgbdata[y * rgblinesize + x * 3 + 2])
					/ 3.0 / 255.0;
			}

		if (persptransformframe(&(data->hc.img)) < 0)
			return NULL;

		if (pthread_mutex_unlock(&(data->hc.framemutex))) {
			perror("pthread_mutex_unlock: ");
			return NULL;
		}
		
		free(rgbdata);
/*
		clock_gettime(CLOCK_REALTIME, &tcur);
		waitfortimestamp(&(data->av.s->streams[data->av.vstreamid]
			->avg_frame_rate), &tstart, &tcur, &timestamp);
*/
//		gtk_widget_queue_draw(data->gui.mainwindow);
	
//		++timestamp;
	}

	return NULL;
}

int initgui(struct gtkdata *guidata, struct hcdata *hcd)
{
	guidata->mainwindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	guidata->drawarea = gtk_drawing_area_new();
	
	gtk_container_add(GTK_CONTAINER(guidata->mainwindow),
		guidata->drawarea);

	g_signal_connect(guidata->mainwindow, "destroy",
		G_CALLBACK(gtk_main_quit), NULL);
	g_signal_connect(G_OBJECT(guidata->drawarea), "draw",
		G_CALLBACK(draw), hcd);
	
	gtk_widget_show_all(guidata->mainwindow);

	gdk_window_set_events(gtk_widget_get_window(guidata->mainwindow), 
		GDK_ALL_EVENTS_MASK);

	return 0;
}

int main(int argc, char **argv)
{
	struct alldata data;
	pthread_t decodethread;
	pthread_t scanthread;
	 
	av_register_all();
	avformat_network_init();

	if (openinput(&(data.av.s), argv[1], &(data.av.vstreamid)) < 0)
		return 1;
	
	if (openvideocodec(data.av.s, &(data.av.vcodecc),
		data.av.vstreamid) < 0)
		return 1;

	data.av.frame = av_frame_alloc();

	if (nd_imgcreate(&(data.hc.img), IMGWIDTH, data.av.vcodecc->height
		* IMGWIDTH / data.av.vcodecc->width, ND_PF_GRAYSCALE) < 0) {
		fprintf(stderr, "nd_imgcreate: %s.\n", nd_strerror(nd_error));
		return 1;
	}

	if (hc_hcascaderead(&(data.hc.hc), argv[2]) < 0) {
		fprintf(stderr, "nd_hcascaderead: %s.\n",
			nd_strerror(nd_error));
		return 1;
	}

	data.hc.scanconf.scalestep = 0.9;
	data.hc.scanconf.winhstep = 1;
	data.hc.scanconf.winwstep = 1;

	data.hc.continuescan = 1;
/*
	gtk_init(&argc, &argv);

	initgui(&(data.gui), &(data.hc));	
*/
	data.outputdir = argv[3];
	
	if (pthread_mutex_init(&(data.hc.framemutex), NULL))
		return 1;

	if (pthread_create(&decodethread, NULL, decodeframe, &data))
		return 1;

	if (pthread_create(&scanthread, NULL, scanframe, &data))
		return 1;
/*
	gtk_window_resize(GTK_WINDOW(data.gui.mainwindow),
		data.hc.img.w, data.hc.img.h);

	gtk_main();
*/
	if (pthread_join(decodethread, NULL))
		return 1;
	
	data.hc.continuescan = 0;

	if (pthread_join(scanthread, NULL))
		return 1;
	
	av_frame_unref(data.av.frame);
	avcodec_free_context(&(data.av.vcodecc));
	avformat_close_input(&(data.av.s));
	
	return 0;
}
