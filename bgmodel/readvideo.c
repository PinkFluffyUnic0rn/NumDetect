#include <stdlib.h>
#include <sys/wait.h>
#include <sys/mman.h>

#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavformat/avio.h>
#include <libavdevice/avdevice.h>

#include <libswscale/swscale.h>

#include "libavutil/channel_layout.h"
#include "libavutil/frame.h"

#include <gtk/gtk.h>
#include <cairo.h>

#include <unistd.h>

#include "nd_image.h"
#include "nd_procsync.h"
#include "nd_error.h"

#include "gui.h"
#include "bgmodel.h"

#define INITFRAMESCOUNT 30
#define INITFRAMESSTEP 25

struct avdata {
	AVFormatContext *s;
	AVCodecContext *vcodecc;
	AVFrame *frame;
	int vstreamid;
};

struct hcdata {
	struct nd_image img;
};

struct playbackstate {
	double framerate;
	struct timespec tstart;
	int64_t timestamp;
};

const int threadcount = 8;

struct nd_image *drawimg;

int openinput(AVFormatContext **s, const char *pathname, int *vstreamid)
{
	int ret;
	int sn;
	
	*s = NULL;

	if ((ret = avformat_open_input(s, pathname,
		NULL, NULL)) < 0) {
//		av_find_input_format("video4linux2"), NULL)) < 0) {
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

	for (sn = 0; sn < (*s)->nb_streams; ++sn)
		if ((*s)->streams[sn]->codec->codec_type == AVMEDIA_TYPE_VIDEO)
			*vstreamid = sn;
	
	return 0;
}

int openvideocodec(AVFormatContext *s, AVCodecContext **vcodecc, int vstreamid)
{
	int ret;
	AVCodec *vcodec;
	
	vcodec = NULL;

	if ((vcodec = avcodec_find_decoder(
		s->streams[vstreamid]->codec->codec_id)) == NULL) {
		fprintf(stderr, "Cannot find decoder.\n");
		return (-1);
	}

	*vcodecc = avcodec_alloc_context3(vcodec);

	avcodec_copy_context(*vcodecc, s->streams[vstreamid]->codec);

	if (vcodec->capabilities & CODEC_CAP_TRUNCATED)
		(*vcodecc)->flags |= CODEC_FLAG_TRUNCATED;

	if ((ret = avcodec_open2(*vcodecc, vcodec, NULL)) < 0) {
		char buf[255];

		av_strerror(ret, buf, 255);
		fprintf(stderr, "%s\n", buf);

		return (-1);
	}

	vcodec = NULL;

	return 0;
}

int initavdata(struct avdata *av, const char *filepath)
{
	av_register_all();
	avdevice_register_all();
	avcodec_register_all();
	avformat_network_init();

	if (openinput(&(av->s), filepath, &(av->vstreamid)) < 0)
		return (-1);

	if (openvideocodec(av->s, &(av->vcodecc),
		av->vstreamid) < 0)
		return (-1);

	av->frame = av_frame_alloc();

	return 0;
}

int imgcreateshared(struct nd_image **img, int w, int h, int format)
{
	size_t imgsz;

	imgsz = sizeof(struct nd_image)
		+ sizeof(double) * w * h * nd_imgchanscount(format);

	if ((*img = mmap(0, imgsz, PROT_READ | PROT_WRITE,
		MAP_ANON | MAP_SHARED, -1, 0)) == MAP_FAILED) {
	
		fprintf(stderr, "Cannot create shared memory mapping\n");
		return (-1);
	}

	(*img)->w = w;
	(*img)->h = h;
	(*img)->format = format;
	(*img)->data = (void *) (*img) + sizeof(struct nd_image);

	return 0;
}

int imgdestroyshared(struct nd_image **img)
{
	if (munmap(*img, sizeof(struct nd_image)
		+ (*img)->w * (*img)->h * sizeof(double)) < 0) {
		fprintf(stderr, "Cannot destroy shared memory mapping\n");
		return (-1);
	}

	*img = NULL;

	return 0;
}

int readframe(AVFormatContext *s, AVCodecContext *vcodecc, int vstreamid,
	AVFrame *frame, int *isdecoded)
{
	AVPacket pack;
	int ret;
	
//	av_init_packet(&pack);	

	memset(&pack, 0, sizeof(AVPacket));
	pack.data = NULL;
	pack.buf = NULL;
	pack.size = 0;

	if ((ret = av_read_frame(s, &pack)) < 0) {
		if (ret != AVERROR_EOF) {
			char buf[255];
			
			av_strerror(ret, buf, 255);
			fprintf(stderr, "%s\n", buf);

			return (-1);
		}

		return 1;
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

int frametorgb(AVFrame *frame, int rgbw, int rgbh, uint8_t **rgbdata,
	int *rgblinesize)
{
	struct SwsContext *swsc;
	
	if ((swsc = sws_getContext(frame->width, frame->height, frame->format,
		rgbw, rgbh, AV_PIX_FMT_RGB24, SWS_BICUBIC, NULL, NULL, NULL))
		== NULL) {
		fprintf(stderr, "Cannot create a scaling context.\n");
		return (-1);
	}
	
	if ((*rgbdata = malloc(sizeof(uint8_t) * rgbw * rgbh * 3)) == NULL) {
		fprintf(stderr, "Cannot allocate memory.\n");
		return (-1);
	}

	*rgblinesize = rgbw * 3;

	if (sws_scale(swsc, (const uint8_t **) frame->data, frame->linesize,
		0, frame->height, rgbdata, rgblinesize) <= 0) {
		fprintf(stderr, "Cannot scale a frame.\n");
		return (-1);
	}

	sws_freeContext(swsc);
	
	return 0;
}

int pbstateinit(struct playbackstate *pbs, double fr)
{
	pbs->framerate = fr;
	pbs->timestamp = 0;

	if (clock_gettime(CLOCK_REALTIME, &(pbs->tstart)) < 0) {
		fprintf(stderr, "Cannot get current time.\n");
		return (-1);
	}

	return 0;
}

int pbnexttimestamp(struct playbackstate *pbs)
{
	struct timespec tcur;
	double tcurd;
	double tframed;

	if (clock_gettime(CLOCK_REALTIME, &tcur) < 0) {
		fprintf(stderr, "Cannot get current time.\n");
		return (-1);
	}
	
	tcurd = 1.0 / pbs->framerate * pbs->timestamp;
	tframed = (tcur.tv_sec - pbs->tstart.tv_sec)
		+ (tcur.tv_nsec - pbs->tstart.tv_nsec) * 1.0e-9;

	if (tcurd > tframed) {
		struct timespec twait;
		double twaitd;

		twaitd = tcurd - tframed;
	
		twait.tv_sec = floor(twaitd);
		twait.tv_nsec = (twaitd - floor(twaitd)) * 1.0e9;
		
		if (nanosleep(&twait, NULL) < 0) {
			fprintf(stderr, "Cannot put process to sleep.\n");
			return (-1);
		}
	}
	
	++(pbs->timestamp);

	return 0;
}

int scanloop(struct nd_image *img)
{
	struct bg_context ctx;
	
	if (bg_initcontext(&ctx, 5, 0.3, 0.05, 640, img->h * 640 / img->w,
		30, 25) < 0)
		return (-1);

	while (1) {
//		struct nd_image imgcopy;
		struct bg_rect *r;
		int rc;

		if (img->data == NULL)
			return 0;

		if (nd_pslock(1) < 0)
			return (-1);

		if (bg_putimgtocontext(&ctx, img, &r, &rc) < 0)
			return (-1);

		if (nd_psunlock(1) < 0)
			return (-1);

		int x, y;
		for (y = 0; y < drawimg->h; ++y)
			for (x = 0; x < drawimg->w; ++x) {
				double *pix;
			
				pix = drawimg->data + (y * drawimg->w + x) * 3;
				memcpy(pix, img->data + (y * drawimg->w + x) * 3,
					sizeof(double) * 3);
			}
	
		if (rc == 0)
			continue;
		
		int i;
		for (i = 0; i < rc; ++i) {
			double *pix;
			int x0, y0, x1, y1;
			int xx, yy;
			
			y0 = drawimg->h * r[i].y0 / ctx.bgm.h;
			y1 = drawimg->h * r[i].y1 / ctx.bgm.h;

			x0 = drawimg->w * r[i].x0 / ctx.bgm.w;
			x1 = drawimg->w * r[i].x1 / ctx.bgm.w;
	
			for (xx = x0; xx < x1; ++xx) {
				pix = drawimg->data + (y0 * drawimg->w + xx) * 3;
				
				pix[0] = 0.0;
				pix[1] = 1.0;
				pix[2] = 0.0;
				
				pix = drawimg->data + ((y1 - 1) * drawimg->w + xx) * 3;
			
				pix[0] = 0.0;
				pix[1] = 1.0;
				pix[2] = 0.0;
			}

			for (yy = y0; yy < y1; ++yy) {
				pix = drawimg->data + (yy * drawimg->w + x0) * 3;
			
				pix[0] = 0.0;
				pix[1] = 1.0;
				pix[2] = 0.0;
			
				pix = drawimg->data + (yy * drawimg->w + (x1 - 1)) * 3;
			
				pix[0] = 0.0;
				pix[1] = 1.0;
				pix[2] = 0.0;

			}
		}

		free(r);

	}

	fprintf(stderr, "Unexpexted loop quit\n");
	return (-1);
}

int rgbtograypix(uint8_t *rgbdata, int rgblinesize, struct nd_image *img,
	int x, int y)
{
	img->data[y * img->w + x] = (rgbdata[y * rgblinesize + x * 3 + 0]
		+ rgbdata[y * rgblinesize + x * 3 + 1]
		+ rgbdata[y * rgblinesize + x * 3 + 2])
		/ 3.0 / 255.0;

	return 0;
}

int rgbtorgbpix(uint8_t *rgbdata, int rgblinesize, struct nd_image *img,
	int x, int y)

{
	img->data[(y * img->w + x)
		* nd_imgchanscount(img->format) + 0]
		= rgbdata[y * rgblinesize + x * 3 + 2]
		/ 255.0;
	img->data[(y * img->w + x)
		* nd_imgchanscount(img->format) + 1]
		= rgbdata[y * rgblinesize + x * 3 + 1]
		/ 255.0;
	img->data[(y * img->w + x)
		* nd_imgchanscount(img->format) + 2]
		= rgbdata[y * rgblinesize + x * 3 + 0]
		/ 255.0;

	return 0;
}

int rgbdatatoimg(uint8_t *rgbdata, int rgblinesize, struct nd_image *img)
{
	int x, y;
	
	for (y = 0; y < img->h; ++y)
		for(x = 0; x < img->w; ++x) {
			if (img->format == ND_PF_GRAYSCALE)
				rgbtograypix(rgbdata, rgblinesize, img, x, y);
			else
				rgbtorgbpix(rgbdata, rgblinesize, img, x, y);
		}			

	return 0;
}

struct decodedata {
	struct avdata *av;
	struct nd_image *img;
	struct playbackstate *pbs;
	struct nd_image *drawimg;
	struct xdata *guidata;
};

int decode(void *userdata)
{
	struct decodedata *ddata;
	int isdecoded;
	uint8_t *rgbdata;
	int rgblinesize;
	int retval;

	ddata = userdata;
	
	isdecoded = 0;

	do {
		if ((retval = readframe(ddata->av->s, ddata->av->vcodecc,
			ddata->av->vstreamid, ddata->av->frame,
			&(isdecoded))) < 0)
			return (-1);

		if (retval > 0)
			return 0;
	} while (!isdecoded);

	if (frametorgb(ddata->av->frame, ddata->img->w, ddata->img->h,
		&rgbdata, &rgblinesize) < 0)
		return (-1);

//	rgbdatatoimg(rgbdata, rgblinesize, ddata->drawimg);
//	windowforceredraw(ddata->guidata->connection, &(ddata->guidata->win));	
	
	if (nd_pstrylock(0)) {
		rgbdatatoimg(rgbdata, rgblinesize, ddata->img);
		
		if (nd_psunlock(0) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}
	}

	free(rgbdata);

	if (pbnexttimestamp(ddata->pbs) < 0)
		return (-1);

	return 0;
}

int draw(cairo_surface_t *sur, void *userdata)
{
	struct decodedata *ddata;
	unsigned char *surdata;
	int surw, surh;
	int i, j;

	ddata = userdata;

	surw = cairo_image_surface_get_width(sur);
	surh = cairo_image_surface_get_height(sur);
	surdata = cairo_image_surface_get_data(sur);
	
	for (i = 0; i < surh; ++i) {
		unsigned char *inputline;
		double *outputline;

		inputline = surdata + i * surw * 4;
		outputline = ddata->drawimg->data + i * surw * 3;

		for (j = 0; j < surw; ++j) {
			inputline[j * 4 + 0] = 255.0 * (outputline[j * 3 + 0]);
			inputline[j * 4 + 1] = 255.0 * (outputline[j * 3 + 1]);
			inputline[j * 4 + 2] = 255.0 * (outputline[j * 3 + 2]);
			inputline[j * 4 + 3] = 255;
		}
	}

	return 0;
}

int decodeloop(struct avdata *av, struct nd_image *img)
{
	struct xdata guidata;
	struct decodedata ddata;
	struct playbackstate pbs;
	AVRational tmpr;

	tmpr = av->s->streams[av->vstreamid]->avg_frame_rate;
	if (pbstateinit(&pbs, (double) tmpr.num / tmpr.den) < 0)
		return 1;

	initgui(&guidata, img->w, img->h);

	ddata.av = av;
	ddata.img = img;
	ddata.pbs = &pbs;
	ddata.drawimg = drawimg;
	ddata.guidata = &guidata;

	guidata.defaultcallback = decode;
	guidata.drawcallback = draw;
	guidata.keypresscallback = NULL;
	guidata.motioncallback = NULL;
	guidata.buttonpresscallback = NULL;
	guidata.buttonreleasecallback = NULL;
	
	mainloop(&guidata, &ddata);

	fprintf(stderr, "Unexpexted loop quit\n");
	return (-1);
}

int main(int argc, char **argv)
{
	int chpid;
	struct avdata av;
	struct nd_image *img;

	if (initavdata(&av, argv[1]) < 0)
		return 1;
	
	if (imgcreateshared(&img, av.vcodecc->width, av.vcodecc->height,
		ND_PF_RGB) < 0)
		return 1;

	if (imgcreateshared(&drawimg, av.vcodecc->width, av.vcodecc->height,
		ND_PF_RGB) < 0)
		return 1;

	if (nd_psinitprefork() < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if ((chpid = fork()) == 0) {
		if (nd_psinitpostfork(1) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}
		
		if (scanloop(img) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}
		
		if (nd_psclose(1) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}
	}
	else {
		int retval;
	
		if (nd_psinitpostfork(0) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}

		if (decodeloop(&av, img) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}

		if (nd_pslock(0) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}
		
		img->data = NULL;
		
		if (nd_psunlock(0) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}	
		
		wait(&retval);
		
		if (nd_psclose(0) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}
	}

	if (imgdestroyshared(&img) < 0)
		return 1;

	av_frame_unref(av.frame);
	avcodec_free_context(&(av.vcodecc));
	avformat_close_input(&(av.s));
	
	return 0;
}
 
