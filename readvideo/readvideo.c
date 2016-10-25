#include <stdlib.h>
#include <sys/wait.h>
#include <sys/mman.h>

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
#include "hc_procsync.h"

#define IMGWIDTH 640

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
};

struct playbackstate {
	double framerate;
	struct timespec tstart;
	int64_t timestamp;
};

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

int initavdata(struct avdata *av, const char *filepath)
{
	av_register_all();
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
		+ w * h * sizeof(double);

	if ((*img = mmap(0, imgsz, PROT_READ | PROT_WRITE,
		MAP_ANON | MAP_SHARED, -1, 0)) == MAP_FAILED)
		return (-1);

	(*img)->w = w;
	(*img)->h = h;
	(*img)->format = ND_PF_GRAYSCALE;
	(*img)->data = (void *) (*img) + sizeof(struct nd_image);

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

int pbstateinit(struct playbackstate *pbs, double fr)
{
	pbs->framerate = fr;
	pbs->timestamp = 0;


	clock_gettime(CLOCK_REALTIME, &(pbs->tstart));

	return 0;
}

int pbnexttimestamp(struct playbackstate *pbs)
{
	struct timespec tcur;
	double tcurd;
	double tframed;

	clock_gettime(CLOCK_REALTIME, &tcur);
	
	tcurd = 1.0 / pbs->framerate * pbs->timestamp;
	tframed = (tcur.tv_sec - pbs->tstart.tv_sec)
		+ (tcur.tv_nsec - pbs->tstart.tv_nsec) * 1.0e-9;

	if (tcurd > tframed) {
		struct timespec twait;
		
		twait.tv_sec = floor(tcurd) - floor(tframed);
		twait.tv_nsec = ((tcurd - floor(tcurd))
			- (tframed - floor(tframed))) * 1.0e9;

		nanosleep(&twait, NULL);
	}

	++(pbs->timestamp);

	return 0;
}

int main(int argc, char **argv)
{
	int chpid;
	struct avdata av;
	struct nd_image *img;

	if (initavdata(&av, argv[1]) < 0)
		return 1;

	if (imgcreateshared(&img, IMGWIDTH, av.vcodecc->height
		* IMGWIDTH / av.vcodecc->width, ND_PF_GRAYSCALE) < 0)
		return 1;

	nd_psinitprefork();

	if ((chpid = fork()) == 0) {
		struct hcdata hcd;
	
		if (hc_hcascaderead(&(hcd.hc), argv[2]) < 0) {
			fprintf(stderr, "nd_hcascaderead: %s.\n",
				nd_strerror(nd_error));
			return 1;
		}

		hcd.scanconf.scalestep = 0.9;
		hcd.scanconf.winhstep = 1;
		hcd.scanconf.winwstep = 1;

		while (1) {
			struct hc_rect *r;
			int rc;

			nd_pslock(1);

			if (nd_imgpyramidscan(&(hcd.hc), img, &r, &rc,
				&(hcd.scanconf)) < 0) {
				fprintf(stderr, "nd_imgpyramidscan: %s.\n",
					nd_strerror(nd_error));
	
				return (-1);
			}

			if (rc) {
				detectedtofile(img, r, rc, argv[3]);
				free(r);
			}
	
			nd_psunlock(1);
		}
		
		nd_psclose();
	}
	else {
		struct playbackstate pbs;
		AVRational tmpr;
		int retval;

		tmpr = av.s->streams[av.vstreamid]->avg_frame_rate;	
		pbstateinit(&pbs, (double) tmpr.num / tmpr.den);

		while (1) {
			int isdecoded;
			uint8_t *rgbdata;
			int rgblinesize;
			int x, y;

			do {
				if (readframe(av.s, av.vcodecc, av.vstreamid,
					av.frame, &(isdecoded)) < 0)
					return 1;
			} while (!isdecoded);

			if (frametorgb(av.frame, img->w, img->h,
				&rgbdata, &rgblinesize) < 0)
				return 1;

			if (nd_pstrylock(0)) {
				for (y = 0; y < img->h; ++y)
					for(x = 0; x < img->w; ++x) {
						img->data[y * img->w + x]
							= (rgbdata[y * rgblinesize + x * 3 + 0]
							+ rgbdata[y * rgblinesize + x * 3 + 1]
							+ rgbdata[y * rgblinesize + x * 3 + 2])
							/ 3.0 / 255.0;
					}
				
				if (persptransformframe(img) < 0)
					return 1;
				
				nd_psunlock(0);
			}

			free(rgbdata);

			pbnexttimestamp(&pbs);

		}

		wait(&retval);

		nd_psclose();
	}

	av_frame_unref(av.frame);
	avcodec_free_context(&(av.vcodecc));
	avformat_close_input(&(av.s));
	
	return 0;
}
