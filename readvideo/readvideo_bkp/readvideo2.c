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
#include "nd_procsync.h"
#include "hc_hcascade.h"
#include "hc_scanimgpyr.h"
#include "ed_findborder.h"

#include "gui.h"

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


const int threadcount = 8;

int openinput(AVFormatContext **s, const char *pathname, int *vstreamid)
{
	int ret;
	int sn;
	
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

int getperspmat(int w, int h, struct nd_matrix3 *m)
{
	double inpoints[8];
	double outpoints[8];

	inpoints[0] = 367; inpoints[1] = 242;
	inpoints[2] = 429; inpoints[3] = 257;
	inpoints[4] = 429; inpoints[5] = 272;
	inpoints[6] = 367; inpoints[7] = 255;

	outpoints[0] = 100 + w / 2 - 4 * 31 / 3;
	outpoints[1] = 100 + h / 2 - 4 * 7 / 3;
	
	outpoints[2] = 100 + w / 2 + 4 * 31 / 3;
	outpoints[3] = 100 + h / 2 - 4 * 7 / 3;
	
	outpoints[4] = 100 + w / 2 + 4 * 31 / 3;
	outpoints[5] = 100 + h / 2 + 4 * 7 / 3;
	
	outpoints[6] = 100 + w / 2 - 4 * 31 / 3;
	outpoints[7] = 100 + h / 2 + 4 * 7 / 3;

	if (nd_getpersptransform(inpoints, outpoints, m) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	return 0;
}

int recttoorig(struct hc_rect *r, struct nd_matrix3 *perspmat,
	double addtop, double addbot, double imgtoorigrel,
	struct nd_vector3 p[4])
{
	p[0].x = (double) r->x0;
	p[0].y = r->y0 - addtop;
	p[0].z = 1.0;
	
	p[1].x = (double) r->x1;
	p[1].y = r->y0 - addtop;
	p[1].z = 1.0;

	p[2].x = (double) r->x1;
	p[2].y = r->y1 + addbot;
	p[2].z = 1.0;

	p[3].x = (double) r->x0;
	p[3].y = r->y1 + addbot;
	p[3].z = 1.0;

	nd_v3m3mult(p + 0, perspmat, p + 0);
	nd_v3m3mult(p + 1, perspmat, p + 1);
	nd_v3m3mult(p + 2, perspmat, p + 2);
	nd_v3m3mult(p + 3, perspmat, p + 3);

	p[0].x /= p[0].z; p[0].y /= p[0].z;
	p[1].x /= p[1].z; p[1].y /= p[1].z;
	p[2].x /= p[2].z; p[2].y /= p[2].z;
	p[3].x /= p[3].z; p[3].y /= p[3].z;

	p[0].x *= imgtoorigrel; p[0].y *= imgtoorigrel;
	p[1].x *= imgtoorigrel; p[1].y *= imgtoorigrel;
	p[2].x *= imgtoorigrel; p[2].y *= imgtoorigrel;
	p[3].x *= imgtoorigrel; p[3].y *= imgtoorigrel;

	return 0;
}

int getorignum(struct nd_image *imgorig, struct nd_vector3 p[4],
	int winw, int winh, struct nd_image *imginwin)
{
	struct nd_matrix3 persporig;
	double inpoints[8];
	double outpoints[8];

	inpoints[0] = p[0].x; inpoints[1] = p[0].y;
	inpoints[2] = p[1].x; inpoints[3] = p[1].y;
	inpoints[4] = p[2].x; inpoints[5] = p[2].y;
	inpoints[6] = p[3].x; inpoints[7] = p[3].y;

	outpoints[0] = 0; 	outpoints[1] = 0;
	outpoints[2] = winw; 	outpoints[3] = 0;
	outpoints[4] = winw; 	outpoints[5] = winh;
	outpoints[6] = 0; 	outpoints[7] = winh;

	if (nd_getpersptransform(inpoints, outpoints, &persporig) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	if (nd_imgapplytransform(imgorig, &persporig, imgorig) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}
	
	if (nd_imgcreate(imginwin, winw, winh, imgorig->format) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	if (nd_imgcrop(imgorig, 0, 0, winw, winh, imginwin) < 0) {
		nd_imgdestroy(imginwin);
		
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	return 0;
}

int normalizenum(struct nd_image *imginwin)
{
	struct nd_matrix3 persp;
	struct nd_image imghsv;
	double inpoints[8];
	double outpoints[8];
	
	nd_imgcopy(imginwin, &imghsv);
	
	if (nd_imghsvval(&imghsv) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}
		
	if (ed_findborder(&imghsv, inpoints) < 0)
		return (-1);

	outpoints[0] = 0.0;		outpoints[1] = 0.0;
	outpoints[2] = imginwin->w;	outpoints[3] = 0.0;
	outpoints[4] = imginwin->w;	outpoints[5] = imginwin->h;
	outpoints[6] = 0.0;		outpoints[7] = imginwin->h;

	if (nd_getpersptransform(inpoints, outpoints, &persp) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	if (nd_imgapplytransform(imginwin, &persp, imginwin) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	return 0;
}

int detectedtofile(struct nd_image *img, struct nd_image *imgorig,
	struct nd_matrix3 *perspmat, struct hc_rect *r, int rc,
	const char *outputdir)
{
	char imgpath[255];
	char imgpathnormalized[255];
	static int foundn = 0;
	int rn;

	if (sprintf(imgpath, "%s/%d.png", outputdir, foundn) < 0) {
		fprintf(stderr, "sprintf: %s.\n", strerror(errno));
		return (-1);
	}

	if (nd_imgwrite(imgorig, imgpath) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	for (rn = 0; rn < rc; ++rn) {	
		struct nd_vector3 p[4];
		int winw, winh;
		double addtop, addbot;

		if (sprintf(imgpath, "%s/%d_%d.png", outputdir, foundn, rn)
			< 0) {
			fprintf(stderr, "sprintf: %s.\n", strerror(errno));
			return (-1);
		}

		if (sprintf(imgpathnormalized, "%s/%d_%dn.png",
			outputdir, foundn, rn) < 0) {
			fprintf(stderr, "sprintf: %s.\n", strerror(errno));
			return (-1);
		}

		addtop = (double) abs(r[rn].y0 - r[rn].y1) * 0.2;
		addbot = (double) abs(r[rn].y0 - r[rn].y1) * 0.4;

		recttoorig(r + rn, perspmat, addtop, addbot,
			(double) imgorig->w / IMGWIDTH, p);
		
		winw = (r[rn].x1 - r[rn].x0)
			* (double) imgorig->w / IMGWIDTH;
		winh = (r[rn].y1 - r[rn].y0 + addtop + addbot)
				* (double) imgorig->w / IMGWIDTH;

		if (winw <= imgorig->w && winh <= imgorig->h) {	
			struct nd_image imginwin;
		
			getorignum(imgorig, p, winw, winh, &imginwin);
	
			if (nd_imgwrite(&imginwin, imgpath) < 0) {
				nd_imgdestroy(&imginwin);
				
				fprintf(stderr, nd_geterrormessage());
				return (-1);
			}
		
			if (normalizenum(&imginwin) < 0)
				return 1;
			
			if (nd_imgwrite(&imginwin, imgpathnormalized) < 0) {
				nd_imgdestroy(&imginwin);
				
				fprintf(stderr, nd_geterrormessage());
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

int scanloop(struct nd_image *img, struct nd_image *imgorig,
	const char *hcpath, const char *outputdir)
{
	struct hcdata hcd;
	struct nd_matrix3 perspmat;

	if (hc_hcascaderead(&(hcd.hc), hcpath) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	hc_confbuild(&(hcd.scanconf),
		hcd.hc.ww, hcd.hc.wh, img->w, img->h, 0.9, 1, 1, threadcount);

	if (getperspmat(img->w, img->h, &perspmat) < 0)
		return (-1);
		
	while (1) {
		struct nd_image imgtr;
		struct hc_rect *r;
		int rc;

		if (nd_pslock(1) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}

		if (img->data == NULL)
			return 0;

struct timespec ts, te;
clock_gettime(CLOCK_REALTIME, &ts);

		if (nd_imgapplytransform(img, &perspmat, &imgtr) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}
/*
static int imgn = 0;
char a[255];
sprintf(a, "res1/%d.png", imgn++);
nd_imgwrite(&imgtr, a);
*/
		if (hc_imgpyramidscan(&(hcd.hc), &imgtr, &r, &rc,
			&(hcd.scanconf)) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}

clock_gettime(CLOCK_REALTIME, &te);
printf("%lf\n", (te.tv_sec - ts.tv_sec) + (te.tv_nsec - ts.tv_nsec) * 1e-9);

		if (rc) {
			printf("Detected.\n");
			if (detectedtofile(&imgtr, imgorig, &perspmat,
				r, rc, outputdir) < 0)
				fprintf(stderr, "Cannot write to file\n");
			free(r);
		}

		nd_imgdestroy(&imgtr);

		if (nd_psunlock(1) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}
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
		= rgbdata[y * rgblinesize + x * 3 + 1]
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
	struct nd_image *imgorig;
	struct playbackstate *pbs;
};

int decode(void *userdata)
{
	struct decodedata *ddata;
	int isdecoded;
	uint8_t *rgbdata;
	uint8_t *rgbdataorig;
	int rgblinesize;
	int rgblinesizeorig;
	int retval;

	ddata = userdata;

	isdecoded = 0;
/*
struct timespec ts, te;
clock_gettime(CLOCK_REALTIME, &ts);
*/

	do {
		if ((retval = readframe(ddata->av->s, ddata->av->vcodecc,
			ddata->av->vstreamid, ddata->av->frame,
			&(isdecoded))) < 0)
			return (-1);

		if (retval > 0)
			return 0;
	} while (!isdecoded);
/*
clock_gettime(CLOCK_REALTIME, &te);
printf("%lf\n", (te.tv_sec - ts.tv_sec) + (te.tv_nsec - ts.tv_nsec) * 1e-9);
*/
	if (frametorgb(ddata->av->frame, ddata->img->w, ddata->img->h,
		&rgbdata, &rgblinesize) < 0)
		return (-1);

	if (frametorgb(ddata->av->frame, ddata->imgorig->w, ddata->imgorig->h,
		&rgbdataorig, &rgblinesizeorig) < 0)
		return (-1);
	
	if (nd_pstrylock(0)) {
		rgbdatatoimg(rgbdata, rgblinesize, ddata->img);
		rgbdatatoimg(rgbdataorig, rgblinesizeorig, ddata->imgorig);
		
		if (nd_psunlock(0) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return (-1);
		}
	}

	free(rgbdata);
	free(rgbdataorig);


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
		outputline = ddata->img->data + i * surw;

		for (j = 0; j < surw; ++j) {
			inputline[j * 4 + 0]
			= inputline[j * 4 + 1]
			= inputline[j * 4 + 2]
			= inputline[j * 4 + 3] = 255.0 * (outputline[j] - 1.0);
		}
	}

	return 0;
}

int decodeloop(struct avdata *av, struct nd_image *img,
	struct nd_image *imgorig)
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
	ddata.imgorig = imgorig;
	ddata.pbs = &pbs;

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
	struct nd_image *imgorig;

	if (initavdata(&av, argv[1]) < 0)
		return 1;
	
	if (imgcreateshared(&img, IMGWIDTH, av.vcodecc->height
		* IMGWIDTH / av.vcodecc->width, ND_PF_GRAYSCALE) < 0)
		return 1;

	if (imgcreateshared(&imgorig, av.vcodecc->width, av.vcodecc->height,
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
		
		if (scanloop(img, imgorig, argv[2], argv[3]) < 0)
			return 1;
		
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

		if (decodeloop(&av, img, imgorig) < 0)
			return 1;
	
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

	if (imgdestroyshared(&imgorig) < 0)
		return 1;

	av_frame_unref(av.frame);
	avcodec_free_context(&(av.vcodecc));
	avformat_close_input(&(av.s));
	
	return 0;
}
