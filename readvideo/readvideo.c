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
#include "bgmodel.h"
#include "hc_hcascade.h"
#include "hc_scanimgpyr.h"
#include "ed_findborder.h"

#include "gui.h"

#define IMGWIDTH 640
#define ADDHEIGHTTOP 0.2
#define ADDHEIGHTBOT 0.4
// #define ADDHEIGHTTOP 0.0
// #define ADDHEIGHTBOT 0.0

struct avdata {
	AVFormatContext *s;
	AVCodecContext *vcodecc;
	AVFrame *frame;
	int vstreamid;
};

struct carnumscanner {
	struct hc_hcascade hc;
	struct nd_image img;
	struct hc_scanconfig scanconf;
	struct nd_image *models;
	int modelscount;
};

struct playbackstate {
	double framerate;
	struct timespec tstart;
	int64_t timestamp;
};

struct haarcascadepath {
	char **hcpath;
	char **modelspath;
	int *modelscount;
	int hccount;
};

const int threadcount = 8;

struct nd_image *drawimg;

static void nd_safefree(void **p)
{
	if (*p != NULL) {
		free(*p);
		*p = NULL;
	}
}

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

	if (vcodec->capabilities & AV_CODEC_CAP_TRUNCATED)
		(*vcodecc)->flags |= AV_CODEC_FLAG_TRUNCATED;

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

/*
	inpoints[0] = 0; inpoints[1] = 0;
	inpoints[2] = 1; inpoints[3] = 0;
	inpoints[4] = 1; inpoints[5] = 1;
	inpoints[6] = 0; inpoints[7] = 1;

	outpoints[0] = 0; outpoints[1] = 0;
	outpoints[2] = 1; outpoints[3] = 0;
	outpoints[4] = 1; outpoints[5] = 1;
	outpoints[6] = 0; outpoints[7] = 1;
*/
	if (nd_getpersptransform(inpoints, outpoints, m) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

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

int rectgetorigpos(struct hc_rect *r, struct nd_matrix3 *perspmat,
	double addtop, double addbot, double origscaling,
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

	p[0].x *= origscaling; p[0].y *= origscaling;
	p[1].x *= origscaling; p[1].y *= origscaling;
	p[2].x *= origscaling; p[2].y *= origscaling;
	p[3].x *= origscaling; p[3].y *= origscaling;

	return 0;
}

int rectgetorignum(struct nd_image *imgorig, struct hc_rect *r, 
	struct nd_matrix3 *perspmat, double origscaling,
	struct nd_image *imginwin)
{
	struct nd_matrix3 persporig;
	double inpoints[8];
	double outpoints[8];

	struct nd_vector3 p[4];
	int winw, winh;
	double addtop, addbot;
	
	addtop = (double) abs(r->y0 - r->y1)
		* ADDHEIGHTTOP;
	addbot = (double) abs(r->y0 - r->y1)
		* ADDHEIGHTBOT;
	
	rectgetorigpos(r, perspmat, addtop, addbot, origscaling, p);

	winw = (r->x1 - r->x0)
		* origscaling;
	winh = (r->y1 - r->y0 + addtop + addbot)
		* origscaling;

	if (winw >= imgorig->w && winh >= imgorig->h)
		return (-1);

	inpoints[0] = p[0].x; inpoints[1] = p[0].y;
	inpoints[2] = p[1].x; inpoints[3] = p[1].y;
	inpoints[4] = p[2].x; inpoints[5] = p[2].y;
	inpoints[6] = p[3].x; inpoints[7] = p[3].y;

	outpoints[0] = 0; 	outpoints[1] = 0;
	outpoints[2] = winw; 	outpoints[3] = 0;
	outpoints[4] = winw; 	outpoints[5] = winh;
	outpoints[6] = 0; 	outpoints[7] = winh;

	if (nd_getpersptransform(inpoints, outpoints, &persporig) < 0)
		return (-1);

	if (nd_imgapplytransform(imgorig, &persporig, imgorig) < 0)
		return (-1);
	
	if (nd_imgcreate(imginwin, winw, winh, imgorig->format) < 0)
		return (-1);

	if (nd_imgcrop(imgorig, 0, 0, winw, winh, imginwin) < 0) {
		nd_imgdestroy(imginwin);
		
		return (-1);
	}

	return 0;
}

int restofile(struct nd_image *imgnum,
	double inpoints[8], const char *outputfile)
{
	struct nd_matrix3 borderpersp;
	double outpoints[8];
	
	outpoints[0] = 0.0;
	outpoints[1] = 0.0;
	outpoints[2] = imgnum->w;
	outpoints[3] = 0.0;
	outpoints[4] = imgnum->w;
	outpoints[5] = imgnum->h;
	outpoints[6] = 0.0;
	outpoints[7] = imgnum->h;

	if (nd_getpersptransform(inpoints, outpoints,
		&borderpersp) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	if (nd_imgapplytransform(imgnum, &borderpersp,
		imgnum) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return (-1);
	}

	nd_imgwrite(imgnum, outputfile);

	return 0;
}

int printelapsed(struct timespec *ts, struct timespec *te, FILE *file)
{
	double h, m, s;
	double td;
		
	td = (te->tv_sec - ts->tv_sec) + (te->tv_nsec - ts->tv_nsec) * 1e-9;

	h = td / 3600.0;
	m = (h - (int) h) * 60.0;
	s = (m - (int) m) * 60.0;
	
	fprintf(file, "detected: %02d:%02d:%02d\n",
		(int) h, (int) m, (int) s);

	return 0;
}

int getcarnum(struct nd_image *imgorig, struct nd_image *imgtr,
	struct nd_matrix3 perspmat, struct bg_rect *fgr,
	struct carnumscanner *hcd, int hccount, struct timespec *ts,
	const char *outputdir)
{
	int hcn;
	struct nd_image imgcropped;
	static int foundn = 0;

	if (nd_imgcreate(&imgcropped, fgr->x1 - fgr->x0,
		fgr->y1 - fgr->y0, ND_PF_RGB) < 0)
		return (-1);

	if (nd_imgcrop(imgtr, fgr->x0, fgr->y0,
		fgr->x1 - fgr->x0, fgr->y1 - fgr->y0, &imgcropped) < 0)
		return (-1);
		
	if (nd_imggrayscale(&imgcropped) < 0)
		return (-1);

	for (hcn = 0; hcn < hccount; ++hcn) {
		struct carnumscanner *curhc;
		struct hc_rect *r;
		int rc;
		int rn;

		curhc = hcd + hcn;

		if (hc_imgpyramidscan(&(curhc->hc), &imgcropped, &r, &rc,
			&(curhc->scanconf)) < 0)
			return (-1);

		for (rn = 0; rn < rc; ++rn) {
			struct nd_image imgnum;
			struct nd_image imghsv;
			double inpoints[8];
			struct timespec te;
			int bm;

			clock_gettime(CLOCK_REALTIME, &te);
			printelapsed(ts, &te, stdout);

			r[rn].x0 += fgr->x0;
			r[rn].y0 += fgr->y0;
			r[rn].x1 += fgr->x0;
			r[rn].y1 += fgr->y0;

			if (rectgetorignum(imgorig, r + rn, &perspmat,
				(double) imgorig->w / IMGWIDTH, &imgnum) < 0)
				return (-1);

///////////////////////////////////////////////////////////////////////////////
			char imgpath[1024];
			sprintf(imgpath, "%s/%d_0.png",
				outputdir, foundn);
			nd_imgwrite(&imgnum, imgpath);		
///////////////////////////////////////////////////////////////////////////////

			nd_imgcopy(&imgnum, &imghsv);

			if (nd_imghsvval(&imghsv) < 0)
				return (-1);

			if (ed_findborder(&imghsv, curhc->models,
				curhc->modelscount, inpoints, &bm) < 0)
				continue;

			if (sprintf(imgpath, "%s/%d_1.png", outputdir, foundn)
				< 0) {
				fprintf(stderr, "sprintf: %s.\n",
					strerror(errno));
				return (-1);
			}

///////////////////////////////////////////////////////////////////////////////
			if (restofile(&imgnum, inpoints, imgpath) < 0)
				return (-1);
///////////////////////////////////////////////////////////////////////////////
			
			++foundn;
		}

		if (rc)
			free(r);
	}

	nd_imgdestroy(&imgcropped);

	return 0;
}

int loadcarnumscanners(struct haarcascadepath *hcpath,
	struct carnumscanner **hcd, int *hccount)
{
	int hcn;
	int modeln;

	if ((*hcd = malloc(sizeof(struct carnumscanner) * hcpath->hccount))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		return (-1);	
	}

	*hccount = hcpath->hccount;

	modeln = 0;
	for (hcn = 0; hcn < hcpath->hccount; ++hcn) {
		struct carnumscanner *curhc;
		int i;

		curhc = *hcd + hcn;
	
		if (hc_hcascaderead(&(curhc->hc), hcpath->hcpath[hcn]) < 0)
			return (-1);

		if ((curhc->models = malloc(sizeof(struct nd_image)
			* hcpath->modelscount[hcn])) == NULL) {
			nd_seterrormessage(ND_MSGALLOCERROR, __func__);
			return (-1);	
		}

		curhc->modelscount = hcpath->modelscount[hcn];
		
		for (i = 0; i < hcpath->modelscount[hcn]; ++i) {
			if (nd_imgread(hcpath->modelspath[modeln],
				curhc->models + i) < 0)
				return (-1);

			if (nd_imggrayscale(curhc->models + i) < 0)
				return (-1);
		
			++modeln;
		}
	}

	return 0;
}

int initcarnumscanners(struct carnumscanner *hcd, int hccount,
	int imgw, int imgh)
{
	int hcn;
	
	for (hcn = 0; hcn < hccount; ++hcn) {
		struct carnumscanner *curhc;
		
		curhc = hcd + hcn;
		
		if (hc_confbuild(&(curhc->scanconf), curhc->hc.ww,
			curhc->hc.wh, imgw, imgh, 0.9, 1, 1, threadcount) < 0)
			return (-1);
	}

	return 0;
}

int nd_imgdrawrect(struct nd_image *img, struct bg_rect *r)
{
	double *pix;
	int xx, yy;
	
	for (xx = r->x0; xx < r->x1; ++xx) {
		pix = img->data + (r->y0 * img->w + xx) * 3;
		
		pix[0] = 0.0;
		pix[1] = 1.0;
		pix[2] = 0.0;
		
		pix = img->data + ((r->y1 - 1) * img->w + xx) * 3;
	
		pix[0] = 0.0;
		pix[1] = 1.0;
		pix[2] = 0.0;
	}

	for (yy = r->y0; yy < r->y1; ++yy) {
		pix = img->data + (yy * img->w + r->x0) * 3;
	
		pix[0] = 0.0;
		pix[1] = 1.0;
		pix[2] = 0.0;
		
		pix = img->data + (yy * img->w + (r->x1 - 1)) * 3;
	
		pix[0] = 0.0;
		pix[1] = 1.0;
		pix[2] = 0.0;
	}

	return 0;
}

int bgmdrawdetected(struct nd_image *imgtr, struct bg_rect *fgr, int fgrc,
	int bgmw, int bgmh, struct nd_image *drawimg)
{
	int i;

	memcpy(drawimg->data, imgtr->data,
		sizeof(double) * 3 * drawimg->w * drawimg->h);

	for (i = 0; i < fgrc; ++i) {
		struct bg_rect newr;

		newr.x0 = drawimg->w * fgr[i].x0 / bgmw;
		newr.y0 = drawimg->h * fgr[i].y0 / bgmh;
		newr.x1 = drawimg->w * fgr[i].x1 / bgmw;
		newr.y1 = drawimg->h * fgr[i].y1 / bgmh;
		
		nd_imgdrawrect(drawimg, &newr);
	}

	return 0;
}

int scanloop(struct nd_image *img, struct haarcascadepath *hcpath,
	const char *outputdir)
{
	struct carnumscanner *hcd;
	int hccount;
	struct bg_context ctx;
	struct nd_matrix3 perspmat;
	struct nd_image imgmini;
	struct nd_image imgorig;
	struct timespec ts;

	if (bg_initcontext(&ctx, 5, 0.3, 0.05, IMGWIDTH,
		img->h * IMGWIDTH / img->w, 30, 1) < 0)
		return (-1);

	if (loadcarnumscanners(hcpath, &hcd, &hccount) < 0)
		return (-1);

	if (initcarnumscanners(hcd, hccount, img->w, img->h) < 0)
		return (-1);

	if (getperspmat(IMGWIDTH, img->h * IMGWIDTH / img->w, &perspmat) < 0)
		return (-1);

	clock_gettime(CLOCK_REALTIME, &ts);
	
	while (1) {
		struct bg_rect *fgr;
		struct nd_image imgtr;
		int fgrc;
		int i;

		if (nd_pslock(1) < 0)
			return (-1);

		if (img->data == NULL)
			return 0;
 struct timespec ts, te;
 clock_gettime(CLOCK_REALTIME, &ts);
		if (nd_imgcopy(img, &imgorig) < 0)
			return (-1);

		if (nd_imgscalebilinear(&imgorig, (double) IMGWIDTH / imgorig.w,
			(double) IMGWIDTH / imgorig.w , &imgmini) < 0)
			return (-1);

		if (nd_imgapplytransform(&imgmini, &perspmat, &imgtr) < 0)
			return (-1);

		if (bg_putimgtocontext(&ctx, &imgtr, &fgr, &fgrc) < 0)
			return (-1);

 clock_gettime(CLOCK_REALTIME, &te);
 printf("%lf\n", (te.tv_sec - ts.tv_sec) + (te.tv_nsec - ts.tv_nsec) * 1e-9);
///////////////////////////////////////////////////////////////////////////////
		bgmdrawdetected(&imgtr, fgr, fgrc,
			ctx.bgm.w, ctx.bgm.h, drawimg);
///////////////////////////////////////////////////////////////////////////////

		if (nd_psunlock(1) < 0)
			return (-1);

		if (fgrc == 0)
			goto endstep;
	
		for (i = 0; i < fgrc; ++i)	
			getcarnum(&imgorig, &imgtr, perspmat, fgr + i,
				hcd, hccount, &ts, outputdir);

		free(fgr);

	endstep:
		nd_imgdestroy(&imgmini);
		nd_imgdestroy(&imgtr);
		nd_imgdestroy(&imgorig);
	}

	free(hcd);
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

/*
struct decodedata {
	struct avdata *av;
	struct nd_image *img;
	struct playbackstate *pbs;
};
*/

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

/*
int decodeloop(struct avdata *av, struct nd_image *img)
{
	struct playbackstate pbs;
	AVRational tmpr;

	tmpr = av->s->streams[av->vstreamid]->avg_frame_rate;
	if (pbstateinit(&pbs, (double) tmpr.num / tmpr.den) < 0)
		return 1;

	while (1) {
		int isdecoded;
		uint8_t *rgbdata;
		int rgblinesize;
		int retval;

		isdecoded = 0;

// struct timespec ts, te;
// clock_gettime(CLOCK_REALTIME, &ts);

		do {
			if ((retval = readframe(av->s, av->vcodecc,
				av->vstreamid, av->frame,
				&(isdecoded))) < 0)
				return (-1);

			if (retval > 0)
				return 0;
		} while (!isdecoded);

// clock_gettime(CLOCK_REALTIME, &te);
// printf("%lf\n", (te.tv_sec - ts.tv_sec) + (te.tv_nsec - ts.tv_nsec) * 1e-9);

		if (frametorgb(av->frame, img->w, img->h,
			&rgbdata, &rgblinesize) < 0)
			return (-1);
		
		if (nd_pstrylock(0)) {
			rgbdatatoimg(rgbdata, rgblinesize, img);
			
			if (nd_psunlock(0) < 0) {
				fprintf(stderr, nd_geterrormessage());
				return (-1);
			}
		}

		free(rgbdata);

		if (pbnexttimestamp(&pbs) < 0)
			return (-1);
	}

	fprintf(stderr, "Unexpexted loop quit\n");
	return (-1);
}
*/

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
	
	for (i = 0; i < MIN(surh, ddata->drawimg->h); ++i) {
		unsigned char *inputline;
		double *outputline;

		inputline = surdata + i * surw * 4;
		outputline = ddata->drawimg->data + i * ddata->drawimg->w * 3;
		
		for (j = 0; j < MIN(surw, ddata->drawimg->w); ++j) {
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

///////////////////////////////////////////////////////////////////////////////
	if (pbstateinit(&pbs, (double) tmpr.num / tmpr.den) < 0)
		return 1;
///////////////////////////////////////////////////////////////////////////////
	
	
	initgui(&guidata, drawimg->w, drawimg->h);

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

int readhcpath(const char *path, struct haarcascadepath *hcpath)
{
	FILE *file;
	int totalmodelscount;
	int i, ii;
	
	if ((file = fopen(path, "r")) == NULL) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		goto fopenerror;
	}

	if (fscanf(file, "%d", &(hcpath->hccount)) < 0) {
		nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
		goto scanhccounterror;
	}
	
	if ((hcpath->hcpath = malloc(sizeof(char *) * hcpath->hccount))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto mallochcpatherror;
	}

	if ((hcpath->modelscount = malloc(sizeof(int) * hcpath->hccount))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto mallocmodelscounterror;
	}

	totalmodelscount = 0;
	for (i = 0; i < hcpath->hccount; ++i) {
		if (fscanf(file, "%d", hcpath->modelscount + i) < 0) {
			nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
			goto scanmodelscounterror;
		}
		
		totalmodelscount += hcpath->modelscount[i];
	}

	if ((hcpath->modelspath = malloc(sizeof(char *) * totalmodelscount))
		== NULL) {
		nd_seterrormessage(ND_MSGALLOCERROR, __func__);
		goto mallocmodelpatherror;
	}

	fscanf(file, "\n");
	for (i = 0; i < hcpath->hccount; ++i)  {
		size_t strsz;
		char **path;

		path = hcpath->hcpath + i;
		strsz = 0;
		if (getline(path, &strsz, file) <= 0) {
			nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
			goto getlineerror;
		}
		
		(*path)[strlen(*path) - 1] = '\0';
	}

	for (i = 0; i < totalmodelscount; ++i) {
		size_t strsz;
		char **path;
		int res;

		path = hcpath->modelspath + i;
		strsz = 0;
		
		res = getline(path, &strsz, file);
		
		if (res <= 0 && res != EOF) {
			strerror(errno);
			nd_seterrormessage(ND_MSGFILEIOERROR, __func__);
			goto getlineerror;
		}
		
		(*path)[strlen(*path) - 1] = '\0';
	}

	return 0;

getlineerror:
	nd_safefree((void **) &(hcpath->hcpath));

	for (ii = 0; ii < i - 1; ++i)
		free(hcpath->modelspath + ii);
	
	nd_safefree((void **) &(hcpath->modelspath));
mallocmodelpatherror:
scanmodelscounterror:
mallocmodelscounterror:
	nd_safefree((void **) &(hcpath->hcpath));
mallochcpatherror:
scanhccounterror:
fopenerror:
	return (-1);
}

int main(int argc, char **argv)
{
	int chpid;
	struct avdata av;
	struct nd_image *img;
	struct haarcascadepath hcpath;
	
	if (initavdata(&av, argv[1]) < 0)
		return 1;
	
	if (imgcreateshared(&img, av.vcodecc->width, av.vcodecc->height,
		ND_PF_RGB) < 0)
		return 1;

	if (imgcreateshared(&drawimg, IMGWIDTH, img->h * IMGWIDTH / img->w,
		ND_PF_RGB) < 0)
		return 1;

	if (readhcpath(argv[2], &hcpath) < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}

	if (nd_psinitprefork() < 0) {
		fprintf(stderr, nd_geterrormessage());
		return 1;
	}
	
	if ((chpid = fork()) == 0) {
		if (nd_psinitpostfork(1) < 0) {
			fprintf(stderr, nd_geterrormessage());
			return 1;
		}
		
		if (scanloop(img, &hcpath, argv[3]) < 0)
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

		if (decodeloop(&av, img) < 0)
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

	av_frame_unref(av.frame);
	avcodec_free_context(&(av.vcodecc));
	avformat_close_input(&(av.s));
	
	return 0;
}
