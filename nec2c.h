/*
 * nec2.h - header file for nec2
 */

#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>

#ifndef	TRUE
#define	TRUE	1
#endif

#ifndef	FALSE
#define	FALSE	0
#endif

/* commonly used complex constants */
#define	CPLX_00	(0.0+0.0fj)
#define	CPLX_01	(0.0+1.0fj)
#define	CPLX_10	(1.0+0.0fj)
#define	CPLX_11	(1.0+1.0fj)

/* common constants */
#define	PI	3.141592654
#define	POT	1.570796327
#define	TP	6.283185308
#define	PTP	.6283185308
#define	TPJ	(0.0+6.283185308fj)
#define PI8	25.13274123
#define PI10	31.41592654
#define	TA	1.745329252E-02
#define	TD	57.29577951
#define	ETA	376.73
#define	CVEL	299.8
#define	RETA	2.654420938E-3
#define	TOSP	1.128379167
#define ACCS	1.E-12
#define	SP	1.772453851
#define	FPI	12.56637062
#define	CCJ	(0.0-0.01666666667fj)
#define	CONST1	(0.0+4.771341189fj)
#define	CONST2	4.771341188
#define	CONST3	(0.0-29.97922085fj)
#define	CONST4	(0.0+188.365fj)
#define	GAMMA	.5772156649
#define C1	-.02457850915
#define C2	.3674669052
#define C3	.7978845608
#define P10	.0703125
#define P20	.1121520996
#define Q10	.125
#define Q20	.0732421875
#define P11	.1171875
#define P21	.1441955566
#define Q11	.375
#define Q21	.1025390625
#define POF	.7853981635
#define MAXH	20
#define CRIT	1.0E-4
#define NM	131072
#define NTS	4
#define	SMIN	1.e-3

/* Replaces the "10000" limit used to */
/* identify segment/patch connections */
#define	PCHCON  100000

/* carriage return and line feed */
#define	CR	0x0d
#define	LF	0x0a

/* max length of a line read from input file */
#define	LINE_LEN	132

/* version of fortran source for the -v option */
#define		version "nec2c v0.1"

/* nec2.c */
int 	main(int argc, char **argv);
void 	arc(int itg, int ns, double rada, double ang1, double ang2, double rad);
double 	atgn2(double x, double y);
void 	blckot(complex double *ar, int nunit, int ix1, int ix2, int nblks, int neof);
void 	blckin(complex double *ar, int nunit, int ix1, int ix2, int nblks, int neof);
void 	cabc(complex double *curx);
double 	cang(complex double z);
void 	cmset(int nrow, complex double *cm, double rkhx, int iexkx);
void 	cmss(int j1, int j2, int im1, int im2,
	complex double *cm, int nrow, int itrp);
void 	cmsw(int j1, int j2, int i1, int i2, complex double *cm,
	complex double *cw, int ncw, int nrow, int itrp);
void 	cmws(int j, int i1, int i2, complex double *cm, int nr,
	complex double *cw, int nw, int itrp);
void 	cmww(int j, int i1, int i2, complex double *cm, int nr,
	complex double *cw, int nw, int itrp);
void 	conect(int ignd);
void 	couple(complex double *cur, double wlam);
void 	datagn(void);
double 	db10(double x);
double 	db20(double x);
void 	efld(double xi, double yi, double zi, double ai, int ij);
void 	eksc(double s, double z, double rh, double xk, int ij,
	complex double *ezs, complex double *ers, complex double *ezc,
	complex double *erc, complex double *ezk, complex double *erk);
void 	ekscx(double bx, double s, double z, double rhx, double xk,
	int ij, int inx1, int inx2, complex double *ezs,
	complex double *ers, complex double *ezc, complex double *erc,
	complex double *ezk, complex double *erk);
void 	etmns(double p1, double p2, double p3, double p4, double p5,
	double p6, int ipr, complex double *e);
void 	factr(int n, complex double *a, int *ip, int ndim);
void 	factrs(int np, int nrow, complex double *a, int *ip);
complex double fbar(complex double p);
void 	fblock(int nrow, int ncol, int imax, int ipsym);
void 	ffld(double thet, double phi,
        complex double *eth, complex double *eph);
void 	fflds(double rox, double roy, double roz, complex double *scur,
	complex double *ex, complex double *ey, complex double *ez);
void 	gf(double zk, double *co, double *si);
void 	gfld(double rho, double phi, double rz, complex double *eth,
	complex double *epi, complex double *erd, complex double ux, int ksymp);
void 	gh(double zk, double *hr, double *hi);
void 	gwave(complex double *erv, complex double *ezv,
	complex double *erh, complex double *ezh, complex double *eph);
void 	gx(double zz, double rh, double xk,
	complex double *gz, complex double *gzp);
void 	gxx(double zz, double rh, double a, double a2, double xk,
	int ira, complex double *g1, complex double *g1p, complex double *g2,
	complex double *g2p, complex double *g3, complex double *gzp);
void 	helix(double s, double hl, double a1, double b1,
	double a2,double b2, double rad, int ns, int itg);
void 	hfk(double el1, double el2, double rhk,
	double zpkx, double *sgr, double *sgi);
void 	hintg(double xi, double yi, double zi);
void 	hsfld(double xi, double yi, double zi, double ai);
void 	hsflx(double s, double rh, double zpx, complex double *hpk,
	complex double *hps, complex double *hpc);
void 	intrp(double x, double y, complex double *f1,
	complex double *f2, complex double *f3, complex double *f4);
void 	intx(double el1, double el2, double b, int ij, double *sgr, double *sgi);
int 	isegno(int itagi, int mx);
void 	lfactr(complex double *a, int nrow, int ix1, int ix2, int *ip);
void 	load(int *ldtyp, int *ldtag, int *ldtagf, int *ldtagt,
	double *zlr, double *zli, double *zlc);
void 	lunscr(complex double *a, int nrow, int nop,
	int *ix, int *ip, int iu2, int iu3, int iu4);
void 	move(double rox, double roy, double roz, double xs,
	double ys, double zs, int its, int nrpt, int itgi);
void 	nefld(double xob, double yob, double zob, complex double *ex,
	complex double *ey, complex double *ez);
void 	netwk(complex double *cm, complex double *cmb, complex double *cmc,
	complex double *cmd, int *ip, complex double *einc);
void 	nfpat(void);
void 	nhfld(double xob, double yob, double zob, complex double *hx,
	complex double *hy, complex double *hz);
void 	patch(int nx, int ny, double ax1, double ay1, double az1,
	double ax2, double ay2, double az2, double ax3, double ay3,
	double az3, double ax4, double ay4, double az4);
void 	subph(int nx, int ny);
void 	pcint(double xi, double yi, double zi, double cabi,
	double sabi, double salpi, complex double *e);
void 	prnt(int in1, int in2, int in3, double fl1, double fl2,
	double fl3, double fl4, double fl5, double fl6, char *ia, int ichar);
void 	qdsrc(int is, complex double v, complex double *e);
void 	rdpat(void);
void 	readgm(char *gm, int *i1, int *i2, double *x1, double *y1,
	double *z1, double *x2, double *y2, double *z2, double *rad);
void 	readmn(char *gm, int *i1, int *i2, int *i3, int *i4, double *f1,
	double *f2, double *f3, double *f4, double *f5, double *f6);
void 	reflc(int ix, int iy, int iz, int itx, int nop);
void 	rom2(double a, double b, complex double *sum, double dmin);
void 	sbf(int i, int is, double *aa, double *bb, double *cc);
void 	sflds(double t, complex double *e);
void 	solgf(complex double *a, complex double *b, complex double *c,
	complex double *d, complex double *xy, int *ip, int np, int n1,
	int n, int mp, int m1, int m, int n1c, int n2c, int n2cz);
void 	solve(int n, complex double *a, int *ip, complex double *b, int ndim);
void 	solves(complex double *a, int *ip, complex double *b, int neq,
	int nrh, int np, int n, int mp, int m);
void 	tbf(int i, int icap);
void 	test(double f1r, double f2r, double *tr, double f1i,
	double f2i, double *ti, double dmin);
void 	trio(int j);
void 	unere(double xob, double yob, double zob);
void 	wire(double xw1, double yw1, double zw1, double xw2, double yw2,
	double zw2, double rad, double rdel, double rrad, int ns, int itg);
complex double zint(double sigl, double rolam);
int 	min(int a, int b);
/* misc.c */
void 	usage(void);
void 	abort_on_error(int why);
complex double cmplx(double a, double j);
void 	secnds(double *x);
int 	stop(int flag);
int 	load_line(char *buff, FILE *pfile);
void	mem_alloc( void **ptr, int req );
void	mem_realloc( void **ptr, int req );
void	free_ptr( void **ptr );
/* somnec.c */
void 	somnec(double epr, double sig, double fmhz);
void 	bessel(complex double z, complex double *j0, complex double *j0p);
void 	evlua(complex double *erv, complex double *ezv,
	complex double *erh, complex double *eph);
void 	gshank(complex double start, complex double dela, complex double *sum,
	int nans, complex double *seed, int ibk, complex double bk, complex double delb);
void 	hankel(complex double z, complex double *h0, complex double *h0p);
void 	lambda(double t, complex double *xlam, complex double *dxlam);
void 	rom1(int n, complex double *sum, int nx);
void 	saoa( double t, complex double *ans);

