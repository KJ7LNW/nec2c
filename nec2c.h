#ifndef	NEC2C_H
#define	NEC2C_H 1

#include <complex.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
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
#define PI		3.141592654
#define	POT		1.570796327
#define	TP		6.283185308
#define	PTP		.6283185308
#define	TPJ		(0.0+6.283185308fj)
#define PI8		25.13274123
#define PI10	31.41592654
#define	TA		1.745329252E-02
#define	TD		57.29577951
#define	ETA		376.73
#define	CVEL	299.8
#define	RETA	2.654420938E-3
#define	TOSP	1.128379167
#define ACCS	1.E-12
#define	SP		1.772453851
#define	FPI		12.56637062
#define	CCJ		(0.0-0.01666666667fj)
#define	CONST1	(0.0+4.771341189fj)
#define	CONST2	4.771341188
#define	CONST3	(0.0-29.97922085fj)
#define	CONST4	(0.0+188.365fj)
#define	GAMMA	.5772156649
#define C1		-.02457850915
#define C2		.3674669052
#define C3		.7978845608
#define P10		.0703125
#define P20		.1121520996
#define Q10		.125
#define Q20		.0732421875
#define P11		.1171875
#define P21		.1441955566
#define Q11		.375
#define Q21		.1025390625
#define POF		.7853981635
#define MAXH	20
#define CRIT	1.0E-4
#define NM		131072
#define NTS		4
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
#define		version "nec2c 0.7"

/*** Structs encapsulating global ("common") variables ***/
/* common  /crnt/ */
typedef struct
{
  long double
	*air,	/* Ai/lambda, real part */
	*aii,	/* Ai/lambda, imaginary part */
	*bir,	/* Bi/lambda, real part */
	*bii,	/* Bi/lambda, imaginary part */
	*cir,	/* Ci/lambda, real part */
	*cii;	/* Ci/lambda, imaginary part */

  complex long double *cur; /* Amplitude of basis function */

} crnt_t;

/* common  /data/ (geometry data) */
typedef struct
{
  int
	n,		/* Number of wire segments */
	np,		/* Number of wire segments in symmetry cell */
	m,		/* Number of surface patches */
	mp,		/* Number of surface patches in symmetry cell */
	npm,	/* = n+m  */
	np2m,	/* = n+2m */
	np3m,	/* = n+3m */
	ipsym,	/* Symmetry flag */
	*icon1, /* Segments end 1 connection */
	*icon2,	/* Segments end 2 connection */
	*itag;	/* Segments tag number */

  /* Wire segment data */
  long double
	*x1, *y1, *z1,	/* End 1 coordinates of wire segments */
	*x2, *y2, *z2,	/* End 2 coordinates of wire segments */
	*x, *y, *z,		/* Coordinates of segment centers */
	*si, *bi,		/* Length and radius of segments  */
	*cab,			/* cos(a)*cos(b) */
	*sab,			/* cos(a)*sin(b) */
	*salp,			/* Z component - sin(a) */

	/* Surface patch data */
	*px, *py, *pz,		/* Coordinates of patch center */
	*t1x, *t1y, *t1z,	/* Coordinates of t1 vector */
	*t2x, *t2y, *t2z,	/* Coordinates of t2 vector */
	*pbi,				/* Patch surface area */
	*psalp,				/* Z component - sin(a) */

	/* Wavelength in meters */
	wlam;

} data_t;

/* common  /dataj/ */
typedef struct
{
  int
	iexk,
	ind1,
	indd1,
	ind2,
	indd2,
	ipgnd;

  long double
	s,
	b,
	xj,
	yj,
	zj,
	cabj,
	sabj,
	salpj,
	rkh,
	t1xj,
	t1yj,
	t1zj,
	t2xj,
	t2yj,
	t2zj;

  complex long double
	exk,
	eyk,
	ezk,
	exs,
	eys,
	ezs,
	exc,
	eyc,
	ezc;

} dataj_t;

/* common  /fpat/ */
typedef struct
{
  int
	near,
	nfeh,
	nrx,
	nry,
	nrz,
	nth,
	nph,
	ipd,
	iavp,
	inor,
	iax,
	ixtyp;

  long double
	thets,
	phis,
	dth,
	dph,
	rfld,
	gnor,
	clt,
	cht,
	epsr2,
	sig2,
	xpr6,
	pinr,
	pnlr,
	ploss,
	xnr,
	ynr,
	znr,
	dxnr,
	dynr,
	dznr;

} fpat_t;

/*common  /ggrid/ */
typedef struct
{
  int
	nxa[3],
	nya[3];

  long double
	dxa[3],
	dya[3],
	xsa[3],
	ysa[3];

  complex long double
	epscf,
	*ar1,
	*ar2,
	*ar3;

} ggrid_t;

/* common  /gnd/ */
typedef struct
{
  int
	ksymp,	/* Ground flag */
	ifar,	/* Int flag in RP card, for far field calculations */
	iperf,	/* Type of ground flag */
	nradl;	/* Number of radials in ground screen */

  long double
	t2,		/* Const for radial wire ground impedance */
	cl,		/* Distance in wavelengths of cliff edge from origin */
	ch,		/* Cliff height in wavelengths */
	scrwl,	/* Wire length in radial ground screen normalized to w/length */
	scrwr;	/* Radius of wires in screen in wavelengths */

  complex long double
	zrati,	/* Ground medium [Er-js/wE0]^-1/2 */
	zrati2,	/* As above for 2nd ground medium */
	t1,		/* Const for radial wire ground impedance */
	frati;	/* (k1^2-k2^2)/(k1^2+k2^2), k1=w(E0Mu0)^1/2, k1=k2/ZRATI */

} gnd_t;

/* common  /gwav/ */
typedef struct
{
  long double
	r1,		/* Distance from current element to point where field is evaluated  */
	r2,		/* Distance from image of element to point where field is evaluated */
	zmh,	/* Z-Z', Z is height of field evaluation point */
	zph;	/* Z+Z', Z' is height of current element */

  complex long double
	u,		/* (Er-jS/WE0)^-1/2 */
	u2,		/* u^2 */
	xx1,	/* G1*exp(jkR1.r[i])  */
	xx2;	/* G2*exp(jkR2.r'[i]) */

} gwav_t;

/* common  /incom/ */
typedef struct
{
  int isnor;

  long double
	xo,
	yo,
	zo,
	sn,
	xsn,
	ysn;

} incom_t;

/* common  /matpar/ (matrix parameters) */
typedef struct
{
  int
	icase,	/* Storage mode of primary matrix */
	npblk,	/* Num of blocks in first (NBLOKS-1) blocks */
	nlast,	/* Num of blocks in last block */
	imat;	/* Storage reserved in CM for primary NGF matrix A */

} matpar_t;

/* common  /netcx/ */
typedef struct
{
  int
	masym,	/* Matrix symmetry flags */
	neq,
	npeq,
	neq2,
	nonet,	/* Number of two-port networks */
	ntsol,	/* "Network equations are solved" flag */
	nprint,	/* Print control flag */
	*iseg1,	/* Num of seg to which port 1 of network is connected */
	*iseg2,	/* Num of seg to which port 2 of network is connected */
	*ntyp;	/* Type of networks */

  long double
	*x11r,	/* Real and imaginary parts of network impedances */
	*x11i,
	*x12r,
	*x12i,
	*x22r,
	*x22i,
	pin,	/* Total input power from sources */
	pnls;	/* Power lost in networks */

  complex long double zped;

} netcx_t;

/* common  /plot/ */
typedef struct
{
  int
	/* Plot control flags */
	iplp1,
	iplp2,
	iplp3,
	iplp4;

} plot_t;

/* common  /save/ */
typedef struct
{
  int *ip;	/* Vector of indices of pivot elements used to factor matrix */

  long double
	epsr,	/* Relative dielectric constant of ground */
	sig,	/* Conductivity of ground */
	scrwlt,	/* Length of radials in ground screen approximation */
	scrwrt,	/* Radius of wires in ground screen approximation */
	fmhz;	/* Frequency in MHz */

} save_t;

/* common  /segj/ */
typedef struct
{
  int
	*jco,	/* Stores connection data */
	jsno,	/* Total number of entries in ax, bx, cx */
	maxcon; /* Max. no. connections */

  long double
	*ax, *bx, *cx;	/* Store constants A, B, C used in current expansion */

} segj_t;

/* common  /smat/ */
typedef struct
{
  int nop; /* My addition */

  complex long double *ssx;

} smat_t;

/* common  /tmi/ */
typedef struct
{
  int ij;

  long double
	zpk,
	rkb2;

} tmi_t;

/*common  /tmh/ */
typedef struct
{
  long double
	zpka,
	rhks;

} tmh_t;

/* common  /vsorc/ */
typedef struct
{
  int
	*isant,	/* Num of segs on which an aplied field source is located */
	*ivqd,	/* Num of segs on which a current-slope discontinuity source is located */
	*iqds,	/* Same as above (?) */
	nsant,	/* Number of applied field voltage sources */
	nvqd,	/* Number of applied current-slope discontinuity sources */
	nqds;	/* Same as above (?) */

  complex long double
	*vqd,	/* Voltage of applied-current slope discontinuity sources */
	*vqds,	/* Same as above (?) */
	*vsant;	/* Voltages of applied field voltage sources */

} vsorc_t;

/* common  /yparm/ */
typedef struct
{
  int
	ncoup,	/* Num of segs between which coupling will be computed */
	icoup,	/* Num of segs in the coupling array that have been excited */
	*nctag,	/* Tag number of segments */
	*ncseg;	/* Num of segs in set of segs that have same tag number */

  complex long double
	*y11a,	/* Self admittance of segments */
	*y12a;	/* Mutual admittances stored in order 1,2 1,3 2,3 2,4 etc */

} yparm_t;

/* common  /zload/ */
typedef struct
{
  int nload;	/* Number of loading networks */

  complex long double *zarray;	/* = Zi/(Di/lambda) */

} zload_t;

/* Returns the complex long double of the arguments */
#define cmplx(r, i) ((r)+(i)*CPLX_01)

/*------------------------------------------------------------------------*/

/* Function prototypes produced by cproto */
/* calculations.c */
void cabc(complex long double *curx);
void couple(complex long double *cur, long double wlam);
void load(int *ldtyp, int *ldtag, int *ldtagf, int *ldtagt, long double *zlr, long double *zli, long double *zlc);
void gf(long double zk, long double *co, long double *si);
long double db10(long double x);
long double db20(long double x);
void intrp(long double x, long double y, complex long double *f1, complex long double *f2, complex long double *f3, complex long double *f4);
void intx(long double el1, long double el2, long double b, int ij, long double *sgr, long double *sgi);
int min(int a, int b);
void test(long double f1r, long double f2r, long double *tr, long double f1i, long double f2i, long double *ti, long double dmin);
void sbf(int i, int is, long double *aa, long double *bb, long double *cc);
void tbf(int i, int icap);
void trio(int j);
void zint(long double sigl, long double rolam, complex long double *zt);
long double cang(complex long double z);
/* fields.c */
void efld(long double xi, long double yi, long double zi, long double ai, int ij);
void eksc(long double s, long double z, long double rh, long double xk, int ij, complex long double *ezs, complex long double *ers, complex long double *ezc, complex long double *erc, complex long double *ezk, complex long double *erk);
void ekscx(long double bx, long double s, long double z, long double rhx, long double xk, int ij, int inx1, int inx2, complex long double *ezs, complex long double *ers, complex long double *ezc, complex long double *erc, complex long double *ezk, complex long double *erk);
void gh(long double zk, long double *hr, long double *hi);
void gwave(complex long double *erv, complex long double *ezv, complex long double *erh, complex long double *ezh, complex long double *eph);
void gx(long double zz, long double rh, long double xk, complex long double *gz, complex long double *gzp);
void gxx(long double zz, long double rh, long double a, long double a2, long double xk, int ira, complex long double *g1, complex long double *g1p, complex long double *g2, complex long double *g2p, complex long double *g3, complex long double *gzp);
void hfk(long double el1, long double el2, long double rhk, long double zpkx, long double *sgr, long double *sgi);
void hintg(long double xi, long double yi, long double zi);
void hsfld(long double xi, long double yi, long double zi, long double ai);
void hsflx(long double s, long double rh, long double zpx, complex long double *hpk, complex long double *hps, complex long double *hpc);
void nefld(long double xob, long double yob, long double zob, complex long double *ex, complex long double *ey, complex long double *ez);
void nfpat(void);
void nhfld(long double xob, long double yob, long double zob, complex long double *hx, complex long double *hy, complex long double *hz);
void pcint(long double xi, long double yi, long double zi, long double cabi, long double sabi, long double salpi, complex long double *e);
void unere(long double xob, long double yob, long double zob);
/* geometry.c */
void arc(int itg, int ns, long double rada, long double ang1, long double ang2, long double rad);
void conect(int ignd);
void datagn(void);
void helix(long double s, long double hl, long double a1, long double b1, long double a2, long double b2, long double rad, int ns, int itg);
int isegno(int itagi, int mx);
void move(long double rox, long double roy, long double roz, long double xs, long double ys, long double zs, int its, int nrpt, int itgi);
void patch(int nx, int ny, long double ax1, long double ay1, long double az1, long double ax2, long double ay2, long double az2, long double ax3, long double ay3, long double az3, long double ax4, long double ay4, long double az4);
void subph(int nx, int ny);
void readgm(char *gm, int *i1, int *i2, long double *x1, long double *y1, long double *z1, long double *x2, long double *y2, long double *z2, long double *rad);
void reflc(int ix, int iy, int iz, int itx, int nop);
void wire(long double xw1, long double yw1, long double zw1, long double xw2, long double yw2, long double zw2, long double rad, long double rdel, long double rrad, int ns, int itg);
/* ground.c */
void rom2(long double a, long double b, complex long double *sum, long double dmin);
void sflds(long double t, complex long double *e);
/* input.c */
void qdsrc(int is, complex long double v, complex long double *e);
void readmn(char *gm, int *i1, int *i2, int *i3, int *i4, long double *f1, long double *f2, long double *f3, long double *f4, long double *f5, long double *f6);
/* main.c */
int main(int argc, char **argv);
void Null_Pointers(void);
void prnt(int in1, int in2, int in3, long double fl1, long double fl2, long double fl3, long double fl4, long double fl5, long double fl6, char *ia, int ichar);
/* matrix.c */
void cmset(int nrow, complex long double *cm, long double rkhx, int iexkx);
void cmss(int j1, int j2, int im1, int im2, complex long double *cm, int nrow, int itrp);
void cmsw(int j1, int j2, int i1, int i2, complex long double *cm, complex long double *cw, int ncw, int nrow, int itrp);
void cmws(int j, int i1, int i2, complex long double *cm, int nr, complex long double *cw, int nw, int itrp);
void cmww(int j, int i1, int i2, complex long double *cm, int nr, complex long double *cw, int nw, int itrp);
void etmns(long double p1, long double p2, long double p3, long double p4, long double p5, long double p6, int ipr, complex long double *e);
void factr(int n, complex long double *a, int *ip, int ndim);
void factrs(int np, int nrow, complex long double *a, int *ip);
void fblock(int nrow, int ncol, int imax, int ipsym);
void solve(int n, complex long double *a, int *ip, complex long double *b, int ndim);
void solves(complex long double *a, int *ip, complex long double *b, int neq, int nrh, int np, int n, int mp, int m);
/* misc.c */
void usage(void);
void abort_on_error(int why);
void secnds(long double *x);
int stop(int flag);
int load_line(char *buff, FILE *pfile);
void mem_alloc(void **ptr, int req);
void mem_realloc(void **ptr, int req);
void free_ptr(void **ptr);
/* network.c */
void netwk(complex long double *cm, int *ip, complex long double *einc);
/* radiation.c */
void ffld(long double thet, long double phi, complex long double *eth, complex long double *eph);
void fflds(long double rox, long double roy, long double roz, complex long double *scur, complex long double *ex, complex long double *ey, complex long double *ez);
void gfld(long double rho, long double phi, long double rz, complex long double *eth, complex long double *epi, complex long double *erd, complex long double ux, int ksymp);
void rdpat(void);
/* somnec.c */
void somnec(long double epr, long double sig, long double fmhz);
void bessel(complex long double z, complex long double *j0, complex long double *j0p);
void evlua(complex long double *erv, complex long double *ezv, complex long double *erh, complex long double *eph);
void fbar(complex long double p, complex long double *r);
void gshank(complex long double start, complex long double dela, complex long double *sum, int nans, complex long double *seed, int ibk, complex long double bk, complex long double delb);
void hankel(complex long double z, complex long double *h0, complex long double *h0p);
void lambda(long double t, complex long double *xlam, complex long double *dxlam);
void rom1(int n, complex long double *sum, int nx);
void saoa(long double t, complex long double *ans);
#endif

