/******* Translated to the C language by N. Kyriazis  20 Aug 2003 *******/
/*									*/
/* Program NEC(input,tape5=input,output,tape11,tape12,tape13,tape14,	*/
/* tape15,tape16,tape20,tape21)						*/
/*									*/
/* Numerical Electromagnetics Code (NEC2)  developed at Lawrence	*/
/* Livermore lab., Livermore, CA.  (contact G. Burke at 415-422-8414	*/
/* for problems with the NEC code. For problems with the vax implem- 	*/
/* entation, contact J. Breakall at 415-422-8196 or E. Domning at 415 	*/
/* 422-5936) 								*/
/* file created 4/11/80. 						*/
/*									*/
/*                ***********Notice********** 				*/
/* This computer code material was prepared as an account of work 	*/
/* sponsored by the United States government.  Neither the United 	*/
/* States nor the United States Department Of Energy, nor any of 	*/
/* their employees, nor any of their contractors, subcontractors, 	*/
/* or their employees, makes any warranty, express or implied, or	*/
/* assumes any legal liability or responsibility for the accuracy, 	*/
/* completeness or usefulness of any information, apparatus, product 	*/
/* or process disclosed, or represents that its use would not infringe 	*/
/* privately-owned rights. 						*/
/*									*/
/************************************************************************/

#include "nec2c.h"

/*** common data are implemented as global variables ***/
/* common  /data/ */
int n, np, m, mp, ipsym, npm, np2m, np3m; /* n+m,n+2m,n+3m */
int *icon1, *icon2, *itag;
double *x, *y, *z, *si, *bi;
double *x2, *y2, *z2, *cab, *sab, *salp;
double *t1x, *t1y, *t1z, *t2x, *t2y, *t2z;
double *px, *py, *pz, *pbi, *psalp;
double wlam;

/* common  /cmb/ */
complex double *cm;

/* common  /matpar/ */
int icase, npblk, nlast;
int imat, nbbx, npbx, nlbx, nbbl, npbl, nlbl;

/* common  /save/ */
int *ip;
double epsr, sig, scrwlt, scrwrt, fmhz;

/* common  /crnt/ */
double *air, *aii, *bir, *bii, *cir, *cii;
complex double *cur;

/* common  /gnd/ */
int ksymp, ifar, iperf, nradl;
double t2, cl, ch, scrwl, scrwr;
complex double zrati, zrati2, t1, frati;

/* common  /zload/ */
int nload;
complex double *zarray;

/* common  /yparm/ */
int ncoup, icoup, *nctag, *ncseg;
complex double *y11a, *y12a;

/* common  /segj/ */
int *jco, jsno, nscon, maxcon; /* Max. no. connections */
double *ax, *bx, *cx;

/* common  /vsorc/ */
int *ivqd, *isant, *iqds, nvqd, nsant, nqds;
complex double *vqd, *vqds, *vsant;

/* common  /netcx/ */
int masym, neq, npeq, neq2, nonet, ntsol, nprint;
int *iseg1, *iseg2, *ntyp;
double *x11r, *x11i, *x12r;
double *x12i, *x22r, *x22i;
double pin, pnls;
complex double zped;

/* common  /fpat/ */
int near, nfeh, nrx, nry, nrz, nth, nph, ipd, iavp, inor, iax, ixtyp;
double thets, phis, dth, dph, rfld, gnor, clt, cht, epsr2, sig2;
double xpr6, pinr, pnlr, ploss, xnr, ynr, znr, dxnr, dynr, dznr;

/*common  /ggrid/ */
extern int nxa[3], nya[3];
extern double dxa[3], dya[3], xsa[3], ysa[3];
extern complex double epscf, *ar1, *ar2, *ar3;

/* common  /gwav/ */
double r1, r2, zmh, zph;
complex double u, u2, xx1, xx2;

/* common  /plot/ */
int iplp1, iplp2, iplp3, iplp4;

/* common  /dataj/ */
int iexk, ind1, indd1, ind2, indd2, ipgnd;
double s, b, xj, yj, zj, cabj, sabj, salpj, rkh;
double t1xj, t1yj, t1zj, t2xj, t2yj, t2zj;
complex double  exk, eyk, ezk, exs, eys, ezs, exc, eyc, ezc;

/* common  /smat/ */
int nop; /* My addition */
complex double *ssx;

/* common  /incom/ */
int isnor;
double xo, yo, zo, sn, xsn, ysn;

/* common  /tmi/ */
int ija; /* changed to ija to avoid conflict */
double zpk, rkb2;

/*common  /tmh/ */
double zpka, rhks;

/* pointers to input/output files */
FILE *input_fp=NULL, *output_fp=NULL, *plot_fp=NULL;

/* signal handler */
static void sig_handler( int signal );

/*-------------------------------------------------------------------*/

int main( int argc, char **argv )
{
  char infile[81] = "", otfile[81] = "";
  char ain[3], line_buf[81];

  /* input card mnemonic list */
  /* "XT" stands for "exit", added for testing */
#define CMD_NUM  20
  char *atst[CMD_NUM] =
  {
    "FR", "LD", "GN", "EX", "NT", "TL", \
    "XQ", "GD", "RP", "NX", "PT", "KH", \
    "NE", "NH", "PQ", "EK", "CP", "PL", \
    "EN", "WG"
  };

  char *hpol[3] = { "LINEAR", "RIGHT", "LEFT" };
  char *pnet[3] = { "        ", "STRAIGHT", " CROSSED" };

  int *ldtyp, *ldtag, *ldtagf, *ldtagt;
  int ifrtmw, ifrtmp, mpcnt, ib11=0, ic11=0, id11=0, ix11, igo, nfrq;
  int iexk, jump, iptflg, iptflq, iped, iflow, itmp1, iresrv;
  int itmp3, itmp2, itmp4, nthi=0, nphi=0, iptag=0, iptagf=0, iptagt=0;
  int iptaq=0, iptaqf=0, iptaqt=0, nphic=0, inc=0;
  int i, j, itmp5, nthic=0, mhz=0, ifrq=0, isave=0;

  int
    igox,       /* used in place of "igo" in freq loop */
    next_job,   /* start next job (next sructure) flag */
    idx,        /* general purpose index    */
    ain_num,    /* ain mnemonic as a number */
    jmp_iloop,  /* jump to input loop flag  */
    jmp_floop=0,/* jump to freq. loop flag  */
    mreq;       /* Size req. for malloc's   */

  double *zlr, *zli, *zlc, *fnorm;
  double *xtemp, *ytemp, *ztemp, *sitemp, *bitemp;
  double fmhz1, rkh, tmp1, delfrq=0., tmp2, tmp3, tmp4, tmp5, tmp6;
  double xpr1=0., xpr2=0., xpr3=0., xpr4=0., xpr5=0.;
  double zpnorm=0., thetis=0., phiss=0., extim;
  double tim1, tim, tim2, etha, fr, fr2, cmag, ph, ethm, ephm, epha;
  complex double eth, eph, curi, ex, ey, ez, epsc;

  /* getopt() variables */
  extern char *optarg;
  extern int optind, opterr, optopt;
  int option;

  /*** signal handler related code ***/
  /* new and old actions for sigaction() */
  struct sigaction sa_new, sa_old;


  /* initialize new actions */
  sa_new.sa_handler = sig_handler;
  sigemptyset( &sa_new.sa_mask );
  sa_new.sa_flags = 0;

  /* register function to handle signals */
  sigaction( SIGINT,  &sa_new, &sa_old );
  sigaction( SIGSEGV, &sa_new, 0 );
  sigaction( SIGFPE,  &sa_new, 0 );
  sigaction( SIGTERM, &sa_new, 0 );
  sigaction( SIGABRT, &sa_new, 0 );

  /*** command line arguments handler ***/
  if( argc == 1 )
  {
    usage();
    exit(-1);
  }

  /* process command line options */
  while( (option = getopt(argc, argv, "i:o:hv") ) != -1 )
  {
    switch( option )
    {
      case 'i' : /* specify input file name */
	if( strlen(optarg) > 75 )
	  abort_on_error(-1);
	strcpy( infile, optarg );
	break;

      case 'o' : /* specify output file name */
	if( strlen(optarg) > 75 )
	  abort_on_error(-2);
	strcpy( otfile, optarg );
	break;

      case 'h' : /* print usage and exit */
	usage();
	exit(0);

      case 'v' : /* print nec2c version */
	puts( version );
	exit(0);

      default: /* print usage and exit */
	usage();
	exit(-1);

    } /* end of switch( option ) */

  } /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */

  /*** open input file ***/
  if( (input_fp = fopen(infile, "r")) == NULL )
  {
    char mesg[88] = "nec2c: ";

    strcat( mesg, infile );
    perror( mesg );
    exit(-1);
  }

  /* make an output file name if not */
  /* specified by user on invocation */
  if( strlen( otfile ) == 0 )
  {
    /* strip file name extension if there is one */
    idx = 0;
    while( (infile[++idx] != '.') && (infile[idx] != '\0') );
    infile[idx] = '\0';

    /* make the output file name */
    strcpy( otfile, infile );
  }

  /* add extension */
  strcat( otfile, ".out" );

  /* open output file */
  if( (output_fp = fopen(otfile, "w")) == NULL )
  {
    char mesg[88] = "nec2c: ";

    strcat( mesg, otfile );
    perror( mesg );
    exit(-1);
  }

  /*** here we had code to read interactively input/output ***/
  /*** file names. this is done non-interactively above.   ***/

  secnds( &extim );

  /* Null buffer pointers */
  /* type int */
  icon1 = icon2 = ncseg = nctag = ivqd = isant = iqds = NULL;
  itag = ip = ldtyp = ldtag = ldtagf = ldtagt = jco = NULL;
  /* type double */
  air = aii = bir = bii = cir = cii = zlr = zli = zlc = fnorm = NULL;
  ax = bx = cx = xtemp = ytemp = ztemp = sitemp = bitemp = NULL;
  x = y = z = si = bi = x2 = y2 = z2 = cab = sab = salp = NULL;
  t1x = t1y = t1z = t2x = t2y = t2z = px = py = pz = pbi = psalp = NULL;
  /* type complex double */
  ar1 = ar2 = ar3 = cur = cm = zarray = NULL;
  y11a = y12a = vqd = vqds = vsant = ssx = NULL;

  /* Allocate some buffers */
  mem_alloc( (void *)&ar1, sizeof(complex double)*11*10*4 );
  mem_alloc( (void *)&ar2, sizeof(complex double)*17*5*4 );
  mem_alloc( (void *)&ar3, sizeof(complex double)*9*8*4 );


  /* l_1: */
  /* main execution loop, exits at various points */
  /* depending on error conditions or end of jobs */
  while( TRUE )
  {
    ifrtmw=0;
    ifrtmp=0;

    /* print the nec2c header to output file */
    fprintf( output_fp,	"\n\n\n"
	"                              "
	" __________________________________________\n"
	"                              "
	"|                                          |\n"
	"                              "
	"|  NUMERICAL ELECTROMAGNETICS CODE (nec2c) |\n"
	"                              "
	"|   Translated to 'C' in Double Precision  |\n"
	"                              "
	"|__________________________________________|\n" );

    /* read a line from input file */
    if( load_line(line_buf, input_fp) == EOF )
      abort_on_error(-3);

    /* separate card's id mnemonic */
    strncpy( ain, line_buf, 2 );
    ain[2] = '\0';

    /* If its an "XT" card, exit (used for debugging) */
    if( strcmp(ain, "XT") == 0 )
    {
      fprintf( stderr,
	  "\nnec2c: Exiting after an \"XT\" command in main()\n" );
      fprintf( output_fp,
	  "\n\n  nec2c: Exiting after an \"XT\" command in main()" );
      stop(0);
    }

    /* if its a "cm" or "ce" card start reading comments */
    if( (strcmp(ain, "CM") == 0) ||
	(strcmp(ain, "CE") == 0) )
    {
      fprintf( output_fp, "\n\n\n"
	  "                               "
	  "---------------- COMMENTS ----------------\n" );

      /* write comment to output file */
      fprintf( output_fp,
	  "                              %s\n",
	  &line_buf[2] );

      /* Keep reading till a non "CM" card */
      while( strcmp(ain, "CM") == 0 )
      {
	/* read a line from input file */
	if( load_line(line_buf, input_fp) == EOF )
	  abort_on_error(-3);

	/* separate card's id mnemonic */
	strncpy( ain, line_buf, 2 );
	ain[2] = '\0';

	/* write comment to output file */
	fprintf( output_fp,
	    "                              %s\n",
	    &line_buf[2] );

      } /* while( strcmp(ain, "CM") == 0 ) */

      /* no "ce" card at end of comments */
      if( strcmp(ain, "CE") != 0 )
      {
	fprintf( output_fp,
	    "\n\n  ERROR: INCORRECT LABEL FOR A COMMENT CARD" );
	abort_on_error(-4);
      }

    } /* if( strcmp(ain, "CM") == 0 ... */
    else
      rewind( input_fp );

    /* Free some buffer pointers.
     * These are allocated by realloc()
     * so they need to be free()'d
     * before reallocation for a new job
     */
    free_ptr( (void *)&itag );
    free_ptr( (void *)&fnorm );
    free_ptr( (void *)&ldtyp );
    free_ptr( (void *)&ldtag );
    free_ptr( (void *)&ldtagf );
    free_ptr( (void *)&ldtagt );
    free_ptr( (void *)&zlr );
    free_ptr( (void *)&zli );
    free_ptr( (void *)&zlc );
    free_ptr( (void *)&jco );
    free_ptr( (void *)&ax );
    free_ptr( (void *)&bx );
    free_ptr( (void *)&cx );
    free_ptr( (void *)&ivqd );
    free_ptr( (void *)&iqds );
    free_ptr( (void *)&vqd );
    free_ptr( (void *)&vqds );
    free_ptr( (void *)&isant );
    free_ptr( (void *)&vsant );
    free_ptr( (void *)&x );
    free_ptr( (void *)&y );
    free_ptr( (void *)&z );
    free_ptr( (void *)&x2 );
    free_ptr( (void *)&y2 );
    free_ptr( (void *)&z2 );
    free_ptr( (void *)&px );
    free_ptr( (void *)&py );
    free_ptr( (void *)&pz );
    free_ptr( (void *)&t1x );
    free_ptr( (void *)&t1y );
    free_ptr( (void *)&t1z );
    free_ptr( (void *)&t2x );
    free_ptr( (void *)&t2y );
    free_ptr( (void *)&t2z );
    free_ptr( (void *)&si );
    free_ptr( (void *)&bi );
    free_ptr( (void *)&cab );
    free_ptr( (void *)&sab );
    free_ptr( (void *)&salp );
    free_ptr( (void *)&pbi );
    free_ptr( (void *)&psalp );

    /* initializations etc from original fortran code */
    mpcnt=0;
    imat=0;

    /* set up geometry data in subroutine datagn */
    datagn();
    iflow=1;

    /* Allocate some buffers */
    mreq = npm * sizeof(double);
    mem_alloc( (void *)&air, mreq );
    mem_alloc( (void *)&aii, mreq );
    mem_alloc( (void *)&bir, mreq );
    mem_alloc( (void *)&bii, mreq );
    mem_alloc( (void *)&cir, mreq );
    mem_alloc( (void *)&cii, mreq );
    mem_alloc( (void *)&xtemp,  mreq );
    mem_alloc( (void *)&ytemp,  mreq );
    mem_alloc( (void *)&ztemp,  mreq );
    mem_alloc( (void *)&sitemp, mreq );
    mem_alloc( (void *)&bitemp, mreq );

    mreq = np2m * sizeof(int);
    mem_alloc( (void *)&ip, mreq );

    mreq = np3m * sizeof( complex double);
    mem_alloc( (void *)&cur, mreq );

    /* Matrix parameters */
    if( imat == 0)
    {
      neq= n+2*m;
      neq2=0;
      ib11=0;
      ic11=0;
      id11=0;
      ix11=0;
    }

    fprintf( output_fp, "\n\n\n" );

    /* default values for input parameters and flags */
    npeq= np+2*mp;
    iplp1=0;
    iplp2=0;
    iplp3=0;
    iplp4=0;
    igo=1;
    nfrq=1;
    rkh=1.;
    iexk=0;
    ixtyp=0;
    nload=0;
    nonet=0;
    near=-1;
    iptflg=-2;
    iptflq=-1;
    ifar=-1;
    zrati=CPLX_10;
    iped=0;
    ncoup=0;
    icoup=0;
    fmhz= CVEL;
    ksymp=1;
    nradl=0;
    iperf=0;

    /* l_14: */

    /* main input section, exits at various points */
    /* depending on error conditions or end of job */
    next_job = FALSE;
    while( ! next_job )
    {
      jmp_iloop = FALSE;

      /* main input section - standard read statement - jumps */
      /* to appropriate section for specific parameter set up */
      readmn( ain, &itmp1, &itmp2, &itmp3, &itmp4,
	  &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6 );

      /* If its an "XT" card, exit */
      if( strcmp(ain, "XT" ) == 0 )
      {
	fprintf( stderr,
	    "\nnec2c: Exiting after an \"XT\" command in main()\n" );
	fprintf( output_fp,
	    "\n\n  nec2c: Exiting after an \"XT\" command in main()" );
	stop(0);
      }

      mpcnt++;
      fprintf( output_fp,
	  "\n  DATA CARD No: %3d "
	  "%s %3d %5d %5d %5d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E",
	  mpcnt, ain, itmp1, itmp2, itmp3, itmp4,
	  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );

      /* identify card id mnemonic (except "ce" and "cm") */
      for( ain_num = 0; ain_num < CMD_NUM; ain_num++ )
	if( strncmp( ain, atst[ain_num], 2) == 0 )
	  break;

      /* take action according to card id mnemonic */
      switch( ain_num )
      {
	case 0: /* "fr" card, frequency parameters */

	  ifrq= itmp1;
	  nfrq= itmp2;
	  if( nfrq == 0)
	    nfrq=1;
	  fmhz= tmp1;
	  delfrq= tmp2;
	  if( iped == 1)
	    zpnorm=0.;
	  igo=1;
	  iflow=1;

	  continue; /* continue card input loop */

	case 1: /* "ld" card, loading parameters */
	  {
	    int idx;

	    if( iflow != 3 )
	    {
	      iflow=3;
	      /* Free loading buffers */
	      nload=0;
	      free_ptr( (void *)&ldtyp );
	      free_ptr( (void *)&ldtag );
	      free_ptr( (void *)&ldtagf );
	      free_ptr( (void *)&ldtagt );
	      free_ptr( (void *)&zlr );
	      free_ptr( (void *)&zli );
	      free_ptr( (void *)&zlc );

	      if( igo > 2 )
		igo=2;
	      if( itmp1 == -1 )
		continue; /* continue card input loop */
	    }

	    /* Reallocate loading buffers */
	    nload++;
	    idx = nload * sizeof(int);
	    mem_realloc( (void *)&ldtyp,  idx );
	    mem_realloc( (void *)&ldtag,  idx );
	    mem_realloc( (void *)&ldtagf, idx );
	    mem_realloc( (void *)&ldtagt, idx );
	    idx = nload * sizeof(double);
	    mem_realloc( (void *)&zlr, idx );
	    mem_realloc( (void *)&zli, idx );
	    mem_realloc( (void *)&zlc, idx );

	    idx = nload-1;
	    ldtyp[idx]= itmp1;
	    ldtag[idx]= itmp2;
	    if( itmp4 == 0)
	      itmp4= itmp3;
	    ldtagf[idx]= itmp3;
	    ldtagt[idx]= itmp4;

	    if( itmp4 < itmp3 )
	    {
	      fprintf( output_fp,
		  "\n\n  DATA FAULT ON LOADING CARD No: %d: ITAG "
		  "STEP1: %d IS GREATER THAN ITAG STEP2: %d",
		  nload, itmp3, itmp4 );
	      stop(-1);
	    }

	    zlr[idx]= tmp1;
	    zli[idx]= tmp2;
	    zlc[idx]= tmp3;
	  }

	  continue; /* continue card input loop */

	case 2: /* "gn" card, ground parameters under the antenna */

	  iflow=4;

	  if( igo > 2)
	    igo=2;

	  if( itmp1 == -1 )
	  {
	    ksymp=1;
	    nradl=0;
	    iperf=0;
	    continue; /* continue card input loop */
	  }

	  iperf= itmp1;
	  nradl= itmp2;
	  ksymp=2;
	  epsr= tmp1;
	  sig= tmp2;

	  if( nradl != 0)
	  {
	    if( iperf == 2)
	    {
	      fprintf( output_fp,
		  "\n\n  RADIAL WIRE G.S. APPROXIMATION MAY "
		  "NOT BE USED WITH SOMMERFELD GROUND OPTION" );
	      stop(-1);
	    }

	    scrwlt= tmp3;
	    scrwrt= tmp4;
	    continue; /* continue card input loop */
	  }

	  epsr2= tmp3;
	  sig2= tmp4;
	  clt= tmp5;
	  cht= tmp6;

	  continue; /* continue card input loop */

	case 3: /* "ex" card, excitation parameters */

	  if( iflow != 5)
	  {
	    /* Free vsource buffers */
	    free_ptr( (void *)&ivqd );
	    free_ptr( (void *)&iqds );
	    free_ptr( (void *)&vqd );
	    free_ptr( (void *)&vqds );
	    free_ptr( (void *)&isant );
	    free_ptr( (void *)&vsant );

	    nsant=0;
	    nvqd=0;
	    iped=0;
	    iflow=5;
	    if( igo > 3)
	      igo=3;
	  }

	  masym= itmp4/10;
	  if( (itmp1 == 0) || (itmp1 == 5) )
	  {
	    ixtyp= itmp1;
	    ntsol=0;

	    if( ixtyp != 0)
	    {
	      nvqd++;
	      mem_realloc( (void *)&ivqd, nvqd * sizeof(int) );
	      mem_realloc( (void *)&iqds, nvqd * sizeof(int) );
	      mem_realloc( (void *)&vqd,  nvqd * sizeof(complex double) );
	      mem_realloc( (void *)&vqds, nvqd * sizeof(complex double) );

	      {
		int indx = nvqd-1;

		ivqd[indx]= isegno( itmp2, itmp3);
		vqd[indx]= cmplx( tmp1, tmp2);
		if( cabs( vqd[indx]) < 1.e-20)
		  vqd[indx] = CPLX_10;

		iped= itmp4- masym*10;
		zpnorm= tmp3;
		if( (iped == 1) && (zpnorm > 0.0) )
		  iped=2;
		continue; /* continue card input loop */
	      }

	    } /* if( ixtyp != 0) */

	    nsant++;
	    mem_realloc( (void *)&isant, nsant * sizeof(int) );
	    mem_realloc( (void *)&vsant, nsant * sizeof(complex double) );

	    {
	      int indx = nsant-1;

	      isant[indx]= isegno( itmp2, itmp3);
	      vsant[indx]= cmplx( tmp1, tmp2);
	      if( cabs( vsant[indx]) < 1.e-20)
		vsant[indx] = CPLX_10;

	      iped= itmp4- masym*10;
	      zpnorm= tmp3;
	      if( (iped == 1) && (zpnorm > 0.0) )
		iped=2;
	      continue; /* continue card input loop */
	    }

	  } /* if( (itmp1 <= 0) || (itmp1 == 5) ) */

	  if( (ixtyp == 0) || (ixtyp == 5) )
	    ntsol=0;

	  ixtyp= itmp1;
	  nthi= itmp2;
	  nphi= itmp3;
	  xpr1= tmp1;
	  xpr2= tmp2;
	  xpr3= tmp3;
	  xpr4= tmp4;
	  xpr5= tmp5;
	  xpr6= tmp6;
	  nsant=0;
	  nvqd=0;
	  thetis= xpr1;
	  phiss= xpr2;

	  continue; /* continue card input loop */

	case 4: case 5: /* "nt" & "tl" cards, network parameters */
	  {
	    int idx;

	    if( iflow != 6)
	    {
	      nonet=0;
	      ntsol=0;
	      iflow=6;

	      /* Free network buffers */
	      free_ptr( (void *)&ntyp );
	      free_ptr( (void *)&iseg1 );
	      free_ptr( (void *)&iseg2 );
	      free_ptr( (void *)&x11r );
	      free_ptr( (void *)&x11i );
	      free_ptr( (void *)&x12r );
	      free_ptr( (void *)&x12i );
	      free_ptr( (void *)&x22r );
	      free_ptr( (void *)&x22i );

	      if( igo > 3)
		igo=3;

	      if( itmp2 == -1 )
		continue; /* continue card input loop */
	    }

	    /* Re-allocate network buffers */
	    nonet++;
	    idx = nonet * sizeof(int);
	    mem_realloc( (void *)&ntyp, idx );
	    mem_realloc( (void *)&iseg1, idx );
	    mem_realloc( (void *)&iseg2, idx );
	    idx = nonet * sizeof(double);
	    mem_realloc( (void *)&x11r, idx );
	    mem_realloc( (void *)&x11i, idx );
	    mem_realloc( (void *)&x12r, idx );
	    mem_realloc( (void *)&x12i, idx );
	    mem_realloc( (void *)&x22r, idx );
	    mem_realloc( (void *)&x22i, idx );

	    idx = nonet-1;
	    if( ain_num == 4 )
	      ntyp[idx]=1;
	    else
	      ntyp[idx]=2;

	    iseg1[idx]= isegno( itmp1, itmp2);
	    iseg2[idx]= isegno( itmp3, itmp4);
	    x11r[idx]= tmp1;
	    x11i[idx]= tmp2;
	    x12r[idx]= tmp3;
	    x12i[idx]= tmp4;
	    x22r[idx]= tmp5;
	    x22i[idx]= tmp6;

	    if( (ntyp[idx] == 1) || (tmp1 > 0.) )
	      continue; /* continue card input loop */

	    ntyp[idx]=3;
	    x11r[idx]=- tmp1;

	    continue; /* continue card input loop */
	  }

	case 6: /* "xq" execute card - calc. including radiated fields */

	  if( ((iflow == 10) && (itmp1 == 0)) ||
	      ((nfrq  ==  1) && (itmp1 == 0) && (iflow > 7)) )
	    continue; /* continue card input loop */

	  if( itmp1 == 0)
	  {
	    if( iflow > 7)
	      iflow=11;
	    else
	      iflow=7;
	  }
	  else
	  {
	    ifar=0;
	    rfld=0.;
	    ipd=0;
	    iavp=0;
	    inor=0;
	    iax=0;
	    nth=91;
	    nph=1;
	    thets=0.;
	    phis=0.;
	    dth=1.0;
	    dph=0.;

	    if( itmp1 == 2)
	      phis=90.;

	    if( itmp1 == 3)
	    {
	      nph=2;
	      dph=90.;
	    }

	  } /* if( itmp1 == 0) */

	  break;

	case 7: /* "gd" card, ground representation */

	  epsr2= tmp1;
	  sig2= tmp2;
	  clt= tmp3;
	  cht= tmp4;
	  iflow=9;

	  continue; /* continue card input loop */

	case 8: /* "rp" card, standard observation angle parameters */

	  ifar= itmp1;
	  nth= itmp2;
	  nph= itmp3;

	  if( nth == 0)
	    nth=1;
	  if( nph == 0)
	    nph=1;

	  ipd= itmp4/10;
	  iavp= itmp4- ipd*10;
	  inor= ipd/10;
	  ipd= ipd- inor*10;
	  iax= inor/10;
	  inor= inor- iax*10;

	  if( iax != 0)
	    iax=1;
	  if( ipd != 0)
	    ipd=1;
	  if( (nth < 2) || (nph < 2) || (ifar == 1) )
	    iavp=0;

	  thets= tmp1;
	  phis= tmp2;
	  dth= tmp3;
	  dph= tmp4;
	  rfld= tmp5;
	  gnor= tmp6;
	  iflow=10;

	  break;

	case 9: /* "nx" card, do next job */
	  next_job = TRUE;
	  continue; /* continue card input loop */

	case 10: /* "pt" card, print control for current */

	  iptflg= itmp1;
	  iptag= itmp2;
	  iptagf= itmp3;
	  iptagt= itmp4;

	  if( (itmp3 == 0) && (iptflg != -1) )
	    iptflg=-2;
	  if( itmp4 == 0)
	    iptagt= iptagf;

	  continue; /* continue card input loop */

	case 11: /* "kh" card, matrix integration limit */

	  rkh= tmp1;
	  if( igo > 2)
	    igo=2;
	  iflow=1;

	  continue; /* continue card input loop */

	case 12: case 13:  /* "ne"/"nh" cards, near field calculation parameters */

	  if( ain_num == 13 )
	    nfeh=1;
	  else
	    nfeh=0;

	  if( (iflow == 8) && (nfrq != 1) )
	  {
	    fprintf( output_fp,
		"\n\n  WHEN MULTIPLE FREQUENCIES ARE REQUESTED, "
		"ONLY ONE NEAR FIELD CARD CAN BE USED -"
		"\n  LAST CARD READ WILL BE USED" );
	  }

	  near= itmp1;
	  nrx= itmp2;
	  nry= itmp3;
	  nrz= itmp4;
	  xnr= tmp1;
	  ynr= tmp2;
	  znr= tmp3;
	  dxnr= tmp4;
	  dynr= tmp5;
	  dznr= tmp6;
	  iflow=8;

	  if( nfrq != 1)
	    continue; /* continue card input loop */

	  break;

	case 14: /* "pq" card, write control for charge */

	  iptflq= itmp1;
	  iptaq= itmp2;
	  iptaqf= itmp3;
	  iptaqt= itmp4;

	  if( (itmp3 == 0) && (iptflq != -1) )
	    iptflq=-2;
	  if( itmp4 == 0)
	    iptaqt= iptaqf;

	  continue; /* continue card input loop */

	case 15: /* "ek" card,  extended thin wire kernel option */

	  iexk=1;
	  if( itmp1 == -1)
	    iexk=0;
	  if( igo > 2)
	    igo=2;
	  iflow=1;

	  continue; /* continue card input loop */

	case 16: /* "cp" card, maximum coupling between antennas */

	  if( iflow != 2)
	  {
	    ncoup=0;
	    free_ptr( (void *)&nctag );
	    free_ptr( (void *)&ncseg );
	    free_ptr( (void *)&y11a );
	    free_ptr( (void *)&y12a );
	  }

	  icoup=0;
	  iflow=2;

	  if( itmp2 == 0)
	    continue; /* continue card input loop */

	  ncoup++;
	  mem_realloc( (void *)&nctag, (ncoup) * sizeof(int) );
	  mem_realloc( (void *)&ncseg, (ncoup) * sizeof(int) );
	  nctag[ncoup-1]= itmp1;
	  ncseg[ncoup-1]= itmp2;

	  if( itmp4 == 0)
	    continue; /* continue card input loop */

	  ncoup++;
	  mem_realloc( (void *)&nctag, (ncoup) * sizeof(int) );
	  mem_realloc( (void *)&ncseg, (ncoup) * sizeof(int) );
	  nctag[ncoup-1]= itmp3;
	  ncseg[ncoup-1]= itmp4;

	  continue; /* continue card input loop */

	case 17: /* "pl" card, plot flags */

	  iplp1= itmp1;
	  iplp2= itmp2;
	  iplp3= itmp3;
	  iplp4= itmp4;

	  if( plot_fp == NULL )
	  {
	    char plotfile[81];

	    /* Make a plot file name */
	    strcpy( plotfile, infile );
	    strcat( plotfile, ".plt" );

	    /* Open plot file */
	    if( (plot_fp = fopen(plotfile, "w")) == NULL )
	    {
	      char mesg[88] = "nec2c: ";

	      strcat( mesg, plotfile );
	      perror( mesg );
	      exit(-1);
	    }
	  }

	  continue; /* continue card input loop */

	case 19: /* "wg" card, not supported */
	  abort_on_error(-5);

	default:
	  if( ain_num != 18 )
	  {
	    fprintf( output_fp,
		"\n\n  FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION" );
	    stop(-1);
	  }

	  /******************************************************
	   *** normal exit of nec2c when all jobs complete ok ***
	   ******************************************************/

	  /* time the process */
	  secnds( &tmp1 );
	  tmp1 -= extim;
	  fprintf( output_fp, "\n\n  TOTAL RUN TIME: %d msec", (int)tmp1 );
	  stop(0);

      } /* switch( ain_num ) */

      /**************************************
       *** end of the main input section. ***
       *** beginning of frequency do loop ***
       **************************************/

      /* Allocate to normalization buffer */
      {
	int mreq1, mreq2;

	mreq1 = mreq2 = 0;
	if( iped )
	  mreq1 = 4*nfrq * sizeof(double);
	if( iptflg >= 2 )
	  mreq2 = nthi*nphi * sizeof(double);

	if( (mreq1 > 0) || (mreq2 > 0) )
	{
	  if( mreq1 > mreq2 )
	    mem_alloc( (void *)&fnorm, mreq1 );
	  else
	    mem_alloc( (void *)&fnorm, mreq2 );
	}
      }

      /* igox is used in place of "igo" in the   */
      /* freq loop. below is a special igox case */
      if( ((ain_num == 6) || (ain_num == 8)) && (igo == 5) )
	igox = 6;
      else
	igox = igo;

      switch( igox )
      {
	case 1: /* label 41 */
	  /* Memory allocation for primary interacton matrix. */
	  iresrv = np2m * (np+2*mp);
	  mem_alloc( (void *)&cm, iresrv * sizeof(complex double) );

	  /* Memory allocation for symmetry array */
	  nop = neq/npeq;
	  mem_alloc( (void *)&ssx, nop*nop * sizeof( complex double) );

	  mhz=1;

	  if( (n != 0) && (ifrtmw != 1) )
	  {
	    ifrtmw=1;
	    for( i = 0; i < n; i++ )
	    {
	      xtemp[i]= x[i];
	      ytemp[i]= y[i];
	      ztemp[i]= z[i];
	      sitemp[i]= si[i];
	      bitemp[i]= bi[i];
	    }
	  }

	  if( (m != 0) && (ifrtmp != 1) )
	  {
	    ifrtmp=1;
	    for( i = 0; i < m; i++ )
	    {
	      j = i+n;
	      xtemp[j]= px[i];
	      ytemp[j]= py[i];
	      ztemp[j]= pz[i];
	      bitemp[j]= pbi[i];
	    }
	  }

	  fmhz1= fmhz;

	  /* irngf is not used (NGF function not implemented) */
	  if( imat == 0)
	    fblock( npeq, neq, iresrv, ipsym);

	  /* label 42 */
	  /* frequency do loop */
	  do
	  {
	    jmp_floop = FALSE;

	    if( mhz != 1)
	    {
	      if( ifrq == 1)
		fmhz *= delfrq;
	      else
		fmhz += delfrq;
	    }

	    fr= fmhz/ CVEL;
	    wlam= CVEL/ fmhz;
	    fprintf( output_fp, "\n\n\n"
		"                               "
		"--------- FREQUENCY --------\n"
		"                                "
		"FREQUENCY :%11.4E MHz\n"
		"                                "
		"WAVELENGTH:%11.4E Mtr", fmhz, wlam );

	    fprintf( output_fp, "\n\n"
		"                        "
		"APPROXIMATE INTEGRATION EMPLOYED FOR SEGMENTS \n"
		"                        "
		"THAT ARE MORE THAN %.3f WAVELENGTHS APART", rkh );

	    if( iexk == 1)
	      fprintf( output_fp, "\n"
		  "                        "
		  "THE EXTENDED THIN WIRE KERNEL WILL BE USED" );

	    /* frequency scaling of geometric parameters */
	    if( n != 0)
	    {
	      for( i = 0; i < n; i++ )
	      {
		x[i]= xtemp[i]* fr;
		y[i]= ytemp[i]* fr;
		z[i]= ztemp[i]* fr;
		si[i]= sitemp[i]* fr;
		bi[i]= bitemp[i]* fr;
	      }
	    }

	    if( m != 0)
	    {
	      fr2= fr* fr;
	      for( i = 0; i < m; i++ )
	      {
		j = i+n;
		px[i]= xtemp[j]* fr;
		py[i]= ytemp[j]* fr;
		pz[i]= ztemp[j]* fr;
		pbi[i]= bitemp[j]* fr2;
	      }
	    }

	    igo = 2;

	    /* label 46 */
	    case 2: /* structure segment loading */
	    fprintf( output_fp, "\n\n\n"
		"                          "
		"------ STRUCTURE IMPEDANCE LOADING ------" );

	    if( nload != 0)
	      load( ldtyp, ldtag, ldtagf, ldtagt, zlr, zli, zlc );

	    if( nload == 0 )
	      fprintf( output_fp, "\n"
		  "                                 "
		  "THIS STRUCTURE IS NOT LOADED" );

	    fprintf( output_fp, "\n\n\n"
		"                            "
		"-------- ANTENNA ENVIRONMENT --------" );

	    if( ksymp != 1)
	    {
	      frati=CPLX_10;

	      if( iperf != 1)
	      {
		if( sig < 0.)
		  sig=- sig/(59.96*wlam);

		epsc= cmplx( epsr, -sig*wlam*59.96);
		zrati=1./ csqrt( epsc);
		u= zrati;
		u2= u* u;

		if( nradl != 0)
		{
		  scrwl= scrwlt/ wlam;
		  scrwr= scrwrt/ wlam;
		  t1= CPLX_01*2367.067/ (double)nradl;
		  t2= scrwr* (double)nradl;

		  fprintf( output_fp, "\n"
		      "                            "
		      "RADIAL WIRE GROUND SCREEN\n"
		      "                            "
		      "%d WIRES\n"
		      "                            "
		      "WIRE LENGTH: %8.2f METERS\n"
		      "                            "
		      "WIRE RADIUS: %10.3E METERS",
		      nradl, scrwlt, scrwrt );

		  fprintf( output_fp, "\n"
		      "                            "
		      "MEDIUM UNDER SCREEN -" );

		} /* if( nradl != 0) */

		if( iperf != 2)
		  fprintf( output_fp, "\n"
		      "                            "
		      "FINITE GROUND - REFLECTION COEFFICIENT APPROXIMATION" );
		else
		{
		  somnec( epsr, sig, fmhz );
		  frati=( epsc-1.)/( epsc+1.);
		  if( cabs(( epscf- epsc)/ epsc) >= 1.0e-3 )
		  {
		    fprintf( output_fp,
			"\n ERROR IN GROUND PARAMETERS -"
			"\n COMPLEX DIELECTRIC CONSTANT FROM FILE IS: %12.5E%+12.5Ej"
			"\n                                REQUESTED: %12.5E%+12.5Ej",
			creal(epscf), cimag(epscf), creal(epsc), cimag(epsc) );
		    stop(-1);
		  }

		  fprintf( output_fp, "\n"
		      "                            "
		      "FINITE GROUND - SOMMERFELD SOLUTION" );

		} /* if( iperf != 2) */

		fprintf( output_fp, "\n"
		    "                            "
		    "RELATIVE DIELECTRIC CONST: %.3f\n"
		    "                            "
		    "CONDUCTIVITY: %10.3E MHOS/METER\n"
		    "                            "
		    "COMPLEX DIELECTRIC CONSTANT: %11.4E%+11.4Ej",
		    epsr, sig, creal(epsc), cimag(epsc) );

	      } /* if( iperf != 1) */
	      else
		fprintf( output_fp, "\n"
		    "                            "
		    "PERFECT GROUND" );

	    } /* if( ksymp != 1) */
	    else
	      fprintf( output_fp, "\n"
		  "                            "
		  "FREE SPACE" );

	    /* label 50 */
	    /* fill and factor primary interaction matrix */
	    secnds( &tim1 );
	    cmset( neq, cm, rkh, iexk );
	    secnds( &tim2 );
	    tim= tim2- tim1;
	    factrs( npeq, neq, cm, ip );
	    secnds( &tim1 );
	    tim2= tim1- tim2;
	    fprintf( output_fp, "\n\n\n"
		"                             "
		"---------- MATRIX TIMING ----------\n"
		"                               "
		"FILL: %d msec  FACTOR: %d msec",
		(int)tim, (int)tim2 );

	    igo=3;
	    ntsol=0;

	    /* label 53 */
	    case 3: /* excitation set up (right hand side, -e inc.) */

	    nthic=1;
	    nphic=1;
	    inc=1;
	    nprint=0;

	    /* l_54 */
	    do
	    {
	      if( (ixtyp != 0) && (ixtyp != 5) )
	      {
		if( (iptflg <= 0) || (ixtyp == 4) )
		  fprintf( output_fp, "\n\n\n"
		      "                             "
		      "---------- EXCITATION ----------" );

		tmp5= TA* xpr5;
		tmp4= TA* xpr4;

		if( ixtyp == 4)
		{
		  tmp1= xpr1/ wlam;
		  tmp2= xpr2/ wlam;
		  tmp3= xpr3/ wlam;
		  tmp6= xpr6/( wlam* wlam);

		  fprintf( output_fp, "\n"
		      "                                  "
		      "    CURRENT SOURCE\n"
		      "                     -- POSITION (METERS) -- "
		      "      ORIENTATION (DEG)\n"
		      "                     X          Y          Z "
		      "      ALPHA        BETA   DIPOLE MOMENT\n"
		      "               %10.5f %10.5f %10.5f "
		      " %7.2f     %7.2f    %8.3f",
		      xpr1, xpr2, xpr3, xpr4, xpr5, xpr6 );
		}
		else
		{
		  tmp1= TA* xpr1;
		  tmp2= TA* xpr2;
		  tmp3= TA* xpr3;
		  tmp6= xpr6;

		  if( iptflg <= 0)
		    fprintf( output_fp,
			"\n  PLANE WAVE - THETA: %7.2f deg, PHI: %7.2f deg,"
			" ETA=%7.2f DEG, TYPE - %s  AXIAL RATIO: %6.3f",
			xpr1, xpr2, xpr3, hpol[ixtyp-1], xpr6 );

		} /* if( ixtyp == 4) */

	      } /* if( (ixtyp  != 0) && (ixtyp <= 4) ) */

	      /* fills e field right-hand matrix */
	      etmns( tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, ixtyp, cur);

	      /* matrix solving  (netwk calls solves) */
	      if( (nonet != 0) && (inc <= 1) )
	      {
		fprintf( output_fp, "\n\n\n"
		    "                                            "
		    "---------- NETWORK DATA ----------" );

		itmp3=0;
		itmp1= ntyp[0];

		for( i = 0; i < 2; i++ )
		{
		  if( itmp1 == 3)
		    itmp1=2;

		  if( itmp1 == 2)
		    fprintf( output_fp, "\n"
			"  -- FROM -  --- TO --      TRANSMISSION LINE       "
			" --------- SHUNT ADMITTANCES (MHOS) ---------   LINE\n"
			"  TAG   SEG  TAG   SEG    IMPEDANCE      LENGTH    "
			" ----- END ONE -----      ----- END TWO -----   TYPE\n"
			"  No:   No:  No:   No:         OHMS      METERS     "
			" REAL      IMAGINARY      REAL      IMAGINARY" );
		  else
		    if( itmp1 == 1)
		      fprintf( output_fp, "\n"
			  "  -- FROM -  --- TO --            --------"
			  " ADMITTANCE MATRIX ELEMENTS (MHOS) ---------\n"
			  "  TAG   SEG  TAG   SEG   ----- (ONE,ONE) ------  "
			  " ----- (ONE,TWO) -----   ----- (TWO,TWO) -------\n"
			  "  No:   No:  No:   No:      REAL      IMAGINARY     "
			  " REAL     IMAGINARY       REAL      IMAGINARY" );

		  for( j = 0; j < nonet; j++)
		  {
		    itmp2= ntyp[j];

		    if( (itmp2/itmp1) != 1 )
		      itmp3 = itmp2;
		    else
		    {
		      int idx4, idx5;

		      itmp4= iseg1[j];
		      itmp5= iseg2[j];
		      idx4 = itmp4-1;
		      idx5 = itmp5-1;

		      if( (itmp2 >= 2) && (x11i[j] <= 0.) )
		      {
			double xx, yy, zz;

			xx = x[idx5]- x[idx4];
			yy = y[idx5]- y[idx4];
			zz = z[idx5]- z[idx4];
			x11i[j]= wlam* sqrt( xx*xx + yy*yy + zz*zz );
		      }

		      fprintf( output_fp, "\n"
			  " %4d %5d %4d %5d  %11.4E %11.4E  %11.4E %11.4E  %11.4E %11.4E %s",
			  itag[idx4], itmp4, itag[idx5], itmp5,
			  x11r[j], x11i[j], x12r[j], x12i[j],
			  x22r[j], x22i[j], pnet[itmp2-1] );

		    } /* if(( itmp2/ itmp1) == 1) */

		  } /* for( j = 0; j < nonet; j++) */

		  if( itmp3 == 0)
		    break;

		  itmp1= itmp3;

		} /* for( j = 0; j < nonet; j++) */

	      } /* if( (nonet != 0) && (inc <= 1) ) */

	      if( (inc > 1) && (iptflg > 0) )
		nprint=1;

	      netwk( cm, &cm[ib11], &cm[ic11], &cm[id11], ip, cur );
	      ntsol=1;

	      if( iped != 0)
	      {
		itmp1= 4*( mhz-1);

		fnorm[itmp1  ]= creal( zped);
		fnorm[itmp1+1]= cimag( zped);
		fnorm[itmp1+2]= cabs( zped);
		fnorm[itmp1+3]= cang( zped);

		if( iped != 2 )
		{
		  if( fnorm[itmp1+2] > zpnorm)
		    zpnorm= fnorm[itmp1+2];
		}

	      } /* if( iped != 0) */

	      /* printing structure currents */
	      if( n != 0)
	      {
		if( iptflg != -1)
		{
		  if( iptflg <= 0)
		  {
		    fprintf( output_fp, "\n\n\n"
			"                           "
			"-------- CURRENTS AND LOCATION --------\n"
			"                                  "
			"DISTANCES IN WAVELENGTHS" );

		    fprintf( output_fp,	"\n\n"
			"   SEG  TAG    COORDINATES OF SEGM CENTER     SEGM"
			"    ------------- CURRENT (AMPS) -------------\n"
			"   No:  No:       X         Y         Z      LENGTH"
			"     REAL      IMAGINARY    MAGN        PHASE" );
		  }
		  else
		  {
		    if( (iptflg != 3) && (inc <= 1) )
		      fprintf( output_fp, "\n\n\n"
			  "             -------- RECEIVING PATTERN PARAMETERS --------\n"
			  "                      ETA: %7.2f DEGREES\n"
			  "                      TYPE: %s\n"
			  "                      AXIAL RATIO: %6.3f\n\n"
			  "            THETA     PHI      ----- CURRENT ----    SEG\n"
			  "            (DEG)    (DEG)     MAGNITUDE    PHASE    No:",
			  xpr3, hpol[ixtyp-1], xpr6 );

		  } /* if( iptflg <= 0) */

		} /* if( iptflg != -1) */

		ploss=0.;
		itmp1=0;
		jump= iptflg+1;

		for( i = 0; i < n; i++ )
		{
		  curi= cur[i]* wlam;
		  cmag= cabs( curi);
		  ph= cang( curi);

		  if( (nload != 0) && (fabs(creal(zarray[i])) >= 1.e-20) )
		    ploss += 0.5* cmag* cmag* creal( zarray[i])* si[i];

		  if( jump == 0)
		    continue;

		  if( jump > 0 )
		  {
		    if( (iptag != 0) && (itag[i] != iptag) )
		      continue;

		    itmp1++;
		    if( (itmp1 < iptagf) || (itmp1 > iptagt) )
		      continue;

		    if( iptflg != 0)
		    {
		      if( iptflg >= 2 )
		      {
			fnorm[inc-1]= cmag;
			isave= (i+1);
		      }

		      if( iptflg != 3)
		      {
			fprintf( output_fp, "\n"
			    "          %7.2f  %7.2f   %11.4E  %7.2f  %5d",
			    xpr1, xpr2, cmag, ph, i+1 );

			continue;
		      }

		    } /* if( iptflg != 0) */

		  } /* if( jump > 0 ) */
		  else
		  {
		    fprintf( output_fp, "\n"
			" %5d %4d %9.4f %9.4f %9.4f %9.5f"
			" %11.4E %11.4E %11.4E %8.3f",
			i+1, itag[i], x[i], y[i], z[i], si[i],
			creal(curi), cimag(curi), cmag, ph );

		    if( iplp1 != 1 )
		      continue;

		    if( iplp2 == 1)
		      fprintf( plot_fp, "%12.4E %12.4E\n", creal(curi), cimag(curi) );
		    else
		      if( iplp2 == 2)
			fprintf( plot_fp, "%12.4E %12.4E\n", cmag, ph );
		  }

		} /* for( i = 0; i < n; i++ ) */

		if( iptflq != -1)
		{
		  fprintf( output_fp, "\n\n\n"
		      "                                  "
		      "------ CHARGE DENSITIES ------\n"
		      "                                  "
		      "   DISTANCES IN WAVELENGTHS\n\n"
		      "   SEG   TAG    COORDINATES OF SEG CENTER     SEG        "
		      "  CHARGE DENSITY (COULOMBS/METER)\n"
		      "   NO:   NO:     X         Y         Z       LENGTH   "
		      "  REAL      IMAGINARY     MAGN        PHASE" );

		  itmp1 = 0;
		  fr = 1.e-6/fmhz;

		  for( i = 0; i < n; i++ )
		  {
		    if( iptflq != -2 )
		    {
		      if( (iptaq != 0) && (itag[i] != iptaq) )
			continue;

		      itmp1++;
		      if( (itmp1 < iptaqf) || (itmp1 > iptaqt) )
			continue;

		    } /* if( iptflq == -2) */

		    curi= fr* cmplx(- bii[i], bir[i]);
		    cmag= cabs( curi);
		    ph= cang( curi);

		    fprintf( output_fp, "\n"
			" %5d %4d %9.4f %9.4f %9.4f %9.5f"
			" %11.4E %11.4E %11.4E %9.3f",
			i+1, itag[i], x[i], y[i], z[i], si[i],
			creal(curi), cimag(curi), cmag, ph );

		  } /* for( i = 0; i < n; i++ ) */

		} /* if( iptflq != -1) */

	      } /* if( n != 0) */

	      if( m != 0)
	      {
		fprintf( output_fp, "\n\n\n"
		    "                                      "
		    " --------- SURFACE PATCH CURRENTS ---------\n"
		    "                                                "
		    " DISTANCE IN WAVELENGTHS\n"
		    "                                                "
		    " CURRENT IN AMPS/METER\n\n"
		    "                                 ---------"
		    " SURFACE COMPONENTS --------   "
		    " ---------------- RECTANGULAR COMPONENTS ----------------\n"
		    "  PCH   --- PATCH CENTER ---     TANGENT VECTOR 1    "
		    " TANGENT VECTOR 2    ------- X ------    ------- Y ------   "
		    " ------- Z ------\n  No:    X       Y       Z       MAG.    "
		    "   PHASE     MAG.       PHASE    REAL   IMAGINARY    REAL  "
		    " IMAGINARY    REAL   IMAGINARY" );

		j= n-3;
		itmp1= -1;

		for( i = 0; i < m; i++ )
		{
		  j += 3;
		  itmp1++;
		  ex= cur[j];
		  ey= cur[j+1];
		  ez= cur[j+2];
		  eth= ex* t1x[itmp1]+ ey* t1y[itmp1]+ ez* t1z[itmp1];
		  eph= ex* t2x[itmp1]+ ey* t2y[itmp1]+ ez* t2z[itmp1];
		  ethm= cabs( eth);
		  etha= cang( eth);
		  ephm= cabs( eph);
		  epha= cang( eph);

		  fprintf( output_fp, "\n"
		      " %4d %7.3f %7.3f %7.3f %11.4E %8.2f %11.4E %8.2f"
		      " %9.2E %9.2E %9.2E %9.2E %9.2E %9.2E",
		      i+1, px[itmp1], py[itmp1], pz[itmp1],
		      ethm, etha, ephm, epha, creal(ex), cimag(ex),
		      creal(ey), cimag(ey), creal(ez), cimag(ez) );

		  if( iplp1 != 1)
		    continue;

		  if( iplp3 == 1)
		    fprintf( plot_fp, "%12.4E %12.4E\n", creal(ex), cimag(ex) );
		  if( iplp3 == 2)
		    fprintf( plot_fp, "%12.4E %12.4E\n", creal(ey), cimag(ey) );
		  if( iplp3 == 3)
		    fprintf( plot_fp, "%12.4E %12.4E\n", creal(ez), cimag(ez) );
		  if( iplp3 == 4)
		    fprintf( plot_fp, "%12.4E %12.4E %12.4E %12.4E %12.4E %12.4E\n",
			creal(ex),cimag(ex),creal(ey),cimag(ey),creal(ez),cimag(ez) );

		} /* for( i=0; i<m; i++ ) */

	      } /* if( m != 0) */

	      if( (ixtyp == 0) || (ixtyp == 5) )
	      {
		tmp1= pin- pnls- ploss;
		tmp2= 100.* tmp1/ pin;

		fprintf( output_fp, "\n\n\n"
		    "                               "
		    "---------- POWER BUDGET ---------\n"
		    "                               "
		    "INPUT POWER   = %11.4E Watts\n"
		    "                               "
		    "RADIATED POWER= %11.4E Watts\n"
		    "                               "
		    "STRUCTURE LOSS= %11.4E Watts\n"
		    "                               "
		    "NETWORK LOSS  = %11.4E Watts\n"
		    "                               "
		    "EFFICIENCY    = %7.2f Percent",
		    pin, tmp1, ploss, pnls, tmp2 );

	      } /* if( (ixtyp == 0) || (ixtyp == 5) ) */

	      igo = 4;

	      if( ncoup > 0)
		couple( cur, wlam );

	      if( iflow == 7)
	      {
		if( (ixtyp > 0) && (ixtyp < 4) )
		{
		  nthic++;
		  inc++;
		  xpr1 += xpr4;

		  if( nthic <= nthi )
		    continue; /* continue excitation loop */

		  nthic=1;
		  xpr1= thetis;
		  xpr2= xpr2+ xpr5;
		  nphic++;

		  if( nphic <= nphi )
		    continue; /* continue excitation loop */

		  break;

		} /* if( (ixtyp >= 1) && (ixtyp <= 3) ) */

		if( nfrq != 1)
		{
		  jmp_floop = TRUE;
		  break; /* continue the freq loop */
		}

		fprintf( output_fp, "\n\n\n" );
		jmp_iloop = TRUE;

		break; /* continue card input loop */

	      } /*if( iflow == 7) */


	      case 4: /* label_71 */
	      igo = 5;

	      /* label_72 */
	      case 5: /* near field calculation */

	      if( near != -1)
	      {
		nfpat();

		if( mhz == nfrq)
		  near=-1;

		if( nfrq == 1)
		{
		  fprintf( output_fp, "\n\n\n" );
		  jmp_iloop = TRUE;
		  break; /* continue card input loop */
		}

	      } /* if( near != -1) */


	      /* label_78 */
	      case 6: /* standard far field calculation */

	      if( ifar != -1)
	      {
		pinr= pin;
		pnlr= pnls;
		rdpat();
	      }

	      if( (ixtyp == 0) || (ixtyp >= 4) )
	      {
		if( mhz == nfrq )
		  ifar=-1;

		if( nfrq != 1)
		{
		  jmp_floop = TRUE;
		  break;
		}

		fprintf( output_fp, "\n\n\n" );
		jmp_iloop = TRUE;
		break;

	      } /* if( (ixtyp == 0) || (ixtyp >= 4) ) */

	      nthic++;
	      inc++;
	      xpr1 += xpr4;

	      if( nthic <= nthi )
		continue; /* continue excitation loop */

	      nthic = 1;
	      xpr1  = thetis;
	      xpr2 += xpr5;
	      nphic++;

	      if( nphic > nphi )
		break;

	    } /* do (l_54) */
	    while( TRUE );

	    /* jump to freq. or input loop */
	    if( jmp_iloop )
	      break;

	    if( jmp_floop )
	      continue;

	    nphic = 1;
	    xpr2  = phiss;

	    /* normalized receiving pattern printed */
	    if( iptflg >= 2)
	    {
	      itmp1= nthi* nphi;

	      tmp1= fnorm[0];
	      for( j = 1; j < itmp1; j++ )
		if( fnorm[j] > tmp1)
		  tmp1= fnorm[j];

	      fprintf( output_fp, "\n\n\n"
		  "                     "
		  "---- NORMALIZED RECEIVING PATTERN ----\n"
		  "                      "
		  "NORMALIZATION FACTOR: %11.4E\n"
		  "                      "
		  "ETA: %7.2f DEGREES\n"
		  "                      "
		  "TYPE: %s\n"
		  "                      AXIAL RATIO: %6.3f\n"
		  "                      SEGMENT No: %d\n\n"
		  "                      "
		  "THETA     PHI       ---- PATTERN ----\n"
		  "                      "
		  "(DEG)    (DEG)       DB     MAGNITUDE",
		  tmp1, xpr3, hpol[ixtyp-1], xpr6, isave );

	      for( j = 0; j < nphi; j++ )
	      {
		itmp2= nthi*j;

		for( i = 0; i < nthi; i++ )
		{
		  itmp3= i + itmp2;

		  if( itmp3 < itmp1)
		  {
		    tmp2= fnorm[itmp3]/ tmp1;
		    tmp3= db20( tmp2);

		    fprintf( output_fp, "\n"
			"                    %7.2f  %7.2f   %7.2f  %11.4E",
			xpr1, xpr2, tmp3, tmp2 );

		    xpr1 += xpr4;
		  }

		} /* for( i = 0; i < nthi; i++ ) */

		xpr1= thetis;
		xpr2 += xpr5;

	      } /* for( j = 0; j < nphi; j++ ) */

	      xpr2= phiss;

	    } /* if( iptflg >= 2) */

	    if( mhz == nfrq)
	      ifar=-1;

	    if( nfrq == 1)
	    {
	      fprintf( output_fp, "\n\n\n" );
	      jmp_iloop = TRUE;
	      break; /* continue card input loop */
	    }

	  } /*** do (frequency loop) (l_42) ***/
	  while( (++mhz <= nfrq) );

	  /* Jump to card input loop */
	  if( jmp_iloop )
	    break;

	  if( iped != 0)
	  {
	    int iss;

	    if( nvqd > 0)
	      iss = ivqd[nvqd-1];
	    else
	      iss = isant[nsant-1];

	    fprintf( output_fp, "\n\n\n"
		"                            "
		" -------- INPUT IMPEDANCE DATA --------\n"
		"                                     "
		" SOURCE SEGMENT No: %d\n"
		"                                  "
		" NORMALIZATION FACTOR:%12.5E\n\n"
		"              ----------- UNNORMALIZED IMPEDANCE ----------  "
		"  ------------ NORMALIZED IMPEDANCE -----------\n"
		"      FREQ    RESISTANCE    REACTANCE    MAGNITUDE    PHASE  "
		"  RESISTANCE    REACTANCE    MAGNITUDE    PHASE\n"
		"       MHz       OHMS         OHMS         OHMS     DEGREES  "
		"     OHMS         OHMS         OHMS     DEGREES",
		iss, zpnorm );

	    itmp1= nfrq;
	    if( ifrq == 0)
	      tmp1= fmhz-( nfrq-1)* delfrq;
	    else
	      if( ifrq == 1)
		tmp1= fmhz/( pow(delfrq, (nfrq-1)) );

	    for( i = 0; i < itmp1; i++ )
	    {
	      itmp2= 4*i;
	      tmp2= fnorm[itmp2  ]/ zpnorm;
	      tmp3= fnorm[itmp2+1]/ zpnorm;
	      tmp4= fnorm[itmp2+2]/ zpnorm;
	      tmp5= fnorm[itmp2+3];

	      fprintf( output_fp, "\n"
		  " %9.3f   %11.4E  %11.4E  %11.4E  %7.2f  "
		  " %11.4E  %11.4E  %11.4E  %7.2f",
		  tmp1, fnorm[itmp2], fnorm[itmp2+1], fnorm[itmp2+2],
		  fnorm[itmp2+3], tmp2, tmp3, tmp4, tmp5 );

	      if( ifrq == 0)
		tmp1 += delfrq;
	      else
		if( ifrq == 1)
		  tmp1 *= delfrq;

	    } /* for( i = 0; i < itmp1; i++ ) */

	    fprintf( output_fp, "\n\n\n" );

	  } /* if( iped != 0) */

	  nfrq=1;
	  mhz=1;

      } /* switch( igox ) */

    } /* while( ! next_job ): Main input section (l_14) */

  } /* while(TRUE): Main execution loop (l_1) */

  return(0);

} /* end of main() */

/*-----------------------------------------------------------------------*/

/* arc generates segment geometry data for an arc of ns segments */
void arc( int itg, int ns, double rada,
    double ang1, double ang2, double rad )
{
  int ist, i, mreq;
  double ang, dang, xs1, xs2, zs1, zs2;

  ist= n;
  n += ns;
  np= n;
  mp= m;
  ipsym=0;

  if( ns < 1)
    return;

  if( fabs( ang2- ang1) < 360.00001)
  {
    /* Reallocate tags buffer */
    mem_realloc( (void *)&itag, (n+m) * sizeof(int) );

    /* Reallocate wire buffers */
    mreq = n * sizeof(double);
    mem_realloc( (void *)&x, mreq );
    mem_realloc( (void *)&y, mreq );
    mem_realloc( (void *)&z, mreq );
    mem_realloc( (void *)&x2, mreq );
    mem_realloc( (void *)&y2, mreq );
    mem_realloc( (void *)&z2, mreq );
    mem_realloc( (void *)&bi, mreq );

    ang= ang1* TA;
    dang=( ang2- ang1)* TA/ ns;
    xs1= rada* cos( ang);
    zs1= rada* sin( ang);

    for( i = ist; i < n; i++ )
    {
      ang += dang;
      xs2= rada* cos( ang);
      zs2= rada* sin( ang);
      x[i]= xs1;
      y[i]=0.;
      z[i]= zs1;
      x2[i]= xs2;
      y2[i]=0.;
      z2[i]= zs2;
      xs1= xs2;
      zs1= zs2;
      bi[i]= rad;
      itag[i]= itg;

    } /* for( i = ist; i < n; i++ ) */

  } /* if( fabs( ang2- ang1) < 360.00001) */
  else
  {
    fprintf( output_fp, "\n  ERROR -- ARC ANGLE EXCEEDS 360 DEGREES");
    stop(-1);
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* atgn2 is arctangent function modified to return 0 when x=y=0. */
double atgn2( double x, double y)
{
  return( atan2(y, x) );
}

/*-----------------------------------------------------------------------*/

/* cabc computes coefficients of the constant (a), sine (b), and */
/* cosine (c) terms in the current interpolation functions for the */
/* current vector cur. */
void cabc( complex double *curx)
{
  int i, is, j, jx, jco1, jco2;
  double ar, ai, sh;
  complex double curd, cs1, cs2;

  if( n != 0)
  {
    for( i = 0; i < n; i++ )
    {
      air[i]=0.;
      aii[i]=0.;
      bir[i]=0.;
      bii[i]=0.;
      cir[i]=0.;
      cii[i]=0.;
    }

    for( i = 0; i < n; i++ )
    {
      ar= creal( curx[i]);
      ai= cimag( curx[i]);
      tbf( i+1, 1 );

      for( jx = 0; jx < jsno; jx++ )
      {
	j= jco[jx]-1;
	air[j] += ax[jx]* ar;
	aii[j] += ax[jx]* ai;
	bir[j] += bx[jx]* ar;
	bii[j] += bx[jx]* ai;
	cir[j] += cx[jx]* ar;
	cii[j] += cx[jx]* ai;
      }

    } /* for( i = 0; i < n; i++ ) */

    if( nqds != 0)
    {
      for( is = 0; is < nqds; is++ )
      {
	i= iqds[is]-1;
	jx= icon1[i];
	icon1[i]=0;
	tbf(i+1,0);
	icon1[i]= jx;
	sh= si[i]*.5;
	curd= CCJ* vqds[is]/( (log(2.* sh/ bi[i])-1.)*
	    (bx[jsno-1]* cos(TP* sh)+ cx[jsno-1]* sin(TP* sh))* wlam );
	ar= creal( curd);
	ai= cimag( curd);

	for( jx = 0; jx < jsno; jx++ )
	{
	  j= jco[jx]-1;
	  air[j]= air[j]+ ax[jx]* ar;
	  aii[j]= aii[j]+ ax[jx]* ai;
	  bir[j]= bir[j]+ bx[jx]* ar;
	  bii[j]= bii[j]+ bx[jx]* ai;
	  cir[j]= cir[j]+ cx[jx]* ar;
	  cii[j]= cii[j]+ cx[jx]* ai;
	}

      } /* for( is = 0; is < nqds; is++ ) */

    } /* if( nqds != 0) */

    for( i = 0; i < n; i++ )
      curx[i]= cmplx( air[i]+cir[i], aii[i]+cii[i] );

  } /* if( n != 0) */

  if( m == 0)
    return;

  /* convert surface currents from */
  /* t1,t2 components to x,y,z components */
  jco1= np2m;
  jco2= jco1+ m;
  for( i = 1; i <= m; i++ )
  {
    jco1 -= 2;
    jco2 -= 3;
    cs1= curx[jco1];
    cs2= curx[jco1+1];
    curx[jco2]  = cs1* t1x[m-i]+ cs2* t2x[m-i];
    curx[jco2+1]= cs1* t1y[m-i]+ cs2* t2y[m-i];
    curx[jco2+2]= cs1* t1z[m-i]+ cs2* t2z[m-i];
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* cang returns the phase angle of a complex number in degrees. */
double cang( complex double z )
{
  return( carg(z)*TD );
}

/*-----------------------------------------------------------------------*/

/* cmset sets up the complex structure matrix in the array cm */
void cmset( int nrow, complex double *cm, double rkhx, int iexkx )
{
  int mp2, neq, npeq, iout, it, i, j, i1, i2, in2;
  int im1, im2, ist, ij, ipr, jss, jm1, jm2, jst, k, ka, kk;
  complex double zaj, deter, *scm = NULL;

  mp2=2* mp;
  npeq= np+ mp2;
  neq= n+2* m;

  rkh= rkhx;
  iexk= iexkx;
  iout=2* npblk* nrow;
  it= nlast;

  for( i = 0; i < nrow; i++ )
    for( j = 0; j < it; j++ )
      cm[i+j*nrow]= CPLX_00;

  i1= 1;
  i2= it;
  in2= i2;

  if( in2 > np)
    in2= np;

  im1= i1- np;
  im2= i2- np;

  if( im1 < 1)
    im1=1;

  ist=1;
  if( i1 <= np)
    ist= np- i1+2;

  /* wire source loop */
  if( n != 0)
  {
    for( j = 1; j <= n; j++ )
    {
      trio(j);
      for( i = 0; i < jsno; i++ )
      {
	ij= jco[i];
	jco[i]=(( ij-1)/ np)* mp2+ ij;
      }

      if( i1 <= in2)
	cmww( j, i1, in2, cm, nrow, cm, nrow,1);

      if( im1 <= im2)
	cmws( j, im1, im2, &cm[(ist-1)*nrow], nrow, cm, nrow, 1);

      /* matrix elements modified by loading */
      if( nload == 0)
	continue;

      if( j > np)
	continue;

      ipr= j;
      if( (ipr < 1) || (ipr > it) )
	continue;

      zaj= zarray[j-1];

      for( i = 0; i < jsno; i++ )
      {
	jss= jco[i];
	cm[(jss-1)+(ipr-1)*nrow] -= ( ax[i]+ cx[i])* zaj;
      }

    } /* for( j = 1; j <= n; j++ ) */

  } /* if( n != 0) */

  if( m != 0)
  {
    /* matrix elements for patch current sources */
    jm1=1- mp;
    jm2=0;
    jst=1- mp2;

    for( i = 0; i < nop; i++ )
    {
      jm1 += mp;
      jm2 += mp;
      jst += npeq;

      if( i1 <= in2)
	cmsw( jm1, jm2, i1, in2, &cm[(jst-1)], cm, 0, nrow, 1);

      if( im1 <= im2)
	cmss( jm1, jm2, im1, im2, &cm[(jst-1)+(ist-1)*nrow], nrow, 1);
    }

  } /* if( m != 0) */

  if( icase == 1)
    return;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, np2m * sizeof(complex double) );

  /* combine elements for symmetry modes */
  for( i = 0; i < it; i++ )
  {
    for( j = 0; j < npeq; j++ )
    {
      for( k = 0; k < nop; k++ )
      {
	ka= j+ k*npeq;
	scm[k]= cm[ka+i*nrow];
      }

      deter= scm[0];

      for( kk = 1; kk < nop; kk++ )
	deter += scm[kk];

      cm[j+i*nrow]= deter;

      for( k = 1; k < nop; k++ )
      {
	ka= j+ k*npeq;
	deter= scm[0];

	for( kk = 1; kk < nop; kk++ )
	{
	  deter += scm[kk]* ssx[k+kk*nop];
	  cm[ka+i*nrow]= deter;
	}

      } /* for( k = 1; k < nop; k++ ) */

    } /* for( j = 0; j < npeq; j++ ) */

  } /* for( i = 0; i < it; i++ ) */

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* cmss computes matrix elements for surface-surface interactions. */
void cmss( int j1, int j2, int im1, int im2,
    complex double *cm, int nrow, int itrp )
{
  int i1, i2, icomp, ii1, i, il, ii2, jj1, j, jl, jj2;
  double t1xi, t1yi, t1zi, t2xi, t2yi, t2zi, xi, yi, zi;
  complex double g11, g12, g21, g22;

  i1=( im1+1)/2;
  i2=( im2+1)/2;
  icomp= i1*2-3;
  ii1=-2;
  if( icomp+2 < im1)
    ii1=-3;

  /* loop over observation patches */
  il = -1;
  for( i = i1; i <= i2; i++ )
  {
    il++;
    icomp += 2;
    ii1 += 2;
    ii2 = ii1+1;

    t1xi= t1x[il]* psalp[il];
    t1yi= t1y[il]* psalp[il];
    t1zi= t1z[il]* psalp[il];
    t2xi= t2x[il]* psalp[il];
    t2yi= t2y[il]* psalp[il];
    t2zi= t2z[il]* psalp[il];
    xi= px[il];
    yi= py[il];
    zi= pz[il];

    /* loop over source patches */
    jj1=-2;
    for( j = j1; j <= j2; j++ )
    {
      jl=j-1;
      jj1 += 2;
      jj2 = jj1+1;

      s= pbi[jl];
      xj= px[jl];
      yj= py[jl];
      zj= pz[jl];
      t1xj= t1x[jl];
      t1yj= t1y[jl];
      t1zj= t1z[jl];
      t2xj= t2x[jl];
      t2yj= t2y[jl];
      t2zj= t2z[jl];

      hintg( xi, yi, zi);

      g11=-( t2xi* exk+ t2yi* eyk+ t2zi* ezk);
      g12=-( t2xi* exs+ t2yi* eys+ t2zi* ezs);
      g21=-( t1xi* exk+ t1yi* eyk+ t1zi* ezk);
      g22=-( t1xi* exs+ t1yi* eys+ t1zi* ezs);

      if( i == j )
      {
	g11 -= .5;
	g22 += .5;
      }

      /* normal fill */
      if( itrp == 0)
      {
	if( icomp >= im1 )
	{
	  cm[ii1+jj1*nrow]= g11;
	  cm[ii1+jj2*nrow]= g12;
	}

	if( icomp >= im2 )
	  continue;

	cm[ii2+jj1*nrow]= g21;
	cm[ii2+jj2*nrow]= g22;
	continue;

      } /* if( itrp == 0) */

      /* transposed fill */
      if( icomp >= im1 )
      {
	cm[jj1+ii1*nrow]= g11;
	cm[jj2+ii1*nrow]= g12;
      }

      if( icomp >= im2 )
	continue;

      cm[jj1+ii2*nrow]= g21;
      cm[jj2+ii2*nrow]= g22;

    } /* for( j = j1; j <= j2; j++ ) */

  } /* for( i = i1; i <= i2; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* computes matrix elements for e along wires due to patch current */
void cmsw( int j1, int j2, int i1, int i2, complex double *cm,
    complex double *cw, int ncw, int nrow, int itrp )
{
  int neqs, k, icgo, i, ipch, jl, j, js, il, ip;
  int jsnox; /* -1 offset to "jsno" for array indexing */
  double xi, yi, zi, cabi, sabi, salpi, fsign=1., pyl, pxl;
  complex double emel[9];

  neqs= np2m;
  jsnox = jsno-1;

  if( itrp >= 0)
  {
    k=-1;
    icgo=0;

    /* observation loop */
    for( i = i1-1; i < i2; i++ )
    {
      k++;
      xi= x[i];
      yi= y[i];
      zi= z[i];
      cabi= cab[i];
      sabi= sab[i];
      salpi= salp[i];
      ipch=0;

      if( icon1[i] >= PCHCON)
      {
	ipch= icon1[i]-PCHCON;
	fsign=-1.;
      }

      if( icon2[i] >= PCHCON)
      {
	ipch= icon2[i]-PCHCON;
	fsign=1.;
      }

      /* source loop */
      jl = -1;
      for( j = j1; j <= j2; j++ )
      {
	jl += 2;
	js = j-1;
	t1xj= t1x[js];
	t1yj= t1y[js];
	t1zj= t1z[js];
	t2xj= t2x[js];
	t2yj= t2y[js];
	t2zj= t2z[js];
	xj= px[js];
	yj= py[js];
	zj= pz[js];
	s= pbi[js];

	/* ground loop */
	for( ip = 1; ip <= ksymp; ip++ )
	{
	  ipgnd= ip;

	  if( ((ipch == j) || (icgo != 0)) && (ip != 2) )
	  {
	    if( icgo <= 0 )
	    {
	      pcint( xi, yi, zi, cabi, sabi, salpi, emel);

	      pyl= PI* si[i]* fsign;
	      pxl= sin( pyl);
	      pyl= cos( pyl);
	      exc= emel[8]* fsign;

	      trio(i+1);

	      il= i-ncw;
	      if( i < np)
		il += (il/np)*2*mp;

	      if( itrp == 0 )
		cw[k+il*nrow] += exc*( ax[jsnox]+ bx[jsnox]* pxl+ cx[jsnox]* pyl);
	      else
		cw[il+k*nrow] += exc*( ax[jsnox]+ bx[jsnox]* pxl+ cx[jsnox]* pyl);

	    } /* if( icgo <= 0 ) */

	    if( itrp == 0)
	    {
	      cm[k+(jl-1)*nrow]= emel[icgo];
	      cm[k+jl*nrow]    = emel[icgo+4];
	    }
	    else
	    {
	      cm[(jl-1)+k*nrow]= emel[icgo];
	      cm[jl+k*nrow]    = emel[icgo+4];
	    }

	    icgo++;
	    if( icgo == 4)
	      icgo=0;

	    continue;

	  } /* if( ((ipch == (j+1)) || (icgo != 0)) && (ip != 2) ) */

	  unere( xi, yi, zi);

	  /* normal fill */
	  if( itrp == 0)
	  {
	    cm[k+(jl-1)*nrow] += exk* cabi+ eyk* sabi+ ezk* salpi;
	    cm[k+jl*nrow]     += exs* cabi+ eys* sabi+ ezs* salpi;
	    continue;
	  }

	  /* transposed fill */
	  cm[(jl-1)+k*nrow] += exk* cabi+ eyk* sabi+ ezk* salpi;
	  cm[jl+k*nrow]     += exs* cabi+ eys* sabi+ ezs* salpi;

	} /* for( ip = 1; ip <= ksymp; ip++ ) */

      } /* for( j = j1; j <= j2; j++ ) */

    } /* for( i = i1-1; i < i2; i++ ) */

  } /* if( itrp >= 0) */

  return;
}

/*-----------------------------------------------------------------------*/

/* cmws computes matrix elements for wire-surface interactions */
void cmws( int j, int i1, int i2, complex double *cm,
    int nr, complex double *cw, int nw, int itrp )
{
  int ipr, i, ipatch, ik, js=0, ij, jx;
  double xi, yi, zi, tx, ty, tz;
  complex double etk, ets, etc;

  j--;
  s= si[j];
  b= bi[j];
  xj= x[j];
  yj= y[j];
  zj= z[j];
  cabj= cab[j];
  sabj= sab[j];
  salpj= salp[j];

  /* observation loop */
  ipr= -1;
  for( i = i1; i <= i2; i++ )
  {
    ipr++;
    ipatch=(i+1)/2;
    ik= i-( i/2)*2;

    if( (ik != 0) || (ipr == 0) )
    {
      js= ipatch-1;
      xi= px[js];
      yi= py[js];
      zi= pz[js];
      hsfld( xi, yi, zi,0.);

      if( ik != 0 )
      {
	tx= t2x[js];
	ty= t2y[js];
	tz= t2z[js];
      }
      else
      {
	tx= t1x[js];
	ty= t1y[js];
	tz= t1z[js];
      }

    } /* if( (ik != 0) || (ipr == 0) ) */
    else
    {
      tx= t1x[js];
      ty= t1y[js];
      tz= t1z[js];

    } /* if( (ik != 0) || (ipr == 0) ) */

    etk=-( exk* tx+ eyk* ty+ ezk* tz)* psalp[js];
    ets=-( exs* tx+ eys* ty+ ezs* tz)* psalp[js];
    etc=-( exc* tx+ eyc* ty+ ezc* tz)* psalp[js];

    /* fill matrix elements.  element locations */
    /* determined by connection data. */

    /* normal fill */
    if( itrp == 0)
    {
      for( ij = 0; ij < jsno; ij++ )
      {
	jx= jco[ij]-1;
	cm[ipr+jx*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      }

      continue;
    } /* if( itrp == 0) */

    /* transposed fill */
    if( itrp != 2)
    {
      for( ij = 0; ij < jsno; ij++ )
      {
	jx= jco[ij]-1;
	cm[jx+ipr*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      }

      continue;
    } /* if( itrp != 2) */

    /* transposed fill - c(ws) and d(ws)prime (=cw) */
    for( ij = 0; ij < jsno; ij++ )
    {
      jx= jco[ij]-1;
      if( jx < nr)
	cm[jx+ipr*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      else
      {
	jx -= nr;
	cw[jx+ipr*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      }
    } /* for( ij = 0; ij < jsno; ij++ ) */

  } /* for( i = i1; i <= i2; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* cmww computes matrix elements for wire-wire interactions */
void cmww( int j, int i1, int i2, complex double *cm,
    int nr, complex double *cw, int nw, int itrp)
{
  int ipr, iprx, i, ij, jx;
  double xi, yi, zi, ai, cabi, sabi, salpi;
  complex double etk, ets, etc;

  /* set source segment parameters */
  jx = j;
  j--;
  s= si[j];
  b= bi[j];
  xj= x[j];
  yj= y[j];
  zj= z[j];
  cabj= cab[j];
  sabj= sab[j];
  salpj= salp[j];

  /* decide whether ext. t.w. approx. can be used */
  if( iexk != 0)
  {
    ipr = icon1[j];
    if( ipr < 0 )
    {
      ipr= -ipr;
      iprx= ipr-1;

      if( -icon1[iprx] != jx )
	ind1=2;
      else
      {
	xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.e-6) )
	  ind1=2;
	else
	  ind1=0;

      } /* if( -icon1[iprx] != jx ) */

    } /* if( ipr < 0 ) */
    else
    {
      iprx = ipr-1;
      if( ipr == 0 )
	ind1=1;
      else
      {
	if( ipr != jx )
	{
	  if( icon2[iprx] != jx )
	    ind1=2;
	  else
	  {
	    xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	    if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.e-6) )
	      ind1=2;
	    else
	      ind1=0;

	  } /* if( icon2[iprx] != jx ) */

	} /* if( ipr != jx ) */
	else
	  if( cabj* cabj+ sabj* sabj > 1.e-8)
	    ind1=2;
	  else
	    ind1=0;

      } /* if( ipr == 0 ) */

    } /* if( ipr < 0 ) */

    ipr = icon2[j];
    if( ipr < 0 )
    {
      ipr= -ipr;
      iprx = ipr-1;
      if( -icon2[iprx] != jx )
	ind2=2;
      else
      {
	xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.e-6) )
	  ind2=2;
	else
	  ind2=0;

      } /* if( -icon1[iprx] != jx ) */

    } /* if( ipr < 0 ) */
    else
    {
      iprx = ipr-1;
      if( ipr == 0 )
	ind2=1;
      else
      {
	if( ipr != jx )
	{
	  if( icon1[iprx] != jx )
	    ind2=2;
	  else
	  {
	    xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	    if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.e-6) )
	      ind2=2;
	    else
	      ind2=0;

	  } /* if( icon2[iprx] != jx ) */

	} /* if( ipr != jx ) */
	else
	  if( cabj* cabj+ sabj* sabj > 1.e-8)
	    ind2=2;
	  else
	    ind2=0;

      } /* if( ipr == 0 ) */

    } /* if( ipr < 0 ) */

  } /* if( iexk != 0) */

  /* observation loop */
  ipr=-1;
  for( i = i1-1; i < i2; i++ )
  {
    ipr++;
    ij= i-j;
    xi= x[i];
    yi= y[i];
    zi= z[i];
    ai= bi[i];
    cabi= cab[i];
    sabi= sab[i];
    salpi= salp[i];

    efld( xi, yi, zi, ai, ij);

    etk= exk* cabi+ eyk* sabi+ ezk* salpi;
    ets= exs* cabi+ eys* sabi+ ezs* salpi;
    etc= exc* cabi+ eyc* sabi+ ezc* salpi;

    /* fill matrix elements. element locations */
    /* determined by connection data. */

    /* normal fill */
    if( itrp == 0)
    {
      for( ij = 0; ij < jsno; ij++ )
      {
	jx = jco[ij]-1;
	cm[ipr+jx*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      }
      continue;
    }

    /* transposed fill */
    if( itrp != 2)
    {
      for( ij = 0; ij < jsno; ij++ )
      {
	jx= jco[ij]-1;
	cm[jx+ipr*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      }
      continue;
    }

    /* trans. fill for c(ww) - test for elements for d(ww)prime.  (=cw) */
    for( ij = 0; ij < jsno; ij++ )
    {
      jx= jco[ij]-1;
      if( jx < nr)
	cm[jx+ipr*nr] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      else
      {
	jx -= nr;
	cw[jx*ipr*nw] += etk* ax[ij]+ ets* bx[ij]+ etc* cx[ij];
      }

    } /* for( ij = 0; ij < jsno; ij++ ) */

  } /* for( i = i1-1; i < i2; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* connect sets up segment connection data in arrays icon1 and */
/* icon2 by searching for segment ends that are in contact. */
void conect( int ignd )
{
  int i, iz, ic, j, jx, ix, ixx, iseg, iend, jend, nsflg, jump, ipf;
  double sep=0., xi1, yi1, zi1, xi2, yi2, zi2;
  double slen, xa, ya, za, xs, ys, zs;

  nscon= -1;
  maxcon = 1;

  if( ignd != 0)
  {
    fprintf( output_fp, "\n\n     GROUND PLANE SPECIFIED." );

    if( ignd > 0)
      fprintf( output_fp,
	  "\n     WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL"
	  " BE INTERPOLATED TO IMAGE IN GROUND PLANE.\n" );

    if( ipsym == 2)
    {
      np=2* np;
      mp=2* mp;
    }

    if( abs( ipsym) > 2 )
    {
      np= n;
      mp= m;
    }

    /*** possibly should be error condition?? **/
    if( np > n)
    {
      fprintf( output_fp,
	  "\n ERROR: NP > N IN CONECT()" );
      stop(-1);
    }

    if( (np == n) && (mp == m) )
      ipsym=0;

  } /* if( ignd != 0) */

  if( n != 0)
  {
    /* Allocate memory to connections */
    mem_alloc( (void *)&icon1, (n+m) * sizeof(int) );
    mem_alloc( (void *)&icon2, (n+m) * sizeof(int) );

    for( i = 0; i < n; i++ )
    {
      iz = i+1;
      xi1= x[i];
      yi1= y[i];
      zi1= z[i];
      xi2= x2[i];
      yi2= y2[i];
      zi2= z2[i];
      slen= sqrt( (xi2- xi1)*(xi2- xi1) + (yi2- yi1) *
	  (yi2- yi1) + (zi2- zi1)*(zi2- zi1) ) * SMIN;

      /* determine connection data for end 1 of segment. */
      jump = FALSE;
      if( ignd > 0)
      {
	if( zi1 <= -slen)
	{
	  fprintf( output_fp,
	      "\n  GEOMETRY DATA ERROR -- SEGMENT"
	      " %d EXTENDS BELOW GROUND", iz );
	  stop(-1);
	}

	if( zi1 <= slen)
	{
	  icon1[i]= iz;
	  z[i]=0.;
	  jump = TRUE;

	} /* if( zi1 <= slen) */

      } /* if( ignd > 0) */

      if( ! jump )
      {
	ic= i;
	for( j = 1; j < n; j++)
	{
	  ic++;
	  if( ic >= n)
	    ic=0;

	  sep= fabs( xi1- x[ic])+ fabs(yi1- y[ic])+ fabs(zi1- z[ic]);
	  if( sep <= slen)
	  {
	    icon1[i]= -(ic+1);
	    break;
	  }

	  sep= fabs( xi1- x2[ic])+ fabs(yi1- y2[ic])+ fabs(zi1- z2[ic]);
	  if( sep <= slen)
	  {
	    icon1[i]= (ic+1);
	    break;
	  }

	} /* for( j = 1; j < n; j++) */

	if( ((iz > 0) || (icon1[i] <= PCHCON)) && (sep > slen) )
	  icon1[i]=0;

      } /* if( ! jump ) */

      /* determine connection data for end 2 of segment. */
      if( (ignd > 0) || jump )
      {
	if( zi2 <= -slen)
	{
	  fprintf( output_fp,
	      "\n  GEOMETRY DATA ERROR -- SEGMENT"
	      " %d EXTENDS BELOW GROUND", iz );
	  stop(-1);
	}

	if( zi2 <= slen)
	{
	  if( icon1[i] == iz )
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR -- SEGMENT"
		" %d LIES IN GROUND PLANE", iz );
	    stop(-1);
	  }

	  icon2[i]= iz;
	  z2[i]=0.;
	  continue;

	} /* if( zi2 <= slen) */

      } /* if( ignd > 0) */

      ic= i;
      for( j = 1; j < n; j++ )
      {
	ic++;
	if( ic >= n)
	  ic=0;

	sep= fabs(xi2- x[ic])+ fabs(yi2- y[ic])+ fabs(zi2- z[ic]);
	if( sep <= slen)
	{
	  icon2[i]= (ic+1);
	  break;
	}

	sep= fabs(xi2- x2[ic])+ fabs(yi2- y2[ic])+ fabs(zi2- z2[ic]);
	if( sep <= slen)
	{
	  icon2[i]= -(ic+1);
	  break;
	}

      } /* for( j = 1; j < n; j++ ) */

      if( ((iz > 0) || (icon2[i] <= PCHCON)) && (sep > slen) )
	icon2[i]=0;

    } /* for( i = 0; i < n; i++ ) */

    /* find wire-surface connections for new patches */
    if( m != 0)
    {
      ix = -1;
      i = 0;
      while( ++i <= m )
      {
	ix++;
	xs= px[ix];
	ys= py[ix];
	zs= pz[ix];

	for( iseg = 0; iseg < n; iseg++ )
	{
	  xi1= x[iseg];
	  yi1= y[iseg];
	  zi1= z[iseg];
	  xi2= x2[iseg];
	  yi2= y2[iseg];
	  zi2= z2[iseg];

	  /* for first end of segment */
	  slen=( fabs(xi2- xi1)+ fabs(yi2- yi1)+ fabs(zi2- zi1))* SMIN;
	  sep= fabs(xi1- xs)+ fabs(yi1- ys)+ fabs(zi1- zs);

	  /* connection - divide patch into 4 patches at present array loc. */
	  if( sep <= slen)
	  {
	    icon1[iseg]=PCHCON+ i;
	    ic=0;
	    subph( i, ic );
	    break;
	  }

	  sep= fabs(xi2- xs)+ fabs(yi2- ys)+ fabs(zi2- zs);
	  if( sep <= slen)
	  {
	    icon2[iseg]=PCHCON+ i;
	    ic=0;
	    subph( i, ic );
	    break;
	  }

	} /* for( iseg = 0; iseg < n; iseg++ ) */

      } /* while( ++i <= m ) */

    } /* if( m != 0) */

  } /* if( n != 0) */

  fprintf( output_fp, "\n\n"
      "     TOTAL SEGMENTS USED: %d   SEGMENTS IN A"
      " SYMMETRIC CELL: %d   SYMMETRY FLAG: %d",
      n, np, ipsym );

  if( m > 0)
    fprintf( output_fp,	"\n"
	"       TOTAL PATCHES USED: %d   PATCHES"
	" IN A SYMMETRIC CELL: %d",  m, mp );

  iseg=( n+ m)/( np+ mp);
  if( iseg != 1)
  {
    /*** may be error condition?? ***/
    if( ipsym == 0 )
    {
      fprintf( output_fp,
	  "\n  ERROR: IPSYM=0 IN CONECT()" );
      stop(-1);
    }

    if( ipsym < 0 )
      fprintf( output_fp,
	  "\n  STRUCTURE HAS %d FOLD ROTATIONAL SYMMETRY\n", iseg );
    else
    {
      ic= iseg/2;
      if( iseg == 8)
	ic=3;
      fprintf( output_fp,
	  "\n  STRUCTURE HAS %d PLANES OF SYMMETRY\n", ic );
    } /* if( ipsym < 0 ) */

  } /* if( iseg == 1) */

  if( n == 0)
    return;

  /* Allocate to connection buffers */
  mem_alloc( (void *)&jco, maxcon * sizeof(int) );

  /* adjust connected seg. ends to exactly coincide.  print junctions */
  /* of 3 or more seg.  also find old seg. connecting to new seg. */
  iseg = 0;
  ipf = FALSE;
  for( j = 0; j < n; j++ )
  {
    jx = j+1;
    iend=-1;
    jend=-1;
    ix= icon1[j];
    ic=1;
    jco[0]= -jx;
    xa= x[j];
    ya= y[j];
    za= z[j];

    while( TRUE )
    {
      if( (ix != 0) && (ix != (j+1)) && (ix <= PCHCON) )
      {
	nsflg=0;

	do
	{
	  if( ix == 0 )
	  {
	    fprintf( output_fp,
		"\n  CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT: %d", ix );
	    stop(-1);
	  }

	  if( ix < 0 )
	    ix= -ix;
	  else
	    jend= -jend;

	  jump = FALSE;

	  if( ix == jx )
	    break;

	  if( ix < jx )
	  {
	    jump = TRUE;
	    break;
	  }

	  /* Record max. no. of connections */
	  ic++;
	  if( ic >= maxcon )
	  {
	    maxcon = ic+1;
	    mem_realloc( (void *)&jco, maxcon * sizeof(int) );
	  }
	  jco[ic-1]= ix* jend;

	  if( ix > 0)
	    nsflg=1;

	  ixx = ix-1;
	  if( jend != 1)
	  {
	    xa= xa+ x[ixx];
	    ya= ya+ y[ixx];
	    za= za+ z[ixx];
	    ix= icon1[ixx];
	    continue;
	  }

	  xa= xa+ x2[ixx];
	  ya= ya+ y2[ixx];
	  za= za+ z2[ixx];
	  ix= icon2[ixx];

	} /* do */
	while( ix != 0 );

	if( jump && (iend == 1) )
	  break;
	else
	  if( jump )
	  {
	    iend=1;
	    jend=1;
	    ix= icon2[j];
	    ic=1;
	    jco[0]= jx;
	    xa= x2[j];
	    ya= y2[j];
	    za= z2[j];
	    continue;
	  }

	sep= (double)ic;
	xa= xa/ sep;
	ya= ya/ sep;
	za= za/ sep;

	for( i = 0; i < ic; i++ )
	{
	  ix= jco[i];
	  if( ix <= 0)
	  {
	    ix=- ix;
	    ixx = ix-1;
	    x[ixx]= xa;
	    y[ixx]= ya;
	    z[ixx]= za;
	    continue;
	  }

	  ixx = ix-1;
	  x2[ixx]= xa;
	  y2[ixx]= ya;
	  z2[ixx]= za;

	} /* for( i = 0; i < ic; i++ ) */

	if( ic >= 3)
	{
	  if( ! ipf )
	  {
	    fprintf( output_fp, "\n\n"
		"    ---------- MULTIPLE WIRE JUNCTIONS ----------\n"
		"    JUNCTION  SEGMENTS (- FOR END 1, + FOR END 2)" );
	    ipf = TRUE;
	  }

	  iseg++;
	  fprintf( output_fp, "\n   %5d      ", iseg );

	  for( i = 1; i <= ic; i++ )
	  {
	    fprintf( output_fp, "%5d", jco[i-1] );
	    if( !(i % 20) )
	      fprintf( output_fp, "\n              " );
	  }

	} /* if( ic >= 3) */

      } /*if( (ix != 0) && (ix != j) && (ix <= PCHCON) ) */

      if( iend == 1)
	break;

      iend=1;
      jend=1;
      ix= icon2[j];
      ic=1;
      jco[0]= jx;
      xa= x2[j];
      ya= y2[j];
      za= z2[j];

    } /* while( TRUE ) */

  } /* for( j = 0; j < n; j++ ) */

  mem_alloc( (void *)&ax, maxcon * sizeof(double) );
  mem_alloc( (void *)&bx, maxcon * sizeof(double) );
  mem_alloc( (void *)&cx, maxcon * sizeof(double) );

  return;
}

/*-----------------------------------------------------------------------*/

/* couple computes the maximum coupling between pairs of segments. */
void couple( complex double *cur, double wlam )
{
  int j, j1, j2, l1, i, k, itt1, itt2, its1, its2, isg1, isg2, npm1;
  double dbc, c, gmax;
  complex double y11, y12, y22, yl, yin, zl, zin, rho;

  if( (nsant != 1) || (nvqd != 0) )
    return;

  j= isegno( nctag[icoup], ncseg[icoup]);
  if( j != isant[0] )
    return;

  zin= vsant[0];
  icoup++;
  mem_realloc( (void *)&y11a, icoup * sizeof( complex double) );
  y11a[icoup-1]= cur[j-1]*wlam/zin;

  l1=(icoup-1)*(ncoup-1);
  for( i = 0; i < ncoup; i++ )
  {
    if( (i+1) == icoup)
      continue;

    l1++;
    mem_realloc( (void *)&y12a, l1 * sizeof( complex double) );
    k= isegno( nctag[i], ncseg[i]);
    y12a[l1-1]= cur[k-1]* wlam/ zin;
  }

  if( icoup < ncoup)
    return;

  fprintf( output_fp, "\n\n\n"
      "                        -----------"
      " ISOLATION DATA -----------\n\n"
      " ------- COUPLING BETWEEN ------     MAXIMUM    "
      " ---------- FOR MAXIMUM COUPLING ----------\n"
      "            SEG              SEG    COUPLING  LOAD"
      " IMPEDANCE (2ND SEG)         INPUT IMPEDANCE \n"
      " TAG  SEG   No:   TAG  SEG   No:      (DB)       "
      " REAL     IMAGINARY         REAL       IMAGINARY" );

  npm1= ncoup-1;

  for( i = 0; i < npm1; i++ )
  {
    itt1= nctag[i];
    its1= ncseg[i];
    isg1= isegno( itt1, its1);
    l1= i+1;

    for( j = l1; j < ncoup; j++ )
    {
      itt2= nctag[j];
      its2= ncseg[j];
      isg2= isegno( itt2, its2);
      j1= j+ i* npm1-1;
      j2= i+ j* npm1;
      y11= y11a[i];
      y22= y11a[j];
      y12=.5*( y12a[j1]+ y12a[j2]);
      yin= y12* y12;
      dbc= cabs( yin);
      c= dbc/(2.* creal( y11)* creal( y22)- creal( yin));

      if( (c >= 0.0) && (c <= 1.0) )
      {
	if( c >= .01 )
	  gmax=(1.- sqrt(1.- c*c))/c;
	else
	  gmax=.5*( c+.25* c* c* c);

	rho= gmax* conj( yin)/ dbc;
	yl=((1.- rho)/(1.+ rho)+1.)* creal( y22)- y22;
	zl=1./ yl;
	yin= y11- yin/( y22+ yl);
	zin=1./ yin;
	dbc= db10( gmax);

	fprintf( output_fp, "\n"
	    " %4d %4d %5d  %4d %4d %5d  %9.3f"
	    "  %12.5E %12.5E  %12.5E %12.5E",
	    itt1, its1, isg1, itt2, its2, isg2, dbc,
	    creal(zl), cimag(zl), creal(zin), cimag(zin) );

	continue;

      } /* if( (c >= 0.0) && (c <= 1.0) ) */

      fprintf( output_fp, "\n"
	  " %4d %4d %5d   %4d %4d %5d  **ERROR** "
	  "COUPLING IS NOT BETWEEN 0 AND 1. (= %12.5E)",
	  itt1, its1, isg1, itt2, its2, isg2, c );

    } /* for( j = l1; j < ncoup; j++ ) */

  } /* for( i = 0; i < npm1; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* datagn is the main routine for input of geometry data. */
void datagn( void )
{
  char gm[3];
  char ifx[2] = {'*', 'X'}, ify[2]={'*','Y'}, ifz[2]={'*','Z'};
  char ipt[4] = { 'P', 'R', 'T', 'Q' };

  /* input card mnemonic list */
  /* "XT" stands for "exit", added for testing */
#define GM_NUM  12
  char *atst[GM_NUM] =
  {
    "GW", "GX", "GR", "GS", "GE", "GM", \
    "SP", "SM", "GA", "SC", "GH", "GF"
  };

  int nwire, isct, iphd, i1, i2, itg, iy, iz, mreq;
  int ix, i, ns, gm_num; /* geometry card id as a number */
  double rad, xs1, xs2, ys1, ys2, zs1, zs2, x4=0, y4=0, z4=0;
  double x3=0, y3=0, z3=0, xw1, xw2, yw1, yw2, zw1, zw2;
  double dummy;

  ipsym=0;
  nwire=0;
  n=0;
  np=0;
  m=0;
  mp=0;
  isct=0;
  iphd = FALSE;

  /* read geometry data card and branch to */
  /* section for operation requested */
  do
  {
    readgm( gm, &itg, &ns, &xw1, &yw1, &zw1, &xw2, &yw2, &zw2, &rad);

    /* identify card id mnemonic */
    for( gm_num = 0; gm_num < GM_NUM; gm_num++ )
      if( strncmp( gm, atst[gm_num], 2) == 0 )
	break;

    if( iphd == FALSE )
    {
      fprintf( output_fp, "\n\n\n"
	  "                               "
	  "-------- STRUCTURE SPECIFICATION --------\n"
	  "                                     "
	  "COORDINATES MUST BE INPUT IN\n"
	  "                                     "
	  "METERS OR BE SCALED TO METERS\n"
	  "                                     "
	  "BEFORE STRUCTURE INPUT IS ENDED\n" );

      fprintf( output_fp, "\n"
	  "  WIRE                                           "
	  "                                      SEG FIRST  LAST  TAG\n"
	  "   No:        X1         Y1         Z1         X2      "
	  "   Y2         Z2       RADIUS   No:   SEG   SEG  No:" );

      iphd=1;
    }

    if( gm_num != 10 )
      isct=0;

    switch( gm_num )
    {
      case 0: /* "gw" card, generate segment data for straight wire. */

	nwire++;
	i1= n+1;
	i2= n+ ns;

	fprintf( output_fp, "\n"
	    " %5d  %10.4f %10.4f %10.4f %10.4f"
	    " %10.4f %10.4f %10.4f %5d %5d %5d %4d",
	    nwire, xw1, yw1, zw1, xw2, yw2, zw2, rad, ns, i1, i2, itg );

	if( rad != 0)
	{
	  xs1=1.;
	  ys1=1.;
	}
	else
	{
	  readgm( gm, &ix, &iy, &xs1, &ys1, &zs1,
	      &dummy, &dummy, &dummy, &dummy);

	  if( strcmp(gm, "GC" ) != 0 )
	  {
	    fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
	    stop(-1);
	  }

	  fprintf( output_fp,
	      "\n  ABOVE WIRE IS TAPERED.  SEGMENT LENGTH RATIO: %9.5f\n"
	      "                                 "
	      "RADIUS FROM: %9.5f TO: %9.5f", xs1, ys1, zs1 );

	  if( (ys1 == 0) || (zs1 == 0) )
	  {
	    fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
	    stop(-1);
	  }

	  rad= ys1;
	  ys1= pow( (zs1/ys1), (1./(ns-1.)) );
	}

	wire( xw1, yw1, zw1, xw2, yw2, zw2, rad, xs1, ys1, ns, itg);

	continue;

	/* reflect structure along x,y, or z */
	/* axes or rotate to form cylinder.  */
      case 1: /* "gx" card */

	iy= ns/10;
	iz= ns- iy*10;
	ix= iy/10;
	iy= iy- ix*10;

	if( ix != 0)
	  ix=1;
	if( iy != 0)
	  iy=1;
	if( iz != 0)
	  iz=1;

	fprintf( output_fp,
	    "\n  STRUCTURE REFLECTED ALONG THE AXES %c %c %c"
	    " - TAGS INCREMENTED BY %d\n",
	    ifx[ix], ify[iy], ifz[iz], itg );

	reflc( ix, iy, iz, itg, ns);

	continue;

      case 2: /* "gr" card */

	fprintf( output_fp,
	    "\n  STRUCTURE ROTATED ABOUT Z-AXIS %d TIMES"
	    " - LABELS INCREMENTED BY %d\n", ns, itg );

	ix=-1;
	iz = 0;
	reflc( ix, iy, iz, itg, ns);

	continue;

      case 3: /* "gs" card, scale structure dimensions by factor xw1. */

	if( n > 0)
	{
	  for( i = 0; i < n; i++ )
	  {
	    x[i]= x[i]* xw1;
	    y[i]= y[i]* xw1;
	    z[i]= z[i]* xw1;
	    x2[i]= x2[i]* xw1;
	    y2[i]= y2[i]* xw1;
	    z2[i]= z2[i]* xw1;
	    bi[i]= bi[i]* xw1;
	  }
	} /* if( n >= n2) */

	if( m > 0)
	{
	  yw1= xw1* xw1;
	  for( i = 0; i < m; i++ )
	  {
	    px[i]= px[i]* xw1;
	    py[i]= py[i]* xw1;
	    pz[i]= pz[i]* xw1;
	    pbi[i]= pbi[i]* yw1;
	  }
	} /* if( m >= m2) */

	fprintf( output_fp,
	    "\n     STRUCTURE SCALED BY FACTOR: %10.5f", xw1 );

	continue;

      case 4: /* "ge" card, terminate structure geometry input. */

	if( ns != 0)
	{
	  iplp1=1;
	  iplp2=1;
	}

	conect( itg);

	if( n != 0)
	{
	  /* Allocate wire buffers */
	  mreq = n * sizeof(double);
	  mem_alloc( (void *)&si, mreq );
	  mem_alloc( (void *)&sab, mreq );
	  mem_alloc( (void *)&cab, mreq );
	  mem_alloc( (void *)&salp, mreq );

	  fprintf( output_fp, "\n\n\n"
	      "                              "
	      " ---------- SEGMENTATION DATA ----------\n"
	      "                                       "
	      " COORDINATES IN METERS\n"
	      "                           "
	      " I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I\n" );

	  fprintf( output_fp, "\n"
	      "   SEG    COORDINATES OF SEGM CENTER     SEGM    ORIENTATION"
	      " ANGLES    WIRE    CONNECTION DATA   TAG\n"
	      "   No:       X         Y         Z      LENGTH     ALPHA     "
	      " BETA    RADIUS    I-     I    I+   NO:" );

	  for( i = 0; i < n; i++ )
	  {
	    xw1= x2[i]- x[i];
	    yw1= y2[i]- y[i];
	    zw1= z2[i]- z[i];
	    x[i]=( x[i]+ x2[i])*.5;
	    y[i]=( y[i]+ y2[i])*.5;
	    z[i]=( z[i]+ z2[i])*.5;
	    xw2= xw1* xw1+ yw1* yw1+ zw1* zw1;
	    yw2= sqrt( xw2);
	    yw2=( xw2/ yw2+ yw2)*.5;
	    si[i]= yw2;
	    cab[i]= xw1/ yw2;
	    sab[i]= yw1/ yw2;
	    xw2= zw1/ yw2;

	    if( xw2 > 1.)
	      xw2=1.;
	    if( xw2 < -1.)
	      xw2=-1.;

	    salp[i]= xw2;
	    xw2= asin( xw2)* TD;
	    yw2= atan2( yw1, xw1)* TD;

	    fprintf( output_fp, "\n"
		" %5d %9.4f %9.4f %9.4f %9.4f"
		" %9.4f %9.4f %9.4f %5d %5d %5d %5d",
		i+1, x[i], y[i], z[i], si[i], xw2, yw2,
		bi[i], icon1[i], i+1, icon2[i], itag[i] );

	    if( iplp1 == 1)
	      fprintf( plot_fp, "%12.4E %12.4E %12.4E "
		  "%12.4E %12.4E %12.4E %12.4E %5d %5d %5d\n",
		  x[i],y[i],z[i],si[i],xw2,yw2,bi[i],icon1[i],i+1,icon2[i] );

	    if( (si[i] <= 1.e-20) || (bi[i] <= 0.) )
	    {
	      fprintf( output_fp, "\n SEGMENT DATA ERROR" );
	      stop(-1);
	    }

	  } /* for( i = 0; i < n; i++ ) */

	} /* if( n != 0) */

	if( m != 0)
	{
	  fprintf( output_fp, "\n\n\n"
	      "                                   "
	      " --------- SURFACE PATCH DATA ---------\n"
	      "                                            "
	      " COORDINATES IN METERS\n\n"
	      " PATCH      COORD. OF PATCH CENTER           UNIT NORMAL VECTOR      "
	      " PATCH           COMPONENTS OF UNIT TANGENT VECTORS\n"
	      "  NO:       X          Y          Z          X        Y        Z      "
	      " AREA         X1       Y1       Z1        X2       Y2      Z2" );

	  for( i = 0; i < m; i++ )
	  {
	    xw1=( t1y[i]* t2z[i]- t1z[i]* t2y[i])* psalp[i];
	    yw1=( t1z[i]* t2x[i]- t1x[i]* t2z[i])* psalp[i];
	    zw1=( t1x[i]* t2y[i]- t1y[i]* t2x[i])* psalp[i];

	    fprintf( output_fp, "\n"
		" %4d %10.5f %10.5f %10.5f  %8.4f %8.4f %8.4f"
		" %10.5f  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f",
		i+1, px[i], py[i], pz[i], xw1, yw1, zw1, pbi[i],
		t1x[i], t1y[i], t1z[i], t2x[i], t2y[i], t2z[i] );

	  } /* for( i = 0; i < m; i++ ) */

	} /* if( m == 0) */

	npm  = n+m;
	np2m = n+2*m;
	np3m = n+3*m;

	return;

	/* "gm" card, move structure or reproduce */
	/* original structure in new positions.   */
      case 5:

	fprintf( output_fp,
	    "\n     THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS:\n"
	    "   %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
	    itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad );

	xw1= xw1* TA;
	yw1= yw1* TA;
	zw1= zw1* TA;

	move( xw1, yw1, zw1, xw2, yw2, zw2, (int)( rad+.5), ns, itg);
	continue;

      case 6: /* "sp" card, generate single new patch */

	i1= m+1;
	ns++;

	if( itg != 0)
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	fprintf( output_fp, "\n"
	    " %5d%c %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
	    i1, ipt[ns-1], xw1, yw1, zw1, xw2, yw2, zw2 );

	if( (ns == 2) || (ns == 4) )
	  isct=1;

	if( ns > 1)
	{
	  readgm( gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);

	  if( (ns == 2) || (itg > 0) )
	  {
	    x4= xw1+ x3- xw2;
	    y4= yw1+ y3- yw2;
	    z4= zw1+ z3- zw2;
	  }

	  fprintf( output_fp, "\n"
	      "      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
	      x3, y3, z3, x4, y4, z4 );

	  if( strcmp(gm, "SC") != 0 )
	  {
	    fprintf( output_fp, "\n  PATCH DATA ERROR" );
	    stop(-1);
	  }

	} /* if( ns > 1) */
	else
	{
	  xw2= xw2* TA;
	  yw2= yw2* TA;
	}

	patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

	continue;

      case 7: /* "sm" card, generate multiple-patch surface */

	i1= m+1;
	fprintf( output_fp, "\n"
	    " %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f"
	    "     SURFACE - %d BY %d PATCHES",
	    i1, ipt[1], xw1, yw1, zw1, xw2, yw2, zw2, itg, ns );

	if( (itg < 1) || (ns < 1) )
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	readgm( gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);

	if( (ns == 2) || (itg > 0) )
	{
	  x4= xw1+ x3- xw2;
	  y4= yw1+ y3- yw2;
	  z4= zw1+ z3- zw2;
	}

	fprintf( output_fp, "\n"
	    "      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
	    x3, y3, z3, x4, y4, z4 );

	if( strcmp(gm, "SC" ) != 0 )
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

	continue;

      case 8: /* "ga" card, generate segment data for wire arc */

	nwire++;
	i1= n+1;
	i2= n+ ns;

	fprintf( output_fp, "\n"
	    " %5d  ARC RADIUS: %9.5f  FROM: %8.3f TO: %8.3f DEGREES"
	    "       %11.5f %5d %5d %5d %4d",
	    nwire, xw1, yw1, zw1, xw2, ns, i1, i2, itg );

	arc( itg, ns, xw1, yw1, zw1, xw2);

	continue;

      case 9: /* "sc" card */

	if( isct == 0)
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	i1= m+1;
	ns++;

	if( (itg != 0) || ((ns != 2) && (ns != 4)) )
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	xs1= x4;
	ys1= y4;
	zs1= z4;
	xs2= x3;
	ys2= y3;
	zs2= z3;
	x3= xw1;
	y3= yw1;
	z3= zw1;

	if( ns == 4)
	{
	  x4= xw2;
	  y4= yw2;
	  z4= zw2;
	}

	xw1= xs1;
	yw1= ys1;
	zw1= zs1;
	xw2= xs2;
	yw2= ys2;
	zw2= zs2;

	if( ns != 4)
	{
	  x4= xw1+ x3- xw2;
	  y4= yw1+ y3- yw2;
	  z4= zw1+ z3- zw2;
	}

	fprintf( output_fp, "\n"
	    " %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
	    i1, ipt[ns-1], xw1, yw1, zw1, xw2, yw2, zw2 );

	fprintf( output_fp, "\n"
	    "      %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f",
	    x3, y3, z3, x4, y4, z4 );

	patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

	continue;

      case 10: /* "gh" card, generate helix */

	nwire++;
	i1= n+1;
	i2= n+ ns;

	fprintf( output_fp, "\n"
	    " %5d HELIX STRUCTURE - SPACING OF TURNS: %8.3f AXIAL"
	    " LENGTH: %8.3f  %8.3f %5d %5d %5d %4d\n      "
	    " RADIUS X1:%8.3f Y1:%8.3f X2:%8.3f Y2:%8.3f ",
	    nwire, xw1, yw1, rad, ns, i1, i2, itg, zw1, xw2, yw2, zw2 );

	helix( xw1, yw1, zw1, xw2, yw2, zw2, rad, ns, itg);

	continue;

      case 11: /* "gf" card, not supported */
	abort_on_error(-5);

      default: /* error message */

	fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
	fprintf( output_fp, "\n"
	    " %2s %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
	    gm, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad );

	stop(-1);

    } /* switch( gm_num ) */

  } /* do */
  while( TRUE );

  return;
}

/*-----------------------------------------------------------------------*/

/* function db10 returns db for magnitude (field) */
double db10( double x )
{
  if( x < 1.e-20 )
    return( -999.99 );

  return( 10. * log10(x) );
}

/*-----------------------------------------------------------------------*/

/* function db20 returns db for mag**2 (power) i */
double db20( double x )
{
  if( x < 1.e-20 )
    return( -999.99 );

  return( 20. * log10(x) );
}

/*-----------------------------------------------------------------------*/

/* compute near e fields of a segment with sine, cosine, and */
/* constant currents.  ground effect included. */
void efld( double xi, double yi, double zi, double ai, int ij )
{
#define	txk	egnd[0]
#define	tyk	egnd[1]
#define	tzk	egnd[2]
#define	txs	egnd[3]
#define	tys	egnd[4]
#define	tzs	egnd[5]
#define	txc	egnd[6]
#define	tyc	egnd[7]
#define	tzc	egnd[8]

  int ip;
  double xij, yij, ijx, rfl, salpr, zij, zp, rhox;
  double rhoy, rhoz, rh, r, rmag, cth, px, py;
  double xymag, xspec, yspec, rhospc, dmin, shaf;
  complex double epx, epy, refs, refps, zrsin, zratx, zscrn;
  complex double tezs, ters, tezc, terc, tezk, terk, egnd[9];

  xij= xi- xj;
  yij= yi- yj;
  ijx= ij;
  rfl=-1.;

  for( ip = 0; ip < ksymp; ip++ )
  {
    if( ip == 1)
      ijx=1;
    rfl=- rfl;
    salpr= salpj* rfl;
    zij= zi- rfl* zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;

    rh= sqrt( rhox* rhox+ rhoy* rhoy+ rhoz* rhoz+ ai* ai);
    if( rh <= 1.e-10)
    {
      rhox=0.;
      rhoy=0.;
      rhoz=0.;
    }
    else
    {
      rhox= rhox/ rh;
      rhoy= rhoy/ rh;
      rhoz= rhoz/ rh;
    }

    /* lumped current element approx. for large separations */
    r= sqrt( zp* zp+ rh* rh);
    if( r >= rkh)
    {
      rmag= TP* r;
      cth= zp/ r;
      px= rh/ r;
      txk= cmplx( cos( rmag),- sin( rmag));
      py= TP* r* r;
      tyk= ETA* cth* txk* cmplx(1.0,-1.0/ rmag)/ py;
      tzk= ETA* px* txk* cmplx(1.0, rmag-1.0/ rmag)/(2.* py);
      tezk= tyk* cth- tzk* px;
      terk= tyk* px+ tzk* cth;
      rmag= sin( PI* s)/ PI;
      tezc= tezk* rmag;
      terc= terk* rmag;
      tezk= tezk* s;
      terk= terk* s;
      txs=CPLX_00;
      tys=CPLX_00;
      tzs=CPLX_00;

    } /* if( r >= rkh) */

    if( r < rkh)
    {
      /* eksc for thin wire approx. or ekscx for extended t.w. approx. */
      if( iexk != 1)
	eksc( s, zp, rh, TP, ijx, &tezs, &ters,
	    &tezc, &terc, &tezk, &terk );
      else
	ekscx( b, s, zp, rh, TP, ijx, ind1, ind2,
	    &tezs, &ters, &tezc, &terc, &tezk, &terk);

      txs= tezs* cabj+ ters* rhox;
      tys= tezs* sabj+ ters* rhoy;
      tzs= tezs* salpr+ ters* rhoz;

    } /* if( r < rkh) */

    txk= tezk* cabj+ terk* rhox;
    tyk= tezk* sabj+ terk* rhoy;
    tzk= tezk* salpr+ terk* rhoz;
    txc= tezc* cabj+ terc* rhox;
    tyc= tezc* sabj+ terc* rhoy;
    tzc= tezc* salpr+ terc* rhoz;

    if( ip == 1)
    {
      if( iperf <= 0)
      {
	zratx= zrati;
	rmag= r;
	xymag= sqrt( xij* xij+ yij* yij);

	/* set parameters for radial wire ground screen. */
	if( nradl != 0)
	{
	  xspec=( xi* zj+ zi* xj)/( zi+ zj);
	  yspec=( yi* zj+ zi* yj)/( zi+ zj);
	  rhospc= sqrt( xspec* xspec+ yspec* yspec+ t2* t2);

	  if( rhospc <= scrwl)
	  {
	    zscrn= t1* rhospc* log( rhospc/ t2);
	    zratx=( zscrn* zrati)/( ETA* zrati+ zscrn);
	  }
	} /* if( nradl != 0) */

	/* calculation of reflection coefficients when ground is specified. */
	if( xymag <= 1.0e-6)
	{
	  px=0.;
	  py=0.;
	  cth=1.;
	  zrsin=CPLX_10;
	}
	else
	{
	  px=- yij/ xymag;
	  py= xij/ xymag;
	  cth= zij/ rmag;
	  zrsin= csqrt(1.0 - zratx*zratx*(1.0 - cth*cth) );

	} /* if( xymag <= 1.0e-6) */

	refs=( cth- zratx* zrsin)/( cth+ zratx* zrsin);
	refps=-( zratx* cth- zrsin)/( zratx* cth+ zrsin);
	refps= refps- refs;
	epy= px* txk+ py* tyk;
	epx= px* epy;
	epy= py* epy;
	txk= refs* txk+ refps* epx;
	tyk= refs* tyk+ refps* epy;
	tzk= refs* tzk;
	epy= px* txs+ py* tys;
	epx= px* epy;
	epy= py* epy;
	txs= refs* txs+ refps* epx;
	tys= refs* tys+ refps* epy;
	tzs= refs* tzs;
	epy= px* txc+ py* tyc;
	epx= px* epy;
	epy= py* epy;
	txc= refs* txc+ refps* epx;
	tyc= refs* tyc+ refps* epy;
	tzc= refs* tzc;

      } /* if( iperf <= 0) */

      exk= exk- txk* frati;
      eyk= eyk- tyk* frati;
      ezk= ezk- tzk* frati;
      exs= exs- txs* frati;
      eys= eys- tys* frati;
      ezs= ezs- tzs* frati;
      exc= exc- txc* frati;
      eyc= eyc- tyc* frati;
      ezc= ezc- tzc* frati;
      continue;

    } /* if( ip == 1) */

    exk= txk;
    eyk= tyk;
    ezk= tzk;
    exs= txs;
    eys= tys;
    ezs= tzs;
    exc= txc;
    eyc= tyc;
    ezc= tzc;

  } /* for( ip = 0; ip < ksymp; ip++ ) */

  if( iperf != 2)
    return;

  /* field due to ground using sommerfeld/norton */
  sn= sqrt( cabj* cabj+ sabj* sabj);
  if( sn >= 1.0e-5)
  {
    xsn= cabj/ sn;
    ysn= sabj/ sn;
  }
  else
  {
    sn=0.;
    xsn=1.;
    ysn=0.;
  }

  /* displace observation point for thin wire approximation */
  zij= zi+ zj;
  salpr=- salpj;
  rhox= sabj* zij- salpr* yij;
  rhoy= salpr* xij- cabj* zij;
  rhoz= cabj* yij- sabj* xij;
  rh= rhox* rhox+ rhoy* rhoy+ rhoz* rhoz;

  if( rh <= 1.e-10)
  {
    xo= xi- ai* ysn;
    yo= yi+ ai* xsn;
    zo= zi;
  }
  else
  {
    rh= ai/ sqrt( rh);
    if( rhoz < 0.)
      rh=- rh;
    xo= xi+ rh* rhox;
    yo= yi+ rh* rhoy;
    zo= zi+ rh* rhoz;

  } /* if( rh <= 1.e-10) */

  r= xij* xij+ yij* yij+ zij* zij;
  if( r <= .95)
  {
    /* field from interpolation is integrated over segment */
    isnor=1;
    dmin= exk* conj( exk)+ eyk* conj( eyk)+ ezk* conj( ezk);
    dmin=.01* sqrt( dmin);
    shaf=.5* s;
    rom2(- shaf, shaf, egnd, dmin);
  }
  else
  {
    /* norton field equations and lumped current element approximation */
    isnor=2;
    sflds(0., egnd);
  } /* if( r <= .95) */

  if( r > .95)
  {
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rh= r- zp* zp;
    if( rh <= 1.e-10)
      dmin=0.;
    else
      dmin= sqrt( rh/( rh+ ai* ai));

    if( dmin <= .95)
    {
      px=1.- dmin;
      terk=( txk* cabj+ tyk* sabj+ tzk* salpr)* px;
      txk= dmin* txk+ terk* cabj;
      tyk= dmin* tyk+ terk* sabj;
      tzk= dmin* tzk+ terk* salpr;
      ters=( txs* cabj+ tys* sabj+ tzs* salpr)* px;
      txs= dmin* txs+ ters* cabj;
      tys= dmin* tys+ ters* sabj;
      tzs= dmin* tzs+ ters* salpr;
      terc=( txc* cabj+ tyc* sabj+ tzc* salpr)* px;
      txc= dmin* txc+ terc* cabj;
      tyc= dmin* tyc+ terc* sabj;
      tzc= dmin* tzc+ terc* salpr;

    } /* if( dmin <= .95) */

  } /* if( r > .95) */

  exk= exk+ txk;
  eyk= eyk+ tyk;
  ezk= ezk+ tzk;
  exs= exs+ txs;
  eys= eys+ tys;
  ezs= ezs+ tzs;
  exc= exc+ txc;
  eyc= eyc+ tyc;
  ezc= ezc+ tzc;

  return;
}

/*-----------------------------------------------------------------------*/

/* compute e field of sine, cosine, and constant */
/* current filaments by thin wire approximation. */
void eksc( double s, double z, double rh, double xk, int ij,
    complex double *ezs, complex double *ers, complex double *ezc,
    complex double *erc, complex double *ezk, complex double *erk )
{
  double rhk, sh, shk, ss, cs, z1a, z2a, cint, sint;
  complex double gz1, gz2, gp1, gp2, gzp1, gzp2;

  ija= ij;
  zpk= xk* z;
  rhk= xk* rh;
  rkb2= rhk* rhk;
  sh=.5* s;
  shk= xk* sh;
  ss= sin( shk);
  cs= cos( shk);
  z2a= sh- z;
  z1a=-( sh+ z);
  gx( z1a, rh, xk, &gz1, &gp1);
  gx( z2a, rh, xk, &gz2, &gp2);
  gzp1= gp1* z1a;
  gzp2= gp2* z2a;
  *ezs=  CONST1*(( gz2- gz1)* cs* xk-( gzp2+ gzp1)* ss);
  *ezc=- CONST1*(( gz2+ gz1)* ss* xk+( gzp2- gzp1)* cs);
  *erk= CONST1*( gp2- gp1)* rh;
  intx(- shk, shk, rhk, ij, &cint, &sint);
  *ezk=- CONST1*( gzp2- gzp1+ xk* xk* cmplx( cint,- sint));
  gzp1= gzp1* z1a;
  gzp2= gzp2* z2a;

  if( rh >= 1.0e-10)
  {
    *ers=- CONST1*(( gzp2+ gzp1+ gz2+ gz1)*
	ss-( z2a* gz2- z1a* gz1)* cs*xk)/ rh;
    *erc=- CONST1*(( gzp2- gzp1+ gz2- gz1)*
	cs+( z2a* gz2+ z1a* gz1)* ss*xk)/ rh;
    return;
  }

  *ers = CPLX_00;
  *erc = CPLX_00;

  return;
}

/*-----------------------------------------------------------------------*/

/* compute e field of sine, cosine, and constant current */
/* filaments by extended thin wire approximation. */
void ekscx( double bx, double s, double z,
    double rhx, double xk, int ij, int inx1, int inx2,
    complex double *ezs, complex double *ers, complex double *ezc,
    complex double *erc, complex double *ezk, complex double *erk )
{
  int ira;
  double b, rh, sh, rhk, shk, ss, cs, z1a;
  double z2a, a2, bk, bk2, cint, sint;
  complex double gz1, gz2, gzp1, gzp2, gr1, gr2;
  complex double grp1, grp2, grk1, grk2, gzz1, gzz2;

  if( rhx >= bx)
  {
    rh= rhx;
    b= bx;
    ira=0;
  }
  else
  {
    rh= bx;
    b= rhx;
    ira=1;
  }

  sh=.5* s;
  ija= ij;
  zpk= xk* z;
  rhk= xk* rh;
  rkb2= rhk* rhk;
  shk= xk* sh;
  ss= sin( shk);
  cs= cos( shk);
  z2a= sh- z;
  z1a=-( sh+ z);
  a2= b* b;

  if( inx1 != 2)
    gxx( z1a, rh, b, a2, xk, ira, &gz1,
	&gzp1, &gr1, &grp1, &grk1, &gzz1);
  else
  {
    gx( z1a, rhx, xk, &gz1, &grk1);
    gzp1= grk1* z1a;
    gr1= gz1/ rhx;
    grp1= gzp1/ rhx;
    grk1= grk1* rhx;
    gzz1= CPLX_00;
  }

  if( inx2 != 2)
    gxx( z2a, rh, b, a2, xk, ira, &gz2,
	&gzp2, &gr2, &grp2, &grk2, &gzz2);
  else
  {
    gx( z2a, rhx, xk, &gz2, &grk2);
    gzp2= grk2* z2a;
    gr2= gz2/ rhx;
    grp2= gzp2/ rhx;
    grk2= grk2* rhx;
    gzz2= CPLX_00;
  }

  *ezs= CONST1*(( gz2- gz1)* cs* xk-( gzp2+ gzp1)* ss);
  *ezc=- CONST1*(( gz2+ gz1)* ss* xk+( gzp2- gzp1)* cs);
  *ers=- CONST1*(( z2a* grp2+ z1a* grp1+ gr2+ gr1)*ss
      -( z2a* gr2- z1a* gr1)* cs* xk);
  *erc=- CONST1*(( z2a* grp2- z1a* grp1+ gr2- gr1)*cs
      +( z2a* gr2+ z1a* gr1)* ss* xk);
  *erk= CONST1*( grk2- grk1);
  intx(- shk, shk, rhk, ij, &cint, &sint);
  bk= b* xk;
  bk2= bk* bk*.25;
  *ezk=- CONST1*( gzp2- gzp1+ xk* xk*(1.- bk2)*
      cmplx( cint,- sint)-bk2*( gzz2- gzz1));

  return;
}

/*-----------------------------------------------------------------------*/

/* etmns fills the array e with the negative of the */
/* electric field incident on the structure. e is the */
/* right hand side of the matrix equation. */
void etmns( double p1, double p2, double p3, double p4,
    double p5, double p6, int ipr, complex double *e )
{
  int i, is, i1, i2=0, neq;
  double cth, sth, cph, sph, cet, set, pxl, pyl, pzl, wx;
  double wy, wz, qx, qy, qz, arg, ds, dsh, rs, r;
  complex double cx, cy, cz, er, et, ezh, erh, rrv, rrh, tt1, tt2;

  neq= n+2*m;
  nqds=0;

  /* applied field of voltage sources for transmitting case */
  if( (ipr <= 0) || (ipr == 5) )
  {
    for( i = 0; i < neq; i++ )
      e[i]=CPLX_00;

    if( nsant != 0)
    {
      for( i = 0; i < nsant; i++ )
      {
	is= isant[i]-1;
	e[is]= -vsant[i]/( si[is]* wlam);
      }
    }

    if( nvqd == 0)
      return;

    for( i = 0; i < nvqd; i++ )
    {
      is= ivqd[i];
      qdsrc( is, vqd[i], e);
    }
    return;

  } /* if( (ipr <= 0) || (ipr == 5) ) */

  /* incident plane wave, linearly polarized. */
  if( ipr <= 3)
  {
    cth= cos( p1);
    sth= sin( p1);
    cph= cos( p2);
    sph= sin( p2);
    cet= cos( p3);
    set= sin( p3);
    pxl= cth* cph* cet- sph* set;
    pyl= cth* sph* cet+ cph* set;
    pzl=- sth* cet;
    wx=- sth* cph;
    wy=- sth* sph;
    wz=- cth;
    qx= wy* pzl- wz* pyl;
    qy= wz* pxl- wx* pzl;
    qz= wx* pyl- wy* pxl;

    if( ksymp != 1)
    {
      if( iperf != 1)
      {
	rrv= csqrt(1.- zrati* zrati* sth* sth);
	rrh= zrati* cth;
	rrh=( rrh- rrv)/( rrh+ rrv);
	rrv= zrati* rrv;
	rrv=-( cth- rrv)/( cth+ rrv);
      }
      else
      {
	rrv=-CPLX_10;
	rrh=-CPLX_10;
      } /* if( iperf != 1) */

    } /* if( ksymp != 1) */

    if( ipr <= 1)
    {
      if( n != 0)
      {
	for( i = 0; i < n; i++ )
	{
	  arg=- TP*( wx* x[i]+ wy* y[i]+ wz* z[i]);
	  e[i]=-( pxl* cab[i]+ pyl* sab[i]+ pzl*
	      salp[i])* cmplx( cos( arg), sin( arg));
	}

	if( ksymp != 1)
	{
	  tt1=( pyl* cph- pxl* sph)*( rrh- rrv);
	  cx= rrv* pxl- tt1* sph;
	  cy= rrv* pyl+ tt1* cph;
	  cz=- rrv* pzl;

	  for( i = 0; i < n; i++ )
	  {
	    arg=- TP*( wx* x[i]+ wy* y[i]- wz* z[i]);
	    e[i]= e[i]-( cx* cab[i]+ cy* sab[i]+
		cz* salp[i])* cmplx(cos( arg), sin( arg));
	  }

	} /* if( ksymp != 1) */

      } /* if( n != 0) */

      if( m == 0)
	return;

      i= -1;
      i1= n-2;
      for( is = 0; is < m; is++ )
      {
	i++;
	i1 += 2;
	i2 = i1+1;
	arg=- TP*( wx* px[i]+ wy* py[i]+ wz* pz[i]);
	tt1= cmplx( cos( arg), sin( arg))* psalp[i]* RETA;
	e[i2]=( qx* t1x[i]+ qy* t1y[i]+ qz* t1z[i])* tt1;
	e[i1]=( qx* t2x[i]+ qy* t2y[i]+ qz* t2z[i])* tt1;
      }

      if( ksymp == 1)
	return;

      tt1=( qy* cph- qx* sph)*( rrv- rrh);
      cx=-( rrh* qx- tt1* sph);
      cy=-( rrh* qy+ tt1* cph);
      cz= rrh* qz;

      i= -1;
      i1= n-2;
      for( is = 0; is < m; is++ )
      {
	i++;
	i1 += 2;
	i2 = i1+1;
	arg=- TP*( wx* px[i]+ wy* py[i]- wz* pz[i]);
	tt1= cmplx( cos( arg), sin( arg))* psalp[i]* RETA;
	e[i2]= e[i2]+( cx* t1x[i]+ cy* t1y[i]+ cz* t1z[i])* tt1;
	e[i1]= e[i1]+( cx* t2x[i]+ cy* t2y[i]+ cz* t2z[i])* tt1;
      }
      return;

    } /* if( ipr <= 1) */

    /* incident plane wave, elliptic polarization. */
    tt1=-(CPLX_01)* p6;
    if( ipr == 3)
      tt1=- tt1;

    if( n != 0)
    {
      cx= pxl+ tt1* qx;
      cy= pyl+ tt1* qy;
      cz= pzl+ tt1* qz;

      for( i = 0; i < n; i++ )
      {
	arg=- TP*( wx* x[i]+ wy* y[i]+ wz* z[i]);
	e[i]=-( cx* cab[i]+ cy* sab[i]+ cz*
	    salp[i])* cmplx( cos( arg), sin( arg));
      }

      if( ksymp != 1)
      {
	tt2=( cy* cph- cx* sph)*( rrh- rrv);
	cx= rrv* cx- tt2* sph;
	cy= rrv* cy+ tt2* cph;
	cz=- rrv* cz;

	for( i = 0; i < n; i++ )
	{
	  arg=- TP*( wx* x[i]+ wy* y[i]- wz* z[i]);
	  e[i]= e[i]-( cx* cab[i]+ cy* sab[i]+
	      cz* salp[i])* cmplx(cos( arg), sin( arg));
	}

      } /* if( ksymp != 1) */

    } /* if( n != 0) */

    if( m == 0)
      return;

    cx= qx- tt1* pxl;
    cy= qy- tt1* pyl;
    cz= qz- tt1* pzl;

    i= -1;
    i1= n-2;
    for( is = 0; is < m; is++ )
    {
      i++;
      i1 += 2;
      i2 = i1+1;
      arg=- TP*( wx* px[i]+ wy* py[i]+ wz* pz[i]);
      tt2= cmplx( cos( arg), sin( arg))* psalp[i]* RETA;
      e[i2]=( cx* t1x[i]+ cy* t1y[i]+ cz* t1z[i])* tt2;
      e[i1]=( cx* t2x[i]+ cy* t2y[i]+ cz* t2z[i])* tt2;
    }

    if( ksymp == 1)
      return;

    tt1=( cy* cph- cx* sph)*( rrv- rrh);
    cx=-( rrh* cx- tt1* sph);
    cy=-( rrh* cy+ tt1* cph);
    cz= rrh* cz;

    i= -1;
    i1= n-2;
    for( is=0; is < m; is++ )
    {
      i++;
      i1 += 2;
      i2 = i1+1;
      arg=- TP*( wx* px[i]+ wy* py[i]- wz* pz[i]);
      tt1= cmplx( cos( arg), sin( arg))* psalp[i]* RETA;
      e[i2]= e[i2]+( cx* t1x[i]+ cy* t1y[i]+ cz* t1z[i])* tt1;
      e[i1]= e[i1]+( cx* t2x[i]+ cy* t2y[i]+ cz* t2z[i])* tt1;
    }

    return;

  } /* if( ipr <= 3) */

  /* incident field of an elementary current source. */
  wz= cos( p4);
  wx= wz* cos( p5);
  wy= wz* sin( p5);
  wz= sin( p4);
  ds= p6*59.958;
  dsh= p6/(2.* TP);

  is= 0;
  i1= n-2;
  for( i = 0; i < npm; i++ )
  {
    if( i >= n )
    {
      i1 += 2;
      i2 = i1+1;
      pxl= px[is]- p1;
      pyl= py[is]- p2;
      pzl= pz[is]- p3;
      is++;
    }

    pxl= x[i]- p1;
    pyl= y[i]- p2;
    pzl= z[i]- p3;

      rs= pxl* pxl+ pyl* pyl+ pzl* pzl;
    if( rs < 1.0e-30)
      continue;

    r= sqrt( rs);
    pxl= pxl/ r;
    pyl= pyl/ r;
    pzl= pzl/ r;
    cth= pxl* wx+ pyl* wy+ pzl* wz;
    sth= sqrt(1.- cth* cth);
    qx= pxl- wx* cth;
    qy= pyl- wy* cth;
    qz= pzl- wz* cth;

    arg= sqrt( qx* qx+ qy* qy+ qz* qz);
    if( arg >= 1.e-30)
    {
      qx= qx/ arg;
      qy= qy/ arg;
      qz= qz/ arg;
    }
    else
    {
      qx=1.;
      qy=0.;
      qz=0.;

    } /* if( arg >= 1.e-30) */

    arg=- TP* r;
    tt1= cmplx( cos( arg), sin( arg));

    if( i < n )
    {
      tt2= cmplx(1.0,-1.0/( r* TP))/ rs;
      er= ds* tt1* tt2* cth;
      et=.5* ds* tt1*((CPLX_01)* TP/ r+ tt2)* sth;
      ezh= er* cth- et* sth;
      erh= er* sth+ et* cth;
      cx= ezh* wx+ erh* qx;
      cy= ezh* wy+ erh* qy;
      cz= ezh* wz+ erh* qz;
      e[i]=-( cx* cab[i]+ cy* sab[i]+ cz* salp[i]);
    }
    else
    {
      pxl= wy* qz- wz* qy;
      pyl= wz* qx- wx* qz;
      pzl= wx* qy- wy* qx;
      tt2= dsh* tt1* cmplx(1./ r, TP)/ r* sth* psalp[is];
      cx= tt2* pxl;
      cy= tt2* pyl;
      cz= tt2* pzl;
      e[i2]= cx* t1x[is]+ cy* t1y[is]+ cz* t1z[is];
      e[i1]= cx* t2x[is]+ cy* t2y[is]+ cz* t2z[is];

    } /* if( i >= n) */

  } /* for( i = 0; i < npm; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine to factor a matrix into a unit lower triangular matrix */
/* and an upper triangular matrix using the gauss-doolittle algorithm */
/* presented on pages 411-416 of a. ralston--a first course in */
/* numerical analysis.  comments below refer to comments in ralstons */
/* text.    (matrix transposed.) */

void factr( int n, complex double *a, int *ip, int ndim)
{
  int r, rm1, rp1, pj, pr, iflg, k, j, jp1, i;
  double dmax, elmag;
  complex double arj, *scm = NULL;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, np2m * sizeof(complex double) );

  /* Un-transpose the matrix for Gauss elimination */
  for( i = 1; i < n; i++ )
    for( j = 0; j < i; j++ )
    {
      arj = a[i+j*ndim];
      a[i+j*ndim] = a[j+i*ndim];
      a[j+i*ndim] = arj;
    }

  iflg=FALSE;
  /* step 1 */
  for( r = 0; r < n; r++ )
  {
    for( k = 0; k < n; k++ )
      scm[k]= a[k+r*ndim];

    /* steps 2 and 3 */
    rm1= r;
    if( rm1 > 0)
    {
      for( j = 0; j < rm1; j++ )
      {
	pj= ip[j]-1;
	arj= scm[pj];
	a[j+r*ndim]= arj;
	scm[pj]= scm[j];
	jp1= j+1;

	for( i = jp1; i < n; i++ )
	  scm[i] -= a[i+j*ndim]* arj;

      } /* for( j = 0; j < rm1; j++ ) */

    } /* if( rm1 >= 0.) */

    /* step 4 */
    dmax= creal( scm[r]*conj(scm[r]) );

    rp1= r+1;
    ip[r]= rp1;
    if( rp1 < n)
    {
      for( i = rp1; i < n; i++ )
      {
	elmag= creal( scm[i]* conj(scm[i]) );
	if( elmag >= dmax)
	{
	  dmax= elmag;
	  ip[r]= i+1;
	}
      }
    } /* if( rp1 < n) */

    if( dmax < 1.e-10)
      iflg=TRUE;

    pr= ip[r]-1;
    a[r+r*ndim]= scm[pr];
    scm[pr]= scm[r];

    /* step 5 */
    if( rp1 < n)
    {
      arj=1./ a[r+r*ndim];

      for( i = rp1; i < n; i++ )
	a[i+r*ndim]= scm[i]* arj;
    }

    if( iflg == TRUE )
    {
      fprintf( output_fp,
	  "\n  PIVOT(%d)= %16.8E", r, dmax );
      iflg=FALSE;
    }

  } /* for( r=0; r < n; r++ ) */

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* factrs, for symmetric structure, transforms submatricies to form */
/* matricies of the symmetric modes and calls routine to factor */
/* matricies.  if no symmetry, the routine is called to factor the */
/* complete matrix. */
void factrs( int np, int nrow, complex double *a, int *ip )
{
  int kk, ka;

  for( kk = 0; kk < nop; kk++ )
  {
    ka= kk* np;
    factr( np, &a[ka], &ip[ka], nrow );
  }
  return;
}

/*-----------------------------------------------------------------------*/

/* fbar is sommerfeld attenuation function for numerical distance p */
complex double  fbar( complex double p )
{
  int i, minus;
  double tms, sms;
  complex double z, zs, sum, pow, term, fbar;

  z= CPLX_01* csqrt( p);
  if( cabs( z) <= 3.)
  {
    /* series expansion */
    zs= z* z;
    sum= z;
    pow= z;

    for( i = 1; i <= 100; i++ )
    {
      pow=- pow* zs/ (double)i;
      term= pow/(2.* i+1.);
      sum= sum+ term;
      tms= creal( term* conj( term));
      sms= creal( sum* conj( sum));
      if( tms/sms < ACCS)
	break;
    }

    fbar=1.-(1.- sum* TOSP)* z* cexp( zs)* SP;
    return( fbar );

  } /* if( cabs( z) <= 3.) */

  /* asymptotic expansion */
  if( creal( z) < 0.)
  {
    minus=1;
    z=- z;
  }
  else
    minus=0;

  zs=.5/( z* z);
  sum=CPLX_00;
  term=CPLX_10;

  for( i = 1; i <= 6; i++ )
  {
    term =- term*(2.*i -1.)* zs;
    sum += term;
  }

  if( minus == 1)
    sum -= 2.* SP* z* cexp( z* z);
  fbar=- sum;

  return( fbar );
}

/*-----------------------------------------------------------------------*/

/* fblock sets parameters for out-of-core */
/* solution for the primary matrix (a) */
void fblock( int nrow, int ncol, int imax, int ipsym )
{
  int i, j, k, ka, kk;
  double phaz, arg;
  complex double deter;

  if( nrow*ncol <= imax)
  {
    npblk= nrow;
    nlast= nrow;
    imat= nrow* ncol;

    if( nrow == ncol)
    {
      icase=1;
      return;
    }
    else
      icase=2;

  } /* if( nrow*ncol <= imax) */

  if( nop*nrow != ncol)
  {
    fprintf( output_fp,
	"\n  SYMMETRY ERROR - NROW: %d NCOL: %d", nrow, ncol );
    stop(-1);
  }

  /* set up ssx matrix for rotational symmetry. */
  if( ipsym <= 0)
  {
    phaz = TP/nop;

    for( i = 1; i < nop; i++ )
    {
      for( j= i; j < nop; j++ )
      {
	arg= phaz* (double)i * (double)j;
	ssx[i+j*nop]= cmplx( cos( arg), sin( arg));
	ssx[j+i*nop]= ssx[i+j*nop];
      }
    }
    return;

  } /* if( ipsym <= 0) */

  /* set up ssx matrix for plane symmetry */
  kk=1;
  ssx[0]=CPLX_10;

  k = 2;
  for( ka = 1; k != nop; ka++ )
    k *= 2;

  for( k = 0; k < ka; k++ )
  {
    for( i = 0; i < kk; i++ )
    {
      for( j = 0; j < kk; j++ )
      {
	deter= ssx[i+j*nop];
	ssx[i+(j+kk)*nop]= deter;
	ssx[i+kk+(j+kk)*nop]=- deter;
	ssx[i+kk+j*nop]= deter;
      }
    }
    kk *= 2;

  } /* for( k = 0; k < ka; k++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* ffld calculates the far zone radiated electric fields, */
/* the factor exp(j*k*r)/(r/lamda) not included */
void ffld( double thet, double phi,
    complex double *eth, complex double *eph )
{
  int k, i, ip, jump;
  double phx, phy, roz, rozs, thx, thy, thz, rox, roy;
  double tthet=0., darg=0., omega, el, sill, top, bot, a;
  double too, boo, b, c, d, rr, ri, arg, dr, rfl, rrz;
  complex double cix, ciy, ciz, exa, ccx, ccy, ccz, cdp;
  complex double zrsin, rrv, rrh, rrv1, rrh1, rrv2, rrh2;
  complex double tix, tiy, tiz, zscrn, ex, ey, ez, gx, gy, gz;

  phx=- sin( phi);
  phy= cos( phi);
  roz= cos( thet);
  rozs= roz;
  thx= roz* phy;
  thy=- roz* phx;
  thz=- sin( thet);
  rox=- thz* phy;
  roy= thz* phx;

  jump = FALSE;
  if( n != 0)
  {
    /* loop for structure image if any */
    /* calculation of reflection coeffecients */
    for( k = 0; k < ksymp; k++ )
    {
      if( k != 0 )
      {
	/* for perfect ground */
	if( iperf == 1)
	{
	  rrv=-CPLX_10;
	  rrh=-CPLX_10;
	}
	else
	{
	  /* for infinite planar ground */
	  zrsin= csqrt(1.- zrati* zrati* thz* thz);
	  rrv=-( roz- zrati* zrsin)/( roz+ zrati* zrsin);
	  rrh=( zrati* roz- zrsin)/( zrati* roz+ zrsin);

	} /* if( iperf == 1) */

	/* for the cliff problem, two reflction coefficients calculated */
	if( ifar > 1)
	{
	  rrv1= rrv;
	  rrh1= rrh;
	  tthet= tan( thet);

	  if( ifar != 4)
	  {
	    zrsin= csqrt(1.- zrati2* zrati2* thz* thz);
	    rrv2=-( roz- zrati2* zrsin)/( roz+ zrati2* zrsin);
	    rrh2=( zrati2* roz- zrsin)/( zrati2* roz+ zrsin);
	    darg=- TP*2.* ch* roz;
	  }
	} /* if( ifar > 1) */

	roz=- roz;
	ccx= cix;
	ccy= ciy;
	ccz= ciz;

      } /* if( k != 0 ) */

      cix=CPLX_00;
      ciy=CPLX_00;
      ciz=CPLX_00;

      /* loop over structure segments */
      for( i = 0; i < n; i++ )
      {
	omega=-( rox* cab[i]+ roy* sab[i]+ roz* salp[i]);
	el= PI* si[i];
	sill= omega* el;
	top= el+ sill;
	bot= el- sill;

	if( fabs( omega) >= 1.0e-7)
	  a=2.* sin( sill)/ omega;
	else
	  a=(2.- omega* omega* el* el/3.)* el;

	if( fabs( top) >= 1.0e-7)
	  too= sin( top)/ top;
	else
	  too=1.- top* top/6.;

	if( fabs( bot) >= 1.0e-7)
	  boo= sin( bot)/ bot;
	else
	  boo=1.- bot* bot/6.;

	b= el*( boo- too);
	c= el*( boo+ too);
	rr= a* air[i]+ b* bii[i]+ c* cir[i];
	ri= a* aii[i]- b* bir[i]+ c* cii[i];
	arg= TP*( x[i]* rox+ y[i]* roy+ z[i]* roz);

	if( (k != 1) || (ifar < 2) )
	{
	  /* summation for far field integral */
	  exa= cmplx( cos( arg), sin( arg))* cmplx( rr, ri);
	  cix= cix+ exa* cab[i];
	  ciy= ciy+ exa* sab[i];
	  ciz= ciz+ exa* salp[i];
	  continue;
	}

	/* calculation of image contribution */
	/* in cliff and ground screen problems */

	/* specular point distance */
	dr= z[i]* tthet;

	d= dr* phy+ x[i];
	if( ifar == 2)
	{
	  if(( cl- d) > 0.)
	  {
	    rrv= rrv1;
	    rrh= rrh1;
	  }
	  else
	  {
	    rrv= rrv2;
	    rrh= rrh2;
	    arg= arg+ darg;
	  }
	} /* if( ifar == 2) */
	else
	{
	  d= sqrt( d*d + (y[i]-dr*phx)*(y[i]-dr*phx) );
	  if( ifar == 3)
	  {
	    if(( cl- d) > 0.)
	    {
	      rrv= rrv1;
	      rrh= rrh1;
	    }
	    else
	    {
	      rrv= rrv2;
	      rrh= rrh2;
	      arg= arg+ darg;
	    }
	  } /* if( ifar == 3) */
	  else
	  {
	    if(( scrwl- d) >= 0.)
	    {
	      /* radial wire ground screen reflection coefficient */
	      d= d+ t2;
	      zscrn= t1* d* log( d/ t2);
	      zscrn=( zscrn* zrati)/( ETA* zrati+ zscrn);
	      zrsin= csqrt(1.- zscrn* zscrn* thz* thz);
	      rrv=( roz+ zscrn* zrsin)/(- roz+ zscrn* zrsin);
	      rrh=( zscrn* roz+ zrsin)/( zscrn* roz- zrsin);
	    } /* if(( scrwl- d) < 0.) */
	    else
	    {
	      if( ifar == 4)
	      {
		rrv= rrv1;
		rrh= rrh1;
	      } /* if( ifar == 4) */
	      else
	      {
		if( ifar == 5)
		  d= dr* phy+ x[i];

		if(( cl- d) > 0.)
		{
		  rrv= rrv1;
		  rrh= rrh1;
		}
		else
		{
		  rrv= rrv2;
		  rrh= rrh2;
		  arg= arg+ darg;
		} /* if(( cl- d) > 0.) */

	      } /* if( ifar == 4) */

	    } /* if(( scrwl- d) < 0.) */

	  } /* if( ifar == 3) */

	} /* if( ifar == 2) */

	/* contribution of each image segment modified by */
	/* reflection coef, for cliff and ground screen problems */
	exa= cmplx( cos( arg), sin( arg))* cmplx( rr, ri);
	tix= exa* cab[i];
	tiy= exa* sab[i];
	tiz= exa* salp[i];
	cdp=( tix* phx+ tiy* phy)*( rrh- rrv);
	cix= cix+ tix* rrv+ cdp* phx;
	ciy= ciy+ tiy* rrv+ cdp* phy;
	ciz= ciz- tiz* rrv;

      } /* for( i = 0; i < n; i++ ) */

      if( k == 0 )
	continue;

      /* calculation of contribution of structure image for infinite ground */
      if( ifar < 2)
      {
	cdp=( cix* phx+ ciy* phy)*( rrh- rrv);
	cix= ccx+ cix* rrv+ cdp* phx;
	ciy= ccy+ ciy* rrv+ cdp* phy;
	ciz= ccz- ciz* rrv;
      }
      else
      {
	cix= cix+ ccx;
	ciy= ciy+ ccy;
	ciz= ciz+ ccz;
      }

    } /* for( k=0; k < ksymp; k++ ) */

    if( m > 0)
      jump = TRUE;
    else
    {
      *eth=( cix* thx+ ciy* thy+ ciz* thz)* CONST3;
      *eph=( cix* phx+ ciy* phy)* CONST3;
      return;
    }

  } /* if( n != 0) */

  if( ! jump )
  {
    cix=CPLX_00;
    ciy=CPLX_00;
    ciz=CPLX_00;
  }

  /* electric field components */
  roz= rozs;
  rfl=-1.;
  for( ip = 0; ip < ksymp; ip++ )
  {
    rfl=- rfl;
    rrz= roz* rfl;
    fflds( rox, roy, rrz, &cur[n], &gx, &gy, &gz);

    if( ip != 1 )
    {
      ex= gx;
      ey= gy;
      ez= gz;
      continue;
    }

    if( iperf == 1)
    {
      gx=- gx;
      gy=- gy;
      gz=- gz;
    }
    else
    {
      rrv= csqrt(1.- zrati* zrati* thz* thz);
      rrh= zrati* roz;
      rrh=( rrh- rrv)/( rrh+ rrv);
      rrv= zrati* rrv;
      rrv=-( roz- rrv)/( roz+ rrv);
      *eth=( gx* phx+ gy* phy)*( rrh- rrv);
      gx= gx* rrv+ *eth* phx;
      gy= gy* rrv+ *eth* phy;
      gz= gz* rrv;

    } /* if( iperf == 1) */

    ex= ex+ gx;
    ey= ey+ gy;
    ez= ez- gz;

  } /* for( ip = 0; ip < ksymp; ip++ ) */

  ex= ex+ cix* CONST3;
  ey= ey+ ciy* CONST3;
  ez= ez+ ciz* CONST3;
  *eth= ex* thx+ ey* thy+ ez* thz;
  *eph= ex* phx+ ey* phy;

  return;
}

/*-----------------------------------------------------------------------*/

/* calculates the xyz components of the electric */
/* field due to surface currents */
void fflds( double rox, double roy, double roz,
    complex double *scur, complex double *ex,
    complex double *ey, complex double *ez )
{
  double *xs, *ys, *zs, *s;
  int j, i, k;
  double arg;
  complex double ct;

  xs = px; ys = py; zs = pz; s = pbi;
  *ex=CPLX_00;
  *ey=CPLX_00;
  *ez=CPLX_00;

  i= -1;
  for( j = 0; j < m; j++ )
  {
    i++;
    arg= TP*( rox* xs[i]+ roy* ys[i]+ roz* zs[i]);
    ct= cmplx( cos( arg)* s[i], sin( arg)* s[i]);
    k=3*(j+1)-1;
    *ex= *ex+ scur[k-2]* ct;
    *ey= *ey+ scur[k-1]* ct;
    *ez= *ez+ scur[k  ]* ct;
  }

  ct= rox* *ex+ roy* *ey+ roz* *ez;
  *ex= CONST4*( ct* rox- *ex);
  *ey= CONST4*( ct* roy- *ey);
  *ez= CONST4*( ct* roz- *ez);

  return;
}

/*-----------------------------------------------------------------------*/

/* gf computes the integrand exp(jkr)/(kr) for numerical integration. */
void gf( double zk, double *co, double *si )
{
  double zdk, rk, rks;

  zdk= zk- zpk;
  rk= sqrt( rkb2+ zdk* zdk);
  *si= sin( rk)/ rk;

  if( ija != 0 )
  {
    *co= cos( rk)/ rk;
    return;
  }

  if( rk >= .2)
  {
    *co=( cos( rk)-1.)/ rk;
    return;
  }

  rks= rk* rk;
  *co=((-1.38888889e-3* rks+4.16666667e-2)* rks-.5)* rk;

  return;
}

/*-----------------------------------------------------------------------*/

/* gfld computes the radiated field including ground wave. */
void gfld( double rho, double phi, double rz,
    complex double *eth, complex double *epi,
    complex double *erd, complex double ux, int ksymp )
{
  int i, k;
  double b, r, thet, arg, phx, phy, rx, ry, dx, dy, dz, rix, riy, rhs, rhp;
  double rhx, rhy, calp, cbet, sbet, cph, sph, el, rfl, riz, thx, thy, thz;
  double rxyz, rnx, rny, rnz, omega, sill, top, bot, a, too, boo, c, rr, ri;
  complex double cix, ciy, ciz, exa, erv;
  complex double ezv, erh, eph, ezh, ex, ey;

  r= sqrt( rho*rho+ rz*rz );
  if( (ksymp == 1) || (cabs(ux) > .5) || (r > 1.e5) )
  {
    /* computation of space wave only */
    if( rz >= 1.0e-20)
      thet= atan( rho/ rz);
    else
      thet= PI*.5;

    ffld( thet, phi, eth, epi);
    arg=- TP* r;
    exa= cmplx( cos( arg), sin( arg))/ r;
    *eth= *eth* exa;
    *epi= *epi* exa;
    *erd=CPLX_00;
    return;
  } /* if( (ksymp == 1) && (cabs(ux) > .5) && (r > 1.e5) ) */

  /* computation of space and ground waves. */
  u= ux;
  u2= u* u;
  phx=- sin( phi);
  phy= cos( phi);
  rx= rho* phy;
  ry=- rho* phx;
  cix=CPLX_00;
  ciy=CPLX_00;
  ciz=CPLX_00;

  /* summation of field from individual segments */
  for( i = 0; i < n; i++ )
  {
    dx= cab[i];
    dy= sab[i];
    dz= salp[i];
    rix= rx- x[i];
    riy= ry- y[i];
    rhs= rix* rix+ riy* riy;
    rhp= sqrt( rhs);

    if( rhp >= 1.0e-6)
    {
      rhx= rix/ rhp;
      rhy= riy/ rhp;
    }
    else
    {
      rhx=1.;
      rhy=0.;
    }

    calp=1.- dz* dz;
    if( calp >= 1.0e-6)
    {
      calp= sqrt( calp);
      cbet= dx/ calp;
      sbet= dy/ calp;
      cph= rhx* cbet+ rhy* sbet;
      sph= rhy* cbet- rhx* sbet;
    }
    else
    {
      cph= rhx;
      sph= rhy;
    }

    el= PI* si[i];
    rfl=-1.;

    /* integration of (current)*(phase factor) over segment and image for */
    /* constant, sine, and cosine current distributions */
    for( k = 0; k < 2; k++ )
    {
      rfl=- rfl;
      riz= rz- z[i]* rfl;
      rxyz= sqrt( rix* rix+ riy* riy+ riz* riz);
      rnx= rix/ rxyz;
      rny= riy/ rxyz;
      rnz= riz/ rxyz;
      omega=-( rnx* dx+ rny* dy+ rnz* dz* rfl);
      sill= omega* el;
      top= el+ sill;
      bot= el- sill;

      if( fabs( omega) >= 1.0e-7)
	a=2.* sin( sill)/ omega;
      else
	a=(2.- omega* omega* el* el/3.)* el;

      if( fabs( top) >= 1.0e-7)
	too= sin( top)/ top;
      else
	too=1.- top* top/6.;

      if( fabs( bot) >= 1.0e-7)
	boo= sin( bot)/ bot;
      else
	boo=1.- bot* bot/6.;

      b= el*( boo- too);
      c= el*( boo+ too);
      rr= a* air[i]+ b* bii[i]+ c* cir[i];
      ri= a* aii[i]- b* bir[i]+ c* cii[i];
      arg= TP*( x[i]* rnx+ y[i]* rny+ z[i]* rnz* rfl);
      exa= cmplx( cos( arg), sin( arg))* cmplx( rr, ri)/ TP;

      if( k != 1 )
      {
	xx1= exa;
	r1= rxyz;
	zmh= riz;
	continue;
      }

      xx2= exa;
      r2= rxyz;
      zph= riz;

    } /* for( k = 0; k < 2; k++ ) */

    /* call subroutine to compute the field */
    /* of segment including ground wave. */
    gwave( &erv, &ezv, &erh, &ezh, &eph);
    erh= erh* cph* calp+ erv* dz;
    eph= eph* sph* calp;
    ezh= ezh* cph* calp+ ezv* dz;
    ex= erh* rhx- eph* rhy;
    ey= erh* rhy+ eph* rhx;
    cix= cix+ ex;
    ciy= ciy+ ey;
    ciz= ciz+ ezh;

  } /* for( i = 0; i < n; i++ ) */

  arg=- TP* r;
  exa= cmplx( cos( arg), sin( arg));
  cix= cix* exa;
  ciy= ciy* exa;
  ciz= ciz* exa;
  rnx= rx/ r;
  rny= ry/ r;
  rnz= rz/ r;
  thx= rnz* phy;
  thy=- rnz* phx;
  thz=- rho/ r;
  *eth= cix* thx+ ciy* thy+ ciz* thz;
  *epi= cix* phx+ ciy* phy;
  *erd= cix* rnx+ ciy* rny+ ciz* rnz;

  return;
}

/*-----------------------------------------------------------------------*/

/* integrand for h field of a wire */
void gh( double zk, double *hr, double *hi)
{
  double rs, r, ckr, skr, rr2, rr3;

  rs= zk- zpka;
  rs= rhks+ rs* rs;
  r= sqrt( rs);
  ckr= cos( r);
  skr= sin( r);
  rr2=1./ rs;
  rr3= rr2/ r;
  *hr= skr* rr2+ ckr* rr3;
  *hi= ckr* rr2- skr* rr3;

  return;
}

/*-----------------------------------------------------------------------*/

/* gwave computes the electric field, including ground wave, of a */
/* current element over a ground plane using formulas of k.a. norton */
/* (proc. ire, sept., 1937, pp.1203,1236.) */

void gwave( complex double *erv, complex double *ezv,
    complex double *erh, complex double *ezh, complex double *eph )
{
  double sppp, sppp2, cppp2, cppp, spp, spp2, cpp2, cpp;
  complex double rk1, rk2, t1, t2, t3, t4, p1, rv;
  complex double omr, w, f, q1, rh, v, g, xr1, xr2;
  complex double x1, x2, x3, x4, x5, x6, x7;

  sppp= zmh/ r1;
  sppp2= sppp* sppp;
  cppp2=1.- sppp2;

  if( cppp2 < 1.0e-20)
    cppp2=1.0e-20;

  cppp= sqrt( cppp2);
  spp= zph/ r2;
  spp2= spp* spp;
  cpp2=1.- spp2;

  if( cpp2 < 1.0e-20)
    cpp2=1.0e-20;

  cpp= sqrt( cpp2);
  rk1=- TPJ* r1;
  rk2=- TPJ* r2;
  t1=1. -u2* cpp2;
  t2= csqrt( t1);
  t3=(1. -1./ rk1)/ rk1;
  t4=(1. -1./ rk2)/ rk2;
  p1= rk2* u2* t1/(2.* cpp2);
  rv=( spp- u* t2)/( spp+ u* t2);
  omr=1.- rv;
  w=1./ omr;
  w=(4.0 + 0.0fj)* p1* w* w;
  f= fbar( w);
  q1= rk2* t1/(2.* u2* cpp2);
  rh=( t2- u* spp)/( t2+ u* spp);
  v=1./(1.+ rh);
  v=(4.0 + 0.0fj)* q1* v* v;
  g= fbar( v);
  xr1= xx1/ r1;
  xr2= xx2/ r2;
  x1= cppp2* xr1;
  x2= rv* cpp2* xr2;
  x3= omr* cpp2* f* xr2;
  x4= u* t2* spp*2.* xr2/ rk2;
  x5= xr1* t3*(1.-3.* sppp2);
  x6= xr2* t4*(1.-3.* spp2);
  *ezv=( x1+ x2+ x3- x4- x5- x6)* (-CONST4);
  x1= sppp* cppp* xr1;
  x2= rv* spp* cpp* xr2;
  x3= cpp* omr* u* t2* f* xr2;
  x4= spp* cpp* omr* xr2/ rk2;
  x5=3.* sppp* cppp* t3* xr1;
  x6= cpp* u* t2* omr* xr2/ rk2*.5;
  x7=3.* spp* cpp* t4* xr2;
  *erv=-( x1+ x2- x3+ x4- x5+ x6- x7)* (-CONST4);
  *ezh=-( x1- x2+ x3- x4- x5- x6+ x7)* (-CONST4);
  x1= sppp2* xr1;
  x2= rv* spp2* xr2;
  x4= u2* t1* omr* f* xr2;
  x5= t3*(1.-3.* cppp2)* xr1;
  x6= t4*(1.-3.* cpp2)*(1.- u2*(1.+ rv)- u2* omr* f)* xr2;
  x7= u2* cpp2* omr*(1.-1./ rk2)*( f*( u2* t1- spp2-1./ rk2)+1./rk2)* xr2;
  *erh=( x1- x2- x4- x5+ x6+ x7)* (-CONST4);
  x1= xr1;
  x2= rh* xr2;
  x3=( rh+1.)* g* xr2;
  x4= t3* xr1;
  x5= t4*(1.- u2*(1.+ rv)- u2* omr* f)* xr2;
  x6=.5* u2* omr*( f*( u2* t1- spp2-1./ rk2)+1./ rk2)* xr2/ rk2;
  *eph=-( x1- x2+ x3- x4+ x5+ x6)* (-CONST4);

  return;
}

/*-----------------------------------------------------------------------*/

/* segment end contributions for thin wire approx. */
void gx( double zz, double rh, double xk,
    complex double *gz, complex double *gzp)
{
  double r, r2, rkz;

  r2= zz* zz+ rh* rh;
  r= sqrt( r2);
  rkz= xk* r;
  *gz= cmplx( cos( rkz),- sin( rkz))/ r;
  *gzp=- cmplx(1.0, rkz)* *gz/ r2;

  return;
}

/*-----------------------------------------------------------------------*/

/* segment end contributions for ext. thin wire approx. */
void gxx( double zz, double rh, double a, double a2, double xk, int ira,
    complex double *g1, complex double *g1p, complex double *g2,
    complex double *g2p, complex double *g3, complex double *gzp )
{
  double r, r2, r4, rk, rk2, rh2, t1, t2;
  complex double  gz, c1, c2, c3;

  r2= zz* zz+ rh* rh;
  r= sqrt( r2);
  r4= r2* r2;
  rk= xk* r;
  rk2= rk* rk;
  rh2= rh* rh;
  t1=.25* a2* rh2/ r4;
  t2=.5* a2/ r2;
  c1= cmplx(1.0, rk);
  c2=3.* c1- rk2;
  c3= cmplx(6.0, rk)* rk2-15.* c1;
  gz= cmplx( cos( rk),- sin( rk))/ r;
  *g2= gz*(1.+ t1* c2);
  *g1= *g2- t2* c1* gz;
  gz= gz/ r2;
  *g2p= gz*( t1* c3- c1);
  *gzp= t2* c2* gz;
  *g3= *g2p+ *gzp;
  *g1p= *g3* zz;

  if( ira != 1)
  {
    *g3=( *g3+ *gzp)* rh;
    *gzp=- zz* c1* gz;

    if( rh <= 1.0e-10)
    {
      *g2=0.;
      *g2p=0.;
      return;
    }

    *g2= *g2/ rh;
    *g2p= *g2p* zz/ rh;
    return;

  } /* if( ira != 1) */

  t2=.5* a;
  *g2=- t2* c1* gz;
  *g2p= t2* gz* c2/ r2;
  *g3= rh2* *g2p- a* gz* c1;
  *g2p= *g2p* zz;
  *gzp=- zz* c1* gz;

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine helix generates segment geometry */
/* data for a helix of ns segments */
void helix( double s, double hl, double a1, double b1,
    double a2, double b2, double rad, int ns, int itg )
{
  int ist, i, mreq;
  double turns, zinc, copy, sangle, hdia, turn, pitch, hmaj, hmin;

  ist= n;
  n += ns;
  np= n;
  mp= m;
  ipsym=0;

  if( ns < 1)
    return;

  turns= fabs( hl/ s);
  zinc= fabs( hl/ ns);

  /* Reallocate tags buffer */
  mem_realloc( (void *)&itag, (n+m) * sizeof(int) );/*????*/

  /* Reallocate wire buffers */
  mreq = n * sizeof(double);
  mem_realloc( (void *)&x, mreq );
  mem_realloc( (void *)&y, mreq );
  mem_realloc( (void *)&z, mreq );
  mem_realloc( (void *)&x2, mreq );
  mem_realloc( (void *)&y2, mreq );
  mem_realloc( (void *)&z2, mreq );
  mem_realloc( (void *)&bi, mreq );

  z[ist]=0.;
  for( i = ist; i < n; i++ )
  {
    bi[i]= rad;
    itag[i]= itg;

    if( i != ist )
      z[i]= z[i-1]+ zinc;

    z2[i]= z[i]+ zinc;

    if( a2 == a1)
    {
      if( b1 == 0.)
	b1= a1;

      x[i]= a1* cos(2.* PI* z[i]/ s);
      y[i]= b1* sin(2.* PI* z[i]/ s);
      x2[i]= a1* cos(2.* PI* z2[i]/ s);
      y2[i]= b1* sin(2.* PI* z2[i]/ s);
    }
    else
    {
      if( b2 == 0.)
	b2= a2;

      x[i]=( a1+( a2- a1)* z[i]/ fabs( hl))* cos(2.* PI* z[i]/ s);
      y[i]=( b1+( b2- b1)* z[i]/ fabs( hl))* sin(2.* PI* z[i]/ s);
      x2[i]=( a1+( a2- a1)* z2[i]/ fabs( hl))* cos(2.* PI* z2[i]/ s);
      y2[i]=( b1+( b2- b1)* z2[i]/ fabs( hl))* sin(2.* PI* z2[i]/ s);

    } /* if( a2 == a1) */

    if( hl > 0.)
      continue;

    copy= x[i];
    x[i]= y[i];
    y[i]= copy;
    copy= x2[i];
    x2[i]= y2[i];
    y2[i]= copy;

  } /* for( i = ist; i < n; i++ ) */

  if( a2 != a1)
  {
    sangle= atan( a2/( fabs( hl)+( fabs( hl)* a1)/( a2- a1)));
    fprintf( output_fp,
	"\n       THE CONE ANGLE OF THE SPIRAL IS %10.4f", sangle );
    return;
  }

  if( a1 == b1)
  {
    hdia=2.* a1;
    turn= hdia* PI;
    pitch= atan( s/( PI* hdia));
    turn= turn/ cos( pitch);
    pitch=180.* pitch/ PI;
  }
  else
  {
    if( a1 >= b1)
    {
      hmaj=2.* a1;
      hmin=2.* b1;
    }
    else
    {
      hmaj=2.* b1;
      hmin=2.* a1;
    }

    hdia= sqrt(( hmaj*hmaj+ hmin*hmin)/2* hmaj);
    turn=2.* PI* hdia;
    pitch=(180./ PI)* atan( s/( PI* hdia));

  } /* if( a1 == b1) */

  fprintf( output_fp, "\n"
      "       THE PITCH ANGLE IS: %.4f    THE LENGTH OF WIRE/TURN IS: %.4f",
      pitch, turn );

  return;
}

/*-----------------------------------------------------------------------*/

/* hfk computes the h field of a uniform current */
/* filament by numerical integration */
void hfk( double el1, double el2, double rhk,
    double zpkx, double *sgr, double *sgi )
{
  int nx = 1, nma = 65536, nts = 4;
  int ns, nt;
  int flag = TRUE;
  double rx = 1.0e-4;
  double z, ze, s, ep, zend, dz=0., zp, dzot=0., t00r, g1r, g5r, t00i;
  double g1i, g5i, t01r, g3r, t01i, g3i, t10r, t10i, te1i, te1r, t02r;
  double g2r, g4r, t02i, g2i, g4i, t11r, t11i, t20r, t20i, te2i, te2r;

  zpka= zpkx;
  rhks= rhk* rhk;
  z= el1;
  ze= el2;
  s= ze- z;
  ep= s/(10.* nma);
  zend= ze- ep;
  *sgr=0.0;
  *sgi=0.0;
  ns= nx;
  nt=0;
  gh( z, &g1r, &g1i);

  while( TRUE )
  {
    if( flag )
    {
      dz= s/ ns;
      zp= z+ dz;

      if( zp > ze )
      {
	dz= ze- z;
	if( fabs(dz) <= ep )
	{
	  *sgr= *sgr* rhk*.5;
	  *sgi= *sgi* rhk*.5;
	  return;
	}
      }

      dzot= dz*.5;
      zp= z+ dzot;
      gh( zp, &g3r, &g3i);
      zp= z+ dz;
      gh( zp, &g5r, &g5i);

    } /* if( flag ) */

    t00r=( g1r+ g5r)* dzot;
    t00i=( g1i+ g5i)* dzot;
    t01r=( t00r+ dz* g3r)*0.5;
    t01i=( t00i+ dz* g3i)*0.5;
    t10r=(4.0* t01r- t00r)/3.0;
    t10i=(4.0* t01i- t00i)/3.0;

    test( t01r, t10r, &te1r, t01i, t10i, &te1i, 0.);
    if( (te1i <= rx) && (te1r <= rx) )
    {
      *sgr= *sgr+ t10r;
      *sgi= *sgi+ t10i;
      nt += 2;

      z += dz;
      if( z >= zend)
      {
	*sgr= *sgr* rhk*.5;
	*sgi= *sgi* rhk*.5;
	return;
      }

      g1r= g5r;
      g1i= g5i;
      if( nt >= nts)
	if( ns > nx)
	{
	  ns= ns/2;
	  nt=1;
	}
      flag = TRUE;
      continue;

    } /* if( (te1i <= rx) && (te1r <= rx) ) */

    zp= z+ dz*0.25;
    gh( zp, &g2r, &g2i);
    zp= z+ dz*0.75;
    gh( zp, &g4r, &g4i);
    t02r=( t01r+ dzot*( g2r+ g4r))*0.5;
    t02i=( t01i+ dzot*( g2i+ g4i))*0.5;
    t11r=(4.0* t02r- t01r)/3.0;
    t11i=(4.0* t02i- t01i)/3.0;
    t20r=(16.0* t11r- t10r)/15.0;
    t20i=(16.0* t11i- t10i)/15.0;

    test( t11r, t20r, &te2r, t11i, t20i, &te2i, 0.);
    if( (te2i > rx) || (te2r > rx) )
    {
      nt=0;
      if( ns >= nma)
	fprintf( output_fp, "\n  STEP SIZE LIMITED AT Z= %10.5f", z );
      else
      {
	ns= ns*2;
	dz= s/ ns;
	dzot= dz*0.5;
	g5r= g3r;
	g5i= g3i;
	g3r= g2r;
	g3i= g2i;

	flag = FALSE;
	continue;
      }

    } /* if( (te2i > rx) || (te2r > rx) ) */

    *sgr= *sgr+ t20r;
    *sgi= *sgi+ t20i;
    nt++;

    z += dz;
    if( z >= zend)
    {
      *sgr= *sgr* rhk*.5;
      *sgi= *sgi* rhk*.5;
      return;
    }

    g1r= g5r;
    g1i= g5i;
    if( nt >= nts)
      if( ns > nx)
      {
	ns= ns/2;
	nt=1;
      }
    flag = TRUE;

  } /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* hintg computes the h field of a patch current */
void hintg( double xi, double yi, double zi )
{
  int ip;
  double rx, ry, rfl, xymag, pxx, pyy, cth;
  double rz, rsq, r, rk, cr, sr, t1zr, t2zr;
  complex double  gam, f1x, f1y, f1z, f2x, f2y, f2z, rrv, rrh;

  rx= xi- xj;
  ry= yi- yj;
  rfl=-1.;
  exk=CPLX_00;
  eyk=CPLX_00;
  ezk=CPLX_00;
  exs=CPLX_00;
  eys=CPLX_00;
  ezs=CPLX_00;

  for( ip = 1; ip <= ksymp; ip++ )
  {
    rfl=- rfl;
    rz= zi- zj* rfl;
    rsq= rx* rx+ ry* ry+ rz* rz;

    if( rsq < 1.0e-20)
      continue;

    r = sqrt( rsq );
    rk= TP* r;
    cr= cos( rk);
    sr= sin( rk);
    gam=-( cmplx(cr,-sr)+rk*cmplx(sr,cr) )/( FPI*rsq*r )* s;
    exc= gam* rx;
    eyc= gam* ry;
    ezc= gam* rz;
    t1zr= t1zj* rfl;
    t2zr= t2zj* rfl;
    f1x= eyc* t1zr- ezc* t1yj;
    f1y= ezc* t1xj- exc* t1zr;
    f1z= exc* t1yj- eyc* t1xj;
    f2x= eyc* t2zr- ezc* t2yj;
    f2y= ezc* t2xj- exc* t2zr;
    f2z= exc* t2yj- eyc* t2xj;

    if( ip != 1)
    {
      if( iperf == 1)
      {
	f1x=- f1x;
	f1y=- f1y;
	f1z=- f1z;
	f2x=- f2x;
	f2y=- f2y;
	f2z=- f2z;
      }
      else
      {
	xymag= sqrt( rx* rx+ ry* ry);
	if( xymag <= 1.0e-6)
	{
	  pxx=0.;
	  pyy=0.;
	  cth=1.;
	  rrv=CPLX_10;
	}
	else
	{
	  pxx=- ry/ xymag;
	  pyy= rx/ xymag;
	  cth= rz/ r;
	  rrv= csqrt(1.- zrati* zrati*(1.- cth* cth));

	} /* if( xymag <= 1.0e-6) */

	rrh= zrati* cth;
	rrh=( rrh- rrv)/( rrh+ rrv);
	rrv= zrati* rrv;
	rrv=-( cth- rrv)/( cth+ rrv);
	gam=( f1x* pxx+ f1y* pyy)*( rrv- rrh);
	f1x= f1x* rrh+ gam* pxx;
	f1y= f1y* rrh+ gam* pyy;
	f1z= f1z* rrh;
	gam=( f2x* pxx+ f2y* pyy)*( rrv- rrh);
	f2x= f2x* rrh+ gam* pxx;
	f2y= f2y* rrh+ gam* pyy;
	f2z= f2z* rrh;

      } /* if( iperf == 1) */

    } /* if( ip != 1) */

    exk += f1x;
    eyk += f1y;
    ezk += f1z;
    exs += f2x;
    eys += f2y;
    ezs += f2z;

  } /* for( ip = 1; ip <= ksymp; ip++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* hsfld computes the h field for constant, sine, and */
/* cosine current on a segment including ground effects. */
void hsfld( double xi, double yi, double zi, double ai )
{
  int ip;
  double xij, yij, rfl, salpr, zij, zp, rhox, rhoy, rhoz, rh, phx;
  double phy, phz, rmag, xymag, xspec, yspec, rhospc, px, py, cth;
  complex double hpk, hps, hpc, qx, qy, qz, rrv, rrh, zratx;

  xij= xi- xj;
  yij= yi- yj;
  rfl=-1.;

  for( ip = 0; ip < ksymp; ip++ )
  {
    rfl=- rfl;
    salpr= salpj* rfl;
    zij= zi- rfl* zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;
    rh= sqrt( rhox* rhox+ rhoy* rhoy+ rhoz* rhoz+ ai* ai);

    if( rh <= 1.0e-10)
    {
      exk=0.;
      eyk=0.;
      ezk=0.;
      exs=0.;
      eys=0.;
      ezs=0.;
      exc=0.;
      eyc=0.;
      ezc=0.;
      continue;
    }

    rhox= rhox/ rh;
    rhoy= rhoy/ rh;
    rhoz= rhoz/ rh;
    phx= sabj* rhoz- salpr* rhoy;
    phy= salpr* rhox- cabj* rhoz;
    phz= cabj* rhoy- sabj* rhox;

    hsflx( s, rh, zp, &hpk, &hps, &hpc);

    if( ip == 1 )
    {
      if( iperf != 1 )
      {
	zratx= zrati;
	rmag= sqrt( zp* zp+ rh* rh);
	xymag= sqrt( xij* xij+ yij* yij);

	/* set parameters for radial wire ground screen. */
	if( nradl != 0)
	{
	  xspec=( xi* zj+ zi* xj)/( zi+ zj);
	  yspec=( yi* zj+ zi* yj)/( zi+ zj);
	  rhospc= sqrt( xspec* xspec+ yspec* yspec+ t2* t2);

	  if( rhospc <= scrwl)
	  {
	    rrv= t1* rhospc* log( rhospc/ t2);
	    zratx=( rrv* zrati)/( ETA* zrati+ rrv);
	  }
	}

	/* calculation of reflection coefficients when ground is specified. */
	if( xymag <= 1.0e-6)
	{
	  px=0.;
	  py=0.;
	  cth=1.;
	  rrv=CPLX_10;
	}
	else
	{
	  px=- yij/ xymag;
	  py= xij/ xymag;
	  cth= zij/ rmag;
	  rrv= csqrt(1.- zratx* zratx*(1.- cth* cth));
	}

	rrh= zratx* cth;
	rrh=-( rrh- rrv)/( rrh+ rrv);
	rrv= zratx* rrv;
	rrv=( cth- rrv)/( cth+ rrv);
	qy=( phx* px+ phy* py)*( rrv- rrh);
	qx= qy* px+ phx* rrh;
	qy= qy* py+ phy* rrh;
	qz= phz* rrh;
	exk= exk- hpk* qx;
	eyk= eyk- hpk* qy;
	ezk= ezk- hpk* qz;
	exs= exs- hps* qx;
	eys= eys- hps* qy;
	ezs= ezs- hps* qz;
	exc= exc- hpc* qx;
	eyc= eyc- hpc* qy;
	ezc= ezc- hpc* qz;
	continue;

      } /* if( iperf != 1 ) */

      exk= exk- hpk* phx;
      eyk= eyk- hpk* phy;
      ezk= ezk- hpk* phz;
      exs= exs- hps* phx;
      eys= eys- hps* phy;
      ezs= ezs- hps* phz;
      exc= exc- hpc* phx;
      eyc= eyc- hpc* phy;
      ezc= ezc- hpc* phz;
      continue;

    } /* if( ip == 1 ) */

    exk= hpk* phx;
    eyk= hpk* phy;
    ezk= hpk* phz;
    exs= hps* phx;
    eys= hps* phy;
    ezs= hps* phz;
    exc= hpc* phx;
    eyc= hpc* phy;
    ezc= hpc* phz;

  } /* for( ip = 0; ip < ksymp; ip++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* calculates h field of sine cosine, and constant current of segment */
void hsflx( double s, double rh, double zpx,
    complex double *hpk, complex double *hps,
    complex double *hpc )
{
  double r1, r2, zp, z2a, hss, dh, z1;
  double rhz, dk, cdk, sdk, hkr, hki, rh2;
  complex double fjk, ekr1, ekr2, t1, t2, cons;

  fjk = -TPJ;
  if( rh >= 1.0e-10)
  {
    if( zpx >= 0.)
    {
      zp= zpx;
      hss=1.;
    }
    else
    {
      zp=- zpx;
      hss=-1.;
    }

    dh=.5* s;
    z1= zp+ dh;
    z2a= zp- dh;
    if( z2a >= 1.0e-7)
      rhz= rh/ z2a;
    else
      rhz=1.;

    dk= TP* dh;
    cdk= cos( dk);
    sdk= sin( dk);
    hfk(- dk, dk, rh* TP, zp* TP, &hkr, &hki);
    *hpk= cmplx( hkr, hki);

    if( rhz >= 1.0e-3)
    {
      rh2= rh* rh;
      r1= sqrt( rh2+ z1* z1);
      r2= sqrt( rh2+ z2a* z2a);
      ekr1= cexp( fjk* r1);
      ekr2= cexp( fjk* r2);
      t1= z1* ekr1/ r1;
      t2= z2a* ekr2/ r2;
      *hps=( cdk*( ekr2- ekr1)- CPLX_01* sdk*( t2+ t1))* hss;
      *hpc=- sdk*( ekr2+ ekr1)- CPLX_01* cdk*( t2- t1);
      cons=- CPLX_01/(2.* TP* rh);
      *hps= cons* *hps;
      *hpc= cons* *hpc;
      return;

    } /* if( rhz >= 1.0e-3) */

    ekr1= cmplx( cdk, sdk)/( z2a* z2a);
    ekr2= cmplx( cdk,- sdk)/( z1* z1);
    t1= TP*(1./ z1-1./ z2a);
    t2= cexp( fjk* zp)* rh/ PI8;
    *hps= t2*( t1+( ekr1+ ekr2)* sdk)* hss;
    *hpc= t2*(- CPLX_01* t1+( ekr1- ekr2)* cdk);
    return;

  } /* if( rh >= 1.0e-10) */

  *hps=CPLX_00;
  *hpc=CPLX_00;
  *hpk=CPLX_00;

  return;
}

/*-----------------------------------------------------------------------*/

/* intrp uses bivariate cubic interpolation to obtain */
/* the values of 4 functions at the point (x,y). */
void intrp( double x, double y, complex double *f1,
    complex double *f2, complex double *f3, complex double *f4 )
{
  static int ix, iy, ixs=-10, iys=-10, igrs=-10, ixeg=0, iyeg=0;
  static int nxm2, nym2, nxms, nyms, nd, ndp;
  int nda[3]={11,17,9}, ndpa[3]={110,85,72};
  int igr, iadd, iadz, i, k, jump;
  static double dx = 1., dy = 1., xs = 0., ys = 0., xz, yz;
  double xx, yy;
  complex double a[4][4], b[4][4], c[4][4], d[4][4];
  complex double p1, p2, p3, p4, fx1, fx2, fx3, fx4;

  jump = TRUE;
  if( (x < xs) || (y < ys) )
    jump = FALSE;
  else
  {
    ix= (int)(( x- xs)/ dx)+1;
    iy= (int)(( y- ys)/ dy)+1;
  }

  /* if point lies in same 4 by 4 point region */
  /* as previous point, old values are reused. */
  if( (ix < ixeg) ||
      (iy < iyeg) ||
      (abs(ix- ixs) >= 2) ||
      (abs(iy- iys) >= 2) ||
      (! jump) )
  {
    /* determine correct grid and grid region */
    if( x <= xsa[1])
      igr=0;
    else
    {
      if( y > ysa[2])
	igr=2;
      else
	igr=1;
    }

    if( igr != igrs)
    {
      igrs= igr;
      dx= dxa[igrs];
      dy= dya[igrs];
      xs= xsa[igrs];
      ys= ysa[igrs];
      nxm2= nxa[igrs]-2;
      nym2= nya[igrs]-2;
      nxms=(( nxm2+1)/3)*3+1;
      nyms=(( nym2+1)/3)*3+1;
      nd= nda[igrs];
      ndp= ndpa[igrs];
      ix= (int)(( x- xs)/ dx)+1;
      iy= (int)(( y- ys)/ dy)+1;

    } /* if( igr != igrs) */

    ixs=(( ix-1)/3)*3+2;
    if( ixs < 2)
      ixs=2;
    ixeg=-10000;

    if( ixs > nxm2)
    {
      ixs= nxm2;
      ixeg= nxms;
    }

    iys=(( iy-1)/3)*3+2;
    if( iys < 2)
      iys=2;
    iyeg=-10000;

    if( iys > nym2)
    {
      iys= nym2;
      iyeg= nyms;
    }

    /* compute coefficients of 4 cubic polynomials in x for */
    /* the 4 grid values of y for each of the 4 functions */
    iadz= ixs+( iys-3)* nd- ndp;
    for( k = 0; k < 4; k++ )
    {
      iadz += ndp;
      iadd = iadz;

      for( i = 0; i < 4; i++ )
      {
	iadd += nd;

	switch( igrs )
	{
	  case 0:
	    p1= ar1[iadd-2];
	    p2= ar1[iadd-1];
	    p3= ar1[iadd];
	    p4= ar1[iadd+1];
	    break;

	  case 1:
	    p1= ar2[iadd-2];
	    p2= ar2[iadd-1];
	    p3= ar2[iadd];
	    p4= ar2[iadd+1];
	    break;

	  case 2:
	    p1= ar3[iadd-2];
	    p2= ar3[iadd-1];
	    p3= ar3[iadd];
	    p4= ar3[iadd+1];
	} /* switch( igrs ) */

	a[i][k]=( p4- p1+3.*( p2- p3))*.1666666667;
	b[i][k]=( p1-2.* p2+ p3)*.5;
	c[i][k]= p3-(2.* p1+3.* p2+ p4)*.1666666667;
	d[i][k]= p2;

      } /* for( i = 0; i < 4; i++ ) */

    } /* for( k = 0; k < 4; k++ ) */

    xz=( ixs-1)* dx+ xs;
    yz=( iys-1)* dy+ ys;

  } /* if( (abs(ix- ixs) >= 2) || */

  /* evaluate polymomials in x and use cubic */
  /* interpolation in y for each of the 4 functions. */
  xx=( x- xz)/ dx;
  yy=( y- yz)/ dy;
  fx1=(( a[0][0]* xx+ b[0][0])* xx+ c[0][0])* xx+ d[0][0];
  fx2=(( a[1][0]* xx+ b[1][0])* xx+ c[1][0])* xx+ d[1][0];
  fx3=(( a[2][0]* xx+ b[2][0])* xx+ c[2][0])* xx+ d[2][0];
  fx4=(( a[3][0]* xx+ b[3][0])* xx+ c[3][0])* xx+ d[3][0];
  p1= fx4- fx1+3.*( fx2- fx3);
  p2=3.*( fx1-2.* fx2+ fx3);
  p3=6.* fx3-2.* fx1-3.* fx2- fx4;
  *f1=(( p1* yy+ p2)* yy+ p3)* yy*.1666666667+ fx2;
  fx1=(( a[0][1]* xx+ b[0][1])* xx+ c[0][1])* xx+ d[0][1];
  fx2=(( a[1][1]* xx+ b[1][1])* xx+ c[1][1])* xx+ d[1][1];
  fx3=(( a[2][1]* xx+ b[2][1])* xx+ c[2][1])* xx+ d[2][1];
  fx4=(( a[3][1]* xx+ b[3][1])* xx+ c[3][1])* xx+ d[3][1];
  p1= fx4- fx1+3.*( fx2- fx3);
  p2=3.*( fx1-2.* fx2+ fx3);
  p3=6.* fx3-2.* fx1-3.* fx2- fx4;
  *f2=(( p1* yy+ p2)* yy+ p3)* yy*.1666666667+ fx2;
  fx1=(( a[0][2]* xx+ b[0][2])* xx+ c[0][2])* xx+ d[0][2];
  fx2=(( a[1][2]* xx+ b[1][2])* xx+ c[1][2])* xx+ d[1][2];
  fx3=(( a[2][2]* xx+ b[2][2])* xx+ c[2][2])* xx+ d[2][2];
  fx4=(( a[3][2]* xx+ b[3][2])* xx+ c[3][2])* xx+ d[3][2];
  p1= fx4- fx1+3.*( fx2- fx3);
  p2=3.*( fx1-2.* fx2+ fx3);
  p3=6.* fx3-2.* fx1-3.* fx2- fx4;
  *f3=(( p1* yy+ p2)* yy+ p3)* yy*.1666666667+ fx2;
  fx1=(( a[0][3]* xx+ b[0][3])* xx+ c[0][3])* xx+ d[0][3];
  fx2=(( a[1][3]* xx+ b[1][3])* xx+ c[1][3])* xx+ d[1][3];
  fx3=(( a[2][3]* xx+ b[2][3])* xx+ c[2][3])* xx+ d[2][3];
  fx4=(( a[3][3]* xx+ b[3][3])* xx+ c[3][3])* xx+ d[3][3];
  p1= fx4- fx1+3.*( fx2- fx3);
  p2=3.*( fx1-2.* fx2+ fx3);
  p3=6.* fx3-2.* fx1-3.* fx2- fx4;
  *f4=(( p1* yy+ p2)* yy+ p3)* yy*.16666666670+ fx2;

  return;
}

/*-----------------------------------------------------------------------*/

/* intx performs numerical integration of exp(jkr)/r by the method of */
/* variable interval width romberg integration.  the integrand value */
/* is supplied by subroutine gf. */
void intx( double el1, double el2, double b, int ij, double *sgr, double *sgi)
{
  int ns, nt;
  int nx = 1, nma = 65536, nts = 4;
  int flag = TRUE;
  double z, s, ze, fnm, ep, zend, fns, dz=0., zp, dzot=0., t00r, g1r, g5r, t00i;
  double g1i, g5i, t01r, g3r, t01i, g3i, t10r, t10i, te1i, te1r, t02r;
  double g2r, g4r, t02i, g2i, g4i, t11r, t11i, t20r, t20i, te2i, te2r;
  double rx = 1.0e-4;

  z= el1;
  ze= el2;
  if( ij == 0)
    ze=0.;
  s= ze- z;
  fnm= nma;
  ep= s/(10.* fnm);
  zend= ze- ep;
  *sgr=0.;
  *sgi=0.;
  ns= nx;
  nt=0;
  gf( z, &g1r, &g1i);

  while( TRUE )
  {
    if( flag )
    {
      fns= ns;
      dz= s/ fns;
      zp= z+ dz;

      if( zp > ze)
      {
	dz= ze- z;
	if( fabs(dz) <= ep)
	{
	  /* add contribution of near singularity for diagonal term */
	  if(ij == 0)
	  {
	    *sgr=2.*( *sgr+ log(( sqrt( b* b+ s* s)+ s)/ b));
	    *sgi=2.* *sgi;
	  }
	  return;
	}

      } /* if( zp > ze) */

      dzot= dz*.5;
      zp= z+ dzot;
      gf( zp, &g3r, &g3i);
      zp= z+ dz;
      gf( zp, &g5r, &g5i);

    } /* if( flag ) */

    t00r=( g1r+ g5r)* dzot;
    t00i=( g1i+ g5i)* dzot;
    t01r=( t00r+ dz* g3r)*0.5;
    t01i=( t00i+ dz* g3i)*0.5;
    t10r=(4.0* t01r- t00r)/3.0;
    t10i=(4.0* t01i- t00i)/3.0;

    /* test convergence of 3 point romberg result. */
    test( t01r, t10r, &te1r, t01i, t10i, &te1i, 0.);
    if( (te1i <= rx) && (te1r <= rx) )
    {
      *sgr= *sgr+ t10r;
      *sgi= *sgi+ t10i;
      nt += 2;

      z += dz;
      if( z >= zend)
      {
	/* add contribution of near singularity for diagonal term */
	if(ij == 0)
	{
	  *sgr=2.*( *sgr+ log(( sqrt( b* b+ s* s)+ s)/ b));
	  *sgi=2.* *sgi;
	}
	return;
      }

      g1r= g5r;
      g1i= g5i;
      if( nt >= nts)
	if( ns > nx)
	{
	  /* Double step size */
	  ns= ns/2;
	  nt=1;
	}
      flag = TRUE;
      continue;

    } /* if( (te1i <= rx) && (te1r <= rx) ) */

    zp= z+ dz*0.25;
    gf( zp, &g2r, &g2i);
    zp= z+ dz*0.75;
    gf( zp, &g4r, &g4i);
    t02r=( t01r+ dzot*( g2r+ g4r))*0.5;
    t02i=( t01i+ dzot*( g2i+ g4i))*0.5;
    t11r=(4.0* t02r- t01r)/3.0;
    t11i=(4.0* t02i- t01i)/3.0;
    t20r=(16.0* t11r- t10r)/15.0;
    t20i=(16.0* t11i- t10i)/15.0;

    /* test convergence of 5 point romberg result. */
    test( t11r, t20r, &te2r, t11i, t20i, &te2i, 0.);
    if( (te2i > rx) || (te2r > rx) )
    {
      nt=0;
      if( ns >= nma)
	fprintf( output_fp, "\n  STEP SIZE LIMITED AT Z= %10.5f", z );
      else
      {
	/* halve step size */
	ns= ns*2;
	fns= ns;
	dz= s/ fns;
	dzot= dz*0.5;
	g5r= g3r;
	g5i= g3i;
	g3r= g2r;
	g3i= g2i;

	flag = FALSE;
	continue;
      }

    } /* if( (te2i > rx) || (te2r > rx) ) */

    *sgr= *sgr+ t20r;
    *sgi= *sgi+ t20i;
    nt++;

    z += dz;
    if( z >= zend)
    {
      /* add contribution of near singularity for diagonal term */
      if(ij == 0)
      {
	*sgr=2.*( *sgr+ log(( sqrt( b* b+ s* s)+ s)/ b));
	*sgi=2.* *sgi;
      }
      return;
    }

    g1r= g5r;
    g1i= g5i;
    if( nt >= nts)
      if( ns > nx)
      {
	/* Double step size */
	ns= ns/2;
	nt=1;
      }
    flag = TRUE;

  } /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* isegno returns the segment number of the mth segment having the */
/* tag number itagi.  if itagi=0 segment number m is returned. */
int isegno( int itagi, int mx)
{
  int icnt, i, iseg;

  if( mx <= 0)
  {
    fprintf( output_fp,
	"\n  CHECK DATA, PARAMETER SPECIFYING SEGMENT"
	" POSITION IN A GROUP OF EQUAL TAGS MUST NOT BE ZERO" );
    stop(-1);
  }

  icnt=0;
  if( itagi == 0)
  {
    iseg = mx;
    return( iseg );
  }

  if( n > 0)
  {
    for( i = 0; i < n; i++ )
    {
      if( itag[i] != itagi )
	continue;

      icnt++;
      if( icnt == mx)
      {
	iseg= i+1;
	return( iseg );
      }

    } /* for( i = 0; i < n; i++ ) */

  } /* if( n > 0) */

  fprintf( output_fp, "\n\n"
      "  NO SEGMENT HAS AN ITAG OF %d",  itagi );
  stop(-1);

  return(0);
}

/*-----------------------------------------------------------------------*/

/* load calculates the impedance of specified */
/* segments for various types of loading */
void load( int *ldtyp, int *ldtag, int *ldtagf, int *ldtagt,
    double *zlr, double *zli, double *zlc )
{
  int i, iwarn, istep, istepx, l1, l2, ldtags, jump, ichk;
  complex double zt, tpcj;

  tpcj = (0.0+1.883698955e+9fj);
  fprintf( output_fp, "\n"
      "  LOCATION        RESISTANCE  INDUCTANCE  CAPACITANCE   "
      "  IMPEDANCE (OHMS)   CONDUCTIVITY  CIRCUIT\n"
      "  ITAG FROM THRU     OHMS       HENRYS      FARADS     "
      "  REAL     IMAGINARY   MHOS/METER      TYPE" );

  /* initialize d array, used for temporary */
  /* storage of loading information. */
  mem_alloc( (void *)&zarray, npm * sizeof(complex double) );
  for( i = 0; i < n; i++ )
    zarray[i]=CPLX_00;

  iwarn=FALSE;
  istep=0;

  /* cycle over loading cards */
  while( TRUE )
  {
    istepx = istep;
    istep++;

    if( istep > nload)
    {
      if( iwarn == TRUE )
	fprintf( output_fp,
	    "\n  NOTE, SOME OF THE ABOVE SEGMENTS "
	    "HAVE BEEN LOADED TWICE - IMPEDANCES ADDED" );

      if( nop == 1)
	return;

      for( i = 0; i < np; i++ )
      {
	zt= zarray[i];
	l1= i;

	for( l2 = 1; l2 < nop; l2++ )
	{
	  l1 += np;
	  zarray[l1]= zt;
	}
      }
      return;

    } /* if( istep > nload) */

    if( ldtyp[istepx] > 5 )
    {
      fprintf( output_fp,
	  "\n  IMPROPER LOAD TYPE CHOSEN,"
	  " REQUESTED TYPE IS %d", ldtyp[istepx] );
      stop(-1);
    }

    /* search segments for proper itags */
    ldtags= ldtag[istepx];
    jump= ldtyp[istepx]+1;
    ichk=0;
    l1= 1;
    l2= n;

    if( ldtags == 0)
    {
      if( (ldtagf[istepx] != 0) || (ldtagt[istepx] != 0) )
      {
	l1= ldtagf[istepx];
	l2= ldtagt[istepx];

      } /* if( (ldtagf[istepx] != 0) || (ldtagt[istepx] != 0) ) */

    } /* if( ldtags == 0) */

    for( i = l1-1; i < l2; i++ )
    {
      if( ldtags != 0)
      {
	if( ldtags != itag[i])
	  continue;

	if( ldtagf[istepx] != 0)
	{
	  ichk++;
	  if( (ichk < ldtagf[istepx]) || (ichk > ldtagt[istepx]) )
	    continue;
	}
	else
	  ichk=1;

      } /* if( ldtags != 0) */
      else
	ichk=1;

      /* calculation of lamda*imped. per unit length, */
      /* jump to appropriate section for loading type */
      switch( jump )
      {
	case 1:
	  zt= zlr[istepx]/ si[i]+ tpcj* zli[istepx]/( si[i]* wlam);
	  if( fabs( zlc[istepx]) > 1.0e-20)
	    zt += wlam/( tpcj* si[i]* zlc[istepx]);
	  break;

	case 2:
	  zt= tpcj* si[i]* zlc[istepx]/ wlam;
	  if( fabs( zli[istepx]) > 1.0e-20)
	    zt += si[i]* wlam/( tpcj* zli[istepx]);
	  if( fabs( zlr[istepx]) > 1.0e-20)
	    zt += si[i]/ zlr[istepx];
	  zt=1./ zt;
	  break;

	case 3:
	  zt= zlr[istepx]* wlam+ tpcj* zli[istepx];
	  if( fabs( zlc[istepx]) > 1.0e-20)
	    zt += 1./( tpcj* si[i]* si[i]* zlc[istepx]);
	  break;

	case 4:
	  zt= tpcj* si[i]* si[i]* zlc[istepx];
	  if( fabs( zli[istepx]) > 1.0e-20)
	    zt += 1./( tpcj* zli[istepx]);
	  if( fabs( zlr[istepx]) > 1.0e-20)
	    zt += 1./( zlr[istepx]* wlam);
	  zt=1./ zt;
	  break;

	case 5:
	  zt= cmplx( zlr[istepx], zli[istepx])/ si[i];
	  break;

	case 6:
	  zt= zint( zlr[istepx]* wlam, bi[i]);

      } /* switch( jump ) */

      if(( fabs( creal( zarray[i]))+ fabs( cimag( zarray[i]))) > 1.0e-20)
	iwarn=TRUE;
      zarray[i] += zt;

    } /* for( i = l1-1; i < l2; i++ ) */

    if( ichk == 0 )
    {
      fprintf( output_fp,
	  "\n  LOADING DATA CARD ERROR,"
	  " NO SEGMENT HAS AN ITAG = %d", ldtags );
      stop(-1);
    }

    /* printing the segment loading data, jump to proper print */
    switch( jump )
    {
      case 1:
	prnt( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
	    zli[istepx], zlc[istepx],0.,0.,0.," SERIES ", 2);
	break;

      case 2:
	prnt( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
	    zli[istepx], zlc[istepx],0.,0.,0.,"PARALLEL",2);
	break;

      case 3:
	prnt( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
	    zli[istepx], zlc[istepx],0.,0.,0., "SERIES (PER METER)", 5);
	break;

      case 4:
	prnt( ldtags, ldtagf[istepx], ldtagt[istepx], zlr[istepx],
	    zli[istepx], zlc[istepx],0.,0.,0.,"PARALLEL (PER METER)",5);
	break;

      case 5:
	prnt( ldtags, ldtagf[istepx], ldtagt[istepx],0.,0.,0.,
	    zlr[istepx], zli[istepx],0.,"FIXED IMPEDANCE ",4);
	break;

      case 6:
	prnt( ldtags, ldtagf[istepx], ldtagt[istepx],
	    0.,0.,0.,0.,0., zlr[istepx],"  WIRE  ",2);

    } /* switch( jump ) */

  } /* while( TRUE ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine move moves the structure with respect to its */
/* coordinate system or reproduces structure in new positions. */
/* structure is rotated about x,y,z axes by rox,roy,roz */
/* respectively, then shifted by xs,ys,zs */
void move( double rox, double roy, double roz, double xs,
    double ys, double zs, int its, int nrpt, int itgi )
{
  int nrp, ix, i1, k, ir, i, ii, mreq;
  double sps, cps, sth, cth, sph, cph, xx, xy;
  double xz, yx, yy, yz, zx, zy, zz, xi, yi, zi;

  if( fabs( rox)+ fabs( roy) > 1.0e-10)
    ipsym= ipsym*3;

  sps= sin( rox);
  cps= cos( rox);
  sth= sin( roy);
  cth= cos( roy);
  sph= sin( roz);
  cph= cos( roz);
  xx= cph* cth;
  xy= cph* sth* sps- sph* cps;
  xz= cph* sth* cps+ sph* sps;
  yx= sph* cth;
  yy= sph* sth* sps+ cph* cps;
  yz= sph* sth* cps- cph* sps;
  zx=- sth;
  zy= cth* sps;
  zz= cth* cps;

  if( nrpt == 0)
    nrp=1;
  else
    nrp= nrpt;

  ix=1;
  if( n > 0)
  {
    i1= isegno( its, 1);
    if( i1 < 1)
      i1= 1;

    ix= i1;
    if( nrpt == 0)
      k= i1-1;
    else
    {
      k= n;
      /* Reallocate tags buffer */
      mreq = n+m + (n+1-i1)*nrpt;
      mem_realloc( (void *)&itag, mreq * sizeof(int) );

      /* Reallocate wire buffers */
      mreq = (n+(n+1-i1)*nrpt) * sizeof(double);
      mem_realloc( (void *)&x, mreq );
      mem_realloc( (void *)&y, mreq );
      mem_realloc( (void *)&z, mreq );
      mem_realloc( (void *)&x2, mreq );
      mem_realloc( (void *)&y2, mreq );
      mem_realloc( (void *)&z2, mreq );
      mem_realloc( (void *)&bi, mreq );
    }

    for( ir = 0; ir < nrp; ir++ )
    {
      for( i = i1-1; i < n; i++ )
      {
	xi= x[i];
	yi= y[i];
	zi= z[i];
	x[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	y[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	z[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	xi= x2[i];
	yi= y2[i];
	zi= z2[i];
	x2[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	y2[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	z2[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	bi[k]= bi[i];
	itag[k]= itag[i];
	if( itag[i] != 0)
	  itag[k]= itag[i]+ itgi;

	k++;

      } /* for( i = i1; i < n; i++ ) */

      i1= n+1;
      n= k;

    } /* for( ir = 0; ir < nrp; ir++ ) */

  } /* if( n >= n2) */

  if( m > 0)
  {
    i1 = 0;
    if( nrpt == 0)
      k= 0;
    else
      k = m;

    /* Reallocate patch buffers */
    mreq = m * (1+nrpt) * sizeof(double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

    for( ii = 0; ii < nrp; ii++ )
    {
      for( i = i1; i < m; i++ )
      {
	xi= px[i];
	yi= py[i];
	zi= pz[i];
	px[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	py[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	pz[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	xi= t1x[i];
	yi= t1y[i];
	zi= t1z[i];
	t1x[k]= xi* xx+ yi* xy+ zi* xz;
	t1y[k]= xi* yx+ yi* yy+ zi* yz;
	t1z[k]= xi* zx+ yi* zy+ zi* zz;
	xi= t2x[i];
	yi= t2y[i];
	zi= t2z[i];
	t2x[k]= xi* xx+ yi* xy+ zi* xz;
	t2y[k]= xi* yx+ yi* yy+ zi* yz;
	t2z[k]= xi* zx+ yi* zy+ zi* zz;
	psalp[k]= psalp[i];
	pbi[k]= pbi[i];
	k++;

      } /* for( i = i1; i < m; i++ ) */

      i1= m;
      m = k;

    } /* for( ii = 0; ii < nrp; ii++ ) */

  } /* if( m >= m2) */

  if( (nrpt == 0) && (ix == 1) )
    return;

  np= n;
  mp= m;
  ipsym=0;

  return;
}

/*-----------------------------------------------------------------------*/

/* nefld computes the near field at specified points in space after */
/* the structure currents have been computed. */
void nefld( double xob, double yob, double zob,
    complex double *ex, complex double *ey, complex double *ez )
{
  int i, ix, ipr, iprx, jc, ipa;
  double zp, xi, ax;
  complex double acx, bcx, ccx;

  *ex=CPLX_00;
  *ey=CPLX_00;
  *ez=CPLX_00;
  ax=0.;

  if( n != 0)
  {
    for( i = 0; i < n; i++ )
    {
      xj= xob- x[i];
      yj= yob- y[i];
      zj= zob- z[i];
      zp= cab[i]* xj+ sab[i]* yj+ salp[i]* zj;

      if( fabs( zp) > 0.5001* si[i])
	continue;

      zp= xj* xj+ yj* yj+ zj* zj- zp* zp;
      xj= bi[i];

      if( zp > 0.9* xj* xj)
	continue;

      ax= xj;
      break;

    } /* for( i = 0; i < n; i++ ) */

    for( i = 0; i < n; i++ )
    {
      ix = i+1;
      s= si[i];
      b= bi[i];
      xj= x[i];
      yj= y[i];
      zj= z[i];
      cabj= cab[i];
      sabj= sab[i];
      salpj= salp[i];

      if( iexk != 0)
      {
	ipr= icon1[i];

	if( ipr < 0 )
	{
	  ipr = -ipr;
	  iprx = ipr-1;

	  if( -icon1[iprx] != ix )
	    ind1=2;
	  else
	  {
	    xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	    if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.0e-6) )
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* if( ipr < 0 ) */
	else
	  if( ipr == 0 )
	    ind1=1;
	  else
	  {
	    iprx = ipr-1;

	    if( ipr != ix )
	    {
	      if( icon2[iprx] != ix )
		ind1=2;
	      else
	      {
		xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
		if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.0e-6) )
		  ind1=2;
		else
		  ind1=0;
	      }
	    } /* if( ipr != ix ) */
	    else
	    {
	      if( cabj* cabj+ sabj* sabj > 1.0e-8)
		ind1=2;
	      else
		ind1=0;
	    }
	  } /* else */

	ipr= icon2[i];

	if( ipr < 0 )
	{
	  ipr = -ipr;
	  iprx = ipr-1;

	  if( -icon2[iprx] != ix )
	    ind1=2;
	  else
	  {
	    xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	    if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.0e-6) )
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* if( ipr < 0 ) */
	else
	  if( ipr == 0 )
	    ind2=1;
	  else
	  {
	    iprx = ipr-1;

	    if( ipr != ix )
	    {
	      if( icon1[iprx] != ix )
		ind2=2;
	      else
	      {
		xi= fabs( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
		if( (xi < 0.999999) || (fabs(bi[iprx]/b-1.) > 1.0e-6) )
		  ind2=2;
		else
		  ind2=0;
	      }
	    } /* if( ipr != (i+1) ) */
	    else
	    {
	      if( cabj* cabj+ sabj* sabj > 1.0e-8)
		ind1=2;
	      else
		ind1=0;
	    }

	  } /* else */

      } /* if( iexk != 0) */

      efld( xob, yob, zob, ax,1);
      acx= cmplx( air[i], aii[i]);
      bcx= cmplx( bir[i], bii[i]);
      ccx= cmplx( cir[i], cii[i]);
      *ex += exk* acx+ exs* bcx+ exc* ccx;
      *ey += eyk* acx+ eys* bcx+ eyc* ccx;
      *ez += ezk* acx+ ezs* bcx+ ezc* ccx;

    } /* for( i = 0; i < n; i++ ) */

    if( m == 0)
      return;

  } /* if( n != 0) */

  jc= n-1;
  for( i = 0; i < m; i++ )
  {
    s= pbi[i];
    xj= px[i];
    yj= py[i];
    zj= pz[i];
    t1xj= t1x[i];
    t1yj= t1y[i];
    t1zj= t1z[i];
    t2xj= t2x[i];
    t2yj= t2y[i];
    t2zj= t2z[i];
    jc += 3;
    acx= t1xj* cur[jc-2]+ t1yj* cur[jc-1]+ t1zj* cur[jc];
    bcx= t2xj* cur[jc-2]+ t2yj* cur[jc-1]+ t2zj* cur[jc];

    for( ipa = 0; ipa < ksymp; ipa++ )
    {
      ipgnd= ipa+1;
      unere( xob, yob, zob);
      *ex= *ex+ acx* exk+ bcx* exs;
      *ey= *ey+ acx* eyk+ bcx* eys;
      *ez= *ez+ acx* ezk+ bcx* ezs;
    }

  } /* for( i = 0; i < m; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine netwk solves for structure currents for a given */
/* excitation including the effect of non-radiating networks if */
/* present. */
void netwk( complex double *cm, complex double *cmb,
    complex double *cmc, complex double *cmd, int *ip,
    complex double *einc )
{
  int *ipnt = NULL, *nteqa = NULL, *ntsca = NULL;
  int jump1, jump2, nteq=0, ntsc=0, nseg2, irow2=0, j, ndimn;
  int neqz2, neqt, irow1=0, i, nseg1, isc1=0, isc2=0;
  double asmx, asa, pwr, y11r, y11i, y12r, y12i, y22r, y22i;
  complex double *vsrc = NULL, *rhs = NULL, *cmn = NULL;
  complex double *rhnt = NULL, *rhnx = NULL, ymit, vlt, cux;

  neqz2= neq2;
  if( neqz2 == 0)
    neqz2=1;

  pin=0.;
  pnls=0.;
  neqt= neq+ neq2;
  ndimn = j = (2*nonet + nsant);

  /* Allocate network buffers */
  if( nonet > 0 )
  {
    mem_alloc( (void *)&rhs, np3m * sizeof(complex double) );

    i = j * sizeof(complex double);
    mem_alloc( (void *)&rhnt, i );
    mem_alloc( (void *)&rhnx, i );
    mem_alloc( (void *)&cmn, i * j );

    i = j * sizeof(int);
    mem_alloc( (void *)&ntsca, i );
    mem_alloc( (void *)&nteqa, i );
    mem_alloc( (void *)&ipnt, i );

    mem_alloc( (void *)&vsrc, nsant * sizeof(complex double) );
  }

  if( ntsol == 0)
  {
    /* compute relative matrix asymmetry */
    if( masym != 0)
    {
      irow1=0;
      if( nonet != 0)
      {
	for( i = 0; i < nonet; i++ )
	{
	  nseg1= iseg1[i];
	  for( isc1 = 0; isc1 < 2; isc1++ )
	  {
	    if( irow1 == 0)
	    {
	      ipnt[irow1]= nseg1;
	      nseg1= iseg2[i];
	      irow1++;
	      continue;
	    }

	    for( j = 0; j < irow1; j++ )
	      if( nseg1 == ipnt[j])
		break;

	    if( j == irow1 )
	    {
	      ipnt[irow1]= nseg1;
	      irow1++;
	    }

	    nseg1= iseg2[i];

	  } /* for( isc1 = 0; isc1 < 2; isc1++ ) */

	} /* for( i = 0; i < nonet; i++ ) */

      } /* if( nonet != 0) */

      if( nsant != 0)
      {
	for( i = 0; i < nsant; i++ )
	{
	  nseg1= isant[i];
	  if( irow1 == 0)
	  {
	    ipnt[irow1]= nseg1;
	    irow1++;
	    continue;
	  }

	  for( j = 0; j < irow1; j++ )
	    if( nseg1 == ipnt[j])
	      break;

	  if( j == irow1 )
	  {
	    ipnt[irow1]= nseg1;
	    irow1++;
	  }

	} /* for( i = 0; i < nsant; i++ ) */

      } /* if( nsant != 0) */

      if( irow1 >= 2)
      {
	for( i = 0; i < irow1; i++ )
	{
	  isc1= ipnt[i]-1;
	  asmx= si[isc1];

	  for( j = 0; j < neqt; j++ )
	    rhs[j] = CPLX_00;

	  rhs[isc1] = CPLX_10;
	  solves( cm, ip, rhs, neq, 1, np, n, mp, m);
	  cabc( rhs);

	  for( j = 0; j < irow1; j++ )
	  {
	    isc1= ipnt[j]-1;
	    cmn[j+i*ndimn]= rhs[isc1]/ asmx;
	  }

	} /* for( i = 0; i < irow1; i++ ) */

	asmx=0.;
	asa=0.;

	for( i = 1; i < irow1; i++ )
	{
	  isc1= i;
	  for( j = 0; j < isc1; j++ )
	  {
	    cux= cmn[i+j*ndimn];
	    pwr= cabs(( cux- cmn[j+i*ndimn])/ cux);
	    asa += pwr* pwr;

	    if( pwr < asmx)
	      continue;

	    asmx= pwr;
	    nteq= ipnt[i];
	    ntsc= ipnt[j];

	  } /* for( j = 0; j < isc1; j++ ) */

	} /* for( i = 1; i < irow1; i++ ) */

	asa= sqrt( asa*2./ (double)( irow1*( irow1-1)));
	fprintf( output_fp, "\n\n"
	    "   MAXIMUM RELATIVE ASYMMETRY OF THE DRIVING POINT ADMITTANCE\n"
	    "   MATRIX IS %10.3E FOR SEGMENTS %d AND %d\n"
	    "   RMS RELATIVE ASYMMETRY IS %10.3E",
	    asmx, nteq, ntsc, asa );

      } /* if( irow1 >= 2) */

    } /* if( masym != 0) */

    /* solution of network equations */
    if( nonet != 0)
    {
      for( i = 0; i < ndimn; i++ )
      {
	rhnx[i]=CPLX_00;
	for( j = 0; j < ndimn; j++ )
	  cmn[j+i*ndimn]=CPLX_00;
      }

      nteq=0;
      ntsc=0;

      /* sort network and source data and */
      /* assign equation numbers to segments */
      for( j = 0; j < nonet; j++ )
      {
	nseg1= iseg1[j];
	nseg2= iseg2[j];

	if( ntyp[j] <= 1)
	{
	  y11r= x11r[j];
	  y11i= x11i[j];
	  y12r= x12r[j];
	  y12i= x12i[j];
	  y22r= x22r[j];
	  y22i= x22i[j];
	}
	else
	{
	  y22r= TP* x11i[j]/ wlam;
	  y12r=0.;
	  y12i=1./( x11r[j]* sin( y22r));
	  y11r= x12r[j];
	  y11i=- y12i* cos( y22r);
	  y22r= x22r[j];
	  y22i= y11i+ x22i[j];
	  y11i= y11i+ x12i[j];

	  if( ntyp[j] != 2)
	  {
	    y12r=- y12r;
	    y12i=- y12i;
	  }

	} /* if( ntyp[j] <= 1) */

	jump1 = FALSE;
	if( nsant != 0)
	{
	  for( i = 0; i < nsant; i++ )
	    if( nseg1 == isant[i])
	    {
	      isc1 = i;
	      jump1 = TRUE;
	      break;
	    }

	} /* if( nsant != 0) */

	jump2 = FALSE;
	if( ! jump1 )
	{
	  isc1=-1;

	  if( nteq != 0)
	  {
	    for( i = 0; i < nteq; i++ )
	      if( nseg1 == nteqa[i])
	      {
		irow1 = i;
		jump2 = TRUE;
		break;
	      }

	  } /* if( nteq != 0) */

	  if( ! jump2 )
	  {
	    irow1= nteq;
	    nteqa[nteq]= nseg1;
	    nteq++;
	  }

	} /* if( ! jump1 ) */
	else
	{
	  if( ntsc != 0)
	  {
	    for( i = 0; i < ntsc; i++ )
	    {
	      if( nseg1 == ntsca[i])
	      {
		irow1 = ndimn- (i+1);
		jump2 = TRUE;
		break;
	      }
	    }

	  } /* if( ntsc != 0) */

	  if( ! jump2 )
	  {
	    irow1= ndimn- (ntsc+1);
	    ntsca[ntsc]= nseg1;
	    vsrc[ntsc]= vsant[isc1];
	    ntsc++;
	  }

	} /* if( ! jump1 ) */

	jump1 = FALSE;
	if( nsant != 0)
	{
	  for( i = 0; i < nsant; i++ )
	  {
	    if( nseg2 == isant[i])
	    {
	      isc2= i;
	      jump1 = TRUE;
	      break;
	    }
	  }

	} /* if( nsant != 0) */

	jump2 = FALSE;
	if( ! jump1 )
	{
	  isc2=-1;

	  if( nteq != 0)
	  {
	    for( i = 0; i < nteq; i++ )
	      if( nseg2 == nteqa[i])
	      {
		irow2= i;
		jump2 = TRUE;
		break;
	      }

	  } /* if( nteq != 0) */

	  if( ! jump2 )
	  {
	    irow2= nteq;
	    nteqa[nteq]= nseg2;
	    nteq++;
	  }

	}  /* if( ! jump1 ) */
	else
	{
	  if( ntsc != 0)
	  {
	    for( i = 0; i < ntsc; i++ )
	      if( nseg2 == ntsca[i])
	      {
		irow2 = ndimn- (i+1);
		jump2 = TRUE;
		break;
	      }

	  } /* if( ntsc != 0) */

	  if( ! jump2 )
	  {
	    irow2= ndimn- (ntsc+1);
	    ntsca[ntsc]= nseg2;
	    vsrc[ntsc]= vsant[isc2];
	    ntsc++;
	  }

	} /* if( ! jump1 ) */

	/* fill network equation matrix and right hand side vector with */
	/* network short-circuit admittance matrix coefficients. */
	if( isc1 == -1)
	{
	  cmn[irow1+irow1*ndimn] -= cmplx( y11r, y11i)* si[nseg1-1];
	  cmn[irow1+irow2*ndimn] -= cmplx( y12r, y12i)* si[nseg1-1];
	}
	else
	{
	  rhnx[irow1] += cmplx( y11r, y11i)* vsant[isc1]/wlam;
	  rhnx[irow2] += cmplx( y12r, y12i)* vsant[isc1]/wlam;
	}

	if( isc2 == -1)
	{
	  cmn[irow2+irow2*ndimn] -= cmplx( y22r, y22i)* si[nseg2-1];
	  cmn[irow2+irow1*ndimn] -= cmplx( y12r, y12i)* si[nseg2-1];
	}
	else
	{
	  rhnx[irow1] += cmplx( y12r, y12i)* vsant[isc2]/wlam;
	  rhnx[irow2] += cmplx( y22r, y22i)* vsant[isc2]/wlam;
	}

      } /* for( j = 0; j < nonet; j++ ) */

      /* add interaction matrix admittance */
      /* elements to network equation matrix */
      for( i = 0; i < nteq; i++ )
      {
	for( j = 0; j < neqt; j++ )
	  rhs[j] = CPLX_00;

	irow1= nteqa[i]-1;
	rhs[irow1]=CPLX_10;
	solves( cm, ip, rhs, neq, 1, np, n, mp, m);
	cabc( rhs);

	for( j = 0; j < nteq; j++ )
	{
	  irow1= nteqa[j]-1;
	  cmn[i+j*ndimn] += rhs[irow1];
	}

      } /* for( i = 0; i < nteq; i++ ) */

      /* factor network equation matrix */
      factr( nteq, cmn, ipnt, ndimn);

    } /* if( nonet != 0) */

  } /* if( ntsol != 0) */

  if( nonet != 0)
  {
    /* add to network equation right hand side */
    /* the terms due to element interactions */
    for( i = 0; i < neqt; i++ )
      rhs[i]= einc[i];

    solves( cm, ip, rhs, neq, 1, np, n, mp, m);
    cabc( rhs);

    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      rhnt[i]= rhnx[i]+ rhs[irow1];
    }

    /* solve network equations */
    solve( nteq, cmn, ipnt, rhnt, ndimn);

    /* add fields due to network voltages to electric fields */
    /* applied to structure and solve for induced current */
    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      einc[irow1] -= rhnt[i];
    }

    solves( cm, ip, einc, neq, 1, np, n, mp, m);
    cabc( einc);

    if( nprint == 0)
    {
      fprintf( output_fp, "\n\n\n"
	  "                          "
	  "--------- STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS --------" );

      fprintf( output_fp, "\n"
	  "  TAG   SEG       VOLTAGE (VOLTS)          CURRENT (AMPS)        "
	  " IMPEDANCE (OHMS)       ADMITTANCE (MHOS)     POWER\n"
	  "  No:   No:     REAL      IMAGINARY     REAL      IMAGINARY    "
	  " REAL      IMAGINARY     REAL      IMAGINARY   (WATTS)" );
    }

    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      vlt= rhnt[i]* si[irow1]* wlam;
      cux= einc[irow1]* wlam;
      ymit= cux/ vlt;
      zped= vlt/ cux;
      irow2= itag[irow1];
      pwr=.5* creal( vlt* conj( cux));
      pnls= pnls- pwr;

      if( nprint == 0)
	fprintf( output_fp, "\n"
	    " %4d %5d %11.4E %11.4E %11.4E %11.4E"
	    " %11.4E %11.4E %11.4E %11.4E %11.4E",
	    irow2, irow1+1, creal(vlt), cimag(vlt), creal(cux), cimag(cux),
	    creal(zped), cimag(zped), creal(ymit), cimag(ymit), pwr );
    }

    if( ntsc != 0)
    {
      for( i = 0; i < ntsc; i++ )
      {
	irow1= ntsca[i]-1;
	vlt= vsrc[i];
	cux= einc[irow1]* wlam;
	ymit= cux/ vlt;
	zped= vlt/ cux;
	irow2= itag[irow1];
	pwr=.5* creal( vlt* conj( cux));
	pnls= pnls- pwr;

	if( nprint == 0)
	  fprintf( output_fp, "\n"
	      " %4d %5d %11.4E %11.4E %11.4E %11.4E"
	      " %11.4E %11.4E %11.4E %11.4E %11.4E",
	      irow2, irow1+1, creal(vlt), cimag(vlt), creal(cux), cimag(cux),
	      creal(zped), cimag(zped), creal(ymit), cimag(ymit), pwr );

      } /* for( i = 0; i < ntsc; i++ ) */

    } /* if( ntsc != 0) */

  } /* if( nonet != 0) */
  else
  {
    /* solve for currents when no networks are present */
    solves( cm, ip, einc, neq, 1, np, n, mp, m);
    cabc( einc);
    ntsc=0;
  }

  if( (nsant+nvqd) == 0)
    return;

  fprintf( output_fp, "\n\n\n"
      "                        "
      "--------- ANTENNA INPUT PARAMETERS ---------" );

  fprintf( output_fp, "\n"
      "  TAG   SEG       VOLTAGE (VOLTS)         "
      "CURRENT (AMPS)         IMPEDANCE (OHMS)    "
      "    ADMITTANCE (MHOS)     POWER\n"
      "  No:   No:     REAL      IMAGINARY"
      "     REAL      IMAGINARY     REAL      "
      "IMAGINARY    REAL       IMAGINARY   (WATTS)" );

  if( nsant != 0)
  {
    for( i = 0; i < nsant; i++ )
    {
      isc1= isant[i]-1;
      vlt= vsant[i];

      if( ntsc == 0)
      {
	cux= einc[isc1]* wlam;
	irow1=0;
      }
      else
      {
	for( j = 0; j < ntsc; j++ )
	  if( ntsca[j] == isc1+1)
	  {
	    irow1= ndimn- (j+1);
	    cux= rhnx[irow1];
	    for( j = 0; j < nteq; j++ )
	      cux -= cmn[j+irow1*ndimn]*rhnt[j];
	    cux=(einc[isc1]+ cux)* wlam;
	    irow1++;
	  }

      } /* if( ntsc == 0) */

      ymit= cux/ vlt;
      zped= vlt/ cux;
      pwr=.5* creal( vlt* conj( cux));
      pin= pin+ pwr;

      if( irow1 != 0)
	pnls= pnls+ pwr;

      irow2= itag[isc1];
      fprintf( output_fp, "\n"
	  " %4d %5d %11.4E %11.4E %11.4E %11.4E"
	  " %11.4E %11.4E %11.4E %11.4E %11.4E",
	  irow2, isc1+1, creal(vlt), cimag(vlt), creal(cux), cimag(cux),
	  creal(zped), cimag(zped), creal(ymit), cimag(ymit), pwr );

    } /* for( i = 0; i < nsant; i++ ) */

  } /* if( nsant != 0) */

  if( nvqd == 0)
    return;

  for( i = 0; i < nvqd; i++ )
  {
    isc1= ivqd[i]-1;
    vlt= vqd[i];
    cux= cmplx( air[isc1], aii[isc1]);
    ymit= cmplx( bir[isc1], bii[isc1]);
    zped= cmplx( cir[isc1], cii[isc1]);
    pwr= si[isc1]* TP*.5;
    cux=( cux- ymit* sin( pwr)+ zped* cos( pwr))* wlam;
    ymit= cux/ vlt;
    zped= vlt/ cux;
    pwr=.5* creal( vlt* conj( cux));
    pin= pin+ pwr;
    irow2= itag[isc1];

    fprintf( output_fp,	"\n"
	" %4d %5d %11.4E %11.4E %11.4E %11.4E"
	" %11.4E %11.4E %11.4E %11.4E %11.4E",
	irow2, isc1+1, creal(vlt), cimag(vlt), creal(cux), cimag(cux),
	creal(zped), cimag(zped), creal(ymit), cimag(ymit), pwr );

  } /* for( i = 0; i < nvqd; i++ ) */

  /* Free network buffers */
  free_ptr( (void *)&ipnt );
  free_ptr( (void *)&nteqa );
  free_ptr( (void *)&ntsca );
  free_ptr( (void *)&vsrc );
  free_ptr( (void *)&rhs );
  free_ptr( (void *)&cmn );
  free_ptr( (void *)&rhnt );
  free_ptr( (void *)&rhnx );

  return;
}

/*-----------------------------------------------------------------------*/

/* compute near e or h fields over a range of points */
void nfpat( void )
{
  int i, j, kk;
  double znrt, cth=0., sth=0., ynrt, cph=0., sph=0., xnrt, xob, yob;
  double zob, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, xxx;
  complex double ex, ey, ez;

  if( nfeh != 1)
  {
    fprintf( output_fp,	"\n\n\n"
	"                             "
	"-------- NEAR ELECTRIC FIELDS --------\n"
	"     ------- LOCATION -------     ------- EX ------    ------- EY ------    ------- EZ ------\n"
	"      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n"
	"    METERS    METERS    METERS     VOLTS/M  DEGREES    VOLTS/M   DEGREES     VOLTS/M  DEGREES" );
  }
  else
  {
    fprintf( output_fp,	"\n\n\n"
	"                                   "
	"-------- NEAR MAGNETIC FIELDS ---------\n\n"
	"     ------- LOCATION -------     ------- HX ------    ------- HY ------    ------- HZ ------\n"
	"      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n"
	"    METERS    METERS    METERS      AMPS/M  DEGREES      AMPS/M  DEGREES      AMPS/M  DEGREES" );
  }

  znrt= znr- dznr;
  for( i = 0; i < nrz; i++ )
  {
    znrt += dznr;
    if( near != 0)
    {
      cth= cos( TA* znrt);
      sth= sin( TA* znrt);
    }

    ynrt= ynr- dynr;
    for( j = 0; j < nry; j++ )
    {
      ynrt += dynr;
      if( near != 0)
      {
	cph= cos( TA* ynrt);
	sph= sin( TA* ynrt);
      }

      xnrt= xnr- dxnr;
      for( kk = 0; kk < nrx; kk++ )
      {
	xnrt += dxnr;
	if( near != 0)
	{
	  xob= xnrt* sth* cph;
	  yob= xnrt* sth* sph;
	  zob= xnrt* cth;
	}
	else
	{
	  xob= xnrt;
	  yob= ynrt;
	  zob= znrt;
	}

	tmp1= xob/ wlam;
	tmp2= yob/ wlam;
	tmp3= zob/ wlam;

	if( nfeh != 1)
	  nefld( tmp1, tmp2, tmp3, &ex, &ey, &ez);
	else
	  nhfld( tmp1, tmp2, tmp3, &ex, &ey, &ez);

	tmp1= cabs( ex);
	tmp2= cang( ex);
	tmp3= cabs( ey);
	tmp4= cang( ey);
	tmp5= cabs( ez);
	tmp6= cang( ez);

	fprintf( output_fp, "\n"
	    " %9.4f %9.4f %9.4f  %11.4E %7.2f  %11.4E %7.2f  %11.4E %7.2f",
	    xob, yob, zob, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );

	if( iplp1 != 2)
	  continue;

	if( iplp4 < 0 )
	  xxx= xob;
	else
	  if( iplp4 == 0 )
	    xxx= yob;
	  else
	    xxx= zob;

	if( iplp2 == 2)
	{
	  switch( iplp3 )
	  {
	    case 1:
	      fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", xxx, tmp1, tmp2 );
	      break;
	    case 2:
	      fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", xxx, tmp3, tmp4 );
	      break;
	    case 3:
	      fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", xxx, tmp5, tmp6 );
	      break;
	    case 4:
	      fprintf( plot_fp, "%12.4E %12.4E %12.4E %12.4E %12.4E %12.4E %12.4E\n",
		  xxx, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );
	  }
	  continue;
	}

	if( iplp2 != 1)
	  continue;

	switch( iplp3 )
	{
	  case 1:
	    fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", xxx, creal(ex), cimag(ex) );
	    break;
	  case 2:
	    fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", xxx, creal(ey), cimag(ey) );
	    break;
	  case 3:
	    fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", xxx, creal(ez), cimag(ez) );
	    break;
	  case 4:
	    fprintf( plot_fp, "%12.4E %12.4E %12.4E %12.4E %12.4E %12.4E %12.4E\n",
		xxx,creal(ex),cimag(ex),creal(ey),cimag(ey),creal(ez),cimag(ez) );
	}
      } /* for( kk = 0; kk < nrx; kk++ ) */

    } /* for( j = 0; j < nry; j++ ) */

  } /* for( i = 0; i < nrz; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* nhfld computes the near field at specified points in space after */
/* the structure currents have been computed. */

void nhfld( double xob, double yob, double zob,
    complex double *hx, complex double *hy, complex double *hz )
{
  int i, jc;
  double ax, zp;
  complex double acx, bcx, ccx;

  *hx=CPLX_00;
  *hy=CPLX_00;
  *hz=CPLX_00;
  ax=0.;

  if( n != 0)
  {
    for( i = 0; i < n; i++ )
    {
      xj= xob- x[i];
      yj= yob- y[i];
      zj= zob- z[i];
      zp= cab[i]* xj+ sab[i]* yj+ salp[i]* zj;

      if( fabs( zp) > 0.5001* si[i])
	continue;

      zp= xj* xj+ yj* yj+ zj* zj- zp* zp;
      xj= bi[i];

      if( zp > 0.9* xj* xj)
	continue;

      ax= xj;
      break;
    }

    for( i = 0; i < n; i++ )
    {
      s= si[i];
      b= bi[i];
      xj= x[i];
      yj= y[i];
      zj= z[i];
      cabj= cab[i];
      sabj= sab[i];
      salpj= salp[i];
      hsfld( xob, yob, zob, ax);
      acx= cmplx( air[i], aii[i]);
      bcx= cmplx( bir[i], bii[i]);
      ccx= cmplx( cir[i], cii[i]);
      *hx += exk* acx+ exs* bcx+ exc* ccx;
      *hy += eyk* acx+ eys* bcx+ eyc* ccx;
      *hz += ezk* acx+ ezs* bcx+ ezc* ccx;
    }

    if( m == 0)
      return;

  } /* if( n != 0) */

  jc= n-1;
  for( i = 0; i < m; i++ )
  {
    s= pbi[i];
    xj= px[i];
    yj= py[i];
    zj= pz[i];
    t1xj= t1x[i];
    t1yj= t1y[i];
    t1zj= t1z[i];
    t2xj= t2x[i];
    t2yj= t2y[i];
    t2zj= t2z[i];
    hintg( xob, yob, zob);
    jc += 3;
    acx= t1xj* cur[jc-2]+ t1yj* cur[jc-1]+ t1zj* cur[jc];
    bcx= t2xj* cur[jc-2]+ t2yj* cur[jc-1]+ t2zj* cur[jc];
    *hx= *hx+ acx* exk+ bcx* exs;
    *hy= *hy+ acx* eyk+ bcx* eys;
    *hz= *hz+ acx* ezk+ bcx* ezs;
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* patch generates and modifies patch geometry data */
void patch( int nx, int ny,
    double ax1, double ay1, double az1,
    double ax2, double ay2, double az2,
    double ax3, double ay3, double az3,
    double ax4, double ay4, double az4 )
{
  int mi, ntp, iy, ix, mreq;
  double s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., xst=0.;
  double znv, xnv, ynv, xa, xn2, yn2, zn2, salpn, xs, ys, zs, xt, yt, zt;

  /* new patches.  for nx=0, ny=1,2,3,4 patch is (respectively) */;
  /* arbitrary, rectagular, triangular, or quadrilateral. */
  /* for nx and ny  > 0 a rectangular surface is produced with */
  /* nx by ny rectangular patches. */

  m++;
  mi= m-1;

  /* Reallocate patch buffers */
  mreq = m * sizeof(double);
  mem_realloc( (void *)&px, mreq );
  mem_realloc( (void *)&py, mreq );
  mem_realloc( (void *)&pz, mreq );
  mem_realloc( (void *)&t1x, mreq );
  mem_realloc( (void *)&t1y, mreq );
  mem_realloc( (void *)&t1z, mreq );
  mem_realloc( (void *)&t2x, mreq );
  mem_realloc( (void *)&t2y, mreq );
  mem_realloc( (void *)&t2z, mreq );
  mem_realloc( (void *)&pbi, mreq );
  mem_realloc( (void *)&psalp, mreq );

  if( nx > 0)
    ntp=2;
  else
    ntp= ny;

  if( ntp <= 1)
  {
    px[mi]= ax1;
    py[mi]= ay1;
    pz[mi]= az1;
    pbi[mi]= az2;
    znv= cos( ax2);
    xnv= znv* cos( ay2);
    ynv= znv* sin( ay2);
    znv= sin( ax2);
    xa= sqrt( xnv* xnv+ ynv* ynv);

    if( xa >= 1.0e-6)
    {
      t1x[mi]=- ynv/ xa;
      t1y[mi]= xnv/ xa;
      t1z[mi]=0.;
    }
    else
    {
      t1x[mi]=1.;
      t1y[mi]=0.;
      t1z[mi]=0.;
    }

  } /* if( ntp <= 1) */
  else
  {
    s1x= ax2- ax1;
    s1y= ay2- ay1;
    s1z= az2- az1;
    s2x= ax3- ax2;
    s2y= ay3- ay2;
    s2z= az3- az2;

    if( nx != 0)
    {
      s1x= s1x/ nx;
      s1y= s1y/ nx;
      s1z= s1z/ nx;
      s2x= s2x/ ny;
      s2y= s2y/ ny;
      s2z= s2z/ ny;
    }

    xnv= s1y* s2z- s1z* s2y;
    ynv= s1z* s2x- s1x* s2z;
    znv= s1x* s2y- s1y* s2x;
    xa= sqrt( xnv* xnv+ ynv* ynv+ znv* znv);
    xnv= xnv/ xa;
    ynv= ynv/ xa;
    znv= znv/ xa;
    xst= sqrt( s1x* s1x+ s1y* s1y+ s1z* s1z);
    t1x[mi]= s1x/ xst;
    t1y[mi]= s1y/ xst;
    t1z[mi]= s1z/ xst;

    if( ntp <= 2)
    {
      px[mi]= ax1+.5*( s1x+ s2x);
      py[mi]= ay1+.5*( s1y+ s2y);
      pz[mi]= az1+.5*( s1z+ s2z);
      pbi[mi]= xa;
    }
    else
    {
      if( ntp != 4)
      {
	px[mi]=( ax1+ ax2+ ax3)/3.;
	py[mi]=( ay1+ ay2+ ay3)/3.;
	pz[mi]=( az1+ az2+ az3)/3.;
	pbi[mi]=.5* xa;
      }
      else
      {
	s1x= ax3- ax1;
	s1y= ay3- ay1;
	s1z= az3- az1;
	s2x= ax4- ax1;
	s2y= ay4- ay1;
	s2z= az4- az1;
	xn2= s1y* s2z- s1z* s2y;
	yn2= s1z* s2x- s1x* s2z;
	zn2= s1x* s2y- s1y* s2x;
	xst= sqrt( xn2* xn2+ yn2* yn2+ zn2* zn2);
	salpn=1./(3.*( xa+ xst));
	px[mi]=( xa*( ax1+ ax2+ ax3)+ xst*( ax1+ ax3+ ax4))* salpn;
	py[mi]=( xa*( ay1+ ay2+ ay3)+ xst*( ay1+ ay3+ ay4))* salpn;
	pz[mi]=( xa*( az1+ az2+ az3)+ xst*( az1+ az3+ az4))* salpn;
	pbi[mi]=.5*( xa+ xst);
	s1x=( xnv* xn2+ ynv* yn2+ znv* zn2)/ xst;

	if( s1x <= 0.9998)
	{
	  fprintf( output_fp,
	      "\n  ERROR -- CORNERS OF QUADRILATERAL"
	      " PATCH DO NOT LIE IN A PLANE" );
	  stop(-1);
	}

      } /* if( ntp != 4) */

    } /* if( ntp <= 2) */

  } /* if( ntp <= 1) */

  t2x[mi]= ynv* t1z[mi]- znv* t1y[mi];
  t2y[mi]= znv* t1x[mi]- xnv* t1z[mi];
  t2z[mi]= xnv* t1y[mi]- ynv* t1x[mi];
  psalp[mi]=1.;

  if( nx != 0)
  {
    m += nx*ny-1;

    /* Reallocate patch buffers */
    mreq = m * sizeof(double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

    xn2= px[mi]- s1x- s2x;
    yn2= py[mi]- s1y- s2y;
    zn2= pz[mi]- s1z- s2z;
    xs= t1x[mi];
    ys= t1y[mi];
    zs= t1z[mi];
    xt= t2x[mi];
    yt= t2y[mi];
    zt= t2z[mi];

    for( iy = 0; iy < ny; iy++ )
    {
      xn2 += s2x;
      yn2 += s2y;
      zn2 += s2z;

      for( ix = 1; ix <= nx; ix++ )
      {
	xst= (double)ix;
	px[mi]= xn2+ xst* s1x;
	py[mi]= yn2+ xst* s1y;
	pz[mi]= zn2+ xst* s1z;
	pbi[mi]= xa;
	psalp[mi]=1.;
	t1x[mi]= xs;
	t1y[mi]= ys;
	t1z[mi]= zs;
	t2x[mi]= xt;
	t2y[mi]= yt;
	t2z[mi]= zt;
	mi++;
      } /* for( ix = 0; ix < nx; ix++ ) */

    } /* for( iy = 0; iy < ny; iy++ ) */

  } /* if( nx != 0) */

  ipsym=0;
  np= n;
  mp= m;

  return;
}

/*-----------------------------------------------------------------------*/

/*** this function was an 'entry point' (part of) 'patch()' ***/
void subph( int nx, int ny )
{
  int mia, ix, iy, mi, mreq;
  double xs, ys, zs, xa, xst, s1x, s1y, s1z, s2x, s2y, s2z, saln, xt, yt;

  /* Shift patches to make room for new ones */
  m += 3;

  /* Reallocate patch buffers */
  mreq = m * sizeof(double);
  mem_realloc( (void *)&px, mreq );
  mem_realloc( (void *)&py, mreq );
  mem_realloc( (void *)&pz, mreq );
  mem_realloc( (void *)&t1x, mreq );
  mem_realloc( (void *)&t1y, mreq );
  mem_realloc( (void *)&t1z, mreq );
  mem_realloc( (void *)&t2x, mreq );
  mem_realloc( (void *)&t2y, mreq );
  mem_realloc( (void *)&t2z, mreq );
  mem_realloc( (void *)&pbi, mreq );
  mem_realloc( (void *)&psalp, mreq );

  if( (ny == 0) && (nx != m) )
  {
    for( iy = m-1; iy >= nx; iy-- )
    {
      px[iy]= px[iy-3];
      py[iy]= py[iy-3];
      pz[iy]= pz[iy-3];
      pbi[iy]= pbi[iy-3];
      psalp[iy]= psalp[iy-3];
      t1x[iy]= t1x[iy-3];
      t1y[iy]= t1y[iy-3];
      t1z[iy]= t1z[iy-3];
      t2x[iy]= t2x[iy-3];
      t2y[iy]= t2y[iy-3];
      t2z[iy]= t2z[iy-3];
    }

  } /* if( (ny <= 0) || (nx != m) ) */

  /* divide patch for connection */
  mi= nx-1;
  xs= px[mi];
  ys= py[mi];
  zs= pz[mi];
  xa= pbi[mi]/4.;
  xst= sqrt( xa)/2.;
  s1x= t1x[mi];
  s1y= t1y[mi];
  s1z= t1z[mi];
  s2x= t2x[mi];
  s2y= t2y[mi];
  s2z= t2z[mi];
  saln= psalp[mi];
  xt= xst;
  yt= xst;

  if( ny == 0)
    mia= mi;
  else
  {
    mp++;
    mia= m-1;
  }

  for( ix = 1; ix <= 4; ix++ )
  {
    px[mia]= xs+ xt* s1x+ yt* s2x;
    py[mia]= ys+ xt* s1y+ yt* s2y;
    pz[mia]= zs+ xt* s1z+ yt* s2z;
    pbi[mia]= xa;
    t1x[mia]= s1x;
    t1y[mia]= s1y;
    t1z[mia]= s1z;
    t2x[mia]= s2x;
    t2y[mia]= s2y;
    t2z[mia]= s2z;
    psalp[mia]= saln;

    if( ix == 2)
      yt=- yt;

    if( (ix == 1) || (ix == 3) )
      xt=- xt;

    mia++;
  }

  if( nx <= mp)
    mp += 3;

  if( ny > 0 )
    pz[mi]=10000.;

  return;
}

/*-----------------------------------------------------------------------*/

/* integrate over patches at wire connection point */
void pcint( double xi, double yi, double zi, double cabi,
    double sabi, double salpi, complex double *e )
{
  int nint, i1, i2;
  double d, ds, da, gcon, fcon, xxj, xyj, xzj, xs, s1;
  double xss, yss, zss, s2x, s2, g1, g2, g3, g4, f2, f1;
  complex double e1, e2, e3, e4, e5, e6, e7, e8, e9;

  nint = 10;
  d= sqrt( s)*.5;
  ds=4.* d/ (double) nint;
  da= ds* ds;
  gcon=1./ s;
  fcon=1./(2.* TP* d);
  xxj= xj;
  xyj= yj;
  xzj= zj;
  xs= s;
  s= da;
  s1= d+ ds*.5;
  xss= xj+ s1*( t1xj+ t2xj);
  yss= yj+ s1*( t1yj+ t2yj);
  zss= zj+ s1*( t1zj+ t2zj);
  s1= s1+ d;
  s2x= s1;
  e1=CPLX_00;
  e2=CPLX_00;
  e3=CPLX_00;
  e4=CPLX_00;
  e5=CPLX_00;
  e6=CPLX_00;
  e7=CPLX_00;
  e8=CPLX_00;
  e9=CPLX_00;

  for( i1 = 0; i1 < nint; i1++ )
  {
    s1= s1- ds;
    s2= s2x;
    xss= xss- ds* t1xj;
    yss= yss- ds* t1yj;
    zss= zss- ds* t1zj;
    xj= xss;
    yj= yss;
    zj= zss;

    for( i2 = 0; i2 < nint; i2++ )
    {
      s2= s2- ds;
      xj= xj- ds* t2xj;
      yj= yj- ds* t2yj;
      zj= zj- ds* t2zj;
      unere( xi, yi, zi);
      exk= exk* cabi+ eyk* sabi+ ezk* salpi;
      exs= exs* cabi+ eys* sabi+ ezs* salpi;
      g1=( d+ s1)*( d+ s2)* gcon;
      g2=( d- s1)*( d+ s2)* gcon;
      g3=( d- s1)*( d- s2)* gcon;
      g4=( d+ s1)*( d- s2)* gcon;
      f2=( s1* s1+ s2* s2)* TP;
      f1= s1/ f2-( g1- g2- g3+ g4)* fcon;
      f2= s2/ f2-( g1+ g2- g3- g4)* fcon;
      e1= e1+ exk* g1;
      e2= e2+ exk* g2;
      e3= e3+ exk* g3;
      e4= e4+ exk* g4;
      e5= e5+ exs* g1;
      e6= e6+ exs* g2;
      e7= e7+ exs* g3;
      e8= e8+ exs* g4;
      e9= e9+ exk* f1+ exs* f2;

    } /* for( i2 = 0; i2 < nint; i2++ ) */

  } /* for( i1 = 0; i1 < nint; i1++ ) */

  e[0]= e1;
  e[1]= e2;
  e[2]= e3;
  e[3]= e4;
  e[4]= e5;
  e[5]= e6;
  e[6]= e7;
  e[7]= e8;
  e[8]= e9;
  xj= xxj;
  yj= xyj;
  zj= xzj;
  s= xs;

  return;
}

/*-----------------------------------------------------------------------*/

/* prnt sets up the print formats for impedance loading */
void prnt( int in1, int in2, int in3, double fl1, double fl2,
    double fl3, double fl4, double fl5, double fl6, char *ia, int ichar )
{
  /* record to be output and buffer used to make it */
  char record[101+ichar*4], buf[15];
  int in[3], i1, i;
  double fl[6];

  in[0]= in1;
  in[1]= in2;
  in[2]= in3;
  fl[0]= fl1;
  fl[1]= fl2;
  fl[2]= fl3;
  fl[3]= fl4;
  fl[4]= fl5;
  fl[5]= fl6;

  /* integer format */
  i1=0;
  strcpy( record, "\n " );

  if( (in1 == 0) && (in2 == 0) && (in3 == 0) )
  {
    strcat( record, " ALL" );
    i1=1;
  }

  for( i = i1; i < 3; i++ )
  {
    if( in[i] == 0)
      strcat( record, "     " );
    else
    {
      sprintf( buf, "%5d", in[i] );
      strcat( record, buf );
    }
  }

  /* floating point format */
  for( i = 0; i < 6; i++ )
  {
    if( fabs( fl[i]) >= 1.0e-20 )
    {
      sprintf( buf, " %11.4E", fl[i] );
      strcat( record, buf );
    }
    else
      strcat( record, "            " );
  }

  strcat( record, "   " );
  strcat( record, ia );
  fprintf( output_fp, "%s", record );

  return;
}

/*-----------------------------------------------------------------------*/

/* fill incident field array for charge discontinuity voltage source */
void qdsrc( int is, complex double v, complex double *e )
{
  int i, jx, j, jp1, ipr, ij, i1;
  double xi, yi, zi, ai, cabi, sabi, salpi, tx, ty, tz;
  complex double curd, etk, ets, etc;

  is--;
  i= icon1[is];
  icon1[is]=0;
  tbf( is+1,0);
  icon1[is]= i;
  s= si[is]*.5;
  curd= CCJ* v/(( log(2.* s/ bi[is])-1.)*( bx[jsno-1]*
	cos( TP* s)+ cx[jsno-1]* sin( TP* s))* wlam);
  vqds[nqds]= v;
  iqds[nqds]= is+1;
  nqds++;

  for( jx = 0; jx < jsno; jx++ )
  {
    j= jco[jx]-1;
    jp1 = j+1;
    s= si[j];
    b= bi[j];
    xj= x[j];
    yj= y[j];
    zj= z[j];
    cabj= cab[j];
    sabj= sab[j];
    salpj= salp[j];

    if( iexk != 0)
    {
      ipr= icon1[j];

      if( ipr < 0 )
      {
	ipr=- ipr;
	ipr--;
	if( -icon1[ipr-1] != jp1 )
	  ind1=2;
	else
	{
	  xi= fabs( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	  if( (xi < 0.999999) || (fabs(bi[ipr]/b-1.) > 1.0e-6) )
	    ind1=2;
	  else
	    ind1=0;
	}
      }  /* if( ipr < 0 ) */
      else
	if( ipr == 0 )
	  ind1=1;
	else /* ipr > 0 */
	{
	  ipr--;
	  if( ipr != j )
	  {
	    if( icon2[ipr] != jp1)
	      ind1=2;
	    else
	    {
	      xi= fabs( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	      if( (xi < 0.999999) || (fabs(bi[ipr]/b-1.) > 1.0e-6) )
		ind1=2;
	      else
		ind1=0;
	    }
	  } /* if( ipr != j ) */
	  else
	  {
	    if( cabj* cabj+ sabj* sabj > 1.0e-8)
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* else */

      ipr= icon2[j];
      if( ipr < 0 )
      {
	ipr = -ipr;
	ipr--;
	if( -icon2[ipr] != jp1 )
	  ind1=2;
	else
	{
	  xi= fabs( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	  if( (xi < 0.999999) || (fabs(bi[ipr]/b-1.) > 1.0e-6) )
	    ind1=2;
	  else
	    ind1=0;
	}
      } /* if( ipr < 0 ) */
      else
	if( ipr == 0 )
	  ind2=1;
	else /* ipr > 0 */
	{
	  ipr--;
	  if( ipr != j )
	  {
	    if( icon1[ipr] != jp1)
	      ind2=2;
	    else
	    {
	      xi= fabs( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	      if( (xi < 0.999999) || (fabs(bi[ipr]/b-1.) > 1.0e-6) )
		ind2=2;
	      else
		ind2=0;
	    }
	  } /* if( ipr != j )*/
	  else
	  {
	    if( cabj* cabj+ sabj* sabj > 1.0e-8)
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* else */

    } /* if( iexk != 0) */

    for( i = 0; i < n; i++ )
    {
      ij= i- j;
      xi= x[i];
      yi= y[i];
      zi= z[i];
      ai= bi[i];
      efld( xi, yi, zi, ai, ij);
      cabi= cab[i];
      sabi= sab[i];
      salpi= salp[i];
      etk= exk* cabi+ eyk* sabi+ ezk* salpi;
      ets= exs* cabi+ eys* sabi+ ezs* salpi;
      etc= exc* cabi+ eyc* sabi+ ezc* salpi;
      e[i]= e[i]-( etk* ax[jx]+ ets* bx[jx]+ etc* cx[jx])* curd;
    }

    if( m != 0)
    {
      i1= n-1;
      for( i = 0; i < m; i++ )
      {
	xi= px[i];
	yi= py[i];
	zi= pz[i];
	hsfld( xi, yi, zi,0.);
	i1++;
	tx= t2x[i];
	ty= t2y[i];
	tz= t2z[i];
	etk= exk* tx+ eyk* ty+ ezk* tz;
	ets= exs* tx+ eys* ty+ ezs* tz;
	etc= exc* tx+ eyc* ty+ ezc* tz;
	e[i1] += ( etk* ax[jx]+ ets* bx[jx]+ etc* cx[jx] )* curd* psalp[i];
	i1++;
	tx= t1x[i];
	ty= t1y[i];
	tz= t1z[i];
	etk= exk* tx+ eyk* ty+ ezk* tz;
	ets= exs* tx+ eys* ty+ ezs* tz;
	etc= exc* tx+ eyc* ty+ ezc* tz;
	e[i1] += ( etk* ax[jx]+ ets* bx[jx]+ etc* cx[jx])* curd* psalp[i];
      }

    } /* if( m != 0) */

    if( nload > 0 )
      e[j] += zarray[j]* curd*(ax[jx]+ cx[jx]);

  } /* for( jx = 0; jx < jsno; jx++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* compute radiation pattern, gain, normalized gain */
void rdpat( void )
{
  char  *hpol[3] = { "LINEAR", "RIGHT ", "LEFT  " };
  char    hcir[] = " CIRCLE";
  char  *igtp[2] = { "----- POWER GAINS ----- ", "--- DIRECTIVE GAINS ---" };
  char  *igax[4] = { " MAJOR", " MINOR", " VERTC", " HORIZ" };
  char *igntp[5] =  { " MAJOR AXIS", "  MINOR AXIS",
    "    VERTICAL", "  HORIZONTAL", "       TOTAL " };

    char *hclif=NULL, *isens;
    int i, j, jump, itmp1, itmp2, kth, kph, itmp3, itmp4;
    double exrm=0., exra=0., prad, gcon, gcop, gmax, pint, tmp1, tmp2;
    double phi, pha, thet, tha, erdm=0., erda=0., ethm2, ethm, *gain = NULL;
    double etha, ephm2, ephm, epha, tilta, emajr2, eminr2, axrat;
    double dfaz, dfaz2, cdfaz, tstor1=0., tstor2, stilta, gnmj;
    double gnmn, gnv, gnh, gtot, tmp3, tmp4, da, tmp5, tmp6;
    complex double  eth, eph, erd;

    /* Allocate memory to gain buffer */
    if( inor > 0 )
      mem_alloc( (void *)&gain, nth*nph * sizeof(double) );

    if( ifar >= 2)
    {
      fprintf( output_fp, "\n\n\n"
	  "                                 "
	  "------ FAR FIELD GROUND PARAMETERS ------\n\n" );

      jump = FALSE;
      if( ifar > 3)
      {
	fprintf( output_fp, "\n"
	    "                                        "
	    "RADIAL WIRE GROUND SCREEN\n"
	    "                                        "
	    "%5d WIRES\n"
	    "                                        "
	    "WIRE LENGTH= %8.2f METERS\n"
	    "                                        "
	    "WIRE RADIUS= %10.3E METERS",
	    nradl, scrwlt, scrwrt );

	if( ifar == 4)
	  jump = TRUE;

      } /* if( ifar > 3) */

      if( ! jump )
      {
	if( (ifar == 2) || (ifar == 5) )
	  hclif= hpol[0];
	if( (ifar == 3) || (ifar == 6) )
	  hclif= hcir;

	cl= clt/ wlam;
	ch= cht/ wlam;
	zrati2= csqrt(1./ cmplx( epsr2,- sig2* wlam*59.96));

	fprintf( output_fp, "\n"
	    "                                        "
	    "%6s CLIFF\n"
	    "                                        "
	    "EDGE DISTANCE= %9.2f METERS\n"
	    "                                        "
	    "HEIGHT= %8.2f METERS\n"
	    "                                        "
	    "SECOND MEDIUM -\n"
	    "                                        "
	    "RELATIVE DIELECTRIC CONST.= %7.3f\n"
	    "                                        "
	    "CONDUCTIVITY= %10.3f MHOS",
	    hclif, clt, cht, epsr2, sig2 );

      } /* if( ! jump ) */

    } /* if( ifar >= 2) */

    if( ifar == 1)
    {
      fprintf( output_fp, "\n\n\n"
	  "                             "
	  "------- RADIATED FIELDS NEAR GROUND --------\n\n"
	  "    ------- LOCATION -------     --- E(THETA) ---    "
	  " ---- E(PHI) ----    --- E(RADIAL) ---\n"
	  "      RHO    PHI        Z           MAG    PHASE     "
	  "    MAG    PHASE        MAG     PHASE\n"
	  "    METERS DEGREES    METERS      VOLTS/M DEGREES   "
	  "   VOLTS/M DEGREES     VOLTS/M  DEGREES" );
    }
    else
    {
      itmp1=2* iax;
      itmp2= itmp1+1;

      fprintf( output_fp, "\n\n\n"
	  "                             "
	  "---------- RADIATION PATTERNS -----------\n" );

      if( rfld >= 1.0e-20)
      {
	exrm=1./ rfld;
	exra= rfld/ wlam;
	exra=-360.*( exra- floor( exra));

	fprintf( output_fp, "\n"
	    "                             "
	    "RANGE: %13.6E METERS\n"
	    "                             "
	    "EXP(-JKR)/R: %12.5E AT PHASE: %7.2f DEGREES\n",
	    rfld, exrm, exra );
      }

      fprintf( output_fp, "\n"
	  " ---- ANGLES -----     %23s      ---- POLARIZATION ----  "
	  " ---- E(THETA) ----    ----- E(PHI) ------\n"
	  "  THETA      PHI      %6s   %6s    TOTAL       AXIAL    "
	  "  TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE\n"
	  " DEGREES   DEGREES        DB       DB       DB       RATIO  "
	  " DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES",
	  igtp[ipd], igax[itmp1], igax[itmp2] );

    } /* if( ifar == 1) */

    if( (ixtyp == 0) || (ixtyp == 5) )
    {
      gcop= wlam* wlam*2.* PI/(376.73* pinr);
      prad= pinr- ploss- pnlr;
      gcon= gcop;
      if( ipd != 0)
	gcon= gcon* pinr/ prad;
    }
    else
      if( ixtyp == 4)
      {
	pinr=394.51* xpr6* xpr6* wlam* wlam;
	gcop= wlam* wlam*2.* PI/(376.73* pinr);
	prad= pinr- ploss- pnlr;
	gcon= gcop;
	if( ipd != 0)
	  gcon= gcon* pinr/ prad;
      }
      else
      {
	prad=0.;
	gcon=4.* PI/(1.+ xpr6* xpr6);
	gcop= gcon;
      }

    i=0;
    gmax=-1.e+10;
    pint=0.;
    tmp1= dph* TA;
    tmp2=.5* dth* TA;
    phi= phis- dph;

    for( kph = 1; kph <= nph; kph++ )
    {
      phi += dph;
      pha= phi* TA;
      thet= thets- dth;

      for( kth = 1; kth <= nth; kth++ )
      {
	thet += dth;
	if( (ksymp == 2) && (thet > 90.01) && (ifar != 1) )
	  continue;

	tha= thet* TA;
	if( ifar != 1)
	  ffld( tha, pha, &eth, &eph);
	else
	{
	  gfld( rfld/wlam, pha, thet/wlam,
	      &eth, &eph, &erd, zrati, ksymp);
	  erdm= cabs( erd);
	  erda= cang( erd);
	}

	ethm2= creal( eth* conj( eth));
	ethm= sqrt( ethm2);
	etha= cang( eth);
	ephm2= creal( eph* conj( eph));
	ephm= sqrt( ephm2);
	epha= cang( eph);

	/* elliptical polarization calc. */
	if( ifar != 1)
	{
	  if( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) )
	  {
	    tilta=0.;
	    emajr2=0.;
	    eminr2=0.;
	    axrat=0.;
	    isens= " ";
	  }
	  else
	  {
	    dfaz= epha- etha;
	    if( epha >= 0.)
	      dfaz2= dfaz-360.;
	    else
	      dfaz2= dfaz+360.;

	    if( fabs(dfaz) > fabs(dfaz2) )
	      dfaz= dfaz2;

	    cdfaz= cos( dfaz* TA);
	    tstor1= ethm2- ephm2;
	    tstor2=2.* ephm* ethm* cdfaz;
	    tilta=.5* atan2( tstor2, tstor1);
	    stilta= sin( tilta);
	    tstor1= tstor1* stilta* stilta;
	    tstor2= tstor2* stilta* cos( tilta);
	    emajr2=- tstor1+ tstor2+ ethm2;
	    eminr2= tstor1- tstor2+ ephm2;
	    if( eminr2 < 0.)
	      eminr2=0.;

	    axrat= sqrt( eminr2/ emajr2);
	    tilta= tilta* TD;
	    if( axrat <= 1.0e-5)
	      isens= hpol[0];
	    else
	      if( dfaz <= 0.)
		isens= hpol[1];
	      else
		isens= hpol[2];

	  } /* if( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) ) */

	  gnmj= db10( gcon* emajr2);
	  gnmn= db10( gcon* eminr2);
	  gnv = db10( gcon* ethm2);
	  gnh = db10( gcon* ephm2);
	  gtot= db10( gcon*(ethm2+ ephm2) );

	  if( inor > 0)
	  {
	    i++;
	    switch( inor )
	    {
	      case 1:
		tstor1= gnmj;
		break;

	      case 2:
		tstor1= gnmn;
		break;

	      case 3:
		tstor1= gnv;
		break;

	      case 4:
		tstor1= gnh;
		break;

	      case 5:
		tstor1= gtot;
	    }

	    gain[i-1]= tstor1;
	    if( tstor1 > gmax)
	      gmax= tstor1;

	  } /* if( inor > 0) */

	  if( iavp != 0)
	  {
	    tstor1= gcop*( ethm2+ ephm2);
	    tmp3= tha- tmp2;
	    tmp4= tha+ tmp2;

	    if( kth == 1)
	      tmp3= tha;
	    else
	      if( kth == nth)
		tmp4= tha;

	    da= fabs( tmp1*( cos( tmp3)- cos( tmp4)));
	    if( (kph == 1) || (kph == nph) )
	      da *=.5;
	    pint += tstor1* da;

	    if( iavp == 2)
	      continue;
	  }

	  if( iax != 1)
	  {
	    tmp5= gnmj;
	    tmp6= gnmn;
	  }
	  else
	  {
	    tmp5= gnv;
	    tmp6= gnh;
	  }

	  ethm= ethm* wlam;
	  ephm= ephm* wlam;

	  if( rfld >= 1.0e-20 )
	  {
	    ethm= ethm* exrm;
	    etha= etha+ exra;
	    ephm= ephm* exrm;
	    epha= epha+ exra;
	  }

	  fprintf( output_fp, "\n"
	      " %7.2f %9.2f  %8.2f %8.2f %8.2f %11.4f"
	      " %9.2f %6s %11.4E %9.2f %11.4E %9.2f",
	      thet, phi, tmp5, tmp6, gtot, axrat,
	      tilta, isens, ethm, etha, ephm, epha );

	  if( iplp1 != 3)
	    continue;

	  if( iplp3 != 0)
	  {
	    if( iplp2 == 1 )
	    {
	      if( iplp3 == 1 )
		fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", thet, ethm, etha );
	      else
		if( iplp3 == 2 )
		  fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", thet, ephm, epha );
	    }

	    if( iplp2 == 2 )
	    {
	      if( iplp3 == 1 )
		fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", phi, ethm, etha );
	      else
		if( iplp3 == 2 )
		  fprintf( plot_fp, "%12.4E %12.4E %12.4E\n", phi, ephm, epha );
	    }
	  }

	  if( iplp4 == 0 )
	    continue;

	  if( iplp2 == 1 )
	  {
	    switch( iplp4 )
	    {
	      case 1:
		fprintf( plot_fp, "%12.4E %12.4E\n", thet, tmp5 );
		break;
	      case 2:
		fprintf( plot_fp, "%12.4E %12.4E\n", thet, tmp6 );
		break;
	      case 3:
		fprintf( plot_fp, "%12.4E %12.4E\n", thet, gtot );
	    }
	  }

	  if( iplp2 == 2 )
	  {
	    switch( iplp4 )
	    {
	      case 1:
		fprintf( plot_fp, "%12.4E %12.4E\n", phi, tmp5 );
		break;
	      case 2:
		fprintf( plot_fp, "%12.4E %12.4E\n", phi, tmp6 );
		break;
	      case 3:
		fprintf( plot_fp, "%12.4E %12.4E\n", phi, gtot );
	    }
	  }

	  continue;
	} /* if( ifar != 1) */

	fprintf( output_fp, "\n"
	    " %9.2f %7.2f %9.2f  %11.4E %7.2f  %11.4E %7.2f  %11.4E %7.2f",
	    rfld, phi, thet, ethm, etha, ephm, epha, erdm, erda );

      } /* for( kth = 1; kth <= nth; kth++ ) */

    } /* for( kph = 1; kph <= nph; kph++ ) */

    if( iavp != 0)
    {
      tmp3= thets* TA;
      tmp4= tmp3+ dth* TA* (double)( nth-1);
      tmp3= fabs( dph* TA* (double)( nph-1)*( cos( tmp3)- cos( tmp4)));
      pint /= tmp3;
      tmp3 /= PI;

      fprintf( output_fp, "\n\n\n"
	  "  AVERAGE POWER GAIN: %11.4E - SOLID ANGLE"
	  " USED IN AVERAGING: (%+7.4f)*PI STERADIANS",
	  pint, tmp3 );
    }

    if( inor == 0)
      return;

    if( fabs( gnor) > 1.0e-20)
      gmax= gnor;
    itmp1=( inor-1);

    fprintf( output_fp,	"\n\n\n"
	"                             "
	" ---------- NORMALIZED GAIN ----------\n"
	"                                      %6s GAIN\n"
	"                                  "
	" NORMALIZATION FACTOR: %.2f db\n\n"
	"    ---- ANGLES ----                ---- ANGLES ----"
	"                ---- ANGLES ----\n"
	"    THETA      PHI        GAIN      THETA      PHI  "
	"      GAIN      THETA      PHI       GAIN\n"
	"   DEGREES   DEGREES        DB     DEGREES   DEGREES "
	"       DB     DEGREES   DEGREES       DB",
	igntp[itmp1], gmax );

    itmp2= nph* nth;
    itmp1=( itmp2+2)/3;
    itmp2= itmp1*3- itmp2;
    itmp3= itmp1;
    itmp4=2* itmp1;

    if( itmp2 == 2)
      itmp4--;

    for( i = 0; i < itmp1; i++ )
    {
      itmp3++;
      itmp4++;
      j= i/ nth;
      tmp1= thets+ (double)( i - j*nth )* dth;
      tmp2= phis+ (double)(j)* dph;
      j=( itmp3-1)/ nth;
      tmp3= thets+ (double)( itmp3- j* nth-1)* dth;
      tmp4= phis+ (double)(j)* dph;
      j=( itmp4-1)/ nth;
      tmp5= thets+ (double)( itmp4- j* nth-1)* dth;
      tmp6= phis+ (double)(j)* dph;
      tstor1= gain[i]- gmax;

      if( ((i+1) == itmp1) && (itmp2 != 0) )
      {
	if( itmp2 != 2)
	{
	  tstor2= gain[itmp3-1]- gmax;
	  fprintf( output_fp, "\n"
	      " %9.2f %9.2f %9.2f   %9.2f %9.2f %9.2f   ",
	      tmp1, tmp2, tstor1, tmp3, tmp4, tstor2 );
	  return;
	}

	fprintf( output_fp, "\n"
	    " %9.2f %9.2f %9.2f   ",
	    tmp1, tmp2, tstor1 );
	return;

      } /* if( ((i+1) == itmp1) && (itmp2 != 0) ) */

      tstor2= gain[itmp3-1]- gmax;
      pint= gain[itmp4-1]- gmax;

      fprintf( output_fp, "\n"
	  " %9.2f %9.2f %9.2f   %9.2f %9.2f %9.2f   %9.2f %9.2f %9.2f",
	  tmp1, tmp2, tstor1, tmp3, tmp4, tstor2, tmp5, tmp6, pint );

    } /* for( i = 0; i < itmp1; i++ ) */

    free_ptr( (void *)&gain );

    return;
}

/*-----------------------------------------------------------------------*/

void readgm( char *gm, int *i1, int *i2, double *x1, double *y1,
    double *z1, double *x2, double *y2, double *z2, double *rad )
{
  char line_buf[134];
  int nlin, i, line_idx;
  int nint = 2, nflt = 7;
  int iarr[2] = { 0, 0 };
  double rarr[7] = { 0., 0., 0., 0., 0., 0., 0. };


  /* read a line from input file */
  load_line( line_buf, input_fp );

  /* get line length */
  nlin= strlen( line_buf );

  /* abort if card's mnemonic too short or missing */
  if( nlin < 2 )
  {
    fprintf( output_fp,
	"\n  GEOMETRY DATA CARD ERROR:"
	"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
    stop(-1);
  }

  /* extract card's mnemonic code */
  strncpy( gm, line_buf, 2 );
  gm[2] = '\0';

  /* Exit if "XT" command read (for testing) */
  if( strcmp( gm, "XT" ) == 0 )
  {
    fprintf( stderr,
	"\nnec2c: Exiting after an \"XT\" command in readgm()\n" );
    fprintf( output_fp,
	"\n\n  nec2c: Exiting after an \"XT\" command in readgm()" );
    stop(0);
  }

  /* Return if only mnemonic on card */
  if( nlin == 2 )
  {
    *i1 = *i2 = 0;
    *x1 = *y1 = *z1 = *x2 = *y2 = *z2 = *rad = 0.;
    return;
  }

  /* read integers from line */
  line_idx = 1;
  for( i = 0; i < nint; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
	    (line_buf[  line_idx] >  '9')) &&
	    (line_buf[  line_idx] != '+')  &&
	    (line_buf[  line_idx] != '-') )
      if( (line_buf[line_idx] == '\0') )
      {
	*i1= iarr[0];
	*i2= iarr[1];
	*x1= rarr[0];
	*y1= rarr[1];
	*z1= rarr[2];
	*x2= rarr[3];
	*y2= rarr[4];
	*z2= rarr[5];
	*rad= rarr[6];
	return;
      }

    /* read an integer from line */
    iarr[i] = atoi( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while( (line_buf[++line_idx] != ' ') &&
	(line_buf[  line_idx] != ',') &&
	(line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if( ((line_buf[line_idx] <  '0')  ||
	   (line_buf[line_idx] >  '9')) &&
	   (line_buf[line_idx] != '+')  &&
	   (line_buf[line_idx] != '-') )
      {
	fprintf( output_fp,
	    "\n  GEOMETRY DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
	    gm, line_buf[line_idx], (line_idx+1)  );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *x1= rarr[0];
      *y1= rarr[1];
      *z1= rarr[2];
      *x2= rarr[3];
      *y2= rarr[4];
      *z2= rarr[5];
      *rad= rarr[6];
      return;
    }

  } /* for( i = 0; i < nint; i++ ) */

  /* read doubles from line */
  for( i = 0; i < nflt; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
	    (line_buf[  line_idx] >  '9')) &&
	    (line_buf[  line_idx] != '+')  &&
	    (line_buf[  line_idx] != '-')  &&
	    (line_buf[  line_idx] != '.') )
      if( (line_buf[line_idx] == '\0') )
      {
	*i1= iarr[0];
	*i2= iarr[1];
	*x1= rarr[0];
	*y1= rarr[1];
	*z1= rarr[2];
	*x2= rarr[3];
	*y2= rarr[4];
	*z2= rarr[5];
	*rad= rarr[6];
	return;
      }

    /* read a double from line */
    rarr[i] = atof( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while( (line_buf[++line_idx] != ' ') &&
	   (line_buf[  line_idx] != ',') &&
	   (line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if( ((line_buf[line_idx] <  '0')  ||
	   (line_buf[line_idx] >  '9')) &&
	   (line_buf[line_idx] != '.')  &&
	   (line_buf[line_idx] != '+')  &&
	   (line_buf[line_idx] != '-')  &&
	   (line_buf[line_idx] != 'E')  &&
	   (line_buf[line_idx] != 'e') )
      {
	fprintf( output_fp,
	    "\n  GEOMETRY DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d.\n",
	    gm, line_buf[line_idx], (line_idx+1) );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *x1= rarr[0];
      *y1= rarr[1];
      *z1= rarr[2];
      *x2= rarr[3];
      *y2= rarr[4];
      *z2= rarr[5];
      *rad= rarr[6];
      return;
    }

  } /* for( i = 0; i < nflt; i++ ) */

  *i1  = iarr[0];
  *i2  = iarr[1];
  *x1  = rarr[0];
  *y1  = rarr[1];
  *z1  = rarr[2];
  *x2  = rarr[3];
  *y2  = rarr[4];
  *z2  = rarr[5];
  *rad = rarr[6];

  return;
}

/*-----------------------------------------------------------------------*/

void readmn( char *gm, int *i1, int *i2, int *i3, int *i4,
    double *f1, double *f2, double *f3,
    double *f4, double *f5, double *f6 )
{
  char line_buf[134];
  int nlin, i, line_idx;
  int nint = 4, nflt = 6;
  int iarr[4] = { 0, 0, 0, 0 };
  double rarr[6] = { 0., 0., 0., 0., 0., 0. };

  /* read a line from input file */
  load_line( line_buf, input_fp );

  /* get line length */
  nlin= strlen( line_buf );

  /* abort if card's mnemonic too short or missing */
  if( nlin < 2 )
  {
    fprintf( output_fp,
	"\n  COMMAND DATA CARD ERROR:"
	"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
    stop(-1);
  }

  /* extract card's mnemonic code */
  strncpy( gm, line_buf, 2 );
  gm[2] = '\0';

  /* Exit if "XT" command read (for testing) */
  if( strcmp( gm, "XT" ) == 0 )
  {
    fprintf( stderr,
	"\nnec2c: Exiting after an \"XT\" command in readgm()\n" );
    fprintf( output_fp,
	"\n\n  nec2c: Exiting after an \"XT\" command in readgm()" );
    stop(0);
  }

  /* Return if only mnemonic on card */
  if( nlin == 2 )
  {
    *i1 = *i2 = *i3 = *i4 = 0;
    *f1 = *f2 = *f3 = *f4 = *f5 = *f6 = 0.0;
    return;
  }

  /* read integers from line */
  line_idx = 1;
  for( i = 0; i < nint; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
	    (line_buf[  line_idx] >  '9')) &&
	    (line_buf[  line_idx] != '+')  &&
	    (line_buf[  line_idx] != '-') )
      if( (line_buf[line_idx] == '\0') )
      {
	*i1= iarr[0];
	*i2= iarr[1];
	*i3= iarr[2];
	*i4= iarr[3];
	*f1= rarr[0];
	*f2= rarr[1];
	*f3= rarr[2];
	*f4= rarr[3];
	*f5= rarr[4];
	*f6= rarr[5];
	return;
      }

    /* read an integer from line */
    iarr[i] = atoi( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while( (line_buf[++line_idx] != ' ') &&
	   (line_buf[  line_idx] != ',') &&
	   (line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if( ((line_buf[line_idx] <  '0')  ||
	   (line_buf[line_idx] >  '9')) &&
	   (line_buf[line_idx] != '+')  &&
	   (line_buf[line_idx] != '-') )
      {
	fprintf( output_fp,
	    "\n  COMMAND DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
	    gm, line_buf[line_idx], (line_idx+1) );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *i3= iarr[2];
      *i4= iarr[3];
      *f1= rarr[0];
      *f2= rarr[1];
      *f3= rarr[2];
      *f4= rarr[3];
      *f5= rarr[4];
      *f6= rarr[5];
      return;
    }

  } /* for( i = 0; i < nint; i++ ) */

  /* read doubles from line */
  for( i = 0; i < nflt; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
	    (line_buf[  line_idx] >  '9')) &&
	    (line_buf[  line_idx] != '+')  &&
	    (line_buf[  line_idx] != '-')  &&
	    (line_buf[  line_idx] != '.') )
      if( (line_buf[line_idx] == '\0') )
      {
	*i1= iarr[0];
	*i2= iarr[1];
	*i3= iarr[2];
	*i4= iarr[3];
	*f1= rarr[0];
	*f2= rarr[1];
	*f3= rarr[2];
	*f4= rarr[3];
	*f5= rarr[4];
	*f6= rarr[5];
	return;
      }

    /* read a double from line */
    rarr[i] = atof( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' */
    line_idx--;
    while( (line_buf[++line_idx] != ' ') &&
	   (line_buf[  line_idx] != ',') &&
	   (line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if( ((line_buf[line_idx] <  '0')  ||
	   (line_buf[line_idx] >  '9')) &&
	   (line_buf[line_idx] != '.')  &&
	   (line_buf[line_idx] != '+')  &&
	   (line_buf[line_idx] != '-')  &&
	   (line_buf[line_idx] != 'E')  &&
	   (line_buf[line_idx] != 'e') )
      {
	fprintf( output_fp,
	    "\n  COMMAND DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d\n",
	    gm, line_buf[line_idx], (line_idx+1) );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *i3= iarr[2];
      *i4= iarr[3];
      *f1= rarr[0];
      *f2= rarr[1];
      *f3= rarr[2];
      *f4= rarr[3];
      *f5= rarr[4];
      *f6= rarr[5];
      return;
    }

  } /* for( i = 0; i < nflt; i++ ) */

  *i1= iarr[0];
  *i2= iarr[1];
  *i3= iarr[2];
  *i4= iarr[3];
  *f1= rarr[0];
  *f2= rarr[1];
  *f3= rarr[2];
  *f4= rarr[3];
  *f5= rarr[4];
  *f6= rarr[5];

  return;
}

/*-----------------------------------------------------------------------*/

/* reflc reflects partial structure along x,y, or z axes or rotates */
/* structure to complete a symmetric structure. */
void reflc( int ix, int iy, int iz, int itx, int nop )
{
  int iti, i, nx, itagi, k, mreq;
  double e1, e2, fnop, sam, cs, ss, xk, yk;

  np= n;
  mp= m;
  ipsym=0;
  iti= itx;

  if( ix >= 0)
  {
    if( nop == 0)
      return;

    ipsym=1;

    /* reflect along z axis */
    if( iz != 0)
    {
      ipsym=2;

      if( n > 0 )
      {
	/* Reallocate tags buffer */
	mem_realloc( (void *)&itag, (2*n+m) * sizeof(int) );

	/* Reallocate wire buffers */
	mreq = 2*n * sizeof(double);
	mem_realloc( (void *)&x, mreq );
	mem_realloc( (void *)&y, mreq );
	mem_realloc( (void *)&z, mreq );
	mem_realloc( (void *)&x2, mreq );
	mem_realloc( (void *)&y2, mreq );
	mem_realloc( (void *)&z2, mreq );
	mem_realloc( (void *)&bi, mreq );

	for( i = 0; i < n; i++ )
	{
	  nx= i+ n;
	  e1= z[i];
	  e2= z2[i];

	  if( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--SEGMENT %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  x[nx]= x[i];
	  y[nx]= y[i];
	  z[nx]=- e1;
	  x2[nx]= x2[i];
	  y2[nx]= y2[i];
	  z2[nx]=- e2;
	  itagi= itag[i];

	  if( itagi == 0)
	    itag[nx]=0;
	  if( itagi != 0)
	    itag[nx]= itagi+ iti;

	  bi[nx]= bi[i];

	} /* for( i = 0; i < n; i++ ) */

	n= n*2;
	iti= iti*2;

      } /* if( n > 0) */

      if( m > 0 )
      {
	/* Reallocate patch buffers */
	mreq = 2*m * sizeof(double);
	mem_realloc( (void *)&px, mreq );
	mem_realloc( (void *)&py, mreq );
	mem_realloc( (void *)&pz, mreq );
	mem_realloc( (void *)&t1x, mreq );
	mem_realloc( (void *)&t1y, mreq );
	mem_realloc( (void *)&t1z, mreq );
	mem_realloc( (void *)&t2x, mreq );
	mem_realloc( (void *)&t2y, mreq );
	mem_realloc( (void *)&t2z, mreq );
	mem_realloc( (void *)&pbi, mreq );
	mem_realloc( (void *)&psalp, mreq );

	for( i = 0; i < m; i++ )
	{
	  nx = i+m;
	  if( fabs(pz[i]) <= 1.0e-10)
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--PATCH %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  px[nx]= px[i];
	  py[nx]= py[i];
	  pz[nx]=- pz[i];
	  t1x[nx]= t1x[i];
	  t1y[nx]= t1y[i];
	  t1z[nx]=- t1z[i];
	  t2x[nx]= t2x[i];
	  t2y[nx]= t2y[i];
	  t2z[nx]=- t2z[i];
	  psalp[nx]=- psalp[i];
	  pbi[nx]= pbi[i];
	}

	m= m*2;

      } /* if( m >= m2) */

    } /* if( iz != 0) */

    /* reflect along y axis */
    if( iy != 0)
    {
      if( n > 0)
      {
	/* Reallocate tags buffer */
	mem_realloc( (void *)&itag, (2*n+m) * sizeof(int) );/*????*/

	/* Reallocate wire buffers */
	mreq = 2*n * sizeof(double);
	mem_realloc( (void *)&x, mreq );
	mem_realloc( (void *)&y, mreq );
	mem_realloc( (void *)&z, mreq );
	mem_realloc( (void *)&x2, mreq );
	mem_realloc( (void *)&y2, mreq );
	mem_realloc( (void *)&z2, mreq );
	mem_realloc( (void *)&bi, mreq );

	for( i = 0; i < n; i++ )
	{
	  nx= i+ n;
	  e1= y[i];
	  e2= y2[i];

	  if( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--SEGMENT %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  x[nx]= x[i];
	  y[nx]=- e1;
	  z[nx]= z[i];
	  x2[nx]= x2[i];
	  y2[nx]=- e2;
	  z2[nx]= z2[i];
	  itagi= itag[i];

	  if( itagi == 0)
	    itag[nx]=0;
	  if( itagi != 0)
	    itag[nx]= itagi+ iti;

	  bi[nx]= bi[i];

	} /* for( i = n2-1; i < n; i++ ) */

	n= n*2;
	iti= iti*2;

      } /* if( n >= n2) */

      if( m > 0 )
      {
	/* Reallocate patch buffers */
	mreq = 2*m * sizeof(double);
	mem_realloc( (void *)&px, mreq );
	mem_realloc( (void *)&py, mreq );
	mem_realloc( (void *)&pz, mreq );
	mem_realloc( (void *)&t1x, mreq );
	mem_realloc( (void *)&t1y, mreq );
	mem_realloc( (void *)&t1z, mreq );
	mem_realloc( (void *)&t2x, mreq );
	mem_realloc( (void *)&t2y, mreq );
	mem_realloc( (void *)&t2z, mreq );
	mem_realloc( (void *)&pbi, mreq );
	mem_realloc( (void *)&psalp, mreq );

	for( i = 0; i < m; i++ )
	{
	  nx= i+m;
	  if( fabs( py[i]) <= 1.0e-10)
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--PATCH %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  px[nx]= px[i];
	  py[nx]=- py[i];
	  pz[nx]= pz[i];
	  t1x[nx]= t1x[i];
	  t1y[nx]=- t1y[i];
	  t1z[nx]= t1z[i];
	  t2x[nx]= t2x[i];
	  t2y[nx]=- t2y[i];
	  t2z[nx]= t2z[i];
	  psalp[nx]=- psalp[i];
	  pbi[nx]= pbi[i];

	} /* for( i = m2; i <= m; i++ ) */

	m= m*2;

      } /* if( m >= m2) */

    } /* if( iy != 0) */

    /* reflect along x axis */
    if( ix == 0 )
      return;

    if( n > 0 )
    {
      /* Reallocate tags buffer */
      mem_realloc( (void *)&itag, (2*n+m) * sizeof(int) );/*????*/

      /* Reallocate wire buffers */
      mreq = 2*n * sizeof(double);
      mem_realloc( (void *)&x, mreq );
      mem_realloc( (void *)&y, mreq );
      mem_realloc( (void *)&z, mreq );
      mem_realloc( (void *)&x2, mreq );
      mem_realloc( (void *)&y2, mreq );
      mem_realloc( (void *)&z2, mreq );
      mem_realloc( (void *)&bi, mreq );

      for( i = 0; i < n; i++ )
      {
	nx= i+ n;
	e1= x[i];
	e2= x2[i];

	if( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
	{
	  fprintf( output_fp,
	      "\n  GEOMETRY DATA ERROR--SEGMENT %d"
	      " LIES IN PLANE OF SYMMETRY", i+1 );
	  stop(-1);
	}

	x[nx]=- e1;
	y[nx]= y[i];
	z[nx]= z[i];
	x2[nx]=- e2;
	y2[nx]= y2[i];
	z2[nx]= z2[i];
	itagi= itag[i];

	if( itagi == 0)
	  itag[nx]=0;
	if( itagi != 0)
	  itag[nx]= itagi+ iti;

	bi[nx]= bi[i];
      }

      n= n*2;

    } /* if( n > 0) */

    if( m == 0 )
      return;

    /* Reallocate patch buffers */
    mreq = 2*m * sizeof(double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

    for( i = 0; i < m; i++ )
    {
      nx= i+m;
      if( fabs( px[i]) <= 1.0e-10)
      {
	fprintf( output_fp,
	    "\n  GEOMETRY DATA ERROR--PATCH %d"
	    " LIES IN PLANE OF SYMMETRY", i+1 );
	stop(-1);
      }

      px[nx]=- px[i];
      py[nx]= py[i];
      pz[nx]= pz[i];
      t1x[nx]=- t1x[i];
      t1y[nx]= t1y[i];
      t1z[nx]= t1z[i];
      t2x[nx]=- t2x[i];
      t2y[nx]= t2y[i];
      t2z[nx]= t2z[i];
      psalp[nx]=- psalp[i];
      pbi[nx]= pbi[i];
    }

    m= m*2;
    return;

  } /* if( ix >= 0) */

  /* reproduce structure with rotation to form cylindrical structure */
  fnop= (double)nop;
  ipsym=-1;
  sam=TP/ fnop;
  cs= cos( sam);
  ss= sin( sam);

  if( n > 0)
  {
    n *= nop;
    nx= np;

    /* Reallocate tags buffer */
    mem_realloc( (void *)&itag, (n+m) * sizeof(int) );/*????*/

    /* Reallocate wire buffers */
    mreq = n * sizeof(double);
    mem_realloc( (void *)&x, mreq );
    mem_realloc( (void *)&y, mreq );
    mem_realloc( (void *)&z, mreq );
    mem_realloc( (void *)&x2, mreq );
    mem_realloc( (void *)&y2, mreq );
    mem_realloc( (void *)&z2, mreq );
    mem_realloc( (void *)&bi, mreq );

    for( i = nx; i < n; i++ )
    {
      k= i- np;
      xk= x[k];
      yk= y[k];
      x[i]= xk* cs- yk* ss;
      y[i]= xk* ss+ yk* cs;
      z[i]= z[k];
      xk= x2[k];
      yk= y2[k];
      x2[i]= xk* cs- yk* ss;
      y2[i]= xk* ss+ yk* cs;
      z2[i]= z2[k];
      bi[i]= bi[k];
      itagi= itag[k];

      if( itagi == 0)
	itag[i]=0;
      if( itagi != 0)
	itag[i]= itagi+ iti;
    }

  } /* if( n >= n2) */

  if( m == 0 )
    return;

  m *= nop;
  nx= mp;

  /* Reallocate patch buffers */
  mreq = m * sizeof(double);
  mem_realloc( (void *)&px, mreq  );
  mem_realloc( (void *)&py, mreq  );
  mem_realloc( (void *)&pz, mreq );
  mem_realloc( (void *)&t1x, mreq );
  mem_realloc( (void *)&t1y, mreq );
  mem_realloc( (void *)&t1z, mreq );
  mem_realloc( (void *)&t2x, mreq );
  mem_realloc( (void *)&t2y, mreq );
  mem_realloc( (void *)&t2z, mreq );
  mem_realloc( (void *)&pbi, mreq );
  mem_realloc( (void *)&psalp, mreq );

  for( i = nx; i < m; i++ )
  {
    k = i-mp;
    xk= px[k];
    yk= py[k];
    px[i]= xk* cs- yk* ss;
    py[i]= xk* ss+ yk* cs;
    pz[i]= pz[k];
    xk= t1x[k];
    yk= t1y[k];
    t1x[i]= xk* cs- yk* ss;
    t1y[i]= xk* ss+ yk* cs;
    t1z[i]= t1z[k];
    xk= t2x[k];
    yk= t2y[k];
    t2x[i]= xk* cs- yk* ss;
    t2y[i]= xk* ss+ yk* cs;
    t2z[i]= t2z[k];
    psalp[i]= psalp[k];
    pbi[i]= pbi[k];

  } /* for( i = nx; i < m; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* for the sommerfeld ground option, rom2 integrates over the source */
/* segment to obtain the total field due to ground.  the method of */
/* variable interval width romberg integration is used.  there are 9 */
/* field components - the x, y, and z components due to constant, */
/* sine, and cosine current distributions. */
void rom2( double a, double b, complex double *sum, double dmin )
{
  int i, ns, nt, flag=TRUE;
  int nts = 4, nx = 1, n = 9;
  double ze, ep, zend, dz=0., dzot=0., tmag1, tmag2, tr, ti;
  double z, s; /***also global***/
  double rx = 1.0e-4;
  complex double g1[9], g2[9], g3[9], g4[9], g5[9];
  complex double t00, t01[9], t10[9], t02, t11, t20[9];

  z= a;
  ze= b;
  s= b- a;

  if( s <= 0.)
  {
    fprintf( output_fp, "\n  ERROR - B LESS THAN A IN ROM2" );
    stop(-1);
  }

  ep= s/(1.e4* npm);
  zend= ze- ep;

  for( i = 0; i < n; i++ )
    sum[i]=CPLX_00;

  ns= nx;
  nt=0;
  sflds( z, g1);

  while( TRUE )
  {
    if( flag )
    {
      dz= s/ ns;
      if( z+ dz > ze)
      {
	dz= ze- z;
	if( dz <= ep)
	  return;
      }

      dzot= dz*.5;
      sflds( z+ dzot, g3);
      sflds( z+ dz, g5);

    } /* if( flag ) */

    tmag1=0.;
    tmag2=0.;

    /* evaluate 3 point romberg result and test convergence. */
    for( i = 0; i < n; i++ )
    {
      t00=( g1[i]+ g5[i])* dzot;
      t01[i]=( t00+ dz* g3[i])*.5;
      t10[i]=(4.* t01[i]- t00)/3.;
      if( i > 2)
	continue;

      tr= creal( t01[i]);
      ti= cimag( t01[i]);
      tmag1= tmag1+ tr* tr+ ti* ti;
      tr= creal( t10[i]);
      ti= cimag( t10[i]);
      tmag2= tmag2+ tr* tr+ ti* ti;

    } /* for( i = 0; i < n; i++ ) */

    tmag1= sqrt( tmag1);
    tmag2= sqrt( tmag2);
    test( tmag1, tmag2, &tr, 0., 0., &ti, dmin);

    if( tr <= rx)
    {
      for( i = 0; i < n; i++ )
	sum[i] += t10[i];
      nt += 2;

      z += dz;
      if( z > zend)
	return;

      for( i = 0; i < n; i++ )
	g1[i]= g5[i];

      if( (nt >= nts) && (ns > nx) )
      {
	ns= ns/2;
	nt=1;
      }
      flag = TRUE;
      continue;

    } /* if( tr <= rx) */

    sflds( z+ dz*.25, g2);
    sflds( z+ dz*.75, g4);
    tmag1=0.;
    tmag2=0.;

    /* evaluate 5 point romberg result and test convergence. */
    for( i = 0; i < n; i++ )
    {
      t02=( t01[i]+ dzot*( g2[i]+ g4[i]))*.5;
      t11=( 4.0 * t02- t01[i] )/3.;
      t20[i]=(16.* t11- t10[i])/15.;
      if( i > 2)
	continue;

      tr= creal( t11);
      ti= cimag( t11);
      tmag1= tmag1+ tr* tr+ ti* ti;
      tr= creal( t20[i]);
      ti= cimag( t20[i]);
      tmag2= tmag2+ tr* tr+ ti* ti;

    } /* for( i = 0; i < n; i++ ) */

    tmag1= sqrt( tmag1);
    tmag2= sqrt( tmag2);
    test( tmag1, tmag2, &tr, 0.,0., &ti, dmin);

    if( tr > rx)
    {
      nt=0;
      if( ns < npm )
      {
	ns= ns*2;
	dz= s/ ns;
	dzot= dz*.5;

	for( i = 0; i < n; i++ )
	{
	  g5[i]= g3[i];
	  g3[i]= g2[i];
	}

	flag=FALSE;
	continue;

      } /* if( ns < npm) */

      fprintf( output_fp,
	  "\n  ROM2 -- STEP SIZE LIMITED AT Z = %12.5E", z );

    } /* if( tr > rx) */

    for( i = 0; i < n; i++ )
      sum[i]= sum[i]+ t20[i];
    nt= nt+1;

    z= z+ dz;
    if( z > zend)
      return;

    for( i = 0; i < n; i++ )
      g1[i]= g5[i];

    flag = TRUE;
    if( (nt < nts) || (ns <= nx) )
      continue;

    ns= ns/2;
    nt=1;

  } /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* compute component of basis function i on segment is. */
void sbf( int i, int is, double *aa, double *bb, double *cc )
{
  int ix, jsno, june, jcox, jcoxx, jend, iend, njun1=0, njun2;
  double d, sig, pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;

  *aa=0.;
  *bb=0.;
  *cc=0.;
  june=0;
  jsno=0;
  pp=0.;
  ix=i-1;

  jcox= icon1[ix];
  if( jcox > PCHCON)
    jcox= i;
  jcoxx = jcox-1;

  jend=-1;
  iend=-1;
  sig=-1.;

  do
  {
    if( jcox != 0 )
    {
      if( jcox < 0 )
	jcox=- jcox;
      else
      {
	sig=- sig;
	jend=- jend;
      }

      jcoxx = jcox-1;
      jsno++;
      d= PI* si[jcoxx];
      sdh= sin( d);
      cdh= cos( d);
      sd=2.* sdh* cdh;

      if( d <= 0.015)
      {
	omc=4.* d* d;
	omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
      }
      else
	omc=1.- cdh* cdh+ sdh* sdh;

      aj=1./( log(1./( PI* bi[jcoxx]))-.577215664);
      pp -= omc/ sd* aj;

      if( jcox == is)
      {
	*aa= aj/ sd* sig;
	*bb= aj/(2.* cdh);
	*cc=- aj/(2.* sdh)* sig;
	june= iend;
      }

      if( jcox != i )
      {
	if( jend != 1)
	  jcox= icon1[jcoxx];
	else
	  jcox= icon2[jcoxx];

	if( abs(jcox) != i )
	{
	  if( jcox == 0 )
	  {
	    fprintf( output_fp,
		"\n  SBF - SEGMENT CONNECTION ERROR FOR SEGMENT %d", i);
	    stop(-1);
	  }
	  else
	    continue;
	}

      } /* if( jcox != i ) */
      else
	if( jcox == is)
	  *bb=- *bb;

      if( iend == 1)
	break;

    } /* if( jcox != 0 ) */

    pm=- pp;
    pp=0.;
    njun1= jsno;

    jcox= icon2[ix];
    if( jcox > PCHCON)
      jcox= i;

    jend=1;
    iend=1;
    sig=-1.;

  } /* do */
  while( jcox != 0 );

  njun2= jsno- njun1;
  d= PI* si[ix];
  sdh= sin( d);
  cdh= cos( d);
  sd=2.* sdh* cdh;
  cd= cdh* cdh- sdh* sdh;

  if( d <= 0.015)
  {
    omc=4.* d* d;
    omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
  }
  else
    omc=1.- cd;

  ap=1./( log(1./( PI* bi[ix])) -.577215664);
  aj= ap;

  if( njun1 == 0)
  {
    if( njun2 == 0)
    {
      *aa =-1.;
      qp= PI* bi[ix];
      xxi= qp* qp;
      xxi= qp*(1.-.5* xxi)/(1.- xxi);
      *cc=1./( cdh- xxi* sdh);
      return;
    }

    qp= PI* bi[ix];
    xxi= qp* qp;
    xxi= qp*(1.-.5* xxi)/(1.- xxi);
    qp=-( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));

    if( june == 1)
    {
      *aa=- *aa* qp;
      *bb=  *bb* qp;
      *cc=- *cc* qp;
      if( i != is)
	return;
    }

    *aa -= 1.;
    d = cd - xxi * sd;
    *bb += (sdh + ap * qp * (cdh - xxi * sdh)) / d;
    *cc += (cdh + ap * qp * (sdh + xxi * cdh)) / d;
    return;

  } /* if( njun1 == 0) */

  if( njun2 == 0)
  {
    qm= PI* bi[ix];
    xxi= qm* qm;
    xxi= qm*(1.-.5* xxi)/(1.- xxi);
    qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));

    if( june == -1)
    {
      *aa= *aa* qm;
      *bb= *bb* qm;
      *cc= *cc* qm;
      if( i != is)
	return;
    }

    *aa -= 1.;
    d= cd- xxi* sd;
    *bb += ( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
    *cc += ( cdh- aj* qm*( sdh+ xxi* cdh))/ d;
    return;

  } /* if( njun2 == 0) */

  qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
  qm=( ap* omc- pp* sd)/ qp;
  qp=-( aj* omc+ pm* sd)/ qp;

  if( june != 0 )
  {
    if( june < 0 )
    {
      *aa= *aa* qm;
      *bb= *bb* qm;
      *cc= *cc* qm;
    }
    else
    {
      *aa=- *aa* qp;
      *bb= *bb* qp;
      *cc=- *cc* qp;
    }

    if( i != is)
      return;

  } /* if( june != 0 ) */

  *aa -= 1.;
  *bb += ( aj* qm+ ap* qp)* sdh/ sd;
  *cc += ( aj* qm- ap* qp)* cdh/ sd;

  return;
}

/*-----------------------------------------------------------------------*/

/* sfldx returns the field due to ground for a current element on */
/* the source segment at t relative to the segment center. */
void sflds( double t, complex double *e )
{
  double xt, yt, zt, rhx, rhy, rhs, rho, phx, phy;
  double cph, sph, zphs, r2s, rk, sfac, thet;
  complex double  erv, ezv, erh, ezh, eph, er, et, hrv, hzv, hrh;

  xt= xj+ t* cabj;
  yt= yj+ t* sabj;
  zt= zj+ t* salpj;
  rhx= xo- xt;
  rhy= yo- yt;
  rhs= rhx* rhx+ rhy* rhy;
  rho= sqrt( rhs);

  if( rho <= 0.)
  {
    rhx=1.;
    rhy=0.;
    phx=0.;
    phy=1.;
  }
  else
  {
    rhx= rhx/ rho;
    rhy= rhy/ rho;
    phx=- rhy;
    phy= rhx;
  }

  cph= rhx* xsn+ rhy* ysn;
  sph= rhy* xsn- rhx* ysn;

  if( fabs( cph) < 1.0e-10)
    cph=0.;
  if( fabs( sph) < 1.0e-10)
    sph=0.;

  zph= zo+ zt;
  zphs= zph* zph;
  r2s= rhs+ zphs;
  r2= sqrt( r2s);
  rk= r2* TP;
  xx2= cmplx( cos( rk),- sin( rk));

  /* use norton approximation for field due to ground.  current is */
  /* lumped at segment center with current moment for constant, sine, */
  /* or cosine distribution. */
  if( isnor != 1)
  {
    zmh=1.;
    r1=1.;
    xx1=0.;
    gwave( &erv, &ezv, &erh, &ezh, &eph);

    et=-CONST1* frati* xx2/( r2s* r2);
    er=2.* et* cmplx(1.0, rk);
    et= et* cmplx(1.0 - rk* rk, rk);
    hrv=( er+ et)* rho* zph/ r2s;
    hzv=( zphs* er- rhs* et)/ r2s;
    hrh=( rhs* er- zphs* et)/ r2s;
    erv= erv- hrv;
    ezv= ezv- hzv;
    erh= erh+ hrh;
    ezh= ezh+ hrv;
    eph= eph+ et;
    erv= erv* salpj;
    ezv= ezv* salpj;
    erh= erh* sn* cph;
    ezh= ezh* sn* cph;
    eph= eph* sn* sph;
    erh= erv+ erh;
    e[0]=( erh* rhx+ eph* phx)* s;
    e[1]=( erh* rhy+ eph* phy)* s;
    e[2]=( ezv+ ezh)* s;
    e[3]=0.;
    e[4]=0.;
    e[5]=0.;
    sfac= PI* s;
    sfac= sin( sfac)/ sfac;
    e[6]= e[0]* sfac;
    e[7]= e[1]* sfac;
    e[8]= e[2]* sfac;

    return;
  } /* if( isnor != 1) */

  /* interpolate in sommerfeld field tables */
  if( rho >= 1.0e-12)
    thet= atan( zph/ rho);
  else
    thet= POT;

  /* combine vertical and horizontal components and convert */
  /* to x,y,z components. multiply by exp(-jkr)/r. */
  intrp( r2, thet, &erv, &ezv, &erh, &eph );
  xx2= xx2/ r2;
  sfac= sn* cph;
  erh= xx2*( salpj* erv+ sfac* erh);
  ezh= xx2*( salpj* ezv- sfac* erv);
  /* x,y,z fields for constant current */
  eph= sn* sph* xx2* eph;
  e[0]= erh* rhx+ eph* phx;
  e[1]= erh* rhy+ eph* phy;
  e[2]= ezh;
  /* x,y,z fields for sine current */
  rk= TP* t;
  sfac= sin( rk);
  e[3]= e[0]* sfac;
  e[4]= e[1]* sfac;
  /* x,y,z fields for cosine current */
  e[5]= e[2]* sfac;
  sfac= cos( rk);
  e[6]= e[0]* sfac;
  e[7]= e[1]* sfac;
  e[8]= e[2]* sfac;

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine to solve the matrix equation lu*x=b where l is a unit */
/* lower triangular matrix and u is an upper triangular matrix both */
/* of which are stored in a.  the rhs vector b is input and the */
/* solution is returned through vector b.   (matrix transposed. */
void solve( int n, complex double *a, int *ip,
    complex double *b, int ndim )
{
  int i, ip1, j, k, pia;
  complex double sum, *scm = NULL;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, np2m * sizeof(complex double) );

  /* forward substitution */
  for( i = 0; i < n; i++ )
  {
    pia= ip[i]-1;
    scm[i]= b[pia];
    b[pia]= b[i];
    ip1= i+1;

    if( ip1 < n)
      for( j = ip1; j < n; j++ )
	b[j] -= a[j+i*ndim]* scm[i];
  }

  /* backward substitution */
  for( k = 0; k < n; k++ )
  {
    i= n-k-1;
    sum=CPLX_00;
    ip1= i+1;

    if( ip1 < n)
      for( j = ip1; j < n; j++ )
	sum += a[i+j*ndim]* b[j];

    b[i]=( scm[i]- sum)/ a[i+i*ndim];
  }

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine solves, for symmetric structures, handles the */
/* transformation of the right hand side vector and solution */
/* of the matrix eq. */
void solves( complex double *a, int *ip, complex double *b,
    int neq, int nrh, int np, int n, int mp, int m)
{
  int  npeq, nrow, ic, i, kk, ia, ib, j, k;
  double fnop, fnorm;
  complex double  sum, *scm = NULL;

  npeq= np+ 2*mp;
  fnop= nop;
  fnorm=1./ fnop;
  nrow= neq;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, np2m * sizeof(complex double) );

  if( nop != 1)
  {
    for( ic = 0; ic < nrh; ic++ )
    {
      if( (n != 0) && (m != 0) )
      {
	for( i = 0; i < neq; i++ )
	  scm[i]= b[i+ic*neq];

	kk=2* mp;
	ia= np-1;
	ib= n-1;
	j= np-1;

	for( k = 0; k < nop; k++ )
	{
	  if( k != 0 )
	  {
	    for( i = 0; i < np; i++ )
	    {
	      ia++;
	      j++;
	      b[j+ic*neq]= scm[ia];
	    }

	    if( k == (nop-1) )
	      continue;

	  } /* if( k != 0 ) */

	  for( i = 0; i < kk; i++ )
	  {
	    ib++;
	    j++;
	    b[j+ic*neq]= scm[ib];
	  }

	} /* for( k = 0; k < nop; k++ ) */

      } /* if( (n != 0) && (m != 0) ) */

      /* transform matrix eq. rhs vector according to symmetry modes */
      for( i = 0; i < npeq; i++ )
      {
	for( k = 0; k < nop; k++ )
	{
	  ia= i+ k* npeq;
	  scm[k]= b[ia+ic*neq];
	}

	sum= scm[0];
	for( k = 1; k < nop; k++ )
	  sum += scm[k];

	b[i+ic*neq]= sum* fnorm;

	for( k = 1; k < nop; k++ )
	{
	  ia= i+ k* npeq;
	  sum= scm[0];

	  for( j = 1; j < nop; j++ )
	    sum += scm[j]* conj( ssx[k+j*nop]);

	  b[ia+ic*neq]= sum* fnorm;
	}

      } /* for( i = 0; i < npeq; i++ ) */

    } /* for( ic = 0; ic < nrh; ic++ ) */

  } /* if( nop != 1) */

  /* solve each mode equation */
  for( kk = 0; kk < nop; kk++ )
  {
    ia= kk* npeq;
    ib= ia;

    for( ic = 0; ic < nrh; ic++ )
      solve( npeq, &a[ib], &ip[ia], &b[ia+ic*neq], nrow );

  } /* for( kk = 0; kk < nop; kk++ ) */

  if( nop == 1)
  {
    free_ptr( (void *)&scm );
    return;
  }

  /* inverse transform the mode solutions */
  for( ic = 0; ic < nrh; ic++ )
  {
    for( i = 0; i < npeq; i++ )
    {
      for( k = 0; k < nop; k++ )
      {
	ia= i+ k* npeq;
	scm[k]= b[ia+ic*neq];
      }

      sum= scm[0];
      for( k = 1; k < nop; k++ )
	sum += scm[k];

      b[i+ic*neq]= sum;
      for( k = 1; k < nop; k++ )
      {
	ia= i+ k* npeq;
	sum= scm[0];

	for( j = 1; j < nop; j++ )
	  sum += scm[j]* ssx[k+j*nop];

	b[ia+ic*neq]= sum;
      }

    } /* for( i = 0; i < npeq; i++ ) */

    if( (n == 0) || (m == 0) )
      continue;

    for( i = 0; i < neq; i++ )
      scm[i]= b[i+ic*neq];

    kk=2* mp;
    ia= np-1;
    ib= n-1;
    j= np-1;

    for( k = 0; k < nop; k++ )
    {
      if( k != 0 )
      {
	for( i = 0; i < np; i++ )
	{
	  ia++;
	  j++;
	  b[ia+ic*neq]= scm[j];
	}

	if( k == nop)
	  continue;

      } /* if( k != 0 ) */

      for( i = 0; i < kk; i++ )
      {
	ib++;
	j++;
	b[ib+ic*neq]= scm[j];
      }

    } /* for( k = 0; k < nop; k++ ) */

  } /* for( ic = 0; ic < nrh; ic++ ) */

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* compute basis function i */
void tbf( int i, int icap )
{
  int ix, jcox, jcoxx, jend, iend, njun1=0, njun2, jsnop, jsnox;
  double pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;
  double d, sig; /*** also global ***/

  jsno=0;
  pp=0.;
  ix = i-1;
  jcox= icon1[ix];

  if( jcox > PCHCON)
    jcox= i;

  jend=-1;
  iend=-1;
  sig=-1.;

  do
  {
    if( jcox != 0 )
    {
      if( jcox < 0 )
	jcox=- jcox;
      else
      {
	sig=- sig;
	jend=- jend;
      }

      jcoxx = jcox-1;
      jsno++;
      jsnox = jsno-1;
      jco[jsnox]= jcox;
      d= PI* si[jcoxx];
      sdh= sin( d);
      cdh= cos( d);
      sd=2.* sdh* cdh;

      if( d <= 0.015)
      {
	omc=4.* d* d;
	omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
      }
      else
	omc=1.- cdh* cdh+ sdh* sdh;

      aj=1./( log(1./( PI* bi[jcoxx]))-.577215664);
      pp= pp- omc/ sd* aj;
      ax[jsnox]= aj/ sd* sig;
      bx[jsnox]= aj/(2.* cdh);
      cx[jsnox]=- aj/(2.* sdh)* sig;

      if( jcox != i)
      {
	if( jend == 1)
	  jcox= icon2[jcoxx];
	else
	  jcox= icon1[jcoxx];

	if( abs(jcox) != i )
	{
	  if( jcox != 0 )
	    continue;
	  else
	  {
	    fprintf( output_fp,
		"\n  TBF - SEGMENT CONNECTION ERROR FOR SEGMENT %5d", i );
	    stop(-1);
	  }
	}

      } /* if( jcox != i) */
      else
	bx[jsnox] =- bx[jsnox];

      if( iend == 1)
	break;

    } /* if( jcox != 0 ) */

    pm=- pp;
    pp=0.;
    njun1= jsno;

    jcox= icon2[ix];
    if( jcox > PCHCON)
      jcox= i;

    jend=1;
    iend=1;
    sig=-1.;

  } /* do */
  while( jcox != 0 );

  njun2= jsno- njun1;
  jsnop= jsno;
  jco[jsnop]= i;
  d= PI* si[ix];
  sdh= sin( d);
  cdh= cos( d);
  sd=2.* sdh* cdh;
  cd= cdh* cdh- sdh* sdh;

  if( d <= 0.015)
  {
    omc=4.* d* d;
    omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
  }
  else
    omc=1.- cd;

  ap=1./( log(1./( PI* bi[ix]))-.577215664);
  aj= ap;

  if( njun1 == 0)
  {
    if( njun2 == 0)
    {
      bx[jsnop]=0.;

      if( icap == 0)
	xxi=0.;
      else
      {
	qp= PI* bi[ix];
	xxi= qp* qp;
	xxi= qp*(1.-.5* xxi)/(1.- xxi);
      }

      cx[jsnop]=1./( cdh- xxi* sdh);
      jsno= jsnop+1;
      ax[jsnop]=-1.;
      return;

    } /* if( njun2 == 0) */

    if( icap == 0)
      xxi=0.;
    else
    {
      qp= PI* bi[ix];
      xxi= qp* qp;
      xxi= qp*(1.-.5* xxi)/(1.- xxi);
    }

    qp=-( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));
    d= cd- xxi* sd;
    bx[jsnop]=( sdh+ ap* qp*( cdh- xxi* sdh))/ d;
    cx[jsnop]=( cdh+ ap* qp*( sdh+ xxi* cdh))/ d;

    for( iend = 0; iend < njun2; iend++ )
    {
      ax[iend]=- ax[iend]* qp;
      bx[iend]= bx[iend]* qp;
      cx[iend]=- cx[iend]* qp;
    }

    jsno= jsnop+1;
    ax[jsnop]=-1.;
    return;

  } /* if( njun1 == 0) */

  if( njun2 == 0)
  {
    if( icap == 0)
      xxi=0.;
    else
    {
      qm= PI* bi[ix];
      xxi= qm* qm;
      xxi= qm*(1.-.5* xxi)/(1.- xxi);
    }

    qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));
    d= cd- xxi* sd;
    bx[jsnop]=( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
    cx[jsnop]=( cdh- aj* qm*( sdh+ xxi* cdh))/ d;

    for( iend = 0; iend < njun1; iend++ )
    {
      ax[iend]= ax[iend]* qm;
      bx[iend]= bx[iend]* qm;
      cx[iend]= cx[iend]* qm;
    }

    jsno= jsnop+1;
    ax[jsnop]=-1.;
    return;

  } /* if( njun2 == 0) */

  qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
  qm=( ap* omc- pp* sd)/ qp;
  qp=-( aj* omc+ pm* sd)/ qp;
  bx[jsnop]=( aj* qm+ ap* qp)* sdh/ sd;
  cx[jsnop]=( aj* qm- ap* qp)* cdh/ sd;

  for( iend = 0; iend < njun1; iend++ )
  {
    ax[iend]= ax[iend]* qm;
    bx[iend]= bx[iend]* qm;
    cx[iend]= cx[iend]* qm;
  }

  jend= njun1;
  for( iend = jend; iend < jsno; iend++ )
  {
    ax[iend]=- ax[iend]* qp;
    bx[iend]= bx[iend]* qp;
    cx[iend]=- cx[iend]* qp;
  }

  jsno= jsnop+1;
  ax[jsnop]=-1.;
}

/*-----------------------------------------------------------------------*/

/* test for convergence in numerical integration */
void test( double f1r, double f2r, double *tr,
    double f1i, double f2i, double *ti, double dmin )
{
  double den;

  den= fabs( f2r);
  *tr= fabs( f2i);

  if( den < *tr)
    den= *tr;
  if( den < dmin)
    den= dmin;

  if( den < 1.0e-37)
  {
    *tr=0.;
    *ti=0.;
    return;
  }

  *tr= fabs(( f1r- f2r)/ den);
  *ti= fabs(( f1i- f2i)/ den);

  return;
}

/*-----------------------------------------------------------------------*/

/* compute the components of all basis functions on segment j */
void trio( int j )
{
  int jcox, jcoxx, jsnox, jx, jend=0, iend=0;

  jsno=0;
  jx = j-1;
  jcox= icon1[jx];
  jcoxx = jcox-1;

  if( jcox <= PCHCON)
  {
    jend=-1;
    iend=-1;
  }

  if( (jcox == 0) || (jcox > PCHCON) )
  {
    jcox= icon2[jx];
    jcoxx = jcox-1;

    if( jcox <= PCHCON)
    {
      jend=1;
      iend=1;
    }

    if( jcox == 0 || (jcox > PCHCON) )
    {
      jsnox = jsno;
      jsno++;

      /* Allocate to connections buffers */
      if( jsno >= maxcon )
      {
	maxcon = jsno +1;
	mem_realloc( (void *)&jco, maxcon * sizeof(int) );
	mem_realloc( (void *) &ax, maxcon * sizeof(double) );
	mem_realloc( (void *) &bx, maxcon * sizeof(double) );
	mem_realloc( (void *) &cx, maxcon * sizeof(double) );
      }

      sbf( j, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
      jco[jsnox]= j;
      return;
    }

  } /* if( (jcox == 0) || (jcox > PCHCON) ) */

  do
  {
    if( jcox < 0 )
      jcox=- jcox;
    else
      jend=- jend;
    jcoxx = jcox-1;

    if( jcox != j)
    {
      jsnox = jsno;
      jsno++;

      /* Allocate to connections buffers */
      if( jsno >= maxcon )
      {
	maxcon = jsno +1;
	mem_realloc( (void *)&jco, maxcon * sizeof(int) );
	mem_realloc( (void *) &ax, maxcon * sizeof(double) );
	mem_realloc( (void *) &bx, maxcon * sizeof(double) );
	mem_realloc( (void *) &cx, maxcon * sizeof(double) );
      }

      sbf( jcox, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
      jco[jsnox]= jcox;

      if( jend != 1)
	jcox= icon1[jcoxx];
      else
	jcox= icon2[jcoxx];

      if( jcox == 0 )
      {
	fprintf( output_fp,
	    "\n  TRIO - SEGMENT CONNENTION ERROR FOR SEGMENT %5d", j );
	stop(-1);
      }
      else
	continue;

    } /* if( jcox != j) */

    if( iend == 1)
      break;

    jcox= icon2[jx];

    if( jcox > PCHCON)
      break;

    jend=1;
    iend=1;

  } /* do */
  while( jcox != 0 );

  jsnox = jsno;
  jsno++;

  /* Allocate to connections buffers */
  if( jsno >= maxcon )
  {
    maxcon = jsno +1;
    mem_realloc( (void *)&jco, maxcon * sizeof(int) );
    mem_realloc( (void *) &ax, maxcon * sizeof(double) );
    mem_realloc( (void *) &bx, maxcon * sizeof(double) );
    mem_realloc( (void *) &cx, maxcon * sizeof(double) );
  }

  sbf( j, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
  jco[jsnox]= j;

  return;

}

/*-----------------------------------------------------------------------*/

/* calculates the electric field due to unit current */
/* in the t1 and t2 directions on a patch */
void unere( double xob, double yob, double zob )
{
  double zr, t1zr, t2zr, rx, ry, rz, r, tt1;
  double tt2, rt, xymag, px, py, cth, r2;
  complex double er, q1, q2, rrv, rrh, edp;

  zr= zj;
  t1zr= t1zj;
  t2zr= t2zj;

  if( ipgnd == 2)
  {
    zr=- zr;
    t1zr=- t1zr;
    t2zr=- t2zr;
  }

  rx= xob- xj;
  ry= yob- yj;
  rz= zob- zr;
  r2= rx* rx+ ry* ry+ rz* rz;

  if( r2 <= 1.0e-20)
  {
    exk=CPLX_00;
    eyk=CPLX_00;
    ezk=CPLX_00;
    exs=CPLX_00;
    eys=CPLX_00;
    ezs=CPLX_00;
    return;
  }

  r= sqrt( r2);
  tt1=- TP* r;
  tt2= tt1* tt1;
  rt= r2* r;
  er= cmplx( sin( tt1),- cos( tt1))*( CONST2* s);
  q1= cmplx( tt2-1., tt1)* er/ rt;
  q2= cmplx(3.- tt2,-3.* tt1)* er/( rt* r2);
  er = q2*( t1xj* rx+ t1yj* ry+ t1zr* rz);
  exk= q1* t1xj+ er* rx;
  eyk= q1* t1yj+ er* ry;
  ezk= q1* t1zr+ er* rz;
  er= q2*( t2xj* rx+ t2yj* ry+ t2zr* rz);
  exs= q1* t2xj+ er* rx;
  eys= q1* t2yj+ er* ry;
  ezs= q1* t2zr+ er* rz;

  if( ipgnd == 1)
    return;

  if( iperf == 1)
  {
    exk=- exk;
    eyk=- eyk;
    ezk=- ezk;
    exs=- exs;
    eys=- eys;
    ezs=- ezs;
    return;
  }

  xymag= sqrt( rx* rx+ ry* ry);
  if( xymag <= 1.0e-6)
  {
    px=0.;
    py=0.;
    cth=1.;
    rrv=CPLX_10;
  }
  else
  {
    px=- ry/ xymag;
    py= rx/ xymag;
    cth= rz/ sqrt( xymag* xymag+ rz* rz);
    rrv= csqrt(1.- zrati* zrati*(1.- cth* cth));
  }

  rrh= zrati* cth;
  rrh=( rrh- rrv)/( rrh+ rrv);
  rrv= zrati* rrv;
  rrv=-( cth- rrv)/( cth+ rrv);
  edp=( exk* px+ eyk* py)*( rrh- rrv);
  exk= exk* rrv+ edp* px;
  eyk= eyk* rrv+ edp* py;
  ezk= ezk* rrv;
  edp=( exs* px+ eys* py)*( rrh- rrv);
  exs= exs* rrv+ edp* px;
  eys= eys* rrv+ edp* py;
  ezs= ezs* rrv;

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine wire generates segment geometry */
/* data for a straight wire of ns segments. */
void wire( double xw1, double yw1, double zw1,
    double xw2, double yw2, double zw2, double rad,
    double rdel, double rrad, int ns, int itg )
{
  int ist, i, mreq;
  double xd, yd, zd, delz, rd, fns, radz;
  double xs1, ys1, zs1, xs2, ys2, zs2;

  ist= n;
  n= n+ ns;
  np= n;
  mp= m;
  ipsym=0;

  if( ns < 1)
    return;

  /* Reallocate tags buffer */
  mem_realloc( (void *)&itag, (n+m) * sizeof(int) );/*????*/

  /* Reallocate wire buffers */
  mreq = n * sizeof(double);
  mem_realloc( (void *)&x, mreq );
  mem_realloc( (void *)&y, mreq );
  mem_realloc( (void *)&z, mreq );
  mem_realloc( (void *)&x2, mreq );
  mem_realloc( (void *)&y2, mreq );
  mem_realloc( (void *)&z2, mreq );
  mem_realloc( (void *)&bi, mreq );

  xd= xw2- xw1;
  yd= yw2- yw1;
  zd= zw2- zw1;

  if( fabs( rdel-1.) >= 1.0e-6)
  {
    delz= sqrt( xd* xd+ yd* yd+ zd* zd);
    xd= xd/ delz;
    yd= yd/ delz;
    zd= zd/ delz;
    delz= delz*(1.- rdel)/(1.- pow(rdel, ns) );
    rd= rdel;
  }
  else
  {
    fns= ns;
    xd= xd/ fns;
    yd= yd/ fns;
    zd= zd/ fns;
    delz=1.;
    rd=1.;
  }

  radz= rad;
  xs1= xw1;
  ys1= yw1;
  zs1= zw1;

  for( i = ist; i < n; i++ )
  {
    itag[i]= itg;
    xs2= xs1+ xd* delz;
    ys2= ys1+ yd* delz;
    zs2= zs1+ zd* delz;
    x[i]= xs1;
    y[i]= ys1;
    z[i]= zs1;
    x2[i]= xs2;
    y2[i]= ys2;
    z2[i]= zs2;
    bi[i]= radz;
    delz= delz* rd;
    radz= radz* rrad;
    xs1= xs2;
    ys1= ys2;
    zs1= zs2;
  }

  x2[n-1]= xw2;
  y2[n-1]= yw2;
  z2[n-1]= zw2;

  return;
}

/*-----------------------------------------------------------------------*/

/* zint computes the internal impedance of a circular wire */
complex double zint( double sigl, double rolam )
{
#define cc1	( 6.0e-7     + 1.9e-6fj)
#define cc2	(-3.4e-6     + 5.1e-6fj)
#define cc3	(-2.52e-5    + 0.fj)
#define cc4	(-9.06e-5    - 9.01e-5fj)
#define cc5	( 0.         - 9.765e-4fj)
#define cc6	(.0110486    - .0110485fj)
#define cc7	( 0.         - .3926991fj)
#define cc8	( 1.6e-6     - 3.2e-6fj)
#define cc9	( 1.17e-5    - 2.4e-6fj)
#define cc10	( 3.46e-5    + 3.38e-5fj)
#define cc11	( 5.0e-7     + 2.452e-4fj)
#define cc12	(-1.3813e-3  + 1.3811e-3fj)
#define cc13	(-6.25001e-2 - 1.0e-7fj)
#define cc14	(.7071068    + .7071068fj)
#define cn	cc14

#define th(d) ( (((((cc1*(d)+cc2)*(d)+cc3)*(d)+cc4)*(d)+cc5)*(d)+cc6)*(d) + cc7 )
#define ph(d) ( (((((cc8*(d)+cc9)*(d)+cc10)*(d)+cc11)*(d)+cc12)*(d)+cc13)*(d)+cc14 )
#define f(d)  ( csqrt(POT/(d))*cexp(-cn*(d)+th(-8./x)) )
#define g(d)  ( cexp(cn*(d)+th(8./x))/csqrt(TP*(d)) )

  double x, y, s, ber, bei;
  double tpcmu = 2.368705e+3;
  double cmotp = 60.00;
  complex double zint, br1, br2;

  x= sqrt( tpcmu* sigl)* rolam;
  if( x <= 110.)
  {
    if( x <= 8.)
    {
      y= x/8.;
      y= y* y;
      s= y* y;

      ber=((((((-9.01e-6* s+1.22552e-3)* s-.08349609)* s+ 2.6419140)*
	      s-32.363456)* s+113.77778)* s-64.)* s+1.;

      bei=((((((1.1346e-4* s-.01103667)* s+.52185615)* s-10.567658)*
	      s+72.817777)* s-113.77778)* s+16.)* y;

      br1= cmplx( ber, bei);

      ber=(((((((-3.94e-6* s+4.5957e-4)* s-.02609253)* s+ .66047849)*
		s-6.0681481)* s+14.222222)* s-4.)* y)* x;

      bei=((((((4.609e-5* s-3.79386e-3)* s+.14677204)* s- 2.3116751)*
	      s+11.377778)* s-10.666667)* s+.5)* x;

      br2= cmplx( ber, bei);
      br1= br1/ br2;
      zint= CPLX_01* sqrt( cmotp/sigl )* br1/ rolam;

      return( zint );

    } /* if( x <= 8.) */

    br2= CPLX_01* f(x)/ PI;
    br1= g( x)+ br2;
    br2= g( x)* ph(8./ x)- br2* ph(-8./ x);
    br1= br1/ br2;
    zint= CPLX_01* sqrt( cmotp/ sigl)* br1/ rolam;

    return( zint );

  } /* if( x <= 110.) */

  br1= cmplx(.70710678,-.70710678);
  zint= CPLX_01* sqrt( cmotp/ sigl)* br1/ rolam;

  return( zint );
}

/*-----------------------------------------------------------------------*/

/* returns smallest of two arguments */
int min( int a, int b )
{
  if( a < b )
    return(a);
  else
    return(b);
}

/*-----------------------------------------------------------------------*/

static void sig_handler( int signal )
{
  switch( signal )
  {
    case SIGINT :
      fprintf( stderr, "\n%s\n", "nec2c: exiting via user interrupt" );
      exit( signal );

    case SIGSEGV :
      fprintf( stderr, "\n%s\n", "nec2c: segmentation fault" );
      exit( signal );

    case SIGFPE :
      fprintf( stderr, "\n%s\n", "nec2c: floating point exception" );
      exit( signal );

    case SIGABRT :
      fprintf( stderr, "\n%s\n", "nec2c: abort signal received" );
      exit( signal );

    case SIGTERM :
      fprintf( stderr, "\n%s\n", "nec2c: termination request received" );

      stop( signal );
  }

} /* end of sig_handler() */

/*------------------------------------------------------------------------*/

