/*** Translated to the C language by N. Kyriazis  20 Aug 2003 ***

  Program NEC(input,tape5=input,output,tape11,tape12,tape13,tape14,
  tape15,tape16,tape20,tape21)

  Numerical Electromagnetics Code (NEC2)  developed at Lawrence
  Livermore lab., Livermore, CA.  (contact G. Burke at 415-422-8414
  for problems with the NEC code. For problems with the vax implem-
  entation, contact J. Breakall at 415-422-8196 or E. Domning at 415
  422-5936)
  file created 4/11/80.

				***********Notice**********
 This computer code material was prepared as an account of work
 sponsored by the United States government.  Neither the United
 States nor the United States Department Of Energy, nor any of
 their employees, nor any of their contractors, subcontractors,
 or their employees, makes any warranty, express or implied, or
 assumes any legal liability or responsibility for the accuracy,
 completeness or usefulness of any information, apparatus, product
 or process disclosed, or represents that its use would not infringe
 privately-owned rights.

 ******************************************************************/

#include "nec2c.h"

/* common  /tmi/ */
extern tmi_t tmi;

/*common  /ggrid/ */
extern ggrid_t ggrid;

/* common  /data/ */
extern data_t data;

/* common  /crnt/ */
extern crnt_t crnt;

/* common  /vsorc/ */
extern vsorc_t vsorc;

/* common  /segj/ */
extern segj_t segj;

/* common  /yparm/ */
extern yparm_t yparm;

/* common  /zload/ */
extern zload_t zload;

/* common  /smat/ */
extern smat_t smat;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/*-----------------------------------------------------------------------*/

/* cabc computes coefficients of the constant (a), sine (b), and */
/* cosine (c) terms in the current interpolation functions for the */
/* current vector cur. */
void cabc( complex long double *curx)
{
  int i, is, j, jx, jco1, jco2;
  long double ar, ai, sh;
  complex long double curd, cs1, cs2;

  if( data.n != 0)
  {
	for( i = 0; i < data.n; i++ )
	{
	  crnt.air[i]=0.;
	  crnt.aii[i]=0.;
	  crnt.bir[i]=0.;
	  crnt.bii[i]=0.;
	  crnt.cir[i]=0.;
	  crnt.cii[i]=0.;
	}

	for( i = 0; i < data.n; i++ )
	{
	  ar= creall( curx[i]);
	  ai= cimagl( curx[i]);
	  tbf( i+1, 1 );

	  for( jx = 0; jx < segj.jsno; jx++ )
	  {
		j= segj.jco[jx]-1;
		crnt.air[j] += segj.ax[jx]* ar;
		crnt.aii[j] += segj.ax[jx]* ai;
		crnt.bir[j] += segj.bx[jx]* ar;
		crnt.bii[j] += segj.bx[jx]* ai;
		crnt.cir[j] += segj.cx[jx]* ar;
		crnt.cii[j] += segj.cx[jx]* ai;
	  }

	} /* for( i = 0; i < n; i++ ) */

	if( vsorc.nqds != 0)
	{
	  for( is = 0; is < vsorc.nqds; is++ )
	  {
		i= vsorc.iqds[is]-1;
		jx= data.icon1[i];
		data.icon1[i]=0;
		tbf(i+1,0);
		data.icon1[i]= jx;
		sh= data.si[i]*.5;
		curd= CCJ* vsorc.vqds[is]/( (logl(2.* sh/ data.bi[i])-1.)*
			(segj.bx[segj.jsno-1]* cosl(TP* sh)+ segj.cx[segj.jsno-1]*
			 sinl(TP* sh))* data.wlam );
		ar= creall( curd);
		ai= cimagl( curd);

		for( jx = 0; jx < segj.jsno; jx++ )
		{
		  j= segj.jco[jx]-1;
		  crnt.air[j]= crnt.air[j]+ segj.ax[jx]* ar;
		  crnt.aii[j]= crnt.aii[j]+ segj.ax[jx]* ai;
		  crnt.bir[j]= crnt.bir[j]+ segj.bx[jx]* ar;
		  crnt.bii[j]= crnt.bii[j]+ segj.bx[jx]* ai;
		  crnt.cir[j]= crnt.cir[j]+ segj.cx[jx]* ar;
		  crnt.cii[j]= crnt.cii[j]+ segj.cx[jx]* ai;
		}

	  } /* for( is = 0; is < vsorc.nqds; is++ ) */

	} /* if( vsorc.nqds != 0) */

	for( i = 0; i < data.n; i++ )
	  curx[i]= cmplx( crnt.air[i]+crnt.cir[i], crnt.aii[i]+crnt.cii[i] );

  } /* if( n != 0) */

  if( data.m == 0)
	return;

  /* convert surface currents from */
  /* t1,t2 components to x,y,z components */
  jco1= data.np2m;
  jco2= jco1+ data.m;
  for( i = 1; i <= data.m; i++ )
  {
	jco1 -= 2;
	jco2 -= 3;
	cs1= curx[jco1];
	cs2= curx[jco1+1];
	curx[jco2]  = cs1* data.t1x[data.m-i]+ cs2* data.t2x[data.m-i];
	curx[jco2+1]= cs1* data.t1y[data.m-i]+ cs2* data.t2y[data.m-i];
	curx[jco2+2]= cs1* data.t1z[data.m-i]+ cs2* data.t2z[data.m-i];
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* couple computes the maximum coupling between pairs of segments. */
void couple( complex long double *cur, long double wlam )
{
  int j, j1, j2, l1, i, k, itt1, itt2, its1, its2, isg1, isg2, npm1;
  long double dbc, c, gmax;
  complex long double y11, y12, y22, yl, yin, zl, zin, rho;

  if( (vsorc.nsant != 1) || (vsorc.nvqd != 0) )
	return;

  j= isegno( yparm.nctag[yparm.icoup], yparm.ncseg[yparm.icoup]);
  if( j != vsorc.isant[0] )
	return;

  zin= vsorc.vsant[0];
  yparm.icoup++;
  mem_realloc( (void *)&yparm.y11a, yparm.icoup * sizeof( complex long double) );
  yparm.y11a[yparm.icoup-1]= cur[j-1]*wlam/zin;

  l1=(yparm.icoup-1)*(yparm.ncoup-1);
  for( i = 0; i < yparm.ncoup; i++ )
  {
	if( (i+1) == yparm.icoup)
	  continue;

	l1++;
	mem_realloc( (void *)&yparm.y12a, l1 * sizeof( complex long double) );
	k= isegno( yparm.nctag[i], yparm.ncseg[i]);
	yparm.y12a[l1-1]= cur[k-1]* wlam/ zin;
  }

  if( yparm.icoup < yparm.ncoup)
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

  npm1= yparm.ncoup-1;

  for( i = 0; i < npm1; i++ )
  {
	itt1= yparm.nctag[i];
	its1= yparm.ncseg[i];
	isg1= isegno( itt1, its1);
	l1= i+1;

	for( j = l1; j < yparm.ncoup; j++ )
	{
	  itt2= yparm.nctag[j];
	  its2= yparm.ncseg[j];
	  isg2= isegno( itt2, its2);
	  j1= j+ i* npm1-1;
	  j2= i+ j* npm1;
	  y11= yparm.y11a[i];
	  y22= yparm.y11a[j];
	  y12=.5*( yparm.y12a[j1]+ yparm.y12a[j2]);
	  yin= y12* y12;
	  dbc= cabsl( yin);
	  c= dbc/(2.* creall( y11)* creall( y22)- creall( yin));

	  if( (c >= 0.0) && (c <= 1.0) )
	  {
		if( c >= .01 )
		  gmax=(1.- sqrtl(1.- c*c))/c;
		else
		  gmax=.5*( c+.25* c* c* c);

		rho= gmax* conjl( yin)/ dbc;
		yl=((1.- rho)/(1.+ rho)+1.)* creall( y22)- y22;
		zl=1./ yl;
		yin= y11- yin/( y22+ yl);
		zin=1./ yin;
		dbc= db10( gmax);

		fprintf( output_fp, "\n"
			" %4d %4d %5d  %4d %4d %5d  %9.3LF"
			"  %12.5LE %12.5LE  %12.5LE %12.5LE",
			itt1, its1, isg1, itt2, its2, isg2, dbc,
			creall(zl), cimagl(zl), creall(zin), cimagl(zin) );

		continue;

	  } /* if( (c >= 0.0) && (c <= 1.0) ) */

	  fprintf( output_fp, "\n"
		  " %4d %4d %5d   %4d %4d %5d  **ERROR** "
		  "COUPLING IS NOT BETWEEN 0 AND 1. (= %12.5LE)",
		  itt1, its1, isg1, itt2, its2, isg2, c );

	} /* for( j = l1; j < yparm.ncoup; j++ ) */

  } /* for( i = 0; i < npm1; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* load calculates the impedance of specified */
/* segments for various types of loading */
void load( int *ldtyp, int *ldtag, int *ldtagf, int *ldtagt,
	long double *zlr, long double *zli, long double *zlc )
{
  int i, iwarn, istep, istepx, l1, l2, ldtags, jump, ichk;
  complex long double zt=CPLX_00, tpcj;

  tpcj = (0.0+1.883698955e+9fj);
  fprintf( output_fp, "\n"
	  "  LOCATION        RESISTANCE  INDUCTANCE  CAPACITANCE   "
	  "  IMPEDANCE (OHMS)   CONDUCTIVITY  CIRCUIT\n"
	  "  ITAG FROM THRU     OHMS       HENRYS      FARADS     "
	  "  REAL     IMAGINARY   MHOS/METER      TYPE" );

  /* initialize d array, used for temporary */
  /* storage of loading information. */
  mem_realloc( (void *)&zload.zarray, data.npm * sizeof(complex long double) );
  for( i = 0; i < data.n; i++ )
	zload.zarray[i]=CPLX_00;

  iwarn=FALSE;
  istep=0;

  /* cycle over loading cards */
  while( TRUE )
  {
	istepx = istep;
	istep++;

	if( istep > zload.nload)
	{
	  if( iwarn == TRUE )
		fprintf( output_fp,
			"\n  NOTE, SOME OF THE ABOVE SEGMENTS "
			"HAVE BEEN LOADED TWICE - IMPEDANCES ADDED" );

	  smat.nop = data.n/data.np;
	  if( smat.nop == 1)
		return;

	  for( i = 0; i < data.np; i++ )
	  {
		zt= zload.zarray[i];
		l1= i;

		for( l2 = 1; l2 < smat.nop; l2++ )
		{
		  l1 += data.np;
		  zload.zarray[l1]= zt;
		}
	  }
	  return;

	} /* if( istep > zload.nload) */

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
	l2= data.n;

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
		if( ldtags != data.itag[i])
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
		  zt= zlr[istepx]/ data.si[i]+ tpcj* zli[istepx]/( data.si[i]* data.wlam);
		  if( fabsl( zlc[istepx]) > 1.0e-20)
			zt += data.wlam/( tpcj* data.si[i]* zlc[istepx]);
		  break;

		case 2:
		  zt= tpcj* data.si[i]* zlc[istepx]/ data.wlam;
		  if( fabsl( zli[istepx]) > 1.0e-20)
			zt += data.si[i]* data.wlam/( tpcj* zli[istepx]);
		  if( fabsl( zlr[istepx]) > 1.0e-20)
			zt += data.si[i]/ zlr[istepx];
		  zt=1./ zt;
		  break;

		case 3:
		  zt= zlr[istepx]* data.wlam+ tpcj* zli[istepx];
		  if( fabsl( zlc[istepx]) > 1.0e-20)
			zt += 1./( tpcj* data.si[i]* data.si[i]* zlc[istepx]);
		  break;

		case 4:
		  zt= tpcj* data.si[i]* data.si[i]* zlc[istepx];
		  if( fabsl( zli[istepx]) > 1.0e-20)
			zt += 1./( tpcj* zli[istepx]);
		  if( fabsl( zlr[istepx]) > 1.0e-20)
			zt += 1./( zlr[istepx]* data.wlam);
		  zt=1./ zt;
		  break;

		case 5:
		  zt= cmplx( zlr[istepx], zli[istepx])/ data.si[i];
		  break;

		case 6:
		  zint( zlr[istepx]* data.wlam, data.bi[i], &zt );

	  } /* switch( jump ) */

	  if(( fabsl( creall( zload.zarray[i]))+ fabsl( cimagl( zload.zarray[i]))) > 1.0e-20)
		iwarn=TRUE;
	  zload.zarray[i] += zt;

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

/* gf computes the integrand exp(jkr)/(kr) for numerical integration. */
void gf( long double zk, long double *co, long double *si )
{
  long double zdk, rk, rks;

  zdk= zk- tmi.zpk;
  rk= sqrtl( tmi.rkb2+ zdk* zdk);
  *si= sinl( rk)/ rk;

  if( tmi.ij != 0 )
  {
	*co= cosl( rk)/ rk;
	return;
  }

  if( rk >= .2)
  {
	*co=( cosl( rk)-1.)/ rk;
	return;
  }

  rks= rk* rk;
  *co=((-1.38888889e-3* rks+4.16666667e-2)* rks-.5)* rk;

  return;
}

/*-----------------------------------------------------------------------*/

/* function db10 returns db for magnitude (field) */
long double db10( long double x )
{
  if( x < 1.e-20 )
	return( -999.99 );

  return( 10. * log10l(x) );
}

/*-----------------------------------------------------------------------*/

/* function db20 returns db for mag**2 (power) i */
long double db20( long double x )
{
  if( x < 1.e-20 )
	return( -999.99 );

  return( 20. * log10l(x) );
}

/*-----------------------------------------------------------------------*/

/* intrp uses bivariate cubic interpolation to obtain */
/* the values of 4 functions at the point (x,y). */
void intrp( long double x, long double y, complex long double *f1,
	complex long double *f2, complex long double *f3, complex long double *f4 )
{
  static int ix, iy, ixs=-10, iys=-10, igrs=-10, ixeg=0, iyeg=0;
  static int nxm2, nym2, nxms, nyms, nd, ndp;
  int nda[3]={11,17,9}, ndpa[3]={110,85,72};
  int igr, iadd, iadz, i, k, jump;
  static long double dx = 1., dy = 1., xs = 0., ys = 0., xz, yz;
  long double xx, yy;
  static complex long double a[4][4], b[4][4], c[4][4], d[4][4];
  complex long double p1=CPLX_00, p2=CPLX_00, p3=CPLX_00, p4=CPLX_00;
  complex long double fx1, fx2, fx3, fx4;

  jump = FALSE;
  if( (x < xs) || (y < ys) )
	jump = TRUE;
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
	  jump )
  {
	/* determine correct grid and grid region */
	if( x <= ggrid.xsa[1])
	  igr=0;
	else
	{
	  if( y > ggrid.ysa[2])
		igr=2;
	  else
		igr=1;
	}

	if( igr != igrs)
	{
	  igrs= igr;
	  dx= ggrid.dxa[igrs];
	  dy= ggrid.dya[igrs];
	  xs= ggrid.xsa[igrs];
	  ys= ggrid.ysa[igrs];
	  nxm2= ggrid.nxa[igrs]-2;
	  nym2= ggrid.nya[igrs]-2;
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
			p1= ggrid.ar1[iadd-2];
			p2= ggrid.ar1[iadd-1];
			p3= ggrid.ar1[iadd];
			p4= ggrid.ar1[iadd+1];
			break;

		  case 1:
			p1= ggrid.ar2[iadd-2];
			p2= ggrid.ar2[iadd-1];
			p3= ggrid.ar2[iadd];
			p4= ggrid.ar2[iadd+1];
			break;

		  case 2:
			p1= ggrid.ar3[iadd-2];
			p2= ggrid.ar3[iadd-1];
			p3= ggrid.ar3[iadd];
			p4= ggrid.ar3[iadd+1];

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
void intx( long double el1, long double el2, long double b,
	int ij, long double *sgr, long double *sgi)
{
  int ns, nt;
  int nx = 1, nma = 65536, nts = 4;
  int flag = TRUE;
  long double z, s, ze, fnm, ep, zend, fns, dz=0., zp, dzot=0., t00r, g1r, g5r, t00i;
  long double g1i, g5i, t01r, g3r, t01i, g3i, t10r, t10i, te1i, te1r, t02r;
  long double g2r, g4r, t02i, g2i, g4i, t11r, t11i, t20r, t20i, te2i, te2r;
  long double rx = 1.0e-4;

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
		if( fabsl(dz) <= ep)
		{
		  /* add contribution of near singularity for diagonal term */
		  if(ij == 0)
		  {
			*sgr=2.*( *sgr+ logl(( sqrtl( b* b+ s* s)+ s)/ b));
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
		  *sgr=2.*( *sgr+ logl(( sqrtl( b* b+ s* s)+ s)/ b));
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
		fprintf( output_fp, "\n  STEP SIZE LIMITED AT Z= %10.5LF", z );
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
		*sgr=2.*( *sgr+ logl(( sqrtl( b* b+ s* s)+ s)/ b));
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

/* returns smallest of two arguments */
int min( int a, int b )
{
  if( a < b )
	return(a);
  else
	return(b);
}

/*-----------------------------------------------------------------------*/

/* test for convergence in numerical integration */
void test( long double f1r, long double f2r, long double *tr,
	long double f1i, long double f2i, long double *ti, long double dmin )
{
  long double den;

  den= fabsl( f2r);
  *tr= fabsl( f2i);

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

  *tr= fabsl(( f1r- f2r)/ den);
  *ti= fabsl(( f1i- f2i)/ den);

  return;
}

/*-----------------------------------------------------------------------*/

/* compute component of basis function i on segment is. */
void sbf( int i, int is, long double *aa, long double *bb, long double *cc )
{
  int ix, jsno, june, jcox, jcoxx, jend, iend, njun1=0, njun2;
  long double d, sig, pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;

  *aa=0.;
  *bb=0.;
  *cc=0.;
  june=0;
  jsno=0;
  pp=0.;
  ix=i-1;

  jcox= data.icon1[ix];
  if( jcox > PCHCON) jcox= i;
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
	  d= PI* data.si[jcoxx];
	  sdh= sinl( d);
	  cdh= cosl( d);
	  sd=2.* sdh* cdh;

	  if( d <= 0.015)
	  {
		omc=4.* d* d;
		omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
	  }
	  else
		omc=1.- cdh* cdh+ sdh* sdh;

	  aj=1./( logl(1./( PI* data.bi[jcoxx]))-.577215664);
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
		  jcox= data.icon1[jcoxx];
		else
		  jcox= data.icon2[jcoxx];

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

	jcox= data.icon2[ix];
	if( jcox > PCHCON) jcox= i;

	jend=1;
	iend=1;
	sig=-1.;

  } /* do */
  while( jcox != 0 );

  njun2= jsno- njun1;
  d= PI* data.si[ix];
  sdh= sinl( d);
  cdh= cosl( d);
  sd=2.* sdh* cdh;
  cd= cdh* cdh- sdh* sdh;

  if( d <= 0.015)
  {
	omc=4.* d* d;
	omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
  }
  else
	omc=1.- cd;

  ap=1./( logl(1./( PI* data.bi[ix])) -.577215664);
  aj= ap;

  if( njun1 == 0)
  {
	if( njun2 == 0)
	{
	  *aa =-1.;
	  qp= PI* data.bi[ix];
	  xxi= qp* qp;
	  xxi= qp*(1.-.5* xxi)/(1.- xxi);
	  *cc=1./( cdh- xxi* sdh);
	  return;
	}

	qp= PI* data.bi[ix];
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
	qm= PI* data.bi[ix];
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

/* compute basis function i */
void tbf( int i, int icap )
{
  int ix, jcox, jcoxx, jend, iend, njun1=0, njun2, jsnop, jsnox;
  long double pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;
  long double d, sig; /*** also global ***/

  segj.jsno=0;
  pp=0.;
  ix = i-1;
  jcox= data.icon1[ix];

  if( jcox > PCHCON) jcox= i;

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
	  segj.jsno++;
	  jsnox = segj.jsno-1;
	  segj.jco[jsnox]= jcox;
	  d= PI* data.si[jcoxx];
	  sdh= sinl( d);
	  cdh= cosl( d);
	  sd=2.* sdh* cdh;

	  if( d <= 0.015)
	  {
		omc=4.* d* d;
		omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
	  }
	  else
		omc=1.- cdh* cdh+ sdh* sdh;

	  aj=1./( logl(1./( PI* data.bi[jcoxx]))-.577215664);
	  pp= pp- omc/ sd* aj;
	  segj.ax[jsnox]= aj/ sd* sig;
	  segj.bx[jsnox]= aj/(2.* cdh);
	  segj.cx[jsnox]=- aj/(2.* sdh)* sig;

	  if( jcox != i)
	  {
		if( jend == 1)
		  jcox= data.icon2[jcoxx];
		else
		  jcox= data.icon1[jcoxx];

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
		segj.bx[jsnox] =- segj.bx[jsnox];

	  if( iend == 1)
		break;

	} /* if( jcox != 0 ) */

	pm=- pp;
	pp=0.;
	njun1= segj.jsno;

	jcox= data.icon2[ix];
	if( jcox > PCHCON) jcox= i;

	jend=1;
	iend=1;
	sig=-1.;

  } /* do */
  while( jcox != 0 );

  njun2= segj.jsno- njun1;
  jsnop= segj.jsno;
  segj.jco[jsnop]= i;
  d= PI* data.si[ix];
  sdh= sinl( d);
  cdh= cosl( d);
  sd=2.* sdh* cdh;
  cd= cdh* cdh- sdh* sdh;

  if( d <= 0.015)
  {
	omc=4.* d* d;
	omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
  }
  else
	omc=1.- cd;

  ap=1./( logl(1./( PI* data.bi[ix]))-.577215664);
  aj= ap;

  if( njun1 == 0)
  {
	if( njun2 == 0)
	{
	  segj.bx[jsnop]=0.;

	  if( icap == 0)
		xxi=0.;
	  else
	  {
		qp= PI* data.bi[ix];
		xxi= qp* qp;
		xxi= qp*(1.-.5* xxi)/(1.- xxi);
	  }

	  segj.cx[jsnop]=1./( cdh- xxi* sdh);
	  segj.jsno= jsnop+1;
	  segj.ax[jsnop]=-1.;
	  return;

	} /* if( njun2 == 0) */

	if( icap == 0)
	  xxi=0.;
	else
	{
	  qp= PI* data.bi[ix];
	  xxi= qp* qp;
	  xxi= qp*(1.-.5* xxi)/(1.- xxi);
	}

	qp=-( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));
	d= cd- xxi* sd;
	segj.bx[jsnop]=( sdh+ ap* qp*( cdh- xxi* sdh))/ d;
	segj.cx[jsnop]=( cdh+ ap* qp*( sdh+ xxi* cdh))/ d;

	for( iend = 0; iend < njun2; iend++ )
	{
	  segj.ax[iend]=- segj.ax[iend]* qp;
	  segj.bx[iend]= segj.bx[iend]* qp;
	  segj.cx[iend]=- segj.cx[iend]* qp;
	}

	segj.jsno= jsnop+1;
	segj.ax[jsnop]=-1.;
	return;

  } /* if( njun1 == 0) */

  if( njun2 == 0)
  {
	if( icap == 0)
	  xxi=0.;
	else
	{
	  qm= PI* data.bi[ix];
	  xxi= qm* qm;
	  xxi= qm*(1.-.5* xxi)/(1.- xxi);
	}

	qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));
	d= cd- xxi* sd;
	segj.bx[jsnop]=( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
	segj.cx[jsnop]=( cdh- aj* qm*( sdh+ xxi* cdh))/ d;

	for( iend = 0; iend < njun1; iend++ )
	{
	  segj.ax[iend]= segj.ax[iend]* qm;
	  segj.bx[iend]= segj.bx[iend]* qm;
	  segj.cx[iend]= segj.cx[iend]* qm;
	}

	segj.jsno= jsnop+1;
	segj.ax[jsnop]=-1.;
	return;

  } /* if( njun2 == 0) */

  qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
  qm=( ap* omc- pp* sd)/ qp;
  qp=-( aj* omc+ pm* sd)/ qp;
  segj.bx[jsnop]=( aj* qm+ ap* qp)* sdh/ sd;
  segj.cx[jsnop]=( aj* qm- ap* qp)* cdh/ sd;

  for( iend = 0; iend < njun1; iend++ )
  {
	segj.ax[iend]= segj.ax[iend]* qm;
	segj.bx[iend]= segj.bx[iend]* qm;
	segj.cx[iend]= segj.cx[iend]* qm;
  }

  jend= njun1;
  for( iend = jend; iend < segj.jsno; iend++ )
  {
	segj.ax[iend]=- segj.ax[iend]* qp;
	segj.bx[iend]= segj.bx[iend]* qp;
	segj.cx[iend]=- segj.cx[iend]* qp;
  }

  segj.jsno= jsnop+1;
  segj.ax[jsnop]=-1.;
}

/*-----------------------------------------------------------------------*/

/* compute the components of all basis functions on segment j */
void trio( int j )
{
  int jcox, jcoxx, jsnox, jx, jend=0, iend=0;

  segj.jsno=0;
  jx = j-1;
  jcox= data.icon1[jx];
  jcoxx = jcox-1;

  if( jcox <= PCHCON)
  {
	jend=-1;
	iend=-1;
  }

  if( (jcox == 0) || (jcox > PCHCON) )
  {
	jcox= data.icon2[jx];
	jcoxx = jcox-1;

	if( jcox <= PCHCON)
	{
	  jend=1;
	  iend=1;
	}

	if( jcox == 0 || (jcox > PCHCON) )
	{
	  jsnox = segj.jsno;
	  segj.jsno++;

	  /* Allocate to connections buffers */
	  if( segj.jsno >= segj.maxcon )
	  {
		segj.maxcon = segj.jsno +1;
		mem_realloc( (void *)&segj.jco, segj.maxcon * sizeof(int) );
		mem_realloc( (void *) &segj.ax, segj.maxcon * sizeof(long double) );
		mem_realloc( (void *) &segj.bx, segj.maxcon * sizeof(long double) );
		mem_realloc( (void *) &segj.cx, segj.maxcon * sizeof(long double) );
	  }

	  sbf( j, j, &segj.ax[jsnox], &segj.bx[jsnox], &segj.cx[jsnox]);
	  segj.jco[jsnox]= j;
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
	  jsnox = segj.jsno;
	  segj.jsno++;

	  /* Allocate to connections buffers */
	  if( segj.jsno >= segj.maxcon )
	  {
		segj.maxcon = segj.jsno +1;
		mem_realloc( (void *)&segj.jco, segj.maxcon * sizeof(int) );
		mem_realloc( (void *) &segj.ax, segj.maxcon * sizeof(long double) );
		mem_realloc( (void *) &segj.bx, segj.maxcon * sizeof(long double) );
		mem_realloc( (void *) &segj.cx, segj.maxcon * sizeof(long double) );
	  }

	  sbf( jcox, j, &segj.ax[jsnox], &segj.bx[jsnox], &segj.cx[jsnox]);
	  segj.jco[jsnox]= jcox;

	  if( jend != 1)
		jcox= data.icon1[jcoxx];
	  else
		jcox= data.icon2[jcoxx];

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

	jcox= data.icon2[jx];

	if( jcox > PCHCON ) break;

	jend=1;
	iend=1;

  } /* do */
  while( jcox != 0 );

  jsnox = segj.jsno;
  segj.jsno++;

  /* Allocate to connections buffers */
  if( segj.jsno >= segj.maxcon )
  {
	segj.maxcon = segj.jsno +1;
	mem_realloc( (void *)&segj.jco, segj.maxcon * sizeof(int) );
	mem_realloc( (void *) &segj.ax, segj.maxcon * sizeof(long double) );
	mem_realloc( (void *) &segj.bx, segj.maxcon * sizeof(long double) );
	mem_realloc( (void *) &segj.cx, segj.maxcon * sizeof(long double) );
  }

  sbf( j, j, &segj.ax[jsnox], &segj.bx[jsnox], &segj.cx[jsnox]);
  segj.jco[jsnox]= j;

  return;

}

/*-----------------------------------------------------------------------*/

/* zint computes the internal impedance of a circular wire */
void zint( long double sigl, long double rolam, complex long double *zint )
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
#define f(d)  ( csqrtl(POT/(d))*cexpl(-cn*(d)+th(-8./x)) )
#define g(d)  ( cexpl(cn*(d)+th(8./x))/csqrtl(TP*(d)) )

  long double x, y, s, ber, bei;
  long double tpcmu = 2.368705e+3;
  long double cmotp = 60.00;
  complex long double br1, br2;

  x= sqrtl( tpcmu* sigl)* rolam;
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
	  *zint= CPLX_01* sqrtl( cmotp/sigl )* br1/ rolam;

  } /* if( x <= 8.) */

	br2= CPLX_01* f(x)/ PI;
	br1= g( x)+ br2;
	br2= g( x)* ph(8./ x)- br2* ph(-8./ x);
	br1= br1/ br2;
	*zint= CPLX_01* sqrtl( cmotp/ sigl)* br1/ rolam;

  } /* if( x <= 110.) */

  br1= cmplx(.70710678,-.70710678);
  *zint= CPLX_01* sqrtl( cmotp/ sigl)* br1/ rolam;
}

/*-----------------------------------------------------------------------*/

/* cang returns the phase angle of a complex number in degrees. */
long double cang( complex long double z )
{
  return( cargl(z)*TD );
}

/*-----------------------------------------------------------------------*/
