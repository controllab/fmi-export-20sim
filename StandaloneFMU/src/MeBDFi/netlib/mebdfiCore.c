/* mebdfi.f -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "../include/f2c.h"

/* Common Block Declarations */

struct {
    doublereal hstpsz[28]	/* was [2][14] */;
} stpsze_;

#define stpsze_1 stpsze_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__5 = 5;
static doublereal c_b318 = 0.;

/************************************************************************/
/* STATICS WHICH MUST BE INITIALIZED BEFORE EACH CALL                   */
/************************************************************************/

/* from mebdfi_ */
static integer g_i1 = 0, g_i2 = 0, g_i3 = 0, g_i4 = 0, g_i5 = 0, g_i6 = 0, g_i7 = 0, g_i8 = 0, g_i9 = 0, g_i10 = 0;

/* from ovdriv_ */
static doublereal g_ovdriv_t = 0.0;
static doublereal g_h__ = 0.0;
static integer g_kflag = 0;
static integer g_nhcut = 0;
static doublereal g_hmin = 0.0, g_hmax = 0.0;
static integer g_jstart = 0;
static doublereal g_crate1 = 0.0, g_crate2 = 0.0;

/* from itrat2_ */
static doublereal g_d1 = 0.0;

/* from stiff_ */
static integer g_imas = 0; /* weird only used, never given a value !!!! */
static integer g_miter = 0;
static integer g_ibnd = 0;
static doublereal g_upbnd = 0.0;
static integer g_nq = 0;
static integer g_l = 0;
static integer g_idoub = 0;
static integer g_kfail = 0;
static doublereal g_stiff_rmax = 0.0;
static doublereal g_rc = 0.0;
static integer g_jsnold = 0;
static logical g_jnewim = 0L;
static doublereal g_tcrat1 = 0.0, g_tcrat2 = 0.0;
static doublereal g_hold = 0.0;
static integer g_mfold = 0;
static logical g_cfail = 0L;
static doublereal g_avnewj = 0.0;
static logical g_sample = 0L;
static integer g_isamp = 0;
static integer g_iemb = 0;
static integer g_meqc1 = 0, g_meqc2 = 0;
static integer g_mq1tmp = 0, g_mq2tmp = 0;
static integer g_lmax = 0;
static integer g_iweval = 0;
static doublereal g_eddn = 0.0;
static doublereal g_eup = 0.0;
static doublereal g_e = 0.0;
static doublereal g_edn = 0.0;
static doublereal g_rh = 0.0;
static integer g_jchang = 0;
static integer g_jsinup = 0;
static integer g_ijus = 0;
static doublereal g_pllfal = 0.0;
static integer g_lmp4 = 0;
static doublereal g_vhold = 0.0;
static doublereal g_dup = 0.0;
static logical g_finish = 0L;
static doublereal g_ffail = 0.0;
static doublereal g_bnd = 0.0;
static doublereal g_qqq = 0.0;

/* from dlamch_
   maybe those really only need to be initialized once !!!
 */
/*          eps   = relative machine precision */
/*          sfmin = safe minimum, such that 1/sfmin does not overflow */
/*          base  = base of the machine */
/*          prec  = eps*base */
/*          t     = number of (base) digits in the mantissa */
/*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/*          emin  = minimum exponent before (gradual) underflow */
/*          rmin  = underflow threshold - base**(emin-1) */
/*          emax  = largest exponent before overflow */
/*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

/* not in use anymore: this goes wrong because there is a one time initialization
   and this will interfere with the get and set functions
   no need to store them globally!!!
static doublereal g_dlamch_rmin, g_dlamch_rmax;
static doublereal g_eps;
static doublereal g_base;
static doublereal g_dlamch_t;
static doublereal g_rnd;
static doublereal g_emin, g_prec, g_emax;
static doublereal g_sfmin;
static doublereal g_rmach;
*/

/************************************************************************/
/*                                                                      */
/************************************************************************/
/*
	34 integers
	25 doublereal
	4 logicals
*/
/************************************************************************/
/*                                                                      */
/************************************************************************/
void RememberMeBDFiStatics(integer *mebdfStaticInts,
							doublereal *mebdfStaticDoubles,
							logical *mebdfStaticLogicals)
{
	/* integers */
	mebdfStaticInts[0] = g_i1;
	mebdfStaticInts[1] = g_i2;
	mebdfStaticInts[2] = g_i3;
	mebdfStaticInts[3] = g_i4;
	mebdfStaticInts[4] = g_i5;
	mebdfStaticInts[5] = g_i6;
	mebdfStaticInts[6] = g_i7;
	mebdfStaticInts[7] = g_i8;
	mebdfStaticInts[8] = g_i9;
	mebdfStaticInts[9] = g_i10;
	mebdfStaticInts[10] = g_kflag;
	mebdfStaticInts[11] = g_nhcut;
	mebdfStaticInts[12] = g_jstart;
	mebdfStaticInts[13] = g_imas;
	mebdfStaticInts[14] = g_miter;
	mebdfStaticInts[15] = g_ibnd;
	mebdfStaticInts[16] = g_nq;
	mebdfStaticInts[17] = g_l;
	mebdfStaticInts[18] = g_idoub;
	mebdfStaticInts[19] = g_kfail;
	mebdfStaticInts[20] = g_jsnold;
	mebdfStaticInts[21] = g_mfold;
	mebdfStaticInts[22] = g_isamp;
	mebdfStaticInts[23] = g_iemb;
	mebdfStaticInts[24] = g_meqc1;
	mebdfStaticInts[25] = g_meqc2;
	mebdfStaticInts[26] = g_mq1tmp;
	mebdfStaticInts[27] = g_mq2tmp;
	mebdfStaticInts[28] = g_lmax;
	mebdfStaticInts[29] = g_iweval;
	mebdfStaticInts[30] = g_jchang;
	mebdfStaticInts[31] = g_jsinup;
	mebdfStaticInts[32] = g_ijus;
	mebdfStaticInts[33] = g_lmp4;
/*
	mebdfStaticInts[34] = g_dlamc1_lbeta;
	mebdfStaticInts[35] = g_dlamch1_lt;
	mebdfStaticInts[36] = g_dlamc2_lt;
	mebdfStaticInts[37] = g_dlamc2_lbeta;
	mebdfStaticInts[38] = g_lemin;
	mebdfStaticInts[39] = g_lemax;
*/

	/* doublereals */
	mebdfStaticDoubles[0] = g_ovdriv_t;
	mebdfStaticDoubles[1] = g_h__;
	mebdfStaticDoubles[2] = g_hmin;
	mebdfStaticDoubles[3] = g_hmax;
	mebdfStaticDoubles[4] = g_crate1;
	mebdfStaticDoubles[5] = g_crate2;
	mebdfStaticDoubles[6] = g_d1;
	mebdfStaticDoubles[7] = g_upbnd;
	mebdfStaticDoubles[8] = g_stiff_rmax;
	mebdfStaticDoubles[9] = g_rc;
	mebdfStaticDoubles[10] = g_tcrat1;
	mebdfStaticDoubles[11] = g_tcrat2;
	mebdfStaticDoubles[12] = g_hold;
	mebdfStaticDoubles[13] = g_avnewj;
	mebdfStaticDoubles[14] = g_eddn;
	mebdfStaticDoubles[15] = g_eup;
	mebdfStaticDoubles[16] = g_e;
	mebdfStaticDoubles[17] = g_edn;
	mebdfStaticDoubles[18] = g_rh;
	mebdfStaticDoubles[19] = g_pllfal;
	mebdfStaticDoubles[20] = g_vhold;
	mebdfStaticDoubles[21] = g_dup;
	mebdfStaticDoubles[22] = g_ffail;
	mebdfStaticDoubles[23] = g_bnd;
	mebdfStaticDoubles[24] = g_qqq;

/*
	mebdfStaticDoubles[25] = g_dlamch_rmin;
	mebdfStaticDoubles[26] = g_dlamch_rmax;
	mebdfStaticDoubles[27] = g_eps;
	mebdfStaticDoubles[28] = g_base;
	mebdfStaticDoubles[29] = g_dlamch_t;
	mebdfStaticDoubles[30] = g_rnd;
	mebdfStaticDoubles[31] = g_emin;
	mebdfStaticDoubles[32] = g_prec;
	mebdfStaticDoubles[33] = g_emax;
	mebdfStaticDoubles[34] = g_sfmin;
	mebdfStaticDoubles[35] = g_rmach;
*/

/*
	mebdfStaticDoubles[25] = g_lrmin;
	mebdfStaticDoubles[26] = g_lrmax;
	mebdfStaticDoubles[27] = g_oldy;
	mebdfStaticDoubles[28] = g_leps;
*/

	/* logicals */
	mebdfStaticLogicals[0] = g_jnewim;
	mebdfStaticLogicals[1] = g_cfail;
	mebdfStaticLogicals[2] = g_sample;
	mebdfStaticLogicals[3] = g_finish;
/*
	mebdfStaticLogicals[4] = g_dlamc1_lrnd;
	mebdfStaticLogicals[5] = g_lieee1;
	mebdfStaticLogicals[6] = g_dlamc2_lrnd;
*/
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
/*
	34 integers
	25 doublereal
	4 logicals
*/
/************************************************************************/
/*                                                                      */
/************************************************************************/
void SetMeBDFiStatics(integer *mebdfStaticInts,
						doublereal *mebdfStaticDoubles,
						logical *mebdfStaticLogicals)
{
	/* integers */
	g_i1 = mebdfStaticInts[0];
	g_i2 = mebdfStaticInts[1];
	g_i3 = mebdfStaticInts[2];
	g_i4 = mebdfStaticInts[3];
	g_i5 = mebdfStaticInts[4];
	g_i6 = mebdfStaticInts[5];
	g_i7 = mebdfStaticInts[6];
	g_i8 = mebdfStaticInts[7];
	g_i9 = mebdfStaticInts[8];
	g_i10 = mebdfStaticInts[9];
	g_kflag = mebdfStaticInts[10];
	g_nhcut = mebdfStaticInts[11];
	g_jstart = mebdfStaticInts[12];
	g_imas = mebdfStaticInts[13];
	g_miter = mebdfStaticInts[14];
	g_ibnd = mebdfStaticInts[15];
	g_nq = mebdfStaticInts[16];
	g_l = mebdfStaticInts[17];
	g_idoub = mebdfStaticInts[18];
	g_kfail = mebdfStaticInts[19];
	g_jsnold = mebdfStaticInts[20];
	g_mfold = mebdfStaticInts[21];
	g_isamp = mebdfStaticInts[22];
	g_iemb = mebdfStaticInts[23];
	g_meqc1 = mebdfStaticInts[24];
	g_meqc2 = mebdfStaticInts[25];
	g_mq1tmp = mebdfStaticInts[26];
	g_mq2tmp = mebdfStaticInts[27];
	g_lmax = mebdfStaticInts[28];
	g_iweval = mebdfStaticInts[29];
	g_jchang = mebdfStaticInts[30];
	g_jsinup = mebdfStaticInts[31];
	g_ijus = mebdfStaticInts[32];
	g_lmp4 = mebdfStaticInts[33];
/*
	g_dlamc1_lbeta = mebdfStaticInts[34];
	g_dlamch1_lt = mebdfStaticInts[35];
	g_dlamc2_lt = mebdfStaticInts[36];
	g_dlamc2_lbeta = mebdfStaticInts[37];
	g_lemin = mebdfStaticInts[38];
	g_lemax = mebdfStaticInts[39];
*/

	/* doublereals */
	g_ovdriv_t = mebdfStaticDoubles[0];
	g_h__ = mebdfStaticDoubles[1];
	g_hmin = mebdfStaticDoubles[2];
	g_hmax = mebdfStaticDoubles[3];
	g_crate1 = mebdfStaticDoubles[4];
	g_crate2 = mebdfStaticDoubles[5];
	g_d1 = mebdfStaticDoubles[6];
	g_upbnd = mebdfStaticDoubles[7];
	g_stiff_rmax = mebdfStaticDoubles[8];
	g_rc = mebdfStaticDoubles[9];
	g_tcrat1 = mebdfStaticDoubles[10];
	g_tcrat2 = mebdfStaticDoubles[11];
	g_hold = mebdfStaticDoubles[12];
	g_avnewj = mebdfStaticDoubles[13];
	g_eddn = mebdfStaticDoubles[14];
	g_eup = mebdfStaticDoubles[15];
	g_e = mebdfStaticDoubles[16];
	g_edn = mebdfStaticDoubles[17];
	g_rh = mebdfStaticDoubles[18];
	g_pllfal = mebdfStaticDoubles[19];
	g_vhold = mebdfStaticDoubles[20];
	g_dup = mebdfStaticDoubles[21];
	g_ffail = mebdfStaticDoubles[22];
	g_bnd = mebdfStaticDoubles[23];
	g_qqq = mebdfStaticDoubles[24];
/*
	g_dlamch_rmin = mebdfStaticDoubles[25];
	g_dlamch_rmax = mebdfStaticDoubles[26];
	g_eps = mebdfStaticDoubles[27];
	g_base = mebdfStaticDoubles[28];
	g_dlamch_t = mebdfStaticDoubles[29];
	g_rnd = mebdfStaticDoubles[30];
	g_emin = mebdfStaticDoubles[31];
	g_prec = mebdfStaticDoubles[32];
	g_emax = mebdfStaticDoubles[33];
	g_sfmin = mebdfStaticDoubles[34];
	g_rmach = mebdfStaticDoubles[35];
*/
/*
	g_lrmin = mebdfStaticDoubles[25];
	g_lrmax = mebdfStaticDoubles[26];
	g_oldy = mebdfStaticDoubles[27];
	g_leps = mebdfStaticDoubles[28];
*/

	/* logicals */
	g_jnewim = mebdfStaticLogicals[0];
	g_cfail = mebdfStaticLogicals[1];
	g_sample = mebdfStaticLogicals[2];
	g_finish = mebdfStaticLogicals[3];
/*
	g_dlamc1_lrnd = mebdfStaticLogicals[4];
	g_lieee1 = mebdfStaticLogicals[5];
	g_dlamc2_lrnd = mebdfStaticLogicals[6];
*/
}




/************************************************************************/
/* END STATICS                                                          */
/************************************************************************/

/*<    >*/
/* Subroutine */ int mebdfi_(integer *n, doublereal *t0, doublereal *ho, 
	doublereal *y0, doublereal *yprime, doublereal *tout, doublereal *
	tend, integer *mf, integer *idid, integer *lout, integer *lwork, 
	doublereal *work, integer *liwork, integer *iwork, integer *mbnd, 
	integer *maxder, integer *itol, doublereal *rtol, doublereal *atol, 
	doublereal *rpar, integer *ipar, U_fp pderv, U_fp resid, integer *
	ierr)
{
    /* Format strings */
    static char fmt_9020[] = "(/,/,\002 ***** ERROR ***** INTEGRATION HALTED\
 IN DRIVER\002,/,/,\002   >>> ILLEGAL VALUE FOR NUMBER OF EQUATIONS <<< \002\
,/,\002                     WITH N = \002,i6)";
    static char fmt_9000[] = "(/,/,\002 ***** ERROR ***** INTEGRATION HALTED\
 IN DRIVER\002,/,/,\002   >>> REAL WORKSPACE IS INSUFFICIENT <<< \002,/,\002\
       WORKSPACE MUST BE AT LEAST \002,i8,\002 ELEMENTS LONG\002)";
    static char fmt_9010[] = "(/,/,\002 ***** ERROR ***** INTEGRATION HALTED\
 IN DRIVER\002,/,/,\002   >>> INTEGER WORKSPACE IS INSUFFICIENT <<< \002,/\
,\002       WORKSPACE MUST BE AT LEAST \002,i6,\002 ELEMENTS LONG\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double sqrt(doublereal);

    /* Local variables */
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int ovdriv_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, U_fp, U_fp, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *);
    
	
	/* epsjac is the machine precision, won't change for other models */
	/* only dependent on the machine */
	/* uround will be the sqrt(epsjac) */
	/* uround will be placed in work[1] */
	/* epsjac won't be placed somewhere */
	static doublereal uround;
	static doublereal l_epsjac;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_9010, 0 };



/* *********************************************************************** */
/* *********************************************************************** */
/*    Written by T.J. Abdulla and J.R. Cash, */
/*    Department of Mathematics, */
/*    Imperial College, */
/*    London SW7 2AZ */
/*    England */

/*    t.abdulla@ic.ac.uk   or  j.cash@ic.ac.uk */

/*    The author would be pleased to receive any comments */
/*          good or bad!! */
/* *********************************************************************** */
/* *********************************************************************** */


/*     THIS IS THE SEPTEMBER 20th 1999 VERSION OF OVDRIV, A PACKAGE FOR */
/*     THE SOLUTION OF THE INITIAL VALUE PROBLEM FOR SYSTEMS OF */
/*     IMPLICIT DIFFERENTIAL ALGEBRAIC EQUATIONS */
/*     G(t,Y,Y')=0, Y=(Y(1),Y(2),Y(3),.....,Y(N)). */
/*     SUBROUTINE OVDRIV IS A DRIVER ROUTINE FOR THIS PACKAGE. */

/*                    REFERENCES */

/*     1.  J. R. CASH, THE INTEGRATION OF STIFF INITIAL VALUE PROBLEMS */
/*         IN O.D.E.S USING MODIFIED EXTENDED BACKWARD DIFFERENTIATION */
/*         FORMULAE, COMP. AND MATHS. WITH APPLICS., 9, 645-657, (1983). */
/*     2.  J.R. CASH AND S. CONSIDINE, AN MEBDF CODE FOR STIFF */
/*         INITIAL VALUE PROBLEMS, ACM TRANS MATH SOFTWARE, 142-158, */
/*         (1992). */
/*     3.  J.R. CASH, STABLE RECURSIONS WITH APPLICATIONS TO THE */
/*         NUMERICAL SOLUTION OF STIFF SYSTEMS, ACADEMIC PRESS,(1979). */
/*     4.  A.C. HINDMARSH, ODEPACK, A SYSTEMISED COLLECTION OF ODE */
/*         SOLVERS, in SCIENTIFIC COMPUTING, R.S. STEPLEMAN et. al. */
/*         (eds) North-Holland, AMSTERDAM, pp55-64 , (1983). */
/*     5.  E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL */
/*         EQUATIONS II, STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS, */
/*         SPRINGER 1996, page 267. */

/*     ---------------------------------------------------------------- */
/*     OVDRIV IS TO BE CALLED ONCE FOR EACH OUTPUT VALUE OF T, AND */
/*     IN TURN MAKES REPEATED CALLS TO THE CORE INTEGRATOR STIFF. */

/*     THE INPUT PARAMETERS ARE .. */
/*     N     =  THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS. */
/*     T0    =  THE INITIAL VALUE OF T, THE INDEPENDENT VARIABLE */
/*              (USED ONLY ON THE FIRST CALL) */
/*     HO    =  THE NEXT STEP SIZE IN T (USED FOR INPUT ONLY ON THE */
/*              FIRST CALL) */
/*     Y0    =  A VECTOR OF LENGTH N CONTAINING THE INITIAL VALUES OF Y */
/*              (USED FOR INPUT ONLY ON FIRST CALL) */
/*     YPRIME   A VECTOR OF LENGTH N CONTAINING THE INITIAL VALUES OF */
/*              DY/DT */
/*     TOUT  =  THE VALUE OF T AT WHICH OUTPUT IS DESIRED NEXT. */
/*              INTEGRATION WILL NORMALLY GO SLIGHTLY BEYOND TOUT */
/*              AND THE PACKAGE WILL INTERPOLATE TO T = TOUT */
/*     TEND  =  END OF THE RANGE OF INTEGRATION. */
/*     MF    =  THE METHOD FLAG.  AT PRESENT MF=21,22,23 OR 24 IS */
/*              ALLOWED. THESE ARE EXTENDED BACKWARD DIFFERENTIATION */
/*              FORMULAE USING THE CHORD METHOD WITH ANALYTIC OR NUMERICAL */
/*              JACOBIAN FOR MF=21,22 RESPECTIVELY. MF=23/24 ARE  THE SAME */
/*              AS FOR 21/22 BUT THE JACOBIAN IS NOW BANDED.   THE USER */
/*              NEEDS TO SPECIFY SUBROUTINE PDERV IF MF=21 OR 23. */
/*     IDID   = THE INTEGER USED ON INPUT TO INDICATE THE TYPE OF CALL. */
/*              1   THIS IS THE FIRST CALL FOR THE PROBLEM. */
/*              0   THIS IS NOT THE FIRST CALL FOR THIS PROBLEM */
/*                  AND INTEGRATION IS TO CONTINUE. */
/*             -1   THIS IS NOT THE FIRST CALL FOR THE PROBLEM, */
/*                  AND THE USER HAS RESET N, RTOL, ATOL,H  AND/OR MF. */
/*              2   SAME AS 0 EXCEPT THAT TOUT HAS TO BE HIT */
/*                  EXACTLY (NO INTERPOLATION IS DONE). */
/*                  ASSUMES TOUT .GE. THE CURRENT T. */
/*              3   SAME AS 0 EXCEPT CONTROL RETURNS TO CALLING */
/*                  PROGRAM AFTER ONE STEP. TOUT IS IGNORED, UNTIL THE */
/*                  INTEGRATION REACHES TOUT OR BEYOND. IF IT PASSES TOUT */
/*                  THE PROGRAM INTERPOLATES THE SOLUTION VALUES AND */
/*                  RETURNS THE SOLUTION VALUE AT TOUT. */
/*              SINCE THE NORMAL OUTPUT VALUE OF IDID IS 0, */
/*              IT NEED NOT BE RESET FOR NORMAL CONTINUATION. */
/*              THE FIRST CALL TO THE DRIVER IS WITH IDID=1 AND FOR */
/*              A SUCCESSFUL STEP THE DRIVER RETURNS WITH IDID=1.THUS */
/*              THE CALL WITH IDID = 1 IS SIMPLY THE FIRST */
/*              INITIALISING STEP FOR THE CODE.  THE USER */
/*              THEN NEEDS TO CONTINUE WITH IDID=0,-1,2 OR 3 AS ABOVE. */
/*     LOUT   = THE LOGICAL OUTPUT CHANNEL FOR MESSAGE PASSING. */
/*     MBND   = AN ARRAY OF DIMENSION 4 FOR USE WHEN THE NEWTON ITERATION */
/*              MATRIX IS BANDED.  IF THIS MATRIX HAS ML DIAGONALS */
/*              BELOW THE MAIN DIAGONAL AND MU DIAGONALS ABOVE THE */
/*              MAIN DIAGONAL THEN: */
/*              MBND(1) = ML */
/*              MBND(2) = MU */
/*              MBND(3) = MU + ML + 1 */
/*              MBND(4) = 2*ML + MU + 1 */
/*     MAXDER=  THE MAXIMUM ORDER IS MAXDER + 1. */
/*              THE VALUE OF MAXDER CANNOT EXCEED 7.  THIS IS THE */
/*              VALUE RECOMMENDED UNLESS IT IS BELIEVED THAT THERE */
/*              ARE SEVERE STABILITY PROBLEMS IN WHICH CASE MAXDER=3 */
/*              OR 4 SHOULD BE TRIED INSTEAD. */
/*     ITOL  =  AN INDICATOR OF THE TYPE OF ERROR CONTROL. SEE */
/*              DESCRIPTION BELOW UNDER ATOL. */
/*     RTOL  =  A RELATIVE ERROR TOLERANCE PARAMETER. CAN BE EITHER A */
/*              SCALAR OR AN ARRAY OF LENGTH N.  SEE DESCRIPTION */
/*              BELOW UNDER ATOL. */
/*     ATOL  =  THE ABSOLUTE ERROR BOUND. */
/*              THE INPUT PARAMETERS ITOL, RTOL AND ATOL DETERMINE */
/*              THE ERROR CONTROL PERFORMED BY THE SOLVER.  THE */
/*              SOLVER WILL CONTROL THE VECTOR e = (e(i)) OF ESTIMATED */
/*              LOCAL ERRORS IN y ACCORDING TO AN INEQUALITY OF THE FORM */
/*                  RMS-NORM OF (e(i)/ewt(i)) .LE. 1 */
/*              THE ROOT MEAN SQUARE NORM IS */
/*                   RMS-NORM(V) = SQRT((SUM v(i)**2)/N).  HERE */
/*                ewt = (ewt(i)) IS A VECTOR OF WEIGHTS WHICH MUST */
/*              ALWAYS BE POSITIVE, AND THE VALUES OF RTOL AND ATOL */
/*              SHOULD BE NON-NEGATIVE. IF ITOL = 1 THEN SINGLE STEP ERROR */
/*              ESTIMATES DIVIDED BY YMAX(I) WILL BE KEPT LESS THAN 1 */
/*              IN ROOT-MEAN-SQUARE NORM.  THE VECTOR YMAX OF WEIGHTS IS */
/*              COMPUTED IN OVDRIV. INITIALLY YMAX(I) IS SET AS */
/*              THE MAXIMUM OF 1 AND ABS(Y(I)).  THEREAFTER YMAX(I) IS */
/*              THE LARGEST VALUE OF ABS(Y(I)) SEEN SO FAR, OR THE */
/*              INITIAL VALUE YMAX(I) IF THAT IS LARGER. */
/*              IF ITOL = 1 THE USER NEEDS TO SET ATOL = RTOL = */
/*              THE PRECISION REQUIRED.  THEN */
/*                         ewt(i) = RTOL(1)*YMAX(i) */
/*              IF ITOL IS GREATER THAN 1 THEN */
/*                 ewt(i) = rtol(i)*abs(y(i)) + atol(i) */
/*              THE FOLLOWING TABLE GIVES THE TYPES (SCALAR/ARRAY) */
/*              OF RTOL AND ATOL, AND THE CORRESPONDING FORM OF ewt(i) */
/*                  ITOL   RTOL      ATOL       ewt(i) */
/*                   2    SCALAR    SCALAR   rtol*abs(y(i))   + atol */
/*                   3    SCALAR    ARRAY    rtol*abs(y(i))   + atol(i) */
/*                   4    ARRAY     SCALAR   rtol(i)*abs(y(i))+ atol */
/*                   5    ARRAY     ARRAY    rtol(i)*abs(y(i))+ atol(i) */
/*              IF EITHER OF THESE PARAMETERS IS A SCALAR, IT NEED */
/*              NOT BE DIMENSIONED IN THE USER'S CALLING PROGRAM. */
/*     NIND1    = THE NUMBER OF VARIABLES OF INDEX 1,2,3 RESPECTIVELY. */
/*     NIND2,NIND3  THESE ARE SET IN IWORK(1),(2),(3). */
/*              THE EQUATIONS MUST BE DEFINED SO THAT THE INDEX 1 */
/*              VARIABLES PRECEDE THE INDEX 2 VARIABLES WHICH IN */
/*              TURN PRECEDE THE INDEX 3 VARIABLES. */
/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) */
/*              WHICH CAN BE USED FOR COMMUNICATION BETWEEN  THE USER'S */
/*              CALLING PROGRAM AND THE F AND PDERV SUBROUTINES. */
/*     IERR     IERR IS AN INTEGER FLAG WHICH IS ALWAYS EQUAL TO ZERO */
/*              ON INPUT.  SUBROUTINES F AND PDERV SHOULD ALTER */
/*              IERR ONLY IF ONE OF THEM ENCOUNTERS AN ILLEGAL OPERATION SUCH */
/*              AS THE SQUARE ROOT OF A NEGATIVE NUMBER OR EXPONENT */
/*              OVERFLOW. THE USER CAN THEN ALTER H AND CALL THE */
/*              SUBROUTINE AGAIN WITH IDID=-1 IF HE WISHES. */

/*     AFTER THE INITIAL CALL, IF A NORMAL RETURN OCCURED AND A NORMAL */
/*     CONTINUATION IS DESIRED, SIMPLY RESET TOUT AND CALL AGAIN. */
/*     ALL OTHER PARAMETERS WILL BE READY FOR THE NEXT CALL. */
/*     A CHANGE OF PARAMETERS WITH IDID = -1 CAN BE MADE AFTER */
/*     EITHER A SUCCESSFUL OR AN UNSUCCESSFUL RETURN. */

/*     THE OUTPUT PARAMETERS ARE.. */
/*     T0    =  THE VALUE OF T WHICH RELATES TO THE CURRENT SOLUTION */
/*              POINT Y0() */
/*     HO    =  THE STEPSIZE H USED LAST, WHETHER SUCCESSFULLY OR NOT. */
/*     Y0    =  THE COMPUTED VALUES OF Y AT T = TOUT */
/*     YPRIME=  THE COMPUTED VALUES OF DY/DT AT T=TOUT. */
/*     TOUT  =  UNCHANGED FROM ITS INPUT VALUE. */
/*     IDID  =  INTEGER USED ON OUTPUT TO INDICATE RESULTS, WITH */
/*              THE FOLLOWING VALUES AND MEANINGS.. */

/*      0   INTEGRATION WAS COMPLETED TO TOUT OR BEYOND. */

/*     -1   THE INTEGRATION WAS HALTED AFTER FAILING TO PASS THE */
/*          ERROR TEST EVEN AFTER REDUCING H BY A FACTOR OF */
/*          1.E10 FROM ITS INITIAL VALUE. */

/*     -2   AFTER SOME INITIAL SUCCESS, THE INTEGRATION WAS */
/*          HALTED EITHER BY REPEATED ERROR TEST FAILURES OR BY */
/*          A TEST ON RTOL/ATOL.  TOO MUCH ACCURACY HAS BEEN REQUESTED. */

/*     -3   THE INTEGRATION WAS HALTED AFTER FAILING TO ACHIEVE */
/*          CORRECTOR CONVERGENCE EVEN AFTER REDUCING H BY A */
/*          FACTOR OF 1.E10 FROM ITS INITIAL VALUE. */

/*     -4   IMMEDIATE HALT BECAUSE OF ILLEGAL VALUES OF INPUT */
/*          PARAMETERS.  SEE PRINTED MESSAGE. */

/*     -5   IDID WAS -1 ON INPUT, BUT THE DESIRED CHANGES OF */
/*          PARAMETERS WERE NOT IMPLEMENTED BECAUSE TOUT */
/*          WAS NOT BEYOND T.  INTERPOLATION AT T = TOUT WAS */
/*          PERFORMED AS ON A NORMAL RETURN.  TO TRY AGAIN, */
/*          SIMPLY CALL AGAIN WITH IDID = -1 AND A NEW TOUT. */

/*     -6   MAXIMUM ALLOWABLE NUMBER OF INTEGRATION STEPS EXCEEDED. */
/*          TO CONTINUE THE USER SHOULD RESET IWORK(14). */


/*     -7   STEPSIZE IS TOO SMALL (LESS THAN SQRT(UROUND)/100) */


/*     -11   INSUFFICIENT REAL WORKSPACE FOR THE INTEGRATION */

/*     -12   INSUFFICIENT INTEGER WORKSPACE FOR THE INTEGRATION */


/*     IN ADDITION TO OVDRIVE, THE FOLLOWING ROUTINES ARE PROVIDED */
/*     IN THE PACKAGE.. */

/*     INTERP( - )   INTERPOLATES TO GET THE OUTPUT VALUES */
/*                   AT T=TOUT FROM THE DATA IN THE Y ARRAY. */
/*     STIFF( - )    IS THE CORE INTEGRATOR ROUTINE.  IT PERFORMS A */
/*                   SINGLE STEP AND ASSOCIATED ERROR CONTROL. */
/*     COSET( - )    SETS COEFFICIENTS FOR BACKWARD DIFFERENTIATION */
/*                   SCHEMES FOR USE IN THE CORE INTEGRATOR. */
/*     PSET( - )     COMPUTES AND PROCESSES THE NEWTON ITERATION */
/*                   MATRIX DG/DY + (1/(H*BETA))DG/DY' */
/*     DEC( - )      PERFORMS AN LU DECOMPOSITION ON A MATRIX. */
/*     SOL( - )      SOLVES LINEAR SYSTEMS A*X = B AFTER DEC */
/*                   HAS BEEN CALLED FOR THE MATRIX A */
/*     DGBFA ( - )   FACTORS A DOUBLE PRECISION BAND MATRIX BY */
/*                   ELIMINATION. */
/*     DGBSL ( - )   SOLVES A BANDED LINEAR SYSTEM A*x=b */

/*                   ALSO SUPPLIED ARE THE BLAS ROUTINES */

/*                   daxpy, dscal, idamax, ddot. */


/*     THE FOLLOWING ROUTINES ARE TO BE SUPPLIED BY THE USER AND */
/*                   SHOULD BE DECLARED AS EXTERNAL. */

/*     PDERV(T,Y,PD,N,YPRIME,MBND(4),CON,IPAR,RPAR,IERR) */
/*                         COMPUTES THE N*N NEWTON ITERATION MATRIX */
/*                         OF PARTIAL DERIVATIVES PD=DG/DY + c(DG/DY') */
/*                         WHERE c IS A SCALAR.  IF A NUMERICAL */
/*                         JACOBIAN IS REQUIRED (MF=22 OR 24) THEN THIS */
/*                         SUBROUTINE CAN BE A DUMMY ONE. THE ITERATION */
/*                         MATRIX IS STORED AS AN N BY N ARRAY IF THE */
/*                         MATRIX IS FULL.  IF THE ITERATION MATRIX IS */
/*                         BANDED THE ARRAY PD IS OF SIZE MBND(4)*N. */
/*                         IF THE ITERATION MATRIX IS FULL, PD(I,J) IS */
/*                         TO BE SET TO THE PARTIAL DERIVATIVE OF THE */
/*                         ith COMPONENT OF G(t,Y,Y') WITH */
/*                         RESPECT TO Y(J).  IF THE JACOBIAN IS BANDED */
/*                         WITH mu DIAGONALS ABOVE THE MAIN DIAGONAL */
/*                         THE PARTIAL DERIVATIVE OF THE ith COMPONENT */
/*                         OF G(t,Y,Y') WITH RESPECT TO Y(J) SHOULD BE */
/*                         PUT IN PD(i-j+mu+1,j). PDERV IS CALLED ONLY IF */
/*                         MITER = 1 OR 3.  OTHERWISE A DUMMY ROUTINE CAN */
/*                         BE SUBSTITUTED. */

/*     THE DIMENSION OF PD  MUST BE AT LEAST N**2 FOR A FULL ITERATION */
/*     MATRIX AND MBND(4)*N FOR A BANDED ITERATION MATRIX */
/*     (MF=23 OR 24) .  THE DIMENSIONS */
/*     OF YMAX,ERROR,SAVE1,SAVE2,IPIV AND THE FIRST DIMENSION */
/*     OF Y SHOULD ALL BE AT LEAST N. */

/*     RESID(N,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR) */
/*                          COMPUTES THE RESIDUAL VECTOR */
/*                          DELTA = G(T,Y,Y'). */


/*     UROUND   THIS IS THE UNIT ROUNDOFF AND HAS TO BE SET AS */
/*                        UROUND = DLAMCH('Epsilon') */
/*     EPSJAC   = sqrt(UROUND). */

/*     HUSED  (=WORK(2))    LAST STEPSIZE SUCCESSFULLY USED BY THE INTEGRATOR */
/*     NQUSED (=IWORK(4))   LAST ORDER SUCCESSFULLY USED */
/*     NSTEP  (=IWORK(5))   NUMBER OF SUCCESSFUL STEPS TAKEN SO FAR */
/*     NFAIL  (=IWORK(6))   NUMBER OF FAILED STEPS */
/*     NRE    (=IWORK(7))   NUMBER OF RESIDUAL EVALUATIONS  SO FAR */
/*     NJE    (=IWORK(8))   NUMBER OF JACOBIAN EVALUATIONS  SO FAR */
/*     NDEC   (=IWORK(9))   NUMBER OF LU DECOMPOSITIONS  SO FAR */
/*     NBSOL  (=IWORK(10))   NUMBER OF 'BACKSOLVES'  SO FAR */
/*     NPSET  (=IWORK(11))   NUMBER OF TIMES A NEW COEFFICIENT MATRIX HAS BEEN */
/*               FORMED SO FAR */
/*     NCOSET (=IWORK(12))   NUMBER OF TIMES THE ORDER OF THE METHOD USED HAS */
/*               BEEN CHANGED SO FAR */
/*     MAXORD (=IWORK(13))   THE MAXIMUM ORDER USED SO FAR IN THE INTEGRATION */
/*     MAXSTP (=IWORK(14))   THE MAXIMUM ALLOWED NUMBER OF STEPS SET BY THE */
/*                           USER. */

/*    IF IT IS ONLY REQUIRED TO CONTROL THE ACCURACY IN THE */
/*    DIFFERENTIAL VARIABLES THEN THE USER SHOULD FIND THE */
/*    STRING 'AMMEND' AND MAKE CHANGES THERE */


/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     THIS SUBROUTINE IS FOR THE PURPOSE               * */
/*     OF SPLITTING UP THE WORK ARRAYS WORK AND IWORK   * */
/*     FOR USE INSIDE THE INTEGRATOR STIFF              * */
/*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER IDID,LIWORK,LOUT,LWORK,MF,N,MAXDER,ITOL >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<    >*/
/*<       INTEGER IWORK(LIWORK), MBND(4) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10 >*/
/*     COMMON BLOCKS */
/*     .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*<       EXTERNAL OVDRIV,PDERV,RESID >*/
/*     .. */
/*     .. SAVE STATEMENT .. */
/*<       SAVE  I1,I2,I3,I4,I5,I6,I7,I8,I9,I10 >*/
/*     .. */
/*<       IF (IDID.EQ.1) THEN >*/
    /* Parameter adjustments */
    --yprime;
    --y0;
    --work;
    --iwork;
    --mbnd;
    --rtol;
    --atol;
    --rpar;
    --ipar;

    /* Function Body */
    if (*idid == 1) {
/*<          IF (N.LE.0) THEN >*/
	if (*n <= 0) {
/*<             WRITE (LOUT,9020) N >*/
	    io___1.ciunit = *lout;
	    s_wsfe(&io___1);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          ELSE >*/
	} else {
/*<             IF(MF.LT.23) MBND(4)=N >*/
	    if (*mf < 23) {
		mbnd[4] = *n;
	    }
/*<             I1 = N*12+ 3 >*/
	    g_i1 = *n * 12 + 3;
/*<             I2 = I1 + N*12 >*/
	    g_i2 = g_i1 + *n * 12;
/*<             I3 = I2 + N*2 >*/
	    g_i3 = g_i2 + (*n << 1);
/*<             I4 = I3 + N >*/
	    g_i4 = g_i3 + *n;
/*<             I5 = I4 + N >*/
	    g_i5 = g_i4 + *n;
/*<             I6 = I5 + N >*/
	    g_i6 = g_i5 + *n;
/*<             I7 = I6 + N >*/
	    g_i7 = g_i6 + *n;
/*<             I8 = I7 + N >*/
	    g_i8 = g_i7 + *n;
/*<             I9 = I8 + N >*/
	    g_i9 = g_i8 + *n;
/*<             I10 = I9 + MBND(4)*N >*/
	    g_i10 = g_i9 + mbnd[4] * *n;
/*<             UROUND = DLAMCH('Epsilon') >*/
	    uround = dlamch_("Epsilon", (ftnlen)7);
/*<             WORK(1) = UROUND >*/
	    work[1] = uround;
/*<             EPSJAC = SQRT(WORK(1)) >*/
	    l_epsjac = sqrt(work[1]);
/*<             IF (LWORK.LT.(I10+1)) THEN >*/
	    if (*lwork < g_i10 + 1) {
/*<                IDID = -11 >*/
		*idid = -11;
/*<                WRITE (LOUT,9000) I10 + 1 >*/
		io___14.ciunit = *lout;
		s_wsfe(&io___14);
		i__1 = g_i10 + 1;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		e_wsfe();
/*<             ENDIF >*/
	    }
/*<             IF (LIWORK.LT.N+14) THEN >*/
	    if (*liwork < *n + 14) {
/*<                IDID = -12 >*/
		*idid = -12;
/*<                WRITE (LOUT,9010) N+14 >*/
		io___15.ciunit = *lout;
		s_wsfe(&io___15);
		i__1 = *n + 14;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		e_wsfe();
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<          IF (IDID.LT.0) RETURN >*/
	if (*idid < 0) {
	    return 0;
	}
/*<       END IF >*/
    }

/*    THE DIMENSION OF THE REAL WORKSPACE, WORK, HAS TO BE AT LEAST */
/*     (32 + MBND(4))*N+2 WHILE THE DIMENSION OF THE INTEGER */
/*    WORKSPACE HAS TO BE AT LEAST N+14. */

/*<    >*/
    ovdriv_(n, t0, ho, &y0[1], &yprime[1], tout, tend, mf, idid, lout, &work[
	    3], &work[g_i1], &work[g_i2], &work[g_i3], &work[g_i4], &work[g_i5], &work[
	    g_i6], &work[g_i7], &work[g_i8], &work[g_i9], &work[g_i10], &iwork[15], &
	    mbnd[1], &iwork[1], &iwork[2], &iwork[3], maxder, itol, &rtol[1], 
	    &atol[1], &rpar[1], &ipar[1], (U_fp)pderv, (U_fp)resid, &iwork[4],
	     &iwork[5], &iwork[6], &iwork[7], &iwork[8], &iwork[9], &iwork[10]
	    , &iwork[11], &iwork[12], &iwork[13], &iwork[14], &work[1], &work[
	    2], &l_epsjac, ierr);
/*     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     WORK() HOUSES THE FOLLOWING ARRAYS */

/*     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N) */
/*     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N) */
/*     PWCOPY(MBND(4)*N) */
/*     IF THE BANDED OPTION IS NOT BEING USED THEN MBND(4)=N. */
/*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*<       RETURN >*/
    return 0;
/*<  9 >*/
/*<  9 >*/
/*<  9 >*/
/*<       END >*/
} /* mebdfi_ */

/* -------------------------------------------------------------------------- */


/*<    >*/
/* Subroutine */ int ovdriv_(integer *n, doublereal *t0, doublereal *ho, 
	doublereal *y0, doublereal *yprime, doublereal *tout, doublereal *
	tend, integer *mf, integer *idid, integer *lout, doublereal *y, 
	doublereal *yhold, doublereal *ynhold, doublereal *ymax, doublereal *
	errors, doublereal *save1, doublereal *save2, doublereal *scale, 
	doublereal *arh, doublereal *pw, doublereal *pwcopy, integer *ipiv, 
	integer *mbnd, integer *nind1, integer *nind2, integer *nind3, 
	integer *maxder, integer *itol, doublereal *rtol, doublereal *atol, 
	doublereal *rpar, integer *ipar, U_fp pderv, U_fp resid, integer *
	nqused, integer *nstep, integer *nfail, integer *nre, integer *nje, 
	integer *ndec, integer *nbsol, integer *npset, integer *ncoset, 
	integer *maxord, integer *maxstp, doublereal *uround, doublereal *
	hused, doublereal *epsjac, integer *ierr)
{
    /* Format strings */
    static char fmt_9160[] = "(/,/,\002STEPSIZE IS TOO SMALL\002)";
    static char fmt_9080[] = "(/,/,\002 IDID = -1 ON INPUT WITH (T-TOUT)*H .\
GE. 0.\002,/,\002 T =\002,e16.8,\002   TOUT =\002,e16.8,\002   H =\002,e16.8\
,/,\002 INTERPOLATION WAS DONE AS ON NORMAL RETURN.\002,/,\002 DESIRED PARAM\
ETER CHANGES WERE NOT MADE.\002)";
    static char fmt_9070[] = "(/,/,\002 ILLEGAL INPUT.. IDID =\002,i5,/,/)";
    static char fmt_9040[] = "(/,/,\002 ILLEGAL INPUT.. RTOL .LE. 0.\002,/,/)"
	    ;
    static char fmt_9045[] = "(/,/,\002 ILLEGAL INPUT.. ATOL .LE. 0.\002,/,/)"
	    ;
    static char fmt_9050[] = "(/,/,\002 ILLEGAL INPUT.. N .LE. 0\002,/,/)";
    static char fmt_9060[] = "(/,/,\002 ILLEGAL INPUT.. (T0-TOUT)*H .GE. 0\
.\002,/,/)";
    static char fmt_9090[] = "(/,/,\002 ILLEGAL INPUT.. METHOD FLAG, MF, =\
 \002,i6,/,\002         ALLOWED VALUES ARE 21 OR 22\002,/)";
    static char fmt_9110[] = "(/,/,\002 ILLEGAL VALUE FOR ITOL\002,/,/)";
    static char fmt_9120[] = "(/,/,\002 ILLEGAL VALUE FOR MAXDER\002,/,/)";
    static char fmt_9140[] = "(/,/,\002BAD INPUT FOR NUMBER OF VARIABLES OF \
INDEX 1,2,3\002,/,/)";
    static char fmt_9000[] = "(\002 WARNING..  T + H = T ON NEXT STEP.\002)";
    static char fmt_9010[] = "(/,/,\002 KFLAG = -2 FROM INTEGRATOR AT T =\
 \002,e16.8,\002  H =\002,e16.8,/,\002  THE REQUESTED ERROR IS SMALLER THAN \
CAN BE HANDLED\002,/,/)";
    static char fmt_9030[] = "(/,/,\002 KFLAG = -3 FROM INTEGRATOR AT T =\
 \002,e16.8,/,\002  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED\002,/)";
    static char fmt_9130[] = "(/,/,\002 NUMBER OF STEPS EXCEEDS MAXIMUM\002,\
/,/)";
    static char fmt_9100[] = "(/,/,\002 PROBLEM APPEARS UNSOLVABLE WITH GIVE\
N INPUT\002,/,\002         HMIN REDUCED BY A FACTOR OF 1.0E10\002,/,/)";

    /* System generated locals */
    integer y_dim1, y_offset, yhold_dim1, yhold_offset, ynhold_dim1, 
	    ynhold_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    
	/* klag is only given a value, it is never used */
	/* so leave it static */
	/**/static/**/ integer klag;

	/* only given a value, never used */
	static doublereal d__;


	/* used as counter in local scope only */
    /*static*/ integer i__;

    extern /* Subroutine */ int stiff_(doublereal *, doublereal *, doublereal 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, U_fp, U_fp, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    
	/* only used locally */
	/*static*/ doublereal vhold;
    /*static*/ integer nn;

    extern /* Subroutine */ int interp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

	/* only used locally */
    /*static*/ integer kgo;
    /*static*/ doublereal ayi;

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 6, 0, fmt_9160, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_9080, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_9070, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_9040, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_9045, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_9040, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_9040, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_9050, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_9060, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_9090, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_9110, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_9120, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_9140, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_9030, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_9130, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_9100, 0 };


/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/

/*     START OF PROGRAM PROPER */

/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER IDID,LOUT,MF,N,MAXDER,ITOL >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<    >*/
/*<       INTEGER IPIV(N), MBND(4) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,KGO,NHCUT >*/
/*     .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*<       EXTERNAL INTERP,STIFF,PDERV,RESID >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*<       INTRINSIC DABS,DMAX1 >*/
/*     .. */
/*     .. COMMON BLOCKS .. */
/*<    >*/
/*<       SAVE T,H,HMIN,HMAX,KFLAG,JSTART >*/
/*     .. */
/*<       IF (IDID.EQ.0) THEN >*/
    /* Parameter adjustments */
    --ipiv;
    --arh;
    --scale;
    --save2;
    --save1;
    --errors;
    --ymax;
    ynhold_dim1 = *n;
    ynhold_offset = 1 + ynhold_dim1 * 1;
    ynhold -= ynhold_offset;
    yhold_dim1 = *n;
    yhold_offset = 1 + yhold_dim1 * 1;
    yhold -= yhold_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    --yprime;
    --y0;
    --pw;
    --pwcopy;
    --mbnd;
    --rtol;
    --atol;
    --rpar;
    --ipar;

    /* Function Body */
    if (*idid == 0) {
/*        I.E. NORMAL CONTINUATION OF INTEGRATION */
/*<          T0=T >*/
	*t0 = g_ovdriv_t;
/*<          HMAX = DABS(TEND-T0)*10.0D+0 >*/
	g_hmax = (d__1 = *tend - *t0, abs(d__1)) * 10.;
/*<          IF ((T-TOUT)*H.GE.0.0D+0) THEN >*/
	if ((g_ovdriv_t - *tout) * g_h__ >= 0.) {
/*           HAVE OVERSHOT THE OUTPUT POINT, SO INTERPOLATE */
/*<             CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
	    interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<             IDID = KFLAG >*/
	    *idid = g_kflag;
/*<             T0 = TOUT >*/
	    *t0 = *tout;
/*<             HO = H >*/
	    *ho = g_h__;
/*<             RETURN >*/
	    return 0;
/*<          END IF >*/
	}
/*<       ELSE IF (IDID.EQ.2) THEN >*/
    } else if (*idid == 2) {
/*        I.E. CONTINUING INTEGRATION BUT WISH TO HIT TOUT */
/*<          T0 = T >*/
	*t0 = g_ovdriv_t;
/*<          HMAX = DABS(TEND-T0)*10.0D+0 >*/
	g_hmax = (d__1 = *tend - *t0, abs(d__1)) * 10.;
/*<          IF (((T+H)-TOUT)*H.GT.0.0D+0) THEN >*/
	if ((g_ovdriv_t + g_h__ - *tout) * g_h__ > 0.) {
/*           WE HAVE ALREADY OVERSHOT THE OUTPUT POINT OR WE WILL */
/*           DO SO ON THE NEXT STEP */
/*<    >*/
	    if ((g_ovdriv_t - *tout) * g_h__ >= 0. || (d__1 = g_ovdriv_t - *tout, abs(d__1)) <= *
		    uround * 100. * g_hmax) {
/*              HAVE OVERSHOT THE OUTPUT POINT, SO INTERPOLATE */
/*<                CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
		interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<                T0 = TOUT >*/
		*t0 = *tout;
/*<                HO = H >*/
		*ho = g_h__;
/*<                IDID = KFLAG >*/
		*idid = g_kflag;
/*<                RETURN >*/
		return 0;
/*<             ELSE >*/
	    } else {
/*              WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE */
/*              SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY' */
/*<                H = (TOUT-T)* (1.0D+0-4.0D+0*UROUND) >*/
		g_h__ = (*tout - g_ovdriv_t) * (1. - *uround * 4.);
/*<                JSTART = -1 >*/
		g_jstart = -1;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF (IDID.EQ.-1) THEN >*/
    } else if (*idid == -1) {
/*        NOT FIRST CALL BUT PARAMETERS RESET */
/*<          H = HO >*/
	g_h__ = *ho;
/*<          IF(H.LT.EPSJAC/100.0D+0) THEN >*/
	if (g_h__ < *epsjac / 100.) {
/*<             WRITE(6,9160) >*/
	    s_wsfe(&io___21);
	    e_wsfe();
/*<             IDID = -7 >*/
	    *idid = -7;
/*<             RETURN >*/
	    return 0;
/*<          ENDIF >*/
	}
/*<          T0 = T >*/
	*t0 = g_ovdriv_t;
/*<          IF ((T-TOUT)*H.GE.0.0D+0) THEN >*/
	if ((g_ovdriv_t - *tout) * g_h__ >= 0.) {
/*           HAVE OVERSHOT TOUT */
/*<             WRITE (LOUT,9080) T,TOUT,H >*/
	    io___22.ciunit = *lout;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&g_ovdriv_t, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g_h__, (ftnlen)sizeof(doublereal));
	    e_wsfe();
/*<             CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
	    interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<             HO = H >*/
	    *ho = g_h__;
/*<             T0 = TOUT >*/
	    *t0 = *tout;
/*<             IDID = -5 >*/
	    *idid = -5;
/*<             RETURN >*/
	    return 0;
/*<          ELSE >*/
	} else {
/*<             JSTART = -1 >*/
	    g_jstart = -1;
/*<          END IF >*/
	}
/*<       ELSE IF (IDID.EQ.3) THEN >*/
    } else if (*idid == 3) {
/*<          T0 = T >*/
	*t0 = g_ovdriv_t;
/*<          IF ((T-TOUT)*H.GE.0.0D+0) THEN >*/
 	if ((g_ovdriv_t - *tout) * g_h__ >= 0.) {
/*           HAVE OVERSHOT,SO INTERPOLATE */
/*<             CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
	    interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<             IDID = KFLAG >*/
	    *idid = g_kflag;
/*<             T0 = TOUT >*/
	    *t0 = *tout;
/*<             HO = H >*/
	    *ho = g_h__;
/*<             RETURN >*/
	    return 0;
/*<          END IF >*/
	}
/*<       ELSE >*/
    } else {
/*        IDID SHOULD BE 1 AND THIS IS THE FIRST CALL FOR THIS PROBLEM */
/*        CHECK THE ARGUMENTS THAT WERE PASSED FOR CORRECTNESS */
/*<          IF (IDID.NE.1) THEN >*/
	if (*idid != 1) {
/*           VALUE OF IDID NOT ALLOWED */
/*<             WRITE (LOUT,9070) IDID >*/
	    io___23.ciunit = *lout;
	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&(*idid), (ftnlen)sizeof(integer));
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          END IF >*/
	}
/*<          NN=N >*/
	nn = *n;
/*<          IF(ITOL.LE.3) NN = 1 >*/
	if (*itol <= 3) {
	    nn = 1;
	}
/*<          DO 1 I=1,NN >*/
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             IF (RTOL(I).LT.0.0D+0) THEN >*/
	    if (rtol[i__] < 0.) {
/*           ILLEGAL VALUE FOR RELATIVE ERROR TOLERENCE */
/*<                WRITE (LOUT,9040) >*/
		io___26.ciunit = *lout;
		s_wsfe(&io___26);
		e_wsfe();
/*<                IDID = -4 >*/
		*idid = -4;
/*<             END IF >*/
	    }
/*<  1       CONTINUE >*/
/* L1: */
	}
/*<          NN=N >*/
	nn = *n;
/*<          IF(ITOL.EQ.1.OR.ITOL.EQ.2.OR.ITOL.EQ.4) NN=1 >*/
	if (*itol == 1 || *itol == 2 || *itol == 4) {
	    nn = 1;
	}
/*<          DO 2 I=1,NN >*/
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             IF (ATOL(I).LT.0.0D+0) THEN >*/
	    if (atol[i__] < 0.) {
/*           ILLEGAL ABSOLUTE ERROR TOLERANCE */
/*<                WRITE(LOUT,9045) >*/
		io___27.ciunit = *lout;
		s_wsfe(&io___27);
		e_wsfe();
/*<                IDID=-4 >*/
		*idid = -4;
/*<             ENDIF >*/
	    }
/*<  2       CONTINUE >*/
/* L2: */
	}
/*<          IF(ITOL.EQ.1.AND.RTOL(1).EQ.0) THEN >*/
	if (*itol == 1 && rtol[1] == 0.) {
/*           ILLEGAL ERROR TOLERANCE */
/*<             WRITE(LOUT,9040) >*/
	    io___28.ciunit = *lout;
	    s_wsfe(&io___28);
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          ENDIF >*/
	}
/*<          IF(ITOL.NE.1) THEN >*/
	if (*itol != 1) {
/*<             VHOLD = 0.0D+0 >*/
	    vhold = 0.;
/*<             DO 3 I=1,N >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                IF(ITOL.EQ.2) THEN >*/
		if (*itol == 2) {
/*<                   VHOLD = DMAX1(RTOL(1),ATOL(1)) >*/
		    vhold = max(rtol[1],atol[1]);
/*<                ELSE IF (ITOL.EQ.3) THEN >*/
		} else if (*itol == 3) {
/*<                   VHOLD = DMAX1(RTOL(1),ATOL(I)) >*/
/* Computing MAX */
		    d__1 = rtol[1], d__2 = atol[i__];
		    vhold = max(d__1,d__2);
/*<                ELSE IF (ITOL.EQ.4) THEN >*/
		} else if (*itol == 4) {
/*<                   VHOLD = DMAX1(RTOL(I),ATOL(1)) >*/
/* Computing MAX */
		    d__1 = rtol[i__];
		    vhold = max(d__1,atol[1]);
/*<                ELSE IF (ITOL.EQ.5) THEN >*/
		} else if (*itol == 5) {
/*<                   VHOLD = DMAX1(RTOL(I),ATOL(I)) >*/
/* Computing MAX */
		    d__1 = rtol[i__], d__2 = atol[i__];
		    vhold = max(d__1,d__2);
/*<                ENDIF >*/
		}
/*<                IF(VHOLD.LE.0.0D+0) THEN >*/
		if (vhold <= 0.) {
/*<                   WRITE(LOUT,9040) >*/
		    io___30.ciunit = *lout;
		    s_wsfe(&io___30);
		    e_wsfe();
/*<                   IDID = -4 >*/
		    *idid = -4;
/*<                ENDIF >*/
		}
/*<  3          CONTINUE >*/
/* L3: */
	    }
/*<          ENDIF >*/
	}
/*<          IF (N.LE.0) THEN >*/
	if (*n <= 0) {
/*           ILLEGAL VALUE FOR THE NUMBER OF EQUATIONS */
/*<             WRITE (LOUT,9050) >*/
	    io___31.ciunit = *lout;
	    s_wsfe(&io___31);
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          END IF >*/
	}
/*<          IF ((T0-TOUT)*HO.GE.0.0D+0) THEN >*/
	if ((*t0 - *tout) * *ho >= 0.) {
/*           PARAMETERS FOR INTEGRATION ARE ILLEGAL */
/*<             WRITE (LOUT,9060) >*/
	    io___32.ciunit = *lout;
	    s_wsfe(&io___32);
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          END IF >*/
	}
/*<    >*/
	if (*mf != 21 && *mf != 22 && *mf != 23 && *mf != 24) {
/*           ILLEGAL VALUE FOR METHOD FLAG */
/*<             WRITE (LOUT,9090) MF >*/
	    io___33.ciunit = *lout;
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&(*mf), (ftnlen)sizeof(integer));
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          END IF >*/
	}
/*<          IF(ITOL.LT.1.OR.ITOL.GT.5) THEN >*/
	if (*itol < 1 || *itol > 5) {
/*           ILLEGAL VALUE FOR ERROR CONTROL PARAMETER */
/*<             WRITE (LOUT,9110) >*/
	    io___34.ciunit = *lout;
	    s_wsfe(&io___34);
	    e_wsfe();
/*<             IDID=-4 >*/
	    *idid = -4;
/*<          ENDIF >*/
	}
/*<          IF(MAXDER.LT.1.OR.MAXDER.GT.7) THEN >*/
	if (*maxder < 1 || *maxder > 7) {
/*        ILLEGAL VALUE FOR MAXIMUM ORDER */
/*<             WRITE(LOUT,9120) >*/
	    io___35.ciunit = *lout;
	    s_wsfe(&io___35);
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          ENDIF >*/
	}
/*<          IF(NIND1.EQ.0)  NIND1=N >*/
	if (*nind1 == 0) {
	    *nind1 = *n;
	}
/*<          IF(NIND1 + NIND2 + NIND3 .NE. N) THEN >*/
	if (*nind1 + *nind2 + *nind3 != *n) {
/*    SUM OF VARIABLES OF DIFFERENT INDEX SHOULD BE N. */
/*<             WRITE(LOUT,9140) >*/
	    io___36.ciunit = *lout;
	    s_wsfe(&io___36);
	    e_wsfe();
/*<             IDID = -4 >*/
	    *idid = -4;
/*<          ENDIF >*/
	}
/*<          IF (IDID.NE.1) THEN >*/
	if (*idid != 1) {
/*<             RETURN >*/
	    return 0;
/*<          ELSE >*/
	} else {
/*           THE INITIAL PARAMETERS ARE O.K. SO INITIALISE EVERYTHING */
/*           ELSE NECESSARY FOR THE INTEGRATION. */
/*           IF VALUES OF YMAX OTHER THAN THOSE SET BELOW ARE DESIRED, */
/*           THEY SHOULD BE SET HERE. ALL YMAX(I) MUST BE POSITIVE. IF */
/*           VALUES FOR HMIN OR HMAX, THE BOUNDS ON DABS(H), OTHER THAN */
/*           THOSE BELOW ARE DESIRED, THEY SHOULD BE SET BELOW. */
/*<             IF(ITOL.EQ.1) THEN >*/
	    if (*itol == 1) {
/*<                DO 10 I = 1,N >*/
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/*<                   YMAX(I) = DABS(Y0(I)) >*/
		    ymax[i__] = (d__1 = y0[i__], abs(d__1));
/*<                   YMAX(I)=DMAX1(YMAX(I),1.0D0) >*/
/* Computing MAX */
		    d__1 = ymax[i__];
		    ymax[i__] = max(d__1,1.);
/*<  10            CONTINUE >*/
/* L10: */
		}
/*<             ENDIF >*/
	    }
/*<             DO 15 I=1,N >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                Y(I,1)=Y0(I) >*/
		y[i__ + y_dim1] = y0[i__];
/*<  15         CONTINUE >*/
/* L15: */
	    }
/*<             T = T0 >*/
	    g_ovdriv_t = *t0;
/*<             H = HO >*/
	    g_h__ = *ho;
/*<             HMIN = DABS(HO) >*/
	    g_hmin = abs(*ho);
/*<             HMAX = DABS(T0-TEND)*10.0D+0 >*/
	    g_hmax = (d__1 = *t0 - *tend, abs(d__1)) * 10.;
/*<             JSTART = 0 >*/
	    g_jstart = 0;
/*<             NHCUT = 0 >*/
	    g_nhcut = 0;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*     <<<<<<<<<<<<<<<<< */
/*     <  TAKE A STEP  > */
/*     <<<<<<<<<<<<<<<<< */
/*<  20   IF ((T+H).EQ.T) THEN >*/
L20:
    if (g_ovdriv_t + g_h__ == g_ovdriv_t) {
/*<          WRITE (LOUT,9000) >*/
	io___39.ciunit = *lout;
	s_wsfe(&io___39);
	e_wsfe();
/*<       END IF >*/
    }
/*<    >*/
    stiff_(&g_h__, &g_hmax, &g_hmin, &g_jstart, &g_kflag, mf, &mbnd[1], nind1, nind2, 
	    nind3, &g_ovdriv_t, tout, tend, &y[y_offset], &yprime[1], n, &ymax[1], &
	    errors[1], &save1[1], &save2[1], &scale[1], &pw[1], &pwcopy[1], &
	    yhold[yhold_offset], &ynhold[ynhold_offset], &arh[1], &ipiv[1], 
	    lout, maxder, itol, &rtol[1], &atol[1], &rpar[1], &ipar[1], (U_fp)
	    pderv, (U_fp)resid, nqused, nstep, nfail, nre, nje, ndec, nbsol, 
	    npset, ncoset, maxord, maxstp, uround, epsjac, hused, ierr);
/*      IF(IERR.NE.0) THEN */
/*      WRITE(LOUT,9150) */
/*      RETURN */
/*      ENDIF */
/*<       KGO = 1 - KFLAG >*/
    kgo = 1 - g_kflag;
/*<       IF (KGO.EQ.1) THEN >*/
    if (kgo == 1) {
/*        NORMAL RETURN FROM STIFF */
/*<          GO TO 30 >*/
	goto L30;
/*<       ELSE IF (KGO.EQ.2) THEN >*/
    } else if (kgo == 2) {
/*        COULD NOT ACHIEVE REQUIRED PRECISION WITH HMIN */
/*        SO CHOP HMIN IF WE HAVEN'T DONE SO 10 TIMES */
/*<          GO TO 60 >*/
	goto L60;
/*<       ELSE IF (KGO.EQ.3) THEN >*/
    } else if (kgo == 3) {
/*        ERROR REQUIREMENT SMALLER THAN CAN BE HANDLED FOR THIS PROBLEM */
/*<          WRITE (LOUT,9010) T,H >*/
	io___41.ciunit = *lout;
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&g_ovdriv_t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g_h__, (ftnlen)sizeof(doublereal));
	e_wsfe();
/*<          GO TO 70 >*/
	goto L70;
/*<       ELSE IF (KGO.EQ.4) THEN >*/
    } else if (kgo == 4) {
/*        COULD NOT ACHIEVE CONVERGENCE WITH HMIN */
/*<          WRITE (LOUT,9030) T >*/
	io___42.ciunit = *lout;
	s_wsfe(&io___42);
	do_fio(&c__1, (char *)&g_ovdriv_t, (ftnlen)sizeof(doublereal));
	e_wsfe();
/*<          GO TO 60 >*/
	goto L60;
/*<       END IF >*/
    }
/*<  30   CONTINUE >*/
L30:
/* --------------------------------------------------------------------- */
/*     NORMAL RETURN FROM THE INTEGRATOR. */

/*     THE WEIGHTS YMAX(I) ARE UPDATED IF ITOL=1. */
/*     IF DIFFERENT VALUES ARE DESIRED, THEY SHOULD BE SET HERE. */

/*     ANY OTHER TESTS OR CALCULATIONS THAT ARE REQUIRED AFTER EVERY */
/*     STEP SHOULD BE INSERTED HERE. */

/*     IF IDID = 3, Y0 IS SET TO THE CURRENT Y VALUES ON RETURN. */
/*     IF IDID = 2, H IS CONTROLLED TO HIT TOUT (WITHIN ROUNDOFF */
/*     ERROR), AND THEN THE CURRENT Y VALUES ARE PUT IN Y0 ON RETURN. */
/*     FOR ANY OTHER VALUE OF IDID, CONTROL RETURNS TO THE INTEGRATOR */
/*     UNLESS TOUT HAS BEEN REACHED.  THEN INTERPOLATED VALUES OF Y ARE */
/*     COMPUTED AND STORED IN Y0 ON RETURN. */
/*     IF INTERPOLATION IS NOT DESIRED, THE CALL TO INTERP SHOULD BE */
/*     REMOVED AND CONTROL TRANSFERRED TO STATEMENT 500 INSTEAD OF 520. */
/* -------------------------------------------------------------------- */
/*<       IF(NSTEP.GT.MAXSTP) THEN >*/
    if (*nstep > *maxstp) {
/*<          KGO=5 >*/
	kgo = 5;
/*<          KLAG=4 >*/
	klag = 4;
/*   TOO MUCH WORK */
/*<          WRITE(LOUT,9130) >*/
	io___44.ciunit = *lout;
	s_wsfe(&io___44);
	e_wsfe();
/*<          IDID = -6 >*/
	*idid = -6;
/*<          GOTO 70 >*/
	goto L70;
/*<       END IF >*/
    }
/*<       IF(ITOL.EQ.1) THEN >*/
    if (*itol == 1) {
/*<          D = 0.0D+0 >*/
	d__ = 0.;
/*<          DO 40 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AYI = DABS(Y(I,1)) >*/
	    ayi = (d__1 = y[i__ + y_dim1], abs(d__1));
/*<             YMAX(I)=DMAX1(YMAX(I),AYI) >*/
/* Computing MAX */
	    d__1 = ymax[i__];
	    ymax[i__] = max(d__1,ayi);
/*<  40      CONTINUE >*/
/* L40: */
	}
/*<       ENDIF >*/
    }
/*<       IF (IDID.EQ.3.OR.IDID.EQ.1) GO TO 70 >*/
    if (*idid == 3 || *idid == 1) {
	goto L70;
    }
/*<       IF (DABS(T-TOUT).LE.DABS(15.0D+0*UROUND*TOUT)) THEN >*/
    if ((d__1 = g_ovdriv_t - *tout, abs(d__1)) <= (d__2 = *uround * 15. * *tout, abs(
	    d__2))) {
/*        EFFECTIVELY WE HAVE HIT TOUT */
/*<          IDID = KFLAG >*/
	*idid = g_kflag;
/*<          T0 = TOUT >*/
	*t0 = *tout;
/*<          DO 50 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             Y0(I) = Y(I,1) >*/
	    y0[i__] = y[i__ + y_dim1];
/*<  50      CONTINUE >*/
/* L50: */
	}
/*<          HO = H >*/
	*ho = g_h__;
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       IF (IDID.EQ.2) THEN >*/
    if (*idid == 2) {
/*        CONTINUING INTEGRATION BUT MUST HIT TOUT EXACTLY */
/*<          IF (((T+H)-TOUT)*H.GT.0.0D+0) THEN >*/
	if ((g_ovdriv_t + g_h__ - *tout) * g_h__ > 0.) {
/*           WE HAVE ALREADY OVERSHOT THE OUTPUT POINT OR WE WILL DO */
/*           SO ON THE NEXT STEP */
/*<    >*/
	    if ((g_ovdriv_t - *tout) * g_h__ >= 0. || (d__1 = g_ovdriv_t - *tout, abs(d__1)) <= *
		    uround * 100. * g_hmax) {
/*              HAVE OVERSHOT, SO INTERPOLATE */
/*<                CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
		interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<                T0 = TOUT >*/
		*t0 = *tout;
/*<                HO = H >*/
		*ho = g_h__;
/*<                IDID = KFLAG >*/
		*idid = g_kflag;
/*<                RETURN >*/
		return 0;
/*<             ELSE >*/
	    } else {
/*              WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE */
/*              SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY' */
/*<                H = (TOUT-T)* (1.0D+0-4.0D+0*UROUND) >*/
		g_h__ = (*tout - g_ovdriv_t) * (1. - *uround * 4.);
/*<                JSTART = -1 >*/
		g_jstart = -1;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF ((T-TOUT)*H.GE.0.0D+0) THEN >*/
    } else if ((g_ovdriv_t - *tout) * g_h__ >= 0.) {
/*        HAVE OVERSHOT, SO INTERPOLATE */
/*<          CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
	interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<          IDID = KFLAG >*/
	*idid = g_kflag;
/*<          HO = H >*/
	*ho = g_h__;
/*<          T0 = TOUT >*/
	*t0 = *tout;
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       GO TO 20 >*/
    goto L20;
/* ------------------------------------------------------------------- */
/*     ON AN ERROR RETURN FROM THE INTEGRATOR, AN IMMEDIATE RETURN OCCURS */
/*     IF KFLAG = -2, AND RECOVERY ATTEMPTS ARE MADE OTHERWISE. */
/*     H AND HMIN ARE REDUCED BY A FACTOR OF .1 UP TO 10 TIMES */
/*     BEFORE GIVING UP. */
/* -------------------------------------------------------------------- */
/*<  60   CONTINUE >*/
L60:
/*<       IF (NHCUT.EQ.10) THEN >*/
    if (g_nhcut == 10) {
/*        HAVE REDUCED H TEN TIMES */
/*<          WRITE (LOUT,9100) >*/
	io___47.ciunit = *lout;
	s_wsfe(&io___47);
	e_wsfe();
/*<          GO TO 70 >*/
	goto L70;
/*<       END IF >*/
    }
/*<       NHCUT = NHCUT + 1 >*/
    ++g_nhcut;
/*<       HMIN = 0.1D+0*HMIN >*/
    g_hmin *= .1;
/*<       H = 0.1D+0*H >*/
    g_h__ *= .1;
/*<       JSTART = -1 >*/
    g_jstart = -1;
/*<       GO TO 20 >*/
    goto L20;
/*<  70   IF(DABS(T-TOUT).GT.1000.0D+0*UROUND) THEN >*/
L70:
    if ((g_ovdriv_t - *tout) <  *uround * 1e3) { /* grn: removed abs */
/*     if ((d__1 = g_ovdriv_t - *tout, abs(d__1)) >  *uround * 1e3) {  */ 
/*<          DO 80 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             Y0(I) = Y(I,1) >*/
	    y0[i__] = y[i__ + y_dim1];
/*<  80      CONTINUE >*/
/* L80: */
	}
/*<          T0 = T >*/
	*t0 = g_ovdriv_t;
/*<       ELSE >*/
    } else {
/*        HAVE PASSED TOUT SO INTERPOLATE */
/*<          CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
	interp_(n, &g_jstart, &g_h__, &g_ovdriv_t, &y[y_offset], tout, &y0[1]);
/*<          T0 = TOUT >*/
	*t0 = *tout;
/*<          IDID = KFLAG >*/
	*idid = g_kflag;
/*<       END IF >*/
    }
/*<       HO = H >*/
    *ho = g_h__;
/*<       IF(KFLAG.NE.0) IDID = KFLAG >*/
    if (g_kflag != 0) {
	*idid = g_kflag;
    }
/*<       RETURN >*/
    return 0;
/* -------------------------- END OF SUBROUTINE OVDRIV ----------------- */
/*<  9000 FORMAT (' WARNING..  T + H = T ON NEXT STEP.') >*/
/*<  9 >*/
/*<  9 >*/
/* L9020: */
/*<  9 >*/
/*<  9040 FORMAT (/,/,' ILLEGAL INPUT.. RTOL .LE. 0.',/,/) >*/
/*<  9045 FORMAT (/,/,' ILLEGAL INPUT.. ATOL .LE. 0.',/,/) >*/
/*<  9050 FORMAT (/,/,' ILLEGAL INPUT.. N .LE. 0',/,/) >*/
/*<  9060 FORMAT (/,/,' ILLEGAL INPUT.. (T0-TOUT)*H .GE. 0.',/,/) >*/
/*<  9070 FORMAT (/,/,' ILLEGAL INPUT.. IDID =',I5,/,/) >*/
/*<  9 >*/
/*<  9 >*/
/*<  9 >*/
/*<  9110 FORMAT (/,/,' ILLEGAL VALUE FOR ITOL',/,/) >*/
/*<  9120 FORMAT (/,/,' ILLEGAL VALUE FOR MAXDER',/,/) >*/
/*<  9130 FORMAT (/,/,' NUMBER OF STEPS EXCEEDS MAXIMUM',/,/) >*/
/*<  9 >*/
/*<  9 >*/
/* L9150: */
/*<  9160 FORMAT (/,/,'STEPSIZE IS TOO SMALL') >*/
/*<       END >*/
} /* ovdriv_ */

/* -------------------------------------------------------------------------- */

/*<       SUBROUTINE INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
/* Subroutine */ int interp_(integer *n, integer *jstart, doublereal *h__, 
	doublereal *t, doublereal *y, doublereal *tout, doublereal *y0)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    /*static*/ integer i__, j, l;
    /*static*/ doublereal s, s1;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER JSTART,N >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  Y(N,12),Y0(N) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,J,L >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*     .. */
/*<       DO 10 I = 1,N >*/
    /* Parameter adjustments */
    --y0;
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y0(I) = Y(I,1) >*/
	y0[i__] = y[i__ + y_dim1];
/*<  10   CONTINUE >*/
/* L10: */
    }
/*<       L = JSTART + 2 >*/
    l = *jstart + 2;
/*<       S = (TOUT-T)/H >*/
    s = (*tout - *t) / *h__;
/*<       S1 = 1.0D+0 >*/
    s1 = 1.;
/*<       DO 30 J = 2,L >*/
    i__1 = l;
    for (j = 2; j <= i__1; ++j) {
/*<          S1 = S1* (S+DBLE(FLOAT(J-2)))/DBLE(FLOAT(J-1)) >*/
	s1 = s1 * (s + (doublereal) (j - 2)) / (doublereal) (j - 1);
/*<          DO 20 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             Y0(I) = Y0(I) + S1*Y(I,J) >*/
	    y0[i__] += s1 * y[i__ + j * y_dim1];
/*<  20      CONTINUE >*/
/* L20: */
	}
/*<  30   CONTINUE >*/
/* L30: */
    }
/*<       RETURN >*/
    return 0;
/* -------------- END OF SUBROUTINE INTERP --------------------------- */
/*<       END >*/
} /* interp_ */

/*<       SUBROUTINE COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD) >*/
/* Subroutine */ int coset_(integer *nq, doublereal *el, doublereal *elst, 
	doublereal *tq, integer *ncoset, integer *maxord)
{
    /* Initialized data */

	/* only used, not assigned */
    static doublereal pertst[24]	/* was [8][3] */ = { 1.,2.,4.5,7.333,
	    10.42,13.7,17.15,20.74,2.,4.5,7.333,10.42,13.7,17.15,20.74,24.46,
	    4.5,7.333,10.42,13.7,17.15,20.74,24.46,1. };

	/* used as a local counter */
    /*static*/ integer k;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/* -------------------------------------------------------------------- */
/*     COSET IS CALLED BY THE INTEGRATOR AND SETS THE COEFFICIENTS USED */
/*     BY THE CONVENTIONAL BACKWARD DIFFERENTIATION SCHEME AND THE */
/*     MODIFIED EXTENDED BACKWARD DIFFERENTIATION SCHEME.  THE VECTOR */
/*     EL OF LENGTH NQ+1 DETERMINES THE BASIC BDF METHOD WHILE THE VECTOR */
/*     ELST OF LENGTH NQ+2 DETERMINES THE MEBDF.  THE VECTOR TQ OF */
/*     LENGTH 4 IS INVOLVED IN ADJUSTING THE STEPSIZE IN RELATION TO THE */
/*     TRUNCATION ERROR.  ITS VALUES ARE GIVEN BY THE PERTST ARRAY.  THE */
/*     VECTORS EL AND TQ BOTH DEPEND ON METH AND NQ.  THE */
/*     COEFFICIENTS IN PERTST NEED TO BE GIVEN TO ONLY ABOUT ONE PERCENT */
/*     ACCURACY.  THE ORDER IN WHICH THE GROUPS APPEAR BELOW IS: */
/*     COEFFICIENTS FOR ORDER NQ-1, COEFFICIENTS FOR ORDER NQ, */
/*     COEFFICIENTS FOR ORDER NQ+1. */
/* ------------------------------------------------------------------- */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER NQ >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  EL(10),ELST(10),TQ(5) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER K >*/
/*     .. */
/*     .. LOCAL ARRAYS .. */
/*<       DIMENSION  PERTST(8,3) >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*     .. */
/*     .. COMMON BLOCKS .. */
/*<       INTEGER NCOSET,MAXORD >*/
/*     .. */
/*     .. DATA STATEMENTS .. */
/*<    >*/
    /* Parameter adjustments */
    --tq;
    --elst;
    --el;

    /* Function Body */
/*<    >*/
/*<    >*/
/*     .. */
/* ------------------------------------------------------------------- */
/*     THE FOLLOWING COEFFICIENTS SHOULD BE DEFINED TO MACHINE ACCURACY. */
/*     THEIR DERIVATION IS GIVEN IN REFERENCE 1. */
/* ------------------------------------------------------------------- */
/*<       IF (NQ.GT.MAXORD) MAXORD = NQ >*/
    if (*nq > *maxord) {
	*maxord = *nq;
    }
/*<       NCOSET = NCOSET + 1 >*/
    ++(*ncoset);
/*<       GO TO (10,20,30,40,50,60,70) NQ >*/
    switch (*nq) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L50;
	case 6:  goto L60;
	case 7:  goto L70;
    }
/*<    10 EL(1) = 1.0D+0 >*/
L10:
    el[1] = 1.;
/*<       ELST(1) = 1.5D+0 >*/
    elst[1] = 1.5;
/*<       ELST(3) = -0.5D+0 >*/
    elst[3] = -.5;
/*<       GO TO 80 >*/
    goto L80;
/*<    20 EL(1) = 6.6666666666667D-01 >*/
L20:
    el[1] = .66666666666667;
/*<       EL(3) = 3.3333333333333D-01 >*/
    el[3] = .33333333333333;
/*<       ELST(1) = 9.5652173913043D-01 >*/
    elst[1] = .95652173913043;
/*<       ELST(3) = 2.1739130434782D-01 >*/
    elst[3] = .21739130434782;
/*<       ELST(4) = -1.7391304347826D-01 >*/
    elst[4] = -.17391304347826;
/*<       GO TO 80 >*/
    goto L80;
/*<    30 EL(1) = 5.4545454545455D-01 >*/
L30:
    el[1] = .54545454545455;
/*<       EL(3) = 4.5454545454545D-01 >*/
    el[3] = .45454545454545;
/*<       EL(4) = 1.8181818181818D-01 >*/
    el[4] = .18181818181818;
/*<       ELST(1) = 7.6142131979695D-01 >*/
    elst[1] = .76142131979695;
/*<       ELST(3) = 3.2994923857868D-01 >*/
    elst[3] = .32994923857868;
/*<       ELST(4) = 8.6294416243654D-02 >*/
    elst[4] = .086294416243654;
/*<       ELST(5) = -9.1370558375634D-02 >*/
    elst[5] = -.091370558375634;
/*<       GO TO 80 >*/
    goto L80;
/*<    40 EL(1) = 0.48D+0 >*/
L40:
    el[1] = .48;
/*<       EL(3) = 0.52D+0 >*/
    el[3] = .52;
/*<       EL(4) = 0.28D+0 >*/
    el[4] = .28;
/*<       EL(5) = 0.12D+0 >*/
    el[5] = .12;
/*<       ELST(1) = 6.5733706517393D-01 >*/
    elst[1] = .65733706517393;
/*<       ELST(3) = 4.0023990403838D-01 >*/
    elst[3] = .40023990403838;
/*<       ELST(4) = 1.5793682526989D-01 >*/
    elst[4] = .15793682526989;
/*<       ELST(5) = 4.4382247101159D-02 >*/
    elst[5] = .044382247101159;
/*<       ELST(6) = -5.7576969212315D-02 >*/
    elst[6] = -.057576969212315;
/*<       GO TO 80 >*/
    goto L80;
/*<    50 EL(1) = 4.3795620437956D-01 >*/
L50:
    el[1] = .43795620437956;
/*<       EL(3) = 5.62043795620436D-01 >*/
    el[3] = .562043795620436;
/*<       EL(4) = 3.43065693430656D-01 >*/
    el[4] = .343065693430656;
/*<       EL(5) = 1.97080291970802D-01 >*/
    el[5] = .197080291970802;
/*<       EL(6) = 8.75912408759123D-02 >*/
    el[6] = .0875912408759123;
/*<       ELST(1) = 5.9119243917152D-01 >*/
    elst[1] = .59119243917152;
/*<       ELST(3) = 4.4902473356122D-01 >*/
    elst[3] = .44902473356122;
/*<       ELST(4) = 2.1375427307460D-01 >*/
    elst[4] = .2137542730746;
/*<       ELST(5) = 9.0421610027481503D-02 >*/
    elst[5] = .090421610027481503;
/*<       ELST(6) = 2.6409276761177D-02 >*/
    elst[6] = .026409276761177;
/*<       ELST(7) = -4.0217172732757D-02 >*/
    elst[7] = -.040217172732757;
/*<       GO TO 80 >*/
    goto L80;
/*<    60 EL(1) = 4.08163265306120D-01 >*/
L60:
    el[1] = .40816326530612;
/*<       EL(3) = 5.91836734693874D-01 >*/
    el[3] = .591836734693874;
/*<       EL(4) = 3.87755102040813D-01 >*/
    el[4] = .387755102040813;
/*<       EL(5) = 2.51700680272107D-01 >*/
    el[5] = .251700680272107;
/*<       EL(6) = 1.49659863945577D-01 >*/
    el[6] = .149659863945577;
/*<       EL(7) = 6.80272108843534D-02 >*/
    el[7] = .0680272108843534;
/*<       ELST(1) = 5.4475876041119D-01 >*/
    elst[1] = .54475876041119;
/*<       ELST(3) = 4.8525549636077D-01 >*/
    elst[3] = .48525549636077;
/*<       ELST(4) = 2.5789750131312D-01 >*/
    elst[4] = .25789750131312;
/*<       ELST(5) = 1.3133738525800D-01 >*/
    elst[5] = .131337385258;
/*<       ELST(6) = 5.7677396763462D-02 >*/
    elst[6] = .057677396763462;
/*<       ELST(7) = 1.7258197643881D-02 >*/
    elst[7] = .017258197643881;
/*<       ELST(8) = -3.0014256771967D-02 >*/
    elst[8] = -.030014256771967;
/*<       GO TO 80 >*/
    goto L80;
/*<    70 EL(1) = 3.85674931129476D-01 >*/
L70:
    el[1] = .385674931129476;
/*<       EL(3) = 6.14325068870521D-01 >*/
    el[3] = .614325068870521;
/*<       EL(4) = 4.21487603305783D-01 >*/
    el[4] = .421487603305783;
/*<       EL(5) = 2.9292929292929D-01 >*/
    el[5] = .29292929292929;
/*<       EL(6) = 1.96510560146923D-01 >*/
    el[6] = .196510560146923;
/*<       EL(7) = 1.19375573921028D-01 >*/
    el[7] = .119375573921028;
/*<       EL(8) = 5.50964187327820D-02 >*/
    el[8] = .055096418732782;
/*<       ELST(1) = 5.0999746293734D-01 >*/
    elst[1] = .50999746293734;
/*<       ELST(3) = 5.1345839935281D-01 >*/
    elst[3] = .51345839935281;
/*<       ELST(4) = 2.9364346131937D-01 >*/
    elst[4] = .29364346131937;
/*<       ELST(5) = 1.6664672120553D-01 >*/
    elst[5] = .16664672120553;
/*<       ELST(6) = 8.8013735242353D-02 >*/
    elst[6] = .088013735242353;
/*<       ELST(7) = 3.9571794884069D-02 >*/
    elst[7] = .039571794884069;
/*<       ELST(8) = 1.2039080338722D-02 >*/
    elst[8] = .012039080338722;
/*<       ELST(9) = -2.3455862290154D-02 >*/
    elst[9] = -.023455862290154;
/*<    80 DO 90 K = 1,3 >*/
L80:
    for (k = 1; k <= 3; ++k) {
/*<         TQ(K) = PERTST(NQ,K) >*/
	tq[k] = pertst[*nq + (k << 3) - 9];
/*<    90 CONTINUE >*/
/* L90: */
    }
/*<       TQ(4) = 0.5D+0*TQ(2)/DBLE(FLOAT(NQ)) >*/
    tq[4] = tq[2] * .5 / (doublereal) (*nq);
/*<       IF(NQ.NE.1) TQ(5)=PERTST(NQ-1,1) >*/
    if (*nq != 1) {
	tq[5] = pertst[*nq - 2];
    }
/*<       RETURN >*/
    return 0;
/* --------------------- END OF SUBROUTINE COSET --------------------- */
/*<       END >*/
} /* coset_ */



/*<    >*/
/* Subroutine */ int pset_(doublereal *y, doublereal *yprime, integer *n, 
	doublereal *h__, doublereal *t, doublereal *uround, doublereal *
	epsjac, doublereal *con, integer *miter, integer *mbnd, integer *
	nind1, integer *nind2, integer *nind3, integer *ier, S_fp pderv, S_fp 
	resid, integer *nrenew, doublereal *ymax, doublereal *save1, 
	doublereal *save2, doublereal *save3, doublereal *pw, doublereal *
	pwcopy, doublereal *wrkspc, integer *ipiv, integer *itol, doublereal *
	rtol, doublereal *atol, integer *npset, integer *nje, integer *nre, 
	integer *ndec, integer *ipar, doublereal *rpar, integer *ierr)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

	/* grn assume all variables here local, since it is in the comment! */

    /* Local variables */
    /*static*/ integer jjkk, itmb;
    extern /* Subroutine */ int dgbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

    /*static*/ doublereal d__;
    /*static*/ integer i__, j;
    /*static*/ doublereal r__;
    /*static*/ integer i1, j1, i2;
    /*static*/ doublereal r0;
    /*static*/ integer ii, jj, ml, mu;
    /*static*/ doublereal yj, yi, yp, tempry;
    extern /* Subroutine */ int dec_(integer *, integer *, doublereal *, 
	    integer *, integer *);
    /*static*/ integer mba;
    /*static*/ doublereal yjj, yjp;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/* ------------------------------------------------------------------- */
/*     PSET IS CALLED BY STIFF TO COMPUTE AND PROCESS THE MATRIX */
/*     PD=DG/DY + (1/CON)DG/DY'. THIS MATRIX IS THEN SUBJECTED TO LU */
/*     DECOMPOSITION IN PREPARATION FOR LATER SOLUTION OF LINEAR SYSTEMS */
/*     OF ALGEBRAIC EQUATIONS WITH LU AS THE COEFFICIENT MATRIX.  THE */
/*     MATRIX PD IS FOUND BY THE USER-SUPPLIED ROUTINE PDERV IF MITER=1 */
/*     OR 3 OR BY FINITE DIFFERENCING IF MITER = 2 OR 4. */
/*     IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION WITH */
/*     PSET USES THE FOLLOWING .. */
/*     EPSJAC = DSQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS. */
/* ******************************************************************* */
/*     THE ARGUMENT NRENEW IS USED TO SIGNAL WHETHER OR NOT */
/*     WE REQUIRE A NEW JACOBIAN TO BE CALCULATED. */

/*        IF NRENEW > 0 THEN WE REQUIRE A NEW J TO BE COMPUTED */
/*                  = 0 THEN USE A COPY OF THE LAST J COMPUTED */
/* ******************************************************************* */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER IER,MITER,N,NRENEW >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<    >*/
/*<       INTEGER IPIV(N), MBND(4) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,J,J1,JJKK,FOUR,FIVE >*/
/*     .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*<       EXTERNAL DEC,PDERV,DGBFA,RESID >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*<       INTRINSIC DABS,DMAX1,DSQRT >*/
/*     .. */
/*     .. COMMON BLOCKS .. */
/*<       INTEGER NDEC,NRE,NJE,NPSET >*/
/*<       NPSET = NPSET + 1 >*/
    /* Parameter adjustments */
    --ipiv;
    --wrkspc;
    --save3;
    --save2;
    --save1;
    --ymax;
    --yprime;
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    --mbnd;
    --pw;
    --pwcopy;
    --rtol;
    --atol;
    --ipar;
    --rpar;

    /* Function Body */
    ++(*npset);
/*<       ML = MBND(1) >*/
    ml = mbnd[1];
/*<       MU = MBND(2) >*/
    mu = mbnd[2];
/*<       IF (NRENEW.EQ.0) THEN >*/
    if (*nrenew == 0) {
/*<          IF (MITER.LT.3) THEN >*/
	if (*miter < 3) {
/*<             DO 10 I = 1,N*N >*/
	    i__1 = *n * *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                PW(I) = PWCOPY(I) >*/
		pw[i__] = pwcopy[i__];
/*<  10         CONTINUE >*/
/* L10: */
	    }
/*<          ELSE >*/
	} else {
/*<             DO 15 I=0,N-1 >*/
	    i__1 = *n - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
/*<                DO 18 J=MBND(1)+1,MBND(4) >*/
		i__2 = mbnd[4];
		for (j = mbnd[1] + 1; j <= i__2; ++j) {
/*<                   PW(I*MBND(4)+J) = PWCOPY(I*MBND(4)+J) >*/
		    pw[i__ * mbnd[4] + j] = pwcopy[i__ * mbnd[4] + j];
/*<  18            CONTINUE >*/
/* L18: */
		}
/*<  15         CONTINUE >*/
/* L15: */
	    }
/*<          ENDIF >*/
	}
/*<          GO TO 70 >*/
	goto L70;
/*<       ENDIF >*/
    }

/*<       IF (MITER.EQ.2.OR.MITER.EQ.4) GO TO 30 >*/
    if (*miter == 2 || *miter == 4) {
	goto L30;
    }

/*<       NJE = NJE + 1 >*/
    ++(*nje);
/*<       IF ( MITER.NE.3 ) THEN >*/
    if (*miter != 3) {
/*<          CALL PDERV(T,Y,PWCOPY,N,YPRIME,N,CON,IPAR,RPAR,IERR) >*/
	(*pderv)(t, &y[y_offset], &pwcopy[1], n, &yprime[1], n, con, &ipar[1],
		 &rpar[1], ierr);
/*<          DO 20 I=1,N*N >*/
	i__1 = *n * *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             PW(I)=PWCOPY(I) >*/
	    pw[i__] = pwcopy[i__];
/*<  20      CONTINUE >*/
/* L20: */
	}
/*<       ELSE >*/
    } else {
/*<          CALL PDERV(T,Y,PWCOPY,N,YPRIME,MBND(4),CON,IPAR,RPAR,IERR) >*/
	(*pderv)(t, &y[y_offset], &pwcopy[1], n, &yprime[1], &mbnd[4], con, &
		ipar[1], &rpar[1], ierr);
/*<          DO 25 I=0,N-1 >*/
	i__1 = *n - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/*<             ITMB=I*MBND(4) >*/
	    itmb = i__ * mbnd[4];
/*<             DO 28 J=MBND(1)+1,MBND(4) >*/
	    i__2 = mbnd[4];
	    for (j = mbnd[1] + 1; j <= i__2; ++j) {
/*<                PW(ITMB+J) = PWCOPY(ITMB+J-MBND(1)) >*/
		pw[itmb + j] = pwcopy[itmb + j - mbnd[1]];
/*<  28         CONTINUE >*/
/* L28: */
	    }
/*<  25      CONTINUE >*/
/* L25: */
	}
/*<       ENDIF >*/
    }
/*<       GOTO 70 >*/
    goto L70;

/*<  30   CONTINUE       >*/
L30:
/*<       NJE = NJE + 1 >*/
    ++(*nje);
/*<       DO 35 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF(ITOL.EQ.2) THEN >*/
	if (*itol == 2) {
/*<             YMAX(I)=Y(I,1)*RTOL(1)+ATOL(1) >*/
	    ymax[i__] = y[i__ + y_dim1] * rtol[1] + atol[1];
/*<          ELSE IF (ITOL.EQ.3) THEN >*/
	} else if (*itol == 3) {
/*<             YMAX(I)=Y(I,1)*RTOL(1)+ATOL(I) >*/
	    ymax[i__] = y[i__ + y_dim1] * rtol[1] + atol[i__];
/*<          ELSE IF(ITOL.EQ.4) THEN >*/
	} else if (*itol == 4) {
/*<             YMAX(I)=Y(I,1)*RTOL(I)+ATOL(1) >*/
	    ymax[i__] = y[i__ + y_dim1] * rtol[i__] + atol[1];
/*<          ELSE IF(ITOL.EQ.5) THEN >*/
	} else if (*itol == 5) {
/*<             YMAX(I)=Y(I,1)*RTOL(I)+ATOL(I) >*/
	    ymax[i__] = y[i__ + y_dim1] * rtol[i__] + atol[i__];
/*<          END IF >*/
	}
/*<  35   CONTINUE >*/
/* L35: */
    }
/*<       D = 0.0D+0 >*/
    d__ = 0.;
/*<       DO 40 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          D = D + (YPRIME(I)*(YMAX(I)))**2 >*/
/* Computing 2nd power */
	d__1 = yprime[i__] * ymax[i__];
	d__ += d__1 * d__1;
/*<  40   CONTINUE >*/
/* L40: */
    }
/*<       IF(ITOL.EQ.1) D=D*RTOL(1)**2 >*/
    if (*itol == 1) {
/* Computing 2nd power */
	d__1 = rtol[1];
	d__ *= d__1 * d__1;
    }
/*<       D = DSQRT(D)/DFLOAT(N) >*/
    d__ = sqrt(d__) / (doublereal) (*n);
/*<       R0 = DABS(H)*D*DFLOAT(N)*(1.0D+03)*UROUND >*/
    r0 = abs(*h__) * d__ * (doublereal) (*n) * 1e3 * *uround;
/*<       IF(R0.EQ.0.0D+0) R0 = epsjac >*/
    if (r0 == 0.) {
	r0 = *epsjac;
    }
/*<       J1 = 0 >*/
    j1 = 0;
/*<       IF(MITER.EQ.4) GOTO 51       >*/
    if (*miter == 4) {
	goto L51;
    }
/*<       CALL RESID(N,T,Y,SAVE2,YPRIME,IPAR,RPAR,IERR)  >*/
    (*resid)(n, t, &y[y_offset], &save2[1], &yprime[1], &ipar[1], &rpar[1], 
	    ierr);
/*<       NRE=NRE+1 >*/
    ++(*nre);
/*<       DO 60 J = 1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          YJ = Y(J,1) >*/
	yj = y[j + y_dim1];
/*<          YP=YPRIME(J)         >*/
	yp = yprime[j];
/*<          IF(ITOL.EQ.1) THEN >*/
	if (*itol == 1) {
/*<             R=DMAX1(EPSJAC*DABS(YJ),R0/(YMAX(J)*RTOL(1))) >*/
/* Computing MAX */
	    d__1 = *epsjac * abs(yj), d__2 = r0 / (ymax[j] * rtol[1]);
	    r__ = max(d__1,d__2);
/*<          ELSE >*/
	} else {
/*<             R=DMAX1(EPSJAC*DABS(YJ),R0/YMAX(J)) >*/
/* Computing MAX */
	    d__1 = *epsjac * abs(yj), d__2 = r0 / ymax[j];
	    r__ = max(d__1,d__2);
/*<          ENDIF >*/
	}
/*<          R = epsjac >*/
	r__ = *epsjac;
/*<          Y(J,1) = Y(J,1) + R >*/
	y[j + y_dim1] += r__;
/*<          YPRIME(J)=YPRIME(J)+R/CON                            >*/
	yprime[j] += r__ / *con;
/*<          CALL RESID(N,T,Y,WRKSPC,YPRIME,IPAR,RPAR,IERR)         >*/
	(*resid)(n, t, &y[y_offset], &wrkspc[1], &yprime[1], &ipar[1], &rpar[
		1], ierr);
/*<          DO 50 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             JJKK = I + J1 >*/
	    jjkk = i__ + j1;
/*<             TEMPRY = (WRKSPC(I)-SAVE2(I)) >*/
	    tempry = wrkspc[i__] - save2[i__];
/*<             PWCOPY(JJKK) = TEMPRY/R >*/
	    pwcopy[jjkk] = tempry / r__;
/*<             PW(JJKK) =  PWCOPY(JJKK)           >*/
	    pw[jjkk] = pwcopy[jjkk];
/*<  50      CONTINUE                 >*/
/* L50: */
	}
/*<          Y(J,1) = YJ >*/
	y[j + y_dim1] = yj;
/*<          YPRIME(J)=YP >*/
	yprime[j] = yp;
/*<          J1 = J1 + N >*/
	j1 += *n;
/*<  60   CONTINUE >*/
/* L60: */
    }
/*<       NRE = NRE + N  >*/
    *nre += *n;
/*<       GOTO 70 >*/
    goto L70;

/*< 51    CONTINUE >*/
L51:
/*<       CALL RESID(N,T,Y,SAVE2,YPRIME,IPAR,RPAR,IERR) >*/
    (*resid)(n, t, &y[y_offset], &save2[1], &yprime[1], &ipar[1], &rpar[1], 
	    ierr);
/*<       NRE = NRE+1 >*/
    ++(*nre);
/*<       MBA = min0(MBND(3),N) >*/
    mba = min(mbnd[3],*n);
/*<       DO 61 J=1,MBA >*/
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
/*<          DO 161 I=J,N,MBND(3) >*/
	i__2 = *n;
	i__3 = mbnd[3];
	for (i__ = j; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
/*<             SAVE1(I) = Y(I,1) >*/
	    save1[i__] = y[i__ + y_dim1];
/*<             SAVE3(I)=YPRIME(I) >*/
	    save3[i__] = yprime[i__];
/*<             YI=Y(I,1) >*/
	    yi = y[i__ + y_dim1];
/*<             YP=YPRIME(I) >*/
	    yp = yprime[i__];
/*<             IF(ITOL.EQ.1) THEN >*/
	    if (*itol == 1) {
/*<                R=DMAX1(EPSJAC*DABS(YI),R0/(YMAX(I)*RTOL(1))) >*/
/* Computing MAX */
		d__1 = *epsjac * abs(yi), d__2 = r0 / (ymax[i__] * rtol[1]);
		r__ = max(d__1,d__2);
/*<             ELSE >*/
	    } else {
/*<                R=DMAX1(EPSJAC*DABS(YI),R0/YMAX(I)) >*/
/* Computing MAX */
		d__1 = *epsjac * abs(yi), d__2 = r0 / ymax[i__];
		r__ = max(d__1,d__2);
/*<             ENDIF >*/
	    }
/*<             Y(I,1)=Y(I,1)+R >*/
	    y[i__ + y_dim1] += r__;
/*<             YPRIME(I)=YPRIME(I)+R/CON >*/
	    yprime[i__] += r__ / *con;
/*<  161     CONTINUE >*/
/* L161: */
	}
/*<          CALL RESID(N,T,Y,WRKSPC,YPRIME,IPAR,RPAR,IERR) >*/
	(*resid)(n, t, &y[y_offset], &wrkspc[1], &yprime[1], &ipar[1], &rpar[
		1], ierr);
/*<          DO 261 JJ=J,N,MBND(3) >*/
	i__3 = *n;
	i__2 = mbnd[3];
	for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
/*<             Y(JJ,1)=SAVE1(JJ) >*/
	    y[jj + y_dim1] = save1[jj];
/*<             YPRIME(JJ)=SAVE3(JJ) >*/
	    yprime[jj] = save3[jj];
/*<             YJJ=Y(JJ,1) >*/
	    yjj = y[jj + y_dim1];
/*<             YJP=YPRIME(JJ) >*/
	    yjp = yprime[jj];
/*<             IF (ITOL.EQ.1) THEN >*/
	    if (*itol == 1) {
/*<                R=DMAX1(EPSJAC*DABS(YJJ),R0/(YMAX(JJ)*RTOL(1))) >*/
/* Computing MAX */
		d__1 = *epsjac * abs(yjj), d__2 = r0 / (ymax[jj] * rtol[1]);
		r__ = max(d__1,d__2);
/*<             ELSE >*/
	    } else {
/*<                R=DMAX1(EPSJAC*DABS(YJJ),R0/YMAX(JJ)) >*/
/* Computing MAX */
		d__1 = *epsjac * abs(yjj), d__2 = r0 / ymax[jj];
		r__ = max(d__1,d__2);
/*<             ENDIF >*/
	    }
/*<             D=CON/R >*/
	    d__ = *con / r__;
/*<             I1 = MAX0(JJ-MU,1) >*/
/* Computing MAX */
	    i__4 = jj - mu;
	    i1 = max(i__4,1);
/*<             I2 = MIN0(JJ+ML,N) >*/
/* Computing MIN */
	    i__4 = jj + ml;
	    i2 = min(i__4,*n);
/*<             II = JJ*(MBND(4)-1)-ML >*/
	    ii = jj * (mbnd[4] - 1) - ml;
/*<             DO 540 I = I1,I2 >*/
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/*<                TEMPRY = WRKSPC(I) - SAVE2(I) >*/
		tempry = wrkspc[i__] - save2[i__];
/*<                PWCOPY(II+I)=TEMPRY/R >*/
		pwcopy[ii + i__] = tempry / r__;
/*<                PW(II+I) = PWCOPY(II+I) >*/
		pw[ii + i__] = pwcopy[ii + i__];
/*<  540        CONTINUE >*/
/* L540: */
	    }
/*<  261     CONTINUE >*/
/* L261: */
	}
/*<  61   CONTINUE >*/
/* L61: */
    }
/*<       NRE=NRE+MBND(3) >*/
    *nre += mbnd[3];

/*<  70   IF (MITER.GT.2) THEN >*/
L70:
    if (*miter > 2) {
/*<          CALL DGBFA(PW,MBND(4),N,ML,MU,IPIV,IER) >*/
	dgbfa_(&pw[1], &mbnd[4], n, &ml, &mu, &ipiv[1], ier);
/*<          NDEC = NDEC + 1 >*/
	++(*ndec);
/*<       ELSE >*/
    } else {
/*<          CALL DEC(N,N,PW,IPIV,IER) >*/
	dec_(n, n, &pw[1], &ipiv[1], ier);
/*<          NDEC = NDEC + 1 >*/
	++(*ndec);
/*<       ENDIF >*/
    }
/*<       RETURN >*/
    return 0;
/* ---------------------- END OF SUBROUTINE PSET --------------------- */
/*<       END >*/
} /* pset_ */




/*<       SUBROUTINE DEC(N,NDIM,A,IP,IER) >*/
/* Subroutine */ int dec_(integer *n, integer *ndim, doublereal *a, integer *
	ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    /*static*/ integer i__, j, k, m;
    /*static*/ doublereal t;
    /*static*/ integer kp1, nm1;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/* ------------------------------------------------------------------- */
/*     MATRIX TRIANGULISATION BY GAUSSIAN ELIMINATION */
/*     INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAY A. */
/*     A = MATRIX TO BE TRIANGULARISED. */
/*     OUTPUT.. */
/*     A(I,J),  I.LE.J = UPPER TRIANGULAR FACTOR, U. */
/*     A(I,J),  I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I-L. */
/*     IP(K), K.LT.N = INDEX OF KTH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0. */
/*     IER = 0 IF MATRIX IS NON-SINGULAR, OR K IF FOUND TO BE SINGULAR */
/*                  AT STAGE K. */
/*     USE SOL TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*     DETERM(A) = IP(N)*A(1,1)*A(2,2)* . . . *A(N,N). */
/*     IF IP(N) = 0, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*     REFERENCE. */
/*     C.B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, C.A.C.M */
/*     15 (1972), P.274. */
/*     ------------------------------------------------------------------ */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER IER,N,NDIM >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  A(NDIM,N) >*/
/*<       INTEGER IP(N) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,J,K,KP1,M,NM1 >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*<       INTRINSIC DABS >*/
/*     .. */
/*     .. COMMON BLOCKS .. */
/*     .. */
/*<       IER = 0 >*/
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
/*<       IP(N) = 1 >*/
    ip[*n] = 1;
/*<       IF (N.EQ.1) GO TO 70 >*/
    if (*n == 1) {
	goto L70;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 60 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<          KP1 = K + 1 >*/
	kp1 = k + 1;
/*<          M = K >*/
	m = k;
/*<          DO 10 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<             IF (DABS(A(I,K)).GT.DABS(A(M,K))) M = I >*/
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/*<  10      CONTINUE >*/
/* L10: */
	}
/*<          IP(K) = M >*/
	ip[k] = m;
/*<          T = A(M,K) >*/
	t = a[m + k * a_dim1];
/*<          IF (M.EQ.K) GO TO 20 >*/
	if (m == k) {
	    goto L20;
	}
/*<          IP(N) = -IP(N) >*/
	ip[*n] = -ip[*n];
/*<          A(M,K) = A(K,K) >*/
	a[m + k * a_dim1] = a[k + k * a_dim1];
/*<          A(K,K) = T >*/
	a[k + k * a_dim1] = t;
/*<  20      IF (T.EQ.0.0D+0) GO TO 80 >*/
L20:
	if (t == 0.) {
	    goto L80;
	}
/*<          T = 1.0D+0/T >*/
	t = 1. / t;
/*<          DO 30 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<             A(I,K) = -A(I,K)*T >*/
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
/*<  30      CONTINUE >*/
/* L30: */
	}
/*<          DO 50 J = KP1,N >*/
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
/*<             T = A(M,J) >*/
	    t = a[m + j * a_dim1];
/*<             A(M,J) = A(K,J) >*/
	    a[m + j * a_dim1] = a[k + j * a_dim1];
/*<             A(K,J) = T >*/
	    a[k + j * a_dim1] = t;
/*<             IF (T.EQ.0.0D+0) GO TO 50 >*/
	    if (t == 0.) {
		goto L50;
	    }
/*<             DO 40 I = KP1,N >*/
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/*<                A(I,J) = A(I,J) + A(I,K)*T >*/
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
/*<  40         CONTINUE >*/
/* L40: */
	    }
/*<  50      CONTINUE >*/
L50:
	    ;
	}
/*<  60   CONTINUE >*/
/* L60: */
    }
/*<  70   K = N >*/
L70:
    k = *n;
/*<       IF (A(N,N).EQ.0.0D+0) GO TO 80 >*/
    if (a[*n + *n * a_dim1] == 0.) {
	goto L80;
    }
/*<       RETURN >*/
    return 0;
/*<  80   IER = K >*/
L80:
    *ier = k;
/*<       IP(N) = 0 >*/
    ip[*n] = 0;
/*<       RETURN >*/
    return 0;
/* --------------------- END OF SUBROUTINE DEC ---------------------- */
/*<       END >*/
} /* dec_ */




/*<       SUBROUTINE SOL(N,NDIM,A,B,IP) >*/
/* Subroutine */ int sol_(integer *n, integer *ndim, doublereal *a, 
	doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    /*static*/ integer i__, k, m;
    /*static*/ doublereal t;
    /*static*/ integer kb, km1, kp1, nm1;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER N,NDIM >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  A(NDIM,N),B(N) >*/
/*<       INTEGER IP(N) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,K,KB,KM1,KP1,M,NM1 >*/
/*     .. */
/*     .. COMMON BLOCKS .. */
/*     .. */
/*     ------------------------------------------------------------------ */
/*     SOLUTION OF LINEAR SYSTEM, A*X = B. */
/*     INPUT .. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF MATRIX A. */
/*     A = TRIANGULARISED MATRIX OBTAINED FROM DEC. */
/*     B = RIGHT HAND SIDE VECTOR. */
/*     IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*     DO NOT USE IF DEC HAS SET IER .NE. 0 */
/*     OUTPUT.. */
/*     B = SOLUTION VECTOR, X. */
/*     ------------------------------------------------------------------ */
/*<       IF (N.EQ.1) GO TO 50 >*/
    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    if (*n == 1) {
	goto L50;
    }
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;
/*<       DO 20 K = 1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<          KP1 = K + 1 >*/
	kp1 = k + 1;
/*<          M = IP(K) >*/
	m = ip[k];
/*<          T = B(M) >*/
	t = b[m];
/*<          B(M) = B(K) >*/
	b[m] = b[k];
/*<          B(K) = T >*/
	b[k] = t;
/*<          DO 10 I = KP1,N >*/
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/*<             B(I) = B(I) + A(I,K)*T >*/
	    b[i__] += a[i__ + k * a_dim1] * t;
/*<  10      CONTINUE >*/
/* L10: */
	}
/*<  20   CONTINUE >*/
/* L20: */
    }
/*<       DO 40 KB = 1,NM1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<          KM1 = N - KB >*/
	km1 = *n - kb;
/*<          K = KM1 + 1 >*/
	k = km1 + 1;
/*<          B(K) = B(K)/A(K,K) >*/
	b[k] /= a[k + k * a_dim1];
/*<          T = -B(K) >*/
	t = -b[k];
/*<          DO 30 I = 1,KM1 >*/
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             B(I) = B(I) + A(I,K)*T >*/
	    b[i__] += a[i__ + k * a_dim1] * t;
/*<  30      CONTINUE >*/
/* L30: */
	}
/*<  40   CONTINUE >*/
/* L40: */
    }
/*<  50   B(1) = B(1)/A(1,1) >*/
L50:
    b[1] /= a[a_dim1 + 1];
/*<       RETURN >*/
    return 0;
/* ------------------------- END OF SUBROUTINE SOL ------------------ */
/*<       END >*/
} /* sol_ */





//#pragma optimize("", off) //dgbfa

/*<       subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info) >*/
/* Subroutine */ int dgbfa_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    /*static*/ integer i__, j, k, l, m;
    /*static*/ doublereal t;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
	extern int daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    /*static*/ integer i0, j0, j1, lm, mm, ju;
    extern integer idamax_(integer *, doublereal *, integer *);
    /*static*/ integer jz, kp1, nm1;

/*<       integer lda,n,ml,mu,ipvt(1),info >*/
/*<       double precision abd(lda,1) >*/

/*     dgbfa factors a double precision band matrix by elimination. */

/*     dgbfa is usually called by dgbco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                contains the matrix in band storage.  the columns */
/*                of the matrix are stored in the columns of  abd  and */
/*                the diagonals of the matrix are stored in rows */
/*                ml+1 through 2*ml+mu+1 of  abd . */
/*                see the comments below for details. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */
/*                lda must be .ge. 2*ml + mu + 1 . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */
/*                0 .le. ml .lt. n . */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */
/*                0 .le. mu .lt. n . */
/*                more efficient if  ml .le. mu . */
/*     on return */

/*        abd     an upper triangular matrix in band storage and */
/*                the multipliers which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that dgbsl will divide by zero if */
/*                     called.  use  rcond  in dgbco for a reliable */
/*                     indication of singularity. */

/*     band storage */

/*           if  a  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ml = (band width below the diagonal) */
/*                   mu = (band width above the diagonal) */
/*                   m = ml + mu + 1 */
/*                   do 20 j = 1, n */
/*                      i1 = max0(1, j-mu) */
/*                      i2 = min0(n, j+ml) */
/*                      do 10 i = i1, i2 */
/*                         k = i - j + m */
/*                         abd(k,j) = a(i,j) */
/*                10    continue */
/*                20 continue */

/*           this uses rows  ml+1  through  2*ml+mu+1  of  abd . */
/*           in addition, the first  ml  rows in  abd  are used for */
/*           elements generated during the triangularization. */
/*           the total number of rows needed in  abd  is  2*ml+mu+1 . */
/*           the  ml+mu by ml+mu  upper left triangle and the */
/*           ml by ml  lower right triangle are not referenced. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,idamax */
/*     fortran max0,min0 */

/*     internal variables */

/*<       double precision t >*/
/*<       integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1 >*/


/*<       m = ml + mu + 1 >*/
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1 * 1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
/*<       info = 0 >*/
    *info = 0;

/*     zero initial fill-in columns */

/*<       j0 = mu + 2 >*/
    j0 = *mu + 2;
/*<       j1 = min0(n,m) - 1 >*/
    j1 = min(*n,m) - 1;
/*<       if (j1 .lt. j0) go to 30 >*/
    if (j1 < j0) {
	goto L30;
    }
/*<       do 20 jz = j0, j1 >*/
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
/*<          i0 = m + 1 - jz >*/
	i0 = m + 1 - jz;
/*<          do 10 i = i0, ml >*/
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
/*<             abd(i,jz) = 0.0d0 >*/
	    abd[i__ + jz * abd_dim1] = 0.;
/*<  10      continue >*/
/* L10: */
	}
/*<  20   continue >*/
/* L20: */
    }
/*<  30   continue >*/
L30:
/*<       jz = j1 >*/
    jz = j1;
/*<       ju = 0 >*/
    ju = 0;

/*     gaussian elimination with partial pivoting */

/*<       nm1 = n - 1 >*/
    nm1 = *n - 1;
/*<       if (nm1 .lt. 1) go to 130 >*/
    if (nm1 < 1) {
	goto L130;
    }
/*<       do 120 k = 1, nm1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<          kp1 = k + 1 >*/
	kp1 = k + 1;

/*        zero next fill-in column */

/*<          jz = jz + 1 >*/
	++jz;
/*<          if (jz .gt. n) go to 50 >*/
	if (jz > *n) {
	    goto L50;
	}
/*<          if (ml .lt. 1) go to 50 >*/
	if (*ml < 1) {
	    goto L50;
	}
/*<          do 40 i = 1, ml >*/
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             abd(i,jz) = 0.0d0 >*/
	    abd[i__ + jz * abd_dim1] = 0.;
/*<  40      continue >*/
/* L40: */
	}
/*<  50      continue >*/
L50:

/*        find l = pivot index */

/*<          lm = min0(ml,n-k) >*/
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
/*<          l = idamax(lm+1,abd(m,k),1) + m - 1 >*/
	i__2 = lm + 1;
	l = idamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
/*<          ipvt(k) = l + k - m >*/
	ipvt[k] = l + k - m;

/*        zero pivot implies this column already triangularized */

/*<          if (abd(l,k) .eq. 0.0d0) go to 100 >*/
	if (abd[l + k * abd_dim1] == 0.) {
	    goto L100;
	}

/*           interchange if necessary */

/*<          if (l .eq. m) go to 60 >*/
	if (l == m) {
	    goto L60;
	}
/*<          t = abd(l,k) >*/
	t = abd[l + k * abd_dim1];
/*<          abd(l,k) = abd(m,k) >*/
	abd[l + k * abd_dim1] = abd[m + k * abd_dim1];
/*<          abd(m,k) = t >*/
	abd[m + k * abd_dim1] = t;
/*<  60      continue >*/
L60:

/*           compute multipliers */

/*<          t = -1.0d0/abd(m,k) >*/
	t = -1. / abd[m + k * abd_dim1];
/*<          call dscal(lm,t,abd(m+1,k),1) >*/
	dscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           row elimination with column indexing */

/*<          ju = min0(max0(ju,mu+ipvt(k)),n) >*/
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
/*<          mm = m >*/
	mm = m;
/*<          if (ju .lt. kp1) go to 90 >*/
	if (ju < kp1) {
	    goto L90;
	}
/*<          do 80 j = kp1, ju >*/
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
/*<             l = l - 1 >*/
	    --l;
/*<             mm = mm - 1 >*/
	    --mm;
/*<             t = abd(l,j) >*/
	    t = abd[l + j * abd_dim1];
/*<             if (l .eq. mm) go to 70 >*/
	    if (l == mm) {
		goto L70;
	    }
/*<             abd(l,j) = abd(mm,j) >*/
	    abd[l + j * abd_dim1] = abd[mm + j * abd_dim1];
/*<             abd(mm,j) = t >*/
	    abd[mm + j * abd_dim1] = t;
/*<  70         continue >*/
L70:
/*<             call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1) >*/
	    daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/*<  80      continue >*/
/* L80: */
	}
/*<  90      continue >*/
L90:
/*<          go to 110 >*/
	goto L110;
/*<  100     continue >*/
L100:
/*<          info = k >*/
	*info = k;
/*<  110     continue >*/
L110:
/*<  120  continue >*/
/* L120: */
	;
    }
/*<  130  continue >*/
L130:
/*<       ipvt(n) = n >*/
    ipvt[*n] = *n;
/*<       if (abd(m,n) .eq. 0.0d0) info = n >*/
    if (abd[m + *n * abd_dim1] == 0.) {
		*info = *n;
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* dgbfa_ */

//#pragma optimize("", on) //2dgbfa

/* -------------------------------------------------------------------------- */
/*<       subroutine daxpy(n,da,dx,incx,dy,incy) >*/
/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    /*static*/ integer i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */

/*<       double precision dx(1),dy(1),da >*/
/*<       integer i,incx,incy,ix,iy,m,mp1,n >*/

/*<       if(n.le.0)return >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       if (da .eq. 0.0d0) return >*/
    if (*da == 0.) {
	return 0;
    }
/*<       if(incx.eq.1.and.incy.eq.1)go to 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

/*<       ix = 1 >*/
    ix = 1;
/*<       iy = 1 >*/
    iy = 1;
/*<       if(incx.lt.0)ix = (-n+1)*incx + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       if(incy.lt.0)iy = (-n+1)*incy + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       do 10 i = 1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          dy(iy) = dy(iy) + da*dx(ix) >*/
	dy[iy] += *da * dx[ix];
/*<          ix = ix + incx >*/
	ix += *incx;
/*<          iy = iy + incy >*/
	iy += *incy;
/*<  10   continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

/*<  20   m = mod(n,4) >*/
L20:
    m = *n % 4;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          dy(i) = dy(i) + da*dx(i) >*/
	dy[i__] += *da * dx[i__];
/*<  30   continue >*/
/* L30: */
    }
/*<       if( n .lt. 4 ) return >*/
    if (*n < 4) {
	return 0;
    }
/*<  40   mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,4 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
/*<          dy(i) = dy(i) + da*dx(i) >*/
	dy[i__] += *da * dx[i__];
/*<          dy(i + 1) = dy(i + 1) + da*dx(i + 1) >*/
	dy[i__ + 1] += *da * dx[i__ + 1];
/*<          dy(i + 2) = dy(i + 2) + da*dx(i + 2) >*/
	dy[i__ + 2] += *da * dx[i__ + 2];
/*<          dy(i + 3) = dy(i + 3) + da*dx(i + 3) >*/
	dy[i__ + 3] += *da * dx[i__ + 3];
/*<  50   continue >*/
/* L50: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* daxpy_ */



/* --------------------------------------------------------------------------- */
/*<       subroutine  dscal(n,da,dx,incx) >*/
/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    /*static*/ integer i__, m, ix, mp1;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified to correct problem with negative increment, 8/21/90. */

/*<       double precision da,dx(1) >*/
/*<       integer i,incx,ix,m,mp1,n >*/

/*<       if(n.le.0)return >*/
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       if(incx.eq.1)go to 20 >*/
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

/*<       ix = 1 >*/
    ix = 1;
/*<       if(incx.lt.0)ix = (-n+1)*incx + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       do 10 i = 1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          dx(ix) = da*dx(ix) >*/
	dx[ix] = *da * dx[ix];
/*<          ix = ix + incx >*/
	ix += *incx;
/*<  10   continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

/*<  20   m = mod(n,5) >*/
L20:
    m = *n % 5;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          dx(i) = da*dx(i) >*/
	dx[i__] = *da * dx[i__];
/*<  30   continue >*/
/* L30: */
    }
/*<       if( n .lt. 5 ) return >*/
    if (*n < 5) {
	return 0;
    }
/*<  40   mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,5 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
/*<          dx(i) = da*dx(i) >*/
	dx[i__] = *da * dx[i__];
/*<          dx(i + 1) = da*dx(i + 1) >*/
	dx[i__ + 1] = *da * dx[i__ + 1];
/*<         dx(i + 2) = da*dx(i + 2) >*/
	dx[i__ + 2] = *da * dx[i__ + 2];
/*<         dx(i + 3) = da*dx(i + 3) >*/
	dx[i__ + 3] = *da * dx[i__ + 3];
/*<         dx(i + 4) = da*dx(i + 4) >*/
	dx[i__ + 4] = *da * dx[i__ + 4];
/*<  50   continue >*/
/* L50: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* dscal_ */



/* -------------------------------------------------------------------------- */
/*<       integer function idamax(n,dx,incx) >*/
integer idamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    /*static*/ doublereal dmax__;
    /*static*/ integer i__, ix;


/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified to correct problem with negative increment, 8/21/90. */

/*<       double precision dx(1),dmax >*/
/*<       integer i,incx,ix,n >*/

/*<       idamax = 0 >*/
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
/*<       if( n .lt. 1 ) return >*/
    if (*n < 1) {
	return ret_val;
    }
/*<       idamax = 1 >*/
    ret_val = 1;
/*<       if(n.eq.1)return >*/
    if (*n == 1) {
	return ret_val;
    }
/*<       if(incx.eq.1)go to 20 >*/
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

/*<       ix = 1 >*/
    ix = 1;
/*<       if(incx.lt.0)ix = (-n+1)*incx + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       dmax = dabs(dx(ix)) >*/
    dmax__ = (d__1 = dx[ix], abs(d__1));
/*<       ix = ix + incx >*/
    ix += *incx;
/*<       do 10 i = 2,n >*/
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          if(dabs(dx(ix)).le.dmax) go to 5 >*/
	if ((d__1 = dx[ix], abs(d__1)) <= dmax__) {
	    goto L5;
	}
/*<          idamax = i >*/
	ret_val = i__;
/*<          dmax = dabs(dx(ix)) >*/
	dmax__ = (d__1 = dx[ix], abs(d__1));
/*<     5    ix = ix + incx >*/
L5:
	ix += *incx;
/*<  10   continue >*/
/* L10: */
    }
/*<       return >*/
    return ret_val;

/*        code for increment equal to 1 */

/*<  20   dmax = dabs(dx(1)) >*/
L20:
    dmax__ = abs(dx[1]);
/*<       do 30 i = 2,n >*/
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          if(dabs(dx(i)).le.dmax) go to 30 >*/
	if ((d__1 = dx[i__], abs(d__1)) <= dmax__) {
	    goto L30;
	}
/*<          idamax = i >*/
	ret_val = i__;
/*<          dmax = dabs(dx(i)) >*/
	dmax__ = (d__1 = dx[i__], abs(d__1));
/*<  30   continue >*/
L30:
	;
    }
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* idamax_ */

/* -------------------------------------------------------------------------- */
/*<       subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job) >*/
/* Subroutine */ int dgbsl_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *b, integer *job)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    /*static*/ integer k, l, m;
    /*static*/ doublereal t;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    /*static*/ integer kb, la, lb, lm, nm1;

/*<       integer lda,n,ml,mu,ipvt(*),job >*/
/*<       double precision abd(lda,*),b(*) >*/

/*     dgbsl solves the double precision band system */
/*     a * x = b  or  trans(a) * x = b */
/*     using the factors computed by dgbco or dgbfa. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                the output from dgbco or dgbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */

/*        ipvt    integer(n) */
/*                the pivot vector from dgbco or dgbfa. */

/*        b       double precision(n) */
/*                the right hand side vector. */

/*        job     integer */
/*                = 0         to solve  a*x = b , */
/*                = nonzero   to solve  trans(a)*x = b , where */
/*                            trans(a)  is the transpose. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of lda .  it will not occur if the subroutines are */
/*        called correctly and if dgbco has set rcond .gt. 0.0 */
/*        or dgbfa has set info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z) */
/*           if (rcond is too small) go to ... */
/*           do 10 j = 1, p */
/*              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */
/*     fortran min0 */

/*     internal variables */

/*<       double precision ddot,t >*/
/*<       integer k,kb,l,la,lb,lm,m,nm1 >*/

/*<       m = mu + ml + 1 >*/
    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1 * 1;
    abd -= abd_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
/*<       nm1 = n - 1 >*/
    nm1 = *n - 1;
/*<       if (job .ne. 0) go to 50 >*/
    if (*job != 0) {
	goto L50;
    }

/*        job = 0 , solve  a * x = b */
/*        first solve l*y = b */

/*<       if (ml .eq. 0) go to 30 >*/
    if (*ml == 0) {
	goto L30;
    }
/*<       if (nm1 .lt. 1) go to 30 >*/
    if (nm1 < 1) {
	goto L30;
    }
/*<       do 20 k = 1, nm1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/*<          lm = min0(ml,n-k) >*/
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
/*<          l = ipvt(k) >*/
	l = ipvt[k];
/*<          t = b(l) >*/
	t = b[l];
/*<          if (l .eq. k) go to 10 >*/
	if (l == k) {
	    goto L10;
	}
/*<          b(l) = b(k) >*/
	b[l] = b[k];
/*<          b(k) = t >*/
	b[k] = t;
/*<  10      continue >*/
L10:
/*<          call daxpy(lm,t,abd(m+1,k),1,b(k+1),1) >*/
	daxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1);
/*<  20   continue >*/
/* L20: */
    }
/*<  30   continue >*/
L30:

/*        now solve  u*x = y */

/*<       do 40 kb = 1, n >*/
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
/*<          k = n + 1 - kb >*/
	k = *n + 1 - kb;
/*<          b(k) = b(k)/abd(m,k) >*/
	b[k] /= abd[m + k * abd_dim1];
/*<          lm = min0(k,m) - 1 >*/
	lm = min(k,m) - 1;
/*<          la = m - lm >*/
	la = m - lm;
/*<          lb = k - lm >*/
	lb = k - lm;
/*<          t = -b(k) >*/
	t = -b[k];
/*<          call daxpy(lm,t,abd(la,k),1,b(lb),1) >*/
	daxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/*<  40   continue >*/
/* L40: */
    }
/*<       go to 100 >*/
    goto L100;
/*<  50   continue >*/
L50:

/*        job = nonzero, solve  trans(a) * x = b */
/*        first solve  trans(u)*y = b */

/*<       do 60 k = 1, n >*/
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*<          lm = min0(k,m) - 1 >*/
	lm = min(k,m) - 1;
/*<          la = m - lm >*/
	la = m - lm;
/*<          lb = k - lm >*/
	lb = k - lm;
/*<          t = ddot(lm,abd(la,k),1,b(lb),1) >*/
	t = ddot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/*<          b(k) = (b(k) - t)/abd(m,k) >*/
	b[k] = (b[k] - t) / abd[m + k * abd_dim1];
/*<  60   continue >*/
/* L60: */
    }

/*        now solve trans(l)*x = y */

/*<       if (ml .eq. 0) go to 90 >*/
    if (*ml == 0) {
	goto L90;
    }
/*<       if (nm1 .lt. 1) go to 90 >*/
    if (nm1 < 1) {
	goto L90;
    }
/*<       do 80 kb = 1, nm1 >*/
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
/*<          k = n - kb >*/
	k = *n - kb;
/*<          lm = min0(ml,n-k) >*/
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
/*<          b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1) >*/
	b[k] += ddot_(&lm, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &
		c__1);
/*<          l = ipvt(k) >*/
	l = ipvt[k];
/*<          if (l .eq. k) go to 70 >*/
	if (l == k) {
	    goto L70;
	}
/*<          t = b(l) >*/
	t = b[l];
/*<          b(l) = b(k) >*/
	b[l] = b[k];
/*<          b(k) = t >*/
	b[k] = t;
/*<  70      continue >*/
L70:
/*<  80   continue >*/
/* L80: */
	;
    }
/*<  90   continue >*/
L90:
/*<  100  continue >*/
L100:
/*<       return >*/
    return 0;
/*<       end >*/
} /* dgbsl_ */

/* --------------------------------------------------------------------------- */
/*<       double precision function ddot(n,dx,incx,dy,incy) >*/
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    /*static*/ integer i__, m;
    /*static*/ doublereal dtemp;
    /*static*/ integer ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */

/*<       double precision dx(1),dy(1),dtemp >*/
/*<       integer i,incx,incy,ix,iy,m,mp1,n >*/

/*<       ddot = 0.0d0 >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
/*<       dtemp = 0.0d0 >*/
    dtemp = 0.;
/*<       if(n.le.0)return >*/
    if (*n <= 0) {
	return ret_val;
    }
/*<       if(incx.eq.1.and.incy.eq.1)go to 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

/*<       ix = 1 >*/
    ix = 1;
/*<       iy = 1 >*/
    iy = 1;
/*<       if(incx.lt.0)ix = (-n+1)*incx + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       if(incy.lt.0)iy = (-n+1)*incy + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       do 10 i = 1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          dtemp = dtemp + dx(ix)*dy(iy) >*/
	dtemp += dx[ix] * dy[iy];
/*<          ix = ix + incx >*/
	ix += *incx;
/*<          iy = iy + incy >*/
	iy += *incy;
/*<  10   continue >*/
/* L10: */
    }
/*<       ddot = dtemp >*/
    ret_val = dtemp;
/*<       return >*/
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

/*<  20   m = mod(n,5) >*/
L20:
    m = *n % 5;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          dtemp = dtemp + dx(i)*dy(i) >*/
	dtemp += dx[i__] * dy[i__];
/*<  30   continue >*/
/* L30: */
    }
/*<       if( n .lt. 5 ) go to 60 >*/
    if (*n < 5) {
	goto L60;
    }
/*<  40   mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,5 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
/*<    >*/
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/*<  50   continue >*/
/* L50: */
    }
/*<  60   ddot = dtemp >*/
L60:
    ret_val = dtemp;
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* ddot_ */

/* --------------------------------------------------------------------------- */
/*<       SUBROUTINE ERRORS(N,TQ,EDN,E,EUP,BND,EDDN) >*/
/* Subroutine */ int errors_(integer *n, doublereal *tq, doublereal *edn, 
	doublereal *e, doublereal *eup, doublereal *bnd, doublereal *eddn)
{
    /*static*/ doublereal sqhol;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     *************************************************** */

/*     THIS ROUTINE CALCULATES ERRORS USED IN TESTS */
/*     IN STIFF . */

/*     *************************************************** */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER N >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  TQ(5) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*     .. */
/*<       SQHOL = DBLE(FLOAT(N)) >*/
    /* Parameter adjustments */
    --tq;

    /* Function Body */
    sqhol = (doublereal) (*n);
/*<       EDN = TQ(1)*TQ(1)*SQHOL >*/
    *edn = tq[1] * tq[1] * sqhol;

/*     ** ERROR ASSOCIATED WITH  METHOD OF ORDER ONE LOWER. */

/*<       E = TQ(2)*TQ(2)*SQHOL >*/
    *e = tq[2] * tq[2] * sqhol;

/*     ** ERROR ASSOCIATED WITH PRESENT ORDER */

/*<       EUP = TQ(3)*TQ(3)*SQHOL >*/
    *eup = tq[3] * tq[3] * sqhol;

/*     ** ERROR ASSOCIATED WITH HIGHER ORDER METHOD */

/*<       BND = TQ(4)*TQ(4)*SQHOL*0.5D+0 >*/
    *bnd = tq[4] * tq[4] * sqhol * .5;
/*<       EDDN=TQ(5)*TQ(5)*SQHOL >*/
    *eddn = tq[5] * tq[5] * sqhol;

/*     ** ERROR ASSOCIATED WITH METHOD OF ORDER TWO LOWER. */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* errors_ */

/* -------------------------------------------------------------------------- */
/*<       SUBROUTINE PRDICT(T,H,Y,L,N,IPAR,RPAR,IERR) >*/
/* Subroutine */ int prdict_(doublereal *t, doublereal *h__, doublereal *y, 
	integer *l, integer *n, integer *ipar, doublereal *rpar, integer *
	ierr)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    /*static*/ integer i__, j2;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/* ********************************************************************** */
/*     PREDICTS A VALUE FOR Y AT (T+H) GIVEN THE HISTORY ARRAY AT T. */
/* ********************************************************************** */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER L,N >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  Y(N,12),IPAR(*),RPAR(*) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,J2 >*/
/*     .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*     .. */
/*     .. COMMON BLOCKS .. */
/*     .. */
/*<       DO 10 I=1,N >*/
    /* Parameter adjustments */
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    --ipar;
    --rpar;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          DO 20 J2 = 2,L          >*/
	i__2 = *l;
	for (j2 = 2; j2 <= i__2; ++j2) {
/*<             Y(I,1) = Y(I,1) + Y(I,J2)           >*/
	    y[i__ + y_dim1] += y[i__ + j2 * y_dim1];
/*<  20      CONTINUE >*/
/* L20: */
	}
/*<  10   CONTINUE         >*/
/* L10: */
    }
/*<       T = T + H >*/
    *t += *h__;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* prdict_ */

/* ------------------------------------------------------------------------ */
/*<    >*/
/* Subroutine */ int itrat2_(doublereal *qqq, doublereal *y, doublereal *
	yprime, integer *n, doublereal *t, doublereal *hbeta, doublereal *
	errbnd, doublereal *arh, doublereal *crate, doublereal *tcrate, 
	integer *m, logical *worked, doublereal *ymax, doublereal *error, 
	doublereal *save1, doublereal *save2, doublereal *scale, doublereal *
	pw, integer *mf, integer *mbnd, integer *nind1, integer *nind2, 
	integer *nind3, integer *ipiv, integer *lmb, integer *itol, 
	doublereal *rtol, doublereal *atol, integer *ipar, doublereal *rpar, 
	doublereal *hused, integer *nbsol, integer *nre, integer *nqused, 
	S_fp resid, integer *ierr)
{
    /* Initialized data */

	/* its ok to leave this on zero */
    static doublereal zero = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    /*static*/ doublereal d__;
    /*static*/ integer i__;
    extern /* Subroutine */ int dgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    /*static*/ doublereal ayi;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER M,N,ITOL >*/
/*<       LOGICAL WORKED >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<    >*/
/*<       INTEGER IPIV(N),MBND(4) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I >*/
/*     .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*<       EXTERNAL SOL,DGBSL,RESID >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*<       INTRINSIC DMAX1,DMIN1 >*/
/*     .. */
/*     .. COMMON BLOCKS .. */
/*<       INTEGER NBSOL,NRE >*/
/*     .. */
/*     .. DATA STATEMENTS .. */
/*<       DATA  ZERO/0.0D+0/ >*/
    /* Parameter adjustments */
    --ipiv;
    --scale;
    --save2;
    --save1;
    --error;
    --ymax;
    --arh;
    --yprime;
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    --pw;
    --mbnd;
    --rtol;
    --atol;
    --ipar;
    --rpar;

    /* Function Body */
/*     .. */
/*<       DO 5 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          AYI = DABS(Y(I,1)) >*/
	ayi = (d__1 = y[i__ + y_dim1], abs(d__1));
/*<          IF(ITOL.EQ.1) THEN >*/
	if (*itol == 1) {
/*<             SCALE(I) = YMAX(I) >*/
	    scale[i__] = ymax[i__];
/*<          ELSE IF(ITOL.EQ.2) THEN >*/
	} else if (*itol == 2) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[1] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.3) THEN >*/
	} else if (*itol == 3) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[1] * ayi + atol[i__];
/*<          ELSE IF(ITOL.EQ.4) THEN >*/
	} else if (*itol == 4) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[i__] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.5) THEN >*/
	} else if (*itol == 5) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[i__] * ayi + atol[i__];
/*<          ENDIF >*/
	}
/*<  5    CONTINUE >*/
/* L5: */
    }
/*<       IF(NIND2.NE.0) THEN >*/
    if (*nind2 != 0) {
/*<          DO 11 I = NIND1+1,NIND2+NIND1 >*/
	i__1 = *nind2 + *nind1;
	for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/HUSED >*/
	    scale[i__] /= *hused;
/*<  11      CONTINUE >*/
/* L11: */
	}
/*<       ENDIF >*/
    }
/*<       IF(NIND3.NE.0) THEN >*/
    if (*nind3 != 0) {
/*<          DO 12 I = NIND1 +NIND2 + 1,NIND3+NIND2+NIND1 >*/
	i__1 = *nind3 + *nind2 + *nind1;
	for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/(HUSED**2) >*/
/* Computing 2nd power */
	    d__1 = *hused;
	    scale[i__] /= d__1 * d__1;
/*<  12      CONTINUE >*/
/* L12: */
	}
/*<       ENDIF >*/
    }

/*<       IF(LMB.EQ.1) GOTO 25 >*/
    if (*lmb == 1) {
	goto L25;
    }

/*<       call resid(n,t,y,save2,yprime,ipar,rpar,ierr) >*/
    (*resid)(n, t, &y[y_offset], &save2[1], &yprime[1], &ipar[1], &rpar[1], 
	    ierr);
/*<       IF(MF.GE.23) THEN >*/
    if (*mf >= 23) {
/*<          CALL DGBSL(PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE2,0) >*/
	dgbsl_(&pw[1], &mbnd[4], n, &mbnd[1], &mbnd[2], &ipiv[1], &save2[1], &
		c__0);
/*<          NBSOL = NBSOL + 1 >*/
	++(*nbsol);
/*<       ELSE >*/
    } else {
/*<          CALL SOL(N,N,PW,SAVE2,IPIV) >*/
	sol_(n, n, &pw[1], &save2[1], &ipiv[1]);
/*<          NBSOL = NBSOL + 1 >*/
	++(*nbsol);
/*<       ENDIF >*/
    }
/*<       D = ZERO >*/
    d__ = zero;
/*<       DO 20 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ERROR(I) = ERROR(I) - SAVE2(I) >*/
	error[i__] -= save2[i__];
/*<          D = D + (SAVE2(I)/(SCALE(I)))**2 >*/
/* Computing 2nd power */
	d__1 = save2[i__] / scale[i__];
	d__ += d__1 * d__1;
/*<          SAVE1(I) = Y(I,1) + ERROR(I) >*/
	save1[i__] = y[i__ + y_dim1] + error[i__];
/*<  20   CONTINUE >*/
/* L20: */
    }
/*<       IF(ITOL.EQ.1) D=D/(RTOL(1)**2) >*/
    if (*itol == 1) {
/* Computing 2nd power */
	d__1 = rtol[1];
	d__ /= d__1 * d__1;
    }
/*<       TCRATE = TCRATE + CRATE >*/
    *tcrate += *crate;
/*<       D1 = D >*/
    g_d1 = d__;
/*<       M = 1 >*/
    *m = 1;
/*<       do 1014 i=1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          yprime(i)=(save1(i)-arh(i))/qqq >*/
	yprime[i__] = (save1[i__] - arh[i__]) / *qqq;
/*<  1014 continue >*/
/* L1014: */
    }
/*<       NRE = NRE + 1 >*/
    ++(*nre);
/*<  25   CONTINUE >*/
L25:
/*<       WORKED = .TRUE. >*/
    *worked = TRUE_;
/*<  30   CONTINUE >*/
L30:
/*<       call resid(n,t,save1,save2,yprime,ipar,rpar,ierr) >*/
    (*resid)(n, t, &save1[1], &save2[1], &yprime[1], &ipar[1], &rpar[1], ierr)
	    ;
/*<       nre=nre+1 >*/
    ++(*nre);

/*     IF WE ARE HERE THEN PARTIALS ARE O.K. */

/*<       IF( MF.GE. 23) THEN >*/
    if (*mf >= 23) {
/*<          CALL DGBSL(PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE2,0) >*/
	dgbsl_(&pw[1], &mbnd[4], n, &mbnd[1], &mbnd[2], &ipiv[1], &save2[1], &
		c__0);
/*<          NBSOL=NBSOL + 1 >*/
	++(*nbsol);
/*<       ELSE >*/
    } else {
/*<          CALL SOL(N,N,PW,SAVE2,IPIV) >*/
	sol_(n, n, &pw[1], &save2[1], &ipiv[1]);
/*<          NBSOL = NBSOL + 1 >*/
	++(*nbsol);
/*<       ENDIF >*/
    }

/*     WE NOW CALCULATE A WEIGHTED RMS TYPE NORM */

/*<       D = ZERO >*/
    d__ = zero;
/*<       DO 50 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ERROR(I) = ERROR(I) - SAVE2(I) >*/
	error[i__] -= save2[i__];
/*<          D = D + (SAVE2(I)/(SCALE(I)))**2 >*/
/* Computing 2nd power */
	d__1 = save2[i__] / scale[i__];
	d__ += d__1 * d__1;
/*<          SAVE1(I) = Y(I,1) + ERROR(I) >*/
	save1[i__] = y[i__ + y_dim1] + error[i__];
/*<  50   CONTINUE >*/
/* L50: */
    }
/*<       IF(ITOL.EQ.1) D=D/(RTOL(1)**2) >*/
    if (*itol == 1) {
/* Computing 2nd power */
	d__1 = rtol[1];
	d__ /= d__1 * d__1;
    }
/* ------------------------------------------------------------------- */
/*     TEST FOR CONVERGENCE.  IF M.GT.0 , AN ESTIMATE OF THE CONVERGENCE */
/*     RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST. */
/* ------------------------------------------------------------------- */
/*<       IF (M.NE.0) THEN >*/
    if (*m != 0) {
/*<          IF (D1.NE.ZERO) CRATE = DMAX1(0.9D+0*CRATE,D/D1) >*/
	if (g_d1 != zero) {
/* Computing MAX */
	    d__1 = *crate * .9, d__2 = d__ / g_d1;
	    *crate = max(d__1,d__2);
	}
/*<       END IF >*/
    }
/*<       TCRATE = TCRATE + CRATE >*/
    *tcrate += *crate;
/*<    >*/
/* Computing MIN */
    d__1 = 1., d__2 = *crate * 2.;
    if (d__ * min(d__1,d__2) < *errbnd / (doublereal) (*nqused)) {
	return 0;
    }
/*<       IF (M.NE.0) THEN >*/
    if (*m != 0) {
/*<          IF (D.GT.D1) THEN >*/
	if (d__ > g_d1) {
/*<             WORKED = .FALSE. >*/
	    *worked = FALSE_;
/*<             RETURN >*/
	    return 0;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       D1 = D >*/
    g_d1 = d__;
/*<       IF (M.EQ.4) THEN >*/
    if (*m == 4) {
/*<          WORKED = .FALSE. >*/
	*worked = FALSE_;
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       M = M + 1 >*/
    ++(*m);
/*<       do 40 i=1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          yprime(i)=(save1(i)-arh(i))/qqq >*/
	yprime[i__] = (save1[i__] - arh[i__]) / *qqq;
/*<  40   continue >*/
/* L40: */
    }
/*<       GO TO 30 >*/
    goto L30;
/*<       END >*/
} /* itrat2_ */

/* -------------------------------------------------------------------------- */
/*<    >*/
/* Subroutine */ int stiff_(doublereal *h__, doublereal *hmax, doublereal *
	hmin, integer *jstart, integer *kflag, integer *mf, integer *mbnd, 
	integer *nind1, integer *nind2, integer *nind3, doublereal *t, 
	doublereal *tout, doublereal *tend, doublereal *y, doublereal *yprime,
	 integer *n, doublereal *ymax, doublereal *error, doublereal *save1, 
	doublereal *save2, doublereal *scale, doublereal *pw, doublereal *
	pwcopy, doublereal *yhold, doublereal *ynhold, doublereal *arh, 
	integer *ipiv, integer *lout, integer *maxder, integer *itol, 
	doublereal *rtol, doublereal *atol, doublereal *rpar, integer *ipar, 
	S_fp pderv, S_fp resid, integer *nqused, integer *nstep, integer *
	nfail, integer *nre, integer *nje, integer *ndec, integer *nbsol, 
	integer *npset, integer *ncoset, integer *maxord, integer *maxstp, 
	doublereal *uround, doublereal *epsjac, doublereal *hused, integer *
	ierr)
{
    /* Initialized data */

    static doublereal el[10] = { 0.0,1. };
    static doublereal elst[10] = { 0.0,1. };
    static doublereal oldlo = 1.;
    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* Format strings */
    static char fmt_1975[] = "(/,/,\002IERR IS NON-ZERO BECAUSE OF AN ILLEGA\
L FUNCTION   CALL\002)";
    static char fmt_9161[] = "(/,/,\002STEPSIZE IS TOO SMALL\002)";

    /* System generated locals */
    integer y_dim1, y_offset, yhold_dim1, yhold_offset, ynhold_dim1, 
	    ynhold_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), e_wsfe(), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle();

    /* Local variables */
    
	/* grn: checked on being really local */
	/*static*/ doublereal demb;
    /*static*/ integer idid;
    /*static*/ integer meth;
    /*static*/ doublereal told;
    /*static*/ integer newq;
    
	extern /* Subroutine */ int pset_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, S_fp, S_fp, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *);

    /*static*/ doublereal vtol;
    /*static*/ doublereal d__;
    /*static*/ integer i__, j;
    /*static*/ doublereal efail;

    extern /* Subroutine */ int dgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *);

    /*static*/ integer iredo;
    /*static*/ integer iiter;
    
	extern /* Subroutine */ int coset_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

    /*static*/ doublereal delst, ddown;

    /*static*/ integer j1, j2, m1, m2;
    /*static*/ doublereal t0, y0;
	
	/* avold2 is only assigned a value, not used */
	/*static*/ doublereal avold2;

    /*static*/ integer iiter2;

	/* avnew2 is only assigned a value, not used */
    /*static*/ doublereal avnew2;

    extern /* Subroutine */ int itrat2_(doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     integer *, integer *, S_fp, integer *);

    /*static*/ integer m3step;
    /*static*/ doublereal ho;
    /*static*/ integer ll;
    /*static*/ doublereal qi;
    /*static*/ doublereal qq, frfail, tq[5], plfail;
    
	extern /* Subroutine */ int hchose_(doublereal *, doublereal *, logical *)
	    , rscale_(integer *, integer *, doublereal *, doublereal *);

    /*static*/ doublereal prfail, avoldj;
    /*static*/ doublereal trange, dddown;
    
	extern /* Subroutine */ int prdict_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *);
    
	/* ovride is given as argument to hchos but this is ok */
	/* in that function it is assigned a value */
	/*static*/ logical ovride;

	/* given as argument to itrat2 */
    /*static*/ logical worked;

    /*static*/ integer newpar, nrenew;

    extern /* Subroutine */ int interp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), cpyary_(
	    integer *, doublereal *, doublereal *);

    /*static*/ integer jm1, jp1;
    
	extern /* Subroutine */ int errors_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    
	/*static*/ doublereal twodwn, pr1, pr2, pr3, pr0;
	/*static*/ doublereal fac;

    /*static*/ integer j2m1;
    /*static*/ doublereal red;
    /*static*/ integer ier;
    /*static*/ doublereal ayi;
    
	extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    
	/* gokje */
	/*static*/ doublereal enq0, enq1, enq2, enq3;

	/*static*/ integer lmp2, lmp3, nqp2;

    /* Fortran I/O blocks */
    static cilist io___245 = { 0, 6, 0, fmt_1975, 0 };
    static cilist io___246 = { 0, 6, 0, fmt_9161, 0 };
    static cilist io___248 = { 0, 0, 0, 0, 0 };


/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     ------------------------------------------------------------------ */
/*     THE SUBROUTINE STIFF PERFORMS ONE STEP OF THE INTEGRATION OF AN */
/*     INITIAL VALUE PROBLEM FOR A SYSTEM OF */
/*     IMPLICIT DIFFERENTIAL ALGEBRAIC EQUATIONS. */
/*     COMMUNICATION WITH STIFF IS DONE WITH THE FOLLOWING VARIABLES.. */
/*     Y      AN N BY LMAX+3 ARRAY CONTAINING THE DEPENDENT VARIABLES */
/*              AND THEIR BACKWARD DIFFERENCES.  MAXDER (=LMAX-1) IS THE */
/*              MAXIMUM ORDER AVAILABLE.  SEE SUBROUTINE COSET. */
/*              Y(I,J+1) CONTAINS THE JTH BACKWARD DIFFERENCE OF Y(I) */
/*     T      THE INDEPENDENT VARIABLE. T IS UPDATED ON EACH STEP TAKEN. */
/*     H      THE STEPSIZE TO BE ATTEMPTED ON THE NEXT STEP. */
/*              H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING */
/*              THE PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE BUT */
/*              ITS SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM. */
/*     HMIN   THE MINIMUM AND MAXIMUM ABSOLUTE VALUE OF THE STEPSIZE */
/*     HMAX   TO BE USED FOR THE STEP.  THESE MAY BE CHANGED AT ANY */
/*              TIME BUT WILL NOT TAKE EFFECT UNTIL THE NEXT H CHANGE. */
/*     RTOL,ATOL  THE ERROR BOUNDS. SEE DESCRIPTION IN OVDRIV. */
/*     N      THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS. */
/*     MF     THE METHOD FLAG.  MUST BE SET TO 21,22,23 OR 24 AT PRESENT */
/*     KFLAG  A COMPLETION FLAG WITH THE FOLLOWING MEANINGS.. */
/*                  0  THE STEP WAS SUCCESSFUL */
/*                 -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED */
/*                       WITH ABS(H) = HMIN. */
/*                 -2  THE REQUESTED ERROR IS SMALLER THAN CAN */
/*                       BE HANDLED FOR THIS PROBLEM. */
/*                 -3  CORRECTOR CONVERGENCE COULD NOT BE */
/*                       ACHIEVED FOR ABS(H)=HMIN. */
/*            ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF T AND */
/*            THE Y ARRAY ARE AS AT THE BEGINNING OF THE LAST */
/*            STEP ATTEMPTED, AND H IS THE LAST STEP SIZE ATTEMPTED. */
/*     JSTART  AN INTEGER USED ON INPUT AND OUTPUT. */
/*          ON INPUT IT HAS THE FOLLOWING VALUES AND MEANINGS.. */
/*              0  PERFORM THE FIRST STEP. */
/*          .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST */
/*          .LT.0  TAKE THE NEXT STEP WITH A NEW VALUE OF H OR N. */
/*          ON EXIT, JSTART IS NQUSED, THE ORDER OF THE METHOD LAST USED. */
/*     YMAX     AN ARRAY OF N ELEMENTS WITH WHICH THE ESTIMATED LOCAL */
/*              ERRORS IN Y ARE COMPARED */
/*     ERROR    AN ARRAY OF N ELEMENTS. */
/*     SAVE1,2  TWO ARRAYS OF WORKING SPACE BOTH OF LENGTH N. */
/*     PW       A BLOCK OF LOCATIONS USED FOR PARTIAL DERIVATIVES */
/*     IPIV     AN INTEGER ARRAY OF LENGTH N USED FOR PIVOT INFORMATION. */
/*     JNEWIM   IS TO INDICATE IF PRESENT ITERATION MATRIX */
/*                WAS FORMED USING A NEW J OR OLD J. */
/*     JSNOLD   KEEPS TRACK OF NO. OF STEPS TAKEN WITH */
/*                PRESENT ITERATION MATRIX (BE IT FORMED BY */
/*                A NEW J OR NOT). */
/*     AVNEWJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION */
/*                MATRIX WAS FORMED BY A NEW J. */
/*     AVOLDJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION */
/*                MATRIX WAS FORMED BY AN OLD J. */
/*     NRENEW   FLAG THAT IS USED IN COMMUNICATION WITH SUBROUTINE PSET. */
/*                IF  NRENEW > 0  THEN FORM A NEW JACOBIAN BEFORE */
/*                                COMPUTING THE COEFFICIENT MATRIX FOR */
/*                                THE NEWTON-RAPHSON ITERATION */
/*                           = 0  FORM THE COEFFICIENT MATRIX USING A */
/*                                COPY OF AN OLD JACOBIAN */
/*     NEWPAR   FLAG USED IN THIS SUBROUTINE TO INDICATE IF A JACOBIAN */
/*              HAS BEEN EVALUATED FOR THE CURRENT STEP */
/* ********************************************************************** */
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER JSTART,KFLAG,LOUT,MF,N,ITOL >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION HSTPSZ(2,14) >*/
/*<    >*/
/*<       INTEGER IPIV(N), MBND(4) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<    >*/
/*<       LOGICAL FINISH,OVRIDE,WORKED >*/
/*     .. */
/*     .. LOCAL ARRAYS .. */
/*<       DIMENSION  EL(10),ELST(10),TQ(5) >*/
/*     .. */
/*     .. EXTERNAL SUBROUTINES .. */
/*<    >*/
/*     .. */
/*     .. INTRINSIC FUNCTIONS .. */
/*<       INTRINSIC DABS,DMAX1,DMIN1 >*/
/*     .. */
/*     .. COMMON BLOCKS .. */
/*<       COMMON/STPSZE/HSTPSZ >*/
/*<    >*/
/*<       LOGICAL CFAIL,JNEWIM,SAMPLE >*/
/*     .. */
/*     .. SAVE STATEMENT .. */
/*<    >*/
/*     .. */
/*     .. DATA STATEMENTS .. */
/*<       DATA  EL(2),ELST(2),OLDLO/3*1.0D+0/ >*/
    /* Parameter adjustments */
    --mbnd;
    --ipiv;
    --arh;
    ynhold_dim1 = *n;
    ynhold_offset = 1 + ynhold_dim1 * 1;
    ynhold -= ynhold_offset;
    yhold_dim1 = *n;
    yhold_offset = 1 + yhold_dim1 * 1;
    yhold -= yhold_offset;
    --scale;
    --save2;
    --save1;
    --error;
    --ymax;
    --yprime;
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    --pw;
    --pwcopy;
    --rtol;
    --atol;
    --rpar;
    --ipar;

    /* Function Body */
/*<       DATA  ZERO,ONE/0.0D+0,1.0D+0/ >*/
/*     .. */
/*<  6000 TOLD = T >*/
L6000:
    told = *t;
/*<       KFLAG = 0 >*/
    *kflag = 0;
/*<       IF (JSTART.GT.0) GO TO 60 >*/
    if (*jstart > 0) {
	goto L60;
    }
/*<       IF (JSTART.NE.0) GO TO 30 >*/
    if (*jstart != 0) {
	goto L30;
    }
/*     ------------------------------------------------------------------ */
/*     ON THE FIRST CALL, THE ORDER IS SET TO 1 AND THE INITIAL DERIVATIVE */
/*     IS CALCULATED OR GIVEN.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE */
/*     INCREASED IN A SINGLE STEP.  RMAX IS SET EQUAL TO 1.D4 INITIALLY */
/*     TO COMPENSATE FOR THE SMALL INITIAL H, BUT THEN IS NORMALLY = 10. */
/*     IF A FAILURE OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), */
/*     RMAX IS SET AT 2. FOR THE NEXT INCREASE. */
/*     ------------------------------------------------------------------ */
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y(I,2) = H*YPRIME(I) >*/
	y[i__ + (y_dim1 << 1)] = *h__ * yprime[i__];
/*<  10   CONTINUE >*/
/* L10: */
    }
/*<       METH = 2 >*/
    meth = 2;
/*<       MITER = MF - 10*METH >*/
    g_miter = *mf - meth * 10;
/*<       IBND=5 >*/
    g_ibnd = 5;
/*<       UPBND=0.2D+0 >*/
    g_upbnd = .2;
/*<       IF(MF.GT.22) UPBND=0.1D+0 >*/
    if (*mf > 22) {
	g_upbnd = .1;
    }
/*<       NQ = 1 >*/
    g_nq = 1;
/*<       NQUSED = NQ >*/
    *nqused = g_nq;
/*<       L = 2 >*/
    g_l = 2;
/*<       IDOUB = 3 >*/
    g_idoub = 3;
/*<       KFAIL = 0 >*/
    g_kfail = 0;
/*<       RMAX = 10000.0D+0 >*/
    g_stiff_rmax = 1e4;
/*<       IER=0 >*/
    ier = 0;
/*<       RC = ZERO >*/
    g_rc = zero;
/*<       CRATE1 = 0.1D+0 >*/
    g_crate1 = .1;
/*<       CRATE2 = 0.1D+0 >*/
    g_crate2 = .1;
/*<       JSNOLD = 0 >*/
    g_jsnold = 0;
/*<       JNEWIM = .TRUE. >*/
    g_jnewim = TRUE_;
/*<       TCRAT1 = ZERO >*/
    g_tcrat1 = zero;
/*<       TCRAT2 = ZERO >*/
    g_tcrat2 = zero;
/*<       VTOL=DMAX1(RTOL(1),ATOL(1))/10.0D+0 >*/
    vtol = max(rtol[1],atol[1]) / 10.;
/*<       DO 15 I=1,12 >*/
    for (i__ = 1; i__ <= 12; ++i__) {
/*<          HSTPSZ(1,I)=1.0D+0 >*/
	stpsze_1.hstpsz[(i__ << 1) - 2] = 1.;
/*<          HSTPSZ(2,I)=VTOL >*/
	stpsze_1.hstpsz[(i__ << 1) - 1] = vtol;
/*<  15   CONTINUE >*/
/* L15: */
    }
/*<       HOLD = H >*/
    g_hold = *h__;
/*<       MFOLD = MF >*/
    g_mfold = *mf;
/*<       NSTEP = 0 >*/
    *nstep = 0;
/*<       NRE = 1 >*/
    *nre = 1;
/*<       NJE = 0 >*/
    *nje = 0;
/*<       NDEC = 0 >*/
    *ndec = 0;
/*<       NPSET = 0 >*/
    *npset = 0;
/*<       NCOSET = 0 >*/
    *ncoset = 0;
/*<       MAXORD = 1 >*/
    *maxord = 1;
/*<       NFAIL = 0 >*/
    *nfail = 0;
/*<       CFAIL = .TRUE. >*/
    g_cfail = TRUE_;
/*<       AVNEWJ = ZERO >*/
    g_avnewj = zero;
/*<       AVOLDJ = ZERO >*/
    avoldj = zero;
/*<       AVNEW2 = ZERO >*/
    avnew2 = zero;
/*<       AVOLD2 = ZERO >*/
    avold2 = zero;
/*<       SAMPLE = .FALSE. >*/
    g_sample = FALSE_;
/*<       ISAMP = 0 >*/
    g_isamp = 0;
/*<       IEMB=0 >*/
    g_iemb = 0;
/*     ************************************************** */
/*     CFAIL=.TRUE. ENSURES THAT WE CALCULATE A NEW */
/*     J ON THE FIRST CALL. */
/*     ************************************************** */
/*<       MEQC1 = 0 >*/
    g_meqc1 = 0;
/*<       MEQC2 = 0 >*/
    g_meqc2 = 0;
/*<       MQ1TMP = 0 >*/
    g_mq1tmp = 0;
/*<       MQ2TMP = 0 >*/
    g_mq2tmp = 0;
/*<       NBSOL = 0 >*/
    *nbsol = 0;
/*<       HUSED = H >*/
    *hused = *h__;
/*     ----------------------------------------------------------------- */
/*     IF THE CALLER HAS CHANGED N , THE CONSTANTS E, EDN, EUP */
/*     AND BND MUST BE RESET.  E IS A COMPARISON FOR ERRORS AT THE */
/*     CURRENT ORDER NQ.  EUP IS TO TEST FOR INCREASING THE ORDER, */
/*     EDN FOR DECREASING THE ORDER.  BND IS USED TO TEST FOR CONVERGENCE */
/*     OF THE CORRECTOR ITERATES.   IF THE CALLER HAS CHANGED H, Y MUST */
/*     BE RE-SCALED.  IF H IS CHANGED, IDOUB IS SET TO L+1 TO PREVENT */
/*     FURTHER CHANGES IN H FOR THAT MANY STEPS. */
/*     ----------------------------------------------------------------- */
/*<       CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD) >*/
    coset_(&g_nq, el, elst, tq, ncoset, maxord);
/*<       LMAX = MAXDER + 1 >*/
    g_lmax = *maxder + 1;
/*<       RC = RC*EL(1)/OLDLO >*/
    g_rc = g_rc * el[0] / oldlo;
/*<       OLDLO = EL(1) >*/
    oldlo = el[0];
/*<       IWEVAL = MITER >*/
    g_iweval = g_miter;
/*<       NRENEW = 1 >*/
    nrenew = 1;
/*<       NEWPAR = 0 >*/
    newpar = 0;
/*     ***************************************************** */
/*     NRENEW AND NEWPAR ARE TO INSTRUCT ROUTINE THAT */
/*     WE WISH A NEW J TO BE CALCULATED FOR THIS STEP. */
/*     ***************************************************** */
/*<       CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN) >*/
    errors_(n, tq, &g_edn, &g_e, &g_eup, &g_bnd, &g_eddn);
/*<       DO 20 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ARH(I) = EL(2)*Y(I,1) >*/
	arh[i__] = el[1] * y[i__ + y_dim1];
/*<  20   CONTINUE >*/
/* L20: */
    }
/*<       CALL CPYARY(N*L,Y,YHOLD) >*/
    i__1 = *n * g_l;
    cpyary_(&i__1, &y[y_offset], &yhold[yhold_offset]);
/*<       QI = H*EL(1) >*/
    qi = *h__ * el[0];
/*<       QQ = ONE/QI >*/
    qq = one / qi;
/*<       CALL PRDICT(T,H,Y,L,N,IPAR,RPAR,IERR)       >*/
    prdict_(t, h__, &y[y_offset], &g_l, n, &ipar[1], &rpar[1], ierr);
/*<       IF(IERR.NE.0) THEN >*/
    if (*ierr != 0) {
/*<          H=H/2 >*/
	*h__ /= 2;
/*<          IERR = 0 >*/
	*ierr = 0;
/*<          GOTO 6000 >*/
	goto L6000;
/*<       ENDIF >*/
    }
/*<       DO 25 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YPRIME(I)=(Y(I,1)-ARH(I))/h >*/
	yprime[i__] = (y[i__ + y_dim1] - arh[i__]) / *h__;
/*<  25   CONTINUE >*/
/* L25: */
    }
/*<       GO TO 110 >*/
    goto L110;
/*     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/*     DIFFERENT PARAMETERS ON THIS CALL        < */
/*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*<  30   CALL CPYARY(N*L,YHOLD,Y) >*/
L30:
    i__1 = *n * g_l;
    cpyary_(&i__1, &yhold[yhold_offset], &y[y_offset]);
/*<       IF (MF.NE.MFOLD) THEN >*/
    if (*mf != g_mfold) {
/*<          METH = MF/10 >*/
	meth = *mf / 10;
/*<          MITER = MF - 10*METH >*/
	g_miter = *mf - meth * 10;
/*<          MFOLD = MF >*/
	g_mfold = *mf;
/*<          IWEVAL = MITER >*/
	g_iweval = g_miter;
/*<       END IF >*/
    }
/*<       IF(NSTEP.GT.0) GOTO 35 >*/
    if (*nstep > 0) {
	goto L35;
    }
/*<       NJE = 0 >*/
    *nje = 0;
/*<       NRE = 1 >*/
    *nre = 1;
/*<       CFAIL = .TRUE. >*/
    g_cfail = TRUE_;
/*<       NEWPAR = 0 >*/
    newpar = 0;
/*<       MQ1TMP = 0 >*/
    g_mq1tmp = 0;
/*<       MQ2TMP = 0 >*/
    g_mq2tmp = 0;
/*<       MEQC1 = 0 >*/
    g_meqc1 = 0;
/*<       MEQC2 = 0 >*/
    g_meqc2 = 0;
/*<       TCRAT1 = 0.0D+0 >*/
    g_tcrat1 = 0.;
/*<       TCRAT2 = 0.0D+0 >*/
    g_tcrat2 = 0.;
/*<       CRATE1 = 1.0D-1 >*/
    g_crate1 = .1;
/*<       CRATE2 = 1.0D-1 >*/
    g_crate2 = .1;
/*<       NSTEP = 0 >*/
    *nstep = 0;
/*<       NBSOL = 0 >*/
    *nbsol = 0;
/*<       NPSET = 0 >*/
    *npset = 0;
/*<       NCOSET = 0 >*/
    *ncoset = 0;
/*<       NDEC = 0 >*/
    *ndec = 0;
/*<  35   CONTINUE >*/
L35:
/*<       IF (H.NE.HOLD) THEN >*/
    if (*h__ != g_hold) {
/*<          RH = H/HOLD >*/
	g_rh = *h__ / g_hold;
/*<          H = HOLD >*/
	*h__ = g_hold;
/*<          IREDO = 3 >*/
	iredo = 3;
/*<          GO TO 50 >*/
	goto L50;
/*<       ELSE >*/
    } else {
/*<          GO TO 60 >*/
	goto L60;
/*<       END IF >*/
    }
/*     ********************************************* */
/*     RE-SCALE Y AFTER A CHANGE OF STEPSIZE   * */
/*     ********************************************* */
/*<  40   RH = DMAX1(RH,HMIN/DABS(H)) >*/
L40:
/* Computing MAX */
    d__1 = g_rh, d__2 = *hmin / abs(*h__);
    g_rh = max(d__1,d__2);
/*<  50   RH = DMIN1(RH,HMAX/DABS(H),RMAX) >*/
L50:
/* Computing MIN */
    d__1 = g_rh, d__2 = *hmax / abs(*h__), d__1 = min(d__1,d__2);
    g_rh = min(d__1,g_stiff_rmax);
/*<       CALL RSCALE(N,L,RH,Y) >*/
    rscale_(n, &g_l, &g_rh, &y[y_offset]);
/*<       RMAX = 10.0D+0 >*/
    g_stiff_rmax = 10.;
/*<       JCHANG = 1 >*/
    g_jchang = 1;
/*<       H = H*RH >*/
    *h__ *= g_rh;
/*<       RC = RC*RH >*/
    g_rc *= g_rh;
/*<       IF (JSNOLD.GT.IBND) THEN >*/
    if (g_jsnold > g_ibnd) {
/*<          CFAIL = .TRUE. >*/
	g_cfail = TRUE_;
/*<          NEWPAR = 0 >*/
	newpar = 0;
/*<          RC = ZERO >*/
	g_rc = zero;
/* ********************************************************************** */
/*        CFAIL=TRUE AND NEWPAR=0 SHOULD FORCE A NEW J TO BE EVALUATED */
/*        AFTER 7 STEPS WITH AN OLD J, IF WE HAVE HAD A FAILURE OF ANY */
/*        KIND ON THE FIRST, SECOND OR THIRD STAGE OF THE CURRENT STEP */
/* ********************************************************************** */
/*<       END IF >*/
    }
/*<       IDOUB = L + 1 >*/
    g_idoub = g_l + 1;
/*<       CALL CPYARY(N*L,Y,YHOLD) >*/
    i__1 = *n * g_l;
    cpyary_(&i__1, &y[y_offset], &yhold[yhold_offset]);
/*<  60   IF (DABS(RC-ONE).GT.UPBND) IWEVAL = MITER >*/
L60:
    if ((d__1 = g_rc - one, abs(d__1)) > g_upbnd) {
	g_iweval = g_miter;
    }
/*<       HUSED = H >*/
    *hused = *h__;
/*     ------------------------------------------------------------------ */
/*     THIS SECTION COMPUTES THE PREDICTED VALUES OF Y */
/*     AND THE RHS, ARH, FOR USE IN THE NEWTON ITERATION SCHEME. */
/*     RC IS THE RATIO OF THE NEW TO OLD VALUES OF THE COEFFICIENT */
/*     H*EL(1). WHEN RC DIFFERS FROM 1 BY MORE THAN 20 PERCENT, IWEVAL IS */
/*     SET TO MITER TO FORCE THE PARTIALS TO BE UPDATED. */
/*     ------------------------------------------------------------------ */
/*<       QI = H*EL(1) >*/
    qi = *h__ * el[0];
/*<       QQ = ONE/QI >*/
    qq = one / qi;
/*<       DO 70 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ARH(I) = EL(2)*Y(I,1) >*/
	arh[i__] = el[1] * y[i__ + y_dim1];
/*<  70   CONTINUE >*/
/* L70: */
    }
/*<       DO 90 J1 = 2,NQ >*/
    i__1 = g_nq;
    for (j1 = 2; j1 <= i__1; ++j1) {
/*<          JP1 =J1+1 >*/
	jp1 = j1 + 1;
/*<          DO 80 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             ARH(I) = ARH(I) + EL(JP1)*Y(I,J1) >*/
	    arh[i__] += el[jp1 - 1] * y[i__ + j1 * y_dim1];
/*<  80      CONTINUE >*/
/* L80: */
	}
/*<  90   CONTINUE >*/
/* L90: */
    }
/*<       IF (JCHANG.EQ.1) THEN >*/
    if (g_jchang == 1) {
/*        IF WE HAVE CHANGED STEPSIZE THEN PREDICT A VALUE FOR Y(T+H) */
/*        AND EVALUATE THE DERIVATIVE THERE (STORED IN SAVE2()) */
/*<          CALL PRDICT(T,H,Y,L,N,IPAR,RPAR,IERR) >*/
	prdict_(t, h__, &y[y_offset], &g_l, n, &ipar[1], &rpar[1], ierr);
/*<          IF(IERR.NE.0) GOTO 8000 >*/
	if (*ierr != 0) {
	    goto L8000;
	}
/*<          DO 95 I=1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             YPRIME(I)=(Y(I,1)-ARH(I))/QI >*/
	    yprime[i__] = (y[i__ + y_dim1] - arh[i__]) / qi;
/*<  95      CONTINUE >*/
/* L95: */
	}
/*<       ELSE >*/
    } else {
/*        ELSE USE THE VALUES COMPUTED FOR THE SECOND BDF FROM THE LAST */
/*        STEP. Y( ,LMAX+2) HOLDS THE VALUE FOR THE DERIVATIVE AT (T+H) */
/*        AND Y( ,LMAX+3) HOLDS THE APPROXIMATION TO Y AT THIS POINT. */
/*<          LMP2=LMAX+2 >*/
	lmp2 = g_lmax + 2;
/*<          LMP3=LMAX+3 >*/
	lmp3 = g_lmax + 3;
/*<          DO 100 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             Y(I,1)=Y(I,LMP3)                       >*/
	    y[i__ + y_dim1] = y[i__ + lmp3 * y_dim1];
/*<             YPRIME(I) = (Y(I,1)-ARH(I))/QI >*/
	    yprime[i__] = (y[i__ + y_dim1] - arh[i__]) / qi;
/*<  100     CONTINUE          >*/
/* L100: */
	}
/*<          T = T + H >*/
	*t += *h__;
/*<       END IF >*/
    }
/*<  110  IF (IWEVAL.LE.0) GO TO 120 >*/
L110:
    if (g_iweval <= 0) {
	goto L120;
    }
/* ------------------------------------------------------------------- */
/*     IF INDICATED, THE MATRIX P = I/(H*EL(2)) - J IS RE-EVALUATED BEFORE */
/*     STARTING THE CORRECTOR ITERATION.  IWEVAL IS SET = 0 TO INDICATE */
/*     THAT THIS HAS BEEN DONE. P IS COMPUTED AND PROCESSED IN PSET. */
/*     THE PROCESSED MATRIX IS STORED IN PW */
/* ------------------------------------------------------------------- */
/*<       IWEVAL = 0 >*/
    g_iweval = 0;
/*<       RC = ONE >*/
    g_rc = one;
/*<       IITER = MEQC1 - MQ1TMP >*/
    iiter = g_meqc1 - g_mq1tmp;
/*<       IITER2 = MEQC2 - MQ2TMP >*/
    iiter2 = g_meqc2 - g_mq2tmp;
/*<       IF (JNEWIM) THEN >*/
    if (g_jnewim) {
/*<          IF (JSNOLD.GE.3) THEN >*/
	if (g_jsnold >= 3) {
/*<             AVNEWJ = TCRAT1/DBLE(FLOAT(IITER)) >*/
	    g_avnewj = g_tcrat1 / (doublereal) iiter;
/*<             AVNEW2 = TCRAT2/DBLE(FLOAT(IITER2)) >*/
	    avnew2 = g_tcrat2 / (doublereal) iiter2;
/*<          ELSE >*/
	} else {
/*<             AVNEWJ = ONE >*/
	    g_avnewj = one;
/*<             AVNEW2 = ONE >*/
	    avnew2 = one;
/*<          END IF >*/
	}
/*<       ELSE >*/
    } else {

/*          MATRIX P WAS FORMED WITH A COPY OF J */

/*<          IF (JSNOLD.GE.3) THEN >*/
	if (g_jsnold >= 3) {
/*<             AVOLDJ = TCRAT1/DBLE(FLOAT(IITER)) >*/
	    avoldj = g_tcrat1 / (doublereal) iiter;
/*<             AVOLD2 = TCRAT2/DBLE(FLOAT(IITER2)) >*/
	    avold2 = g_tcrat2 / (doublereal) iiter2;
/*<             IF (AVOLDJ.LT.AVNEWJ) THEN >*/
	    if (avoldj < g_avnewj) {
/*<                AVNEWJ = AVOLDJ >*/
		g_avnewj = avoldj;
/*<    >*/
	    } else if ((d__1 = avoldj - g_avnewj, abs(d__1)) > .3 || avoldj > 
		    .85 && avoldj != one) {

/*              SINCE IN CERTAIN INSTANCES AVOLDJ WILL */
/*              BE 1.0 AND THERE WILL BE NO NEED TO */
/*              UPDATE J. */

/*<                CFAIL = .TRUE. >*/
		g_cfail = TRUE_;
/*<                CRATE1 = 0.1D+0 >*/
		g_crate1 = .1;
/*<                CRATE2 = 0.1D+0 >*/
		g_crate2 = .1;
/*<             END IF >*/
	    }
/*<          ELSE >*/
	} else {
/*<             CFAIL = .TRUE. >*/
	    g_cfail = TRUE_;
/*<             CRATE1 = 0.1D+0 >*/
	    g_crate1 = .1;
/*<             CRATE2 = 0.1D+0 >*/
	    g_crate2 = .1;

/*           ********************************************************* */
/*           IF WE HAVE REACHED HERE THINGS MUST HAVE GONE WRONG */
/*           ********************************************************* */

/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       TCRAT1 = ZERO >*/
    g_tcrat1 = zero;
/*<       TCRAT2 = ZERO >*/
    g_tcrat2 = zero;
/*<       IF (CFAIL) THEN >*/
    if (g_cfail) {
/*<          NRENEW = 1 >*/
	nrenew = 1;
/*<          NEWPAR = 1 >*/
	newpar = 1;
/*<          JSINUP = -1 >*/
	g_jsinup = -1;
/*<          JNEWIM = .TRUE. >*/
	g_jnewim = TRUE_;
/*<       END IF >*/
    }
/*<       CFAIL = .FALSE. >*/
    g_cfail = FALSE_;
/*<       JSNOLD = 0 >*/
    g_jsnold = 0;
/*<       MQ1TMP = MEQC1 >*/
    g_mq1tmp = g_meqc1;
/*<       MQ2TMP = MEQC2       >*/
    g_mq2tmp = g_meqc2;
/*<    >*/
    pset_(&y[y_offset], &yprime[1], n, h__, t, uround, epsjac, &qi, &g_miter, &
	    mbnd[1], nind1, nind2, nind3, &ier, (S_fp)pderv, (S_fp)resid, &
	    nrenew, &ymax[1], &save1[1], &save2[1], &scale[1], &pw[1], &
	    pwcopy[1], &error[1], &ipiv[1], itol, &rtol[1], &atol[1], npset, 
	    nje, nre, ndec, &ipar[1], &rpar[1], ierr);
/*<       IF(IERR.NE.0) GOTO 8000 >*/
    if (*ierr != 0) {
	goto L8000;
    }
/*<       QQQ=QI >*/
    g_qqq = qi;

/*     NOTE THAT ERROR() IS JUST BEING USED AS A WORKSPACE BY PSET */
/*<       IF (IER.NE.0) THEN >*/
    if (ier != 0) {
/*     IF IER>0 THEN WE HAVE HAD A SINGULARITY IN THE ITERATION MATRIX */
/*<          IJUS=1 >*/
	g_ijus = 1;
/*<          RED=0.5D+0 >*/
	red = .5;
/*<          NFAIL = NFAIL + 1 >*/
	++(*nfail);
/*<          GO TO 450 >*/
	goto L450;
/*<       END IF >*/
    }
/*<  120  DO 130 I = 1,N >*/
L120:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SAVE1(I) = Y(I,1) >*/
	save1[i__] = y[i__ + y_dim1];
/*<          ERROR(I) = ZERO >*/
	error[i__] = zero;
/*<  130  CONTINUE >*/
/* L130: */
    }
/*<       M1 = 0 >*/
    m1 = 0;
/* ********************************************************************** */
/*     UP TO 4 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS MADE */
/*     ON THE R.M.S. NORM OF EACH CORRECTION ,USING BND, WHICH DEPENDS */
/*     ON ATOL AND RTOL.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE */
/*     VECTOR  ERROR(I).  THE Y ARRAY IS NOT ALTERED IN THE CORRECTOR */
/*     LOOP. THE UPDATED Y VECTOR IS STORED TEMPORARILY IN SAVE1. */
/* ********************************************************************** */
/*<       IF (.NOT.SAMPLE) THEN >*/
    if (! g_sample) {
/*<    >*/
	itrat2_(&g_qqq, &y[y_offset], &yprime[1], n, t, &qi, &g_bnd, &arh[1], &
		g_crate1, &g_tcrat1, &m1, &worked, &ymax[1], &error[1], &save1[1],
		 &save2[1], &scale[1], &pw[1], mf, &mbnd[1], nind1, nind2, 
		nind3, &ipiv[1], &c__1, itol, &rtol[1], &atol[1], &ipar[1], &
		rpar[1], hused, nbsol, nre, nqused, (S_fp)resid, ierr);
/*<          IF(IERR.NE.0) GOTO 8000 >*/
	if (*ierr != 0) {
	    goto L8000;
	}
/*<       ELSE >*/
    } else {
/*<    >*/
	itrat2_(&g_qqq, &y[y_offset], &yprime[1], n, t, &qi, &g_bnd, &arh[1], &
		g_crate1, &g_tcrat1, &m1, &worked, &ymax[1], &error[1], &save1[1],
		 &save2[1], &scale[1], &pw[1], mf, &mbnd[1], nind1, nind2, 
		nind3, &ipiv[1], &c__0, itol, &rtol[1], &atol[1], &ipar[1], &
		rpar[1], hused, nbsol, nre, nqused, (S_fp)resid, ierr);
/*<          IF(IERR.NE.0) GOTO 8000 >*/
	if (*ierr != 0) {
	    goto L8000;
	}
/*<       END IF >*/
    }
/*<       MEQC1 = MEQC1 + M1 + 1 >*/
    g_meqc1 = g_meqc1 + m1 + 1;

/*       NOW TEST TO SEE IF IT WAS SUCCESSFUL OR NOT */


/*<       IF (.NOT.WORKED) THEN >*/
    if (! worked) {
/*<          NFAIL = NFAIL + 1 >*/
	++(*nfail);
/* ********************************************************************** */
/*        THE CORRECTOR ITERATION FAILED TO CONVERGE IN 4 TRIES. IF */
/*        PARTIALS ARE NOT UP TO DATE, THEY ARE RE-EVALUATED FOR THE */
/*        NEXT TRY. OTHERWISE THE Y ARRAY IS REPLACED BY ITS VALUES */
/*        BEFORE PREDICTION AND H IS REDUCED IF POSSIBLE. IF NOT A */
/*        NON-CONVERGENCE EXIT IS TAKEN */
/* ********************************************************************** */
/*<          IF (IWEVAL.EQ.-1) THEN >*/
	if (g_iweval == -1) {
/*           HAVE BEEN USING OLD PARTIALS, UPDATE THEM AND TRY AGAIN */
/*<             IWEVAL = MITER >*/
	    g_iweval = g_miter;
/*<             CFAIL = .TRUE. >*/
	    g_cfail = TRUE_;
/*<             do 135 i=1,n >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                yprime(i)=(y(I,1)-arh(i))*qq >*/
		yprime[i__] = (y[i__ + y_dim1] - arh[i__]) * qq;
/*<  135        continue >*/
/* L135: */
	    }
/*<             GO TO 110 >*/
	    goto L110;
/*<          END IF >*/
	}
/*<          IJUS=0 >*/
	g_ijus = 0;
/*<          RED=0.5D+0 >*/
	red = .5;
/*    ***    failed at step 1 because of Newton */
/*<          GO TO 450 >*/
	goto L450;
/*<       END IF >*/
    }
/*<       IWEVAL = -1 >*/
    g_iweval = -1;
/*<       HUSED = H >*/
    *hused = *h__;
/*<       NQUSED = NQ >*/
    *nqused = g_nq;
/*<       DO 140 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y(I,1) = (SAVE1(I)-ARH(I)) >*/
	y[i__ + y_dim1] = save1[i__] - arh[i__];
/*<  140  CONTINUE >*/
/* L140: */
    }
/*<       DO 145 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SAVE2(I) = Y(I,1)*QQ >*/
	save2[i__] = y[i__ + y_dim1] * qq;
/*<          Y(I,1) = SAVE1(I) >*/
	y[i__ + y_dim1] = save1[i__];
/*<  145  CONTINUE >*/
/* L145: */
    }

/*     UPDATE THE DIFFERENCES AT N+1 */

/*<       DO 160 J = 2,L >*/
    i__1 = g_l;
    for (j = 2; j <= i__1; ++j) {
/*<          JM1 = J-1 >*/
	jm1 = j - 1;
/*<          DO 150 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             Y(I,J) = Y(I,JM1) - YHOLD(I,JM1) >*/
	    y[i__ + j * y_dim1] = y[i__ + jm1 * y_dim1] - yhold[i__ + jm1 * 
		    yhold_dim1];
/*<  150     CONTINUE >*/
/* L150: */
	}
/*<  160  CONTINUE >*/
/* L160: */
    }

/*     COMPUTE ERROR IN THE SOLUTION */

/*<       DO 161 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          AYI = DABS(Y(I,1)) >*/
	ayi = (d__1 = y[i__ + y_dim1], abs(d__1));
/*<          IF(ITOL.EQ.1) THEN >*/
	if (*itol == 1) {
/*<             SCALE(I) = YMAX(I) >*/
	    scale[i__] = ymax[i__];
/*<          ELSE IF(ITOL.EQ.2) THEN >*/
	} else if (*itol == 2) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[1] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.3) THEN >*/
	} else if (*itol == 3) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[1] * ayi + atol[i__];
/*<          ELSE IF(ITOL.EQ.4) THEN >*/
	} else if (*itol == 4) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[i__] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.5) THEN >*/
	} else if (*itol == 5) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[i__] * ayi + atol[i__];
/*<          ENDIF >*/
	}
/*<  161  CONTINUE >*/
/* L161: */
    }
/*<       IF(NIND2.NE.0) THEN >*/
    if (*nind2 != 0) {
/*<          DO 162 I = NIND1+1,NIND2+NIND1 >*/
	i__1 = *nind2 + *nind1;
	for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/HUSED >*/
	    scale[i__] /= *hused;
/*<  162     CONTINUE >*/
/* L162: */
	}
/*<       ENDIF >*/
    }
/*<       IF(NIND3.NE.0) THEN >*/
    if (*nind3 != 0) {
/*<          DO 163 I = NIND1 +NIND2 + 1,NIND1+NIND2+NIND3 >*/
	i__1 = *nind1 + *nind2 + *nind3;
	for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/(HUSED**2) >*/
/* Computing 2nd power */
	    d__1 = *hused;
	    scale[i__] /= d__1 * d__1;
/*<  163     CONTINUE >*/
/* L163: */
	}
/*<       ENDIF >*/
    }
/* **** */
/* ****  AMMEND */
/* ****  CHANGE 1,N BELOW TO 1,NVARS */
/* **** */
/*<       D = ZERO >*/
    d__ = zero;
/*<       DO 170 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          D = D + ((Y(I,L)-YHOLD(I,L))/SCALE(I))**2 >*/
/* Computing 2nd power */
	d__1 = (y[i__ + g_l * y_dim1] - yhold[i__ + g_l * yhold_dim1]) / scale[
		i__];
	d__ += d__1 * d__1;
/*<  170  CONTINUE >*/
/* L170: */
    }

/*    STORING Y FROM FIRST STEP FOR USE IN THIRD STEP. */

/*<       IF(ITOL .EQ. 1) D = D/(RTOL(1)**2) >*/
    if (*itol == 1) {
/* Computing 2nd power */
	d__1 = rtol[1];
	d__ /= d__1 * d__1;
    }
/*<       IF(D.GT.E) GOTO 330 >*/
    if (d__ > g_e) {
	goto L330;
    }
/*<       DO 180 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YNHOLD(I,1) = Y(I,1) >*/
	ynhold[i__ + ynhold_dim1] = y[i__ + y_dim1];
/*<          YNHOLD(I,2) = SAVE2(I) >*/
	ynhold[i__ + (ynhold_dim1 << 1)] = save2[i__];
/*<  180  CONTINUE >*/
/* L180: */
    }
/*<       KFAIL = 0 >*/
    g_kfail = 0;
/*<       IREDO = 0 >*/
    iredo = 0;
/* ---------------------------------------------------------------------- */
/*<       DO 190 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         ARH(I) = EL(2)*Y(I,1) >*/
	arh[i__] = el[1] * y[i__ + y_dim1];
/*<  190  CONTINUE >*/
/* L190: */
    }
/*<       DO 210 J1 = 2,NQ >*/
    i__1 = g_nq;
    for (j1 = 2; j1 <= i__1; ++j1) {
/*<          JP1 = J1+1 >*/
	jp1 = j1 + 1;
/*<          DO 200 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             ARH(I) = ARH(I) + EL(JP1)*Y(I,J1) >*/
	    arh[i__] += el[jp1 - 1] * y[i__ + j1 * y_dim1];
/*<  200     CONTINUE >*/
/* L200: */
	}
/*<  210  CONTINUE >*/
/* L210: */
    }
/*<       CALL PRDICT(T,H,Y,L,N,IPAR,RPAR,IERR)       >*/
    prdict_(t, h__, &y[y_offset], &g_l, n, &ipar[1], &rpar[1], ierr);
/*<       IF(IERR.NE.0) GOTO 8000 >*/
    if (*ierr != 0) {
	goto L8000;
    }
/*<       DO 215 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YPRIME(I)=(Y(I,1)-ARH(I))/QQQ >*/
	yprime[i__] = (y[i__ + y_dim1] - arh[i__]) / g_qqq;
/*<  215  CONTINUE       >*/
/* L215: */
    }
/*<       DO 220 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SAVE1(I) = Y(I,1) >*/
	save1[i__] = y[i__ + y_dim1];
/*<          ERROR(I) = ZERO >*/
	error[i__] = zero;
/*<  220  CONTINUE >*/
/* L220: */
    }
/*<       M2 = 0 >*/
    m2 = 0;

/*     FOR NOW WILL ASSUME THAT WE DO NOT WISH TO SAMPLE */
/*     AT THE N+2 STEP POINT */

/*<    >*/
    itrat2_(&g_qqq, &y[y_offset], &yprime[1], n, t, &qi, &g_bnd, &arh[1], &g_crate2,
	     &g_tcrat2, &m2, &worked, &ymax[1], &error[1], &save1[1], &save2[1],
	     &scale[1], &pw[1], mf, &mbnd[1], nind1, nind2, nind3, &ipiv[1], &
	    c__1, itol, &rtol[1], &atol[1], &ipar[1], &rpar[1], hused, nbsol, 
	    nre, nqused, (S_fp)resid, ierr);
/*<       IF(IERR.NE.0) GOTO 8000 >*/
    if (*ierr != 0) {
	goto L8000;
    }
/*<       MEQC2 = MEQC2 + M2 + 1 >*/
    g_meqc2 = g_meqc2 + m2 + 1;

/*       NOW CHECK TO SEE IF IT WAS SUCCESSFUL OR NOT */

/*<       IF (.NOT.WORKED) THEN >*/
    if (! worked) {
/*<          NFAIL = NFAIL + 1 >*/
	++(*nfail);
/*<          IJUS=0 >*/
	g_ijus = 0;
/*<          RED=0.5D+0 >*/
	red = .5;
/* ***have failed on step 2 */
/*<          GOTO 450 >*/
	goto L450;
/*<       END IF >*/
    }

/*        IF WE ARE DOWN TO HERE THEN THINGS MUST HAVE CONVERGED */

/*<       LMP2=LMAX+2 >*/
    lmp2 = g_lmax + 2;
/*<       LMP3=LMAX+3 >*/
    lmp3 = g_lmax + 3;
/*<       DO 230 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y(I,LMP3) = (SAVE1(I)-ARH(I)) >*/
	y[i__ + lmp3 * y_dim1] = save1[i__] - arh[i__];
/*<  230  CONTINUE >*/
/* L230: */
    }
/*<       DO 233 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y(I,LMP2) = Y(I,LMP3)*QQ >*/
	y[i__ + lmp2 * y_dim1] = y[i__ + lmp3 * y_dim1] * qq;
/*<          Y(I,LMP3) = SAVE1(I) >*/
	y[i__ + lmp3 * y_dim1] = save1[i__];
/*<  233  CONTINUE >*/
/* L233: */
    }

/*     WE ARE NOW COMPUTING THE THIRD STAGE */

/*<       LL = L + 1 >*/
    ll = g_l + 1;
/*<       T = TOLD + H >*/
    *t = told + *h__;
/*<       DELST = ELST(1)-EL(1) >*/
    delst = elst[0] - el[0];
/*<       NQP2 = NQ+2 >*/
    nqp2 = g_nq + 2;
/*<       LMP2=LMAX+2 >*/
    lmp2 = g_lmax + 2;
/*<       DO 280 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          ARH(I) = H*(ELST(NQP2)*Y(I,LMP2)+DELST*YNHOLD(I,2)) >*/
	arh[i__] = *h__ * (elst[nqp2 - 1] * y[i__ + lmp2 * y_dim1] + delst * 
		ynhold[i__ + (ynhold_dim1 << 1)]);
/*<          DO 270 J1 = 1,NQ >*/
	i__2 = g_nq;
	for (j1 = 1; j1 <= i__2; ++j1) {
/*<             ARH(I) = ARH(I) + ELST(J1+1)*YHOLD(I,J1) >*/
	    arh[i__] += elst[j1] * yhold[i__ + j1 * yhold_dim1];
/*<  270     CONTINUE >*/
/* L270: */
	}
/*<  280  CONTINUE >*/
/* L280: */
    }
/*<       DO 290 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          SAVE2(I) = YNHOLD(I,2) >*/
	save2[i__] = ynhold[i__ + (ynhold_dim1 << 1)];
/*<          Y(I,1) = YNHOLD(I,1) >*/
	y[i__ + y_dim1] = ynhold[i__ + ynhold_dim1];
/*<  290  CONTINUE >*/
/* L290: */
    }
/*<       M3STEP = 0 >*/
    m3step = 0;
/*<  300  CONTINUE >*/
L300:
/*<       DO 310 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          YPRIME(I)=(Y(I,1)-ARH(I))/QQQ >*/
	yprime[i__] = (y[i__ + y_dim1] - arh[i__]) / g_qqq;
/*<  310  CONTINUE >*/
/* L310: */
    }
/*<       CALL RESID(N,T,Y,SAVE1,YPRIME,IPAR,RPAR,IERR) >*/
    (*resid)(n, t, &y[y_offset], &save1[1], &yprime[1], &ipar[1], &rpar[1], 
	    ierr);
/*<       if (ierr .ne. 0) go to 8000 >*/
    if (*ierr != 0) {
	goto L8000;
    }
/*<       NRE=NRE+1 >*/
    ++(*nre);

/*<       IF (MF.GE. 23) THEN >*/
    if (*mf >= 23) {
/*<          CALL DGBSL(PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE1,0) >*/
	dgbsl_(&pw[1], &mbnd[4], n, &mbnd[1], &mbnd[2], &ipiv[1], &save1[1], &
		c__0);
/*<          NBSOL=NBSOL+1 >*/
	++(*nbsol);
/*<       ELSE >*/
    } else {
/*<          CALL SOL(N,N,PW,SAVE1,IPIV) >*/
	sol_(n, n, &pw[1], &save1[1], &ipiv[1]);
/*<          NBSOL = NBSOL + 1 >*/
	++(*nbsol);
/*<       ENDIF >*/
    }
/*<       DO 321 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          AYI = DABS(Y(I,1)) >*/
	ayi = (d__1 = y[i__ + y_dim1], abs(d__1));
/*<          IF(ITOL.EQ.1) THEN >*/
	if (*itol == 1) {
/*<             SCALE(I) = YMAX(I) >*/
	    scale[i__] = ymax[i__];
/*<          ELSE IF(ITOL.EQ.2) THEN >*/
	} else if (*itol == 2) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[1] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.3) THEN >*/
	} else if (*itol == 3) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[1] * ayi + atol[i__];
/*<          ELSE IF(ITOL.EQ.4) THEN >*/
	} else if (*itol == 4) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[i__] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.5) THEN >*/
	} else if (*itol == 5) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[i__] * ayi + atol[i__];
/*<          ENDIF >*/
	}
/*<  321  CONTINUE >*/
/* L321: */
    }
/*<       IF(NIND2.NE.0) THEN >*/
    if (*nind2 != 0) {
/*<          DO 322 I = NIND1+1,NIND2+NIND1 >*/
	i__1 = *nind2 + *nind1;
	for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/HUSED >*/
	    scale[i__] /= *hused;
/*<  322     CONTINUE >*/
/* L322: */
	}
/*<       ENDIF >*/
    }
/*<       IF(NIND3.NE.0) THEN >*/
    if (*nind3 != 0) {
/*<          DO 133 I = NIND1+NIND2 + 1,NIND1+NIND2+NIND3 >*/
	i__1 = *nind1 + *nind2 + *nind3;
	for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/(HUSED**2) >*/
/* Computing 2nd power */
	    d__1 = *hused;
	    scale[i__] /= d__1 * d__1;
/*<  133     CONTINUE >*/
/* L133: */
	}
/*<       ENDIF >*/
    }
/*<       D = ZERO >*/
    d__ = zero;
/*<       DO 320 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         D = D + (SAVE1(I)/SCALE(I))**2 >*/
/* Computing 2nd power */
	d__1 = save1[i__] / scale[i__];
	d__ += d__1 * d__1;
/*<         Y(I,1) = Y(I,1) - SAVE1(I) >*/
	y[i__ + y_dim1] -= save1[i__];
/*<  320  CONTINUE >*/
/* L320: */
    }
/*<       IF(ITOL .EQ. 1) D = D/(RTOL(1)**2) >*/
    if (*itol == 1) {
/* Computing 2nd power */
	d__1 = rtol[1];
	d__ /= d__1 * d__1;
    }
/*<       IF ((D*DMIN1(ONE,2.0D+0*CRATE1)).LE.BND) GO TO 360 >*/
/* Computing MIN */
    d__1 = one, d__2 = g_crate1 * 2.;
    if (d__ * min(d__1,d__2) <= g_bnd) {
	goto L360;
    }
/*<       IF (M3STEP.EQ.4) THEN >*/
    if (m3step == 4) {
/*         WRITE (LOUT,9000) */
/*<          IJUS=1 >*/
	g_ijus = 1;
/*<          RED=0.5D+0 >*/
	red = .5;
/*    ****  step 3 fails */
/*<          NFAIL = NFAIL + 1 >*/
	++(*nfail);
/*<          GO TO 450 >*/
	goto L450;
/*<       END IF >*/
    }
/*<       M3STEP = M3STEP + 1 >*/
    ++m3step;
/*      IF(IERR.NE.0) GOTO 8000 */
/*<       GO TO 300 >*/
    goto L300;
/*<  330  KFAIL = KFAIL - 1 >*/
L330:
    --g_kfail;
/* ********************************************************************** */
/*     THE ERROR TEST FAILED. KFAIL KEEPS TRACK OF MULTIPLE FAILURES. */
/*     RESTORE T AND THE Y ARRAY TO THEIR PREVIOUS VALUES AND PREPARE TO */
/*     TRY THE STEP AGAIN. COMPUTE THE OPTIMAL STEP SIZE FOR THIS ORDER */
/*     AND ONE ORDER LOWER. */
/* ********************************************************************** */
/*     ***  failed on step 1 because of accuracy */
/*     COMPUTE ERROR IN THE SOLUTION */

/*<       NFAIL = NFAIL + 1 >*/
    ++(*nfail);
/*<       DO 561 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          AYI = DABS(Y(I,1)) >*/
	ayi = (d__1 = y[i__ + y_dim1], abs(d__1));
/*<          IF(ITOL.EQ.1) THEN >*/
	if (*itol == 1) {
/*<             SCALE(I) = YMAX(I) >*/
	    scale[i__] = ymax[i__];
/*<          ELSE IF(ITOL.EQ.2) THEN >*/
	} else if (*itol == 2) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[1] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.3) THEN >*/
	} else if (*itol == 3) {
/*<             SCALE(I) = RTOL(1)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[1] * ayi + atol[i__];
/*<          ELSE IF(ITOL.EQ.4) THEN >*/
	} else if (*itol == 4) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(1) >*/
	    scale[i__] = rtol[i__] * ayi + atol[1];
/*<          ELSE IF(ITOL.EQ.5) THEN >*/
	} else if (*itol == 5) {
/*<             SCALE(I) = RTOL(I)*AYI + ATOL(I) >*/
	    scale[i__] = rtol[i__] * ayi + atol[i__];
/*<          ENDIF >*/
	}
/*<  561  CONTINUE >*/
/* L561: */
    }
/*<       IF(NIND2.NE.0) THEN >*/
    if (*nind2 != 0) {
/*<          DO 562 I = NIND1+1,NIND2+NIND1 >*/
	i__1 = *nind2 + *nind1;
	for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/HUSED >*/
	    scale[i__] /= *hused;
/*<  562     CONTINUE >*/
/* L562: */
	}
/*<       ENDIF >*/
    }
/*<       IF(NIND3.NE.0) THEN >*/
    if (*nind3 != 0) {
/*<          DO 563 I = NIND1+NIND2 + 1,NIND1+NIND2+NIND3 >*/
	i__1 = *nind1 + *nind2 + *nind3;
	for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
/*<             SCALE(I)=SCALE(I)/(HUSED**2) >*/
/* Computing 2nd power */
	    d__1 = *hused;
	    scale[i__] /= d__1 * d__1;
/*<  563     CONTINUE >*/
/* L563: */
	}
/*<       ENDIF >*/
    }
/*<       DDOWN = ZERO >*/
    ddown = zero;
/*<       TWODWN = ZERO >*/
    twodwn = zero;
/*<       DO 1700 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          DDOWN = DDOWN + ((Y(I,L))/SCALE(I))**2 >*/
/* Computing 2nd power */
	d__1 = y[i__ + g_l * y_dim1] / scale[i__];
	ddown += d__1 * d__1;
/*<          TWODWN = TWODWN + ((Y(I,l-1))/SCALE(I))**2 >*/
/* Computing 2nd power */
	d__1 = y[i__ + (g_l - 1) * y_dim1] / scale[i__];
	twodwn += d__1 * d__1;
/*<  1700 CONTINUE >*/
/* L1700: */
    }
/*<       IF(ITOL .EQ. 1) D = D/(RTOL(1)**2) >*/
    if (*itol == 1) {
/* Computing 2nd power */
	d__1 = rtol[1];
	d__ /= d__1 * d__1;
    }
/*<       T = TOLD >*/
    *t = told;
/*<       HOLD = H >*/
    g_hold = *h__;
/*<       IF(NQ.GT.1) FFAIL = 0.5D+0/DBLE(FLOAT(NQ)) >*/
    if (g_nq > 1) {
	g_ffail = .5 / (doublereal) g_nq;
    }
/*<       IF(NQ.GT.2) FRFAIL = 0.5D+0/DBLE(FLOAT(NQ-1)) >*/
    if (g_nq > 2) {
	frfail = .5 / (doublereal) (g_nq - 1);
    }
/*<       EFAIL = 0.5D+0/DBLE(FLOAT(L)) >*/
    efail = .5 / (doublereal) g_l;
/*<       CALL CPYARY(N*L,YHOLD,Y) >*/
    i__1 = *n * g_l;
    cpyary_(&i__1, &yhold[yhold_offset], &y[y_offset]);
/*<       RMAX = 2.0D+0 >*/
    g_stiff_rmax = 2.;
/*<       IF (DABS(H).LE.HMIN*1.00001D+0) THEN >*/
    if (abs(*h__) <= *hmin * 1.00001) {

/*        REQUESTED ERROR NOT POSSIBLE WITH GIVEN HMIN */

/*<          KFLAG = -1 >*/
	*kflag = -1;
/*<          HOLD = H >*/
	g_hold = *h__;
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<       IF (KFAIL.LE.-3) GO TO 340 >*/
    if (g_kfail <= -3) {
	goto L340;
    }
/*<       IREDO = 2 >*/
    iredo = 2;

/*     PREDICTING A NEW H AFTER INSUFFICIENT ACCURACY */

/*<       PRFAIL = ((D/(0.2D+0*E))**EFAIL)*1.5D+0 + 1.6D-6 >*/
    d__1 = d__ / (g_e * .2);
    prfail = pow_dd(&d__1, &efail) * 1.5 + 1.6e-6;
/*<       PLFAIL = ((DDOWN/(0.2D+0*EDN))**FFAIL)*1.5D+0+1.7D-6 >*/
    d__1 = ddown / (g_edn * .2);
    plfail = pow_dd(&d__1, &g_ffail) * 1.5 + 1.7e-6;
/*<    >*/
    if (g_nq > 2) {
	d__1 = twodwn / (g_eddn * .2);
	g_pllfal = pow_dd(&d__1, &frfail) * 1.5 + 1.7e-6;
    }
/*<       IF(PLLFAL.GT.PLFAIL) PLFAIL=PLLFAL >*/
    if (g_pllfal > plfail) {
	plfail = g_pllfal;
    }
/*<       IF(PLFAIL.LT.PRFAIL.AND.NQ.NE.1) THEN >*/
    if (plfail < prfail && g_nq != 1) {
/*<          NEWQ=NQ-1 >*/
	newq = g_nq - 1;
/*<          NQ=NEWQ >*/
	g_nq = newq;
/*<          RH=ONE/(PLFAIL*DBLE(FLOAT(-KFAIL))) >*/
	g_rh = one / (plfail * (doublereal) (-g_kfail));
/*<          L=NQ+1 >*/
	g_l = g_nq + 1;
/*<          CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD) >*/
	coset_(&g_nq, el, elst, tq, ncoset, maxord);
/*<          RC=RC*EL(1)/OLDLO >*/
	g_rc = g_rc * el[0] / oldlo;
/*<          OLDLO=EL(1) >*/
	oldlo = el[0];
/*<          CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN) >*/
	errors_(n, tq, &g_edn, &g_e, &g_eup, &g_bnd, &g_eddn);
/*<       ELSE >*/
    } else {
/*<          NEWQ = NQ >*/
	newq = g_nq;
/*<          RH = ONE/ (PRFAIL*DBLE(FLOAT(-KFAIL))) >*/
	g_rh = one / (prfail * (doublereal) (-g_kfail));
/*<       ENDIF >*/
    }
/*<       GO TO 40 >*/
    goto L40;
/* ********************************************************************** */
/*     CONTROL REACHES THIS STAGE IF 3 OR MORE FAILURES HAVE OCCURED. */
/*     IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE Y */
/*     ARRAY HAVE ERRORS OF THE WRONG ORDER. HENCE THE FIRST DERIVATIVE */
/*     IS RE-COMPUTED, AND THE ORDER IS SET TO 1. THEN H IS REDUCED BY A */
/*     FACTOR OF 10, AND THE STEP IS RETRIED. AFTER A TOTAL OF 7 */
/*     FAILURES AN EXIT IS TAKEN WITH KFLAG=-2. */
/* ********************************************************************** */
/*<  340  IF (KFAIL.EQ.-7) THEN >*/
L340:
    if (g_kfail == -7) {
/*     ERROR SMALLER THAN CAN BE HANDLED FOR PROBLEM */
/*<          KFLAG = -2 >*/
	*kflag = -2;
/*<          HOLD = H >*/
	g_hold = *h__;
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*     ********************************* */
/*     START FROM ORDER 1 AGAIN    * */
/*     ********************************* */
/*<       JCHANG = 1 >*/
    g_jchang = 1;
/*<       RH = DMAX1(HMIN/DABS(H),0.1D+0) >*/
/* Computing MAX */
    d__1 = *hmin / abs(*h__);
    g_rh = max(d__1,.1);
/*<       CALL HCHOSE(RH,H,OVRIDE) >*/
    hchose_(&g_rh, h__, &ovride);
/*<       H = H*RH >*/
    *h__ *= g_rh;
/*<       DO 350 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y(I,1) = YHOLD(I,1) >*/
	y[i__ + y_dim1] = yhold[i__ + yhold_dim1];
/*<          Y(I,2) = YHOLD(I,2) >*/
	y[i__ + (y_dim1 << 1)] = yhold[i__ + (yhold_dim1 << 1)];
/*<  350  CONTINUE >*/
/* L350: */
    }
/*<       IWEVAL = MITER >*/
    g_iweval = g_miter;
/*<       CFAIL = .TRUE. >*/
    g_cfail = TRUE_;
/*     SINCE WE HAVE HAD PROBLEMS PROCEED WITH THIS ORDER */
/*     FOR 10 STEPS (IF WE CAN) */
/*<       IDOUB = 10 >*/
    g_idoub = 10;
/*<       IF (NQ.EQ.1) GO TO 60 >*/
    if (g_nq == 1) {
	goto L60;
    }
/*<       NQ = 1 >*/
    g_nq = 1;
/*<       L = 2 >*/
    g_l = 2;
/*     RESET ORDER, RECALCULATE ERROR BOUNDS */
/*<       CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD) >*/
    coset_(&g_nq, el, elst, tq, ncoset, maxord);
/*<       LMAX = MAXDER + 1 >*/
    g_lmax = *maxder + 1;
/*<       RC = RC*EL(1)/OLDLO >*/
    g_rc =g_rc * el[0] / oldlo;
/*<       OLDLO = EL(1) >*/
    oldlo = el[0];
/*<       CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN) >*/
    errors_(n, tq, &g_edn, &g_e, &g_eup, &g_bnd, &g_eddn);
/*     NOW JUMP TO NORMAL CONTINUATION POINT */
/*<       GO TO 60 >*/
    goto L60;
/* ********************************************************************** */
/*     THE ITERATION FOR THE CORRECTED SOLUTION HAS CONVERGED. */
/*     UPDATE THE Y ARRAY. */
/* ********************************************************************** */
/*<  360  CONTINUE >*/
L360:
/*   **** */
/*   **** AMMEND **** */
/*   **** CHANGE 1,N BELOW TO 1,NVARS */
/*   **** */
/*<       DEMB=0.0D+0 >*/
    demb = 0.;
/*<       DO 361 I=1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          DEMB=DEMB+((Y(I,1)-YNHOLD(I,1))/SCALE(I))**2 >*/
/* Computing 2nd power */
	d__1 = (y[i__ + y_dim1] - ynhold[i__ + ynhold_dim1]) / scale[i__];
	demb += d__1 * d__1;
/*<  361  CONTINUE >*/
/* L361: */
    }
/*<       IF(DEMB.GT.4.0D+0*DBLE(FLOAT(N))) THEN >*/
    if (demb > (doublereal) (*n) * 4.) {
/*<          IEMB=1 >*/
	g_iemb = 1;
/*<          IJUS=1 >*/
	g_ijus = 1;
/*<          RED=0.5D+0 >*/
	red = .5;
/* ***  failed because of embedded error estimate */
/*<          NFAIL = NFAIL + 1 >*/
	++(*nfail);
/*<          GOTO 450 >*/
	goto L450;
/*<       ENDIF >*/
    }
/*<       DO 380 J2 = 2,LL >*/
    i__1 = ll;
    for (j2 = 2; j2 <= i__1; ++j2) {
/*<          J2M1=J2-1 >*/
	j2m1 = j2 - 1;
/*<          DO 370 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             Y(I,J2) = Y(I,J2M1) - YHOLD(I,J2M1) >*/
	    y[i__ + j2 * y_dim1] = y[i__ + j2m1 * y_dim1] - yhold[i__ + j2m1 *
		     yhold_dim1];
/*<  370     CONTINUE >*/
/* L370: */
	}
/*<  380  CONTINUE >*/
/* L380: */
    }
/*<       do 385 i=1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          yprime(i)=(y(i,1)-arh(i))/qqq >*/
	yprime[i__] = (y[i__ + y_dim1] - arh[i__]) / g_qqq;
/*<  385  continue >*/
/* L385: */
    }
/* --------------------------------------------------------------------- */
/*     IF THE COUNTER IDOUB EQUALS 2 AND WE ARE NOT ALREADY USING THE */
/*     MAXIMUM ALLOWABLE ORDER , STORE Y(I,LMAX+4) WHICH IS USED IN */
/*     ASSESSING THE POSSIBILITY OF INCREASING THE ORDER. IF IDOUB = 0 */
/*     CONTROL PASSES TO 480 WHERE AN ATTEMPT TO CHANGE THE STEPSIZE AND */
/*     ORDER IS MADE. */
/* ---------------------------------------------------------------------- */
/*<       IF (IDOUB.EQ.2.AND.L.NE.LMAX) THEN >*/
    if (g_idoub == 2 && g_l != g_lmax) {
/*<          LMP4=LMAX+4 >*/
	g_lmp4 = g_lmax + 4;
/*<          DO 390 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             Y(I,LMP4) = Y(I,LL) >*/
	    y[i__ + g_lmp4 * y_dim1] = y[i__ + ll * y_dim1];
/*<  390     CONTINUE >*/
/* L390: */
	}
/*<       END IF >*/
    }
/*<       IDOUB = IDOUB - 1 >*/
    --g_idoub;
/*<       TRANGE=(TEND-TOLD-H)*H >*/
    trange = (*tend - told - *h__) * *h__;
/*<       IF(TRANGE.LT.0.0D+0) THEN >*/
    if (trange < 0.) {
/*<          IDOUB = IDOUB + 2 >*/
	g_idoub += 2;
/*<          GOTO 440 >*/
	goto L440;
/*<       ENDIF >*/
    }
/*<       JCHANG = 0 >*/
    g_jchang = 0;
/*<       IF (IDOUB.EQ.0) THEN >*/
    if (g_idoub == 0) {
/*<          SAMPLE = .FALSE. >*/
	g_sample = FALSE_;
/*<          ISAMP = ISAMP + 1 >*/
	++g_isamp;
/*<          IF (ISAMP.EQ.4) THEN >*/
	if (g_isamp == 4) {
/*<             SAMPLE = .TRUE. >*/
	    g_sample = TRUE_;
/*<             ISAMP = 0 >*/
	    g_isamp = 0;
/*<          END IF >*/
	}
/* ********************************************************************** */
/*        NOW COMPUTE THE FACTORS PR1, PR2 AND PR3, BY WHICH */
/*        H COULD BE DIVIDED AT ORDER NQ-1, ORDER NQ AND ORDER NQ+1 */
/*        RESPECTIVELY. THE SMALLEST OF THESE IS DETERMINED AND THE NEW */
/*        ORDER CHOSEN ACCORDINGLY. IF THE ORDER IS TO BE INCREASED WE */
/*        MUST COMPUTE ONE MORE BACKWARD DIFFERENCE. */
/* ********************************************************************** */
/*<          PR3 = 1.D+20 >*/
	pr3 = 1e20;
/*<          FAC = 1.5D+0 >*/
	fac = 1.5;
/*<          IF(IEMB.EQ.1) FAC = 1.8D+0 >*/
	if (g_iemb == 1) {
	    fac = 1.8;
	}
/*<          DO 400 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             AYI = DABS(Y(I,1)) >*/
	    ayi = (d__1 = y[i__ + y_dim1], abs(d__1));
/*<             IF(ITOL.EQ.1) THEN >*/
	    if (*itol == 1) {
/*<                VHOLD = YMAX(I) >*/
		g_vhold = ymax[i__];
/*<             ELSE IF(ITOL.EQ.2) THEN >*/
	    } else if (*itol == 2) {
/*<                VHOLD = RTOL(1)*AYI + ATOL(1) >*/
		g_vhold = rtol[1] * ayi + atol[1];
/*<             ELSE IF(ITOL.EQ.3) THEN >*/
	    } else if (*itol == 3) {
/*<                VHOLD = RTOL(1)*AYI + ATOL(I) >*/
		g_vhold = rtol[1] * ayi + atol[i__];
/*<             ELSE IF(ITOL.EQ.4) THEN >*/
	    } else if (*itol == 4) {
/*<                VHOLD = RTOL(I)*AYI + ATOL(1) >*/
		g_vhold = rtol[i__] * ayi + atol[1];
/*<             ELSE IF(ITOL.EQ.5) THEN >*/
	    } else if (*itol == 5) {
/*<                VHOLD = RTOL(I)*AYI + ATOL(I) >*/
		g_vhold = rtol[i__] * ayi + atol[i__];
/*<             ENDIF >*/
	    }
/*<             SCALE(I)=VHOLD >*/
	    scale[i__] = g_vhold;
/*<  400     CONTINUE >*/
/* L400: */
	}
/*<          IF(NIND2.NE.0) THEN >*/
	if (*nind2 != 0) {
/*<             DO 4461 I=NIND1+1,NIND1+NIND2 >*/
	    i__1 = *nind1 + *nind2;
	    for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
/*<                SCALE(I) = SCALE(I)/HUSED >*/
		scale[i__] /= *hused;
/*<  4461       CONTINUE >*/
/* L4461: */
	    }
/*<          ENDIF >*/
	}
/*<          IF(NIND3.NE.0) THEN >*/
	if (*nind3 != 0) {
/*<             DO 4462 I=NIND1+NIND2+1,N >*/
	    i__1 = *n;
	    for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
/*<                SCALE(I) = SCALE(I)/(HUSED**2) >*/
/* Computing 2nd power */
		d__1 = *hused;
		scale[i__] /= d__1 * d__1;
/*<  4462       CONTINUE >*/
/* L4462: */
	    }
/*<          ENDIF >*/
	}
/*<          IF(L.NE.LMAX) THEN >*/
	if (g_l != g_lmax) {
/*<             LMP4 = LMAX + 4 >*/
	    g_lmp4 = g_lmax + 4;
/*<             DUP = ZERO >*/
	    g_dup = zero;
/*<             DO 401 I=1,N >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                DUP = DUP + ((Y(I,LL)-Y(I,LMP4))/SCALE(I))**2 >*/
/* Computing 2nd power */
		d__1 = (y[i__ + ll * y_dim1] - y[i__ + g_lmp4 * y_dim1]) / 
			scale[i__];
		g_dup += d__1 * d__1;
/*<  401        CONTINUE >*/
/* L401: */
	    }
/*<             IF(ITOL .EQ. 1) DUP = DUP/(RTOL(1)**2) >*/
	    if (*itol == 1) {
/* Computing 2nd power */
		d__1 = rtol[1];
		g_dup /= d__1 * d__1;
	    }
/*<             ENQ3 = 0.5D+0/DBLE(FLOAT(L+1)) >*/
	    enq3 = .5 / (doublereal) (g_l + 1);
/*<             PR3 = ((DUP/EUP)**ENQ3)*(FAC+0.2D+0) + 1.8D-6 >*/
	    d__1 = g_dup / g_eup;
	    pr3 = pow_dd(&d__1, &enq3) * (fac + .2) + 1.8e-6;
/*<          END IF >*/
	}
/*<          ENQ2 = 0.5D+0/DBLE(FLOAT(L)) >*/
	enq2 = .5 / (doublereal) g_l;
/*<          D = ZERO >*/
	d__ = zero;
/*<          DDOWN=ZERO >*/
	ddown = zero;
/*<          DO 410 I = 1,N >*/
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             D = D + (Y(I,LL)/SCALE(I))**2 >*/
/* Computing 2nd power */
	    d__1 = y[i__ + ll * y_dim1] / scale[i__];
	    d__ += d__1 * d__1;
/*<             DDOWN = DDOWN + (Y(I,LMP4)/SCALE(I))**2 >*/
/* Computing 2nd power */
	    d__1 = y[i__ + g_lmp4 * y_dim1] / scale[i__];
	    ddown += d__1 * d__1;
/*<  410     CONTINUE >*/
/* L410: */
	}
/*<          IF(ITOL .EQ.1) D = D/(RTOL(1)**2) >*/
	if (*itol == 1) {
/* Computing 2nd power */
	    d__1 = rtol[1];
	    d__ /= d__1 * d__1;
	}
/*<          PR2 = ((D/E)**ENQ2)*FAC + 1.6D-6 >*/
	d__1 = d__ / g_e;
	pr2 = pow_dd(&d__1, &enq2) * fac + 1.6e-6;
/*<          PR1 = 1.D+20 >*/
	pr1 = 1e20;
/*<          IF (NQ.GT.1) THEN >*/
	if (g_nq > 1) {
/*<             DDDOWN=ZERO >*/
	    dddown = zero;
/*<             DDOWN = ZERO >*/
	    ddown = zero;
/*<             DO 57420 I = 1,N >*/
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                DDOWN = DDOWN + (Y(I,L)/SCALE(I))**2 >*/
/* Computing 2nd power */
		d__1 = y[i__ + g_l * y_dim1] / scale[i__];
		ddown += d__1 * d__1;
/*<                DDDOWN = DDDOWN + (Y(I,L-1)/SCALE(I))**2 >*/
/* Computing 2nd power */
		d__1 = y[i__ + (g_l - 1) * y_dim1] / scale[i__];
		dddown += d__1 * d__1;
/*< 57420       CONTINUE >*/
/* L57420: */
	    }
/*<             IF(ITOL .EQ. 1) DDOWN = DDOWN/(RTOL(1)**2) >*/
	    if (*itol == 1) {
/* Computing 2nd power */
		d__1 = rtol[1];
		ddown /= d__1 * d__1;
	    }
/*<             ENQ1 = 0.5D+0/DBLE(FLOAT(NQ)) >*/
	    enq1 = .5 / (doublereal) g_nq;
/*<             PR1 = ((DDOWN/EDN)**ENQ1)*(FAC+0.1D+0) + 1.7D-6 >*/
	    d__1 = ddown / g_edn;
	    pr1 = pow_dd(&d__1, &enq1) * (fac + .1) + 1.7e-6;
/*<             IF(NQ.GT.2) THEN >*/
	    if (g_nq > 2) {
/*<                ENQ0 = 0.5D+0/DBLE(FLOAT(NQ-1)) >*/
		enq0 = .5 / (doublereal) (g_nq - 1);
/*<                PR0 = ((DDDOWN/EDDN)**ENQ0)*(FAC+0.1D+0) + 1.7D-6 >*/
		d__1 = dddown / g_eddn;
		pr0 = pow_dd(&d__1, &enq0) * (fac + .1) + 1.7e-6;
/*<                IF(PR0.GT.PR1) PR1 = PR0 >*/
		if (pr0 > pr1) {
		    pr1 = pr0;
		}
/*<                IF(DDDOWN.LT.DDOWN) DDOWN = DDDOWN >*/
		if (dddown < ddown) {
		    ddown = dddown;
		}
/*<             ENDIF >*/
	    }
/*<          END IF >*/
	}
/*<          IF(L.EQ.LMAX) DUP = 0.0D+0 >*/
	if (g_l == g_lmax) {
	    g_dup = 0.;
	}
/*<          IF(NQ.LE.1) GOTO 6578 >*/
	if (g_nq <= 1) {
	    goto L6578;
	}
/*<          IF(DUP.GT.D.AND.D.GT.DDOWN) THEN >*/
	if (g_dup > d__ && d__ > ddown) {
/*<             PR2=1.0D+30 >*/
	    pr2 = 1e30;
/*<             PR3=1.0D+30 >*/
	    pr3 = 1e30;
/*<          ENDIF >*/
	}
/*<  6578    CONTINUE >*/
L6578:
/*<          IF (PR2.LE.PR3) THEN >*/
	if (pr2 <= pr3) {
/*<             IF (PR2.GT.PR1) THEN >*/
	    if (pr2 > pr1) {
/*<                NEWQ = NQ - 1 >*/
		newq = g_nq - 1;
/*<                RH = 1.0D+0/PR1 >*/
		g_rh = 1. / pr1;
/*<             ELSE >*/
	    } else {
/*<                NEWQ = NQ >*/
		newq = g_nq;
/*<                RH = 1.0D+0/PR2 >*/
		g_rh = 1. / pr2;
/*<             END IF >*/
	    }
/*<          ELSE IF (PR3.LT.PR1) THEN >*/
	} else if (pr3 < pr1) {
/*<             NEWQ = L >*/
	    newq = g_l;
/*<             RH = 1.0D+0/PR3 >*/
	    g_rh = 1. / pr3;
/*<          ELSE >*/
	} else {
/*<             NEWQ = NQ - 1 >*/
	    newq = g_nq - 1;
/*<             RH = 1.0D+0/PR1 >*/
	    g_rh = 1. / pr1;
/*<          END IF >*/
	}
/*<          IEMB=0 >*/
	g_iemb = 0;
/*<          IF(RH.GT.1.0D+0.AND.RH.LT.1.1D+0) THEN >*/
	if (g_rh > 1. && g_rh < 1.1) {
/*<             IDOUB=10 >*/
	    g_idoub = 10;
/*<             NQ=NQUSED >*/
	    g_nq = *nqused;
/*<             L=NQ+1 >*/
	    g_l = g_nq + 1;
/*<             GOTO 440 >*/
	    goto L440;
/*<          ENDIF >*/
	}
/*<          RH = DMIN1(RH,RMAX) >*/
	g_rh = min(g_rh, g_stiff_rmax);
/*<          CALL HCHOSE(RH,H,OVRIDE) >*/
	hchose_(&g_rh, h__, &ovride);
/*<          IF ((JSINUP.LE.20).AND.(KFLAG.EQ.0).AND.(RH.LT.1.1D+0)) THEN >*/
	if (g_jsinup <= 20 && *kflag == 0 && g_rh < 1.1) {
/*           WE HAVE RUN INTO PROBLEMS */
/*<             IDOUB = 10 >*/
	    g_idoub = 10;
/*<             NQ = NQUSED >*/
	    g_nq = *nqused;
/*<             L = NQ + 1 >*/
	    g_l = g_nq + 1;
/*<             GO TO 440 >*/
	    goto L440;
/*<          END IF >*/
	}
/* ********************************************************************** */
/*        IF THERE IS A CHANGE IN ORDER, RESET NQ, L AND THE */
/*        COEFFICIENTS. IN ANY CASE H IS RESET  AND THE */
/*        Y ARRAY IS RE-SCALED */
/* ********************************************************************** */
/*<          IF(IMAS.NE.0.and.NIND3.NE.0) THEN >*/
	if (g_imas != 0 && *nind3 != 0) {
/*<             IF(NQ.LE.2.AND.PR3.LT.1.0D+0) NEWQ=NQ+1 >*/
	    if (g_nq <= 2 && pr3 < 1.) {
		newq = g_nq + 1;
	    }
/*<          ENDIF >*/
	}
/*<          IF (NEWQ.NE.NQ) THEN >*/
	if (newq != g_nq) {
/*<             IF (NEWQ.GT.NQ) THEN >*/
	    if (newq > g_nq) {
/*              ADD AN EXTRA TERM TO THE HISTORY ARRAY */
/*<                DO 430 I = 1,N >*/
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/*<                   Y(I,LL) = Y(I,L) - YHOLD(I,L) >*/
		    y[i__ + ll * y_dim1] = y[i__ + g_l * y_dim1] - yhold[i__ + 
			    g_l * yhold_dim1];
/*<  430           CONTINUE >*/
/* L430: */
		}
/*<             END IF >*/
	    }
/*<             NQ = NEWQ >*/
	    g_nq = newq;
/*<             L = NQ + 1 >*/
	    g_l = g_nq + 1;
/*           RESET ORDER,RECALCULATE ERROR BOUNDS */
/*<             CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD) >*/
	    coset_(&g_nq, el, elst, tq, ncoset, maxord);
/*<             LMAX = MAXDER + 1 >*/
	    g_lmax = *maxder + 1;
/*<             RC = RC*EL(1)/OLDLO >*/
	    g_rc = g_rc * el[0] / oldlo;
/*<             OLDLO = EL(1) >*/
	    oldlo = el[0];
/*<             CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN) >*/
	    errors_(n, tq, &g_edn, &g_e, &g_eup, &g_bnd, &g_eddn);
/*<          END IF >*/
	}
/*<          RH = DMAX1(RH,HMIN/DABS(H)) >*/
/* Computing MAX */
	d__1 = g_rh, d__2 = *hmin / abs(*h__);
	g_rh = max(d__1,d__2);
/*<          RH = DMIN1(RH,HMAX/DABS(H),RMAX) >*/
/* Computing MIN */
	d__1 = g_rh, d__2 = *hmax / abs(*h__), d__1 = min(d__1,d__2);
	g_rh = min(d__1, g_stiff_rmax);
/*<          CALL RSCALE(N,L,RH,Y) >*/
	rscale_(n, &g_l, &g_rh, &y[y_offset]);
/*<          RMAX = 10.0D+0 >*/
	g_stiff_rmax = 10.;
/*<          JCHANG = 1 >*/
	g_jchang = 1;
/*<          H = H*RH >*/
	*h__ *= g_rh;
/*<          RC = RC*RH >*/
	g_rc *= g_rh;
/*<          IF(JSNOLD.GT.IBND) RC=ZERO >*/
	if (g_jsnold > g_ibnd) {
	    g_rc = zero;
	}
/*<          IDOUB = L + 1 >*/
	g_idoub = g_l + 1;
/*<       END IF >*/
    }
/*<  440  CONTINUE >*/
L440:
/* ---------------------------------------------------------------------- */
/*     STORE THE Y ARRAY IN THE MATRIX YHOLD.  STORE IN THE Y ARRAY THE */
/*     INFORMATION NECESSARY TO PERFORM AN INTERPOLATION TO FIND THE */
/*     SOLUTION AT THE SPECIFIED OUTPUT POINT IF APPROPRIATE. */
/* ---------------------------------------------------------------------- */
/*<       CALL CPYARY(N*L,Y,YHOLD) >*/
    i__1 = *n * g_l;
    cpyary_(&i__1, &y[y_offset], &yhold[yhold_offset]);
/*<       NSTEP = NSTEP + 1 >*/
    ++(*nstep);
/*<       JSINUP = JSINUP + 1 >*/
    ++g_jsinup;
/*<       JSNOLD = JSNOLD + 1 >*/
    ++g_jsnold;
/*<       JSTART = NQUSED >*/
    *jstart = *nqused;
/*<       T = TOLD + HUSED >*/
    *t = told + *hused;
/*<       HOLD = H >*/
    g_hold = *h__;
/*<       KFAIL = 0 >*/
    g_kfail = 0;
/*<       NEWPAR = 0 >*/
    newpar = 0;
/*<       CFAIL = .FALSE. >*/
    g_cfail = FALSE_;
/*<       RETURN >*/
    return 0;
/*<  450  CONTINUE >*/
L450:
/*<       FINISH = .FALSE. >*/
    g_finish = FALSE_;
/*<       T=TOLD >*/
    *t = told;
/*<       RMAX=2.0D+0 >*/
    g_stiff_rmax = 2.;
/*<       DO 460 J1=1,L >*/
    i__1 = g_l;
    for (j1 = 1; j1 <= i__1; ++j1) {
/*<          DO 460 I=1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             Y(I,J1)=YHOLD(I,J1) >*/
	    y[i__ + j1 * y_dim1] = yhold[i__ + j1 * yhold_dim1];
/*<  460  CONTINUE >*/
/* L460: */
	}
    }
/*<       IF(DABS(H).LE.HMIN*1.00001D+0) THEN >*/
    if (abs(*h__) <= *hmin * 1.00001) {

/*   CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED */

/*<          IF(NSTEP.EQ.0) THEN >*/
	if (*nstep == 0) {
/*<             KFLAG=-1 >*/
	    *kflag = -1;
/*<          ELSE >*/
	} else {
/*<             KFLAG=-3 >*/
	    *kflag = -3;
/*<          END IF >*/
	}

/*    TO SUPPRESS ERROR MESSAGES AT START AS H MAY */
/*    HAVE BEEN TOO LARGE ON THE FIRST STEP. */

/*<          HOLD=H >*/
	g_hold = *h__;
/*<          FINISH = .TRUE. >*/
	g_finish = TRUE_;
/*<       END IF >*/
    }
/*<       RH = RED >*/
    g_rh = red;
/*<       IREDO=1 >*/
    iredo = 1;

/*     TRY AGAIN WITH UPDATED PARTIALS */

/*<  8000 if (ierr .ne. 0) then >*/
L8000:
    if (*ierr != 0) {
/*<          write(*,1975) >*/
	s_wsfe(&io___245);
	e_wsfe();
/*<  1 >*/
/*<          h= h/2          >*/
	*h__ /= 2;
/*<          IF(H.LT.EPSJAC/100.0D+0) THEN >*/
	if (*h__ < *epsjac / 100.) {
/*<             WRITE(6,9161) >*/
	    s_wsfe(&io___246);
	    e_wsfe();
/*<  9161       FORMAT (/,/,'STEPSIZE IS TOO SMALL') >*/
/*<             IDID = -7 >*/
	    idid = -7;
/*<             RETURN >*/
	    return 0;
/*<          ENDIF >*/
	}
/*<          T= TOLD >*/
	*t = told;
/*<          IF ((T-TOUT)*H.GE.0.0D+0) THEN >*/
	if ((*t - *tout) * *h__ >= 0.) {
/*           HAVE OVERSHOT TOUT */
/*<             WRITE (LOUT,*) T,TOUT,H >*/
	    io___248.ciunit = *lout;
	    s_wsle(&io___248);
	    do_lio(&c__5, &c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
	    do_lio(&c__5, &c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal))
		    ;
	    do_lio(&c__5, &c__1, (char *)&(*h__), (ftnlen)sizeof(doublereal));
	    e_wsle();
/*<             CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0) >*/
	    interp_(n, jstart, h__, t, &y[y_offset], tout, &y0);
/*<             HO = H >*/
	    ho = *h__;
/*<             T0 = TOUT >*/
	    t0 = *tout;
/*<             IDID = -5 >*/
	    idid = -5;
/*<             RETURN >*/
	    return 0;
/*<          ENDIF >*/
	}
/*<          IERR = 0 >*/
	*ierr = 0;
/*<          jstart = -1 >*/
	*jstart = -1;
/*<          goto 30 >*/
	goto L30;
/*<       endif >*/
    }

/*<       IF(IJUS.EQ.0) CALL HCHOSE(RH,H,OVRIDE) >*/
    if (g_ijus == 0) {
	hchose_(&g_rh, h__, &ovride);
    }
/*<       IF(.NOT.FINISH) THEN >*/
    if (! g_finish) {
/*<          GO TO 40 >*/
	goto L40;
/*<       ELSE >*/
    } else {
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }
/*<  9000 FORMAT (1X,' CORRECTOR HAS NOT CONVERGED') >*/
/* L9000: */
/*<       END >*/
    return 0;
} /* stiff_ */

/* ------------------- END OF SUBROUTINE STIFF -------------------------- */
/*<       SUBROUTINE RSCALE(N,L,RH,Y) >*/
/* Subroutine */ int rscale_(integer *n, integer *l, doublereal *rh, 
	doublereal *y)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2, i__3;

    /* Local variables */
    /*static*/ integer i__, j, j1;
    /*static*/ doublereal di[64]	/* was [8][8] */, ta, tb, tc, td, te, tf, zz;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER L,N >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  Y(N,12) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I,J,J1 >*/
/*     .. */
/*     .. LOCAL ARRAYS .. */
/*<       DIMENSION  DI(8,8) >*/
/*     .. */
/*     .. DATA STATEMENTS .. */
/*     SUBROUTINE IS FOR RESCALING THE HISTORY ARRAY AFTER A CHANGE IN */
/*     STEPSIZE */

/*     N      ORDER OF THE PROBLEM */
/*     L      NUMBER OF TERMS IN THE HISTORY ARRAY TO BE RESCALED */
/*     RH     RATIO OF THE STEPSIZE CHANGE (I.E. RH = HNEW/HOLD) */
/*     Y()    THE HISTORY ARRAY */

/*<       DATA  ZERO/0.0D+0/ >*/
    /* Parameter adjustments */
    y_dim1 = *n;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;

    /* Function Body */
/*     .. */
/*<       DI(2,2) = RH >*/
    di[9] = *rh;
/*<       IF (L.GT.2) THEN >*/
    if (*l > 2) {
/*<          TA = RH*RH >*/
	ta = *rh * *rh;
/*<          DI(2,3) = RH* (1.0D+0-RH)/2.0D+0 >*/
	di[17] = *rh * (1. - *rh) / 2.;
/*<          DI(3,3) = TA >*/
	di[18] = ta;
/*<          IF (L.GT.3) THEN >*/
	if (*l > 3) {
/*<             TB = TA*RH >*/
	    tb = ta * *rh;
/*<             DI(2,4) = RH* ((RH-3.0D+0)*RH+2.0D+0)/6.0D+0 >*/
	    di[25] = *rh * ((*rh - 3.) * *rh + 2.) / 6.;
/*<             DI(3,4) = TA* (1.0D+0-RH) >*/
	    di[26] = ta * (1. - *rh);
/*<             DI(4,4) = TB >*/
	    di[27] = tb;
/*<             IF (L.GT.4) THEN >*/
	    if (*l > 4) {
/*<                TC = TB*RH >*/
		tc = tb * *rh;
/*<    >*/
		di[33] = -(((*rh - 6.) * *rh + 11.) * *rh - 6.) * *rh / 24.;
/*<                DI(3,5) = TA* ((7.0D+0*RH-18.0D+0)*RH+11.0D+0)/12.0D+0 >*/
		di[34] = ta * ((*rh * 7. - 18.) * *rh + 11.) / 12.;
/*<                DI(4,5) = 1.5D+0*TB* (1.0D+0-RH) >*/
		di[35] = tb * 1.5 * (1. - *rh);
/*<                DI(5,5) = TC >*/
		di[36] = tc;
/*<                IF (L.GT.5) THEN >*/
		if (*l > 5) {
/*<                   TD = TC*RH >*/
		    td = tc * *rh;
/*<    >*/
		    di[41] = ((((*rh - 10.) * *rh + 35.) * *rh - 50.) * *rh + 
			    24.) * *rh / 120.;
/*<    >*/
		    di[42] = -(((*rh * 3. - 14.) * *rh + 21.) * *rh - 10.) * 
			    ta / 12.;
/*<                   DI(4,6) = ((5.0D+0*RH-12.0D+0)*RH+7.0D+0)*TB/4.0D+0 >*/
		    di[43] = ((*rh * 5. - 12.) * *rh + 7.) * tb / 4.;
/*<                   DI(5,6) = 2.0D+0*TC* (1.0D+0-RH) >*/
		    di[44] = tc * 2. * (1. - *rh);
/*<                   DI(6,6) = TD >*/
		    di[45] = td;
/*<                   IF (L.GT.6) THEN >*/
		    if (*l > 6) {
/*<                      TE = TD*RH >*/
			te = td * *rh;
/*<    >*/
			di[49] = -(*rh) * (*rh - 1.) * (*rh - 2.) * (*rh - 3.)
				 * (*rh - 4.) * (*rh - 5.) / 720.;
/*<    >*/
			di[50] = ta * ((((*rh * 62. - 450.) * *rh + 1190.) * *
				rh - 1350.) * *rh + 548.) / 720.;
/*<    >*/
			di[51] = tb * (((*rh * -18. + 75.) * *rh - 102.) * *
				rh + 45.) / 24.;
/*<    >*/
			di[52] = tc * ((*rh * 13. - 30.) * *rh + 17.) / 6.;
/*<                      DI(6,7) = 2.5D+0*TD* (1.0D+0-RH) >*/
			di[53] = td * 2.5 * (1. - *rh);
/*<                      DI(7,7) = TE >*/
			di[54] = te;
/*<                      IF (L.GT.7) THEN >*/
			if (*l > 7) {
/*<                         TF = TE*RH >*/
			    tf = te * *rh;
/*<    >*/
			    di[57] = *rh * (*rh - 1.) * (*rh - 2.) * (*rh - 
				    3.) * (*rh - 4.) * (*rh - 5.) * (*rh - 6.)
				     / 5040.;
/*<    >*/
			    di[58] = ta * (((((*rh * -126. + 1302.) * *rh - 
				    5250.) * *rh + 10290.) * *rh - 9744.) * *
				    rh + 3528.) / 5040.;
/*<    >*/
			    di[59] = tb * ((((*rh * 43. - 270.) * *rh + 625.) 
				    * *rh - 630.) * *rh + 232.) / 120.;
/*<    >*/
			    di[60] = tc * (((*rh * -10. + 39.) * *rh - 50.) * 
				    *rh + 21.) / 6.;
/*<    >*/
			    di[61] = td * ((*rh * 20. - 45.) * *rh + 25.) / 
				    6.;
/*<                         DI(7,8) = 3.0D+0*TE* (1.0D+0-RH) >*/
			    di[62] = te * 3. * (1. - *rh);
/*<                         DI(8,8) = TF >*/
			    di[63] = tf;
/*<                      END IF >*/
			}
/*<                   END IF >*/
		    }
/*<                END IF >*/
		}
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       DO 30 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DO 20 J = 2,L >*/
	i__2 = *l;
	for (j = 2; j <= i__2; ++j) {
/*<           ZZ = ZERO >*/
	    zz = zero;
/*<           DO 10 J1 = J,L >*/
	    i__3 = *l;
	    for (j1 = j; j1 <= i__3; ++j1) {
/*<             ZZ = ZZ + DI(J,J1)*Y(I,J1) >*/
		zz += di[j + (j1 << 3) - 9] * y[i__ + j1 * y_dim1];
/*<    10     CONTINUE >*/
/* L10: */
	    }
/*<           Y(I,J) = ZZ >*/
	    y[i__ + j * y_dim1] = zz;
/*<    20   CONTINUE >*/
/* L20: */
	}
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* rscale_ */

/* --------------------------------------------------------------------------- */
/*<       SUBROUTINE CPYARY(NELEM,SOURCE,TARGET) >*/
/* Subroutine */ int cpyary_(integer *nelem, doublereal *source, doublereal *
	target)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    /*static*/ integer i__;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/

/*     COPIES THE ARRAY SOURCE() INTO THE ARRAY TARGET() */

/*     THIS SUBROUTINE COULD BE REPLACED BY THE BLAS ROUTINE SCOPY */
/*     (AFTER CHANGING THE ARGUMENT LIST APPROPRIATELY) */

/*     .. SCALAR ARGUMENTS .. */
/*<       INTEGER NELEM >*/
/*     .. */
/*     .. ARRAY ARGUMENTS .. */
/*<       DIMENSION  SOURCE(NELEM),TARGET(NELEM) >*/
/*     .. */
/*     .. LOCAL SCALARS .. */
/*<       INTEGER I >*/
/*     .. */
/*<       DO 10 I = 1,NELEM >*/
    /* Parameter adjustments */
    --target;
    --source;

    /* Function Body */
    i__1 = *nelem;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          TARGET(I) = SOURCE(I) >*/
	target[i__] = source[i__];
/*<  10   CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* cpyary_ */

/* ---------------------------------------------------------------------------- */
/*<       SUBROUTINE HCHOSE(RH,H,OVRIDE) >*/
/* Subroutine */ int hchose_(doublereal *rh, doublereal *h__, logical *ovride)
{
    /*static*/ integer i__, i2;

/*<       IMPLICIT DOUBLE PRECISION(A-H,O-Z) >*/
/*<       COMMON / STPSZE / HSTPSZ(2,14) >*/
/*<       LOGICAL OVRIDE >*/

/*     FIRST MOVE ALL ELEMENTS DOWN ONE PLACE */

/*<       IF (H.NE.HSTPSZ(2,1)) THEN >*/
    if (*h__ != stpsze_1.hstpsz[1]) {
/*<          DO 5 I=12,2,-1 >*/
	for (i__ = 12; i__ >= 2; --i__) {
/*<             I2=I-1 >*/
	    i2 = i__ - 1;
/*<             HSTPSZ(1,I)=HSTPSZ(1,I2) >*/
	    stpsze_1.hstpsz[(i__ << 1) - 2] = stpsze_1.hstpsz[(i2 << 1) - 2];
/*<  5       HSTPSZ(2,I)=HSTPSZ(2,I2) >*/
/* L5: */
	    stpsze_1.hstpsz[(i__ << 1) - 1] = stpsze_1.hstpsz[(i2 << 1) - 1];
	}

/*          NOW INSERT VALUE OF H USED BEFORE THIS CALL */

/*<          HSTPSZ(1,2)=H/HSTPSZ(2,1) >*/
	stpsze_1.hstpsz[2] = *h__ / stpsze_1.hstpsz[1];
/*<          HSTPSZ(2,1)=H >*/
	stpsze_1.hstpsz[1] = *h__;
/*<       END IF >*/
    }

/*     NOW DECIDE ON THE NEW CHANGE */

/*<       IF (RH.GT.1.0D+0) THEN >*/
    if (*rh > 1.) {
/*<          OVRIDE=.FALSE. >*/
	*ovride = FALSE_;
/*<       ELSE IF (HSTPSZ(1,2).LE.1.0D+0) THEN >*/
    } else if (stpsze_1.hstpsz[2] <= 1.) {
/*<          OVRIDE=.FALSE. >*/
	*ovride = FALSE_;
/*<       ELSE IF ((RH*H).LE.HSTPSZ(2,2)) THEN >*/
    } else if (*rh * *h__ <= stpsze_1.hstpsz[3]) {
/*<          OVRIDE=.FALSE. >*/
	*ovride = FALSE_;
/*<       ELSE >*/
    } else {
/*<          RH=HSTPSZ(2,2)/H >*/
	*rh = stpsze_1.hstpsz[3] / *h__;
/*<          OVRIDE=.TRUE. >*/
	*ovride = TRUE_;
/*<       END IF >*/
    }
/*<       HSTPSZ(1,1)=RH >*/
    stpsze_1.hstpsz[0] = *rh;
/*<       RETURN >*/
    return 0;
/*<       END       >*/
} /* hchose_ */


/*  ************************************************************ */

/*<       DOUBLE PRECISION FUNCTION DLAMCH( CMACH ) >*/
doublereal dlamch_(char *cmach, ftnlen cmach_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal base;
    static integer beta;
    static doublereal emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static doublereal rmin, rmax, t, rmach;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small, sfmin;
    extern /* Subroutine */ int dlamc2_(integer *, integer *, logical *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *);
    static integer it;
    static doublereal rnd, eps;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER          CMACH >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMCH determines double precision machine parameters. */

/*  Arguments */
/*  ========= */

/*  CMACH   (input) CHARACTER*1 */
/*          Specifies the value to be returned by DLAMCH: */
/*          = 'E' or 'e',   DLAMCH := eps */
/*          = 'S' or 's ,   DLAMCH := sfmin */
/*          = 'B' or 'b',   DLAMCH := base */
/*          = 'P' or 'p',   DLAMCH := eps*base */
/*          = 'N' or 'n',   DLAMCH := t */
/*          = 'R' or 'r',   DLAMCH := rnd */
/*          = 'M' or 'm',   DLAMCH := emin */
/*          = 'U' or 'u',   DLAMCH := rmin */
/*          = 'L' or 'l',   DLAMCH := emax */
/*          = 'O' or 'o',   DLAMCH := rmax */

/*          where */

/*          eps   = relative machine precision */
/*          sfmin = safe minimum, such that 1/sfmin does not overflow */
/*          base  = base of the machine */
/*          prec  = eps*base */
/*          t     = number of (base) digits in the mantissa */
/*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/*          emin  = minimum exponent before (gradual) underflow */
/*          rmin  = underflow threshold - base**(emin-1) */
/*          emax  = largest exponent before overflow */
/*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

/* ===================================================================== */

/*     .. Parameters .. */
/*<       DOUBLE PRECISION   ONE, ZERO >*/
/*<       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 ) >*/
/*     .. */
/*     .. Local Scalars .. */
/*<       LOGICAL            FIRST, LRND >*/
/*<       INTEGER            BETA, IMAX, IMIN, IT >*/
/*<    >*/
/*     .. */
/*     .. External Functions .. */
/*<       LOGICAL            LSAME >*/
/*<       EXTERNAL           LSAME >*/
/*     .. */
/*     .. External Subroutines .. */
/*<       EXTERNAL           DLAMC2 >*/
/*     .. */
/*     .. Save statement .. */
/*<    >*/
/*     .. */
/*     .. Data statements .. */
/*<       DATA               FIRST / .TRUE. / >*/
/*     .. */
/*     .. Executable Statements .. */

/*<       IF( FIRST ) THEN >*/
    if (first) {
/*<          FIRST = .FALSE. >*/
	first = FALSE_;
/*<          CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX ) >*/
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
/*<          BASE = BETA >*/
	base = (doublereal) beta;
/*<          T = IT >*/
	t = (doublereal) it;
/*<          IF( LRND ) THEN >*/
	if (lrnd) {
/*<             RND = ONE >*/
	    rnd = 1.;
/*<             EPS = ( BASE**( 1-IT ) ) / 2 >*/
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1) / 2;
/*<          ELSE >*/
	} else {
/*<             RND = ZERO >*/
	    rnd = 0.;
/*<             EPS = BASE**( 1-IT ) >*/
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1);
/*<          END IF >*/
	}
/*<          PREC = EPS*BASE >*/
	prec = eps * base;
/*<          EMIN = IMIN >*/
	emin = (doublereal) imin;
/*<          EMAX = IMAX >*/
	emax = (doublereal) imax;
/*<          SFMIN = RMIN >*/
	sfmin = rmin;
/*<          SMALL = ONE / RMAX >*/
	small = 1. / rmax;
/*<          IF( SMALL.GE.SFMIN ) THEN >*/
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rounding */
/*           causing overflow when computing  1/sfmin. */

/*<             SFMIN = SMALL*( ONE+EPS ) >*/
	    sfmin = small * (eps + 1.);
/*<          END IF >*/
	}
/*<       END IF >*/
    }

/*<       IF( LSAME( CMACH, 'E' ) ) THEN >*/
    if (lsame_(cmach, "E", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = EPS >*/
	rmach = eps;
/*<       ELSE IF( LSAME( CMACH, 'S' ) ) THEN >*/
    } else if (lsame_(cmach, "S", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = SFMIN >*/
	rmach = sfmin;
/*<       ELSE IF( LSAME( CMACH, 'B' ) ) THEN >*/
    } else if (lsame_(cmach, "B", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = BASE >*/
	rmach = base;
/*<       ELSE IF( LSAME( CMACH, 'P' ) ) THEN >*/
    } else if (lsame_(cmach, "P", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = PREC >*/
	rmach = prec;
/*<       ELSE IF( LSAME( CMACH, 'N' ) ) THEN >*/
    } else if (lsame_(cmach, "N", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = T >*/
	rmach = t;
/*<       ELSE IF( LSAME( CMACH, 'R' ) ) THEN >*/
    } else if (lsame_(cmach, "R", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = RND >*/
	rmach = rnd;
/*<       ELSE IF( LSAME( CMACH, 'M' ) ) THEN >*/
    } else if (lsame_(cmach, "M", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = EMIN >*/
	rmach = emin;
/*<       ELSE IF( LSAME( CMACH, 'U' ) ) THEN >*/
    } else if (lsame_(cmach, "U", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = RMIN >*/
	rmach = rmin;
/*<       ELSE IF( LSAME( CMACH, 'L' ) ) THEN >*/
    } else if (lsame_(cmach, "L", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = EMAX >*/
	rmach = emax;
/*<       ELSE IF( LSAME( CMACH, 'O' ) ) THEN >*/
    } else if (lsame_(cmach, "O", (ftnlen)1, (ftnlen)1)) {
/*<          RMACH = RMAX >*/
	rmach = rmax;
/*<       END IF >*/
    }

/*<       DLAMCH = RMACH >*/
    ret_val = rmach;
/*<       RETURN >*/
    return ret_val;

/*     End of DLAMCH */

/*<       END >*/
} /* dlamch_ */









/* *********************************************************************** */

/*<       SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 ) >*/
/* Subroutine */ int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */

	/* grn: these MUST be static because
	  they won't change on new calls because the machine
	  won't change between calls.
	  The boolean "first" determines if everything is 
	  initialized properly */
    static doublereal a, b, c__, f;
    static doublereal savec;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static doublereal t1, t2;
    static doublereal one, qtr;

	/* placed back from global to local static */
	static integer g_dlamc1_lbeta;
	static logical g_dlamc1_lrnd;
	static integer g_dlamch1_lt;
	static logical g_lieee1;

/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*<       LOGICAL            IEEE1, RND >*/
/*<       INTEGER            BETA, T >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC1 determines the machine parameters given by BETA, T, RND, and */
/*  IEEE1. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) INTEGER */
/*          The base of the machine. */

/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */

/*  IEEE1   (output) LOGICAL */
/*          Specifies whether rounding appears to be done in the IEEE */
/*          'round to nearest' style. */

/*  Further Details */
/*  =============== */

/*  The routine is based on the routine  ENVRON  by Malcolm and */
/*  incorporates suggestions by Gentleman and Marovich. See */

/*     Malcolm M. A. (1972) Algorithms to reveal properties of */
/*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

/*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
/*        that reveal properties of floating point arithmetic units. */
/*        Comms. of the ACM, 17, 276-277. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*<       LOGICAL            FIRST, LIEEE1, LRND >*/
/*<       INTEGER            LBETA, LT >*/
/*<       DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2 >*/
/*     .. */
/*     .. External Functions .. */
/*<       DOUBLE PRECISION   DLAMC3 >*/
/*<       EXTERNAL           DLAMC3 >*/
/*     .. */
/*     .. Save statement .. */
/*<       SAVE               FIRST, LIEEE1, LBETA, LRND, LT >*/
/*     .. */
/*     .. Data statements .. */
/*<       DATA               FIRST / .TRUE. / >*/
/*     .. */
/*     .. Executable Statements .. */

/*<       IF( FIRST ) THEN >*/
    if (first) {
/*<          FIRST = .FALSE. >*/
	first = FALSE_;
/*<          ONE = 1 >*/
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
/*        IEEE1, T and RND. */

/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are  stored and not held in registers,  or */
/*        are not affected by optimizers. */

/*        Compute  a = 2.0**m  with the  smallest positive integer m such */
/*        that */

/*           fl( a + 1.0 ) = a. */

/*<          A = 1 >*/
	a = 1.;
/*<          C = 1 >*/
	c__ = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
/*<    10    CONTINUE >*/
L10:
/*<          IF( C.EQ.ONE ) THEN >*/
	if (c__ == one) {
/*<             A = 2*A >*/
	    a *= 2;
/*<             C = DLAMC3( A, ONE ) >*/
	    c__ = dlamc3_(&a, &one);
/*<             C = DLAMC3( C, -A ) >*/
	    d__1 = -a;
	    c__ = dlamc3_(&c__, &d__1);
/*<             GO TO 10 >*/
	    goto L10;
/*<          END IF >*/
	}
/* +       END WHILE */

/*        Now compute  b = 2.0**m  with the smallest positive integer m */
/*        such that */

/*           fl( a + b ) .gt. a. */

/*<          B = 1 >*/
	b = 1.;
/*<          C = DLAMC3( A, B ) >*/
	c__ = dlamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
/*<    20    CONTINUE >*/
L20:
/*<          IF( C.EQ.A ) THEN >*/
	if (c__ == a) {
/*<             B = 2*B >*/
	    b *= 2;
/*<             C = DLAMC3( A, B ) >*/
	    c__ = dlamc3_(&a, &b);
/*<             GO TO 20 >*/
	    goto L20;
/*<          END IF >*/
	}
/* +       END WHILE */

/*        Now compute the base.  a and c  are neighbouring floating point */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
/*        their difference is beta. Adding 0.25 to c is to ensure that it */
/*        is truncated to beta and not ( beta - 1 ). */

/*<          QTR = ONE / 4 >*/
	qtr = one / 4;
/*<          SAVEC = C >*/
	savec = c__;
/*<          C = DLAMC3( C, -A ) >*/
	d__1 = -a;
	c__ = dlamc3_(&c__, &d__1);
/*<          LBETA = C + QTR >*/
	g_dlamc1_lbeta = (integer) (c__ + qtr);

/*        Now determine whether rounding or chopping occurs,  by adding a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

/*<          B = LBETA >*/
	b = (doublereal) g_dlamc1_lbeta;
/*<          F = DLAMC3( B / 2, -B / 100 ) >*/
	d__1 = b / 2;
	d__2 = -b / 100;
	f = dlamc3_(&d__1, &d__2);
/*<          C = DLAMC3( F, A ) >*/
	c__ = dlamc3_(&f, &a);
/*<          IF( C.EQ.A ) THEN >*/
	if (c__ == a) {
/*<             LRND = .TRUE. >*/
	    g_dlamc1_lrnd = TRUE_;
/*<          ELSE >*/
	} else {
/*<             LRND = .FALSE. >*/
	    g_dlamc1_lrnd = FALSE_;
/*<          END IF >*/
	}
/*<          F = DLAMC3( B / 2, B / 100 ) >*/
	d__1 = b / 2;
	d__2 = b / 100;
	f = dlamc3_(&d__1, &d__2);
/*<          C = DLAMC3( F, A ) >*/
	c__ = dlamc3_(&f, &a);
/*<    >*/
	if (g_dlamc1_lrnd && c__ == a) {
	    g_dlamc1_lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round to */
/*        nearest' style. B/2 is half a unit in the last place of the two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
/*        A, but adding B/2 to SAVEC should change SAVEC. */

/*<          T1 = DLAMC3( B / 2, A ) >*/
	d__1 = b / 2;
	t1 = dlamc3_(&d__1, &a);
/*<          T2 = DLAMC3( B / 2, SAVEC ) >*/
	d__1 = b / 2;
	t2 = dlamc3_(&d__1, &savec);
/*<          LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND >*/
	g_lieee1 = t1 == a && t2 > savec && g_dlamc1_lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part of */
/*        log to the base beta of a,  however it is safer to determine  t */
/*        by powering.  So we find t as the smallest positive integer for */
/*        which */

/*           fl( beta**t + 1.0 ) = 1.0. */

/*<          LT = 0 >*/
	g_dlamch1_lt = 0;
/*<          A = 1 >*/
	a = 1.;
/*<          C = 1 >*/
	c__ = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
/*<    30    CONTINUE >*/
L30:
/*<          IF( C.EQ.ONE ) THEN >*/
	if (c__ == one) {
/*<             LT = LT + 1 >*/
	    ++g_dlamch1_lt;
/*<             A = A*LBETA >*/
	    a *= g_dlamc1_lbeta;
/*<             C = DLAMC3( A, ONE ) >*/
	    c__ = dlamc3_(&a, &one);
/*<             C = DLAMC3( C, -A ) >*/
	    d__1 = -a;
	    c__ = dlamc3_(&c__, &d__1);
/*<             GO TO 30 >*/
	    goto L30;
/*<          END IF >*/
	}
/* +       END WHILE */

/*<       END IF >*/
    }

/*<       BETA = LBETA >*/
    *beta = g_dlamc1_lbeta;
/*<       T = LT >*/
    *t = g_dlamch1_lt;
/*<       RND = LRND >*/
    *rnd = g_dlamc1_lrnd;
/*<       IEEE1 = LIEEE1 >*/
    *ieee1 = g_lieee1;
/*<       RETURN >*/
    return 0;

/*     End of DLAMC1 */

/*<       END >*/
} /* dlamc1_ */


/* *********************************************************************** */

/*<       SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX ) >*/
/* Subroutine */ int dlamc2_(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax)
{
    /* Initialized data */

    static logical first = TRUE_;
    /**/static/**/ logical iwarn = FALSE_;

    /* Format strings */
    static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre\
ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the value EMIN loo\
ks\002,\002 acceptable please comment out \002,/\002 the IF block as marked \
within the code of routine\002,\002 DLAMC2,\002,/\002 otherwise supply EMIN \
explicitly.\002,/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
	
	/* grn: these MUST be static because
	  they won't change on new calls because the machine
	  won't change between calls.
	  The boolean "first" determines if everything is 
	  initialized properly */

    static logical ieee;
    static doublereal half;
    static doublereal zero, a, b, c__;
    static integer i__;
    static doublereal rbase;
    static integer gnmin;
    static doublereal small;
    static integer gpmin;
    static doublereal third, sixth;

	/* placed back from global to local static */
	static logical g_dlamc2_lrnd;
	static integer g_dlamc2_lt; /* carefull also in dlamc1_ */
	static integer g_dlamc2_lbeta;
	static doublereal g_leps;
	static integer g_lemin, g_lemax;
	static doublereal g_lrmin, g_lrmax;


    extern /* Subroutine */ int dlamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal dlamc3_(doublereal *, doublereal *);

    static logical lieee1;
    
	extern /* Subroutine */ int dlamc4_(integer *, doublereal *, integer *), 
	    dlamc5_(integer *, integer *, integer *, logical *, integer *, 
	    doublereal *);
    static integer ngnmin, ngpmin;
    static doublereal one, two;

    /* Fortran I/O blocks */
    static cilist io___324 = { 0, 6, 0, fmt_9999, 0 };

/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*<       LOGICAL            RND >*/
/*<       INTEGER            BETA, EMAX, EMIN, T >*/
/*<       DOUBLE PRECISION   EPS, RMAX, RMIN >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC2 determines the machine parameters specified in its argument */
/*  list. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) INTEGER */
/*          The base of the machine. */

/*  T       (output) INTEGER */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) LOGICAL */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */

/*  EPS     (output) DOUBLE PRECISION */
/*          The smallest positive number such that */

/*             fl( 1.0 - EPS ) .LT. 1.0, */

/*          where fl denotes the computed value. */

/*  EMIN    (output) INTEGER */
/*          The minimum exponent before (gradual) underflow occurs. */

/*  RMIN    (output) DOUBLE PRECISION */
/*          The smallest normalized number for the machine, given by */
/*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
/*          of BETA. */

/*  EMAX    (output) INTEGER */
/*          The maximum exponent before overflow occurs. */

/*  RMAX    (output) DOUBLE PRECISION */
/*          The largest positive number for the machine, given by */
/*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
/*          value of BETA. */

/*  Further Details */
/*  =============== */

/*  The computation of  EPS  is based on a routine PARANOIA by */
/*  W. Kahan of the University of California at Berkeley. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*<       LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND >*/
/*<    >*/
/*<    >*/
/*     .. */
/*     .. External Functions .. */
/*<       DOUBLE PRECISION   DLAMC3 >*/
/*<       EXTERNAL           DLAMC3 >*/
/*     .. */
/*     .. External Subroutines .. */
/*<       EXTERNAL           DLAMC1, DLAMC4, DLAMC5 >*/
/*     .. */
/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          ABS, MAX, MIN >*/
/*     .. */
/*     .. Save statement .. */
/*<    >*/
/*     .. */
/*     .. Data statements .. */
/*<       DATA               FIRST / .TRUE. / , IWARN / .FALSE. / >*/
/*     .. */
/*     .. Executable Statements .. */

/*<       IF( FIRST ) THEN >*/
    if (first) {
/*<          FIRST = .FALSE. >*/
	first = FALSE_;
/*<          ZERO = 0 >*/
	zero = 0.;
/*<          ONE = 1 >*/
	one = 1.;
/*<          TWO = 2 >*/
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */

/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are stored  and not held in registers,  or */
/*        are not affected by optimizers. */

/*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

/*<          CALL DLAMC1( LBETA, LT, LRND, LIEEE1 ) >*/
	dlamc1_(&g_dlamc2_lbeta, &g_dlamc2_lt, &g_dlamc2_lrnd, &lieee1);

/*        Start to find EPS. */

/*<          B = LBETA >*/
	b = (doublereal) g_dlamc2_lbeta;
/*<          A = B**( -LT ) >*/
	i__1 = -g_dlamc2_lt;
	a = pow_di(&b, &i__1);
/*<          LEPS = A >*/
	g_leps = a;

/*        Try some tricks to see whether or not this is the correct  EPS. */

/*<          B = TWO / 3 >*/
	b = two / 3;
/*<          HALF = ONE / 2 >*/
	half = one / 2;
/*<          SIXTH = DLAMC3( B, -HALF ) >*/
	d__1 = -half;
	sixth = dlamc3_(&b, &d__1);
/*<          THIRD = DLAMC3( SIXTH, SIXTH ) >*/
	third = dlamc3_(&sixth, &sixth);
/*<          B = DLAMC3( THIRD, -HALF ) >*/
	d__1 = -half;
	b = dlamc3_(&third, &d__1);
/*<          B = DLAMC3( B, SIXTH ) >*/
	b = dlamc3_(&b, &sixth);
/*<          B = ABS( B ) >*/
	b = abs(b);
/*<    >*/
	if (b < g_leps) {
	    b = g_leps;
	}

/*<          LEPS = 1 >*/
	g_leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
/*<    10    CONTINUE >*/
L10:
/*<          IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN >*/
	if (g_leps > b && b > zero) {
/*<             LEPS = B >*/
	    g_leps = b;
/*<             C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) ) >*/
	    d__1 = half * g_leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = g_leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c__ = dlamc3_(&d__1, &d__2);
/*<             C = DLAMC3( HALF, -C ) >*/
	    d__1 = -c__;
	    c__ = dlamc3_(&half, &d__1);
/*<             B = DLAMC3( HALF, C ) >*/
	    b = dlamc3_(&half, &c__);
/*<             C = DLAMC3( HALF, -B ) >*/
	    d__1 = -b;
	    c__ = dlamc3_(&half, &d__1);
/*<             B = DLAMC3( HALF, C ) >*/
	    b = dlamc3_(&half, &c__);
/*<             GO TO 10 >*/
	    goto L10;
/*<          END IF >*/
	}
/* +       END WHILE */

/*<    >*/
	if (a < g_leps) {
	    g_leps = a;
	}

/*        Computation of EPS complete. */

/*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). */
/*        Keep dividing  A by BETA until (gradual) underflow occurs. This */
/*        is detected when we cannot recover the previous A. */

/*<          RBASE = ONE / LBETA >*/
	rbase = one / g_dlamc2_lbeta;
/*<          SMALL = ONE >*/
	small = one;
/*<          DO 20 I = 1, 3 >*/
	for (i__ = 1; i__ <= 3; ++i__) {
/*<             SMALL = DLAMC3( SMALL*RBASE, ZERO ) >*/
	    d__1 = small * rbase;
	    small = dlamc3_(&d__1, &zero);
/*<    20    CONTINUE >*/
/* L20: */
	}
/*<          A = DLAMC3( ONE, SMALL ) >*/
	a = dlamc3_(&one, &small);
/*<          CALL DLAMC4( NGPMIN, ONE, LBETA ) >*/
	dlamc4_(&ngpmin, &one, &g_dlamc2_lbeta);
/*<          CALL DLAMC4( NGNMIN, -ONE, LBETA ) >*/
	d__1 = -one;
	dlamc4_(&ngnmin, &d__1, &g_dlamc2_lbeta);
/*<          CALL DLAMC4( GPMIN, A, LBETA ) >*/
	dlamc4_(&gpmin, &a, &g_dlamc2_lbeta);
/*<          CALL DLAMC4( GNMIN, -A, LBETA ) >*/
	d__1 = -a;
	dlamc4_(&gnmin, &d__1, &g_dlamc2_lbeta);
/*<          IEEE = .FALSE. >*/
	ieee = FALSE_;

/*<          IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN >*/
	if (ngpmin == ngnmin && gpmin == gnmin) {
/*<             IF( NGPMIN.EQ.GPMIN ) THEN >*/
	    if (ngpmin == gpmin) {
/*<                LEMIN = NGPMIN >*/
		g_lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual underflow; */
/*              e.g.,  VAX ) */
/*<             ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN >*/
	    } else if (gpmin - ngpmin == 3) {
/*<                LEMIN = NGPMIN - 1 + LT >*/
		g_lemin = ngpmin - 1 + g_dlamc2_lt;
/*<                IEEE = .TRUE. >*/
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual underflow; */
/*              e.g., IEEE standard followers ) */
/*<             ELSE >*/
	    } else {
/*<                LEMIN = MIN( NGPMIN, GPMIN ) >*/
		g_lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
/*<                IWARN = .TRUE. >*/
		iwarn = TRUE_;
/*<             END IF >*/
	    }

/*<          ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN >*/
	} else if (ngpmin == gpmin && ngnmin == gnmin) {
/*<             IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN >*/
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
/*<                LEMIN = MAX( NGPMIN, NGNMIN ) >*/
		g_lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow; */
/*              e.g., CYBER 205 ) */
/*<             ELSE >*/
	    } else {
/*<                LEMIN = MIN( NGPMIN, NGNMIN ) >*/
		g_lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
/*<                IWARN = .TRUE. >*/
		iwarn = TRUE_;
/*<             END IF >*/
	    }

/*<    >*/
	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
/*<             IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN >*/
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
/*<                LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT >*/
		g_lemin = max(ngpmin,ngnmin) - 1 + g_dlamc2_lt;
/*            ( Twos-complement machines with gradual underflow; */
/*              no known machine ) */
/*<             ELSE >*/
	    } else {
/*<                LEMIN = MIN( NGPMIN, NGNMIN ) >*/
		g_lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
/*<                IWARN = .TRUE. >*/
		iwarn = TRUE_;
/*<             END IF >*/
	    }

/*<          ELSE >*/
	} else {
/*<             LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN ) >*/
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    g_lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
/*<             IWARN = .TRUE. >*/
	    iwarn = TRUE_;
/*<          END IF >*/
	}
/* ** */
/* Comment out this if block if EMIN is ok */
/*<          IF( IWARN ) THEN >*/
	if (iwarn) {
/*<             FIRST = .TRUE. >*/
	    first = TRUE_;
/*<             WRITE( 6, FMT = 9999 )LEMIN >*/
	    s_wsfe(&io___324);
	    do_fio(&c__1, (char *)&g_lemin, (ftnlen)sizeof(integer));
	    e_wsfe();
/*<          END IF >*/
	}
/* ** */

/*        Assume IEEE arithmetic if we found denormalised  numbers above, */
/*        or if arithmetic seems to round in the  IEEE style,  determined */
/*        in routine DLAMC1. A true IEEE machine should have both  things */
/*        true; however, faulty machines may have one or the other. */

/*<          IEEE = IEEE .OR. LIEEE1 >*/
	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could compute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
/*        this computation. */

/*<          LRMIN = 1 >*/
	g_lrmin = 1.;
/*<          DO 30 I = 1, 1 - LEMIN >*/
	i__1 = 1 - g_lemin;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             LRMIN = DLAMC3( LRMIN*RBASE, ZERO ) >*/
	    d__1 = g_lrmin * rbase;
	    g_lrmin = dlamc3_(&d__1, &zero);
/*<    30    CONTINUE >*/
/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

/*<          CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX ) >*/
	dlamc5_(&g_dlamc2_lbeta, &g_dlamc2_lt, &g_lemin, &ieee, &g_lemax, &g_lrmax);
/*<       END IF >*/
    }

/*<       BETA = LBETA >*/
    *beta = g_dlamc2_lbeta;
/*<       T = LT >*/
    *t = g_dlamc2_lt;
/*<       RND = LRND >*/
    *rnd = g_dlamc2_lrnd;
/*<       EPS = LEPS >*/
    *eps = g_leps;
/*<       EMIN = LEMIN >*/
    *emin = g_lemin;
/*<       RMIN = LRMIN >*/
    *rmin = g_lrmin;
/*<       EMAX = LEMAX >*/
    *emax = g_lemax;
/*<       RMAX = LRMAX >*/
    *rmax = g_lrmax;

/*<       RETURN >*/
    return 0;

/*<  9 >*/

/*     End of DLAMC2 */

/*<       END >*/
} /* dlamc2_ */


/* *********************************************************************** */

/*<       DOUBLE PRECISION FUNCTION DLAMC3( A, B ) >*/
doublereal dlamc3_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*<       DOUBLE PRECISION   A, B >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing */
/*  the addition of  A  and  B ,  for use in situations where optimizers */
/*  might hold one of these in a register. */

/*  Arguments */
/*  ========= */

/*  A, B    (input) DOUBLE PRECISION */
/*          The values A and B. */

/* ===================================================================== */

/*     .. Executable Statements .. */

/*<       DLAMC3 = A + B >*/
    ret_val = *a + *b;

/*<       RETURN >*/
    return ret_val;

/*     End of DLAMC3 */

/*<       END >*/
} /* dlamc3_ */


/* *********************************************************************** */

/*<       SUBROUTINE DLAMC4( EMIN, START, BASE ) >*/
/* Subroutine */ int dlamc4_(integer *emin, doublereal *start, integer *base)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
	
    /*static*/ doublereal zero, a;
    /*static*/ integer i__;
    /*static*/ doublereal rbase, b1, b2, c1, c2, d1, d2;

    extern doublereal dlamc3_(doublereal *, doublereal *);
    
	/*static*/ doublereal one;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*<       INTEGER            BASE, EMIN >*/
/*<       DOUBLE PRECISION   START >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC4 is a service routine for DLAMC2. */

/*  Arguments */
/*  ========= */

/*  EMIN    (output) EMIN */
/*          The minimum exponent before (gradual) underflow, computed by */
/*          setting A = START and dividing by BASE until the previous A */
/*          can not be recovered. */

/*  START   (input) DOUBLE PRECISION */
/*          The starting point for determining EMIN. */

/*  BASE    (input) INTEGER */
/*          The base of the machine. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*<       INTEGER            I >*/
/*<       DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO >*/
/*     .. */
/*     .. External Functions .. */
/*<       DOUBLE PRECISION   DLAMC3 >*/
/*<       EXTERNAL           DLAMC3 >*/
/*     .. */
/*     .. Executable Statements .. */

/*<       A = START >*/
    a = *start;
/*<       ONE = 1 >*/
    one = 1.;
/*<       RBASE = ONE / BASE >*/
    rbase = one / *base;
/*<       ZERO = 0 >*/
    zero = 0.;
/*<       EMIN = 1 >*/
    *emin = 1;
/*<       B1 = DLAMC3( A*RBASE, ZERO ) >*/
    d__1 = a * rbase;
    b1 = dlamc3_(&d__1, &zero);
/*<       C1 = A >*/
    c1 = a;
/*<       C2 = A >*/
    c2 = a;
/*<       D1 = A >*/
    d1 = a;
/*<       D2 = A >*/
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
/*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
/*<    10 CONTINUE >*/
L10:
/*<    >*/
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
/*<          EMIN = EMIN - 1 >*/
	--(*emin);
/*<          A = B1 >*/
	a = b1;
/*<          B1 = DLAMC3( A / BASE, ZERO ) >*/
	d__1 = a / *base;
	b1 = dlamc3_(&d__1, &zero);
/*<          C1 = DLAMC3( B1*BASE, ZERO ) >*/
	d__1 = b1 * *base;
	c1 = dlamc3_(&d__1, &zero);
/*<          D1 = ZERO >*/
	d1 = zero;
/*<          DO 20 I = 1, BASE >*/
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             D1 = D1 + B1 >*/
	    d1 += b1;
/*<    20    CONTINUE >*/
/* L20: */
	}
/*<          B2 = DLAMC3( A*RBASE, ZERO ) >*/
	d__1 = a * rbase;
	b2 = dlamc3_(&d__1, &zero);
/*<          C2 = DLAMC3( B2 / RBASE, ZERO ) >*/
	d__1 = b2 / rbase;
	c2 = dlamc3_(&d__1, &zero);
/*<          D2 = ZERO >*/
	d2 = zero;
/*<          DO 30 I = 1, BASE >*/
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             D2 = D2 + B2 >*/
	    d2 += b2;
/*<    30    CONTINUE >*/
/* L30: */
	}
/*<          GO TO 10 >*/
	goto L10;
/*<       END IF >*/
    }
/* +    END WHILE */

/*<       RETURN >*/
    return 0;

/*     End of DLAMC4 */

/*<       END >*/
} /* dlamc4_ */


/* *********************************************************************** */

/*<       SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX ) >*/
/* Subroutine */ int dlamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    /*static*/ integer lexp;
    /*static*/ integer uexp, i__;
    /*static*/ doublereal y, z__;
    /*static*/ integer nbits;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    /*static*/ doublereal recbas;
    /*static*/ integer exbits, expsum, try__;

	/* placed back from global to local static */
	static doublereal g_oldy;

/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*<       LOGICAL            IEEE >*/
/*<       INTEGER            BETA, EMAX, EMIN, P >*/
/*<       DOUBLE PRECISION   RMAX >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAMC5 attempts to compute RMAX, the largest machine floating-point */
/*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
/*  approximately to a power of 2.  It will fail on machines where this */
/*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, */
/*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
/*  too large (i.e. too close to zero), probably with overflow. */

/*  Arguments */
/*  ========= */

/*  BETA    (input) INTEGER */
/*          The base of floating-point arithmetic. */

/*  P       (input) INTEGER */
/*          The number of base BETA digits in the mantissa of a */
/*          floating-point value. */

/*  EMIN    (input) INTEGER */
/*          The minimum exponent before (gradual) underflow. */

/*  IEEE    (input) LOGICAL */
/*          A logical flag specifying whether or not the arithmetic */
/*          system is thought to comply with the IEEE standard. */

/*  EMAX    (output) INTEGER */
/*          The largest exponent before overflow */

/*  RMAX    (output) DOUBLE PRECISION */
/*          The largest machine floating-point number. */

/* ===================================================================== */

/*     .. Parameters .. */
/*<       DOUBLE PRECISION   ZERO, ONE >*/
/*<       PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 ) >*/
/*     .. */
/*     .. Local Scalars .. */
/*<       INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP >*/
/*<       DOUBLE PRECISION   OLDY, RECBAS, Y, Z >*/
/*     .. */
/*     .. External Functions .. */
/*<       DOUBLE PRECISION   DLAMC3 >*/
/*<       EXTERNAL           DLAMC3 >*/
/*     .. */
/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          MOD >*/
/*     .. */
/*     .. Executable Statements .. */

/*     First compute LEXP and UEXP, two powers of 2 that bound */
/*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
/*     approximately to the bound that is closest to abs(EMIN). */
/*     (EMAX is the exponent of the required number RMAX). */

/*<       LEXP = 1 >*/
    lexp = 1;
/*<       EXBITS = 1 >*/
    exbits = 1;
/*<    10 CONTINUE >*/
L10:
/*<       TRY = LEXP*2 >*/
    try__ = lexp << 1;
/*<       IF( TRY.LE.( -EMIN ) ) THEN >*/
    if (try__ <= -(*emin)) {
/*<          LEXP = TRY >*/
	lexp = try__;
/*<          EXBITS = EXBITS + 1 >*/
	++exbits;
/*<          GO TO 10 >*/
	goto L10;
/*<       END IF >*/
    }
/*<       IF( LEXP.EQ.-EMIN ) THEN >*/
    if (lexp == -(*emin)) {
/*<          UEXP = LEXP >*/
	uexp = lexp;
/*<       ELSE >*/
    } else {
/*<          UEXP = TRY >*/
	uexp = try__;
/*<          EXBITS = EXBITS + 1 >*/
	++exbits;
/*<       END IF >*/
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
/*     than or equal to EMIN. EXBITS is the number of bits needed to */
/*     store the exponent. */

/*<       IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN >*/
    if (uexp + *emin > -lexp - *emin) {
/*<          EXPSUM = 2*LEXP >*/
	expsum = lexp << 1;
/*<       ELSE >*/
    } else {
/*<          EXPSUM = 2*UEXP >*/
	expsum = uexp << 1;
/*<       END IF >*/
    }

/*     EXPSUM is the exponent range, approximately equal to */
/*     EMAX - EMIN + 1 . */

/*<       EMAX = EXPSUM + EMIN - 1 >*/
    *emax = expsum + *emin - 1;
/*<       NBITS = 1 + EXBITS + P >*/
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a */
/*     floating-point number. */

/*<       IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN >*/
    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a */
/*        floating-point number, which is unlikely, or some bits are */
/*        not used in the representation of numbers, which is possible, */
/*        (e.g. Cray machines) or the mantissa has an implicit bit, */
/*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the */
/*        most likely. We have to assume the last alternative. */
/*        If this is true, then we need to reduce EMAX by one because */
/*        there must be some way of representing zero in an implicit-bit */
/*        system. On machines like Cray, we are reducing EMAX by one */
/*        unnecessarily. */

/*<          EMAX = EMAX - 1 >*/
	--(*emax);
/*<       END IF >*/
    }

/*<       IF( IEEE ) THEN >*/
    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent */
/*        for infinity and NaN. */

/*<          EMAX = EMAX - 1 >*/
	--(*emax);
/*<       END IF >*/
    }

/*     Now create RMAX, the largest machine number, which should */
/*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

/*     First compute 1.0 - BETA**(-P), being careful that the */
/*     result is less than 1.0 . */

/*<       RECBAS = ONE / BETA >*/
    recbas = 1. / *beta;
/*<       Z = BETA - ONE >*/
    z__ = *beta - 1.;
/*<       Y = ZERO >*/
    y = 0.;
/*<       DO 20 I = 1, P >*/
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Z = Z*RECBAS >*/
	z__ *= recbas;
/*<    >*/
	if (y < 1.) {
	    g_oldy = y;
	}
/*<          Y = DLAMC3( Y, Z ) >*/
	y = dlamc3_(&y, &z__);
/*<    20 CONTINUE >*/
/* L20: */
    }
/*<    >*/
    if (y >= 1.) {
	y = g_oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

/*<       DO 30 I = 1, EMAX >*/
    i__1 = *emax;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          Y = DLAMC3( Y*BETA, ZERO ) >*/
	d__1 = y * *beta;
	y = dlamc3_(&d__1, &c_b318);
/*<    30 CONTINUE >*/
/* L30: */
    }

/*<       RMAX = Y >*/
    *rmax = y;
/*<       RETURN >*/
    return 0;

/*     End of DLAMC5 */

/*<       END >*/
} /* dlamc5_ */

/*<       LOGICAL          FUNCTION LSAME( CA, CB ) >*/
logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    /*static*/ integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER          CA, CB >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/* ===================================================================== */

/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          ICHAR >*/
/*     .. */
/*     .. Local Scalars .. */
/*<       INTEGER            INTA, INTB, ZCODE >*/
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */

/*<       LSAME = CA.EQ.CB >*/
    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
/*<    >*/
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

/*<       ZCODE = ICHAR( 'Z' ) >*/
    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

/*<       INTA = ICHAR( CA ) >*/
    inta = *(unsigned char *)ca;
/*<       INTB = ICHAR( CB ) >*/
    intb = *(unsigned char *)cb;

/*<       IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN >*/
    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
/*        upper case 'Z'. */

/*<          IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32 >*/
	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
/*<          IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32 >*/
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

/*<       ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN >*/
    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
/*        upper case 'Z'. */

/*<    >*/
	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
/*<    >*/
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

/*<       ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN >*/
    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
/*        plus 128 of either lower or upper case 'Z'. */

/*<          IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32 >*/
	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
/*<          IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32 >*/
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
/*<       END IF >*/
    }
/*<       LSAME = INTA.EQ.INTB >*/
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

/*<       END >*/
    return ret_val;
} /* lsame_ */


#ifdef __cplusplus
	}
#endif
