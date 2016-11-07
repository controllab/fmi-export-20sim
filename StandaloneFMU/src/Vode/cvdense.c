/******************************************************************
 *                                                                *
 * File          : cvdense.c                                      *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the CVODE dense linear     *
 * solver, CVDENSE.                                               *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "cvdense.h"
#include "cvode.h"
#include "dense.h"
#include "cvnltyps.h"
#include "cvvector.h"
#include "cvnlmath.h"


/* Error Messages */

#define CVDENSE_INIT  "CVDenseInit-- "

#define MSG_MEM_FAIL  CVDENSE_INIT "A memory request failed.\n\n"


/* Other Constants */

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)


/******************************************************************
 *                                                                *           
 * Types : CVDenseMemRec, CVDenseMemRec*                              *
 *----------------------------------------------------------------*
 * The type CVDenseMem is pointer to a CVDenseMemRec. This        *
 * structure contains CVDense solver-specific data.               *
 *                                                                *
 ******************************************************************/

typedef struct {

    CVDenseJacFn d_jac; /* jac = Jacobian routine to be called    */

    DenseMat *d_M;       /* M = I - gamma J, gamma = h / l1        */

    integer *d_pivots;  /* pivots = pivot array for PM = LU       */

    DenseMat *d_savedJ;  /* savedJ = old Jacobian                  */

    long int  d_nstlj;  /* nstlj = nst at last Jacobian eval.     */
    
    long int d_nje;     /* nje = no. of calls to jac              */

    void *d_J_data;     /* J_data is passed to jac                */

} CVDenseMemRec;//, *CVDenseMem;


/* CVDENSE linit, lsetup, lsolve, and lfree routines */
 
static int  CVDenseInit(CVodeMemData* cv_mem, bool *setupNonNull);

static int  CVDenseSetup(CVodeMemData* cv_mem, int convfail, iN_Vector* ypred,
			 iN_Vector* fpred, bool *jcurPtr, iN_Vector* vtemp1,
			 iN_Vector* vtemp2, iN_Vector* vtemp3);

static int  CVDenseSolve(CVodeMemData* cv_mem, iN_Vector* b, iN_Vector* ycur,
			 iN_Vector* fcur);

static void CVDenseFree(CVodeMemData* cv_mem);


/*************** CVDenseDQJac ****************************************

 This routine generates a dense difference quotient approximation to
 the Jacobian of f(t,y). It assumes that a dense matrix of type
 DenseMat is stored column-wise, and that elements within each column
 are contiguous. The address of the jth column of J is obtained via
 the macro DENSE_COL and an N_Vector with the jth column as the
 component array is created using N_VMAKE and N_VDATA. Finally, the
 actual computation of the jth column of the Jacobian is done with a
 call to N_VLinearSum.

**********************************************************************/
/* int g_jac_counter = 0; */
void CVDenseDQJac(integer N, DenseMat *J, RhsFn f, void *f_data, real tn,
		  iN_Vector* y, iN_Vector* fy, iN_Vector* ewt, real h, real uround,
		  void *jac_data, long int *nfePtr, iN_Vector* vtemp1,
		  iN_Vector* vtemp2, iN_Vector* vtemp3)
{
  real fnorm, minInc, inc, inc_inv, yjsaved, srur;
  real *y_data, *ewt_data;
  iN_Vector* ftemp, *jthCol;
  integer j;

  ftemp = vtemp1; /* Rename work vector for use as f vector value */

  /* Obtain pointers to the data for ewt, y */
  ewt_data   = N_VDATA(ewt);
  y_data     = N_VDATA(y);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * N * fnorm) : ONE;

  N_VMAKE(jthCol, NULL, N);

  /* this is the only for loop for 0..N-1 in CVODE */
  for (j=0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */
    
    N_VDATA(jthCol) = DENSE_COL(J,j);
    yjsaved = y_data[j];
    inc = MAX(srur*ABS(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;
    f(N, tn, y, ftemp, f_data);
    
	inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);
    y_data[j] = yjsaved;
  }

  N_VDISPOSE(jthCol);

  /* Increment counter nfe = *nfePtr */
  *nfePtr += N;
}

/* Readability Replacements */
/* grn: 2012-03-07, just put them in the code, it makes it more readable that way */
/*                  else you don't see if it is a local variable, or a member of the struct!!!!! */

/*#define N         (cv_mem->cv_N) */
/*#define lmm       (cv_mem->cv_lmm)*/
/*#define f         (cv_mem->cv_f)*/
/*#define f_data    (cv_mem->cv_f_data)*/
/*#define uround    (cv_mem->cv_uround)*/
/*#define nst       (cv_mem->cv_nst)*/
/*#define tn        (cv_mem->cv_tn)*/
/*#define h         (cv_mem->cv_h)*/
/*#define gamma     (cv_mem->cv_gamma)*/
/*#define gammap    (cv_mem->cv_gammap)*/
/*#define gamrat    (cv_mem->cv_gamrat)*/
/*#define ewt       (cv_mem->cv_ewt)*/
/*#define nfe       (cv_mem->cv_nfe)*/
/*#define errfp     (cv_mem->cv_errfp)*/
/*#define iopt      (cv_mem->cv_iopt)*/
/*#define linit     (cv_mem->cv_linit)*/
/*#define lsetup    (cv_mem->cv_lsetup)*/
/*#define lsolve    (cv_mem->cv_lsolve)*/
/*#define lfree     (cv_mem->cv_lfree)*/
/*#define lmem      (cv_mem->cv_lmem)*/

/*#define jac       (cvdense_mem->d_jac)*/
/*#define M         (cvdense_mem->d_M)*/
/*#define pivots    (cvdense_mem->d_pivots)*/
/*#define savedJ    (cvdense_mem->d_savedJ)*/
/*#define nstlj     (cvdense_mem->d_nstlj)*/
/*#define nje       (cvdense_mem->d_nje)*/
/*#define J_data    (cvdense_mem->d_J_data)*/

                  
/*************** CVDense *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the dense linear solver module. CVDense sets the
 cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 to be CVDenseInit, CVDenseSetup, CVDenseSolve, and CVDenseFree,
 respectively. It allocates memory for a structure of type
 CVDenseMemRec and sets the cv_lmem field in (*cvode_mem) to the
 address of this structure. Finally, it sets d_J_data field in the            
 CVDenseMemRec structure to be the input parameter jac_data and the
 d_jac field to be:                                         
                                                                
 (1) the input parameter djac if djac != NULL or                
                                                               
 (2) CVDenseDQJac if djac == NULL.                             

**********************************************************************/

void CVDense(void *cvode_mem, CVDenseJacFn djac, void *jac_data)
{
  CVodeMemData* cv_mem;
  CVDenseMemRec* cvdense_mem;

  /* Return immediately if cvode_mem is NULL */
  cv_mem = (CVodeMemData*) cvode_mem;
  if (cv_mem == NULL) return;  /* CVode reports this error */

  /* Set four main function fields in cv_mem */
  cv_mem->cv_linit  = CVDenseInit;
  cv_mem->cv_lsetup = CVDenseSetup;
  cv_mem->cv_lsolve = CVDenseSolve;
  cv_mem->cv_lfree  = CVDenseFree;

  /* Get memory for CVDenseMemRec */
  cv_mem->cv_lmem = cvdense_mem = (CVDenseMemRec*) malloc(sizeof(CVDenseMemRec));
  if (cvdense_mem == NULL) return;  /* CVDenseInit reports this error */

  // grn: 2012-03-05, make the whole thing NULL, this will be
  // good when we're going to "Remember" values
  memset(cvdense_mem, 0, sizeof(CVDenseMemRec));
  /* Set Jacobian routine field to user's djac or CVDenseDQJac */
  if (djac == NULL) {
    cvdense_mem->d_jac = CVDenseDQJac;
  } else {
	cvdense_mem->d_jac = djac;
  }
  cvdense_mem->d_J_data = jac_data;
}

/*************** CVDenseInit *****************************************

 This routine initializes remaining memory specific to the dense
 linear solver.  If any memory request fails, all memory previously
 allocated is freed, and an error message printed, before returning.

**********************************************************************/

static int CVDenseInit(CVodeMemData* cv_mem, bool *setupNonNull)
{
  CVDenseMemRec* cvdense_mem;
  
  cvdense_mem = (CVDenseMemRec*) (cv_mem->cv_lmem);

  /* Print error message and return if cvdense_mem is NULL */  
  if (cvdense_mem == NULL) {
//    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }

  /* Set flag setupNonNull = TRUE */  
  *setupNonNull = TRUE;
  
  /* Allocate memory for M, savedJ, and pivot array */
  
  (cvdense_mem->d_M) = DenseAllocMat((cv_mem->cv_N));
  if ((cvdense_mem->d_M) == NULL) {
//    fprintf(errfp, MSG_MEM_FAIL);
    return(LINIT_ERR);
  }
  cvdense_mem->d_savedJ = DenseAllocMat((cv_mem->cv_N));
  if ((cvdense_mem->d_savedJ) == NULL) {
//    fprintf(errfp, MSG_MEM_FAIL);
    DenseFreeMat(cvdense_mem->d_M);
    return(LINIT_ERR);
  }
  cvdense_mem->d_pivots = DenseAllocPiv(cv_mem->cv_N);
  if ((cvdense_mem->d_pivots) == NULL) {
//    fprintf(errfp, MSG_MEM_FAIL);
    DenseFreeMat(cvdense_mem->d_M);
    DenseFreeMat(cvdense_mem->d_savedJ);
    return(LINIT_ERR);
  }
  
  /* Initialize nje and nstlj, and set workspace lengths */
   
  cvdense_mem->d_nje = 0;
  if ((cv_mem->cv_iopt) != NULL) {
   (cv_mem->cv_iopt)[DENSE_NJE] = cvdense_mem->d_nje;
   (cv_mem->cv_iopt)[DENSE_LRW] = 2*(cv_mem->cv_N)*(cv_mem->cv_N);
   (cv_mem->cv_iopt)[DENSE_LIW] = cv_mem->cv_N;
  }
  cvdense_mem->d_nstlj = 0;
  
  return(LINIT_OK);
}

/*************** CVDenseSetup ****************************************

 This routine does the setup operations for the dense linear solver.
 It makes a decision whether or not to call the Jacobian evaluation
 routine based on various state variables, and if not it uses the 
 saved copy.  In any case, it constructs the Newton matrix 
 M = I - gamma*J, updates counters, and calls the dense LU 
 factorization routine.

**********************************************************************/

static int CVDenseSetup(CVodeMemData* cv_mem, int convfail, iN_Vector* ypred,
			iN_Vector* fpred, bool *jcurPtr, iN_Vector* vtemp1,
			iN_Vector* vtemp2, iN_Vector* vtemp3)
{
  bool jbad, jok;
  real dgamma;
  integer ier;
  CVDenseMemRec* cvdense_mem;
  
  cvdense_mem = (CVDenseMemRec*) (cv_mem->cv_lmem);
 
  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
 
  dgamma = ABS(((cv_mem->cv_gamma)/(cv_mem->cv_gammap)) - ONE);
  jbad = ((cv_mem->cv_nst) == 0) || ((cv_mem->cv_nst) > (cvdense_mem->d_nstlj) + CVD_MSBJ) ||
         ((convfail == FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (convfail == FAIL_OTHER);
  jok = !jbad;
 
  if (jok) {
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    DenseCopy(cvdense_mem->d_savedJ, cvdense_mem->d_M);
  } else {
    /* If jok = FALSE, call jac routine for new J value */
    cvdense_mem->d_nje++;
    if ((cv_mem->cv_iopt) != NULL)
		(cv_mem->cv_iopt)[DENSE_NJE] = cvdense_mem->d_nje;
    cvdense_mem->d_nstlj = cv_mem->cv_nst;
    *jcurPtr = TRUE;
    DenseZero(cvdense_mem->d_M); 
    (cvdense_mem->d_jac)(cv_mem->cv_N,
						 cvdense_mem->d_M,
						 cv_mem->cv_f,
						 cv_mem->cv_f_data,
						 cv_mem->cv_tn,
						 ypred, fpred,
						 cv_mem->cv_ewt,
						 cv_mem->cv_h,
						cv_mem->cv_uround,
						cvdense_mem->d_J_data,
						&(cv_mem->cv_nfe),
						vtemp1, vtemp2, vtemp3);
    DenseCopy(cvdense_mem->d_M, cvdense_mem->d_savedJ);
  }
  
  /* Scale and add I to get M = I - gamma*J */
  DenseScale(-(cv_mem->cv_gamma), cvdense_mem->d_M);
  DenseAddI(cvdense_mem->d_M);

  /* Do LU factorization of M */
  ier = DenseFactor(cvdense_mem->d_M, cvdense_mem->d_pivots); 
  
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0)
	  return(1);
  return(0);
}

/*************** CVDenseSolve ****************************************

 This routine handles the solve operation for the dense linear solver
 by calling the dense backsolve routine.  The returned value is 0.

**********************************************************************/

static int CVDenseSolve(CVodeMemData* cv_mem, iN_Vector* b, iN_Vector* ycur,
			iN_Vector* fcur)
{
  CVDenseMemRec* cvdense_mem;
  
  cvdense_mem = (CVDenseMemRec*) (cv_mem->cv_lmem);
  
  DenseBacksolve(cvdense_mem->d_M, cvdense_mem->d_pivots, b);

  /* If BDF, scale the correction to account for change in gamma */
  if (((cv_mem->cv_lmm) == BDF) && ((cv_mem->cv_gamrat) != ONE)) {
    N_VScale(TWO/(ONE + (cv_mem->cv_gamrat)), b, b);
  }

  return(0);
}

/*************** CVDenseFree *****************************************

 This routine frees memory specific to the dense linear solver.

**********************************************************************/

static void CVDenseFree(CVodeMemData* cv_mem)
{
  CVDenseMemRec*  cvdense_mem;

  cvdense_mem = (CVDenseMemRec*) (cv_mem->cv_lmem);
  
  DenseFreeMat(cvdense_mem->d_M);
  DenseFreeMat(cvdense_mem->d_savedJ);
  DenseFreePiv(cvdense_mem->d_pivots);
  free(cv_mem->cv_lmem);
}


/**********************************
 Some routines to remember some memory
 needed to make a step back in time!!!
 grn: 2007-04-09
 *********************************/
void *MyAllocateDenseMemory(CVodeMemData* cv_mem)
{
	CVDenseMemRec *denseStruct = malloc(sizeof(CVDenseMemRec));

	/* grn: 2012-03-07 initialize on NULL */
	memset(denseStruct, 0, sizeof(CVDenseMemRec));

	denseStruct->d_M = DenseAllocMat(cv_mem->cv_N);
	denseStruct->d_savedJ = DenseAllocMat(cv_mem->cv_N);
	denseStruct->d_pivots = DenseAllocPiv(cv_mem->cv_N);

	// the data field is not used in our case

	return denseStruct;
}

void MyFreeDenseMemory(void *givenMemory)
{
	CVDenseMemRec *denseStruct = (CVDenseMemRec*)givenMemory;
	DenseFreeMat(denseStruct->d_M);
	DenseFreeMat(denseStruct->d_savedJ);
	DenseFreePiv(denseStruct->d_pivots);

	free(givenMemory);
}

void MyRememberDenseState(void *givenMemory, CVodeMemData* cv_mem)
{
	CVDenseMemRec *denseStruct = (CVDenseMemRec*)givenMemory;
	CVDenseMemRec *sourceStruct = (CVDenseMemRec*)cv_mem->cv_lmem;
	if( sourceStruct == NULL )
	{
		// nothing to remember yet, maybe later
		return;
	}
	DenseCopy(sourceStruct->d_M, denseStruct->d_M);
	DenseCopy(sourceStruct->d_savedJ, denseStruct->d_savedJ);

	if( sourceStruct->d_pivots != NULL )
	{
		if( denseStruct->d_pivots == NULL )
		{
			// do allocation of the pivot memory
			denseStruct->d_pivots = DenseAllocPiv(cv_mem->cv_N);
		}

		memcpy(denseStruct->d_pivots, sourceStruct->d_pivots, cv_mem->cv_N * sizeof(int));
	}
	else
	{
		// release some memory
		DenseFreePiv(denseStruct->d_pivots);
		denseStruct->d_pivots = NULL;
	}
}

void CopyiNVector(iN_Vector *dest, iN_Vector *src, int size)
{
	if( dest == NULL || src == NULL )
	{
		return;
	}

	if( dest->length != src->length )
		return;

	memcpy(dest->data, src->data, dest->length * sizeof(real));
}

void MyRememberVodeState(CVodeMemData *prevMemory, void *prevDenseMemory, CVodeMemData *cv_mem)
{
	int maxord;
	int i;
	int N;
	prevMemory->cv_uround = cv_mem->cv_uround;		/* machine unit roundoff */
	prevMemory->cv_N = cv_mem->cv_N;				/* ODE system size             */
	prevMemory->cv_f = cv_mem->cv_f;				/* y' = f(t,y(t))              */
	prevMemory->cv_f_data = cv_mem->cv_f_data;		/* user pointer passed to f    */
	prevMemory->cv_lmm = cv_mem->cv_lmm;			/* lmm = ADAMS or BDF          */
	prevMemory->cv_iter = cv_mem->cv_iter;			/* iter = FUNCTIONAL or NEWTON */
	prevMemory->cv_itol = cv_mem->cv_itol;			/* itol = SS or SV             */
	prevMemory->cv_reltol = cv_mem->cv_reltol;		/* ptr to relative tolerance   */
	prevMemory->cv_abstol = cv_mem->cv_abstol;		/* ptr to absolute tolerance   */

	// determine the lenght of the cv_zn
	maxord = maxord = (cv_mem->cv_lmm == ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

	/* Nordsieck History Array */

	N = cv_mem->cv_N;
	for( i = 0; i < maxord; i++)
	{
		CopyiNVector(prevMemory->cv_zn[i], cv_mem->cv_zn[i], N);
	}

	/* Vectors of length N */
	CopyiNVector(prevMemory->cv_ewt, cv_mem->cv_ewt, N);		/* error weight vector                          */
	CopyiNVector(prevMemory->cv_y, cv_mem->cv_y, N);			/* y is used as temporary storage by the solver */
																/* The memory is provided by the user to CVode  */
																/* where the vector is named yout.              */
	CopyiNVector(prevMemory->cv_acor, cv_mem->cv_acor, N);		/* In the context of the solution of the        */
																/* nonlinear equation, acor = y_n(m) - y_n(0).  */
																/* On return, this vector is scaled to give     */
																/* the estimated local error in y.              */
	CopyiNVector(prevMemory->cv_tempv, cv_mem->cv_tempv, N);	/* temporary storage vector                     */
	CopyiNVector(prevMemory->cv_ftemp, cv_mem->cv_ftemp, N);	/* temporary storage vector                     */



  /* Step Data */

	prevMemory->cv_q = cv_mem->cv_q;				/* current order                           */
	prevMemory->cv_qprime = cv_mem->cv_qprime;		/* order to be used on the next step       */ 
													/* = q-1, q, or q+1                        */
	prevMemory->cv_qwait = cv_mem->cv_qwait;		/* number of internal steps to wait before */
													/* considering a change in q               */
	prevMemory->cv_L = cv_mem->cv_L;				/* L = q + 1                               */
	prevMemory->cv_h = cv_mem->cv_h;				/* current step size                     */
	prevMemory->cv_hprime = cv_mem->cv_hprime;		/* step size to be used on the next step */ 
	prevMemory->cv_eta = cv_mem->cv_eta;			/* eta = hprime / h                      */
	prevMemory->cv_hscale = cv_mem->cv_hscale;		/* value of h used in zn                 */
	prevMemory->cv_tn = cv_mem->cv_tn;				/* current internal value of t           */

	memcpy(prevMemory->cv_tau, cv_mem->cv_tau, L_MAX * sizeof(real));	/* vector of previous q+1 successful step    */
																		/* sizes indexed from 1 to q+1               */
	memcpy(prevMemory->cv_tq, cv_mem->cv_tq, NUM_TESTS * sizeof(real));	/* vector of test quantities indexed from    */
																		/* 1 to NUM_TESTS(=5)                        */
  
	memcpy(prevMemory->cv_l, cv_mem->cv_l, L_MAX * sizeof(real));		/* coefficients of l(x) (degree q poly)      */
     

	prevMemory->cv_rl1 = cv_mem->cv_rl1;			/* 1 / l[1]                     */
	prevMemory->cv_gamma = cv_mem->cv_gamma;		/* gamma = h * rl1              */
	prevMemory->cv_gammap = cv_mem->cv_gammap;		/* gamma at the last setup call */
	prevMemory->cv_gamrat = cv_mem->cv_gamrat;		/* gamma / gammap               */

	prevMemory->cv_crate = cv_mem->cv_crate;		/* estimated corrector convergence rate */
	prevMemory->cv_acnrm = cv_mem->cv_acnrm;		/* | acor | wrms                        */
	prevMemory->cv_mnewt = cv_mem->cv_mnewt;		/* Newton iteration counter             */

	/* Limits */

	prevMemory->cv_qmax = cv_mem->cv_qmax;			/* q <= qmax                                          */
	prevMemory->cv_mxstep = cv_mem->cv_mxstep;		/* maximum number of internal steps for one user call */
	prevMemory->cv_maxcor = cv_mem->cv_maxcor;		/* maximum number of corrector iterations for the     */
													/* solution of the nonlinear equation                 */
	prevMemory->cv_mxhnil = cv_mem->cv_mxhnil;		/* maximum number of warning messages issued to the   */
													/* user that t + h == t for the next internal step    */

	prevMemory->cv_hmin = cv_mem->cv_hmin;			/* |h| >= hmin       */
	prevMemory->cv_hmax_inv = cv_mem->cv_hmax_inv;	/* |h| <= 1/hmax_inv */
	prevMemory->cv_etamax = cv_mem->cv_etamax;		/* eta <= etamax     */

  /* Counters */

	prevMemory->cv_nst = cv_mem->cv_nst;			/* number of internal steps taken             */
	prevMemory->cv_nfe = cv_mem->cv_nfe;			/* number of f calls                          */
	prevMemory->cv_ncfn = cv_mem->cv_ncfn;			/* number of corrector convergence failures   */
	prevMemory->cv_netf = cv_mem->cv_netf;			/* number of error test failures              */
	prevMemory->cv_nni = cv_mem->cv_nni;			/* number of Newton iterations performed      */
	prevMemory->cv_nsetups = cv_mem->cv_nsetups;	/* number of setup calls                      */
	prevMemory->cv_nhnil = cv_mem->cv_nhnil;		/* number of messages issued to the user that */
													/* t + h == t for the next iternal step       */
	prevMemory->cv_lrw = cv_mem->cv_lrw;			/* number of real words in CVODE work vectors */
	prevMemory->cv_liw = cv_mem->cv_liw;			/* no. of integer words in CVODE work vectors */

	/* Linear Solver Data */

	/* Linear Solver functions to be called */

	/* Linear Solver specific memory */
	MyRememberDenseState(prevDenseMemory, cv_mem);

	/* Flag to indicate successful cv_linit call */

	prevMemory->cv_linitOK = cv_mem->cv_linitOK;

	/* Saved Values */

	prevMemory->cv_qu = cv_mem->cv_qu;						/* last successful q value used   */
	prevMemory->cv_nstlp = cv_mem->cv_nstlp;				/* step number of last setup call */
	prevMemory->cv_hu = cv_mem->cv_hu;						/* last successful h value used   */
	prevMemory->cv_saved_tq5 = cv_mem->cv_saved_tq5;		/* saved value of tq[5]           */
	prevMemory->cv_imxer = cv_mem->cv_imxer;				/* index of max value of          */
															/* |acor[i]|*ewt[i]               */
	prevMemory->cv_jcur = cv_mem->cv_jcur;					/* Is the Jacobian info used by   */
															/* linear solver current?         */
	prevMemory->cv_tolsf = cv_mem->cv_tolsf;				/* tolerance scale factor         */
	prevMemory->cv_setupNonNull = cv_mem->cv_setupNonNull;	/* Does setup do something?       */

	/* Arrays for Optional Input and Optional Output */

	/*cv_iopt; */  /* long int optional input, output */
	/*cv_ropt; */ /* real optional input, output     */

	/* Error File */

	/*cv_errfp; */      /* CVODE error messages are sent to errfp */

	/* Pointer to Machine Environment-Specific Information */

	/* cv_machenv */;
}


void MySetBackRememberedDenseState(void *givenMemory, CVodeMemData* cv_mem)
{
	CVDenseMemRec *denseStruct = (CVDenseMemRec*)givenMemory;
	CVDenseMemRec *destStruct = (CVDenseMemRec*)cv_mem->cv_lmem;

	if( destStruct == NULL )
	{
		// nothing to remember yet, maybe later
		return;
	}

	// if the destination is NULL, it is fine, apparently it was not necessary
	DenseCopy(denseStruct->d_M, destStruct->d_M);
	DenseCopy(denseStruct->d_savedJ, destStruct->d_savedJ);

#if 0
	if( denseStruct->d_pivots != NULL && destStruct->d_pivots != NULL )
		memcpy(destStruct->d_pivots, denseStruct->d_pivots, cv_mem->cv_N * sizeof(int));
#else
	if( denseStruct->d_pivots != NULL )
	{
		if( destStruct->d_pivots == NULL )
		{
			// do allocation of the pivot memory
			destStruct->d_pivots = DenseAllocPiv(cv_mem->cv_N);
		}

		memcpy(destStruct->d_pivots, denseStruct->d_pivots, cv_mem->cv_N * sizeof(int));
	}
	else
	{
		// release some memory
		DenseFreePiv(destStruct->d_pivots);
		destStruct->d_pivots = NULL;
	}

#endif
}

void MySetbackVodeState(CVodeMemData *prevMemory, void *prevDenseMemory, CVodeMemData *cv_mem)
{
	int maxord;
	int i;
	int N;
	cv_mem->cv_uround = prevMemory->cv_uround;		/* machine unit roundoff */
	cv_mem->cv_N = prevMemory->cv_N;				/* ODE system size             */
	cv_mem->cv_f = prevMemory->cv_f;				/* y' = f(t,y(t))              */
	cv_mem->cv_f_data = prevMemory->cv_f_data;		/* user pointer passed to f    */
	cv_mem->cv_lmm = prevMemory->cv_lmm;			/* lmm = ADAMS or BDF          */
	cv_mem->cv_iter = prevMemory->cv_iter;			/* iter = FUNCTIONAL or NEWTON */
	cv_mem->cv_itol = prevMemory->cv_itol;			/* itol = SS or SV             */
	cv_mem->cv_reltol = prevMemory->cv_reltol;		/* ptr to relative tolerance   */
	cv_mem->cv_abstol = prevMemory->cv_abstol;		/* ptr to absolute tolerance   */

	// determine the lenght of the cv_zn
	maxord = maxord = (cv_mem->cv_lmm == ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

	N = cv_mem->cv_N;

	/* Nordsieck History Array */
	for( i = 0; i < maxord; i++)
	{
		CopyiNVector(cv_mem->cv_zn[i], prevMemory->cv_zn[i], N);
	}

	/* Vectors of length N */
	CopyiNVector(cv_mem->cv_ewt, prevMemory->cv_ewt, N);		/* error weight vector                          */
	CopyiNVector(cv_mem->cv_y, prevMemory->cv_y, N);			/* y is used as temporary storage by the solver */
																/* The memory is provided by the user to CVode  */
																/* where the vector is named yout.              */
	CopyiNVector(cv_mem->cv_acor, prevMemory->cv_acor, N);		/* In the context of the solution of the        */
																/* nonlinear equation, acor = y_n(m) - y_n(0).  */
																/* On return, this vector is scaled to give     */
																/* the estimated local error in y.              */
	CopyiNVector(cv_mem->cv_tempv, prevMemory->cv_tempv, N);	/* temporary storage vector                     */
	CopyiNVector(cv_mem->cv_ftemp, prevMemory->cv_ftemp, N);	/* temporary storage vector                     */



  /* Step Data */

	cv_mem->cv_q = prevMemory->cv_q;				/* current order                           */
	cv_mem->cv_qprime = prevMemory->cv_qprime;		/* order to be used on the next step       */ 
													/* = q-1, q, or q+1                        */
	cv_mem->cv_qwait = prevMemory->cv_qwait;		/* number of internal steps to wait before */
													/* considering a change in q               */
	cv_mem->cv_L = prevMemory->cv_L;				/* L = q + 1                               */
	cv_mem->cv_h = prevMemory->cv_h;				/* current step size                     */
	cv_mem->cv_hprime = prevMemory->cv_hprime;		/* step size to be used on the next step */ 
	cv_mem->cv_eta = prevMemory->cv_eta;			/* eta = hprime / h                      */
	cv_mem->cv_hscale = prevMemory->cv_hscale;		/* value of h used in zn                 */
	cv_mem->cv_tn = prevMemory->cv_tn;				/* current internal value of t           */

	memcpy(cv_mem->cv_tau, prevMemory->cv_tau, L_MAX * sizeof(real));	/* vector of previous q+1 successful step    */
																		/* sizes indexed from 1 to q+1               */
	memcpy(cv_mem->cv_tq, prevMemory->cv_tq, NUM_TESTS * sizeof(real));	/* vector of test quantities indexed from    */
																		/* 1 to NUM_TESTS(=5)                        */
  
	memcpy(cv_mem->cv_l, prevMemory->cv_l, L_MAX * sizeof(real));		/* coefficients of l(x) (degree q poly)      */
     

	cv_mem->cv_rl1 = prevMemory->cv_rl1;			/* 1 / l[1]                     */
	cv_mem->cv_gamma = prevMemory->cv_gamma;		/* gamma = h * rl1              */
	cv_mem->cv_gammap = prevMemory->cv_gammap;		/* gamma at the last setup call */
	cv_mem->cv_gamrat = prevMemory->cv_gamrat;		/* gamma / gammap               */

	cv_mem->cv_crate = prevMemory->cv_crate;		/* estimated corrector convergence rate */
	cv_mem->cv_acnrm = prevMemory->cv_acnrm;		/* | acor | wrms                        */
	cv_mem->cv_mnewt = prevMemory->cv_mnewt;		/* Newton iteration counter             */

	/* Limits */

	cv_mem->cv_qmax = prevMemory->cv_qmax;			/* q <= qmax                                          */
	cv_mem->cv_mxstep = prevMemory->cv_mxstep;		/* maximum number of internal steps for one user call */
	cv_mem->cv_maxcor = prevMemory->cv_maxcor;		/* maximum number of corrector iterations for the     */
													/* solution of the nonlinear equation                 */
	cv_mem->cv_mxhnil = prevMemory->cv_mxhnil;		/* maximum number of warning messages issued to the   */
													/* user that t + h == t for the next internal step    */

	cv_mem->cv_hmin = prevMemory->cv_hmin;			/* |h| >= hmin       */
	cv_mem->cv_hmax_inv = prevMemory->cv_hmax_inv;	/* |h| <= 1/hmax_inv */
	cv_mem->cv_etamax = prevMemory->cv_etamax;		/* eta <= etamax     */

  /* Counters */

	cv_mem->cv_nst = prevMemory->cv_nst;			/* number of internal steps taken             */
	cv_mem->cv_nfe = prevMemory->cv_nfe;			/* number of f calls                          */
	cv_mem->cv_ncfn = prevMemory->cv_ncfn;			/* number of corrector convergence failures   */
	cv_mem->cv_netf = prevMemory->cv_netf;			/* number of error test failures              */
	cv_mem->cv_nni = prevMemory->cv_nni;			/* number of Newton iterations performed      */
	cv_mem->cv_nsetups = prevMemory->cv_nsetups;	/* number of setup calls                      */
	cv_mem->cv_nhnil = prevMemory->cv_nhnil;		/* number of messages issued to the user that */
													/* t + h == t for the next iternal step       */
	cv_mem->cv_lrw = prevMemory->cv_lrw;			/* number of real words in CVODE work vectors */
	cv_mem->cv_liw = prevMemory->cv_liw;			/* no. of integer words in CVODE work vectors */

	/* Linear Solver Data */

	/* Linear Solver functions to be called */


	/* Linear Solver specific memory */
	MySetBackRememberedDenseState(prevDenseMemory, cv_mem);

	/* Flag to indicate successful cv_linit call */

	cv_mem->cv_linitOK = prevMemory->cv_linitOK;

	/* Saved Values */

	cv_mem->cv_qu = prevMemory->cv_qu;						/* last successful q value used   */
	cv_mem->cv_nstlp = prevMemory->cv_nstlp;				/* step number of last setup call */
	cv_mem->cv_hu = prevMemory->cv_hu;						/* last successful h value used   */
	cv_mem->cv_saved_tq5 = prevMemory->cv_saved_tq5;		/* saved value of tq[5]           */
	cv_mem->cv_imxer = prevMemory->cv_imxer;				/* index of max value of          */
															/* |acor[i]|*ewt[i]               */
	cv_mem->cv_jcur = prevMemory->cv_jcur;					/* Is the Jacobian info used by   */
															/* linear solver current?         */
	cv_mem->cv_tolsf = prevMemory->cv_tolsf;				/* tolerance scale factor         */
	cv_mem->cv_setupNonNull = prevMemory->cv_setupNonNull;	/* Does setup do something?       */

	/* Arrays for Optional Input and Optional Output */

	/*cv_iopt; */  /* long int optional input, output */
	/*cv_ropt; */ /* real optional input, output     */

	/* Error File */

	/*cv_errfp; */      /* CVODE error messages are sent to errfp */

	/* Pointer to Machine Environment-Specific Information */

	/* cv_machenv */;
}
