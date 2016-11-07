/******************************************************************
 *                                                                *
 * File          : cvode.c                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the main CVODE integrator. *
 * It is independent of the CVODE linear solver in use.           *
 *                                                                *
 ******************************************************************/


/************************************************************/
/******************* BEGIN Imports **************************/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "cvode.h"
#include "cvnltyps.h"
#include "cvvector.h"
#include "cvnlmath.h"

/************************************************************/
/******************** END Imports ***************************/
/************************************************************/


/***************************************************************/
/*********************** BEGIN Macros **************************/
/***************************************************************/

/* Macro: loop */

#define loop for(;;)

/***************************************************************/
/************************ END Macros ***************************/
/***************************************************************/



/************************************************************/
/************** BEGIN CVODE Private Constants ***************/
/************************************************************/

#define HALF   RCONST(0.5)  /* real 0.5   */
#define ZERO   RCONST(0.0)  /* real 0.0   */
#define ONE    RCONST(1.0)  /* real 1.0   */
#define TWO    RCONST(2.0)  /* real 2.0   */
#define TWELVE RCONST(12.0) /* real 12.0  */

/***************************************************************/
/************** BEGIN Default Constants ************************/
/***************************************************************/

#define HMIN_DEFAULT     ZERO    /* hmin default value     */
#define HMAX_INV_DEFAULT ZERO    /* hmax_inv default value */
#define MXHNIL_DEFAULT   10      /* mxhnil default value   */
#define MXSTEP_DEFAULT   500     /* mxstep default value   */


/***************************************************************/
/*************** END Default Constants *************************/
/***************************************************************/


/***************************************************************/
/************ BEGIN Routine-Specific Constants *****************/
/***************************************************************/

/* CVodeDky */

#define FUZZ_FACTOR RCONST(100.0)

/* CVHin */

#define HLB_FACTOR RCONST(100.0)
#define HUB_FACTOR RCONST(0.1)
#define H_BIAS     HALF
#define MAX_ITERS  4

/* CVSet */

#define CORTES RCONST(0.1)

/* CVStep return values */

#define SUCCESS_STEP      0
#define REP_ERR_FAIL     -1
#define REP_CONV_FAIL    -2
#define SETUP_FAILED     -3
#define SOLVE_FAILED     -4

/* CVStep control constants */

#define PREDICT_AGAIN    -5
#define DO_ERROR_TEST     1

/* CVStep */

#define THRESH RCONST(1.5)
#define ETAMX1 RCONST(10000.0) 
#define ETAMX2 RCONST(10.0)
#define ETAMX3 RCONST(10.0)
#define ETAMXF RCONST(0.2)
#define ETAMIN RCONST(0.1)
#define ETACF  RCONST(0.25)
#define ADDON  RCONST(0.000001)
#define BIAS1  RCONST(6.0)
#define BIAS2  RCONST(6.0)
#define BIAS3  RCONST(10.0)
#define ONEPSM RCONST(1.000001)

#define SMALL_NST    10   /* nst > SMALL_NST => use ETAMX3          */
#define MXNCF        10   /* max no. of convergence failures during */
		          /* one step try                           */
#define MXNEF         7   /* max no. of error test failures during  */
		          /* one step try                           */
#define MXNEF1        3   /* max no. of error test failures before  */
		          /* forcing a reduction of order           */
#define SMALL_NEF     2   /* if an error failure occurs and         */
                          /* SMALL_NEF <= nef <= MXNEF1, then       */
                          /* reset eta =  MIN(eta, ETAMXF)          */
#define LONG_WAIT    10   /* number of steps to wait before         */
                          /* considering an order change when       */
                          /* q==1 and MXNEF1 error test failures    */
                          /* have occurred                          */

/* CVnls return values */

#define SOLVED            0
#define CONV_FAIL        -1 
#define SETUP_FAIL_UNREC -2
#define SOLVE_FAIL_UNREC -3

/* CVnls input flags */

#define FIRST_CALL      0
#define PREV_CONV_FAIL -1
#define PREV_ERR_FAIL  -2

/* CVnls other constants */

#define FUNC_MAXCOR 3  /* maximum no. of corrector iterations   */
                       /* for iter == FUNCTIONAL                */
#define NEWT_MAXCOR 3  /* maximum no. of corrector iterations   */
                       /* for iter == NEWTON                    */

#define CRDOWN RCONST(0.3) /* constant used in the estimation of the   */
                           /* convergence rate (crate) of the          */
                           /* iterates for the nonlinear equation      */
#define DGMAX  RCONST(0.3) /* iter == NEWTON, |gamma/gammap-1| > DGMAX */
			   /* => call lsetup                           */

#define RDIV      TWO  /* declare divergence if ratio del/delp > RDIV  */
#define MSBP       20  /* max no. of steps between lsetup calls        */

#define TRY_AGAIN  99  /* control constant for CVnlsNewton - should be */
		       /* distinct from CVnls return values            */


/***************************************************************/
/*************** END Routine-Specific Constants  ***************/
/***************************************************************/


/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* CVodeMalloc Error Messages */

#define CVM             "CVodeMalloc-- "

#define MSG_Y0_NULL     CVM "y0=NULL illegal.\n\n"

#define MSG_BAD_N       CVM "N=%ld < 1 illegal.\n\n"

#define MSG_BAD_LMM_1   CVM "lmm=%d illegal.\n"
#define MSG_BAD_LMM_2   "The legal values are ADAMS=%d and BDF=%d.\n\n"
#define MSG_BAD_LMM     MSG_BAD_LMM_1 MSG_BAD_LMM_2

#define MSG_BAD_ITER_1  CVM "iter=%d illegal.\n"
#define MSG_BAD_ITER_2  "The legal values are FUNCTIONAL=%d "
#define MSG_BAD_ITER_3  "and NEWTON=%d.\n\n"
#define MSG_BAD_ITER    MSG_BAD_ITER_1 MSG_BAD_ITER_2 MSG_BAD_ITER_3

#define MSG_BAD_ITOL_1  CVM "itol=%d illegal.\n"
#define MSG_BAD_ITOL_2  "The legal values are SS=%d and SV=%d.\n\n"
#define MSG_BAD_ITOL    MSG_BAD_ITOL_1 MSG_BAD_ITOL_2

#define MSG_F_NULL       CVM "f=NULL illegal.\n\n"

#define MSG_RELTOL_NULL  CVM "reltol=NULL illegal.\n\n"
 
#define MSG_BAD_RELTOL   CVM "*reltol=%g < 0 illegal.\n\n"

#define MSG_ABSTOL_NULL  CVM "abstol=NULL illegal.\n\n"

#define MSG_BAD_ABSTOL   CVM "Some abstol component < 0.0 illegal.\n\n"

#define MSG_BAD_OPTIN_1  CVM "optIn=%d illegal.\n"
#define MSG_BAD_OPTIN_2  "The legal values are FALSE=%d and TRUE=%d.\n\n"
#define MSG_BAD_OPTIN    MSG_BAD_OPTIN_1 MSG_BAD_OPTIN_2

#define MSG_BAD_OPT     CVM "optIn=TRUE, but iopt=ropt=NULL.\n\n"

#define MSG_BAD_HMIN_HMAX_1 CVM "Inconsistent step size limits:\n"
#define MSG_BAD_HMIN_HMAX_2 "ropt[HMIN]=%g > ropt[HMAX]=%g.\n\n"
#define MSG_BAD_HMIN_HMAX   MSG_BAD_HMIN_HMAX_1 MSG_BAD_HMIN_HMAX_2

#define MSG_MEM_FAIL    CVM "A memory request failed.\n\n"

#define MSG_BAD_EWT     CVM "Some initial ewt component = 0.0 illegal.\n\n"


/* CVode error messages */

#define CVODE            "CVode-- "

#define NO_MEM           "cvode_mem=NULL illegal.\n\n"

#define MSG_CVODE_NO_MEM CVODE NO_MEM
 
#define MSG_LINIT_NULL   CVODE "The linear solver's init routine is NULL.\n\n"

#define MSG_LSETUP_NULL  CVODE "The linear solver's setup routine is NULL.\n\n"

#define MSG_LSOLVE_NULL  CVODE "The linear solver's solve routine is NULL.\n\n"

#define MSG_LFREE_NULL   CVODE "The linear solver's free routine is NULL.\n\n"

#define MSG_LINIT_FAIL   CVODE "The linear solver's init routine failed.\n\n"

#define MSG_YOUT_NULL    CVODE "yout=NULL illegal.\n\n"

#define MSG_T_NULL       CVODE "t=NULL illegal.\n\n"

#define MSG_BAD_ITASK_1   CVODE "itask=%d illegal.\nThe legal values are"
#define MSG_BAD_ITASK_2   " NORMAL=%d and ONE_STEP=%d.\n\n"
#define MSG_BAD_ITASK     MSG_BAD_ITASK_1 MSG_BAD_ITASK_2

#define MSG_BAD_H0        CVODE "h0=%g and tout-t0=%g inconsistent.\n\n"

#define MSG_BAD_TOUT_1    CVODE "Trouble interpolating at tout = %g.\n"
#define MSG_BAD_TOUT_2    "tout too far back in direction of integration.\n\n"
#define MSG_BAD_TOUT      MSG_BAD_TOUT_1 MSG_BAD_TOUT_2

#define MSG_MAX_STEPS_1   CVODE "At t=%g, mxstep=%d steps taken on "
#define MSG_MAX_STEPS_2   "this call before\nreaching tout=%g.\n\n"
#define MSG_MAX_STEPS     MSG_MAX_STEPS_1 MSG_MAX_STEPS_2

#define MSG_EWT_NOW_BAD_1  CVODE "At t=%g, "
#define MSG_EWT_NOW_BAD_2  "some ewt component has become <= 0.0.\n\n"
#define MSG_EWT_NOW_BAD    MSG_EWT_NOW_BAD_1 MSG_EWT_NOW_BAD_2

#define MSG_TOO_MUCH_ACC  CVODE "At t=%g, too much accuracy requested.\n\n"

#define MSG_HNIL_1  CVODE "Warning.. internal t=%g and step size h=%g\n"
#define MSG_HNIL_2  "are such that t + h == t on the next step.\n"
#define MSG_HNIL_3  "The solver will continue anyway.\n\n"
#define MSG_HNIL    MSG_HNIL_1 MSG_HNIL_2 MSG_HNIL_3

#define MSG_HNIL_DONE_1   CVODE "The above warning has been issued %d times "
#define MSG_HNIL_DONE_2   "and will not be\nissued again for this problem.\n\n"
#define MSG_HNIL_DONE     MSG_HNIL_DONE_1 MSG_HNIL_DONE_2

#define MSG_ERR_FAILS_1   CVODE "At t=%g and step size h=%g, the error test\n"
#define MSG_ERR_FAILS_2   "failed repeatedly or with |h| = hmin.\n\n"
#define MSG_ERR_FAILS     MSG_ERR_FAILS_1 MSG_ERR_FAILS_2

#define MSG_CONV_FAILS_1  CVODE "At t=%g and step size h=%g, the corrector\n"
#define MSG_CONV_FAILS_2  "convergence failed repeatedly or "
#define MSG_CONV_FAILS_3  "with |h| = hmin.\n\n"
#define MSG_CONV_FAILS    MSG_CONV_FAILS_1 MSG_CONV_FAILS_2 MSG_CONV_FAILS_3

#define MSG_SETUP_FAILED_1 CVODE "At t=%g, the setup routine failed in an "
#define MSG_SETUP_FAILED_2 "unrecoverable manner.\n\n"
#define MSG_SETUP_FAILED   MSG_SETUP_FAILED_1 MSG_SETUP_FAILED_2

#define MSG_SOLVE_FAILED_1 CVODE "At t=%g, the solve routine failed in an "
#define MSG_SOLVE_FAILED_2 "unrecoverable manner.\n\n"
#define MSG_SOLVE_FAILED   MSG_SOLVE_FAILED_1 MSG_SOLVE_FAILED_2

#define MSG_TOO_CLOSE_1    CVODE "tout=%g too close to t0=%g to start"
#define MSG_TOO_CLOSE_2    " integration.\n\n"
#define MSG_TOO_CLOSE      MSG_TOO_CLOSE_1 MSG_TOO_CLOSE_2


/* CVodeDky Error Messages */

#define DKY         "CVodeDky-- "

#define MSG_DKY_NO_MEM  DKY NO_MEM

#define MSG_BAD_K   DKY "k=%d illegal.\n\n"

#define MSG_BAD_T_1 DKY "t=%g illegal.\n"
#define MSG_BAD_T_2 "t not in interval tcur-hu=%g to tcur=%g.\n\n"
#define MSG_BAD_T   MSG_BAD_T_1 MSG_BAD_T_2

#define MSG_BAD_DKY DKY "dky=NULL illegal.\n\n"

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/


/************************************************************/
/*************** END CVODE Private Constants ****************/
/************************************************************/


/**************************************************************/
/********* BEGIN Private Helper Functions Prototypes **********/
/**************************************************************/

static bool CVAllocVectors(CVodeMemData* cv_mem, integer neq, int maxord,
			   void *machEnv);
static void CVFreeVectors(CVodeMemData* cv_mem, int maxord);

static bool CVEwtSet(CVodeMemData* cv_mem, real *rtol, void *atol, int tol_type,
		     iN_Vector* ycur, iN_Vector* ewtvec, integer neq);
static bool CVEwtSetSS(CVodeMemData* cv_mem, real *rtol, real *atol,
		       iN_Vector* ycur, iN_Vector* ewtvec, integer neq);
static bool CVEwtSetSV(CVodeMemData* cv_mem, real *rtol, iN_Vector* atol,
		       iN_Vector* ycur, iN_Vector* ewtvec, integer neq);

static bool CVHin(CVodeMemData* cv_mem, real tout);
static real CVUpperBoundH0(CVodeMemData* cv_mem, real tdist);
static real CVYddNorm(CVodeMemData* cv_mem, real hg);

static int  CVStep(CVodeMemData* cv_mem);

static void CVAdjustParams(CVodeMemData* cv_mem);
static void CVAdjustOrder(CVodeMemData* cv_mem, int deltaq);
static void CVAdjustAdams(CVodeMemData* cv_mem, int deltaq);
static void CVAdjustBDF(CVodeMemData* cv_mem, int deltaq);
static void CVIncreaseBDF(CVodeMemData* cv_mem);
static void CVDecreaseBDF(CVodeMemData* cv_mem);

static void CVRescale(CVodeMemData* cv_mem);

static void CVPredict(CVodeMemData* cv_mem);

static void CVSet(CVodeMemData* cv_mem);
static void CVSetAdams(CVodeMemData* cv_mem);
static real CVAdamsStart(CVodeMemData* cv_mem, real m[]);
static void CVAdamsFinish(CVodeMemData* cv_mem, real m[], real M[], real hsum);
static real CVAltSum(int iend, real a[], int k);
static void CVSetBDF(CVodeMemData* cv_mem);
static void CVSetTqBDF(CVodeMemData* cv_mem, real hsum, real alpha0,
		       real alpha0_hat, real xi_inv, real xistar_inv);

static int CVnls(CVodeMemData* cv_mem, int nflag);
static int CVnlsFunctional(CVodeMemData* cv_mem);
static int CVnlsNewton(CVodeMemData* cv_mem, int nflag);
static int CVNewtonIteration(CVodeMemData* cv_mem);

static int  CVHandleNFlag(CVodeMemData* cv_mem, int *nflagPtr, real saved_t,
			  int *ncfPtr);

static void CVRestore(CVodeMemData* cv_mem, real saved_t);

static bool CVDoErrorTest(CVodeMemData* cv_mem, int *nflagPtr, int *kflagPtr,
			  real saved_t, int *nefPtr, real *dsmPtr);

static void CVCompleteStep(CVodeMemData* cv_mem);

static void CVPrepareNextStep(CVodeMemData* cv_mem, real dsm);
static void CVSetEta(CVodeMemData* cv_mem);
static real CVComputeEtaqm1(CVodeMemData* cv_mem);
static real CVComputeEtaqp1(CVodeMemData* cv_mem);
static void CVChooseEta(CVodeMemData* cv_mem,real etaqm1, real etaq, real etaqp1);

static int  CVHandleFailure(CVodeMemData* cv_mem,int kflag);


/**************************************************************/
/********** END Private Helper Functions Prototypes ***********/
/**************************************************************/


/**************************************************************/
/**************** BEGIN Readability Constants *****************/
/**************************************************************/


#define cv_mem_uround (cv_mem->cv_uround)  
#define cv_mem_zn     (cv_mem->cv_zn) 
#define cv_mem_ewt    (cv_mem->cv_ewt)  
#define cv_mem_y      (cv_mem->cv_y)
#define cv_mem_acor   (cv_mem->cv_acor)
#define cv_mem_tempv  (cv_mem->cv_tempv)
#define cv_mem_ftemp  (cv_mem->cv_ftemp) 
#define cv_mem_q      (cv_mem->cv_q)
#define cv_mem_qprime (cv_mem->cv_qprime)
#define cv_mem_qwait  (cv_mem->cv_qwait)
#define cv_mem_L      (cv_mem->cv_L)
#define cv_mem_h      (cv_mem->cv_h)
#define cv_mem_hprime (cv_mem->cv_hprime)
#define cv_mem_eta    (cv_mem-> cv_eta) 
#define cv_mem_hscale (cv_mem->cv_hscale) 
#define cv_mem_tn     (cv_mem->cv_tn)
#define cv_mem_tau    (cv_mem->cv_tau)
#define cv_mem_tq     (cv_mem->cv_tq)
#define cv_mem_l      (cv_mem->cv_l)
#define cv_mem_rl1    (cv_mem->cv_rl1)
#define cv_mem_gamma  (cv_mem->cv_gamma) 
#define cv_mem_gammap (cv_mem->cv_gammap) 
#define cv_mem_gamrat (cv_mem->cv_gamrat)
#define cv_mem_crate  (cv_mem->cv_crate)
#define cv_mem_acnrm  (cv_mem->cv_acnrm)
#define cv_mem_mnewt  (cv_mem->cv_mnewt)
#define cv_mem_qmax   (cv_mem->cv_qmax) 
#define cv_mem_mxstep (cv_mem->cv_mxstep)
#define cv_mem_maxcor (cv_mem->cv_maxcor)
#define cv_mem_mxhnil (cv_mem->cv_mxhnil)
#define cv_mem_hmin   (cv_mem->cv_hmin)
#define cv_mem_hmax_inv (cv_mem->cv_hmax_inv)
#define cv_mem_etamax (cv_mem->cv_etamax)
#define cv_mem_nst    (cv_mem->cv_nst)
#define cv_mem_nfe    (cv_mem->cv_nfe)
#define cv_mem_ncfn   (cv_mem->cv_ncfn)
#define cv_mem_netf   (cv_mem->cv_netf)
#define cv_mem_nni    (cv_mem-> cv_nni)
#define cv_mem_nsetups (cv_mem->cv_nsetups)
#define cv_mem_nhnil  (cv_mem->cv_nhnil)
#define cv_mem_lrw    (cv_mem->cv_lrw)
#define cv_mem_liw    (cv_mem->cv_liw)
#define cv_mem_linit  (cv_mem->cv_linit)
#define cv_mem_lsetup (cv_mem->cv_lsetup)
#define cv_mem_lsolve (cv_mem->cv_lsolve) 
#define cv_mem_lfree  (cv_mem->cv_lfree) 
#define cv_mem_lmem   (cv_mem->cv_lmem) 
#define cv_mem_linitOK (cv_mem->cv_linitOK)
#define cv_mem_qu     (cv_mem->cv_qu)          
#define cv_mem_nstlp  (cv_mem->cv_nstlp)  
#define cv_mem_hu     (cv_mem->cv_hu)         
#define cv_mem_saved_tq5 (cv_mem->cv_saved_tq5)  
#define cv_mem_jcur   (cv_mem->cv_jcur)         
#define cv_mem_tolsf  (cv_mem->cv_tolsf)      
#define cv_mem_setupNonNull (cv_mem->cv_setupNonNull) 
#define cv_mem_machenv (cv_mem->cv_machenv)

/**************************************************************/
/***************** END Readability Constants ******************/
/**************************************************************/


/***************************************************************/
/************* BEGIN CVODE Implementation **********************/
/***************************************************************/


/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/


/******************** CVodeMalloc *******************************

 CVode Malloc allocates and initializes memory for a problem. All
 problem specification inputs are checked for errors. If any
 error occurs during initialization, it is reported to the file
 whose file pointer is errfp and NULL is returned. Otherwise, the
 pointer to successfully initialized problem memory is returned.
 
*****************************************************************/

void *CVodeMalloc(integer N, RhsFn f, real t0, iN_Vector* y0, int lmm, int iter,
		  int itol, real *reltol, void *abstol, void *f_data,
		  FILE *errfp, bool optIn, long int iopt[], real ropt[],
		  void *machEnv)
{
  bool    allocOK, ioptExists, roptExists, neg_abstol, ewtsetOK;
  int     maxord;
  CVodeMemData* cv_mem;
  FILE *fp;
  
  /* Check for legal input parameters */
  
  fp = (errfp == NULL) ? stdout : errfp;

  if (y0==NULL) {
//    fprintf(fp, MSG_Y0_NULL);
    return(NULL);
  }
  
  if (N <= 0) {
//    fprintf(fp, MSG_BAD_N, N);
    return(NULL);
  }

  if ((lmm != ADAMS) && (lmm != BDF)) {
//    fprintf(fp, MSG_BAD_LMM, lmm, ADAMS, BDF);
    return(NULL);
  }

  if ((iter != FUNCTIONAL) && (iter != NEWTON)) {
//    fprintf(fp, MSG_BAD_ITER, iter, FUNCTIONAL, NEWTON);
    return(NULL);
  }

  if ((itol != SS) && (itol != SV)) {
//    fprintf(fp, MSG_BAD_ITOL, itol, SS, SV);
    return(NULL);
  }

  if (f == NULL) {
//    fprintf(fp, MSG_F_NULL);
    return(NULL);
  }

  if (reltol == NULL) {
//    fprintf(fp, MSG_RELTOL_NULL);
    return(NULL);
  }

  if (*reltol < ZERO) {
//    fprintf(fp, MSG_BAD_RELTOL, *reltol);
    return(NULL);
  }
   
  if (abstol == NULL) {
//    fprintf(fp, MSG_ABSTOL_NULL);
    return(NULL);
  }

  if (itol == SS) {
    neg_abstol = (*((real *)abstol) < ZERO);
  } else {
    neg_abstol = (N_VMin((iN_Vector*)abstol) < ZERO);
  }
  if (neg_abstol) {
//    fprintf(fp, MSG_BAD_ABSTOL);
    return(NULL);
  }

  if ((optIn != FALSE) && (optIn != TRUE)) {
//    fprintf(fp, MSG_BAD_OPTIN, optIn, FALSE, TRUE);
    return(NULL);
  }

  if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
//    fprintf(fp, MSG_BAD_OPT);
    return(NULL);
  } 

  ioptExists = (iopt != NULL);
  roptExists = (ropt != NULL);

  if (optIn && roptExists) {
    if ((ropt[HMAX] > ZERO) && (ropt[HMIN] > ropt[HMAX])) {
//      fprintf(fp, MSG_BAD_HMIN_HMAX, ropt[HMIN], ropt[HMAX]);
      return(NULL);
    }
  }

  /* compute maxord */

  maxord = (lmm == ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

  if (optIn && ioptExists) {
    if (iopt[MAXORD] > 0) 
		maxord = MIN(maxord, iopt[MAXORD]);
  }

  cv_mem = (CVodeMemData*) malloc(sizeof(struct CVodeMemRec));
  if (cv_mem == NULL) {
//    fprintf(fp, MSG_MEM_FAIL);
    return(NULL);
  }

  /* initialize some on zero, is good for BoundsChecker */
  memset(cv_mem, 0, sizeof(struct CVodeMemRec));
 
  /* Allocate the vectors */

  allocOK = CVAllocVectors(cv_mem, N, maxord, machEnv);
  if (!allocOK) {
//    fprintf(fp, MSG_MEM_FAIL);
    free(cv_mem);
    return(NULL);
  }
 
  /* Set the ewt vector */

  ewtsetOK = CVEwtSet(cv_mem, reltol, abstol, itol, y0, cv_mem_ewt, N);
  if (!ewtsetOK) {
//    fprintf(fp, MSG_BAD_EWT);
    CVFreeVectors(cv_mem, maxord);
    free(cv_mem);
    return(NULL);
  }
  
  /* All error checking is complete at this point */
  
  /* Copy the input parameters into CVODE state */

  cv_mem->cv_N = N;  /* readability constants defined below CVodeMalloc */
  cv_mem->cv_f = f;
  cv_mem->cv_f_data = f_data;
  cv_mem->cv_lmm = lmm;    
  cv_mem->cv_iter = iter;
  cv_mem->cv_itol = itol;
  cv_mem->cv_reltol = reltol;      
  cv_mem->cv_abstol = abstol;
  cv_mem->cv_iopt = iopt;
  cv_mem->cv_ropt = ropt;
  cv_mem->cv_errfp = fp;
  cv_mem_tn = t0;
  cv_mem_machenv = machEnv;

  /* Set step parameters */

  cv_mem_q = 1;
  cv_mem_L = 2;
  cv_mem_qwait = cv_mem_L;
  cv_mem_qmax = maxord;
  cv_mem_etamax = ETAMX1;

  /* Set uround */

  cv_mem_uround = UnitRoundoff();

  /* Set the linear solver addresses to NULL, linitOK to FALSE */

  cv_mem_linit = NULL;
  cv_mem_lsetup = NULL;
  cv_mem_lsolve = NULL;
  cv_mem_lfree = NULL;
  cv_mem_lmem = NULL;
  /* We check != NULL later, in CVode and linit, if using NEWTON */
  cv_mem_linitOK = FALSE;

  /* Initialize the history array zn */
  
  N_VScale(ONE, y0, cv_mem_zn[0]); 
  f(N, t0, y0, cv_mem_zn[1], f_data); 
  cv_mem_nfe = 1;
 
  /* Handle the remaining optional inputs */

  cv_mem_hmin = HMIN_DEFAULT;
  cv_mem_hmax_inv = HMAX_INV_DEFAULT;
  if (optIn && roptExists) {
    if (ropt[HMIN] > ZERO)
		cv_mem_hmin = ropt[HMIN];
    if (ropt[HMAX] > ZERO)
		cv_mem_hmax_inv = ONE/ropt[HMAX];
  }

  cv_mem_mxhnil = MXHNIL_DEFAULT;
  cv_mem_mxstep = MXSTEP_DEFAULT;
  if (optIn && ioptExists) {
    if (iopt[MXHNIL] > 0)
		cv_mem_mxhnil = iopt[MXHNIL];
    if (iopt[MXSTEP] > 0)
		cv_mem_mxstep = iopt[MXSTEP];
  }
 
  if ((!optIn) && roptExists) ropt[H0] = ZERO;
 
  /* Set maxcor */

  cv_mem_maxcor = (iter==NEWTON) ? NEWT_MAXCOR : FUNC_MAXCOR;
  
  /* Initialize all the counters */
 
  cv_mem_nst = cv_mem_ncfn = cv_mem_netf = cv_mem_nni = cv_mem_nsetups = cv_mem_nhnil = cv_mem_nstlp = 0;
  
  /* Initialize all other vars corresponding to optional outputs */
  
  cv_mem_qu = 0;
  cv_mem_hu = ZERO;
  cv_mem_tolsf = ONE;

  /* Initialize optional output locations in iopt, ropt */

  if (ioptExists) {
    iopt[NST] = iopt[NFE] = iopt[NSETUPS] = iopt[NNI] = 0;
    iopt[NCFN] = iopt[NETF] = 0;
    iopt[QU] = cv_mem_qu;
    iopt[QCUR] = 0;
    iopt[LENRW] = cv_mem_lrw;
    iopt[LENIW] = cv_mem_liw;
  }
  
  if (roptExists) {
    ropt[HU] = cv_mem_hu;
    ropt[HCUR] = ZERO;
    ropt[TCUR] = t0;
    ropt[TOLSF] = cv_mem_tolsf;
  }
      
  /* Problem has been successfully initialized */
  return((void *)cv_mem);
}


/**************************************************************/
/************** BEGIN More Readability Constants **************/
/**************************************************************/

#define cv_mem_N      (cv_mem->cv_N)
#define cv_mem_f      (cv_mem->cv_f)      
#define cv_mem_f_data (cv_mem->cv_f_data)    
#define cv_mem_lmm    (cv_mem->cv_lmm) 
#define cv_mem_iter   (cv_mem->cv_iter)        
#define cv_mem_itol   (cv_mem->cv_itol)         
#define cv_mem_reltol (cv_mem->cv_reltol)       
#define cv_mem_abstol (cv_mem->cv_abstol)     
#define cv_mem_iopt   (cv_mem->cv_iopt)
#define cv_mem_ropt   (cv_mem->cv_ropt)
#define cv_mem_errfp  (cv_mem->cv_errfp)

/**************************************************************/
/*************** END More Readability Constants ***************/
/**************************************************************/


/********************* CVode ****************************************

 This routine is the main driver of the CVODE package. 

 It integrates over a time interval defined by the user, by calling
 CVStep to do internal time steps.

 The first time that CVode is called for a successfully initialized
 problem, it computes a tentative initial step size h.

 CVode supports two modes, specified by itask: NORMAL and ONE_STEP.
 In the NORMAL mode, the solver steps until it reaches or passes tout
 and then interpolates to obtain y(tout).
 In the ONE_STEP mode, it takes one internal step and returns.

********************************************************************/

int CVode(void *cvode_mem, real tout, iN_Vector* yout, real *t, int itask)
{
  int nstloc, kflag, istate, next_q, ier;
  real rh, next_h;
  bool hOK, ewtsetOK;
  CVodeMemData* cv_mem;

  /* Check for legal inputs in all cases */

  cv_mem = (CVodeMemData*) cvode_mem;
  if (cvode_mem == NULL) {
//    fprintf(stdout, MSG_CVODE_NO_MEM);
    return(CVODE_NO_MEM);
  }
  
  if ((cv_mem_y = yout) == NULL) {
//    fprintf(errfp, MSG_YOUT_NULL);       
    return(ILL_INPUT);
  }
  
  if (t == NULL) {
//    fprintf(errfp, MSG_T_NULL);
    return(ILL_INPUT);
  }
  *t = cv_mem_tn;

  if ((itask != NORMAL) && (itask != ONE_STEP)) {
//    fprintf(errfp, MSG_BAD_ITASK, itask, NORMAL, ONE_STEP);
    return(ILL_INPUT);
  }

  /* On first call, check solver functions and call linit function */
  
  if (cv_mem_nst == 0) {
    if (cv_mem_iter == NEWTON) {
      if (cv_mem_linit == NULL) {
//	fprintf(errfp, MSG_LINIT_NULL);
	return(ILL_INPUT);
      }
      if (cv_mem_lsetup == NULL) {
//	fprintf(errfp, MSG_LSETUP_NULL);
	return(ILL_INPUT);
      }
      if (cv_mem_lsolve == NULL) {
//	fprintf(errfp, MSG_LSOLVE_NULL);
	return(ILL_INPUT);
      }
      if (cv_mem_lfree == NULL) {
//	fprintf(errfp, MSG_LFREE_NULL);
	return(ILL_INPUT);
      }
      cv_mem_linitOK = (cv_mem_linit(cv_mem, &(cv_mem_setupNonNull)) == LINIT_OK);
      if (!cv_mem_linitOK) {
//	fprintf(errfp, MSG_LINIT_FAIL);
	return(ILL_INPUT);
      }
    }

    /* On first call, set initial h (from H0 or CVHin) and scale zn[1] */
    
    cv_mem_h = ZERO;
    if (cv_mem_ropt != NULL) cv_mem_h = cv_mem_ropt[H0];
    if ( (cv_mem_h != ZERO) && ((tout-cv_mem_tn)*cv_mem_h < ZERO) ) {
//      fprintf(errfp, MSG_BAD_H0, h, tout-tn);
      return(ILL_INPUT);
    }
    if (cv_mem_h == ZERO) {
      hOK = CVHin(cv_mem, tout);
      if (!hOK) {
//	fprintf(errfp, MSG_TOO_CLOSE, tout, tn);
	return(ILL_INPUT);
      }
    }
    rh = ABS(cv_mem_h)*cv_mem_hmax_inv;
    if (rh > ONE) cv_mem_h /= rh;
    if (ABS(cv_mem_h) < cv_mem_hmin) cv_mem_h *= cv_mem_hmin/ABS(cv_mem_h);
    cv_mem_hscale = cv_mem_h; 
    N_VScale(cv_mem_h, cv_mem_zn[1], cv_mem_zn[1]);
  }

  /* If not the first call, check if tout already reached */

  if ( (itask == NORMAL) && (cv_mem_nst > 0) && ((cv_mem_tn-tout)*cv_mem_h >= ZERO) ) {
    *t = tout;
    ier =  CVodeDky(cv_mem, tout, 0, yout);
    if (ier != OKAY) {  /* ier must be == BAD_T */
//      fprintf(errfp, MSG_BAD_TOUT, tout);
      return(ILL_INPUT);
    }
    return(SUCCESS);
  }

  /* Looping point for internal steps */

  nstloc = 0;
  loop {
   
    next_h = cv_mem_h;
    next_q = cv_mem_q;
    
    /* Reset and check ewt */

    if (cv_mem_nst > 0) {
      ewtsetOK = CVEwtSet(cv_mem, cv_mem_reltol, cv_mem_abstol, cv_mem_itol, cv_mem_zn[0], cv_mem_ewt, cv_mem_N);
      if (!ewtsetOK) {
//	fprintf(errfp, MSG_EWT_NOW_BAD, tn);
	istate = ILL_INPUT;
	*t = cv_mem_tn;
	N_VScale(ONE, cv_mem_zn[0], yout);
	break;
      }
    }

    /* Check for too many steps */
    
    if (nstloc >= cv_mem_mxstep) {
//      fprintf(errfp, MSG_MAX_STEPS, tn, mxstep, tout);
      istate = TOO_MUCH_WORK;
      *t = cv_mem_tn;
      N_VScale(ONE, cv_mem_zn[0], yout);
      break;
    }

    /* Check for too much accuracy requested */

    if ((cv_mem_tolsf = cv_mem_uround * N_VWrmsNorm(cv_mem_zn[0], cv_mem_ewt)) > ONE) {
//      fprintf(errfp, MSG_TOO_MUCH_ACC, tn);
      istate = TOO_MUCH_ACC;
      *t = cv_mem_tn;
      N_VScale(ONE, cv_mem_zn[0], yout);
      cv_mem_tolsf *= TWO;
      break;
    }

    /* Check for h below roundoff level in tn */

    if (cv_mem_tn + cv_mem_h == cv_mem_tn) {
      cv_mem_nhnil++;
      if (cv_mem_nhnil <= cv_mem_mxhnil) fprintf(cv_mem_errfp, MSG_HNIL, cv_mem_tn, cv_mem_h);
      if (cv_mem_nhnil == cv_mem_mxhnil) fprintf(cv_mem_errfp, MSG_HNIL_DONE, cv_mem_mxhnil);
    }

    /* Call CVStep to take a step */

    kflag = CVStep(cv_mem);

    /* Process failed step cases, and exit loop */
   
    if (kflag != SUCCESS_STEP) {
      istate = CVHandleFailure(cv_mem, kflag);
      *t = cv_mem_tn;
      N_VScale(ONE, cv_mem_zn[0], yout);
      break;
    }
    
    nstloc++;

    /* Check if in one-step mode, and if so copy y and exit loop */
    
    if (itask == ONE_STEP) {
      istate = SUCCESS;


	  /* added by grn 2009-04-16 for V4.1 */
#if 1
		if ((cv_mem_tn-tout)*cv_mem_h >= ZERO) {
			*t = tout;
			(void) CVodeDky(cv_mem, tout, 0, yout);
		}
		else
		{
			/* the original code */
			*t = cv_mem_tn;
			N_VScale(ONE, cv_mem_zn[0], yout);
		}
#else
		*t = tn;
		N_VScale(ONE, zn[0], yout);
#endif

	  next_q = cv_mem_qprime;
      next_h = cv_mem_hprime;

	  break;
    }

    /* Check if tout reached, and if so interpolate and exit loop */

    if ((cv_mem_tn-tout)*cv_mem_h >= ZERO) {
      istate = SUCCESS;
      *t = tout;
      (void) CVodeDky(cv_mem, tout, 0, yout);
      next_q = cv_mem_qprime;
      next_h = cv_mem_hprime;
      break;
    }
  }

  /* End of step loop; load optional outputs and return */

  if (cv_mem_iopt != NULL) {
    cv_mem_iopt[NST] = cv_mem_nst;
    cv_mem_iopt[NFE] = cv_mem_nfe;
    cv_mem_iopt[NSETUPS] = cv_mem_nsetups;
    cv_mem_iopt[NNI] = cv_mem_nni;
    cv_mem_iopt[NCFN] = cv_mem_ncfn;
    cv_mem_iopt[NETF] = cv_mem_netf;
    cv_mem_iopt[QU] = cv_mem_q;
    cv_mem_iopt[QCUR] = next_q;
  }
  
  if (cv_mem_ropt != NULL) {
    cv_mem_ropt[HU] = cv_mem_h;
    cv_mem_ropt[HCUR] = next_h;
    cv_mem_ropt[TCUR] = cv_mem_tn;
    cv_mem_ropt[TOLSF] = cv_mem_tolsf;
  }
  
  return(istate);	
}

/*************** CVodeDky ********************************************

 This routine computes the k-th derivative of the interpolating
 polynomial at the time t and stores the result in the vector dky.
 The formula is:
          q 
   dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] , 
         j=k 
 where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 zn[j] is the j-th column of the Nordsieck history array.

 This function is called by CVode with k = 0 and t = tout, but
 may also be called directly by the user.

**********************************************************************/

int CVodeDky(void *cvode_mem, real t, int k, iN_Vector* dky)
{
  real s, c, r;
  real tfuzz, tp, tn1;
  int i, j;
  CVodeMemData* cv_mem;
  
  cv_mem = (CVodeMemData*) cvode_mem;

  /* Check all inputs for legality */
 
  if (cvode_mem == NULL) {
//    fprintf(stdout, MSG_DKY_NO_MEM);
    return(DKY_NO_MEM);
  }
  
  if (dky == NULL) {
//    fprintf(stdout, MSG_BAD_DKY);
    return(BAD_DKY);
  }

  if ((k < 0) || (k > cv_mem_q)) {
//    fprintf(errfp, MSG_BAD_K, k);
    return(BAD_K);
  }
  
  tfuzz = FUZZ_FACTOR * cv_mem_uround * (cv_mem_tn + cv_mem_hu);
  tp = cv_mem_tn - cv_mem_hu - tfuzz;
  tn1 = cv_mem_tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
//    fprintf(errfp, MSG_BAD_T, t, tn-hu, tn);
    return(BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */

  s = (t - cv_mem_tn) / cv_mem_h;
  for (j=cv_mem_q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == cv_mem_q) {
      N_VScale(c, cv_mem_zn[cv_mem_q], dky);
    } else {
      N_VLinearSum(c, cv_mem_zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(OKAY);
  r = RPowerI(cv_mem_h,-k);
  N_VScale(r, dky, dky);
  return(OKAY);
}
 
/********************* CVodeFree **********************************

 This routine frees the problem memory allocated by CVodeMalloc.
 Such memory includes all the vectors allocated by CVAllocVectors,
 and the memory lmem for the linear solver (deallocated by a call
 to lfree).

*******************************************************************/

void CVodeFree(void *cvode_mem)
{
  CVodeMemData* cv_mem;

  cv_mem = (CVodeMemData*) cvode_mem;
  
  if (cvode_mem == NULL) return;

  CVFreeVectors(cv_mem, cv_mem_qmax);
  if ((cv_mem_iter == NEWTON) && cv_mem_linitOK)
	  cv_mem_lfree(cv_mem);
  free(cv_mem);
}


/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/


/*******************************************************************/
/******** BEGIN Private Helper Functions Implementation ************/
/*******************************************************************/
 
/****************** CVAllocVectors ***********************************

 This routine allocates the CVODE vectors ewt, acor, tempv, ftemp, and
 zn[0], ..., zn[maxord]. The length of the vectors is the input
 parameter neq and the maximum order (needed to allocate zn) is the
 input parameter maxord. If all memory allocations are successful,
 CVAllocVectors returns TRUE. Otherwise all allocated memory is freed
 and CVAllocVectors returns FALSE.
 This routine also sets the optional outputs lrw and liw, which are
 (respectively) the lengths of the real and integer work spaces
 allocated here.

**********************************************************************/

static bool CVAllocVectors(CVodeMemData* cv_mem, integer neq, int maxord,
			   void *machEnv)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */
  
  cv_mem_ewt = N_VNew(neq, cv_mem_machenv);
  if (cv_mem_ewt == NULL)
		return(FALSE);
  
  cv_mem_acor = N_VNew(neq, machEnv);
  if (cv_mem_acor == NULL)
  {
    N_VFree(cv_mem_ewt);
    return(FALSE);
  }
  
  cv_mem_tempv = N_VNew(neq, machEnv);
  if (cv_mem_tempv == NULL)
  {
    N_VFree(cv_mem_ewt);
    N_VFree(cv_mem_acor);
    return(FALSE);
  }
  
  cv_mem_ftemp = N_VNew(neq, machEnv);
  if (cv_mem_ftemp == NULL)
  {
    N_VFree(cv_mem_tempv);
    N_VFree(cv_mem_ewt);
    N_VFree(cv_mem_acor);
    return(FALSE);
  }

  /* Allocate zn[0] ... zn[maxord] */

  for (j=0; j <= maxord; j++)
  {
    cv_mem_zn[j] = N_VNew(neq, machEnv);
    
	if (cv_mem_zn[j] == NULL)
	{
      N_VFree(cv_mem_ewt);
      N_VFree(cv_mem_acor);
      N_VFree(cv_mem_tempv);
      N_VFree(cv_mem_ftemp);
      for (i=0; i < j; i++)
		  N_VFree(cv_mem_zn[i]);
      return(FALSE);
    }
  }

  /* Set solver workspace lengths  */

  cv_mem_lrw = (maxord + 5)*neq;
  cv_mem_liw = 0;

  return(TRUE);
}

/***************** CVFreeVectors *********************************
  
 This routine frees the CVODE vectors allocated in CVAllocVectors.

******************************************************************/

static void CVFreeVectors(CVodeMemData* cv_mem, int maxord)
{
  int j;
  
  N_VFree(cv_mem_ewt);
  N_VFree(cv_mem_acor);
  N_VFree(cv_mem_tempv);
  N_VFree(cv_mem_ftemp);
  for(j=0; j <= maxord; j++)
	  N_VFree(cv_mem_zn[j]);
}

/*********************** CVEwtSet **************************************
  
 This routine is responsible for setting the error weight vector
 ewtvec, according to tol_type, as follows:

 (1) ewtvec[i] = 1 / (*rtol * ABS(ycur[i]) + *atol), i=0,...,neq-1
     if tol_type = SS
 (2) ewtvec[i] = 1 / (*rtol * ABS(ycur[i]) + atol[i]), i=0,...,neq-1
     if tol_type = SV

  CVEwtSet returns TRUE if ewtvec is successfully set as above to a
  positive vector and FALSE otherwise. In the latter case, ewtvec is
  considered undefined after the FALSE return from CVEwtSet.

  All the real work is done in the routines CVEwtSetSS, CVEwtSetSV.
 
***********************************************************************/

static bool CVEwtSet(CVodeMemData* cv_mem, real *rtol, void *atol, int tol_type, 
		     iN_Vector* ycur, iN_Vector* ewtvec, integer neq)
{
  switch(tol_type) {
  case SS: return(CVEwtSetSS(cv_mem, rtol, (real *)atol, ycur, ewtvec, neq));
  case SV: return(CVEwtSetSV(cv_mem, rtol, (iN_Vector*)atol, ycur, ewtvec, neq));
  }
  return FALSE;
}

/*********************** CVEwtSetSS *********************************

 This routine sets ewtvec as decribed above in the case tol_type=SS.
 It tests for non-positive components before inverting. CVEwtSetSS
 returns TRUE if ewtvec is successfully set to a positive vector
 and FALSE otherwise. In the latter case, ewtvec is considered
 undefined after the FALSE return from CVEwtSetSS.

********************************************************************/

static bool CVEwtSetSS(CVodeMemData* cv_mem, real *rtol, real *atol,
		       iN_Vector* ycur, iN_Vector* ewtvec, integer neq)
{
  real rtoli, atoli;
  
  rtoli = *rtol;
  atoli = *atol;
  N_VAbs(ycur, cv_mem_tempv);
  N_VScale(rtoli, cv_mem_tempv, cv_mem_tempv);
  N_VAddConst(cv_mem_tempv, atoli, cv_mem_tempv);
  if (N_VMin(cv_mem_tempv) <= ZERO)
	  return(FALSE);
  N_VInv(cv_mem_tempv, ewtvec);
  return(TRUE);
}

/*********************** CVEwtSetSV *********************************

 This routine sets ewtvec as decribed above in the case tol_type=SV.
 It tests for non-positive components before inverting. CVEwtSetSV
 returns TRUE if ewtvec is successfully set to a positive vector
 and FALSE otherwise. In the latter case, ewtvec is considered
 undefined after the FALSE return from CVEwtSetSV.

********************************************************************/

static bool CVEwtSetSV(CVodeMemData* cv_mem, real *rtol, iN_Vector* atol,
		       iN_Vector* ycur, iN_Vector* ewtvec, integer neq)
{
  real rtoli;
  
  rtoli = *rtol;
  N_VAbs(ycur, cv_mem_tempv);
  N_VLinearSum(rtoli, cv_mem_tempv, ONE, atol, cv_mem_tempv);
  if (N_VMin(cv_mem_tempv) <= ZERO)
	  return(FALSE);
  N_VInv(cv_mem_tempv, ewtvec);
  return(TRUE);
}

/******************* CVHin ***************************************

 This routine computes a tentative initial step size h0. 
 If tout is too close to tn (= t0), then CVHin returns FALSE and
 h remains uninitialized. Otherwise, CVHin sets h to the chosen 
 value h0 and returns TRUE.

 The algorithm used seeks to find h0 as a solution of
       (WRMS norm of (h0^2 ydd / 2)) = 1, 
 where ydd = estimated second derivative of y.

*****************************************************************/

static bool CVHin(CVodeMemData* cv_mem, real tout)
{
  int sign, count;
  real tdiff, tdist, tround, hlb, hub;
  real hg, hgs, hnew, hrat, h0, yddnrm;

  /* Test for tout too close to tn */
  
  if ((tdiff = tout-cv_mem_tn) == ZERO)
	  return(FALSE);
  
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = ABS(tdiff);
  tround = cv_mem_uround * MAX(ABS(cv_mem_tn), ABS(tout));
  if (tdist < TWO*tround)
	  return(FALSE);
  
  /* Set lower and upper bounds on h0, and take geometric mean 
     Exit with this value if the bounds cross each other       */

  hlb = HLB_FACTOR * tround;
  hub = CVUpperBoundH0(cv_mem, tdist);
  hg  = RSqrt(hlb*hub);
  if (hub < hlb) {
    if (sign == -1) hg = -hg;
    cv_mem_h = hg;
    return(TRUE);
  }
  
  /* Loop up to MAX_ITERS times to find h0.
     Stop if new and previous values differ by a factor < 2.
     Stop if hnew/hg > 2 after one iteration, as this probably means
     that the ydd value is bad because of cancellation error.        */

  count = 0;
  loop {
    hgs = hg*sign;
    yddnrm = CVYddNorm(cv_mem, hgs);
    hnew =  (yddnrm*hub*hub > TWO) ? RSqrt(TWO/yddnrm) : RSqrt(hg*hub);
    count++;
    if (count >= MAX_ITERS) break;
    hrat = hnew/hg;
    if ((hrat > HALF) && (hrat < TWO)) break;
    if ((count >= 2) && (hrat > TWO)) {
      hnew = hg;
      break;
    }
    hg = hnew;
  }
  
  /* Apply bounds, bias factor, and attach sign */

  h0 = H_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  cv_mem_h = h0;
  return(TRUE);
}

/******************** CVUpperBoundH0 ******************************

 This routine sets an upper bound on abs(h0) based on
 tdist = tn - t0 and the values of y[i]/ydot[i].

******************************************************************/

static real CVUpperBoundH0(CVodeMemData* cv_mem, real tdist)
{
  real atoli, hub_inv, hub;
  bool vectorAtol;
  iN_Vector* temp1, *temp2;

  vectorAtol = (cv_mem_itol == SV);
  if (!vectorAtol) atoli = *((real *) cv_mem_abstol);
  temp1 = cv_mem_tempv;
  temp2 = cv_mem_acor;
  N_VAbs(cv_mem_zn[0], temp1);
  N_VAbs(cv_mem_zn[1], temp2);
  if (vectorAtol) {
    N_VLinearSum(HUB_FACTOR, temp1, ONE, (iN_Vector*)cv_mem_abstol, temp1);
  } else {
    N_VScale(HUB_FACTOR, temp1, temp1);
    N_VAddConst(temp1, atoli, temp1);
  }
  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);
  hub = HUB_FACTOR*tdist;
  if (hub*hub_inv > ONE) hub = ONE/hub_inv;
  return(hub);
}

/****************** CVYddNorm *************************************

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.

******************************************************************/

static real CVYddNorm(CVodeMemData* cv_mem, real hg)
{
  real yddnrm;
  
  N_VLinearSum(hg, cv_mem_zn[1], ONE, cv_mem_zn[0], cv_mem_y);
  cv_mem_f(cv_mem_N, cv_mem_tn+hg, cv_mem_y, cv_mem_tempv, cv_mem_f_data);
  cv_mem_nfe++;
  N_VLinearSum(ONE, cv_mem_tempv, -ONE, cv_mem_zn[1], cv_mem_tempv);
  N_VScale(ONE/hg, cv_mem_tempv, cv_mem_tempv);

  yddnrm = N_VWrmsNorm(cv_mem_tempv, cv_mem_ewt);
  return(yddnrm);
}

/********************* CVStep **************************************
 
 This routine performs one internal cvode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
  * preliminary adjustments if a new step size was chosen;
  * prediction of the Nordsieck history array zn at tn + h;
  * setting of multistep method coefficients and test quantities;
  * solution of the nonlinear system;
  * testing the local error;
  * updating zn and other state data if successful;
  * resetting stepsize and order for the next step.

 On a failure in the nonlinear system solution or error test, the
 step may be reattempted, depending on the nature of the failure.

********************************************************************/

int CVStep(CVodeMemData* cv_mem)
{
  real saved_t, dsm;
  int ncf, nef, nflag, kflag;
  bool passed;
  
  saved_t = cv_mem_tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;
  
  if ((cv_mem_nst > 0) && (cv_mem_hprime != cv_mem_h)) CVAdjustParams(cv_mem);
  
  /* Looping point for attempts to take a step */
  loop {  
    CVPredict(cv_mem);  
    CVSet(cv_mem);

    nflag = CVnls(cv_mem, nflag);
    kflag = CVHandleNFlag(cv_mem, &nflag, saved_t, &ncf);
    if (kflag == PREDICT_AGAIN) continue;
    if (kflag != DO_ERROR_TEST) return(kflag);
    /* Return if nonlinear solve failed and recovery not possible. */

    passed = CVDoErrorTest(cv_mem, &nflag, &kflag, saved_t, &nef, &dsm);
    if ((!passed) && (kflag == REP_ERR_FAIL)) return(kflag);
    /* Return if error test failed and recovery not possible. */
    if (passed) break;
    /* Retry step if error test failed, nflag == PREV_ERR_FAIL */
  }

  /* Nonlinear system solve and error test were both successful;
     update data, and consider change of step and/or order       */

  CVCompleteStep(cv_mem);
  CVPrepareNextStep(cv_mem, dsm);

  return(SUCCESS_STEP);
}

/********************* CVAdjustParams ********************************

 This routine is called when a change in step size was decided upon,
 and it handles the required adjustments to the history array zn.
 If there is to be a change in order, we call CVAdjustOrder and reset
 q, L = q+1, and qwait.  Then in any case, we call CVRescale, which
 resets h and rescales the Nordsieck array.

**********************************************************************/

static void CVAdjustParams(CVodeMemData* cv_mem)
{
  if (cv_mem_qprime != cv_mem_q) {
    CVAdjustOrder(cv_mem, cv_mem_qprime-cv_mem_q);
    cv_mem_q = cv_mem_qprime;
    cv_mem_L = cv_mem_q+1;
    cv_mem_qwait = cv_mem_L;
  }
  CVRescale(cv_mem);
}

/********************* CVAdjustOrder *****************************

  This routine is a high level routine which handles an order
  change by an amount deltaq (= +1 or -1). If a decrease in order
  is requested and q==2, then the routine returns immediately.
  Otherwise CVAdjustAdams or CVAdjustBDF is called to handle the
  order change (depending on the value of lmm).

******************************************************************/

static void CVAdjustOrder(CVodeMemData* cv_mem, int deltaq)
{
  if ((cv_mem_q==2) && (deltaq != 1)) return;
  
  switch(cv_mem_lmm){
    case ADAMS: CVAdjustAdams(cv_mem, deltaq);
                break;
    case BDF:   CVAdjustBDF(cv_mem, deltaq);
                break;
  }
}

/*************** CVAdjustAdams ***********************************

 This routine adjusts the history array on a change of order q by
 deltaq, in the case that lmm == ADAMS.

*****************************************************************/

static void CVAdjustAdams(CVodeMemData* cv_mem, int deltaq)
{
  int i, j;
  real xi, hsum;

  /* On an order increase, set new column of zn to zero and return */
  
  if (deltaq==1) {
    N_VConst(ZERO, cv_mem_zn[cv_mem_L]);
    return;
  }

  /* On an order decrease, each zn[j] is adjusted by a multiple
     of zn[q].  The coefficients in the adjustment are the 
     coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
     integrated, where xi_j = [t_n - t_(n-j)]/h.               */

  for (i=0; i <= cv_mem_qmax; i++) cv_mem_l[i] = ZERO;
  cv_mem_l[1] = ONE;
  hsum = ZERO;
  for (j=1; j <= cv_mem_q-2; j++) {
    hsum += cv_mem_tau[j];
    xi = hsum / cv_mem_hscale;
    for (i=j+1; i >= 1; i--)
		cv_mem_l[i] = cv_mem_l[i]*xi + cv_mem_l[i-1];
  }
  
  for (j=1; j <= cv_mem_q-2; j++)
	  cv_mem_l[j+1] = cv_mem_q * (cv_mem_l[j] / (j+1));
  
  for (j=2; j < cv_mem_q; j++)
    N_VLinearSum(-cv_mem_l[j], cv_mem_zn[cv_mem_q], ONE, cv_mem_zn[j], cv_mem_zn[j]);
}

/********************** CVAdjustBDF *******************************

 This is a high level routine which handles adjustments to the
 history array on a change of order by deltaq in the case that 
 lmm == BDF.  CVAdjustBDF calls CVIncreaseBDF if deltaq = +1 and 
 CVDecreaseBDF if deltaq = -1 to do the actual work.

******************************************************************/

static void CVAdjustBDF(CVodeMemData* cv_mem, int deltaq)
{
  switch(deltaq) {
    case 1 : CVIncreaseBDF(cv_mem);
             return;
    case -1: CVDecreaseBDF(cv_mem);
             return;
  }
}

/******************** CVIncreaseBDF **********************************

 This routine adjusts the history array on an increase in the 
 order q in the case that lmm == BDF.  
 A new column zn[q+1] is set equal to a multiple of the saved 
 vector (= acor) in zn[qmax].  Then each zn[j] is adjusted by
 a multiple of zn[q+1].  The coefficients in the adjustment are the 
 coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
 where xi_j = [t_n - t_(n-j)]/h.

*********************************************************************/

static void CVIncreaseBDF(CVodeMemData* cv_mem)
{
  real alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;
  
  for (i=0; i <= cv_mem_qmax; i++)
	  cv_mem_l[i] = ZERO;
  cv_mem_l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = cv_mem_hscale;
  if (cv_mem_q > 1) {
    for (j=1; j < cv_mem_q; j++) {
      hsum += cv_mem_tau[j+1];
      xi = hsum / cv_mem_hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--)
		  cv_mem_l[i] = cv_mem_l[i]*xiold + cv_mem_l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;
  N_VScale(A1, cv_mem_zn[cv_mem_qmax], cv_mem_zn[cv_mem_L]);
  for (j=2; j <= cv_mem_q; j++) {
    N_VLinearSum(cv_mem_l[j], cv_mem_zn[cv_mem_L], ONE, cv_mem_zn[j], cv_mem_zn[j]);
  }  
}

/********************* CVDecreaseBDF ******************************

 This routine adjusts the history array on a decrease in the 
 order q in the case that lmm == BDF.  
 Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 in the adjustment are the coefficients of the polynomial
 x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h.

******************************************************************/

static void CVDecreaseBDF(CVodeMemData* cv_mem)
{
  real hsum, xi;
  int i, j;
  
  for (i=0; i <= cv_mem_qmax; i++)
	  cv_mem_l[i] = ZERO;
  cv_mem_l[2] = ONE;
  hsum = ZERO;
  for(j=1; j <= cv_mem_q-2; j++) {
    hsum += cv_mem_tau[j];
    xi = hsum /cv_mem_hscale;
    for (i=j+2; i >= 2; i--)
		cv_mem_l[i] = cv_mem_l[i]*xi + cv_mem_l[i-1];
  }
  
  for(j=2; j < cv_mem_q; j++)
    N_VLinearSum(-cv_mem_l[j], cv_mem_zn[cv_mem_q], ONE, cv_mem_zn[j], cv_mem_zn[j]);
}

/**************** CVRescale ***********************************

  This routine rescales the Nordsieck array by multiplying the
  jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
  h is rescaled by eta, and hscale is reset to h.

***************************************************************/

static void CVRescale(CVodeMemData* cv_mem)
{
  int j;
  real factor;
  
  factor = cv_mem_eta;
  for (j=1; j <= cv_mem_q; j++) {
    N_VScale(factor, cv_mem_zn[j], cv_mem_zn[j]);
    factor *= cv_mem_eta;
  }
  cv_mem_h = cv_mem_hscale * cv_mem_eta;
  cv_mem_hscale = cv_mem_h;
}

/********************* CVPredict *************************************

 This routine advances tn by the tentative step size h, and computes
 the predicted array z_n(0), which is overwritten on zn.  The
 prediction of zn is done by repeated additions.

*********************************************************************/

static void CVPredict(CVodeMemData* cv_mem)
{
  int j, k;
  
  cv_mem_tn += cv_mem_h;
  for (k = 1; k <= cv_mem_q; k++)
    for (j = cv_mem_q; j >= k; j--) 
      N_VLinearSum(ONE, cv_mem_zn[j-1], ONE, cv_mem_zn[j], cv_mem_zn[j-1]); 
}

/************************** CVSet *********************************

 This routine is a high level routine which calls CVSetAdams or
 CVSetBDF to set the polynomial l, the test quantity array tq, 
 and the related variables  rl1, gamma, and gamrat.

******************************************************************/

static void CVSet(CVodeMemData* cv_mem)
{
  switch(cv_mem_lmm) {
    case ADAMS: CVSetAdams(cv_mem);
                break;
    case BDF  : CVSetBDF(cv_mem);
                break;
  }
  cv_mem_rl1 = ONE / cv_mem_l[1];
  cv_mem_gamma = cv_mem_h * cv_mem_rl1;
  if (cv_mem_nst == 0)
	  cv_mem_gammap = cv_mem_gamma;
  cv_mem_gamrat = (cv_mem_nst > 0) ? cv_mem_gamma / cv_mem_gammap : ONE;  /* protect x / x != 1.0 */
}

/******************** CVSetAdams *********************************

 This routine handles the computation of l and tq for the
 case lmm == ADAMS.

 The components of the vector l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                          q-1
 (d/dx) Lambda(x) = c * PRODUCT (1 + x / xi_i) , where
                          i=1
 Lambda(-1) = 0, Lambda(0) = 1, and c is a normalization factor.
 Here xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.

*****************************************************************/

static void CVSetAdams(CVodeMemData* cv_mem)
{
  real m[L_MAX], M[3], hsum;
  
  if (cv_mem_q == 1) {
    cv_mem_l[0] = cv_mem_l[1] = cv_mem_tq[1] = cv_mem_tq[5] = ONE;
    cv_mem_tq[2] = TWO;
    cv_mem_tq[3] = TWELVE;
    cv_mem_tq[4] = CORTES * cv_mem_tq[2];       /* = 0.1 * tq[2] */
    return;
  }
  
  hsum = CVAdamsStart(cv_mem, m);
  
  M[0] = CVAltSum(cv_mem_q-1, m, 1);
  M[1] = CVAltSum(cv_mem_q-1, m, 2);
  
  CVAdamsFinish(cv_mem, m, M, hsum);
}

/****************** CVAdamsStart ********************************

 This routine generates in m[] the coefficients of the product
 polynomial needed for the Adams l and tq coefficients for q > 1.
  
******************************************************************/

static real CVAdamsStart(CVodeMemData* cv_mem, real m[])
{
  real hsum, xi_inv, sum;
  int i, j;
  
  hsum = cv_mem_h;
  m[0] = ONE;
  for (i=1; i <= cv_mem_q; i++) m[i] = ZERO;
  for (j=1; j < cv_mem_q; j++) {
    if ((j==cv_mem_q-1) && (cv_mem_qwait == 1)) {
      sum = CVAltSum(cv_mem_q-2, m, 2);
      cv_mem_tq[1] = m[cv_mem_q-2] / (cv_mem_q * sum);
    }
    xi_inv = cv_mem_h / hsum;
    for (i=j; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    hsum += cv_mem_tau[j];
    /* The m[i] are coefficients of product(1 to j) (1 + x/xi_i) */
  }
  return(hsum);
}

/****************** CVAdamsFinish  *******************************

 This routine completes the calculation of the Adams l and tq.

******************************************************************/

static void CVAdamsFinish(CVodeMemData* cv_mem, real m[], real M[], real hsum)
{
  int i;
  real M0_inv, xi, xi_inv;
  
  M0_inv = ONE / M[0];
  
  cv_mem_l[0] = ONE;
  for (i=1; i <= cv_mem_q; i++)
	  cv_mem_l[i] = M0_inv * (m[i-1] / i);
  xi = hsum / cv_mem_h;
  xi_inv = ONE / xi;
  
  cv_mem_tq[2] = xi * M[0] / M[1];
  cv_mem_tq[5] = xi / cv_mem_l[cv_mem_q];

  if (cv_mem_qwait == 1) {
    for (i=cv_mem_q; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    M[2] = CVAltSum(cv_mem_q, m, 2);
    cv_mem_tq[3] = cv_mem_L * M[0] / M[2];
  }

  cv_mem_tq[4] = CORTES * cv_mem_tq[2];
}

/****************** CVAltSum **************************************
  
 CVAltSum returns the value of the alternating sum
   sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 If iend < 0 then CVAltSum returns 0.
 This operation is needed to compute the integral, from -1 to 0,
 of a polynomial x^(k-1) M(x) given the coefficients of M(x).

******************************************************************/

static real CVAltSum(int iend, real a[], int k)
{
  int i, sign;
  real sum;
  
  if (iend < 0) return(ZERO);
  
  sum = ZERO;
  sign = 1;
  for (i=0; i <= iend; i++) {
    sum += sign * (a[i] / (i+k));
    sign = -sign;
  }
  return(sum);
}

/***************** CVSetBDF **************************************

 This routine computes the coefficients l and tq in the case
 lmm == BDF.  CVSetBDF calls CVSetTqBDF to set the test
 quantity vector tq. 

 The components of the vector l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                                 q-1
 Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
                                 i=1
 xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.


*****************************************************************/

static void CVSetBDF(CVodeMemData* cv_mem)
{
  real alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;
  
  cv_mem_l[0] = cv_mem_l[1] = xi_inv = xistar_inv = ONE;
  for (i=2; i <= cv_mem_q; i++)
	  cv_mem_l[i] = ZERO;
  alpha0 = alpha0_hat = -ONE;
  hsum = cv_mem_h;
  if (cv_mem_q > 1) {
    for (j=2; j < cv_mem_q; j++) {
      hsum += cv_mem_tau[j-1];
      xi_inv = cv_mem_h / hsum;
      alpha0 -= ONE / j;
      for(i=j; i >= 1; i--)
		  cv_mem_l[i] += cv_mem_l[i-1]*xi_inv;
      /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }
    
    /* j = q */
    alpha0 -= ONE / cv_mem_q;
    xistar_inv = -cv_mem_l[1] - alpha0;
    hsum += cv_mem_tau[cv_mem_q-1];
    xi_inv = cv_mem_h / hsum;
    alpha0_hat = -cv_mem_l[1] - xi_inv;
    for (i=cv_mem_q; i >= 1; i--)
		cv_mem_l[i] += cv_mem_l[i-1]*xistar_inv;
  }

  CVSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/****************** CVSetTqBDF ************************************

 This routine sets the test quantity vector tq in the case
 lmm == BDF.

******************************************************************/

static void CVSetTqBDF(CVodeMemData* cv_mem, real hsum, real alpha0,
		       real alpha0_hat, real xi_inv, real xistar_inv)
{
  real A1, A2, A3, A4, A5, A6;
  real C, CPrime, CPrimePrime;
  
  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE + cv_mem_q * A1;
  cv_mem_tq[2] = ABS(alpha0 * (A2 / A1));
  cv_mem_tq[5] = ABS((A2) / (cv_mem_l[cv_mem_q] * xi_inv/xistar_inv));
  if (cv_mem_qwait == 1) {
    C = xistar_inv / cv_mem_l[cv_mem_q];
    A3 = alpha0 + ONE / cv_mem_q;
    A4 = alpha0_hat + xi_inv;
    CPrime = A3 / (ONE - A4 + A3);
    cv_mem_tq[1] = ABS(CPrime / C);
    hsum += cv_mem_tau[cv_mem_q];
    xi_inv = cv_mem_h / hsum;
    A5 = alpha0 - (ONE / (cv_mem_q+1));
    A6 = alpha0_hat - xi_inv;
    CPrimePrime = A2 / (ONE - A6 + A5);
    cv_mem_tq[3] = ABS(CPrimePrime * xi_inv * (cv_mem_q+2) * A5);
  }
  cv_mem_tq[4] = CORTES * cv_mem_tq[2];
}

/****************** CVnls *****************************************

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 Depending on iter, it calls CVnlsFunctional or CVnlsNewton
 to do the work.

******************************************************************/

static int CVnls(CVodeMemData* cv_mem, int nflag)
{
  switch(cv_mem_iter) {
    case FUNCTIONAL : return(CVnlsFunctional(cv_mem));
    case NEWTON     : return(CVnlsNewton(cv_mem, nflag));
    default :
        return(CVnlsNewton(cv_mem, nflag));
  }
}

/***************** CVnlsFunctional ********************************

 This routine attempts to solve the nonlinear system using 
 functional iteration (no matrices involved).

******************************************************************/

static int CVnlsFunctional(CVodeMemData* cv_mem)
{
  int m;
  real del, delp, dcon;

  /* Initialize counter and evaluate f at predicted y */
  
  cv_mem_crate = ONE;
  m = 0;
  cv_mem_f(cv_mem_N, cv_mem_tn, cv_mem_zn[0], cv_mem_tempv, cv_mem_f_data);
  cv_mem_nfe++;
  N_VConst(ZERO, cv_mem_acor);

  /* Loop until convergence; accumulate corrections in acor */

  loop {
    /* Correct y directly from the last f value */
    N_VLinearSum(cv_mem_h, cv_mem_tempv, -ONE, cv_mem_zn[1], cv_mem_tempv);
    N_VScale(cv_mem_rl1, cv_mem_tempv, cv_mem_tempv);
    N_VLinearSum(ONE, cv_mem_zn[0], ONE, cv_mem_tempv, cv_mem_y);
    /* Get WRMS norm of current correction to use in convergence test */
    N_VLinearSum(ONE, cv_mem_tempv, -ONE, cv_mem_acor, cv_mem_acor);
    del = N_VWrmsNorm(cv_mem_acor, cv_mem_ewt);
    N_VScale(ONE, cv_mem_tempv, cv_mem_acor);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0)
		cv_mem_crate = MAX(CRDOWN * cv_mem_crate, del / delp);
    dcon = del * MIN(ONE, cv_mem_crate) / cv_mem_tq[4];
    if (dcon <= ONE) {
      cv_mem_acnrm = (m == 0) ? del : N_VWrmsNorm(cv_mem_acor, cv_mem_ewt);
      return(SOLVED);  /* Convergence achieved */
    }

    /* Stop at maxcor iterations or if iter. seems to be diverging */
    m++;
    if ((m==cv_mem_maxcor) || ((m >= 2) && (del > RDIV * delp)))
      return(CONV_FAIL);
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    cv_mem_f(cv_mem_N, cv_mem_tn, cv_mem_y, cv_mem_tempv, cv_mem_f_data);
    cv_mem_nfe++;
  }
}

/*********************** CVnlsNewton **********************************

 This routine handles the Newton iteration. It calls lsetup if 
 indicated, calls CVNewtonIteration to perform the iteration, and 
 retries a failed attempt at Newton iteration if that is indicated.
 See return values at top of this file.

**********************************************************************/

static int CVnlsNewton(CVodeMemData* cv_mem, int nflag)
{
  iN_Vector* vtemp1, *vtemp2, *vtemp3;
  int convfail, ier;
  bool callSetup;
  
  vtemp1 = cv_mem_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = cv_mem_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = cv_mem_tempv; /* rename tempv as vtemp3 for readability */
  
  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    NO_FAILURES : FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (cv_mem_setupNonNull) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (cv_mem_nst == 0) || (cv_mem_nst >= cv_mem_nstlp + MSBP) || (ABS(cv_mem_gamrat-ONE) > DGMAX);
  } else {  
    cv_mem_crate = ONE;
    callSetup = FALSE;
  }
  
  /* Looping point for the solution of the nonlinear system.
     Evaluate f at the predicted y, call lsetup if indicated, and
     call CVNewtonIteration for the Newton iteration itself.      */
  
  loop {

    cv_mem_f(cv_mem_N, cv_mem_tn, cv_mem_zn[0], cv_mem_ftemp, cv_mem_f_data);
    cv_mem_nfe++; 
    
    if (callSetup) {
      ier = cv_mem_lsetup(cv_mem, convfail, cv_mem_zn[0], cv_mem_ftemp, &cv_mem_jcur, 
		   vtemp1, vtemp2, vtemp3);
      cv_mem_nsetups++;
      callSetup = FALSE;
      cv_mem_gamrat = cv_mem_crate = ONE; 
      cv_mem_gammap = cv_mem_gamma;
      cv_mem_nstlp = cv_mem_nst;
      /* Return if lsetup failed */
      if (ier < 0)
		  return(SETUP_FAIL_UNREC);
      if (ier > 0)
		  return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, cv_mem_acor);
    N_VScale(ONE, cv_mem_zn[0], cv_mem_y);

    /* Do the Newton iteration */
    ier = CVNewtonIteration(cv_mem);

    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=FAIL_BAD_J.  Otherwise return.                 */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = TRUE;
    convfail = FAIL_BAD_J;
  }
}

/********************** CVNewtonIteration ****************************

 This routine performs the Newton iteration. If the iteration succeeds,
 it returns the value SOLVED. If not, it may signal the CVnlsNewton 
 routine to call lsetup again and reattempt the iteration, by
 returning the value TRY_AGAIN. (In this case, CVnlsNewton must set 
 convfail to FAIL_BAD_J before calling setup again). 
 Otherwise, this routine returns one of the appropriate values 
 SOLVE_FAIL_UNREC or CONV_FAIL back to CVnlsNewton.

*********************************************************************/

static int CVNewtonIteration(CVodeMemData* cv_mem)
{
  int m, ret;
  real del, delp, dcon;
  iN_Vector* b;
  
  
  cv_mem_mnewt = m = 0;

  /* Looping point for Newton iteration */
  loop {

    /* Evaluate the residual of the nonlinear system*/
    N_VLinearSum(cv_mem_rl1, cv_mem_zn[1], ONE, cv_mem_acor, cv_mem_tempv);
    N_VLinearSum(cv_mem_gamma, cv_mem_ftemp, -ONE, cv_mem_tempv, cv_mem_tempv);

    /* Call the lsolve function */
    b = cv_mem_tempv;
    ret = cv_mem_lsolve(cv_mem, b, cv_mem_y, cv_mem_ftemp); 
    cv_mem_nni++;
    
    if (ret < 0)
		return(SOLVE_FAIL_UNREC);
    
    /* If lsolve had a recoverable failure and Jacobian data is
       not current, signal to try the solution again            */
    if (ret > 0) { 
      if ((!cv_mem_jcur) && (cv_mem_setupNonNull))
		  return(TRY_AGAIN);
      return(CONV_FAIL);
    }
/* *************** */
    /* Get WRMS norm of correction; add correction to acor and y */
    del = N_VWrmsNorm(b, cv_mem_ewt);
    N_VLinearSum(ONE, cv_mem_acor, ONE, b, cv_mem_acor);
    N_VLinearSum(ONE, cv_mem_zn[0], ONE, cv_mem_acor, cv_mem_y);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) {
      cv_mem_crate = MAX(CRDOWN * cv_mem_crate, del/delp);
    }
    dcon = del * MIN(ONE, cv_mem_crate) / cv_mem_tq[4];
    
    if (dcon <= ONE) {
      cv_mem_acnrm = (m==0) ? del : N_VWrmsNorm(cv_mem_acor, cv_mem_ewt);
      cv_mem_jcur = FALSE;
      return(SOLVED); /* Nonlinear system was solved successfully */
    }
    
    cv_mem_mnewt = ++m;
    
    /* Stop at maxcor iterations or if iter. seems to be diverging.
       If still not converged and Jacobian data is not current, 
       signal to try the solution again                            */
    if ((m == cv_mem_maxcor) || ((m >= 2) && (del > RDIV*delp))) {
      if ((!cv_mem_jcur) && (cv_mem_setupNonNull))
		  return(TRY_AGAIN);
      return(CONV_FAIL);
    }
    
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    cv_mem_f(cv_mem_N, cv_mem_tn, cv_mem_y, cv_mem_ftemp, cv_mem_f_data);
    cv_mem_nfe++;
  }
}

/********************** CVHandleNFlag *******************************

 This routine takes action on the return value nflag = *nflagPtr
 returned by CVnls, as follows:
 
 If CVnls succeeded in solving the nonlinear system, then
 CVHandleNFlag returns the constant DO_ERROR_TEST, which tells CVStep
 to perform the error test.

 If the nonlinear system was not solved successfully, then ncfn and
 ncf = *ncfPtr are incremented and Nordsieck array zn is restored.

 If the solution of the nonlinear system failed due to an
 unrecoverable failure by setup, we return the value SETUP_FAILED.

 If it failed due to an unrecoverable failure in solve, then we return
 the value SOLVE_FAILED.

 Otherwise, a recoverable failure occurred when solving the 
 nonlinear system (CVnls returned nflag == CONV_FAIL). 
   In this case, we return the value REP_CONV_FAIL if ncf is now
   equal to MXNCF or |h| = hmin. 
   If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
   PREDICT_AGAIN, telling CVStep to reattempt the step.

*********************************************************************/

static int CVHandleNFlag(CVodeMemData* cv_mem, int *nflagPtr, real saved_t,
			 int *ncfPtr)
{
  int nflag;
  
  nflag = *nflagPtr;
  
  if (nflag == SOLVED) return(DO_ERROR_TEST);

  /* The nonlinear soln. failed; increment ncfn and restore zn */
  cv_mem_ncfn++;
  CVRestore(cv_mem, saved_t);
  
  /* Return if lsetup or lsolve failed unrecoverably */
  if (nflag == SETUP_FAIL_UNREC) return(SETUP_FAILED);
  if (nflag == SOLVE_FAIL_UNREC) return(SOLVE_FAILED);
  
  /* At this point, nflag == CONV_FAIL; increment ncf */
  
  (*ncfPtr)++;
  cv_mem_etamax = ONE;
  /* If we had MXNCF failures or |h| = hmin, return REP_CONV_FAIL */
  if ((ABS(cv_mem_h) <= cv_mem_hmin*ONEPSM) || (*ncfPtr == MXNCF))
    return(REP_CONV_FAIL);

  /* Reduce step size; return to reattempt the step */
  cv_mem_eta = MAX(ETACF, cv_mem_hmin / ABS(cv_mem_h));
  *nflagPtr = PREV_CONV_FAIL;
  CVRescale(cv_mem);
  return(PREDICT_AGAIN);
}

/********************** CVRestore ************************************

 This routine restores the value of tn to saved_t and undoes the
 prediction.  After execution of CVRestore, the Nordsieck array zn has
 the same values as before the call to CVPredict.

********************************************************************/

static void CVRestore(CVodeMemData* cv_mem, real saved_t)
{
  int j, k;
  
  cv_mem_tn = saved_t;
  for (k = 1; k <= cv_mem_q; k++)
    for (j = cv_mem_q; j >= k; j--)
      N_VLinearSum(ONE, cv_mem_zn[j-1], -ONE, cv_mem_zn[j], cv_mem_zn[j-1]);
}

/******************* CVDoErrorTest ********************************

 This routine performs the local error test. 
 The weighted local error norm dsm is loaded into *dsmPtr, and 
 the test dsm ?<= 1 is made.

 If the test passes, CVDoErrorTest returns TRUE. 

 If the test fails, we undo the step just taken (call CVRestore), 
 set *nflagPtr to PREV_ERR_FAIL, and return FALSE. 

 If MXNEF error test failures have occurred or if ABS(h) = hmin,
 we set *kflagPtr = REP_ERR_FAIL. (Otherwise *kflagPtr has the
 value last returned by CVHandleNflag.)

 If more than MXNEF1 error test failures have occurred, an order
 reduction is forced.

******************************************************************/

static bool CVDoErrorTest(CVodeMemData* cv_mem, int *nflagPtr, int *kflagPtr,
			  real saved_t, int *nefPtr, real *dsmPtr)
{
  real dsm;
  
  dsm = cv_mem_acnrm / cv_mem_tq[2];

  /* If est. local error norm dsm passes test, return TRUE */  
  *dsmPtr = dsm; 
  if (dsm <= ONE) return(TRUE);
  
  /* Test failed; increment counters, set nflag, and restore zn array */
  (*nefPtr)++;
  cv_mem_netf++;
  *nflagPtr = PREV_ERR_FAIL;
  CVRestore(cv_mem, saved_t);

  /* At MXNEF failures or |h| = hmin, return with kflag = REP_ERR_FAIL */
  if ((ABS(cv_mem_h) <= cv_mem_hmin*ONEPSM) || (*nefPtr == MXNEF)) {
    *kflagPtr = REP_ERR_FAIL;
    return(FALSE);
  }

  /* Set etamax = 1 to prevent step size increase at end of this step */
  cv_mem_etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  if (*nefPtr <= MXNEF1) {
    cv_mem_eta = ONE / (RPowerR(BIAS2*dsm,ONE/cv_mem_L) + ADDON);
    cv_mem_eta = MAX(ETAMIN, MAX(cv_mem_eta, cv_mem_hmin / ABS(cv_mem_h)));
    if (*nefPtr >= SMALL_NEF)
		cv_mem_eta = MIN(cv_mem_eta, ETAMXF);
    CVRescale(cv_mem);
    return(FALSE);
  }
  
  /* After MXNEF1 failures, force an order reduction and retry step */
  if (cv_mem_q > 1) {
    cv_mem_eta = MAX(ETAMIN, cv_mem_hmin / ABS(cv_mem_h));
    CVAdjustOrder(cv_mem,-1);
    cv_mem_L = cv_mem_q;
    cv_mem_q--;
    cv_mem_qwait = cv_mem_L;
    CVRescale(cv_mem);
    return(FALSE);
  }

  /* If already at order 1, restart: reload zn from scratch */
  cv_mem_eta = MAX(ETAMIN, cv_mem_hmin / ABS(cv_mem_h));
  cv_mem_h *= cv_mem_eta;
  cv_mem_hscale = cv_mem_h;
  cv_mem_qwait = LONG_WAIT;
  cv_mem_f(cv_mem_N, cv_mem_tn, cv_mem_zn[0], cv_mem_tempv, cv_mem_f_data);
  cv_mem_nfe++;
  N_VScale(cv_mem_h, cv_mem_tempv, cv_mem_zn[1]);
  return(FALSE);
}

/*************** CVCompleteStep **********************************

 This routine performs various update operations when the solution
 to the nonlinear system has passed the local error test. 
 We increment the step counter nst, record the values hu and qu,
 update the tau array, and apply the corrections to the zn array.
 The tau[i] are the last q values of h, with tau[1] the most recent.
 The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 we save acor and tq[5] for a possible order increase.

******************************************************************/

static void CVCompleteStep(CVodeMemData* cv_mem)
{
  int i, j;
  
  cv_mem_nst++;
  cv_mem_hu = cv_mem_h;
  cv_mem_qu = cv_mem_q;

  for (i=cv_mem_q; i >= 2; i--) 
	  cv_mem_tau[i] = cv_mem_tau[i-1];
  if ((cv_mem_q==1) && (cv_mem_nst > 1))
	  cv_mem_tau[2] = cv_mem_tau[1];
  cv_mem_tau[1] = cv_mem_h;

  for (j=0; j <= cv_mem_q; j++) 
    N_VLinearSum(cv_mem_l[j], cv_mem_acor, ONE, cv_mem_zn[j], cv_mem_zn[j]);
  cv_mem_qwait--;
  if ((cv_mem_qwait == 1) && (cv_mem_q != cv_mem_qmax)) {
    N_VScale(ONE, cv_mem_acor, cv_mem_zn[cv_mem_qmax]);
    cv_mem_saved_tq5 = cv_mem_tq[5];
  }
}

/************* CVPrepareNextStep **********************************

 This routine handles the setting of stepsize and order for the
 next step -- hprime and qprime.  Along with hprime, it sets the
 ratio eta = hprime/h.  It also updates other state variables 
 related to a change of step size or order.  Finally, we rescale
 the acor array to be the estimated local error vector.

******************************************************************/

static void CVPrepareNextStep(CVodeMemData* cv_mem, real dsm)
{
  real etaqm1, etaq, etaqp1;
  
  /* If etamax = 1, defer step size or order changes */
  if (cv_mem_etamax == ONE) {
    cv_mem_qwait = MAX(cv_mem_qwait, 2);
    cv_mem_qprime = cv_mem_q;
    cv_mem_hprime = cv_mem_h;
    cv_mem_eta = ONE;
    cv_mem_etamax = (cv_mem_nst <= SMALL_NST) ? ETAMX2 : ETAMX3;
    N_VScale(ONE/cv_mem_tq[2], cv_mem_acor, cv_mem_acor);
    return;
  }

  /* etaq is the ratio of new to old h at the current order */  
  etaq = ONE /(RPowerR(BIAS2*dsm,ONE/cv_mem_L) + ADDON);
  
  /* If no order change, adjust eta and acor in CVSetEta and return */
  if (cv_mem_qwait != 0) {
    cv_mem_eta = etaq;
    cv_mem_qprime = cv_mem_q;
    CVSetEta(cv_mem);
    return;
  }
  
  /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are 
     the ratios of new to old h at orders q-1 and q+1, respectively.
     CVChooseEta selects the largest; CVSetEta adjusts eta and acor */
  cv_mem_qwait = 2; 
  etaqm1 = CVComputeEtaqm1(cv_mem);
  etaqp1 = CVComputeEtaqp1(cv_mem);
  CVChooseEta(cv_mem, etaqm1, etaq, etaqp1);
  CVSetEta(cv_mem);
}

/***************** CVSetEta ***************************************

 This routine adjusts the value of eta according to the various
 heuristic limits and the optional input hmax.  It also resets
 etamax and rescales acor to be the estimated local error vector.

*******************************************************************/

static void CVSetEta(CVodeMemData* cv_mem)
{

  /* If eta below the threshhold THRESH, reject a change of step size */
  if (cv_mem_eta < THRESH) {
    cv_mem_eta = ONE;
    cv_mem_hprime = cv_mem_h;
  } else {
    /* Limit eta by etamax and hmax, then set hprime */
    cv_mem_eta = MIN(cv_mem_eta, cv_mem_etamax);
    cv_mem_eta /= MAX(ONE, ABS(cv_mem_h)*cv_mem_hmax_inv*cv_mem_eta);
    cv_mem_hprime = cv_mem_h * cv_mem_eta;
  }

  /* Reset etamx for the next step size change, and scale acor */
  cv_mem_etamax = (cv_mem_nst <= SMALL_NST) ? ETAMX2 : ETAMX3;
  N_VScale(ONE/cv_mem_tq[2], cv_mem_acor, cv_mem_acor);
}

/*************** CVComputeEtaqm1 **********************************

 This routine computes and returns the value of etaqm1 for a
 possible decrease in order by 1.

******************************************************************/

static real CVComputeEtaqm1(CVodeMemData* cv_mem)
{
  real etaqm1, ddn;
  
  etaqm1 = ZERO;
  if (cv_mem_q > 1) {
    ddn = N_VWrmsNorm(cv_mem_zn[cv_mem_q], cv_mem_ewt) / cv_mem_tq[1];
    etaqm1 = ONE/(RPowerR(BIAS1*ddn, ONE/cv_mem_q) + ADDON);
  }
  return(etaqm1);
}

/*************** CVComputeEtaqp1 **********************************

 This routine computes and returns the value of etaqp1 for a
 possible increase in order by 1.

******************************************************************/

static real CVComputeEtaqp1(CVodeMemData* cv_mem)
{
  real etaqp1, dup, cquot;
  
  etaqp1 = ZERO;
  if (cv_mem_q != cv_mem_qmax) {
    cquot = (cv_mem_tq[5] / cv_mem_saved_tq5) * RPowerI(cv_mem_h/cv_mem_tau[2], cv_mem_L);
    N_VLinearSum(-cquot, cv_mem_zn[cv_mem_qmax], ONE, cv_mem_acor, cv_mem_tempv);
    dup = N_VWrmsNorm(cv_mem_tempv, cv_mem_ewt) /cv_mem_tq[3];
    etaqp1 = ONE / (RPowerR(BIAS3*dup, ONE/(cv_mem_L+1)) + ADDON);
  }
  return(etaqp1);
}

/******************* CVChooseEta **********************************

 Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 q - 1, q, or q + 1, respectively), this routine chooses the 
 maximum eta value, sets eta to that value, and sets qprime to the
 corresponding value of q.  If there is a tie, the preference
 order is to (1) keep the same order, then (2) decrease the order,
 and finally (3) increase the order.  If the maximum eta value
 is below the threshhold THRESH, the order is kept unchanged and
 eta is set to 1.

******************************************************************/

static void CVChooseEta(CVodeMemData* cv_mem, real etaqm1, real etaq, real etaqp1)
{
  real etam;
  
  etam = MAX(etaqm1, MAX(etaq, etaqp1));
  
  if (etam < THRESH) {
    cv_mem_eta = ONE;
    cv_mem_qprime = cv_mem_q;
    return;
  }

  if (etam == etaq) {
    cv_mem_eta = etaq;
    cv_mem_qprime = cv_mem_q;
  } else if (etam == etaqm1) {
    cv_mem_eta = etaqm1;
    cv_mem_qprime = cv_mem_q - 1;
  } else {
    cv_mem_eta = etaqp1;
    cv_mem_qprime = cv_mem_q + 1;
    N_VScale(ONE, cv_mem_acor, cv_mem_zn[cv_mem_qmax]);
  }
}

/****************** CVHandleFailure ******************************

 This routine prints error messages for all cases of failure by
 CVStep. It returns to CVode the value that CVode is to return to
 the user.

*****************************************************************/

static int CVHandleFailure(CVodeMemData* cv_mem, int kflag)
{

  /* Set imxer to the index of maximum weighted local error */
  N_VProd(cv_mem_acor, cv_mem_ewt, cv_mem_tempv);
  N_VAbs(cv_mem_tempv, cv_mem_tempv);

  /* Depending on kflag, print error message and return error flag */
  switch (kflag) {
    case REP_ERR_FAIL:  //fprintf(errfp, MSG_ERR_FAILS, tn, h);
                        return(ERR_FAILURE);
    case REP_CONV_FAIL: //fprintf(errfp, MSG_CONV_FAILS, tn, h);
                        return(CONV_FAILURE);
    case SETUP_FAILED:  //fprintf(errfp, MSG_SETUP_FAILED, tn);
                        return(SETUP_FAILURE);
    case SOLVE_FAILED:  //fprintf(errfp, MSG_SOLVE_FAILED, tn);
                        return(SOLVE_FAILURE);
    default:
//        fprintf(errfp, MSG_SOLVE_FAILED, tn);
                        return(SOLVE_FAILURE);
  }
}

/*******************************************************************/
/********* END Private Helper Functions Implementation *************/
/*******************************************************************/


/***************************************************************/
/************** END CVODE Implementation ***********************/
/***************************************************************/
