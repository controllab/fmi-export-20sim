/* HYBRD.F -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "../include/f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;

/*<    >*/
/* Subroutine */ int hybrd_(S_fp fcn, integer *n, doublereal *x, doublereal *
	fvec, doublereal *xtol, integer *maxfev, integer *ml, integer *mu, 
	doublereal *epsfcn, doublereal *diag, integer *mode, doublereal *
	factor, integer *nprint, integer *info, integer *nfev, doublereal *
	fjac, integer *ldfjac, doublereal *r__, integer *lr, doublereal *qtf, 
	doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *wa4, void *user_data)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p1 = .1;
    static doublereal p5 = .5;
    static doublereal p001 = .001;
    static doublereal p0001 = 1e-4;
    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    logical sing;
    integer iter;
    doublereal temp;
    integer msum, i__, j, l, iflag;
    doublereal delta;
    extern /* Subroutine */ int qrfac_(integer *, integer *, doublereal *, 
	    integer *, logical *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    logical jeval;
    integer ncsuc;
    doublereal ratio;
    extern doublereal enorm_(integer *, doublereal *);
    doublereal fnorm;
    extern /* Subroutine */ int qform_(integer *, integer *, doublereal *, 
	    integer *, doublereal *), fdjac1_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, void *);
    doublereal pnorm, xnorm, fnorm1;
    extern /* Subroutine */ int r1updt_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, logical *);
    integer nslow1, nslow2;
    extern /* Subroutine */ int r1mpyq_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);
    integer ncfail;
    extern /* Subroutine */ int dogleg_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    doublereal actred, epsmch, prered;
    extern doublereal dpmpar_(integer *);
    integer jm1, iwa[1];
    doublereal sum;

/*<       integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr >*/
/*<       double precision xtol,epsfcn,factor >*/
/*<    >*/
/*<       external fcn >*/
/*     ********** */

/*     subroutine hybrd */

/*     the purpose of hybrd is to find a zero of a system of */
/*     n nonlinear functions in n variables by a modification */
/*     of the powell hybrid method. the user must provide a */
/*     subroutine which calculates the functions. the jacobian is */
/*     then calculated by a forward-difference approximation. */

/*     the subroutine statement is */

/*       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn, */
/*                        diag,mode,factor,nprint,info,nfev,fjac, */
/*                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4) */

/*     where */

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions. fcn must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         subroutine fcn(n,x,fvec,iflag) */
/*         integer n,iflag */
/*         double precision x(n),fvec(n) */
/*         ---------- */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/*         --------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of hybrd. */
/*         in this case set iflag to a negative integer. */

/*       n is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       x is an array of length n. on input x must contain */
/*         an initial estimate of the solution vector. on output x */
/*         contains the final estimate of the solution vector. */

/*       fvec is an output array of length n which contains */
/*         the functions evaluated at the output x. */

/*       xtol is a nonnegative input variable. termination */
/*         occurs when the relative error between two consecutive */
/*         iterates is at most xtol. */

/*       maxfev is a positive integer input variable. termination */
/*         occurs when the number of calls to fcn is at least maxfev */
/*         by the end of an iteration. */

/*       ml is a nonnegative integer input variable which specifies */
/*         the number of subdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         ml to at least n - 1. */

/*       mu is a nonnegative integer input variable which specifies */
/*         the number of superdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         mu to at least n - 1. */

/*       epsfcn is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. this */
/*         approximation assumes that the relative errors in the */
/*         functions are of the order of epsfcn. if epsfcn is less */
/*         than the machine precision, it is assumed that the relative */
/*         errors in the functions are of the order of the machine */
/*         precision. */

/*       diag is an array of length n. if mode = 1 (see */
/*         below), diag is internally set. if mode = 2, diag */
/*         must contain positive entries that serve as */
/*         multiplicative scale factors for the variables. */

/*       mode is an integer input variable. if mode = 1, the */
/*         variables will be scaled internally. if mode = 2, */
/*         the scaling is specified by the input diag. other */
/*         values of mode are equivalent to mode = 1. */

/*       factor is a positive input variable used in determining the */
/*         initial step bound. this bound is set to the product of */
/*         factor and the euclidean norm of diag*x if nonzero, or else */
/*         to factor itself. in most cases factor should lie in the */
/*         interval (.1,100.). 100. is a generally recommended value. */

/*       nprint is an integer input variable that enables controlled */
/*         printing of iterates if it is positive. in this case, */
/*         fcn is called with iflag = 0 at the beginning of the first */
/*         iteration and every nprint iterations thereafter and */
/*         immediately prior to return, with x and fvec available */
/*         for printing. if nprint is not positive, no special calls */
/*         of fcn with iflag = 0 are made. */

/*       info is an integer output variable. if the user has */
/*         terminated execution, info is set to the (negative) */
/*         value of iflag. see description of fcn. otherwise, */
/*         info is set as follows. */

/*         info = 0   improper input parameters. */

/*         info = 1   relative error between two consecutive iterates */
/*                    is at most xtol. */

/*         info = 2   number of calls to fcn has reached or exceeded */
/*                    maxfev. */

/*         info = 3   xtol is too small. no further improvement in */
/*                    the approximate solution x is possible. */

/*         info = 4   iteration is not making good progress, as */
/*                    measured by the improvement from the last */
/*                    five jacobian evaluations. */

/*         info = 5   iteration is not making good progress, as */
/*                    measured by the improvement from the last */
/*                    ten iterations. */

/*       nfev is an integer output variable set to the number of */
/*         calls to fcn. */

/*       fjac is an output n by n array which contains the */
/*         orthogonal matrix q produced by the qr factorization */
/*         of the final approximate jacobian. */

/*       ldfjac is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array fjac. */

/*       r is an output array of length lr which contains the */
/*         upper triangular matrix produced by the qr factorization */
/*         of the final approximate jacobian, stored rowwise. */

/*       lr is a positive integer input variable not less than */
/*         (n*(n+1))/2. */

/*       qtf is an output array of length n which contains */
/*         the vector (q transpose)*fvec. */

/*       wa1, wa2, wa3, and wa4 are work arrays of length n. */

/*     subprograms called */

/*       user-supplied ...... fcn */

/*       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1, */
/*                            qform,qrfac,r1mpyq,r1updt */

/*       fortran-supplied ... dabs,dmax1,dmin1,min0,mod */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
/*<       integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2 >*/
/*<       integer iwa(1) >*/
/*<       logical jeval,sing >*/
/*<    >*/
/*<       double precision dpmpar,enorm >*/
/*<    >*/
    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --diag;
    --fvec;
    --x;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;
    --r__;

    /* Function Body */

	/* grn: 2007-03-21. xnorm initialiation, prevent compiler warning about use before initialization */
	xnorm = 0.0;

/*     epsmch is the machine precision. */

/*<       epsmch = dpmpar(1) >*/
    epsmch = dpmpar_(&c__1);

/*<       info = 0 >*/
    *info = 0;
/*<       iflag = 0 >*/
    iflag = 0;
/*<       nfev = 0 >*/
    *nfev = 0;

/*     check the input parameters for errors. */

/*<    >*/
    if (*n <= 0 || *xtol < zero || *maxfev <= 0 || *ml < 0 || *mu < 0 || *
	    factor <= zero || *ldfjac < *n || *lr < *n * (*n + 1) / 2) {
	goto L300;
    }
/*<       if (mode .ne. 2) go to 20 >*/
    if (*mode != 2) {
	goto L20;
    }
/*<       do 10 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          if (diag(j) .le. zero) go to 300 >*/
	if (diag[j] <= zero) {
	    goto L300;
	}
/*<    10    continue >*/
/* L10: */
    }
/*<    20 continue >*/
L20:

/*     evaluate the function at the starting point */
/*     and calculate its norm. */

/*<       iflag = 1 >*/
    iflag = 1;
/*<       call fcn(n,x,fvec,iflag) >*/
    (*fcn)(n, &x[1], &fvec[1], &iflag, user_data);
/*<       nfev = 1 >*/
    *nfev = 1;
/*<       if (iflag .lt. 0) go to 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<       fnorm = enorm(n,fvec) >*/
    fnorm = enorm_(n, &fvec[1]);

/*     determine the number of calls to fcn needed to compute */
/*     the jacobian matrix. */

/*<       msum = min0(ml+mu+1,n) >*/
/* Computing MIN */
    i__1 = *ml + *mu + 1;
    msum = min(i__1,*n);

/*     initialize iteration counter and monitors. */

/*<       iter = 1 >*/
    iter = 1;
/*<       ncsuc = 0 >*/
    ncsuc = 0;
/*<       ncfail = 0 >*/
    ncfail = 0;
/*<       nslow1 = 0 >*/
    nslow1 = 0;
/*<       nslow2 = 0 >*/
    nslow2 = 0;

/*     beginning of the outer loop. */

/*<    30 continue >*/
L30:
/*<          jeval = .true. >*/
    jeval = TRUE_;

/*        calculate the jacobian matrix. */

/*<          iflag = 2 >*/
    iflag = 2;
/*<    >*/
    fdjac1_((S_fp)fcn, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag,
	     ml, mu, epsfcn, &wa1[1], &wa2[1], user_data);
/*<          nfev = nfev + msum >*/
    *nfev += msum;
/*<          if (iflag .lt. 0) go to 300 >*/
    if (iflag < 0) {
	goto L300;
    }

/*        compute the qr factorization of the jacobian. */

/*<          call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3) >*/
    qrfac_(n, n, &fjac[fjac_offset], ldfjac, &c_false, iwa, &c__1, &wa1[1], &
	    wa2[1], &wa3[1]);

/*        on the first iteration and if mode is 1, scale according */
/*        to the norms of the columns of the initial jacobian. */

/*<          if (iter .ne. 1) go to 70 >*/
    if (iter != 1) {
	goto L70;
    }
/*<          if (mode .eq. 2) go to 50 >*/
    if (*mode == 2) {
	goto L50;
    }
/*<          do 40 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             diag(j) = wa2(j) >*/
	diag[j] = wa2[j];
/*<             if (wa2(j) .eq. zero) diag(j) = one >*/
	if (wa2[j] == zero) {
	    diag[j] = one;
	}
/*<    40       continue >*/
/* L40: */
    }
/*<    50    continue >*/
L50:

/*        on the first iteration, calculate the norm of the scaled x */
/*        and initialize the step bound delta. */

/*<          do 60 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             wa3(j) = diag(j)*x(j) >*/
	wa3[j] = diag[j] * x[j];
/*<    60       continue >*/
/* L60: */
    }
/*<          xnorm = enorm(n,wa3) >*/
    xnorm = enorm_(n, &wa3[1]);
/*<          delta = factor*xnorm >*/
    delta = *factor * xnorm;
/*<          if (delta .eq. zero) delta = factor >*/
    if (delta == zero) {
	delta = *factor;
    }
/*<    70    continue >*/
L70:

/*        form (q transpose)*fvec and store in qtf. */

/*<          do 80 i = 1, n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<             qtf(i) = fvec(i) >*/
	qtf[i__] = fvec[i__];
/*<    80       continue >*/
/* L80: */
    }
/*<          do 120 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             if (fjac(j,j) .eq. zero) go to 110 >*/
	if (fjac[j + j * fjac_dim1] == zero) {
	    goto L110;
	}
/*<             sum = zero >*/
	sum = zero;
/*<             do 90 i = j, n >*/
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<                sum = sum + fjac(i,j)*qtf(i) >*/
	    sum += fjac[i__ + j * fjac_dim1] * qtf[i__];
/*<    90          continue >*/
/* L90: */
	}
/*<             temp = -sum/fjac(j,j) >*/
	temp = -sum / fjac[j + j * fjac_dim1];
/*<             do 100 i = j, n >*/
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<                qtf(i) = qtf(i) + fjac(i,j)*temp >*/
	    qtf[i__] += fjac[i__ + j * fjac_dim1] * temp;
/*<   100          continue >*/
/* L100: */
	}
/*<   110       continue >*/
L110:
/*<   120       continue >*/
/* L120: */
	;
    }

/*        copy the triangular factor of the qr factorization into r. */

/*<          sing = .false. >*/
    sing = FALSE_;
/*<          do 150 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             l = j >*/
	l = j;
/*<             jm1 = j - 1 >*/
	jm1 = j - 1;
/*<             if (jm1 .lt. 1) go to 140 >*/
	if (jm1 < 1) {
	    goto L140;
	}
/*<             do 130 i = 1, jm1 >*/
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<                r(l) = fjac(i,j) >*/
	    r__[l] = fjac[i__ + j * fjac_dim1];
/*<                l = l + n - i >*/
	    l = l + *n - i__;
/*<   130          continue >*/
/* L130: */
	}
/*<   140       continue >*/
L140:
/*<             r(l) = wa1(j) >*/
	r__[l] = wa1[j];
/*<             if (wa1(j) .eq. zero) sing = .true. >*/
	if (wa1[j] == zero) {
	    sing = TRUE_;
	}
/*<   150       continue >*/
/* L150: */
    }

/*        accumulate the orthogonal factor in fjac. */

/*<          call qform(n,n,fjac,ldfjac,wa1) >*/
    qform_(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);

/*        rescale if necessary. */

/*<          if (mode .eq. 2) go to 170 >*/
    if (*mode == 2) {
	goto L170;
    }
/*<          do 160 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             diag(j) = dmax1(diag(j),wa2(j)) >*/
/* Computing MAX */
	d__1 = diag[j], d__2 = wa2[j];
	diag[j] = max(d__1,d__2);
/*<   160       continue >*/
/* L160: */
    }
/*<   170    continue >*/
L170:

/*        beginning of the inner loop. */

/*<   180    continue >*/
L180:

/*           if requested, call fcn to enable printing of iterates. */

/*<             if (nprint .le. 0) go to 190 >*/
    if (*nprint <= 0) {
	goto L190;
    }
/*<             iflag = 0 >*/
    iflag = 0;
/*<             if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag) >*/
    if ((iter - 1) % *nprint == 0) {
	(*fcn)(n, &x[1], &fvec[1], &iflag, user_data);
    }
/*<             if (iflag .lt. 0) go to 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<   190       continue >*/
L190:

/*           determine the direction p. */

/*<             call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3) >*/
    dogleg_(n, &r__[1], lr, &diag[1], &qtf[1], &delta, &wa1[1], &wa2[1], &wa3[
	    1]);

/*           store the direction p and x + p. calculate the norm of p. */

/*<             do 200 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<                wa1(j) = -wa1(j) >*/
	wa1[j] = -wa1[j];
/*<                wa2(j) = x(j) + wa1(j) >*/
	wa2[j] = x[j] + wa1[j];
/*<                wa3(j) = diag(j)*wa1(j) >*/
	wa3[j] = diag[j] * wa1[j];
/*<   200          continue >*/
/* L200: */
    }
/*<             pnorm = enorm(n,wa3) >*/
    pnorm = enorm_(n, &wa3[1]);

/*           on the first iteration, adjust the initial step bound. */

/*<             if (iter .eq. 1) delta = dmin1(delta,pnorm) >*/
    if (iter == 1) {
	delta = min(delta,pnorm);
    }

/*           evaluate the function at x + p and calculate its norm. */

/*<             iflag = 1 >*/
    iflag = 1;
/*<             call fcn(n,wa2,wa4,iflag) >*/
    (*fcn)(n, &wa2[1], &wa4[1], &iflag, user_data);
/*<             nfev = nfev + 1 >*/
    ++(*nfev);
/*<             if (iflag .lt. 0) go to 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<             fnorm1 = enorm(n,wa4) >*/
    fnorm1 = enorm_(n, &wa4[1]);

/*           compute the scaled actual reduction. */

/*<             actred = -one >*/
    actred = -one;
/*<             if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2 >*/
    if (fnorm1 < fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / fnorm;
	actred = one - d__1 * d__1;
    }

/*           compute the scaled predicted reduction. */

/*<             l = 1 >*/
    l = 1;
/*<             do 220 i = 1, n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<                sum = zero >*/
	sum = zero;
/*<                do 210 j = i, n >*/
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
/*<                   sum = sum + r(l)*wa1(j) >*/
	    sum += r__[l] * wa1[j];
/*<                   l = l + 1 >*/
	    ++l;
/*<   210             continue >*/
/* L210: */
	}
/*<                wa3(i) = qtf(i) + sum >*/
	wa3[i__] = qtf[i__] + sum;
/*<   220          continue >*/
/* L220: */
    }
/*<             temp = enorm(n,wa3) >*/
    temp = enorm_(n, &wa3[1]);
/*<             prered = zero >*/
    prered = zero;
/*<             if (temp .lt. fnorm) prered = one - (temp/fnorm)**2 >*/
    if (temp < fnorm) {
/* Computing 2nd power */
	d__1 = temp / fnorm;
	prered = one - d__1 * d__1;
    }

/*           compute the ratio of the actual to the predicted */
/*           reduction. */

/*<             ratio = zero >*/
    ratio = zero;
/*<             if (prered .gt. zero) ratio = actred/prered >*/
    if (prered > zero) {
	ratio = actred / prered;
    }

/*           update the step bound. */

/*<             if (ratio .ge. p1) go to 230 >*/
    if (ratio >= p1) {
	goto L230;
    }
/*<                ncsuc = 0 >*/
    ncsuc = 0;
/*<                ncfail = ncfail + 1 >*/
    ++ncfail;
/*<                delta = p5*delta >*/
    delta = p5 * delta;
/*<                go to 240 >*/
    goto L240;
/*<   230       continue >*/
L230:
/*<                ncfail = 0 >*/
    ncfail = 0;
/*<                ncsuc = ncsuc + 1 >*/
    ++ncsuc;
/*<    >*/
    if (ratio >= p5 || ncsuc > 1) {
/* Computing MAX */
	d__1 = delta, d__2 = pnorm / p5;
	delta = max(d__1,d__2);
    }
/*<                if (dabs(ratio-one) .le. p1) delta = pnorm/p5 >*/
    if ((d__1 = ratio - one, abs(d__1)) <= p1) {
	delta = pnorm / p5;
    }
/*<   240       continue >*/
L240:

/*           test for successful iteration. */

/*<             if (ratio .lt. p0001) go to 260 >*/
    if (ratio < p0001) {
	goto L260;
    }

/*           successful iteration. update x, fvec, and their norms. */

/*<             do 250 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<                x(j) = wa2(j) >*/
	x[j] = wa2[j];
/*<                wa2(j) = diag(j)*x(j) >*/
	wa2[j] = diag[j] * x[j];
/*<                fvec(j) = wa4(j) >*/
	fvec[j] = wa4[j];
/*<   250          continue >*/
/* L250: */
    }
/*<             xnorm = enorm(n,wa2) >*/
    xnorm = enorm_(n, &wa2[1]);
/*<             fnorm = fnorm1 >*/
    fnorm = fnorm1;
/*<             iter = iter + 1 >*/
    ++iter;
/*<   260       continue >*/
L260:

/*           determine the progress of the iteration. */

/*<             nslow1 = nslow1 + 1 >*/
    ++nslow1;
/*<             if (actred .ge. p001) nslow1 = 0 >*/
    if (actred >= p001) {
	nslow1 = 0;
    }
/*<             if (jeval) nslow2 = nslow2 + 1 >*/
    if (jeval) {
	++nslow2;
    }
/*<             if (actred .ge. p1) nslow2 = 0 >*/
    if (actred >= p1) {
	nslow2 = 0;
    }

/*           test for convergence. */

/*<             if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1 >*/
    if (delta <= *xtol * xnorm || fnorm == zero) {
	*info = 1;
    }
/*<             if (info .ne. 0) go to 300 >*/
    if (*info != 0) {
	goto L300;
    }

/*           tests for termination and stringent tolerances. */

/*<             if (nfev .ge. maxfev) info = 2 >*/
    if (*nfev >= *maxfev) {
	*info = 2;
    }
/*<             if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3 >*/
/* Computing MAX */
    d__1 = p1 * delta;
    if (p1 * max(d__1,pnorm) <= epsmch * xnorm) {
	*info = 3;
    }
/*<             if (nslow2 .eq. 5) info = 4 >*/
    if (nslow2 == 500) { /* changed by grn from 5 to 500 */
	*info = 4;
    }
/*<             if (nslow1 .eq. 10) info = 5 >*/
    if (nslow1 == 1000) { /* changed by grn from 5 to 500 */
	*info = 5;
    }
/*<             if (info .ne. 0) go to 300 >*/
    if (*info != 0) {
	goto L300;
    }

/*           criterion for recalculating jacobian approximation */
/*           by forward differences. */

/*<             if (ncfail .eq. 2) go to 290 >*/
    if (ncfail == 2) {
	goto L290;
    }

/*           calculate the rank one modification to the jacobian */
/*           and update qtf if necessary. */

/*<             do 280 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<                sum = zero >*/
	sum = zero;
/*<                do 270 i = 1, n >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<                   sum = sum + fjac(i,j)*wa4(i) >*/
	    sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
/*<   270             continue >*/
/* L270: */
	}
/*<                wa2(j) = (sum - wa3(j))/pnorm >*/
	wa2[j] = (sum - wa3[j]) / pnorm;
/*<                wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm) >*/
	wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
/*<                if (ratio .ge. p0001) qtf(j) = sum >*/
	if (ratio >= p0001) {
	    qtf[j] = sum;
	}
/*<   280          continue >*/
/* L280: */
    }

/*           compute the qr factorization of the updated jacobian. */

/*<             call r1updt(n,n,r,lr,wa1,wa2,wa3,sing) >*/
    r1updt_(n, n, &r__[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing);
/*<             call r1mpyq(n,n,fjac,ldfjac,wa2,wa3) >*/
    r1mpyq_(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
/*<             call r1mpyq(1,n,qtf,1,wa2,wa3) >*/
    r1mpyq_(&c__1, n, &qtf[1], &c__1, &wa2[1], &wa3[1]);

/*           end of the inner loop. */

/*<             jeval = .false. >*/
    jeval = FALSE_;
/*<             go to 180 >*/
    goto L180;
/*<   290    continue >*/
L290:

/*        end of the outer loop. */

/*<          go to 30 >*/
    goto L30;
/*<   300 continue >*/
L300:

/*     termination, either normal or user imposed. */

/*<       if (iflag .lt. 0) info = iflag >*/
    if (iflag < 0) {
	*info = iflag;
    }
/*<       iflag = 0 >*/
    iflag = 0;
/*<       if (nprint .gt. 0) call fcn(n,x,fvec,iflag) >*/
    if (*nprint > 0) {
	(*fcn)(n, &x[1], &fvec[1], &iflag, user_data);
    }
/*<       return >*/
    return 0;

/*     last card of subroutine hybrd. */

/*<       end >*/
} /* hybrd_ */

#ifdef __cplusplus
	}
#endif

