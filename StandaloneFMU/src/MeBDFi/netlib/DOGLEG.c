/* DOGLEG.F -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "../include/f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/*<       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2) >*/
/* Subroutine */ int dogleg_(integer *n, doublereal *r__, integer *lr, 
	doublereal *diag, doublereal *qtb, doublereal *delta, doublereal *x, 
	doublereal *wa1, doublereal *wa2)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal temp;
    integer i__, j, k, l;
    doublereal alpha, bnorm;
    extern doublereal enorm_(integer *, doublereal *);
    doublereal gnorm, qnorm;
    integer jj;
    doublereal epsmch;
    extern doublereal dpmpar_(integer *);
    doublereal sgnorm;
    integer jp1;
    doublereal sum;

/*<       integer n,lr >*/
/*<       double precision delta >*/
/*<       double precision r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n) >*/
/*     ********** */

/*     subroutine dogleg */

/*     given an m by n matrix a, an n by n nonsingular diagonal */
/*     matrix d, an m-vector b, and a positive number delta, the */
/*     problem is to determine the convex combination x of the */
/*     gauss-newton and scaled gradient directions that minimizes */
/*     (a*x - b) in the least squares sense, subject to the */
/*     restriction that the euclidean norm of d*x be at most delta. */

/*     this subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     qr factorization of a. that is, if a = q*r, where q has */
/*     orthogonal columns and r is an upper triangular matrix, */
/*     then dogleg expects the full upper triangle of r and */
/*     the first n components of (q transpose)*b. */

/*     the subroutine statement is */

/*       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2) */

/*     where */

/*       n is a positive integer input variable set to the order of r. */

/*       r is an input array of length lr which must contain the upper */
/*         triangular matrix r stored by rows. */

/*       lr is a positive integer input variable not less than */
/*         (n*(n+1))/2. */

/*       diag is an input array of length n which must contain the */
/*         diagonal elements of the matrix d. */

/*       qtb is an input array of length n which must contain the first */
/*         n elements of the vector (q transpose)*b. */

/*       delta is a positive input variable which specifies an upper */
/*         bound on the euclidean norm of d*x. */

/*       x is an output array of length n which contains the desired */
/*         convex combination of the gauss-newton direction and the */
/*         scaled gradient direction. */

/*       wa1 and wa2 are work arrays of length n. */

/*     subprograms called */

/*       minpack-supplied ... dpmpar,enorm */

/*       fortran-supplied ... dabs,dmax1,dmin1,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
/*<       integer i,j,jj,jp1,k,l >*/
/*<    >*/
/*<       double precision dpmpar,enorm >*/
/*<       data one,zero /1.0d0,0.0d0/ >*/
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --x;
    --qtb;
    --diag;
    --r__;

    /* Function Body */

/*     epsmch is the machine precision. */

/*<       epsmch = dpmpar(1) >*/
    epsmch = dpmpar_(&c__1);

/*     first, calculate the gauss-newton direction. */

/*<       jj = (n*(n + 1))/2 + 1 >*/
    jj = *n * (*n + 1) / 2 + 1;
/*<       do 50 k = 1, n >*/
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*<          j = n - k + 1 >*/
	j = *n - k + 1;
/*<          jp1 = j + 1 >*/
	jp1 = j + 1;
/*<          jj = jj - k >*/
	jj -= k;
/*<          l = jj + 1 >*/
	l = jj + 1;
/*<          sum = zero >*/
	sum = zero;
/*<          if (n .lt. jp1) go to 20 >*/
	if (*n < jp1) {
	    goto L20;
	}
/*<          do 10 i = jp1, n >*/
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
/*<             sum = sum + r(l)*x(i) >*/
	    sum += r__[l] * x[i__];
/*<             l = l + 1 >*/
	    ++l;
/*<    10       continue >*/
/* L10: */
	}
/*<    20    continue >*/
L20:
/*<          temp = r(jj) >*/
	temp = r__[jj];
/*<          if (temp .ne. zero) go to 40 >*/
	if (temp != zero) {
	    goto L40;
	}
/*<          l = j >*/
	l = j;
/*<          do 30 i = 1, j >*/
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             temp = dmax1(temp,dabs(r(l))) >*/
/* Computing MAX */
	    d__2 = temp, d__3 = (d__1 = r__[l], abs(d__1));
	    temp = max(d__2,d__3);
/*<             l = l + n - i >*/
	    l = l + *n - i__;
/*<    30       continue >*/
/* L30: */
	}
/*<          temp = epsmch*temp >*/
	temp = epsmch * temp;
/*<          if (temp .eq. zero) temp = epsmch >*/
	if (temp == zero) {
	    temp = epsmch;
	}
/*<    40    continue >*/
L40:
/*<          x(j) = (qtb(j) - sum)/temp >*/
	x[j] = (qtb[j] - sum) / temp;
/*<    50    continue >*/
/* L50: */
    }

/*     test whether the gauss-newton direction is acceptable. */

/*<       do 60 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          wa1(j) = zero >*/
	wa1[j] = zero;
/*<          wa2(j) = diag(j)*x(j) >*/
	wa2[j] = diag[j] * x[j];
/*<    60    continue >*/
/* L60: */
    }
/*<       qnorm = enorm(n,wa2) >*/
    qnorm = enorm_(n, &wa2[1]);
/*<       if (qnorm .le. delta) go to 140 >*/
    if (qnorm <= *delta) {
	goto L140;
    }

/*     the gauss-newton direction is not acceptable. */
/*     next, calculate the scaled gradient direction. */

/*<       l = 1 >*/
    l = 1;
/*<       do 80 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          temp = qtb(j) >*/
	temp = qtb[j];
/*<          do 70 i = j, n >*/
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<             wa1(i) = wa1(i) + r(l)*temp >*/
	    wa1[i__] += r__[l] * temp;
/*<             l = l + 1 >*/
	    ++l;
/*<    70       continue >*/
/* L70: */
	}
/*<          wa1(j) = wa1(j)/diag(j) >*/
	wa1[j] /= diag[j];
/*<    80    continue >*/
/* L80: */
    }

/*     calculate the norm of the scaled gradient and test for */
/*     the special case in which the scaled gradient is zero. */

/*<       gnorm = enorm(n,wa1) >*/
    gnorm = enorm_(n, &wa1[1]);
/*<       sgnorm = zero >*/
    sgnorm = zero;
/*<       alpha = delta/qnorm >*/
    alpha = *delta / qnorm;
/*<       if (gnorm .eq. zero) go to 120 >*/
    if (gnorm == zero) {
	goto L120;
    }

/*     calculate the point along the scaled gradient */
/*     at which the quadratic is minimized. */

/*<       do 90 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          wa1(j) = (wa1(j)/gnorm)/diag(j) >*/
	wa1[j] = wa1[j] / gnorm / diag[j];
/*<    90    continue >*/
/* L90: */
    }
/*<       l = 1 >*/
    l = 1;
/*<       do 110 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          sum = zero >*/
	sum = zero;
/*<          do 100 i = j, n >*/
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<             sum = sum + r(l)*wa1(i) >*/
	    sum += r__[l] * wa1[i__];
/*<             l = l + 1 >*/
	    ++l;
/*<   100       continue >*/
/* L100: */
	}
/*<          wa2(j) = sum >*/
	wa2[j] = sum;
/*<   110    continue >*/
/* L110: */
    }
/*<       temp = enorm(n,wa2) >*/
    temp = enorm_(n, &wa2[1]);
/*<       sgnorm = (gnorm/temp)/temp >*/
    sgnorm = gnorm / temp / temp;

/*     test whether the scaled gradient direction is acceptable. */

/*<       alpha = zero >*/
    alpha = zero;
/*<       if (sgnorm .ge. delta) go to 120 >*/
    if (sgnorm >= *delta) {
	goto L120;
    }

/*     the scaled gradient direction is not acceptable. */
/*     finally, calculate the point along the dogleg */
/*     at which the quadratic is minimized. */

/*<       bnorm = enorm(n,qtb) >*/
    bnorm = enorm_(n, &qtb[1]);
/*<       temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta) >*/
    temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / *delta);
/*<    >*/
/* Computing 2nd power */
    d__1 = sgnorm / *delta;
/* Computing 2nd power */
    d__2 = temp - *delta / qnorm;
/* Computing 2nd power */
    d__3 = *delta / qnorm;
/* Computing 2nd power */
    d__4 = sgnorm / *delta;
    temp = temp - *delta / qnorm * (d__1 * d__1) + sqrt(d__2 * d__2 + (one - 
	    d__3 * d__3) * (one - d__4 * d__4));
/*<       alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp >*/
/* Computing 2nd power */
    d__1 = sgnorm / *delta;
    alpha = *delta / qnorm * (one - d__1 * d__1) / temp;
/*<   120 continue >*/
L120:

/*     form appropriate convex combination of the gauss-newton */
/*     direction and the scaled gradient direction. */

/*<       temp = (one - alpha)*dmin1(sgnorm,delta) >*/
    temp = (one - alpha) * min(sgnorm,*delta);
/*<       do 130 j = 1, n >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          x(j) = temp*wa1(j) + alpha*x(j) >*/
	x[j] = temp * wa1[j] + alpha * x[j];
/*<   130    continue >*/
/* L130: */
    }
/*<   140 continue >*/
L140:
/*<       return >*/
    return 0;

/*     last card of subroutine dogleg. */

/*<       end >*/
} /* dogleg_ */

#ifdef __cplusplus
	}
#endif
