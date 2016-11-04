/* QFORM.F -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "../include/f2c.h"

/*<       subroutine qform(m,n,q,ldq,wa) >*/
/* Subroutine */ int qform_(integer *m, integer *n, doublereal *q, integer *
	ldq, doublereal *wa)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;

    /* Local variables */
    doublereal temp;
    integer i__, j, k, l, minmn, jm1, np1;
    doublereal sum;

/*<       integer m,n,ldq >*/
/*<       double precision q(ldq,m),wa(m) >*/
/*     ********** */

/*     subroutine qform */

/*     this subroutine proceeds from the computed qr factorization of */
/*     an m by n matrix a to accumulate the m by m orthogonal matrix */
/*     q from its factored form. */

/*     the subroutine statement is */

/*       subroutine qform(m,n,q,ldq,wa) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of a and the order of q. */

/*       n is a positive integer input variable set to the number */
/*         of columns of a. */

/*       q is an m by m array. on input the full lower trapezoid in */
/*         the first min(m,n) columns of q contains the factored form. */
/*         on output q has been accumulated into a square matrix. */

/*       ldq is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array q. */

/*       wa is a work array of length m. */

/*     subprograms called */

/*       fortran-supplied ... min0 */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
/*<       integer i,j,jm1,k,l,minmn,np1 >*/
/*<       double precision one,sum,temp,zero >*/
/*<       data one,zero /1.0d0,0.0d0/ >*/
    /* Parameter adjustments */
    --wa;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;

    /* Function Body */

/*     zero out upper triangle of q in the first min(m,n) columns. */

/*<       minmn = min0(m,n) >*/
    minmn = min(*m,*n);
/*<       if (minmn .lt. 2) go to 30 >*/
    if (minmn < 2) {
	goto L30;
    }
/*<       do 20 j = 2, minmn >*/
    i__1 = minmn;
    for (j = 2; j <= i__1; ++j) {
/*<          jm1 = j - 1 >*/
	jm1 = j - 1;
/*<          do 10 i = 1, jm1 >*/
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             q(i,j) = zero >*/
	    q[i__ + j * q_dim1] = zero;
/*<    10       continue >*/
/* L10: */
	}
/*<    20    continue >*/
/* L20: */
    }
/*<    30 continue >*/
L30:

/*     initialize remaining columns to those of the identity matrix. */

/*<       np1 = n + 1 >*/
    np1 = *n + 1;
/*<       if (m .lt. np1) go to 60 >*/
    if (*m < np1) {
	goto L60;
    }
/*<       do 50 j = np1, m >*/
    i__1 = *m;
    for (j = np1; j <= i__1; ++j) {
/*<          do 40 i = 1, m >*/
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             q(i,j) = zero >*/
	    q[i__ + j * q_dim1] = zero;
/*<    40       continue >*/
/* L40: */
	}
/*<          q(j,j) = one >*/
	q[j + j * q_dim1] = one;
/*<    50    continue >*/
/* L50: */
    }
/*<    60 continue >*/
L60:

/*     accumulate q from its factored form. */

/*<       do 120 l = 1, minmn >*/
    i__1 = minmn;
    for (l = 1; l <= i__1; ++l) {
/*<          k = minmn - l + 1 >*/
	k = minmn - l + 1;
/*<          do 70 i = k, m >*/
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
/*<             wa(i) = q(i,k) >*/
	    wa[i__] = q[i__ + k * q_dim1];
/*<             q(i,k) = zero >*/
	    q[i__ + k * q_dim1] = zero;
/*<    70       continue >*/
/* L70: */
	}
/*<          q(k,k) = one >*/
	q[k + k * q_dim1] = one;
/*<          if (wa(k) .eq. zero) go to 110 >*/
	if (wa[k] == zero) {
	    goto L110;
	}
/*<          do 100 j = k, m >*/
	i__2 = *m;
	for (j = k; j <= i__2; ++j) {
/*<             sum = zero >*/
	    sum = zero;
/*<             do 80 i = k, m >*/
	    i__3 = *m;
	    for (i__ = k; i__ <= i__3; ++i__) {
/*<                sum = sum + q(i,j)*wa(i) >*/
		sum += q[i__ + j * q_dim1] * wa[i__];
/*<    80          continue >*/
/* L80: */
	    }
/*<             temp = sum/wa(k) >*/
	    temp = sum / wa[k];
/*<             do 90 i = k, m >*/
	    i__3 = *m;
	    for (i__ = k; i__ <= i__3; ++i__) {
/*<                q(i,j) = q(i,j) - temp*wa(i) >*/
		q[i__ + j * q_dim1] -= temp * wa[i__];
/*<    90          continue >*/
/* L90: */
	    }
/*<   100       continue >*/
/* L100: */
	}
/*<   110    continue >*/
L110:
/*<   120    continue >*/
/* L120: */
	;
    }
/*<       return >*/
    return 0;

/*     last card of subroutine qform. */

/*<       end >*/
} /* qform_ */

#ifdef __cplusplus
	}
#endif
