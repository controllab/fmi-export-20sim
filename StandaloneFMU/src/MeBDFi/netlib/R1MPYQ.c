/* R1MPYQ.F -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "../include/f2c.h"

/*<       subroutine r1mpyq(m,n,a,lda,v,w) >*/
/* Subroutine */ int r1mpyq_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *v, doublereal *w)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal temp;
    integer i__, j, nm1, nmj;
    doublereal cos__, sin__;

/*<       integer m,n,lda >*/
/*<       double precision a(lda,n),v(n),w(n) >*/
/*     ********** */

/*     subroutine r1mpyq */

/*     given an m by n matrix a, this subroutine computes a*q where */
/*     q is the product of 2*(n - 1) transformations */

/*           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) */

/*     and gv(i), gw(i) are givens rotations in the (i,n) plane which */
/*     eliminate elements in the i-th and n-th planes, respectively. */
/*     q itself is not given, rather the information to recover the */
/*     gv, gw rotations is supplied. */

/*     the subroutine statement is */

/*       subroutine r1mpyq(m,n,a,lda,v,w) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of a. */

/*       n is a positive integer input variable set to the number */
/*         of columns of a. */

/*       a is an m by n array. on input a must contain the matrix */
/*         to be postmultiplied by the orthogonal matrix q */
/*         described above. on output a*q has replaced a. */

/*       lda is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array a. */

/*       v is an input array of length n. v(i) must contain the */
/*         information necessary to recover the givens rotation gv(i) */
/*         described above. */

/*       w is an input array of length n. w(i) must contain the */
/*         information necessary to recover the givens rotation gw(i) */
/*         described above. */

/*     subroutines called */

/*       fortran-supplied ... dabs,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
/*<       integer i,j,nmj,nm1 >*/
/*<       double precision cos,one,sin,temp >*/
/*<       data one /1.0d0/ >*/
    /* Parameter adjustments */
    --w;
    --v;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
	cos__ = 0.0; /* 2007-03-21. grn, potentially uninitialized variable used some lines lower */
    sin__ = 0.0; /* 2007-03-21. grn, potentially uninitialized variable used some lines lower */
	
	/* Function Body */

/*     apply the first set of givens rotations to a. */

/*<       nm1 = n - 1 >*/
    nm1 = *n - 1;
/*<       if (nm1 .lt. 1) go to 50 >*/
    if (nm1 < 1) {
	goto L50;
    }
/*<       do 20 nmj = 1, nm1 >*/
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
/*<          j = n - nmj >*/
	j = *n - nmj;
/*<          if (dabs(v(j)) .gt. one) cos = one/v(j) >*/
	if ((d__1 = v[j], abs(d__1)) > one) {
	    cos__ = one / v[j];
	}
/*<          if (dabs(v(j)) .gt. one) sin = dsqrt(one-cos**2) >*/
	if ((d__1 = v[j], abs(d__1)) > one) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(one - d__2 * d__2);
	}
/*<          if (dabs(v(j)) .le. one) sin = v(j) >*/
	if ((d__1 = v[j], abs(d__1)) <= one) {
	    sin__ = v[j];
	}
/*<          if (dabs(v(j)) .le. one) cos = dsqrt(one-sin**2) >*/
	if ((d__1 = v[j], abs(d__1)) <= one) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(one - d__2 * d__2);
	}
/*<          do 10 i = 1, m >*/
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             temp = cos*a(i,j) - sin*a(i,n) >*/
	    temp = cos__ * a[i__ + j * a_dim1] - sin__ * a[i__ + *n * a_dim1];
/*<             a(i,n) = sin*a(i,j) + cos*a(i,n) >*/
	    a[i__ + *n * a_dim1] = sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + *n * a_dim1];
/*<             a(i,j) = temp >*/
	    a[i__ + j * a_dim1] = temp;
/*<    10       continue >*/
/* L10: */
	}
/*<    20    continue >*/
/* L20: */
    }

/*     apply the second set of givens rotations to a. */

/*<       do 40 j = 1, nm1 >*/
    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          if (dabs(w(j)) .gt. one) cos = one/w(j) >*/
	if ((d__1 = w[j], abs(d__1)) > one) {
	    cos__ = one / w[j];
	}
/*<          if (dabs(w(j)) .gt. one) sin = dsqrt(one-cos**2) >*/
	if ((d__1 = w[j], abs(d__1)) > one) {
/* Computing 2nd power */
	    d__2 = cos__;
	    sin__ = sqrt(one - d__2 * d__2);
	}
/*<          if (dabs(w(j)) .le. one) sin = w(j) >*/
	if ((d__1 = w[j], abs(d__1)) <= one) {
	    sin__ = w[j];
	}
/*<          if (dabs(w(j)) .le. one) cos = dsqrt(one-sin**2) >*/
	if ((d__1 = w[j], abs(d__1)) <= one) {
/* Computing 2nd power */
	    d__2 = sin__;
	    cos__ = sqrt(one - d__2 * d__2);
	}
/*<          do 30 i = 1, m >*/
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             temp = cos*a(i,j) + sin*a(i,n) >*/
	    temp = cos__ * a[i__ + j * a_dim1] + sin__ * a[i__ + *n * a_dim1];
/*<             a(i,n) = -sin*a(i,j) + cos*a(i,n) >*/
	    a[i__ + *n * a_dim1] = -sin__ * a[i__ + j * a_dim1] + cos__ * a[
		    i__ + *n * a_dim1];
/*<             a(i,j) = temp >*/
	    a[i__ + j * a_dim1] = temp;
/*<    30       continue >*/
/* L30: */
	}
/*<    40    continue >*/
/* L40: */
    }
/*<    50 continue >*/
L50:
/*<       return >*/
    return 0;

/*     last card of subroutine r1mpyq. */

/*<       end >*/
} /* r1mpyq_ */

#ifdef __cplusplus
	}
#endif
