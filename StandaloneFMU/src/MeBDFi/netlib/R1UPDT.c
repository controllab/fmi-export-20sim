/* R1UPDT.F -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "../include/f2c.h"

/* Table of constant values */

static integer c__3 = 3;

/*<       subroutine r1updt(m,n,s,ls,u,v,w,sing) >*/
/* Subroutine */ int r1updt_(integer *m, integer *n, doublereal *s, integer *
	ls, doublereal *u, doublereal *v, doublereal *w, logical *sing)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p5 = .5;
    static doublereal p25 = .25;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal temp;
    integer i__, j, l;
    doublereal giant, cotan;
    integer jj;
    extern doublereal dpmpar_(integer *);
    integer nm1;
    doublereal tan__;
    integer nmj;
    doublereal cos__, sin__, tau;

/*<       integer m,n,ls >*/
/*<       logical sing >*/
/*<       double precision s(ls),u(m),v(n),w(m) >*/
/*     ********** */

/*     subroutine r1updt */

/*     given an m by n lower trapezoidal matrix s, an m-vector u, */
/*     and an n-vector v, the problem is to determine an */
/*     orthogonal matrix q such that */

/*                   t */
/*           (s + u*v )*q */

/*     is again lower trapezoidal. */

/*     this subroutine determines q as the product of 2*(n - 1) */
/*     transformations */

/*           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) */

/*     where gv(i), gw(i) are givens rotations in the (i,n) plane */
/*     which eliminate elements in the i-th and n-th planes, */
/*     respectively. q itself is not accumulated, rather the */
/*     information to recover the gv, gw rotations is returned. */

/*     the subroutine statement is */

/*       subroutine r1updt(m,n,s,ls,u,v,w,sing) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of s. */

/*       n is a positive integer input variable set to the number */
/*         of columns of s. n must not exceed m. */

/*       s is an array of length ls. on input s must contain the lower */
/*         trapezoidal matrix s stored by columns. on output s contains */
/*         the lower trapezoidal matrix produced as described above. */

/*       ls is a positive integer input variable not less than */
/*         (n*(2*m-n+1))/2. */

/*       u is an input array of length m which must contain the */
/*         vector u. */

/*       v is an array of length n. on input v must contain the vector */
/*         v. on output v(i) contains the information necessary to */
/*         recover the givens rotation gv(i) described above. */

/*       w is an output array of length m. w(i) contains information */
/*         necessary to recover the givens rotation gw(i) described */
/*         above. */

/*       sing is a logical output variable. sing is set true if any */
/*         of the diagonal elements of the output s are zero. otherwise */
/*         sing is set false. */

/*     subprograms called */

/*       minpack-supplied ... dpmpar */

/*       fortran-supplied ... dabs,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more, */
/*     john l. nazareth */

/*     ********** */
/*<       integer i,j,jj,l,nmj,nm1 >*/
/*<    >*/
/*<       double precision dpmpar >*/
/*<       data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/ >*/
    /* Parameter adjustments */
    --w;
    --u;
    --v;
    --s;

    /* Function Body */

/*     giant is the largest magnitude. */

/*<       giant = dpmpar(3) >*/
    giant = dpmpar_(&c__3);

/*     initialize the diagonal element pointer. */

/*<       jj = (n*(2*m - n + 1))/2 - (m - n) >*/
    jj = *n * ((*m << 1) - *n + 1) / 2 - (*m - *n);

/*     move the nontrivial part of the last column of s into w. */

/*<       l = jj >*/
    l = jj;
/*<       do 10 i = n, m >*/
    i__1 = *m;
    for (i__ = *n; i__ <= i__1; ++i__) {
/*<          w(i) = s(l) >*/
	w[i__] = s[l];
/*<          l = l + 1 >*/
	++l;
/*<    10    continue >*/
/* L10: */
    }

/*     rotate the vector v into a multiple of the n-th unit vector */
/*     in such a way that a spike is introduced into w. */

/*<       nm1 = n - 1 >*/
    nm1 = *n - 1;
/*<       if (nm1 .lt. 1) go to 70 >*/
    if (nm1 < 1) {
	goto L70;
    }
/*<       do 60 nmj = 1, nm1 >*/
    i__1 = nm1;
    for (nmj = 1; nmj <= i__1; ++nmj) {
/*<          j = n - nmj >*/
	j = *n - nmj;
/*<          jj = jj - (m - j + 1) >*/
	jj -= *m - j + 1;
/*<          w(j) = zero >*/
	w[j] = zero;
/*<          if (v(j) .eq. zero) go to 50 >*/
	if (v[j] == zero) {
	    goto L50;
	}

/*        determine a givens rotation which eliminates the */
/*        j-th element of v. */

/*<          if (dabs(v(n)) .ge. dabs(v(j))) go to 20 >*/
	if ((d__1 = v[*n], abs(d__1)) >= (d__2 = v[j], abs(d__2))) {
	    goto L20;
	}
/*<             cotan = v(n)/v(j) >*/
	cotan = v[*n] / v[j];
/*<             sin = p5/dsqrt(p25+p25*cotan**2) >*/
/* Computing 2nd power */
	d__1 = cotan;
	sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<             cos = sin*cotan >*/
	cos__ = sin__ * cotan;
/*<             tau = one >*/
	tau = one;
/*<             if (dabs(cos)*giant .gt. one) tau = one/cos >*/
	if (abs(cos__) * giant > one) {
	    tau = one / cos__;
	}
/*<             go to 30 >*/
	goto L30;
/*<    20    continue >*/
L20:
/*<             tan = v(j)/v(n) >*/
	tan__ = v[j] / v[*n];
/*<             cos = p5/dsqrt(p25+p25*tan**2) >*/
/* Computing 2nd power */
	d__1 = tan__;
	cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<             sin = cos*tan >*/
	sin__ = cos__ * tan__;
/*<             tau = sin >*/
	tau = sin__;
/*<    30    continue >*/
L30:

/*        apply the transformation to v and store the information */
/*        necessary to recover the givens rotation. */

/*<          v(n) = sin*v(j) + cos*v(n) >*/
	v[*n] = sin__ * v[j] + cos__ * v[*n];
/*<          v(j) = tau >*/
	v[j] = tau;

/*        apply the transformation to s and extend the spike in w. */

/*<          l = jj >*/
	l = jj;
/*<          do 40 i = j, m >*/
	i__2 = *m;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<             temp = cos*s(l) - sin*w(i) >*/
	    temp = cos__ * s[l] - sin__ * w[i__];
/*<             w(i) = sin*s(l) + cos*w(i) >*/
	    w[i__] = sin__ * s[l] + cos__ * w[i__];
/*<             s(l) = temp >*/
	    s[l] = temp;
/*<             l = l + 1 >*/
	    ++l;
/*<    40       continue >*/
/* L40: */
	}
/*<    50    continue >*/
L50:
/*<    60    continue >*/
/* L60: */
	;
    }
/*<    70 continue >*/
L70:

/*     add the spike from the rank 1 update to w. */

/*<       do 80 i = 1, m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          w(i) = w(i) + v(n)*u(i) >*/
	w[i__] += v[*n] * u[i__];
/*<    80    continue >*/
/* L80: */
    }

/*     eliminate the spike. */

/*<       sing = .false. >*/
    *sing = FALSE_;
/*<       if (nm1 .lt. 1) go to 140 >*/
    if (nm1 < 1) {
	goto L140;
    }
/*<       do 130 j = 1, nm1 >*/
    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
/*<          if (w(j) .eq. zero) go to 120 >*/
	if (w[j] == zero) {
	    goto L120;
	}

/*        determine a givens rotation which eliminates the */
/*        j-th element of the spike. */

/*<          if (dabs(s(jj)) .ge. dabs(w(j))) go to 90 >*/
	if ((d__1 = s[jj], abs(d__1)) >= (d__2 = w[j], abs(d__2))) {
	    goto L90;
	}
/*<             cotan = s(jj)/w(j) >*/
	cotan = s[jj] / w[j];
/*<             sin = p5/dsqrt(p25+p25*cotan**2) >*/
/* Computing 2nd power */
	d__1 = cotan;
	sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<             cos = sin*cotan >*/
	cos__ = sin__ * cotan;
/*<             tau = one >*/
	tau = one;
/*<             if (dabs(cos)*giant .gt. one) tau = one/cos >*/
	if (abs(cos__) * giant > one) {
	    tau = one / cos__;
	}
/*<             go to 100 >*/
	goto L100;
/*<    90    continue >*/
L90:
/*<             tan = w(j)/s(jj) >*/
	tan__ = w[j] / s[jj];
/*<             cos = p5/dsqrt(p25+p25*tan**2) >*/
/* Computing 2nd power */
	d__1 = tan__;
	cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<             sin = cos*tan >*/
	sin__ = cos__ * tan__;
/*<             tau = sin >*/
	tau = sin__;
/*<   100    continue >*/
L100:

/*        apply the transformation to s and reduce the spike in w. */

/*<          l = jj >*/
	l = jj;
/*<          do 110 i = j, m >*/
	i__2 = *m;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<             temp = cos*s(l) + sin*w(i) >*/
	    temp = cos__ * s[l] + sin__ * w[i__];
/*<             w(i) = -sin*s(l) + cos*w(i) >*/
	    w[i__] = -sin__ * s[l] + cos__ * w[i__];
/*<             s(l) = temp >*/
	    s[l] = temp;
/*<             l = l + 1 >*/
	    ++l;
/*<   110       continue >*/
/* L110: */
	}

/*        store the information necessary to recover the */
/*        givens rotation. */

/*<          w(j) = tau >*/
	w[j] = tau;
/*<   120    continue >*/
L120:

/*        test for zero diagonal elements in the output s. */

/*<          if (s(jj) .eq. zero) sing = .true. >*/
	if (s[jj] == zero) {
	    *sing = TRUE_;
	}
/*<          jj = jj + (m - j + 1) >*/
	jj += *m - j + 1;
/*<   130    continue >*/
/* L130: */
    }
/*<   140 continue >*/
L140:

/*     move w back into the last column of the output s. */

/*<       l = jj >*/
    l = jj;
/*<       do 150 i = n, m >*/
    i__1 = *m;
    for (i__ = *n; i__ <= i__1; ++i__) {
/*<          s(l) = w(i) >*/
	s[l] = w[i__];
/*<          l = l + 1 >*/
	++l;
/*<   150    continue >*/
/* L150: */
    }
/*<       if (s(jj) .eq. zero) sing = .true. >*/
    if (s[jj] == zero) {
	*sing = TRUE_;
    }
/*<       return >*/
    return 0;

/*     last card of subroutine r1updt. */

/*<       end >*/
} /* r1updt_ */

#ifdef __cplusplus
	}
#endif
