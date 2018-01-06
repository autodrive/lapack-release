#line 1 "dlagtf.f"
/* dlagtf.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "dlagtf.f"
/* > \brief \b DLAGTF computes an LU factorization of a matrix T-λI, where T is a general tridiagonal matrix,
 and λ a scalar, using partial pivoting with row interchanges. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAGTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagtf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagtf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagtf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   LAMBDA, TOL */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IN( * ) */
/*       DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n */
/* > tridiagonal matrix and lambda is a scalar, as */
/* > */
/* >    T - lambda*I = PLU, */
/* > */
/* > where P is a permutation matrix, L is a unit lower tridiagonal matrix */
/* > with at most one non-zero sub-diagonal elements per column and U is */
/* > an upper triangular matrix with at most two non-zero super-diagonal */
/* > elements per column. */
/* > */
/* > The factorization is obtained by Gaussian elimination with partial */
/* > pivoting and implicit row scaling. */
/* > */
/* > The parameter LAMBDA is included in the routine so that DLAGTF may */
/* > be used, in conjunction with DLAGTS, to obtain eigenvectors of T by */
/* > inverse iteration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, A must contain the diagonal elements of T. */
/* > */
/* >          On exit, A is overwritten by the n diagonal elements of the */
/* >          upper triangular matrix U of the factorization of T. */
/* > \endverbatim */
/* > */
/* > \param[in] LAMBDA */
/* > \verbatim */
/* >          LAMBDA is DOUBLE PRECISION */
/* >          On entry, the scalar lambda. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, B must contain the (n-1) super-diagonal elements of */
/* >          T. */
/* > */
/* >          On exit, B is overwritten by the (n-1) super-diagonal */
/* >          elements of the matrix U of the factorization of T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, C must contain the (n-1) sub-diagonal elements of */
/* >          T. */
/* > */
/* >          On exit, C is overwritten by the (n-1) sub-diagonal elements */
/* >          of the matrix L of the factorization of T. */
/* > \endverbatim */
/* > */
/* > \param[in] TOL */
/* > \verbatim */
/* >          TOL is DOUBLE PRECISION */
/* >          On entry, a relative tolerance used to indicate whether or */
/* >          not the matrix (T - lambda*I) is nearly singular. TOL should */
/* >          normally be chose as approximately the largest relative error */
/* >          in the elements of T. For example, if the elements of T are */
/* >          correct to about 4 significant figures, then TOL should be */
/* >          set to about 5*10**(-4). If TOL is supplied as less than eps, */
/* >          where eps is the relative machine precision, then the value */
/* >          eps is used in place of TOL. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N-2) */
/* >          On exit, D is overwritten by the (n-2) second super-diagonal */
/* >          elements of the matrix U of the factorization of T. */
/* > \endverbatim */
/* > */
/* > \param[out] IN */
/* > \verbatim */
/* >          IN is INTEGER array, dimension (N) */
/* >          On exit, IN contains details of the permutation matrix P. If */
/* >          an interchange occurred at the kth step of the elimination, */
/* >          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n) */
/* >          returns the smallest positive integer j such that */
/* > */
/* >             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL, */
/* > */
/* >          where norm( A(j) ) denotes the sum of the absolute values of */
/* >          the jth row of the matrix A. If no such j exists then IN(n) */
/* >          is returned as zero. If IN(n) is returned as positive, then a */
/* >          diagonal element of U is small, indicating that */
/* >          (T - lambda*I) is singular or nearly singular, */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0   : successful exit */
/* >          .lt. 0: if INFO = -k, the kth argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlagtf_(integer *n, doublereal *a, doublereal *lambda, 
	doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__, 
	integer *in, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer k;
    static doublereal tl, eps, piv1, piv2, temp, mult, scale1, scale2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 195 "dlagtf.f"
    /* Parameter adjustments */
#line 195 "dlagtf.f"
    --in;
#line 195 "dlagtf.f"
    --d__;
#line 195 "dlagtf.f"
    --c__;
#line 195 "dlagtf.f"
    --b;
#line 195 "dlagtf.f"
    --a;
#line 195 "dlagtf.f"

#line 195 "dlagtf.f"
    /* Function Body */
#line 195 "dlagtf.f"
    *info = 0;
#line 196 "dlagtf.f"
    if (*n < 0) {
#line 197 "dlagtf.f"
	*info = -1;
#line 198 "dlagtf.f"
	i__1 = -(*info);
#line 198 "dlagtf.f"
	xerbla_("DLAGTF", &i__1, (ftnlen)6);
#line 199 "dlagtf.f"
	return 0;
#line 200 "dlagtf.f"
    }

#line 202 "dlagtf.f"
    if (*n == 0) {
#line 202 "dlagtf.f"
	return 0;
#line 202 "dlagtf.f"
    }

#line 205 "dlagtf.f"
    a[1] -= *lambda;
#line 206 "dlagtf.f"
    in[*n] = 0;
#line 207 "dlagtf.f"
    if (*n == 1) {
#line 208 "dlagtf.f"
	if (a[1] == 0.) {
#line 208 "dlagtf.f"
	    in[1] = 1;
#line 208 "dlagtf.f"
	}
#line 210 "dlagtf.f"
	return 0;
#line 211 "dlagtf.f"
    }

#line 213 "dlagtf.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

#line 215 "dlagtf.f"
    tl = max(*tol,eps);
#line 216 "dlagtf.f"
    scale1 = abs(a[1]) + abs(b[1]);
#line 217 "dlagtf.f"
    i__1 = *n - 1;
#line 217 "dlagtf.f"
    for (k = 1; k <= i__1; ++k) {
#line 218 "dlagtf.f"
	a[k + 1] -= *lambda;
#line 219 "dlagtf.f"
	scale2 = (d__1 = c__[k], abs(d__1)) + (d__2 = a[k + 1], abs(d__2));
#line 220 "dlagtf.f"
	if (k < *n - 1) {
#line 220 "dlagtf.f"
	    scale2 += (d__1 = b[k + 1], abs(d__1));
#line 220 "dlagtf.f"
	}
#line 222 "dlagtf.f"
	if (a[k] == 0.) {
#line 223 "dlagtf.f"
	    piv1 = 0.;
#line 224 "dlagtf.f"
	} else {
#line 225 "dlagtf.f"
	    piv1 = (d__1 = a[k], abs(d__1)) / scale1;
#line 226 "dlagtf.f"
	}
#line 227 "dlagtf.f"
	if (c__[k] == 0.) {
#line 228 "dlagtf.f"
	    in[k] = 0;
#line 229 "dlagtf.f"
	    piv2 = 0.;
#line 230 "dlagtf.f"
	    scale1 = scale2;
#line 231 "dlagtf.f"
	    if (k < *n - 1) {
#line 231 "dlagtf.f"
		d__[k] = 0.;
#line 231 "dlagtf.f"
	    }
#line 233 "dlagtf.f"
	} else {
#line 234 "dlagtf.f"
	    piv2 = (d__1 = c__[k], abs(d__1)) / scale2;
#line 235 "dlagtf.f"
	    if (piv2 <= piv1) {
#line 236 "dlagtf.f"
		in[k] = 0;
#line 237 "dlagtf.f"
		scale1 = scale2;
#line 238 "dlagtf.f"
		c__[k] /= a[k];
#line 239 "dlagtf.f"
		a[k + 1] -= c__[k] * b[k];
#line 240 "dlagtf.f"
		if (k < *n - 1) {
#line 240 "dlagtf.f"
		    d__[k] = 0.;
#line 240 "dlagtf.f"
		}
#line 242 "dlagtf.f"
	    } else {
#line 243 "dlagtf.f"
		in[k] = 1;
#line 244 "dlagtf.f"
		mult = a[k] / c__[k];
#line 245 "dlagtf.f"
		a[k] = c__[k];
#line 246 "dlagtf.f"
		temp = a[k + 1];
#line 247 "dlagtf.f"
		a[k + 1] = b[k] - mult * temp;
#line 248 "dlagtf.f"
		if (k < *n - 1) {
#line 249 "dlagtf.f"
		    d__[k] = b[k + 1];
#line 250 "dlagtf.f"
		    b[k + 1] = -mult * d__[k];
#line 251 "dlagtf.f"
		}
#line 252 "dlagtf.f"
		b[k] = temp;
#line 253 "dlagtf.f"
		c__[k] = mult;
#line 254 "dlagtf.f"
	    }
#line 255 "dlagtf.f"
	}
#line 256 "dlagtf.f"
	if (max(piv1,piv2) <= tl && in[*n] == 0) {
#line 256 "dlagtf.f"
	    in[*n] = k;
#line 256 "dlagtf.f"
	}
#line 258 "dlagtf.f"
/* L10: */
#line 258 "dlagtf.f"
    }
#line 259 "dlagtf.f"
    if ((d__1 = a[*n], abs(d__1)) <= scale1 * tl && in[*n] == 0) {
#line 259 "dlagtf.f"
	in[*n] = *n;
#line 259 "dlagtf.f"
    }

#line 262 "dlagtf.f"
    return 0;

/*     End of DLAGTF */

} /* dlagtf_ */

