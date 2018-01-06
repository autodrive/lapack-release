#line 1 "zgbequ.f"
/* zgbequ.f -- translated by f2c (version 20100827).
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

#line 1 "zgbequ.f"
/* > \brief \b ZGBEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                          AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       DOUBLE PRECISION   AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), R( * ) */
/*       COMPLEX*16         AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBEQU computes row and column scalings intended to equilibrate an */
/* > M-by-N band matrix A and reduce its condition number.  R returns the */
/* > row scale factors and C the column scale factors, chosen to try to */
/* > make the largest element in each row and column of the matrix B with */
/* > elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1. */
/* > */
/* > R(i) and C(j) are restricted to be between SMLNUM = smallest safe */
/* > number and BIGNUM = largest safe number.  Use of these scaling */
/* > factors is not guaranteed to reduce the condition number of A but */
/* > works well in practice. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th */
/* >          column of A is stored in the j-th column of the array AB as */
/* >          follows: */
/* >          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is DOUBLE PRECISION array, dimension (M) */
/* >          If INFO = 0, or INFO > M, R contains the row scale factors */
/* >          for A. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, C contains the column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] ROWCND */
/* > \verbatim */
/* >          ROWCND is DOUBLE PRECISION */
/* >          If INFO = 0 or INFO > M, ROWCND contains the ratio of the */
/* >          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and */
/* >          AMAX is neither too large nor too small, it is not worth */
/* >          scaling by R. */
/* > \endverbatim */
/* > */
/* > \param[out] COLCND */
/* > \verbatim */
/* >          COLCND is DOUBLE PRECISION */
/* >          If INFO = 0, COLCND contains the ratio of the smallest */
/* >          C(i) to the largest C(i).  If COLCND >= 0.1, it is not */
/* >          worth scaling by C. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is DOUBLE PRECISION */
/* >          Absolute value of largest matrix element.  If AMAX is very */
/* >          close to overflow or very close to underflow, the matrix */
/* >          should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, and i is */
/* >                <= M:  the i-th row of A is exactly zero */
/* >                >  M:  the (i-M)-th column of A is exactly zero */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16GBcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, kd;
    static doublereal rcmin, rcmax;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 202 "zgbequ.f"
    /* Parameter adjustments */
#line 202 "zgbequ.f"
    ab_dim1 = *ldab;
#line 202 "zgbequ.f"
    ab_offset = 1 + ab_dim1;
#line 202 "zgbequ.f"
    ab -= ab_offset;
#line 202 "zgbequ.f"
    --r__;
#line 202 "zgbequ.f"
    --c__;
#line 202 "zgbequ.f"

#line 202 "zgbequ.f"
    /* Function Body */
#line 202 "zgbequ.f"
    *info = 0;
#line 203 "zgbequ.f"
    if (*m < 0) {
#line 204 "zgbequ.f"
	*info = -1;
#line 205 "zgbequ.f"
    } else if (*n < 0) {
#line 206 "zgbequ.f"
	*info = -2;
#line 207 "zgbequ.f"
    } else if (*kl < 0) {
#line 208 "zgbequ.f"
	*info = -3;
#line 209 "zgbequ.f"
    } else if (*ku < 0) {
#line 210 "zgbequ.f"
	*info = -4;
#line 211 "zgbequ.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 212 "zgbequ.f"
	*info = -6;
#line 213 "zgbequ.f"
    }
#line 214 "zgbequ.f"
    if (*info != 0) {
#line 215 "zgbequ.f"
	i__1 = -(*info);
#line 215 "zgbequ.f"
	xerbla_("ZGBEQU", &i__1, (ftnlen)6);
#line 216 "zgbequ.f"
	return 0;
#line 217 "zgbequ.f"
    }

/*     Quick return if possible */

#line 221 "zgbequ.f"
    if (*m == 0 || *n == 0) {
#line 222 "zgbequ.f"
	*rowcnd = 1.;
#line 223 "zgbequ.f"
	*colcnd = 1.;
#line 224 "zgbequ.f"
	*amax = 0.;
#line 225 "zgbequ.f"
	return 0;
#line 226 "zgbequ.f"
    }

/*     Get machine constants. */

#line 230 "zgbequ.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 231 "zgbequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 235 "zgbequ.f"
    i__1 = *m;
#line 235 "zgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "zgbequ.f"
	r__[i__] = 0.;
#line 237 "zgbequ.f"
/* L10: */
#line 237 "zgbequ.f"
    }

/*     Find the maximum element in each row. */

#line 241 "zgbequ.f"
    kd = *ku + 1;
#line 242 "zgbequ.f"
    i__1 = *n;
#line 242 "zgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 243 "zgbequ.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 243 "zgbequ.f"
	i__4 = j + *kl;
#line 243 "zgbequ.f"
	i__3 = min(i__4,*m);
#line 243 "zgbequ.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 244 "zgbequ.f"
	    i__2 = kd + i__ - j + j * ab_dim1;
#line 244 "zgbequ.f"
	    d__3 = r__[i__], d__4 = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2));
#line 244 "zgbequ.f"
	    r__[i__] = max(d__3,d__4);
#line 245 "zgbequ.f"
/* L20: */
#line 245 "zgbequ.f"
	}
#line 246 "zgbequ.f"
/* L30: */
#line 246 "zgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 250 "zgbequ.f"
    rcmin = bignum;
#line 251 "zgbequ.f"
    rcmax = 0.;
#line 252 "zgbequ.f"
    i__1 = *m;
#line 252 "zgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 253 "zgbequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 253 "zgbequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 254 "zgbequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 254 "zgbequ.f"
	rcmin = min(d__1,d__2);
#line 255 "zgbequ.f"
/* L40: */
#line 255 "zgbequ.f"
    }
#line 256 "zgbequ.f"
    *amax = rcmax;

#line 258 "zgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 262 "zgbequ.f"
	i__1 = *m;
#line 262 "zgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 263 "zgbequ.f"
	    if (r__[i__] == 0.) {
#line 264 "zgbequ.f"
		*info = i__;
#line 265 "zgbequ.f"
		return 0;
#line 266 "zgbequ.f"
	    }
#line 267 "zgbequ.f"
/* L50: */
#line 267 "zgbequ.f"
	}
#line 268 "zgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 272 "zgbequ.f"
	i__1 = *m;
#line 272 "zgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 273 "zgbequ.f"
	    d__2 = r__[i__];
#line 273 "zgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 273 "zgbequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 274 "zgbequ.f"
/* L60: */
#line 274 "zgbequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 278 "zgbequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 279 "zgbequ.f"
    }

/*     Compute column scale factors */

#line 283 "zgbequ.f"
    i__1 = *n;
#line 283 "zgbequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 284 "zgbequ.f"
	c__[j] = 0.;
#line 285 "zgbequ.f"
/* L70: */
#line 285 "zgbequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 290 "zgbequ.f"
    kd = *ku + 1;
#line 291 "zgbequ.f"
    i__1 = *n;
#line 291 "zgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 292 "zgbequ.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 292 "zgbequ.f"
	i__4 = j + *kl;
#line 292 "zgbequ.f"
	i__2 = min(i__4,*m);
#line 292 "zgbequ.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 293 "zgbequ.f"
	    i__3 = kd + i__ - j + j * ab_dim1;
#line 293 "zgbequ.f"
	    d__3 = c__[j], d__4 = ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2))) * 
		    r__[i__];
#line 293 "zgbequ.f"
	    c__[j] = max(d__3,d__4);
#line 294 "zgbequ.f"
/* L80: */
#line 294 "zgbequ.f"
	}
#line 295 "zgbequ.f"
/* L90: */
#line 295 "zgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 299 "zgbequ.f"
    rcmin = bignum;
#line 300 "zgbequ.f"
    rcmax = 0.;
#line 301 "zgbequ.f"
    i__1 = *n;
#line 301 "zgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 302 "zgbequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 302 "zgbequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 303 "zgbequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 303 "zgbequ.f"
	rcmax = max(d__1,d__2);
#line 304 "zgbequ.f"
/* L100: */
#line 304 "zgbequ.f"
    }

#line 306 "zgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 310 "zgbequ.f"
	i__1 = *n;
#line 310 "zgbequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 311 "zgbequ.f"
	    if (c__[j] == 0.) {
#line 312 "zgbequ.f"
		*info = *m + j;
#line 313 "zgbequ.f"
		return 0;
#line 314 "zgbequ.f"
	    }
#line 315 "zgbequ.f"
/* L110: */
#line 315 "zgbequ.f"
	}
#line 316 "zgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 320 "zgbequ.f"
	i__1 = *n;
#line 320 "zgbequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 321 "zgbequ.f"
	    d__2 = c__[j];
#line 321 "zgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 321 "zgbequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 322 "zgbequ.f"
/* L120: */
#line 322 "zgbequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 326 "zgbequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 327 "zgbequ.f"
    }

#line 329 "zgbequ.f"
    return 0;

/*     End of ZGBEQU */

} /* zgbequ_ */

