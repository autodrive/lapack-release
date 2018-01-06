#line 1 "cgbequ.f"
/* cgbequ.f -- translated by f2c (version 20100827).
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

#line 1 "cgbequ.f"
/* > \brief \b CGBEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                          AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ), R( * ) */
/*       COMPLEX            AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBEQU computes row and column scalings intended to equilibrate an */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          R is REAL array, dimension (M) */
/* >          If INFO = 0, or INFO > M, R contains the row scale factors */
/* >          for A. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >          If INFO = 0, C contains the column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] ROWCND */
/* > \verbatim */
/* >          ROWCND is REAL */
/* >          If INFO = 0 or INFO > M, ROWCND contains the ratio of the */
/* >          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and */
/* >          AMAX is neither too large nor too small, it is not worth */
/* >          scaling by R. */
/* > \endverbatim */
/* > */
/* > \param[out] COLCND */
/* > \verbatim */
/* >          COLCND is REAL */
/* >          If INFO = 0, COLCND contains the ratio of the smallest */
/* >          C(i) to the largest C(i).  If COLCND >= 0.1, it is not */
/* >          worth scaling by C. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is REAL */
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

/* > \date December 2016 */

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgbequ_(integer *m, integer *n, integer *kl, integer *ku,
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
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 202 "cgbequ.f"
    /* Parameter adjustments */
#line 202 "cgbequ.f"
    ab_dim1 = *ldab;
#line 202 "cgbequ.f"
    ab_offset = 1 + ab_dim1;
#line 202 "cgbequ.f"
    ab -= ab_offset;
#line 202 "cgbequ.f"
    --r__;
#line 202 "cgbequ.f"
    --c__;
#line 202 "cgbequ.f"

#line 202 "cgbequ.f"
    /* Function Body */
#line 202 "cgbequ.f"
    *info = 0;
#line 203 "cgbequ.f"
    if (*m < 0) {
#line 204 "cgbequ.f"
	*info = -1;
#line 205 "cgbequ.f"
    } else if (*n < 0) {
#line 206 "cgbequ.f"
	*info = -2;
#line 207 "cgbequ.f"
    } else if (*kl < 0) {
#line 208 "cgbequ.f"
	*info = -3;
#line 209 "cgbequ.f"
    } else if (*ku < 0) {
#line 210 "cgbequ.f"
	*info = -4;
#line 211 "cgbequ.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 212 "cgbequ.f"
	*info = -6;
#line 213 "cgbequ.f"
    }
#line 214 "cgbequ.f"
    if (*info != 0) {
#line 215 "cgbequ.f"
	i__1 = -(*info);
#line 215 "cgbequ.f"
	xerbla_("CGBEQU", &i__1, (ftnlen)6);
#line 216 "cgbequ.f"
	return 0;
#line 217 "cgbequ.f"
    }

/*     Quick return if possible */

#line 221 "cgbequ.f"
    if (*m == 0 || *n == 0) {
#line 222 "cgbequ.f"
	*rowcnd = 1.;
#line 223 "cgbequ.f"
	*colcnd = 1.;
#line 224 "cgbequ.f"
	*amax = 0.;
#line 225 "cgbequ.f"
	return 0;
#line 226 "cgbequ.f"
    }

/*     Get machine constants. */

#line 230 "cgbequ.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 231 "cgbequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 235 "cgbequ.f"
    i__1 = *m;
#line 235 "cgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "cgbequ.f"
	r__[i__] = 0.;
#line 237 "cgbequ.f"
/* L10: */
#line 237 "cgbequ.f"
    }

/*     Find the maximum element in each row. */

#line 241 "cgbequ.f"
    kd = *ku + 1;
#line 242 "cgbequ.f"
    i__1 = *n;
#line 242 "cgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 243 "cgbequ.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 243 "cgbequ.f"
	i__4 = j + *kl;
#line 243 "cgbequ.f"
	i__3 = min(i__4,*m);
#line 243 "cgbequ.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 244 "cgbequ.f"
	    i__2 = kd + i__ - j + j * ab_dim1;
#line 244 "cgbequ.f"
	    d__3 = r__[i__], d__4 = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2));
#line 244 "cgbequ.f"
	    r__[i__] = max(d__3,d__4);
#line 245 "cgbequ.f"
/* L20: */
#line 245 "cgbequ.f"
	}
#line 246 "cgbequ.f"
/* L30: */
#line 246 "cgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 250 "cgbequ.f"
    rcmin = bignum;
#line 251 "cgbequ.f"
    rcmax = 0.;
#line 252 "cgbequ.f"
    i__1 = *m;
#line 252 "cgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 253 "cgbequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 253 "cgbequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 254 "cgbequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 254 "cgbequ.f"
	rcmin = min(d__1,d__2);
#line 255 "cgbequ.f"
/* L40: */
#line 255 "cgbequ.f"
    }
#line 256 "cgbequ.f"
    *amax = rcmax;

#line 258 "cgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 262 "cgbequ.f"
	i__1 = *m;
#line 262 "cgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 263 "cgbequ.f"
	    if (r__[i__] == 0.) {
#line 264 "cgbequ.f"
		*info = i__;
#line 265 "cgbequ.f"
		return 0;
#line 266 "cgbequ.f"
	    }
#line 267 "cgbequ.f"
/* L50: */
#line 267 "cgbequ.f"
	}
#line 268 "cgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 272 "cgbequ.f"
	i__1 = *m;
#line 272 "cgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 273 "cgbequ.f"
	    d__2 = r__[i__];
#line 273 "cgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 273 "cgbequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 274 "cgbequ.f"
/* L60: */
#line 274 "cgbequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 278 "cgbequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 279 "cgbequ.f"
    }

/*     Compute column scale factors */

#line 283 "cgbequ.f"
    i__1 = *n;
#line 283 "cgbequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 284 "cgbequ.f"
	c__[j] = 0.;
#line 285 "cgbequ.f"
/* L70: */
#line 285 "cgbequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 290 "cgbequ.f"
    kd = *ku + 1;
#line 291 "cgbequ.f"
    i__1 = *n;
#line 291 "cgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 292 "cgbequ.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 292 "cgbequ.f"
	i__4 = j + *kl;
#line 292 "cgbequ.f"
	i__2 = min(i__4,*m);
#line 292 "cgbequ.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 293 "cgbequ.f"
	    i__3 = kd + i__ - j + j * ab_dim1;
#line 293 "cgbequ.f"
	    d__3 = c__[j], d__4 = ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2))) * 
		    r__[i__];
#line 293 "cgbequ.f"
	    c__[j] = max(d__3,d__4);
#line 294 "cgbequ.f"
/* L80: */
#line 294 "cgbequ.f"
	}
#line 295 "cgbequ.f"
/* L90: */
#line 295 "cgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 299 "cgbequ.f"
    rcmin = bignum;
#line 300 "cgbequ.f"
    rcmax = 0.;
#line 301 "cgbequ.f"
    i__1 = *n;
#line 301 "cgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 302 "cgbequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 302 "cgbequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 303 "cgbequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 303 "cgbequ.f"
	rcmax = max(d__1,d__2);
#line 304 "cgbequ.f"
/* L100: */
#line 304 "cgbequ.f"
    }

#line 306 "cgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 310 "cgbequ.f"
	i__1 = *n;
#line 310 "cgbequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 311 "cgbequ.f"
	    if (c__[j] == 0.) {
#line 312 "cgbequ.f"
		*info = *m + j;
#line 313 "cgbequ.f"
		return 0;
#line 314 "cgbequ.f"
	    }
#line 315 "cgbequ.f"
/* L110: */
#line 315 "cgbequ.f"
	}
#line 316 "cgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 320 "cgbequ.f"
	i__1 = *n;
#line 320 "cgbequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 321 "cgbequ.f"
	    d__2 = c__[j];
#line 321 "cgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 321 "cgbequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 322 "cgbequ.f"
/* L120: */
#line 322 "cgbequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 326 "cgbequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 327 "cgbequ.f"
    }

#line 329 "cgbequ.f"
    return 0;

/*     End of CGBEQU */

} /* cgbequ_ */

