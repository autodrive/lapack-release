#line 1 "sgbequ.f"
/* sgbequ.f -- translated by f2c (version 20100827).
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

#line 1 "sgbequ.f"
/* > \brief \b SGBEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                          AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBEQU computes row and column scalings intended to equilibrate an */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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

/* > \date November 2011 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, kd;
    static doublereal rcmin, rcmax;
    extern doublereal slamch_(char *, ftnlen);
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
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 193 "sgbequ.f"
    /* Parameter adjustments */
#line 193 "sgbequ.f"
    ab_dim1 = *ldab;
#line 193 "sgbequ.f"
    ab_offset = 1 + ab_dim1;
#line 193 "sgbequ.f"
    ab -= ab_offset;
#line 193 "sgbequ.f"
    --r__;
#line 193 "sgbequ.f"
    --c__;
#line 193 "sgbequ.f"

#line 193 "sgbequ.f"
    /* Function Body */
#line 193 "sgbequ.f"
    *info = 0;
#line 194 "sgbequ.f"
    if (*m < 0) {
#line 195 "sgbequ.f"
	*info = -1;
#line 196 "sgbequ.f"
    } else if (*n < 0) {
#line 197 "sgbequ.f"
	*info = -2;
#line 198 "sgbequ.f"
    } else if (*kl < 0) {
#line 199 "sgbequ.f"
	*info = -3;
#line 200 "sgbequ.f"
    } else if (*ku < 0) {
#line 201 "sgbequ.f"
	*info = -4;
#line 202 "sgbequ.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 203 "sgbequ.f"
	*info = -6;
#line 204 "sgbequ.f"
    }
#line 205 "sgbequ.f"
    if (*info != 0) {
#line 206 "sgbequ.f"
	i__1 = -(*info);
#line 206 "sgbequ.f"
	xerbla_("SGBEQU", &i__1, (ftnlen)6);
#line 207 "sgbequ.f"
	return 0;
#line 208 "sgbequ.f"
    }

/*     Quick return if possible */

#line 212 "sgbequ.f"
    if (*m == 0 || *n == 0) {
#line 213 "sgbequ.f"
	*rowcnd = 1.;
#line 214 "sgbequ.f"
	*colcnd = 1.;
#line 215 "sgbequ.f"
	*amax = 0.;
#line 216 "sgbequ.f"
	return 0;
#line 217 "sgbequ.f"
    }

/*     Get machine constants. */

#line 221 "sgbequ.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 222 "sgbequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 226 "sgbequ.f"
    i__1 = *m;
#line 226 "sgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 227 "sgbequ.f"
	r__[i__] = 0.;
#line 228 "sgbequ.f"
/* L10: */
#line 228 "sgbequ.f"
    }

/*     Find the maximum element in each row. */

#line 232 "sgbequ.f"
    kd = *ku + 1;
#line 233 "sgbequ.f"
    i__1 = *n;
#line 233 "sgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 234 "sgbequ.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 234 "sgbequ.f"
	i__4 = j + *kl;
#line 234 "sgbequ.f"
	i__3 = min(i__4,*m);
#line 234 "sgbequ.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 235 "sgbequ.f"
	    d__2 = r__[i__], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], 
		    abs(d__1));
#line 235 "sgbequ.f"
	    r__[i__] = max(d__2,d__3);
#line 236 "sgbequ.f"
/* L20: */
#line 236 "sgbequ.f"
	}
#line 237 "sgbequ.f"
/* L30: */
#line 237 "sgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 241 "sgbequ.f"
    rcmin = bignum;
#line 242 "sgbequ.f"
    rcmax = 0.;
#line 243 "sgbequ.f"
    i__1 = *m;
#line 243 "sgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 244 "sgbequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 244 "sgbequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 245 "sgbequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 245 "sgbequ.f"
	rcmin = min(d__1,d__2);
#line 246 "sgbequ.f"
/* L40: */
#line 246 "sgbequ.f"
    }
#line 247 "sgbequ.f"
    *amax = rcmax;

#line 249 "sgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 253 "sgbequ.f"
	i__1 = *m;
#line 253 "sgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 254 "sgbequ.f"
	    if (r__[i__] == 0.) {
#line 255 "sgbequ.f"
		*info = i__;
#line 256 "sgbequ.f"
		return 0;
#line 257 "sgbequ.f"
	    }
#line 258 "sgbequ.f"
/* L50: */
#line 258 "sgbequ.f"
	}
#line 259 "sgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 263 "sgbequ.f"
	i__1 = *m;
#line 263 "sgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 264 "sgbequ.f"
	    d__2 = r__[i__];
#line 264 "sgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 264 "sgbequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 265 "sgbequ.f"
/* L60: */
#line 265 "sgbequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 269 "sgbequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 270 "sgbequ.f"
    }

/*     Compute column scale factors */

#line 274 "sgbequ.f"
    i__1 = *n;
#line 274 "sgbequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 275 "sgbequ.f"
	c__[j] = 0.;
#line 276 "sgbequ.f"
/* L70: */
#line 276 "sgbequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 281 "sgbequ.f"
    kd = *ku + 1;
#line 282 "sgbequ.f"
    i__1 = *n;
#line 282 "sgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 283 "sgbequ.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 283 "sgbequ.f"
	i__4 = j + *kl;
#line 283 "sgbequ.f"
	i__2 = min(i__4,*m);
#line 283 "sgbequ.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 284 "sgbequ.f"
	    d__2 = c__[j], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
		    d__1)) * r__[i__];
#line 284 "sgbequ.f"
	    c__[j] = max(d__2,d__3);
#line 285 "sgbequ.f"
/* L80: */
#line 285 "sgbequ.f"
	}
#line 286 "sgbequ.f"
/* L90: */
#line 286 "sgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 290 "sgbequ.f"
    rcmin = bignum;
#line 291 "sgbequ.f"
    rcmax = 0.;
#line 292 "sgbequ.f"
    i__1 = *n;
#line 292 "sgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 293 "sgbequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 293 "sgbequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 294 "sgbequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 294 "sgbequ.f"
	rcmax = max(d__1,d__2);
#line 295 "sgbequ.f"
/* L100: */
#line 295 "sgbequ.f"
    }

#line 297 "sgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 301 "sgbequ.f"
	i__1 = *n;
#line 301 "sgbequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 302 "sgbequ.f"
	    if (c__[j] == 0.) {
#line 303 "sgbequ.f"
		*info = *m + j;
#line 304 "sgbequ.f"
		return 0;
#line 305 "sgbequ.f"
	    }
#line 306 "sgbequ.f"
/* L110: */
#line 306 "sgbequ.f"
	}
#line 307 "sgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 311 "sgbequ.f"
	i__1 = *n;
#line 311 "sgbequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 312 "sgbequ.f"
	    d__2 = c__[j];
#line 312 "sgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 312 "sgbequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 313 "sgbequ.f"
/* L120: */
#line 313 "sgbequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 317 "sgbequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 318 "sgbequ.f"
    }

#line 320 "sgbequ.f"
    return 0;

/*     End of SGBEQU */

} /* sgbequ_ */

