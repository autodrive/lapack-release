#line 1 "dgbequ.f"
/* dgbequ.f -- translated by f2c (version 20100827).
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

#line 1 "dgbequ.f"
/* > \brief \b DGBEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                          AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, M, N */
/*       DOUBLE PRECISION   AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBEQU computes row and column scalings intended to equilibrate an */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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

/* > \date December 2016 */

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgbequ_(integer *m, integer *n, integer *kl, integer *ku,
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
    extern doublereal dlamch_(char *, ftnlen);
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
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 193 "dgbequ.f"
    /* Parameter adjustments */
#line 193 "dgbequ.f"
    ab_dim1 = *ldab;
#line 193 "dgbequ.f"
    ab_offset = 1 + ab_dim1;
#line 193 "dgbequ.f"
    ab -= ab_offset;
#line 193 "dgbequ.f"
    --r__;
#line 193 "dgbequ.f"
    --c__;
#line 193 "dgbequ.f"

#line 193 "dgbequ.f"
    /* Function Body */
#line 193 "dgbequ.f"
    *info = 0;
#line 194 "dgbequ.f"
    if (*m < 0) {
#line 195 "dgbequ.f"
	*info = -1;
#line 196 "dgbequ.f"
    } else if (*n < 0) {
#line 197 "dgbequ.f"
	*info = -2;
#line 198 "dgbequ.f"
    } else if (*kl < 0) {
#line 199 "dgbequ.f"
	*info = -3;
#line 200 "dgbequ.f"
    } else if (*ku < 0) {
#line 201 "dgbequ.f"
	*info = -4;
#line 202 "dgbequ.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 203 "dgbequ.f"
	*info = -6;
#line 204 "dgbequ.f"
    }
#line 205 "dgbequ.f"
    if (*info != 0) {
#line 206 "dgbequ.f"
	i__1 = -(*info);
#line 206 "dgbequ.f"
	xerbla_("DGBEQU", &i__1, (ftnlen)6);
#line 207 "dgbequ.f"
	return 0;
#line 208 "dgbequ.f"
    }

/*     Quick return if possible */

#line 212 "dgbequ.f"
    if (*m == 0 || *n == 0) {
#line 213 "dgbequ.f"
	*rowcnd = 1.;
#line 214 "dgbequ.f"
	*colcnd = 1.;
#line 215 "dgbequ.f"
	*amax = 0.;
#line 216 "dgbequ.f"
	return 0;
#line 217 "dgbequ.f"
    }

/*     Get machine constants. */

#line 221 "dgbequ.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 222 "dgbequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 226 "dgbequ.f"
    i__1 = *m;
#line 226 "dgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 227 "dgbequ.f"
	r__[i__] = 0.;
#line 228 "dgbequ.f"
/* L10: */
#line 228 "dgbequ.f"
    }

/*     Find the maximum element in each row. */

#line 232 "dgbequ.f"
    kd = *ku + 1;
#line 233 "dgbequ.f"
    i__1 = *n;
#line 233 "dgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 234 "dgbequ.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 234 "dgbequ.f"
	i__4 = j + *kl;
#line 234 "dgbequ.f"
	i__3 = min(i__4,*m);
#line 234 "dgbequ.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 235 "dgbequ.f"
	    d__2 = r__[i__], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], 
		    abs(d__1));
#line 235 "dgbequ.f"
	    r__[i__] = max(d__2,d__3);
#line 236 "dgbequ.f"
/* L20: */
#line 236 "dgbequ.f"
	}
#line 237 "dgbequ.f"
/* L30: */
#line 237 "dgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 241 "dgbequ.f"
    rcmin = bignum;
#line 242 "dgbequ.f"
    rcmax = 0.;
#line 243 "dgbequ.f"
    i__1 = *m;
#line 243 "dgbequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 244 "dgbequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 244 "dgbequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 245 "dgbequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 245 "dgbequ.f"
	rcmin = min(d__1,d__2);
#line 246 "dgbequ.f"
/* L40: */
#line 246 "dgbequ.f"
    }
#line 247 "dgbequ.f"
    *amax = rcmax;

#line 249 "dgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 253 "dgbequ.f"
	i__1 = *m;
#line 253 "dgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 254 "dgbequ.f"
	    if (r__[i__] == 0.) {
#line 255 "dgbequ.f"
		*info = i__;
#line 256 "dgbequ.f"
		return 0;
#line 257 "dgbequ.f"
	    }
#line 258 "dgbequ.f"
/* L50: */
#line 258 "dgbequ.f"
	}
#line 259 "dgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 263 "dgbequ.f"
	i__1 = *m;
#line 263 "dgbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 264 "dgbequ.f"
	    d__2 = r__[i__];
#line 264 "dgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 264 "dgbequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 265 "dgbequ.f"
/* L60: */
#line 265 "dgbequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 269 "dgbequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 270 "dgbequ.f"
    }

/*     Compute column scale factors */

#line 274 "dgbequ.f"
    i__1 = *n;
#line 274 "dgbequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 275 "dgbequ.f"
	c__[j] = 0.;
#line 276 "dgbequ.f"
/* L70: */
#line 276 "dgbequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 281 "dgbequ.f"
    kd = *ku + 1;
#line 282 "dgbequ.f"
    i__1 = *n;
#line 282 "dgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 283 "dgbequ.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 283 "dgbequ.f"
	i__4 = j + *kl;
#line 283 "dgbequ.f"
	i__2 = min(i__4,*m);
#line 283 "dgbequ.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 284 "dgbequ.f"
	    d__2 = c__[j], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
		    d__1)) * r__[i__];
#line 284 "dgbequ.f"
	    c__[j] = max(d__2,d__3);
#line 285 "dgbequ.f"
/* L80: */
#line 285 "dgbequ.f"
	}
#line 286 "dgbequ.f"
/* L90: */
#line 286 "dgbequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 290 "dgbequ.f"
    rcmin = bignum;
#line 291 "dgbequ.f"
    rcmax = 0.;
#line 292 "dgbequ.f"
    i__1 = *n;
#line 292 "dgbequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 293 "dgbequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 293 "dgbequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 294 "dgbequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 294 "dgbequ.f"
	rcmax = max(d__1,d__2);
#line 295 "dgbequ.f"
/* L100: */
#line 295 "dgbequ.f"
    }

#line 297 "dgbequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 301 "dgbequ.f"
	i__1 = *n;
#line 301 "dgbequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 302 "dgbequ.f"
	    if (c__[j] == 0.) {
#line 303 "dgbequ.f"
		*info = *m + j;
#line 304 "dgbequ.f"
		return 0;
#line 305 "dgbequ.f"
	    }
#line 306 "dgbequ.f"
/* L110: */
#line 306 "dgbequ.f"
	}
#line 307 "dgbequ.f"
    } else {

/*        Invert the scale factors. */

#line 311 "dgbequ.f"
	i__1 = *n;
#line 311 "dgbequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 312 "dgbequ.f"
	    d__2 = c__[j];
#line 312 "dgbequ.f"
	    d__1 = max(d__2,smlnum);
#line 312 "dgbequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 313 "dgbequ.f"
/* L120: */
#line 313 "dgbequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 317 "dgbequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 318 "dgbequ.f"
    }

#line 320 "dgbequ.f"
    return 0;

/*     End of DGBEQU */

} /* dgbequ_ */

