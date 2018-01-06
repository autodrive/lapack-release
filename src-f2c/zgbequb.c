#line 1 "zgbequb.f"
/* zgbequb.f -- translated by f2c (version 20100827).
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

#line 1 "zgbequb.f"
/* > \brief \b ZGBEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                           AMAX, INFO ) */

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
/* > ZGBEQUB computes row and column scalings intended to equilibrate an */
/* > M-by-N matrix A and reduce its condition number.  R returns the row */
/* > scale factors and C the column scale factors, chosen to try to make */
/* > the largest element in each row and column of the matrix B with */
/* > elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most */
/* > the radix. */
/* > */
/* > R(i) and C(j) are restricted to be a power of the radix between */
/* > SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use */
/* > of these scaling factors is not guaranteed to reduce the condition */
/* > number of A but works well in practice. */
/* > */
/* > This routine differs from ZGEEQU by restricting the scaling factors */
/* > to a power of the radix.  Baring over- and underflow, scaling by */
/* > these factors introduces no additional rounding errors.  However, the */
/* > scaled entries' magnitured are no longer approximately 1 but lie */
/* > between sqrt(radix) and 1/sqrt(radix). */
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
/* >          On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl) */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array A.  LDAB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is DOUBLE PRECISION array, dimension (M) */
/* >          If INFO = 0 or INFO > M, R contains the row scale factors */
/* >          for A. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0,  C contains the column scale factors for A. */
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
/* >          > 0:  if INFO = i,  and i is */
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
/* Subroutine */ int zgbequb_(integer *m, integer *n, integer *kl, integer *
	ku, doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *
	c__, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, 
	integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double log(doublereal), d_imag(doublecomplex *), pow_di(doublereal *, 
	    integer *);

    /* Local variables */
    static integer i__, j, kd;
    static doublereal radix, rcmin, rcmax;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, logrdx, smlnum;


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

/*     Test the input parameters. */

#line 210 "zgbequb.f"
    /* Parameter adjustments */
#line 210 "zgbequb.f"
    ab_dim1 = *ldab;
#line 210 "zgbequb.f"
    ab_offset = 1 + ab_dim1;
#line 210 "zgbequb.f"
    ab -= ab_offset;
#line 210 "zgbequb.f"
    --r__;
#line 210 "zgbequb.f"
    --c__;
#line 210 "zgbequb.f"

#line 210 "zgbequb.f"
    /* Function Body */
#line 210 "zgbequb.f"
    *info = 0;
#line 211 "zgbequb.f"
    if (*m < 0) {
#line 212 "zgbequb.f"
	*info = -1;
#line 213 "zgbequb.f"
    } else if (*n < 0) {
#line 214 "zgbequb.f"
	*info = -2;
#line 215 "zgbequb.f"
    } else if (*kl < 0) {
#line 216 "zgbequb.f"
	*info = -3;
#line 217 "zgbequb.f"
    } else if (*ku < 0) {
#line 218 "zgbequb.f"
	*info = -4;
#line 219 "zgbequb.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 220 "zgbequb.f"
	*info = -6;
#line 221 "zgbequb.f"
    }
#line 222 "zgbequb.f"
    if (*info != 0) {
#line 223 "zgbequb.f"
	i__1 = -(*info);
#line 223 "zgbequb.f"
	xerbla_("ZGBEQUB", &i__1, (ftnlen)7);
#line 224 "zgbequb.f"
	return 0;
#line 225 "zgbequb.f"
    }

/*     Quick return if possible. */

#line 229 "zgbequb.f"
    if (*m == 0 || *n == 0) {
#line 230 "zgbequb.f"
	*rowcnd = 1.;
#line 231 "zgbequb.f"
	*colcnd = 1.;
#line 232 "zgbequb.f"
	*amax = 0.;
#line 233 "zgbequb.f"
	return 0;
#line 234 "zgbequb.f"
    }

/*     Get machine constants.  Assume SMLNUM is a power of the radix. */

#line 238 "zgbequb.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 239 "zgbequb.f"
    bignum = 1. / smlnum;
#line 240 "zgbequb.f"
    radix = dlamch_("B", (ftnlen)1);
#line 241 "zgbequb.f"
    logrdx = log(radix);

/*     Compute row scale factors. */

#line 245 "zgbequb.f"
    i__1 = *m;
#line 245 "zgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "zgbequb.f"
	r__[i__] = 0.;
#line 247 "zgbequb.f"
/* L10: */
#line 247 "zgbequb.f"
    }

/*     Find the maximum element in each row. */

#line 251 "zgbequb.f"
    kd = *ku + 1;
#line 252 "zgbequb.f"
    i__1 = *n;
#line 252 "zgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 253 "zgbequb.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 253 "zgbequb.f"
	i__4 = j + *kl;
#line 253 "zgbequb.f"
	i__3 = min(i__4,*m);
#line 253 "zgbequb.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 254 "zgbequb.f"
	    i__2 = kd + i__ - j + j * ab_dim1;
#line 254 "zgbequb.f"
	    d__3 = r__[i__], d__4 = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2));
#line 254 "zgbequb.f"
	    r__[i__] = max(d__3,d__4);
#line 255 "zgbequb.f"
/* L20: */
#line 255 "zgbequb.f"
	}
#line 256 "zgbequb.f"
/* L30: */
#line 256 "zgbequb.f"
    }
#line 257 "zgbequb.f"
    i__1 = *m;
#line 257 "zgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "zgbequb.f"
	if (r__[i__] > 0.) {
#line 259 "zgbequb.f"
	    i__3 = (integer) (log(r__[i__]) / logrdx);
#line 259 "zgbequb.f"
	    r__[i__] = pow_di(&radix, &i__3);
#line 260 "zgbequb.f"
	}
#line 261 "zgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 265 "zgbequb.f"
    rcmin = bignum;
#line 266 "zgbequb.f"
    rcmax = 0.;
#line 267 "zgbequb.f"
    i__1 = *m;
#line 267 "zgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 268 "zgbequb.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 268 "zgbequb.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 269 "zgbequb.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 269 "zgbequb.f"
	rcmin = min(d__1,d__2);
#line 270 "zgbequb.f"
/* L40: */
#line 270 "zgbequb.f"
    }
#line 271 "zgbequb.f"
    *amax = rcmax;

#line 273 "zgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 277 "zgbequb.f"
	i__1 = *m;
#line 277 "zgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "zgbequb.f"
	    if (r__[i__] == 0.) {
#line 279 "zgbequb.f"
		*info = i__;
#line 280 "zgbequb.f"
		return 0;
#line 281 "zgbequb.f"
	    }
#line 282 "zgbequb.f"
/* L50: */
#line 282 "zgbequb.f"
	}
#line 283 "zgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 287 "zgbequb.f"
	i__1 = *m;
#line 287 "zgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 288 "zgbequb.f"
	    d__2 = r__[i__];
#line 288 "zgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 288 "zgbequb.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 289 "zgbequb.f"
/* L60: */
#line 289 "zgbequb.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)). */

#line 293 "zgbequb.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 294 "zgbequb.f"
    }

/*     Compute column scale factors. */

#line 298 "zgbequb.f"
    i__1 = *n;
#line 298 "zgbequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 299 "zgbequb.f"
	c__[j] = 0.;
#line 300 "zgbequb.f"
/* L70: */
#line 300 "zgbequb.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 305 "zgbequb.f"
    i__1 = *n;
#line 305 "zgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 306 "zgbequb.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 306 "zgbequb.f"
	i__4 = j + *kl;
#line 306 "zgbequb.f"
	i__2 = min(i__4,*m);
#line 306 "zgbequb.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 307 "zgbequb.f"
	    i__3 = kd + i__ - j + j * ab_dim1;
#line 307 "zgbequb.f"
	    d__3 = c__[j], d__4 = ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2))) * 
		    r__[i__];
#line 307 "zgbequb.f"
	    c__[j] = max(d__3,d__4);
#line 308 "zgbequb.f"
/* L80: */
#line 308 "zgbequb.f"
	}
#line 309 "zgbequb.f"
	if (c__[j] > 0.) {
#line 310 "zgbequb.f"
	    i__2 = (integer) (log(c__[j]) / logrdx);
#line 310 "zgbequb.f"
	    c__[j] = pow_di(&radix, &i__2);
#line 311 "zgbequb.f"
	}
#line 312 "zgbequb.f"
/* L90: */
#line 312 "zgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 316 "zgbequb.f"
    rcmin = bignum;
#line 317 "zgbequb.f"
    rcmax = 0.;
#line 318 "zgbequb.f"
    i__1 = *n;
#line 318 "zgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 319 "zgbequb.f"
	d__1 = rcmin, d__2 = c__[j];
#line 319 "zgbequb.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 320 "zgbequb.f"
	d__1 = rcmax, d__2 = c__[j];
#line 320 "zgbequb.f"
	rcmax = max(d__1,d__2);
#line 321 "zgbequb.f"
/* L100: */
#line 321 "zgbequb.f"
    }

#line 323 "zgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 327 "zgbequb.f"
	i__1 = *n;
#line 327 "zgbequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 328 "zgbequb.f"
	    if (c__[j] == 0.) {
#line 329 "zgbequb.f"
		*info = *m + j;
#line 330 "zgbequb.f"
		return 0;
#line 331 "zgbequb.f"
	    }
#line 332 "zgbequb.f"
/* L110: */
#line 332 "zgbequb.f"
	}
#line 333 "zgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 337 "zgbequb.f"
	i__1 = *n;
#line 337 "zgbequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 338 "zgbequb.f"
	    d__2 = c__[j];
#line 338 "zgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 338 "zgbequb.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 339 "zgbequb.f"
/* L120: */
#line 339 "zgbequb.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)). */

#line 343 "zgbequb.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 344 "zgbequb.f"
    }

#line 346 "zgbequb.f"
    return 0;

/*     End of ZGBEQUB */

} /* zgbequb_ */

