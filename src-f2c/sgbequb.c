#line 1 "sgbequb.f"
/* sgbequb.f -- translated by f2c (version 20100827).
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

#line 1 "sgbequb.f"
/* > \brief \b SGBEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                           AMAX, INFO ) */

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
/* > SGBEQUB computes row and column scalings intended to equilibrate an */
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
/* > This routine differs from SGEEQU by restricting the scaling factors */
/* > to a power of the radix.  Barring over- and underflow, scaling by */
/* > these factors introduces no additional rounding errors.  However, the */
/* > scaled entries' magnitudes are no longer approximately 1 but lie */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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
/* >          R is REAL array, dimension (M) */
/* >          If INFO = 0 or INFO > M, R contains the row scale factors */
/* >          for A. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >          If INFO = 0,  C contains the column scale factors for A. */
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

/* > \date June 2016 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgbequb_(integer *m, integer *n, integer *kl, integer *
	ku, doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, kd;
    static doublereal radix, rcmin, rcmax;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, logrdx, smlnum;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

/*     Test the input parameters. */

#line 200 "sgbequb.f"
    /* Parameter adjustments */
#line 200 "sgbequb.f"
    ab_dim1 = *ldab;
#line 200 "sgbequb.f"
    ab_offset = 1 + ab_dim1;
#line 200 "sgbequb.f"
    ab -= ab_offset;
#line 200 "sgbequb.f"
    --r__;
#line 200 "sgbequb.f"
    --c__;
#line 200 "sgbequb.f"

#line 200 "sgbequb.f"
    /* Function Body */
#line 200 "sgbequb.f"
    *info = 0;
#line 201 "sgbequb.f"
    if (*m < 0) {
#line 202 "sgbequb.f"
	*info = -1;
#line 203 "sgbequb.f"
    } else if (*n < 0) {
#line 204 "sgbequb.f"
	*info = -2;
#line 205 "sgbequb.f"
    } else if (*kl < 0) {
#line 206 "sgbequb.f"
	*info = -3;
#line 207 "sgbequb.f"
    } else if (*ku < 0) {
#line 208 "sgbequb.f"
	*info = -4;
#line 209 "sgbequb.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 210 "sgbequb.f"
	*info = -6;
#line 211 "sgbequb.f"
    }
#line 212 "sgbequb.f"
    if (*info != 0) {
#line 213 "sgbequb.f"
	i__1 = -(*info);
#line 213 "sgbequb.f"
	xerbla_("SGBEQUB", &i__1, (ftnlen)7);
#line 214 "sgbequb.f"
	return 0;
#line 215 "sgbequb.f"
    }

/*     Quick return if possible. */

#line 219 "sgbequb.f"
    if (*m == 0 || *n == 0) {
#line 220 "sgbequb.f"
	*rowcnd = 1.;
#line 221 "sgbequb.f"
	*colcnd = 1.;
#line 222 "sgbequb.f"
	*amax = 0.;
#line 223 "sgbequb.f"
	return 0;
#line 224 "sgbequb.f"
    }

/*     Get machine constants.  Assume SMLNUM is a power of the radix. */

#line 228 "sgbequb.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 229 "sgbequb.f"
    bignum = 1. / smlnum;
#line 230 "sgbequb.f"
    radix = slamch_("B", (ftnlen)1);
#line 231 "sgbequb.f"
    logrdx = log(radix);

/*     Compute row scale factors. */

#line 235 "sgbequb.f"
    i__1 = *m;
#line 235 "sgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "sgbequb.f"
	r__[i__] = 0.;
#line 237 "sgbequb.f"
/* L10: */
#line 237 "sgbequb.f"
    }

/*     Find the maximum element in each row. */

#line 241 "sgbequb.f"
    kd = *ku + 1;
#line 242 "sgbequb.f"
    i__1 = *n;
#line 242 "sgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 243 "sgbequb.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 243 "sgbequb.f"
	i__4 = j + *kl;
#line 243 "sgbequb.f"
	i__3 = min(i__4,*m);
#line 243 "sgbequb.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 244 "sgbequb.f"
	    d__2 = r__[i__], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], 
		    abs(d__1));
#line 244 "sgbequb.f"
	    r__[i__] = max(d__2,d__3);
#line 245 "sgbequb.f"
/* L20: */
#line 245 "sgbequb.f"
	}
#line 246 "sgbequb.f"
/* L30: */
#line 246 "sgbequb.f"
    }
#line 247 "sgbequb.f"
    i__1 = *m;
#line 247 "sgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 248 "sgbequb.f"
	if (r__[i__] > 0.) {
#line 249 "sgbequb.f"
	    i__3 = (integer) (log(r__[i__]) / logrdx);
#line 249 "sgbequb.f"
	    r__[i__] = pow_di(&radix, &i__3);
#line 250 "sgbequb.f"
	}
#line 251 "sgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 255 "sgbequb.f"
    rcmin = bignum;
#line 256 "sgbequb.f"
    rcmax = 0.;
#line 257 "sgbequb.f"
    i__1 = *m;
#line 257 "sgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 258 "sgbequb.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 258 "sgbequb.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 259 "sgbequb.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 259 "sgbequb.f"
	rcmin = min(d__1,d__2);
#line 260 "sgbequb.f"
/* L40: */
#line 260 "sgbequb.f"
    }
#line 261 "sgbequb.f"
    *amax = rcmax;

#line 263 "sgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 267 "sgbequb.f"
	i__1 = *m;
#line 267 "sgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 268 "sgbequb.f"
	    if (r__[i__] == 0.) {
#line 269 "sgbequb.f"
		*info = i__;
#line 270 "sgbequb.f"
		return 0;
#line 271 "sgbequb.f"
	    }
#line 272 "sgbequb.f"
/* L50: */
#line 272 "sgbequb.f"
	}
#line 273 "sgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 277 "sgbequb.f"
	i__1 = *m;
#line 277 "sgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 278 "sgbequb.f"
	    d__2 = r__[i__];
#line 278 "sgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 278 "sgbequb.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 279 "sgbequb.f"
/* L60: */
#line 279 "sgbequb.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)). */

#line 283 "sgbequb.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 284 "sgbequb.f"
    }

/*     Compute column scale factors. */

#line 288 "sgbequb.f"
    i__1 = *n;
#line 288 "sgbequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 289 "sgbequb.f"
	c__[j] = 0.;
#line 290 "sgbequb.f"
/* L70: */
#line 290 "sgbequb.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 295 "sgbequb.f"
    i__1 = *n;
#line 295 "sgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 296 "sgbequb.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 296 "sgbequb.f"
	i__4 = j + *kl;
#line 296 "sgbequb.f"
	i__2 = min(i__4,*m);
#line 296 "sgbequb.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 297 "sgbequb.f"
	    d__2 = c__[j], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
		    d__1)) * r__[i__];
#line 297 "sgbequb.f"
	    c__[j] = max(d__2,d__3);
#line 298 "sgbequb.f"
/* L80: */
#line 298 "sgbequb.f"
	}
#line 299 "sgbequb.f"
	if (c__[j] > 0.) {
#line 300 "sgbequb.f"
	    i__2 = (integer) (log(c__[j]) / logrdx);
#line 300 "sgbequb.f"
	    c__[j] = pow_di(&radix, &i__2);
#line 301 "sgbequb.f"
	}
#line 302 "sgbequb.f"
/* L90: */
#line 302 "sgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 306 "sgbequb.f"
    rcmin = bignum;
#line 307 "sgbequb.f"
    rcmax = 0.;
#line 308 "sgbequb.f"
    i__1 = *n;
#line 308 "sgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 309 "sgbequb.f"
	d__1 = rcmin, d__2 = c__[j];
#line 309 "sgbequb.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 310 "sgbequb.f"
	d__1 = rcmax, d__2 = c__[j];
#line 310 "sgbequb.f"
	rcmax = max(d__1,d__2);
#line 311 "sgbequb.f"
/* L100: */
#line 311 "sgbequb.f"
    }

#line 313 "sgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 317 "sgbequb.f"
	i__1 = *n;
#line 317 "sgbequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 318 "sgbequb.f"
	    if (c__[j] == 0.) {
#line 319 "sgbequb.f"
		*info = *m + j;
#line 320 "sgbequb.f"
		return 0;
#line 321 "sgbequb.f"
	    }
#line 322 "sgbequb.f"
/* L110: */
#line 322 "sgbequb.f"
	}
#line 323 "sgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 327 "sgbequb.f"
	i__1 = *n;
#line 327 "sgbequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 328 "sgbequb.f"
	    d__2 = c__[j];
#line 328 "sgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 328 "sgbequb.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 329 "sgbequb.f"
/* L120: */
#line 329 "sgbequb.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)). */

#line 333 "sgbequb.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 334 "sgbequb.f"
    }

#line 336 "sgbequb.f"
    return 0;

/*     End of SGBEQUB */

} /* sgbequb_ */

