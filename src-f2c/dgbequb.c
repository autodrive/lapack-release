#line 1 "dgbequb.f"
/* dgbequb.f -- translated by f2c (version 20100827).
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

#line 1 "dgbequb.f"
/* > \brief \b DGBEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                           AMAX, INFO ) */

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
/* > DGBEQUB computes row and column scalings intended to equilibrate an */
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
/* > This routine differs from DGEEQU by restricting the scaling factors */
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

/* > \date December 2016 */

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgbequb_(integer *m, integer *n, integer *kl, integer *
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
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, logrdx, smlnum;


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

/*     Test the input parameters. */

#line 200 "dgbequb.f"
    /* Parameter adjustments */
#line 200 "dgbequb.f"
    ab_dim1 = *ldab;
#line 200 "dgbequb.f"
    ab_offset = 1 + ab_dim1;
#line 200 "dgbequb.f"
    ab -= ab_offset;
#line 200 "dgbequb.f"
    --r__;
#line 200 "dgbequb.f"
    --c__;
#line 200 "dgbequb.f"

#line 200 "dgbequb.f"
    /* Function Body */
#line 200 "dgbequb.f"
    *info = 0;
#line 201 "dgbequb.f"
    if (*m < 0) {
#line 202 "dgbequb.f"
	*info = -1;
#line 203 "dgbequb.f"
    } else if (*n < 0) {
#line 204 "dgbequb.f"
	*info = -2;
#line 205 "dgbequb.f"
    } else if (*kl < 0) {
#line 206 "dgbequb.f"
	*info = -3;
#line 207 "dgbequb.f"
    } else if (*ku < 0) {
#line 208 "dgbequb.f"
	*info = -4;
#line 209 "dgbequb.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 210 "dgbequb.f"
	*info = -6;
#line 211 "dgbequb.f"
    }
#line 212 "dgbequb.f"
    if (*info != 0) {
#line 213 "dgbequb.f"
	i__1 = -(*info);
#line 213 "dgbequb.f"
	xerbla_("DGBEQUB", &i__1, (ftnlen)7);
#line 214 "dgbequb.f"
	return 0;
#line 215 "dgbequb.f"
    }

/*     Quick return if possible. */

#line 219 "dgbequb.f"
    if (*m == 0 || *n == 0) {
#line 220 "dgbequb.f"
	*rowcnd = 1.;
#line 221 "dgbequb.f"
	*colcnd = 1.;
#line 222 "dgbequb.f"
	*amax = 0.;
#line 223 "dgbequb.f"
	return 0;
#line 224 "dgbequb.f"
    }

/*     Get machine constants.  Assume SMLNUM is a power of the radix. */

#line 228 "dgbequb.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 229 "dgbequb.f"
    bignum = 1. / smlnum;
#line 230 "dgbequb.f"
    radix = dlamch_("B", (ftnlen)1);
#line 231 "dgbequb.f"
    logrdx = log(radix);

/*     Compute row scale factors. */

#line 235 "dgbequb.f"
    i__1 = *m;
#line 235 "dgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "dgbequb.f"
	r__[i__] = 0.;
#line 237 "dgbequb.f"
/* L10: */
#line 237 "dgbequb.f"
    }

/*     Find the maximum element in each row. */

#line 241 "dgbequb.f"
    kd = *ku + 1;
#line 242 "dgbequb.f"
    i__1 = *n;
#line 242 "dgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 243 "dgbequb.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 243 "dgbequb.f"
	i__4 = j + *kl;
#line 243 "dgbequb.f"
	i__3 = min(i__4,*m);
#line 243 "dgbequb.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 244 "dgbequb.f"
	    d__2 = r__[i__], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], 
		    abs(d__1));
#line 244 "dgbequb.f"
	    r__[i__] = max(d__2,d__3);
#line 245 "dgbequb.f"
/* L20: */
#line 245 "dgbequb.f"
	}
#line 246 "dgbequb.f"
/* L30: */
#line 246 "dgbequb.f"
    }
#line 247 "dgbequb.f"
    i__1 = *m;
#line 247 "dgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 248 "dgbequb.f"
	if (r__[i__] > 0.) {
#line 249 "dgbequb.f"
	    i__3 = (integer) (log(r__[i__]) / logrdx);
#line 249 "dgbequb.f"
	    r__[i__] = pow_di(&radix, &i__3);
#line 250 "dgbequb.f"
	}
#line 251 "dgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 255 "dgbequb.f"
    rcmin = bignum;
#line 256 "dgbequb.f"
    rcmax = 0.;
#line 257 "dgbequb.f"
    i__1 = *m;
#line 257 "dgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 258 "dgbequb.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 258 "dgbequb.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 259 "dgbequb.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 259 "dgbequb.f"
	rcmin = min(d__1,d__2);
#line 260 "dgbequb.f"
/* L40: */
#line 260 "dgbequb.f"
    }
#line 261 "dgbequb.f"
    *amax = rcmax;

#line 263 "dgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 267 "dgbequb.f"
	i__1 = *m;
#line 267 "dgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 268 "dgbequb.f"
	    if (r__[i__] == 0.) {
#line 269 "dgbequb.f"
		*info = i__;
#line 270 "dgbequb.f"
		return 0;
#line 271 "dgbequb.f"
	    }
#line 272 "dgbequb.f"
/* L50: */
#line 272 "dgbequb.f"
	}
#line 273 "dgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 277 "dgbequb.f"
	i__1 = *m;
#line 277 "dgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 278 "dgbequb.f"
	    d__2 = r__[i__];
#line 278 "dgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 278 "dgbequb.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 279 "dgbequb.f"
/* L60: */
#line 279 "dgbequb.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)). */

#line 283 "dgbequb.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 284 "dgbequb.f"
    }

/*     Compute column scale factors. */

#line 288 "dgbequb.f"
    i__1 = *n;
#line 288 "dgbequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 289 "dgbequb.f"
	c__[j] = 0.;
#line 290 "dgbequb.f"
/* L70: */
#line 290 "dgbequb.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 295 "dgbequb.f"
    i__1 = *n;
#line 295 "dgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 296 "dgbequb.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 296 "dgbequb.f"
	i__4 = j + *kl;
#line 296 "dgbequb.f"
	i__2 = min(i__4,*m);
#line 296 "dgbequb.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 297 "dgbequb.f"
	    d__2 = c__[j], d__3 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
		    d__1)) * r__[i__];
#line 297 "dgbequb.f"
	    c__[j] = max(d__2,d__3);
#line 298 "dgbequb.f"
/* L80: */
#line 298 "dgbequb.f"
	}
#line 299 "dgbequb.f"
	if (c__[j] > 0.) {
#line 300 "dgbequb.f"
	    i__2 = (integer) (log(c__[j]) / logrdx);
#line 300 "dgbequb.f"
	    c__[j] = pow_di(&radix, &i__2);
#line 301 "dgbequb.f"
	}
#line 302 "dgbequb.f"
/* L90: */
#line 302 "dgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 306 "dgbequb.f"
    rcmin = bignum;
#line 307 "dgbequb.f"
    rcmax = 0.;
#line 308 "dgbequb.f"
    i__1 = *n;
#line 308 "dgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 309 "dgbequb.f"
	d__1 = rcmin, d__2 = c__[j];
#line 309 "dgbequb.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 310 "dgbequb.f"
	d__1 = rcmax, d__2 = c__[j];
#line 310 "dgbequb.f"
	rcmax = max(d__1,d__2);
#line 311 "dgbequb.f"
/* L100: */
#line 311 "dgbequb.f"
    }

#line 313 "dgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 317 "dgbequb.f"
	i__1 = *n;
#line 317 "dgbequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 318 "dgbequb.f"
	    if (c__[j] == 0.) {
#line 319 "dgbequb.f"
		*info = *m + j;
#line 320 "dgbequb.f"
		return 0;
#line 321 "dgbequb.f"
	    }
#line 322 "dgbequb.f"
/* L110: */
#line 322 "dgbequb.f"
	}
#line 323 "dgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 327 "dgbequb.f"
	i__1 = *n;
#line 327 "dgbequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 328 "dgbequb.f"
	    d__2 = c__[j];
#line 328 "dgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 328 "dgbequb.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 329 "dgbequb.f"
/* L120: */
#line 329 "dgbequb.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)). */

#line 333 "dgbequb.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 334 "dgbequb.f"
    }

#line 336 "dgbequb.f"
    return 0;

/*     End of DGBEQUB */

} /* dgbequb_ */

