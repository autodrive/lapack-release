#line 1 "cgbequb.f"
/* cgbequb.f -- translated by f2c (version 20100827).
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

#line 1 "cgbequb.f"
/* > \brief \b CGBEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/*                           AMAX, INFO ) */

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
/* > CGBEQUB computes row and column scalings intended to equilibrate an */
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
/* > This routine differs from CGEEQU by restricting the scaling factors */
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

/* > \date November 2011 */

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgbequb_(integer *m, integer *n, integer *kl, integer *
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
    extern doublereal slamch_(char *, ftnlen);
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

#line 210 "cgbequb.f"
    /* Parameter adjustments */
#line 210 "cgbequb.f"
    ab_dim1 = *ldab;
#line 210 "cgbequb.f"
    ab_offset = 1 + ab_dim1;
#line 210 "cgbequb.f"
    ab -= ab_offset;
#line 210 "cgbequb.f"
    --r__;
#line 210 "cgbequb.f"
    --c__;
#line 210 "cgbequb.f"

#line 210 "cgbequb.f"
    /* Function Body */
#line 210 "cgbequb.f"
    *info = 0;
#line 211 "cgbequb.f"
    if (*m < 0) {
#line 212 "cgbequb.f"
	*info = -1;
#line 213 "cgbequb.f"
    } else if (*n < 0) {
#line 214 "cgbequb.f"
	*info = -2;
#line 215 "cgbequb.f"
    } else if (*kl < 0) {
#line 216 "cgbequb.f"
	*info = -3;
#line 217 "cgbequb.f"
    } else if (*ku < 0) {
#line 218 "cgbequb.f"
	*info = -4;
#line 219 "cgbequb.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 220 "cgbequb.f"
	*info = -6;
#line 221 "cgbequb.f"
    }
#line 222 "cgbequb.f"
    if (*info != 0) {
#line 223 "cgbequb.f"
	i__1 = -(*info);
#line 223 "cgbequb.f"
	xerbla_("CGBEQUB", &i__1, (ftnlen)7);
#line 224 "cgbequb.f"
	return 0;
#line 225 "cgbequb.f"
    }

/*     Quick return if possible. */

#line 229 "cgbequb.f"
    if (*m == 0 || *n == 0) {
#line 230 "cgbequb.f"
	*rowcnd = 1.;
#line 231 "cgbequb.f"
	*colcnd = 1.;
#line 232 "cgbequb.f"
	*amax = 0.;
#line 233 "cgbequb.f"
	return 0;
#line 234 "cgbequb.f"
    }

/*     Get machine constants.  Assume SMLNUM is a power of the radix. */

#line 238 "cgbequb.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 239 "cgbequb.f"
    bignum = 1. / smlnum;
#line 240 "cgbequb.f"
    radix = slamch_("B", (ftnlen)1);
#line 241 "cgbequb.f"
    logrdx = log(radix);

/*     Compute row scale factors. */

#line 245 "cgbequb.f"
    i__1 = *m;
#line 245 "cgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "cgbequb.f"
	r__[i__] = 0.;
#line 247 "cgbequb.f"
/* L10: */
#line 247 "cgbequb.f"
    }

/*     Find the maximum element in each row. */

#line 251 "cgbequb.f"
    kd = *ku + 1;
#line 252 "cgbequb.f"
    i__1 = *n;
#line 252 "cgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 253 "cgbequb.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 253 "cgbequb.f"
	i__4 = j + *kl;
#line 253 "cgbequb.f"
	i__3 = min(i__4,*m);
#line 253 "cgbequb.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 254 "cgbequb.f"
	    i__2 = kd + i__ - j + j * ab_dim1;
#line 254 "cgbequb.f"
	    d__3 = r__[i__], d__4 = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2));
#line 254 "cgbequb.f"
	    r__[i__] = max(d__3,d__4);
#line 255 "cgbequb.f"
/* L20: */
#line 255 "cgbequb.f"
	}
#line 256 "cgbequb.f"
/* L30: */
#line 256 "cgbequb.f"
    }
#line 257 "cgbequb.f"
    i__1 = *m;
#line 257 "cgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "cgbequb.f"
	if (r__[i__] > 0.) {
#line 259 "cgbequb.f"
	    i__3 = (integer) (log(r__[i__]) / logrdx);
#line 259 "cgbequb.f"
	    r__[i__] = pow_di(&radix, &i__3);
#line 260 "cgbequb.f"
	}
#line 261 "cgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 265 "cgbequb.f"
    rcmin = bignum;
#line 266 "cgbequb.f"
    rcmax = 0.;
#line 267 "cgbequb.f"
    i__1 = *m;
#line 267 "cgbequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 268 "cgbequb.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 268 "cgbequb.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 269 "cgbequb.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 269 "cgbequb.f"
	rcmin = min(d__1,d__2);
#line 270 "cgbequb.f"
/* L40: */
#line 270 "cgbequb.f"
    }
#line 271 "cgbequb.f"
    *amax = rcmax;

#line 273 "cgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 277 "cgbequb.f"
	i__1 = *m;
#line 277 "cgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "cgbequb.f"
	    if (r__[i__] == 0.) {
#line 279 "cgbequb.f"
		*info = i__;
#line 280 "cgbequb.f"
		return 0;
#line 281 "cgbequb.f"
	    }
#line 282 "cgbequb.f"
/* L50: */
#line 282 "cgbequb.f"
	}
#line 283 "cgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 287 "cgbequb.f"
	i__1 = *m;
#line 287 "cgbequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 288 "cgbequb.f"
	    d__2 = r__[i__];
#line 288 "cgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 288 "cgbequb.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 289 "cgbequb.f"
/* L60: */
#line 289 "cgbequb.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)). */

#line 293 "cgbequb.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 294 "cgbequb.f"
    }

/*     Compute column scale factors. */

#line 298 "cgbequb.f"
    i__1 = *n;
#line 298 "cgbequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 299 "cgbequb.f"
	c__[j] = 0.;
#line 300 "cgbequb.f"
/* L70: */
#line 300 "cgbequb.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 305 "cgbequb.f"
    i__1 = *n;
#line 305 "cgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 306 "cgbequb.f"
	i__3 = j - *ku;
/* Computing MIN */
#line 306 "cgbequb.f"
	i__4 = j + *kl;
#line 306 "cgbequb.f"
	i__2 = min(i__4,*m);
#line 306 "cgbequb.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 307 "cgbequb.f"
	    i__3 = kd + i__ - j + j * ab_dim1;
#line 307 "cgbequb.f"
	    d__3 = c__[j], d__4 = ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(d__2))) * 
		    r__[i__];
#line 307 "cgbequb.f"
	    c__[j] = max(d__3,d__4);
#line 308 "cgbequb.f"
/* L80: */
#line 308 "cgbequb.f"
	}
#line 309 "cgbequb.f"
	if (c__[j] > 0.) {
#line 310 "cgbequb.f"
	    i__2 = (integer) (log(c__[j]) / logrdx);
#line 310 "cgbequb.f"
	    c__[j] = pow_di(&radix, &i__2);
#line 311 "cgbequb.f"
	}
#line 312 "cgbequb.f"
/* L90: */
#line 312 "cgbequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 316 "cgbequb.f"
    rcmin = bignum;
#line 317 "cgbequb.f"
    rcmax = 0.;
#line 318 "cgbequb.f"
    i__1 = *n;
#line 318 "cgbequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 319 "cgbequb.f"
	d__1 = rcmin, d__2 = c__[j];
#line 319 "cgbequb.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 320 "cgbequb.f"
	d__1 = rcmax, d__2 = c__[j];
#line 320 "cgbequb.f"
	rcmax = max(d__1,d__2);
#line 321 "cgbequb.f"
/* L100: */
#line 321 "cgbequb.f"
    }

#line 323 "cgbequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 327 "cgbequb.f"
	i__1 = *n;
#line 327 "cgbequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 328 "cgbequb.f"
	    if (c__[j] == 0.) {
#line 329 "cgbequb.f"
		*info = *m + j;
#line 330 "cgbequb.f"
		return 0;
#line 331 "cgbequb.f"
	    }
#line 332 "cgbequb.f"
/* L110: */
#line 332 "cgbequb.f"
	}
#line 333 "cgbequb.f"
    } else {

/*        Invert the scale factors. */

#line 337 "cgbequb.f"
	i__1 = *n;
#line 337 "cgbequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 338 "cgbequb.f"
	    d__2 = c__[j];
#line 338 "cgbequb.f"
	    d__1 = max(d__2,smlnum);
#line 338 "cgbequb.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 339 "cgbequb.f"
/* L120: */
#line 339 "cgbequb.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)). */

#line 343 "cgbequb.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 344 "cgbequb.f"
    }

#line 346 "cgbequb.f"
    return 0;

/*     End of CGBEQUB */

} /* cgbequb_ */

