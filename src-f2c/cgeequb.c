#line 1 "cgeequb.f"
/* cgeequb.f -- translated by f2c (version 20100827).
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

#line 1 "cgeequb.f"
/* > \brief \b CGEEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                           INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ), R( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEEQUB computes row and column scalings intended to equilibrate an */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The M-by-N matrix whose equilibration factors are */
/* >          to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
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

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgeequb_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double log(doublereal), d_imag(doublecomplex *), pow_di(doublereal *, 
	    integer *);

    /* Local variables */
    static integer i__, j;
    static doublereal radix, rcmin, rcmax;
    extern doublereal slamch_(char *, ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 195 "cgeequb.f"
    /* Parameter adjustments */
#line 195 "cgeequb.f"
    a_dim1 = *lda;
#line 195 "cgeequb.f"
    a_offset = 1 + a_dim1;
#line 195 "cgeequb.f"
    a -= a_offset;
#line 195 "cgeequb.f"
    --r__;
#line 195 "cgeequb.f"
    --c__;
#line 195 "cgeequb.f"

#line 195 "cgeequb.f"
    /* Function Body */
#line 195 "cgeequb.f"
    *info = 0;
#line 196 "cgeequb.f"
    if (*m < 0) {
#line 197 "cgeequb.f"
	*info = -1;
#line 198 "cgeequb.f"
    } else if (*n < 0) {
#line 199 "cgeequb.f"
	*info = -2;
#line 200 "cgeequb.f"
    } else if (*lda < max(1,*m)) {
#line 201 "cgeequb.f"
	*info = -4;
#line 202 "cgeequb.f"
    }
#line 203 "cgeequb.f"
    if (*info != 0) {
#line 204 "cgeequb.f"
	i__1 = -(*info);
#line 204 "cgeequb.f"
	xerbla_("CGEEQUB", &i__1, (ftnlen)7);
#line 205 "cgeequb.f"
	return 0;
#line 206 "cgeequb.f"
    }

/*     Quick return if possible. */

#line 210 "cgeequb.f"
    if (*m == 0 || *n == 0) {
#line 211 "cgeequb.f"
	*rowcnd = 1.;
#line 212 "cgeequb.f"
	*colcnd = 1.;
#line 213 "cgeequb.f"
	*amax = 0.;
#line 214 "cgeequb.f"
	return 0;
#line 215 "cgeequb.f"
    }

/*     Get machine constants.  Assume SMLNUM is a power of the radix. */

#line 219 "cgeequb.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 220 "cgeequb.f"
    bignum = 1. / smlnum;
#line 221 "cgeequb.f"
    radix = slamch_("B", (ftnlen)1);
#line 222 "cgeequb.f"
    logrdx = log(radix);

/*     Compute row scale factors. */

#line 226 "cgeequb.f"
    i__1 = *m;
#line 226 "cgeequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 227 "cgeequb.f"
	r__[i__] = 0.;
#line 228 "cgeequb.f"
/* L10: */
#line 228 "cgeequb.f"
    }

/*     Find the maximum element in each row. */

#line 232 "cgeequb.f"
    i__1 = *n;
#line 232 "cgeequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 233 "cgeequb.f"
	i__2 = *m;
#line 233 "cgeequb.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 234 "cgeequb.f"
	    i__3 = i__ + j * a_dim1;
#line 234 "cgeequb.f"
	    d__3 = r__[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 234 "cgeequb.f"
	    r__[i__] = max(d__3,d__4);
#line 235 "cgeequb.f"
/* L20: */
#line 235 "cgeequb.f"
	}
#line 236 "cgeequb.f"
/* L30: */
#line 236 "cgeequb.f"
    }
#line 237 "cgeequb.f"
    i__1 = *m;
#line 237 "cgeequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 238 "cgeequb.f"
	if (r__[i__] > 0.) {
#line 239 "cgeequb.f"
	    i__2 = (integer) (log(r__[i__]) / logrdx);
#line 239 "cgeequb.f"
	    r__[i__] = pow_di(&radix, &i__2);
#line 240 "cgeequb.f"
	}
#line 241 "cgeequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 245 "cgeequb.f"
    rcmin = bignum;
#line 246 "cgeequb.f"
    rcmax = 0.;
#line 247 "cgeequb.f"
    i__1 = *m;
#line 247 "cgeequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 248 "cgeequb.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 248 "cgeequb.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 249 "cgeequb.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 249 "cgeequb.f"
	rcmin = min(d__1,d__2);
#line 250 "cgeequb.f"
/* L40: */
#line 250 "cgeequb.f"
    }
#line 251 "cgeequb.f"
    *amax = rcmax;

#line 253 "cgeequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 257 "cgeequb.f"
	i__1 = *m;
#line 257 "cgeequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "cgeequb.f"
	    if (r__[i__] == 0.) {
#line 259 "cgeequb.f"
		*info = i__;
#line 260 "cgeequb.f"
		return 0;
#line 261 "cgeequb.f"
	    }
#line 262 "cgeequb.f"
/* L50: */
#line 262 "cgeequb.f"
	}
#line 263 "cgeequb.f"
    } else {

/*        Invert the scale factors. */

#line 267 "cgeequb.f"
	i__1 = *m;
#line 267 "cgeequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 268 "cgeequb.f"
	    d__2 = r__[i__];
#line 268 "cgeequb.f"
	    d__1 = max(d__2,smlnum);
#line 268 "cgeequb.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 269 "cgeequb.f"
/* L60: */
#line 269 "cgeequb.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)). */

#line 273 "cgeequb.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 274 "cgeequb.f"
    }

/*     Compute column scale factors. */

#line 278 "cgeequb.f"
    i__1 = *n;
#line 278 "cgeequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 279 "cgeequb.f"
	c__[j] = 0.;
#line 280 "cgeequb.f"
/* L70: */
#line 280 "cgeequb.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 285 "cgeequb.f"
    i__1 = *n;
#line 285 "cgeequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 286 "cgeequb.f"
	i__2 = *m;
#line 286 "cgeequb.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 287 "cgeequb.f"
	    i__3 = i__ + j * a_dim1;
#line 287 "cgeequb.f"
	    d__3 = c__[j], d__4 = ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2))) * r__[i__];
#line 287 "cgeequb.f"
	    c__[j] = max(d__3,d__4);
#line 288 "cgeequb.f"
/* L80: */
#line 288 "cgeequb.f"
	}
#line 289 "cgeequb.f"
	if (c__[j] > 0.) {
#line 290 "cgeequb.f"
	    i__2 = (integer) (log(c__[j]) / logrdx);
#line 290 "cgeequb.f"
	    c__[j] = pow_di(&radix, &i__2);
#line 291 "cgeequb.f"
	}
#line 292 "cgeequb.f"
/* L90: */
#line 292 "cgeequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 296 "cgeequb.f"
    rcmin = bignum;
#line 297 "cgeequb.f"
    rcmax = 0.;
#line 298 "cgeequb.f"
    i__1 = *n;
#line 298 "cgeequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 299 "cgeequb.f"
	d__1 = rcmin, d__2 = c__[j];
#line 299 "cgeequb.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 300 "cgeequb.f"
	d__1 = rcmax, d__2 = c__[j];
#line 300 "cgeequb.f"
	rcmax = max(d__1,d__2);
#line 301 "cgeequb.f"
/* L100: */
#line 301 "cgeequb.f"
    }

#line 303 "cgeequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 307 "cgeequb.f"
	i__1 = *n;
#line 307 "cgeequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 308 "cgeequb.f"
	    if (c__[j] == 0.) {
#line 309 "cgeequb.f"
		*info = *m + j;
#line 310 "cgeequb.f"
		return 0;
#line 311 "cgeequb.f"
	    }
#line 312 "cgeequb.f"
/* L110: */
#line 312 "cgeequb.f"
	}
#line 313 "cgeequb.f"
    } else {

/*        Invert the scale factors. */

#line 317 "cgeequb.f"
	i__1 = *n;
#line 317 "cgeequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 318 "cgeequb.f"
	    d__2 = c__[j];
#line 318 "cgeequb.f"
	    d__1 = max(d__2,smlnum);
#line 318 "cgeequb.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 319 "cgeequb.f"
/* L120: */
#line 319 "cgeequb.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)). */

#line 323 "cgeequb.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 324 "cgeequb.f"
    }

#line 326 "cgeequb.f"
    return 0;

/*     End of CGEEQUB */

} /* cgeequb_ */

