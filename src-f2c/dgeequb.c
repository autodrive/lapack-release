#line 1 "dgeequb.f"
/* dgeequb.f -- translated by f2c (version 20100827).
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

#line 1 "dgeequb.f"
/* > \brief \b DGEEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                           INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       DOUBLE PRECISION   AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEEQUB computes row and column scalings intended to equilibrate an */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgeequb_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j;
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

#line 186 "dgeequb.f"
    /* Parameter adjustments */
#line 186 "dgeequb.f"
    a_dim1 = *lda;
#line 186 "dgeequb.f"
    a_offset = 1 + a_dim1;
#line 186 "dgeequb.f"
    a -= a_offset;
#line 186 "dgeequb.f"
    --r__;
#line 186 "dgeequb.f"
    --c__;
#line 186 "dgeequb.f"

#line 186 "dgeequb.f"
    /* Function Body */
#line 186 "dgeequb.f"
    *info = 0;
#line 187 "dgeequb.f"
    if (*m < 0) {
#line 188 "dgeequb.f"
	*info = -1;
#line 189 "dgeequb.f"
    } else if (*n < 0) {
#line 190 "dgeequb.f"
	*info = -2;
#line 191 "dgeequb.f"
    } else if (*lda < max(1,*m)) {
#line 192 "dgeequb.f"
	*info = -4;
#line 193 "dgeequb.f"
    }
#line 194 "dgeequb.f"
    if (*info != 0) {
#line 195 "dgeequb.f"
	i__1 = -(*info);
#line 195 "dgeequb.f"
	xerbla_("DGEEQUB", &i__1, (ftnlen)7);
#line 196 "dgeequb.f"
	return 0;
#line 197 "dgeequb.f"
    }

/*     Quick return if possible. */

#line 201 "dgeequb.f"
    if (*m == 0 || *n == 0) {
#line 202 "dgeequb.f"
	*rowcnd = 1.;
#line 203 "dgeequb.f"
	*colcnd = 1.;
#line 204 "dgeequb.f"
	*amax = 0.;
#line 205 "dgeequb.f"
	return 0;
#line 206 "dgeequb.f"
    }

/*     Get machine constants.  Assume SMLNUM is a power of the radix. */

#line 210 "dgeequb.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 211 "dgeequb.f"
    bignum = 1. / smlnum;
#line 212 "dgeequb.f"
    radix = dlamch_("B", (ftnlen)1);
#line 213 "dgeequb.f"
    logrdx = log(radix);

/*     Compute row scale factors. */

#line 217 "dgeequb.f"
    i__1 = *m;
#line 217 "dgeequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "dgeequb.f"
	r__[i__] = 0.;
#line 219 "dgeequb.f"
/* L10: */
#line 219 "dgeequb.f"
    }

/*     Find the maximum element in each row. */

#line 223 "dgeequb.f"
    i__1 = *n;
#line 223 "dgeequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 224 "dgeequb.f"
	i__2 = *m;
#line 224 "dgeequb.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 225 "dgeequb.f"
	    d__2 = r__[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 225 "dgeequb.f"
	    r__[i__] = max(d__2,d__3);
#line 226 "dgeequb.f"
/* L20: */
#line 226 "dgeequb.f"
	}
#line 227 "dgeequb.f"
/* L30: */
#line 227 "dgeequb.f"
    }
#line 228 "dgeequb.f"
    i__1 = *m;
#line 228 "dgeequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "dgeequb.f"
	if (r__[i__] > 0.) {
#line 230 "dgeequb.f"
	    i__2 = (integer) (log(r__[i__]) / logrdx);
#line 230 "dgeequb.f"
	    r__[i__] = pow_di(&radix, &i__2);
#line 231 "dgeequb.f"
	}
#line 232 "dgeequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 236 "dgeequb.f"
    rcmin = bignum;
#line 237 "dgeequb.f"
    rcmax = 0.;
#line 238 "dgeequb.f"
    i__1 = *m;
#line 238 "dgeequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 239 "dgeequb.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 239 "dgeequb.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 240 "dgeequb.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 240 "dgeequb.f"
	rcmin = min(d__1,d__2);
#line 241 "dgeequb.f"
/* L40: */
#line 241 "dgeequb.f"
    }
#line 242 "dgeequb.f"
    *amax = rcmax;

#line 244 "dgeequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 248 "dgeequb.f"
	i__1 = *m;
#line 248 "dgeequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 249 "dgeequb.f"
	    if (r__[i__] == 0.) {
#line 250 "dgeequb.f"
		*info = i__;
#line 251 "dgeequb.f"
		return 0;
#line 252 "dgeequb.f"
	    }
#line 253 "dgeequb.f"
/* L50: */
#line 253 "dgeequb.f"
	}
#line 254 "dgeequb.f"
    } else {

/*        Invert the scale factors. */

#line 258 "dgeequb.f"
	i__1 = *m;
#line 258 "dgeequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 259 "dgeequb.f"
	    d__2 = r__[i__];
#line 259 "dgeequb.f"
	    d__1 = max(d__2,smlnum);
#line 259 "dgeequb.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 260 "dgeequb.f"
/* L60: */
#line 260 "dgeequb.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)). */

#line 264 "dgeequb.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 265 "dgeequb.f"
    }

/*     Compute column scale factors */

#line 269 "dgeequb.f"
    i__1 = *n;
#line 269 "dgeequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 270 "dgeequb.f"
	c__[j] = 0.;
#line 271 "dgeequb.f"
/* L70: */
#line 271 "dgeequb.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 276 "dgeequb.f"
    i__1 = *n;
#line 276 "dgeequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 277 "dgeequb.f"
	i__2 = *m;
#line 277 "dgeequb.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 278 "dgeequb.f"
	    d__2 = c__[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1)) * 
		    r__[i__];
#line 278 "dgeequb.f"
	    c__[j] = max(d__2,d__3);
#line 279 "dgeequb.f"
/* L80: */
#line 279 "dgeequb.f"
	}
#line 280 "dgeequb.f"
	if (c__[j] > 0.) {
#line 281 "dgeequb.f"
	    i__2 = (integer) (log(c__[j]) / logrdx);
#line 281 "dgeequb.f"
	    c__[j] = pow_di(&radix, &i__2);
#line 282 "dgeequb.f"
	}
#line 283 "dgeequb.f"
/* L90: */
#line 283 "dgeequb.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 287 "dgeequb.f"
    rcmin = bignum;
#line 288 "dgeequb.f"
    rcmax = 0.;
#line 289 "dgeequb.f"
    i__1 = *n;
#line 289 "dgeequb.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 290 "dgeequb.f"
	d__1 = rcmin, d__2 = c__[j];
#line 290 "dgeequb.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 291 "dgeequb.f"
	d__1 = rcmax, d__2 = c__[j];
#line 291 "dgeequb.f"
	rcmax = max(d__1,d__2);
#line 292 "dgeequb.f"
/* L100: */
#line 292 "dgeequb.f"
    }

#line 294 "dgeequb.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 298 "dgeequb.f"
	i__1 = *n;
#line 298 "dgeequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 299 "dgeequb.f"
	    if (c__[j] == 0.) {
#line 300 "dgeequb.f"
		*info = *m + j;
#line 301 "dgeequb.f"
		return 0;
#line 302 "dgeequb.f"
	    }
#line 303 "dgeequb.f"
/* L110: */
#line 303 "dgeequb.f"
	}
#line 304 "dgeequb.f"
    } else {

/*        Invert the scale factors. */

#line 308 "dgeequb.f"
	i__1 = *n;
#line 308 "dgeequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 309 "dgeequb.f"
	    d__2 = c__[j];
#line 309 "dgeequb.f"
	    d__1 = max(d__2,smlnum);
#line 309 "dgeequb.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 310 "dgeequb.f"
/* L120: */
#line 310 "dgeequb.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)). */

#line 314 "dgeequb.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 315 "dgeequb.f"
    }

#line 317 "dgeequb.f"
    return 0;

/*     End of DGEEQUB */

} /* dgeequb_ */

