#line 1 "zgeequ.f"
/* zgeequ.f -- translated by f2c (version 20100827).
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

#line 1 "zgeequ.f"
/* > \brief \b ZGEEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       DOUBLE PRECISION   AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), R( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEEQU computes row and column scalings intended to equilibrate an */
/* > M-by-N matrix A and reduce its condition number.  R returns the row */
/* > scale factors and C the column scale factors, chosen to try to make */
/* > the largest element in each row and column of the matrix B with */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgeequ_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 188 "zgeequ.f"
    /* Parameter adjustments */
#line 188 "zgeequ.f"
    a_dim1 = *lda;
#line 188 "zgeequ.f"
    a_offset = 1 + a_dim1;
#line 188 "zgeequ.f"
    a -= a_offset;
#line 188 "zgeequ.f"
    --r__;
#line 188 "zgeequ.f"
    --c__;
#line 188 "zgeequ.f"

#line 188 "zgeequ.f"
    /* Function Body */
#line 188 "zgeequ.f"
    *info = 0;
#line 189 "zgeequ.f"
    if (*m < 0) {
#line 190 "zgeequ.f"
	*info = -1;
#line 191 "zgeequ.f"
    } else if (*n < 0) {
#line 192 "zgeequ.f"
	*info = -2;
#line 193 "zgeequ.f"
    } else if (*lda < max(1,*m)) {
#line 194 "zgeequ.f"
	*info = -4;
#line 195 "zgeequ.f"
    }
#line 196 "zgeequ.f"
    if (*info != 0) {
#line 197 "zgeequ.f"
	i__1 = -(*info);
#line 197 "zgeequ.f"
	xerbla_("ZGEEQU", &i__1, (ftnlen)6);
#line 198 "zgeequ.f"
	return 0;
#line 199 "zgeequ.f"
    }

/*     Quick return if possible */

#line 203 "zgeequ.f"
    if (*m == 0 || *n == 0) {
#line 204 "zgeequ.f"
	*rowcnd = 1.;
#line 205 "zgeequ.f"
	*colcnd = 1.;
#line 206 "zgeequ.f"
	*amax = 0.;
#line 207 "zgeequ.f"
	return 0;
#line 208 "zgeequ.f"
    }

/*     Get machine constants. */

#line 212 "zgeequ.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 213 "zgeequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 217 "zgeequ.f"
    i__1 = *m;
#line 217 "zgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "zgeequ.f"
	r__[i__] = 0.;
#line 219 "zgeequ.f"
/* L10: */
#line 219 "zgeequ.f"
    }

/*     Find the maximum element in each row. */

#line 223 "zgeequ.f"
    i__1 = *n;
#line 223 "zgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 224 "zgeequ.f"
	i__2 = *m;
#line 224 "zgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 225 "zgeequ.f"
	    i__3 = i__ + j * a_dim1;
#line 225 "zgeequ.f"
	    d__3 = r__[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 225 "zgeequ.f"
	    r__[i__] = max(d__3,d__4);
#line 226 "zgeequ.f"
/* L20: */
#line 226 "zgeequ.f"
	}
#line 227 "zgeequ.f"
/* L30: */
#line 227 "zgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 231 "zgeequ.f"
    rcmin = bignum;
#line 232 "zgeequ.f"
    rcmax = 0.;
#line 233 "zgeequ.f"
    i__1 = *m;
#line 233 "zgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 234 "zgeequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 234 "zgeequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 235 "zgeequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 235 "zgeequ.f"
	rcmin = min(d__1,d__2);
#line 236 "zgeequ.f"
/* L40: */
#line 236 "zgeequ.f"
    }
#line 237 "zgeequ.f"
    *amax = rcmax;

#line 239 "zgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 243 "zgeequ.f"
	i__1 = *m;
#line 243 "zgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "zgeequ.f"
	    if (r__[i__] == 0.) {
#line 245 "zgeequ.f"
		*info = i__;
#line 246 "zgeequ.f"
		return 0;
#line 247 "zgeequ.f"
	    }
#line 248 "zgeequ.f"
/* L50: */
#line 248 "zgeequ.f"
	}
#line 249 "zgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 253 "zgeequ.f"
	i__1 = *m;
#line 253 "zgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 254 "zgeequ.f"
	    d__2 = r__[i__];
#line 254 "zgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 254 "zgeequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 255 "zgeequ.f"
/* L60: */
#line 255 "zgeequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 259 "zgeequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 260 "zgeequ.f"
    }

/*     Compute column scale factors */

#line 264 "zgeequ.f"
    i__1 = *n;
#line 264 "zgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 265 "zgeequ.f"
	c__[j] = 0.;
#line 266 "zgeequ.f"
/* L70: */
#line 266 "zgeequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 271 "zgeequ.f"
    i__1 = *n;
#line 271 "zgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 272 "zgeequ.f"
	i__2 = *m;
#line 272 "zgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 273 "zgeequ.f"
	    i__3 = i__ + j * a_dim1;
#line 273 "zgeequ.f"
	    d__3 = c__[j], d__4 = ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2))) * r__[i__];
#line 273 "zgeequ.f"
	    c__[j] = max(d__3,d__4);
#line 274 "zgeequ.f"
/* L80: */
#line 274 "zgeequ.f"
	}
#line 275 "zgeequ.f"
/* L90: */
#line 275 "zgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 279 "zgeequ.f"
    rcmin = bignum;
#line 280 "zgeequ.f"
    rcmax = 0.;
#line 281 "zgeequ.f"
    i__1 = *n;
#line 281 "zgeequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 282 "zgeequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 282 "zgeequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 283 "zgeequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 283 "zgeequ.f"
	rcmax = max(d__1,d__2);
#line 284 "zgeequ.f"
/* L100: */
#line 284 "zgeequ.f"
    }

#line 286 "zgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 290 "zgeequ.f"
	i__1 = *n;
#line 290 "zgeequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 291 "zgeequ.f"
	    if (c__[j] == 0.) {
#line 292 "zgeequ.f"
		*info = *m + j;
#line 293 "zgeequ.f"
		return 0;
#line 294 "zgeequ.f"
	    }
#line 295 "zgeequ.f"
/* L110: */
#line 295 "zgeequ.f"
	}
#line 296 "zgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 300 "zgeequ.f"
	i__1 = *n;
#line 300 "zgeequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 301 "zgeequ.f"
	    d__2 = c__[j];
#line 301 "zgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 301 "zgeequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 302 "zgeequ.f"
/* L120: */
#line 302 "zgeequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 306 "zgeequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 307 "zgeequ.f"
    }

#line 309 "zgeequ.f"
    return 0;

/*     End of ZGEEQU */

} /* zgeequ_ */

