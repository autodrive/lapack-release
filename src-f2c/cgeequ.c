#line 1 "cgeequ.f"
/* cgeequ.f -- translated by f2c (version 20100827).
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

#line 1 "cgeequ.f"
/* > \brief \b CGEEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          INFO ) */

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
/* > CGEEQU computes row and column scalings intended to equilibrate an */
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
/* Subroutine */ int cgeequ_(integer *m, integer *n, doublecomplex *a, 
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

/*     Test the input parameters. */

#line 188 "cgeequ.f"
    /* Parameter adjustments */
#line 188 "cgeequ.f"
    a_dim1 = *lda;
#line 188 "cgeequ.f"
    a_offset = 1 + a_dim1;
#line 188 "cgeequ.f"
    a -= a_offset;
#line 188 "cgeequ.f"
    --r__;
#line 188 "cgeequ.f"
    --c__;
#line 188 "cgeequ.f"

#line 188 "cgeequ.f"
    /* Function Body */
#line 188 "cgeequ.f"
    *info = 0;
#line 189 "cgeequ.f"
    if (*m < 0) {
#line 190 "cgeequ.f"
	*info = -1;
#line 191 "cgeequ.f"
    } else if (*n < 0) {
#line 192 "cgeequ.f"
	*info = -2;
#line 193 "cgeequ.f"
    } else if (*lda < max(1,*m)) {
#line 194 "cgeequ.f"
	*info = -4;
#line 195 "cgeequ.f"
    }
#line 196 "cgeequ.f"
    if (*info != 0) {
#line 197 "cgeequ.f"
	i__1 = -(*info);
#line 197 "cgeequ.f"
	xerbla_("CGEEQU", &i__1, (ftnlen)6);
#line 198 "cgeequ.f"
	return 0;
#line 199 "cgeequ.f"
    }

/*     Quick return if possible */

#line 203 "cgeequ.f"
    if (*m == 0 || *n == 0) {
#line 204 "cgeequ.f"
	*rowcnd = 1.;
#line 205 "cgeequ.f"
	*colcnd = 1.;
#line 206 "cgeequ.f"
	*amax = 0.;
#line 207 "cgeequ.f"
	return 0;
#line 208 "cgeequ.f"
    }

/*     Get machine constants. */

#line 212 "cgeequ.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 213 "cgeequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 217 "cgeequ.f"
    i__1 = *m;
#line 217 "cgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "cgeequ.f"
	r__[i__] = 0.;
#line 219 "cgeequ.f"
/* L10: */
#line 219 "cgeequ.f"
    }

/*     Find the maximum element in each row. */

#line 223 "cgeequ.f"
    i__1 = *n;
#line 223 "cgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 224 "cgeequ.f"
	i__2 = *m;
#line 224 "cgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 225 "cgeequ.f"
	    i__3 = i__ + j * a_dim1;
#line 225 "cgeequ.f"
	    d__3 = r__[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 225 "cgeequ.f"
	    r__[i__] = max(d__3,d__4);
#line 226 "cgeequ.f"
/* L20: */
#line 226 "cgeequ.f"
	}
#line 227 "cgeequ.f"
/* L30: */
#line 227 "cgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 231 "cgeequ.f"
    rcmin = bignum;
#line 232 "cgeequ.f"
    rcmax = 0.;
#line 233 "cgeequ.f"
    i__1 = *m;
#line 233 "cgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 234 "cgeequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 234 "cgeequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 235 "cgeequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 235 "cgeequ.f"
	rcmin = min(d__1,d__2);
#line 236 "cgeequ.f"
/* L40: */
#line 236 "cgeequ.f"
    }
#line 237 "cgeequ.f"
    *amax = rcmax;

#line 239 "cgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 243 "cgeequ.f"
	i__1 = *m;
#line 243 "cgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "cgeequ.f"
	    if (r__[i__] == 0.) {
#line 245 "cgeequ.f"
		*info = i__;
#line 246 "cgeequ.f"
		return 0;
#line 247 "cgeequ.f"
	    }
#line 248 "cgeequ.f"
/* L50: */
#line 248 "cgeequ.f"
	}
#line 249 "cgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 253 "cgeequ.f"
	i__1 = *m;
#line 253 "cgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 254 "cgeequ.f"
	    d__2 = r__[i__];
#line 254 "cgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 254 "cgeequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 255 "cgeequ.f"
/* L60: */
#line 255 "cgeequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 259 "cgeequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 260 "cgeequ.f"
    }

/*     Compute column scale factors */

#line 264 "cgeequ.f"
    i__1 = *n;
#line 264 "cgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 265 "cgeequ.f"
	c__[j] = 0.;
#line 266 "cgeequ.f"
/* L70: */
#line 266 "cgeequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 271 "cgeequ.f"
    i__1 = *n;
#line 271 "cgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 272 "cgeequ.f"
	i__2 = *m;
#line 272 "cgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 273 "cgeequ.f"
	    i__3 = i__ + j * a_dim1;
#line 273 "cgeequ.f"
	    d__3 = c__[j], d__4 = ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2))) * r__[i__];
#line 273 "cgeequ.f"
	    c__[j] = max(d__3,d__4);
#line 274 "cgeequ.f"
/* L80: */
#line 274 "cgeequ.f"
	}
#line 275 "cgeequ.f"
/* L90: */
#line 275 "cgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 279 "cgeequ.f"
    rcmin = bignum;
#line 280 "cgeequ.f"
    rcmax = 0.;
#line 281 "cgeequ.f"
    i__1 = *n;
#line 281 "cgeequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 282 "cgeequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 282 "cgeequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 283 "cgeequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 283 "cgeequ.f"
	rcmax = max(d__1,d__2);
#line 284 "cgeequ.f"
/* L100: */
#line 284 "cgeequ.f"
    }

#line 286 "cgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 290 "cgeequ.f"
	i__1 = *n;
#line 290 "cgeequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 291 "cgeequ.f"
	    if (c__[j] == 0.) {
#line 292 "cgeequ.f"
		*info = *m + j;
#line 293 "cgeequ.f"
		return 0;
#line 294 "cgeequ.f"
	    }
#line 295 "cgeequ.f"
/* L110: */
#line 295 "cgeequ.f"
	}
#line 296 "cgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 300 "cgeequ.f"
	i__1 = *n;
#line 300 "cgeequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 301 "cgeequ.f"
	    d__2 = c__[j];
#line 301 "cgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 301 "cgeequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 302 "cgeequ.f"
/* L120: */
#line 302 "cgeequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 306 "cgeequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 307 "cgeequ.f"
    }

#line 309 "cgeequ.f"
    return 0;

/*     End of CGEEQU */

} /* cgeequ_ */

