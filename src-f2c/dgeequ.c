#line 1 "dgeequ.f"
/* dgeequ.f -- translated by f2c (version 20100827).
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

#line 1 "dgeequ.f"
/* > \brief \b DGEEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          INFO ) */

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
/* > DGEEQU computes row and column scalings intended to equilibrate an */
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
/* Subroutine */ int dgeequ_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 179 "dgeequ.f"
    /* Parameter adjustments */
#line 179 "dgeequ.f"
    a_dim1 = *lda;
#line 179 "dgeequ.f"
    a_offset = 1 + a_dim1;
#line 179 "dgeequ.f"
    a -= a_offset;
#line 179 "dgeequ.f"
    --r__;
#line 179 "dgeequ.f"
    --c__;
#line 179 "dgeequ.f"

#line 179 "dgeequ.f"
    /* Function Body */
#line 179 "dgeequ.f"
    *info = 0;
#line 180 "dgeequ.f"
    if (*m < 0) {
#line 181 "dgeequ.f"
	*info = -1;
#line 182 "dgeequ.f"
    } else if (*n < 0) {
#line 183 "dgeequ.f"
	*info = -2;
#line 184 "dgeequ.f"
    } else if (*lda < max(1,*m)) {
#line 185 "dgeequ.f"
	*info = -4;
#line 186 "dgeequ.f"
    }
#line 187 "dgeequ.f"
    if (*info != 0) {
#line 188 "dgeequ.f"
	i__1 = -(*info);
#line 188 "dgeequ.f"
	xerbla_("DGEEQU", &i__1, (ftnlen)6);
#line 189 "dgeequ.f"
	return 0;
#line 190 "dgeequ.f"
    }

/*     Quick return if possible */

#line 194 "dgeequ.f"
    if (*m == 0 || *n == 0) {
#line 195 "dgeequ.f"
	*rowcnd = 1.;
#line 196 "dgeequ.f"
	*colcnd = 1.;
#line 197 "dgeequ.f"
	*amax = 0.;
#line 198 "dgeequ.f"
	return 0;
#line 199 "dgeequ.f"
    }

/*     Get machine constants. */

#line 203 "dgeequ.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 204 "dgeequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 208 "dgeequ.f"
    i__1 = *m;
#line 208 "dgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "dgeequ.f"
	r__[i__] = 0.;
#line 210 "dgeequ.f"
/* L10: */
#line 210 "dgeequ.f"
    }

/*     Find the maximum element in each row. */

#line 214 "dgeequ.f"
    i__1 = *n;
#line 214 "dgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 215 "dgeequ.f"
	i__2 = *m;
#line 215 "dgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 216 "dgeequ.f"
	    d__2 = r__[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 216 "dgeequ.f"
	    r__[i__] = max(d__2,d__3);
#line 217 "dgeequ.f"
/* L20: */
#line 217 "dgeequ.f"
	}
#line 218 "dgeequ.f"
/* L30: */
#line 218 "dgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 222 "dgeequ.f"
    rcmin = bignum;
#line 223 "dgeequ.f"
    rcmax = 0.;
#line 224 "dgeequ.f"
    i__1 = *m;
#line 224 "dgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 225 "dgeequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 225 "dgeequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 226 "dgeequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 226 "dgeequ.f"
	rcmin = min(d__1,d__2);
#line 227 "dgeequ.f"
/* L40: */
#line 227 "dgeequ.f"
    }
#line 228 "dgeequ.f"
    *amax = rcmax;

#line 230 "dgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 234 "dgeequ.f"
	i__1 = *m;
#line 234 "dgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "dgeequ.f"
	    if (r__[i__] == 0.) {
#line 236 "dgeequ.f"
		*info = i__;
#line 237 "dgeequ.f"
		return 0;
#line 238 "dgeequ.f"
	    }
#line 239 "dgeequ.f"
/* L50: */
#line 239 "dgeequ.f"
	}
#line 240 "dgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 244 "dgeequ.f"
	i__1 = *m;
#line 244 "dgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 245 "dgeequ.f"
	    d__2 = r__[i__];
#line 245 "dgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 245 "dgeequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 246 "dgeequ.f"
/* L60: */
#line 246 "dgeequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 250 "dgeequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 251 "dgeequ.f"
    }

/*     Compute column scale factors */

#line 255 "dgeequ.f"
    i__1 = *n;
#line 255 "dgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 256 "dgeequ.f"
	c__[j] = 0.;
#line 257 "dgeequ.f"
/* L70: */
#line 257 "dgeequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 262 "dgeequ.f"
    i__1 = *n;
#line 262 "dgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 263 "dgeequ.f"
	i__2 = *m;
#line 263 "dgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 264 "dgeequ.f"
	    d__2 = c__[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1)) * 
		    r__[i__];
#line 264 "dgeequ.f"
	    c__[j] = max(d__2,d__3);
#line 265 "dgeequ.f"
/* L80: */
#line 265 "dgeequ.f"
	}
#line 266 "dgeequ.f"
/* L90: */
#line 266 "dgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 270 "dgeequ.f"
    rcmin = bignum;
#line 271 "dgeequ.f"
    rcmax = 0.;
#line 272 "dgeequ.f"
    i__1 = *n;
#line 272 "dgeequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 273 "dgeequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 273 "dgeequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 274 "dgeequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 274 "dgeequ.f"
	rcmax = max(d__1,d__2);
#line 275 "dgeequ.f"
/* L100: */
#line 275 "dgeequ.f"
    }

#line 277 "dgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 281 "dgeequ.f"
	i__1 = *n;
#line 281 "dgeequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 282 "dgeequ.f"
	    if (c__[j] == 0.) {
#line 283 "dgeequ.f"
		*info = *m + j;
#line 284 "dgeequ.f"
		return 0;
#line 285 "dgeequ.f"
	    }
#line 286 "dgeequ.f"
/* L110: */
#line 286 "dgeequ.f"
	}
#line 287 "dgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 291 "dgeequ.f"
	i__1 = *n;
#line 291 "dgeequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 292 "dgeequ.f"
	    d__2 = c__[j];
#line 292 "dgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 292 "dgeequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 293 "dgeequ.f"
/* L120: */
#line 293 "dgeequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 297 "dgeequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 298 "dgeequ.f"
    }

#line 300 "dgeequ.f"
    return 0;

/*     End of DGEEQU */

} /* dgeequ_ */

