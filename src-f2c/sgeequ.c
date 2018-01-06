#line 1 "sgeequ.f"
/* sgeequ.f -- translated by f2c (version 20100827).
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

#line 1 "sgeequ.f"
/* > \brief \b SGEEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       REAL               AMAX, COLCND, ROWCND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( * ), R( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEEQU computes row and column scalings intended to equilibrate an */
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
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgeequ_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 179 "sgeequ.f"
    /* Parameter adjustments */
#line 179 "sgeequ.f"
    a_dim1 = *lda;
#line 179 "sgeequ.f"
    a_offset = 1 + a_dim1;
#line 179 "sgeequ.f"
    a -= a_offset;
#line 179 "sgeequ.f"
    --r__;
#line 179 "sgeequ.f"
    --c__;
#line 179 "sgeequ.f"

#line 179 "sgeequ.f"
    /* Function Body */
#line 179 "sgeequ.f"
    *info = 0;
#line 180 "sgeequ.f"
    if (*m < 0) {
#line 181 "sgeequ.f"
	*info = -1;
#line 182 "sgeequ.f"
    } else if (*n < 0) {
#line 183 "sgeequ.f"
	*info = -2;
#line 184 "sgeequ.f"
    } else if (*lda < max(1,*m)) {
#line 185 "sgeequ.f"
	*info = -4;
#line 186 "sgeequ.f"
    }
#line 187 "sgeequ.f"
    if (*info != 0) {
#line 188 "sgeequ.f"
	i__1 = -(*info);
#line 188 "sgeequ.f"
	xerbla_("SGEEQU", &i__1, (ftnlen)6);
#line 189 "sgeequ.f"
	return 0;
#line 190 "sgeequ.f"
    }

/*     Quick return if possible */

#line 194 "sgeequ.f"
    if (*m == 0 || *n == 0) {
#line 195 "sgeequ.f"
	*rowcnd = 1.;
#line 196 "sgeequ.f"
	*colcnd = 1.;
#line 197 "sgeequ.f"
	*amax = 0.;
#line 198 "sgeequ.f"
	return 0;
#line 199 "sgeequ.f"
    }

/*     Get machine constants. */

#line 203 "sgeequ.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 204 "sgeequ.f"
    bignum = 1. / smlnum;

/*     Compute row scale factors. */

#line 208 "sgeequ.f"
    i__1 = *m;
#line 208 "sgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "sgeequ.f"
	r__[i__] = 0.;
#line 210 "sgeequ.f"
/* L10: */
#line 210 "sgeequ.f"
    }

/*     Find the maximum element in each row. */

#line 214 "sgeequ.f"
    i__1 = *n;
#line 214 "sgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 215 "sgeequ.f"
	i__2 = *m;
#line 215 "sgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 216 "sgeequ.f"
	    d__2 = r__[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 216 "sgeequ.f"
	    r__[i__] = max(d__2,d__3);
#line 217 "sgeequ.f"
/* L20: */
#line 217 "sgeequ.f"
	}
#line 218 "sgeequ.f"
/* L30: */
#line 218 "sgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 222 "sgeequ.f"
    rcmin = bignum;
#line 223 "sgeequ.f"
    rcmax = 0.;
#line 224 "sgeequ.f"
    i__1 = *m;
#line 224 "sgeequ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 225 "sgeequ.f"
	d__1 = rcmax, d__2 = r__[i__];
#line 225 "sgeequ.f"
	rcmax = max(d__1,d__2);
/* Computing MIN */
#line 226 "sgeequ.f"
	d__1 = rcmin, d__2 = r__[i__];
#line 226 "sgeequ.f"
	rcmin = min(d__1,d__2);
#line 227 "sgeequ.f"
/* L40: */
#line 227 "sgeequ.f"
    }
#line 228 "sgeequ.f"
    *amax = rcmax;

#line 230 "sgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 234 "sgeequ.f"
	i__1 = *m;
#line 234 "sgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "sgeequ.f"
	    if (r__[i__] == 0.) {
#line 236 "sgeequ.f"
		*info = i__;
#line 237 "sgeequ.f"
		return 0;
#line 238 "sgeequ.f"
	    }
#line 239 "sgeequ.f"
/* L50: */
#line 239 "sgeequ.f"
	}
#line 240 "sgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 244 "sgeequ.f"
	i__1 = *m;
#line 244 "sgeequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
/* Computing MAX */
#line 245 "sgeequ.f"
	    d__2 = r__[i__];
#line 245 "sgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 245 "sgeequ.f"
	    r__[i__] = 1. / min(d__1,bignum);
#line 246 "sgeequ.f"
/* L60: */
#line 246 "sgeequ.f"
	}

/*        Compute ROWCND = min(R(I)) / max(R(I)) */

#line 250 "sgeequ.f"
	*rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 251 "sgeequ.f"
    }

/*     Compute column scale factors */

#line 255 "sgeequ.f"
    i__1 = *n;
#line 255 "sgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 256 "sgeequ.f"
	c__[j] = 0.;
#line 257 "sgeequ.f"
/* L70: */
#line 257 "sgeequ.f"
    }

/*     Find the maximum element in each column, */
/*     assuming the row scaling computed above. */

#line 262 "sgeequ.f"
    i__1 = *n;
#line 262 "sgeequ.f"
    for (j = 1; j <= i__1; ++j) {
#line 263 "sgeequ.f"
	i__2 = *m;
#line 263 "sgeequ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 264 "sgeequ.f"
	    d__2 = c__[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1)) * 
		    r__[i__];
#line 264 "sgeequ.f"
	    c__[j] = max(d__2,d__3);
#line 265 "sgeequ.f"
/* L80: */
#line 265 "sgeequ.f"
	}
#line 266 "sgeequ.f"
/* L90: */
#line 266 "sgeequ.f"
    }

/*     Find the maximum and minimum scale factors. */

#line 270 "sgeequ.f"
    rcmin = bignum;
#line 271 "sgeequ.f"
    rcmax = 0.;
#line 272 "sgeequ.f"
    i__1 = *n;
#line 272 "sgeequ.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 273 "sgeequ.f"
	d__1 = rcmin, d__2 = c__[j];
#line 273 "sgeequ.f"
	rcmin = min(d__1,d__2);
/* Computing MAX */
#line 274 "sgeequ.f"
	d__1 = rcmax, d__2 = c__[j];
#line 274 "sgeequ.f"
	rcmax = max(d__1,d__2);
#line 275 "sgeequ.f"
/* L100: */
#line 275 "sgeequ.f"
    }

#line 277 "sgeequ.f"
    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. */

#line 281 "sgeequ.f"
	i__1 = *n;
#line 281 "sgeequ.f"
	for (j = 1; j <= i__1; ++j) {
#line 282 "sgeequ.f"
	    if (c__[j] == 0.) {
#line 283 "sgeequ.f"
		*info = *m + j;
#line 284 "sgeequ.f"
		return 0;
#line 285 "sgeequ.f"
	    }
#line 286 "sgeequ.f"
/* L110: */
#line 286 "sgeequ.f"
	}
#line 287 "sgeequ.f"
    } else {

/*        Invert the scale factors. */

#line 291 "sgeequ.f"
	i__1 = *n;
#line 291 "sgeequ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
/* Computing MAX */
#line 292 "sgeequ.f"
	    d__2 = c__[j];
#line 292 "sgeequ.f"
	    d__1 = max(d__2,smlnum);
#line 292 "sgeequ.f"
	    c__[j] = 1. / min(d__1,bignum);
#line 293 "sgeequ.f"
/* L120: */
#line 293 "sgeequ.f"
	}

/*        Compute COLCND = min(C(J)) / max(C(J)) */

#line 297 "sgeequ.f"
	*colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 298 "sgeequ.f"
    }

#line 300 "sgeequ.f"
    return 0;

/*     End of SGEEQU */

} /* sgeequ_ */

