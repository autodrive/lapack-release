#line 1 "dsyequb.f"
/* dsyequb.f -- translated by f2c (version 20100827).
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

#line 1 "dsyequb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSYEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       CHARACTER          UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEQUB computes row and column scalings intended to equilibrate a */
/* > symmetric matrix A and reduce its condition number */
/* > (with respect to the two-norm).  S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The N-by-N symmetric matrix whose scaling */
/* >          factors are to be computed.  Only the diagonal elements of A */
/* >          are referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* >          SCOND is DOUBLE PRECISION */
/* >          If INFO = 0, S contains the ratio of the smallest S(i) to */
/* >          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the i-th diagonal element is nonpositive. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleSYcomputational */

/* > \par References: */
/*  ================ */
/* > */
/* >  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* >  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* >  DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* >  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsyequb_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *s, doublereal *scond, doublereal *amax, doublereal *
	work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal d__;
    static integer i__, j;
    static doublereal t, u, c0, c1, c2, si;
    static logical up;
    static doublereal avg, std, tol, base;
    static integer iter;
    static doublereal smin, smax, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal sumsq;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal smlnum;


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
/*     .. Executable Statements .. */

/*     Test input parameters. */

#line 181 "dsyequb.f"
    /* Parameter adjustments */
#line 181 "dsyequb.f"
    a_dim1 = *lda;
#line 181 "dsyequb.f"
    a_offset = 1 + a_dim1;
#line 181 "dsyequb.f"
    a -= a_offset;
#line 181 "dsyequb.f"
    --s;
#line 181 "dsyequb.f"
    --work;
#line 181 "dsyequb.f"

#line 181 "dsyequb.f"
    /* Function Body */
#line 181 "dsyequb.f"
    *info = 0;
#line 182 "dsyequb.f"
    if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1))) {
#line 183 "dsyequb.f"
	*info = -1;
#line 184 "dsyequb.f"
    } else if (*n < 0) {
#line 185 "dsyequb.f"
	*info = -2;
#line 186 "dsyequb.f"
    } else if (*lda < max(1,*n)) {
#line 187 "dsyequb.f"
	*info = -4;
#line 188 "dsyequb.f"
    }
#line 189 "dsyequb.f"
    if (*info != 0) {
#line 190 "dsyequb.f"
	i__1 = -(*info);
#line 190 "dsyequb.f"
	xerbla_("DSYEQUB", &i__1, (ftnlen)7);
#line 191 "dsyequb.f"
	return 0;
#line 192 "dsyequb.f"
    }
#line 194 "dsyequb.f"
    up = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 195 "dsyequb.f"
    *amax = 0.;

/*     Quick return if possible. */

#line 199 "dsyequb.f"
    if (*n == 0) {
#line 200 "dsyequb.f"
	*scond = 1.;
#line 201 "dsyequb.f"
	return 0;
#line 202 "dsyequb.f"
    }
#line 204 "dsyequb.f"
    i__1 = *n;
#line 204 "dsyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 205 "dsyequb.f"
	s[i__] = 0.;
#line 206 "dsyequb.f"
    }
#line 208 "dsyequb.f"
    *amax = 0.;
#line 209 "dsyequb.f"
    if (up) {
#line 210 "dsyequb.f"
	i__1 = *n;
#line 210 "dsyequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "dsyequb.f"
	    i__2 = j - 1;
#line 211 "dsyequb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 212 "dsyequb.f"
		d__2 = s[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 212 "dsyequb.f"
		s[i__] = max(d__2,d__3);
/* Computing MAX */
#line 213 "dsyequb.f"
		d__2 = s[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 213 "dsyequb.f"
		s[j] = max(d__2,d__3);
/* Computing MAX */
#line 214 "dsyequb.f"
		d__2 = *amax, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 214 "dsyequb.f"
		*amax = max(d__2,d__3);
#line 215 "dsyequb.f"
	    }
/* Computing MAX */
#line 216 "dsyequb.f"
	    d__2 = s[j], d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 216 "dsyequb.f"
	    s[j] = max(d__2,d__3);
/* Computing MAX */
#line 217 "dsyequb.f"
	    d__2 = *amax, d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 217 "dsyequb.f"
	    *amax = max(d__2,d__3);
#line 218 "dsyequb.f"
	}
#line 219 "dsyequb.f"
    } else {
#line 220 "dsyequb.f"
	i__1 = *n;
#line 220 "dsyequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 221 "dsyequb.f"
	    d__2 = s[j], d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 221 "dsyequb.f"
	    s[j] = max(d__2,d__3);
/* Computing MAX */
#line 222 "dsyequb.f"
	    d__2 = *amax, d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 222 "dsyequb.f"
	    *amax = max(d__2,d__3);
#line 223 "dsyequb.f"
	    i__2 = *n;
#line 223 "dsyequb.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 224 "dsyequb.f"
		d__2 = s[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 224 "dsyequb.f"
		s[i__] = max(d__2,d__3);
/* Computing MAX */
#line 225 "dsyequb.f"
		d__2 = s[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 225 "dsyequb.f"
		s[j] = max(d__2,d__3);
/* Computing MAX */
#line 226 "dsyequb.f"
		d__2 = *amax, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 226 "dsyequb.f"
		*amax = max(d__2,d__3);
#line 227 "dsyequb.f"
	    }
#line 228 "dsyequb.f"
	}
#line 229 "dsyequb.f"
    }
#line 230 "dsyequb.f"
    i__1 = *n;
#line 230 "dsyequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 231 "dsyequb.f"
	s[j] = 1. / s[j];
#line 232 "dsyequb.f"
    }
#line 234 "dsyequb.f"
    tol = 1. / sqrt(*n * 2.);
#line 236 "dsyequb.f"
    for (iter = 1; iter <= 100; ++iter) {
#line 237 "dsyequb.f"
	scale = 0.;
#line 238 "dsyequb.f"
	sumsq = 0.;
/*       BETA = |A|S */
#line 240 "dsyequb.f"
	i__1 = *n;
#line 240 "dsyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "dsyequb.f"
	    work[i__] = 0.;
#line 242 "dsyequb.f"
	}
#line 243 "dsyequb.f"
	if (up) {
#line 244 "dsyequb.f"
	    i__1 = *n;
#line 244 "dsyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "dsyequb.f"
		i__2 = j - 1;
#line 245 "dsyequb.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 246 "dsyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 247 "dsyequb.f"
		    work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    j];
#line 248 "dsyequb.f"
		    work[j] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    i__];
#line 249 "dsyequb.f"
		}
#line 250 "dsyequb.f"
		work[j] += (d__1 = a[j + j * a_dim1], abs(d__1)) * s[j];
#line 251 "dsyequb.f"
	    }
#line 252 "dsyequb.f"
	} else {
#line 253 "dsyequb.f"
	    i__1 = *n;
#line 253 "dsyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 254 "dsyequb.f"
		work[j] += (d__1 = a[j + j * a_dim1], abs(d__1)) * s[j];
#line 255 "dsyequb.f"
		i__2 = *n;
#line 255 "dsyequb.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 256 "dsyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 257 "dsyequb.f"
		    work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    j];
#line 258 "dsyequb.f"
		    work[j] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    i__];
#line 259 "dsyequb.f"
		}
#line 260 "dsyequb.f"
	    }
#line 261 "dsyequb.f"
	}
/*       avg = s^T beta / n */
#line 264 "dsyequb.f"
	avg = 0.;
#line 265 "dsyequb.f"
	i__1 = *n;
#line 265 "dsyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "dsyequb.f"
	    avg += s[i__] * work[i__];
#line 267 "dsyequb.f"
	}
#line 268 "dsyequb.f"
	avg /= *n;
#line 270 "dsyequb.f"
	std = 0.;
#line 271 "dsyequb.f"
	i__1 = *n * 3;
#line 271 "dsyequb.f"
	for (i__ = (*n << 1) + 1; i__ <= i__1; ++i__) {
#line 272 "dsyequb.f"
	    work[i__] = s[i__ - (*n << 1)] * work[i__ - (*n << 1)] - avg;
#line 273 "dsyequb.f"
	}
#line 274 "dsyequb.f"
	dlassq_(n, &work[(*n << 1) + 1], &c__1, &scale, &sumsq);
#line 275 "dsyequb.f"
	std = scale * sqrt(sumsq / *n);
#line 277 "dsyequb.f"
	if (std < tol * avg) {
#line 277 "dsyequb.f"
	    goto L999;
#line 277 "dsyequb.f"
	}
#line 279 "dsyequb.f"
	i__1 = *n;
#line 279 "dsyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "dsyequb.f"
	    t = (d__1 = a[i__ + i__ * a_dim1], abs(d__1));
#line 281 "dsyequb.f"
	    si = s[i__];
#line 282 "dsyequb.f"
	    c2 = (*n - 1) * t;
#line 283 "dsyequb.f"
	    c1 = (*n - 2) * (work[i__] - t * si);
#line 284 "dsyequb.f"
	    c0 = -(t * si) * si + work[i__] * 2 * si - *n * avg;
#line 285 "dsyequb.f"
	    d__ = c1 * c1 - c0 * 4 * c2;
#line 287 "dsyequb.f"
	    if (d__ <= 0.) {
#line 288 "dsyequb.f"
		*info = -1;
#line 289 "dsyequb.f"
		return 0;
#line 290 "dsyequb.f"
	    }
#line 291 "dsyequb.f"
	    si = c0 * -2 / (c1 + sqrt(d__));
#line 293 "dsyequb.f"
	    d__ = si - s[i__];
#line 294 "dsyequb.f"
	    u = 0.;
#line 295 "dsyequb.f"
	    if (up) {
#line 296 "dsyequb.f"
		i__2 = i__;
#line 296 "dsyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 297 "dsyequb.f"
		    t = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 298 "dsyequb.f"
		    u += s[j] * t;
#line 299 "dsyequb.f"
		    work[j] += d__ * t;
#line 300 "dsyequb.f"
		}
#line 301 "dsyequb.f"
		i__2 = *n;
#line 301 "dsyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 302 "dsyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 303 "dsyequb.f"
		    u += s[j] * t;
#line 304 "dsyequb.f"
		    work[j] += d__ * t;
#line 305 "dsyequb.f"
		}
#line 306 "dsyequb.f"
	    } else {
#line 307 "dsyequb.f"
		i__2 = i__;
#line 307 "dsyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 308 "dsyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 309 "dsyequb.f"
		    u += s[j] * t;
#line 310 "dsyequb.f"
		    work[j] += d__ * t;
#line 311 "dsyequb.f"
		}
#line 312 "dsyequb.f"
		i__2 = *n;
#line 312 "dsyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 313 "dsyequb.f"
		    t = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 314 "dsyequb.f"
		    u += s[j] * t;
#line 315 "dsyequb.f"
		    work[j] += d__ * t;
#line 316 "dsyequb.f"
		}
#line 317 "dsyequb.f"
	    }
#line 319 "dsyequb.f"
	    avg += (u + work[i__]) * d__ / *n;
#line 320 "dsyequb.f"
	    s[i__] = si;
#line 322 "dsyequb.f"
	}
#line 324 "dsyequb.f"
    }
#line 326 "dsyequb.f"
L999:
#line 328 "dsyequb.f"
    smlnum = dlamch_("SAFEMIN", (ftnlen)7);
#line 329 "dsyequb.f"
    bignum = 1. / smlnum;
#line 330 "dsyequb.f"
    smin = bignum;
#line 331 "dsyequb.f"
    smax = 0.;
#line 332 "dsyequb.f"
    t = 1. / sqrt(avg);
#line 333 "dsyequb.f"
    base = dlamch_("B", (ftnlen)1);
#line 334 "dsyequb.f"
    u = 1. / log(base);
#line 335 "dsyequb.f"
    i__1 = *n;
#line 335 "dsyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 336 "dsyequb.f"
	i__2 = (integer) (u * log(s[i__] * t));
#line 336 "dsyequb.f"
	s[i__] = pow_di(&base, &i__2);
/* Computing MIN */
#line 337 "dsyequb.f"
	d__1 = smin, d__2 = s[i__];
#line 337 "dsyequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 338 "dsyequb.f"
	d__1 = smax, d__2 = s[i__];
#line 338 "dsyequb.f"
	smax = max(d__1,d__2);
#line 339 "dsyequb.f"
    }
#line 340 "dsyequb.f"
    *scond = max(smin,smlnum) / min(smax,bignum);

#line 342 "dsyequb.f"
    return 0;
} /* dsyequb_ */

