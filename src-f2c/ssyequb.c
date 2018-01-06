#line 1 "ssyequb.f"
/* ssyequb.f -- translated by f2c (version 20100827).
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

#line 1 "ssyequb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SSYEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       REAL               AMAX, SCOND */
/*       CHARACTER          UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEQUB computes row and column scalings intended to equilibrate a */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          S is REAL array, dimension (N) */
/* >          If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* >          SCOND is REAL */
/* >          If INFO = 0, S contains the ratio of the smallest S(i) to */
/* >          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (3*N) */
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

/* > \ingroup realSYcomputational */

/* > \par References: */
/*  ================ */
/* > */
/* >  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* >  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* >  DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* >  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ssyequb_(char *uplo, integer *n, doublereal *a, integer *
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
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
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

#line 181 "ssyequb.f"
    /* Parameter adjustments */
#line 181 "ssyequb.f"
    a_dim1 = *lda;
#line 181 "ssyequb.f"
    a_offset = 1 + a_dim1;
#line 181 "ssyequb.f"
    a -= a_offset;
#line 181 "ssyequb.f"
    --s;
#line 181 "ssyequb.f"
    --work;
#line 181 "ssyequb.f"

#line 181 "ssyequb.f"
    /* Function Body */
#line 181 "ssyequb.f"
    *info = 0;
#line 182 "ssyequb.f"
    if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1))) {
#line 183 "ssyequb.f"
	*info = -1;
#line 184 "ssyequb.f"
    } else if (*n < 0) {
#line 185 "ssyequb.f"
	*info = -2;
#line 186 "ssyequb.f"
    } else if (*lda < max(1,*n)) {
#line 187 "ssyequb.f"
	*info = -4;
#line 188 "ssyequb.f"
    }
#line 189 "ssyequb.f"
    if (*info != 0) {
#line 190 "ssyequb.f"
	i__1 = -(*info);
#line 190 "ssyequb.f"
	xerbla_("SSYEQUB", &i__1, (ftnlen)7);
#line 191 "ssyequb.f"
	return 0;
#line 192 "ssyequb.f"
    }
#line 194 "ssyequb.f"
    up = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 195 "ssyequb.f"
    *amax = 0.;

/*     Quick return if possible. */

#line 199 "ssyequb.f"
    if (*n == 0) {
#line 200 "ssyequb.f"
	*scond = 1.;
#line 201 "ssyequb.f"
	return 0;
#line 202 "ssyequb.f"
    }
#line 204 "ssyequb.f"
    i__1 = *n;
#line 204 "ssyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 205 "ssyequb.f"
	s[i__] = 0.;
#line 206 "ssyequb.f"
    }
#line 208 "ssyequb.f"
    *amax = 0.;
#line 209 "ssyequb.f"
    if (up) {
#line 210 "ssyequb.f"
	i__1 = *n;
#line 210 "ssyequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "ssyequb.f"
	    i__2 = j - 1;
#line 211 "ssyequb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 212 "ssyequb.f"
		d__2 = s[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 212 "ssyequb.f"
		s[i__] = max(d__2,d__3);
/* Computing MAX */
#line 213 "ssyequb.f"
		d__2 = s[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 213 "ssyequb.f"
		s[j] = max(d__2,d__3);
/* Computing MAX */
#line 214 "ssyequb.f"
		d__2 = *amax, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 214 "ssyequb.f"
		*amax = max(d__2,d__3);
#line 215 "ssyequb.f"
	    }
/* Computing MAX */
#line 216 "ssyequb.f"
	    d__2 = s[j], d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 216 "ssyequb.f"
	    s[j] = max(d__2,d__3);
/* Computing MAX */
#line 217 "ssyequb.f"
	    d__2 = *amax, d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 217 "ssyequb.f"
	    *amax = max(d__2,d__3);
#line 218 "ssyequb.f"
	}
#line 219 "ssyequb.f"
    } else {
#line 220 "ssyequb.f"
	i__1 = *n;
#line 220 "ssyequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 221 "ssyequb.f"
	    d__2 = s[j], d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 221 "ssyequb.f"
	    s[j] = max(d__2,d__3);
/* Computing MAX */
#line 222 "ssyequb.f"
	    d__2 = *amax, d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 222 "ssyequb.f"
	    *amax = max(d__2,d__3);
#line 223 "ssyequb.f"
	    i__2 = *n;
#line 223 "ssyequb.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 224 "ssyequb.f"
		d__2 = s[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 224 "ssyequb.f"
		s[i__] = max(d__2,d__3);
/* Computing MAX */
#line 225 "ssyequb.f"
		d__2 = s[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 225 "ssyequb.f"
		s[j] = max(d__2,d__3);
/* Computing MAX */
#line 226 "ssyequb.f"
		d__2 = *amax, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 226 "ssyequb.f"
		*amax = max(d__2,d__3);
#line 227 "ssyequb.f"
	    }
#line 228 "ssyequb.f"
	}
#line 229 "ssyequb.f"
    }
#line 230 "ssyequb.f"
    i__1 = *n;
#line 230 "ssyequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 231 "ssyequb.f"
	s[j] = 1. / s[j];
#line 232 "ssyequb.f"
    }
#line 234 "ssyequb.f"
    tol = 1. / sqrt(*n * 2.);
#line 236 "ssyequb.f"
    for (iter = 1; iter <= 100; ++iter) {
#line 237 "ssyequb.f"
	scale = 0.;
#line 238 "ssyequb.f"
	sumsq = 0.;
/*       BETA = |A|S */
#line 240 "ssyequb.f"
	i__1 = *n;
#line 240 "ssyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "ssyequb.f"
	    work[i__] = 0.;
#line 242 "ssyequb.f"
	}
#line 243 "ssyequb.f"
	if (up) {
#line 244 "ssyequb.f"
	    i__1 = *n;
#line 244 "ssyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "ssyequb.f"
		i__2 = j - 1;
#line 245 "ssyequb.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 246 "ssyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 247 "ssyequb.f"
		    work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    j];
#line 248 "ssyequb.f"
		    work[j] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    i__];
#line 249 "ssyequb.f"
		}
#line 250 "ssyequb.f"
		work[j] += (d__1 = a[j + j * a_dim1], abs(d__1)) * s[j];
#line 251 "ssyequb.f"
	    }
#line 252 "ssyequb.f"
	} else {
#line 253 "ssyequb.f"
	    i__1 = *n;
#line 253 "ssyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 254 "ssyequb.f"
		work[j] += (d__1 = a[j + j * a_dim1], abs(d__1)) * s[j];
#line 255 "ssyequb.f"
		i__2 = *n;
#line 255 "ssyequb.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 256 "ssyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 257 "ssyequb.f"
		    work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    j];
#line 258 "ssyequb.f"
		    work[j] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    i__];
#line 259 "ssyequb.f"
		}
#line 260 "ssyequb.f"
	    }
#line 261 "ssyequb.f"
	}
/*       avg = s^T beta / n */
#line 264 "ssyequb.f"
	avg = 0.;
#line 265 "ssyequb.f"
	i__1 = *n;
#line 265 "ssyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 266 "ssyequb.f"
	    avg += s[i__] * work[i__];
#line 267 "ssyequb.f"
	}
#line 268 "ssyequb.f"
	avg /= *n;
#line 270 "ssyequb.f"
	std = 0.;
#line 271 "ssyequb.f"
	i__1 = *n * 3;
#line 271 "ssyequb.f"
	for (i__ = (*n << 1) + 1; i__ <= i__1; ++i__) {
#line 272 "ssyequb.f"
	    work[i__] = s[i__ - (*n << 1)] * work[i__ - (*n << 1)] - avg;
#line 273 "ssyequb.f"
	}
#line 274 "ssyequb.f"
	slassq_(n, &work[(*n << 1) + 1], &c__1, &scale, &sumsq);
#line 275 "ssyequb.f"
	std = scale * sqrt(sumsq / *n);
#line 277 "ssyequb.f"
	if (std < tol * avg) {
#line 277 "ssyequb.f"
	    goto L999;
#line 277 "ssyequb.f"
	}
#line 279 "ssyequb.f"
	i__1 = *n;
#line 279 "ssyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "ssyequb.f"
	    t = (d__1 = a[i__ + i__ * a_dim1], abs(d__1));
#line 281 "ssyequb.f"
	    si = s[i__];
#line 282 "ssyequb.f"
	    c2 = (*n - 1) * t;
#line 283 "ssyequb.f"
	    c1 = (*n - 2) * (work[i__] - t * si);
#line 284 "ssyequb.f"
	    c0 = -(t * si) * si + work[i__] * 2 * si - *n * avg;
#line 285 "ssyequb.f"
	    d__ = c1 * c1 - c0 * 4 * c2;
#line 287 "ssyequb.f"
	    if (d__ <= 0.) {
#line 288 "ssyequb.f"
		*info = -1;
#line 289 "ssyequb.f"
		return 0;
#line 290 "ssyequb.f"
	    }
#line 291 "ssyequb.f"
	    si = c0 * -2 / (c1 + sqrt(d__));
#line 293 "ssyequb.f"
	    d__ = si - s[i__];
#line 294 "ssyequb.f"
	    u = 0.;
#line 295 "ssyequb.f"
	    if (up) {
#line 296 "ssyequb.f"
		i__2 = i__;
#line 296 "ssyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 297 "ssyequb.f"
		    t = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 298 "ssyequb.f"
		    u += s[j] * t;
#line 299 "ssyequb.f"
		    work[j] += d__ * t;
#line 300 "ssyequb.f"
		}
#line 301 "ssyequb.f"
		i__2 = *n;
#line 301 "ssyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 302 "ssyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 303 "ssyequb.f"
		    u += s[j] * t;
#line 304 "ssyequb.f"
		    work[j] += d__ * t;
#line 305 "ssyequb.f"
		}
#line 306 "ssyequb.f"
	    } else {
#line 307 "ssyequb.f"
		i__2 = i__;
#line 307 "ssyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 308 "ssyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 309 "ssyequb.f"
		    u += s[j] * t;
#line 310 "ssyequb.f"
		    work[j] += d__ * t;
#line 311 "ssyequb.f"
		}
#line 312 "ssyequb.f"
		i__2 = *n;
#line 312 "ssyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 313 "ssyequb.f"
		    t = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 314 "ssyequb.f"
		    u += s[j] * t;
#line 315 "ssyequb.f"
		    work[j] += d__ * t;
#line 316 "ssyequb.f"
		}
#line 317 "ssyequb.f"
	    }
#line 319 "ssyequb.f"
	    avg += (u + work[i__]) * d__ / *n;
#line 320 "ssyequb.f"
	    s[i__] = si;
#line 322 "ssyequb.f"
	}
#line 324 "ssyequb.f"
    }
#line 326 "ssyequb.f"
L999:
#line 328 "ssyequb.f"
    smlnum = slamch_("SAFEMIN", (ftnlen)7);
#line 329 "ssyequb.f"
    bignum = 1. / smlnum;
#line 330 "ssyequb.f"
    smin = bignum;
#line 331 "ssyequb.f"
    smax = 0.;
#line 332 "ssyequb.f"
    t = 1. / sqrt(avg);
#line 333 "ssyequb.f"
    base = slamch_("B", (ftnlen)1);
#line 334 "ssyequb.f"
    u = 1. / log(base);
#line 335 "ssyequb.f"
    i__1 = *n;
#line 335 "ssyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 336 "ssyequb.f"
	i__2 = (integer) (u * log(s[i__] * t));
#line 336 "ssyequb.f"
	s[i__] = pow_di(&base, &i__2);
/* Computing MIN */
#line 337 "ssyequb.f"
	d__1 = smin, d__2 = s[i__];
#line 337 "ssyequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 338 "ssyequb.f"
	d__1 = smax, d__2 = s[i__];
#line 338 "ssyequb.f"
	smax = max(d__1,d__2);
#line 339 "ssyequb.f"
    }
#line 340 "ssyequb.f"
    *scond = max(smin,smlnum) / min(smax,bignum);

#line 342 "ssyequb.f"
    return 0;
} /* ssyequb_ */

