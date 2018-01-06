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
/* > symmetric matrix A (with respect to the Euclidean norm) and reduce */
/* > its condition number. The scale factors S are computed by the BIN */
/* > algorithm (see references) so that the scaled matrix B with elements */
/* > B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of */
/* > the smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The N-by-N symmetric matrix whose scaling factors are to be */
/* >          computed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
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
/* >          the largest S(i). If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is DOUBLE PRECISION */
/* >          Largest absolute value of any matrix element. If AMAX is */
/* >          very close to overflow or very close to underflow, the */
/* >          matrix should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \date November 2017 */

/* > \ingroup doubleSYcomputational */

/* > \par References: */
/*  ================ */
/* > */
/* >  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* >  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* >  DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* >  Tech report version: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.1679 */
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


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 177 "dsyequb.f"
    /* Parameter adjustments */
#line 177 "dsyequb.f"
    a_dim1 = *lda;
#line 177 "dsyequb.f"
    a_offset = 1 + a_dim1;
#line 177 "dsyequb.f"
    a -= a_offset;
#line 177 "dsyequb.f"
    --s;
#line 177 "dsyequb.f"
    --work;
#line 177 "dsyequb.f"

#line 177 "dsyequb.f"
    /* Function Body */
#line 177 "dsyequb.f"
    *info = 0;
#line 178 "dsyequb.f"
    if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1))) {
#line 179 "dsyequb.f"
	*info = -1;
#line 180 "dsyequb.f"
    } else if (*n < 0) {
#line 181 "dsyequb.f"
	*info = -2;
#line 182 "dsyequb.f"
    } else if (*lda < max(1,*n)) {
#line 183 "dsyequb.f"
	*info = -4;
#line 184 "dsyequb.f"
    }
#line 185 "dsyequb.f"
    if (*info != 0) {
#line 186 "dsyequb.f"
	i__1 = -(*info);
#line 186 "dsyequb.f"
	xerbla_("DSYEQUB", &i__1, (ftnlen)7);
#line 187 "dsyequb.f"
	return 0;
#line 188 "dsyequb.f"
    }
#line 190 "dsyequb.f"
    up = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 191 "dsyequb.f"
    *amax = 0.;

/*     Quick return if possible. */

#line 195 "dsyequb.f"
    if (*n == 0) {
#line 196 "dsyequb.f"
	*scond = 1.;
#line 197 "dsyequb.f"
	return 0;
#line 198 "dsyequb.f"
    }
#line 200 "dsyequb.f"
    i__1 = *n;
#line 200 "dsyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 201 "dsyequb.f"
	s[i__] = 0.;
#line 202 "dsyequb.f"
    }
#line 204 "dsyequb.f"
    *amax = 0.;
#line 205 "dsyequb.f"
    if (up) {
#line 206 "dsyequb.f"
	i__1 = *n;
#line 206 "dsyequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 207 "dsyequb.f"
	    i__2 = j - 1;
#line 207 "dsyequb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 208 "dsyequb.f"
		d__2 = s[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 208 "dsyequb.f"
		s[i__] = max(d__2,d__3);
/* Computing MAX */
#line 209 "dsyequb.f"
		d__2 = s[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 209 "dsyequb.f"
		s[j] = max(d__2,d__3);
/* Computing MAX */
#line 210 "dsyequb.f"
		d__2 = *amax, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 210 "dsyequb.f"
		*amax = max(d__2,d__3);
#line 211 "dsyequb.f"
	    }
/* Computing MAX */
#line 212 "dsyequb.f"
	    d__2 = s[j], d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 212 "dsyequb.f"
	    s[j] = max(d__2,d__3);
/* Computing MAX */
#line 213 "dsyequb.f"
	    d__2 = *amax, d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 213 "dsyequb.f"
	    *amax = max(d__2,d__3);
#line 214 "dsyequb.f"
	}
#line 215 "dsyequb.f"
    } else {
#line 216 "dsyequb.f"
	i__1 = *n;
#line 216 "dsyequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 217 "dsyequb.f"
	    d__2 = s[j], d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 217 "dsyequb.f"
	    s[j] = max(d__2,d__3);
/* Computing MAX */
#line 218 "dsyequb.f"
	    d__2 = *amax, d__3 = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 218 "dsyequb.f"
	    *amax = max(d__2,d__3);
#line 219 "dsyequb.f"
	    i__2 = *n;
#line 219 "dsyequb.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 220 "dsyequb.f"
		d__2 = s[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 220 "dsyequb.f"
		s[i__] = max(d__2,d__3);
/* Computing MAX */
#line 221 "dsyequb.f"
		d__2 = s[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 221 "dsyequb.f"
		s[j] = max(d__2,d__3);
/* Computing MAX */
#line 222 "dsyequb.f"
		d__2 = *amax, d__3 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 222 "dsyequb.f"
		*amax = max(d__2,d__3);
#line 223 "dsyequb.f"
	    }
#line 224 "dsyequb.f"
	}
#line 225 "dsyequb.f"
    }
#line 226 "dsyequb.f"
    i__1 = *n;
#line 226 "dsyequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 227 "dsyequb.f"
	s[j] = 1. / s[j];
#line 228 "dsyequb.f"
    }
#line 230 "dsyequb.f"
    tol = 1. / sqrt(*n * 2.);
#line 232 "dsyequb.f"
    for (iter = 1; iter <= 100; ++iter) {
#line 233 "dsyequb.f"
	scale = 0.;
#line 234 "dsyequb.f"
	sumsq = 0.;
/*        beta = |A|s */
#line 236 "dsyequb.f"
	i__1 = *n;
#line 236 "dsyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 237 "dsyequb.f"
	    work[i__] = 0.;
#line 238 "dsyequb.f"
	}
#line 239 "dsyequb.f"
	if (up) {
#line 240 "dsyequb.f"
	    i__1 = *n;
#line 240 "dsyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 241 "dsyequb.f"
		i__2 = j - 1;
#line 241 "dsyequb.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 242 "dsyequb.f"
		    work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    j];
#line 243 "dsyequb.f"
		    work[j] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    i__];
#line 244 "dsyequb.f"
		}
#line 245 "dsyequb.f"
		work[j] += (d__1 = a[j + j * a_dim1], abs(d__1)) * s[j];
#line 246 "dsyequb.f"
	    }
#line 247 "dsyequb.f"
	} else {
#line 248 "dsyequb.f"
	    i__1 = *n;
#line 248 "dsyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 249 "dsyequb.f"
		work[j] += (d__1 = a[j + j * a_dim1], abs(d__1)) * s[j];
#line 250 "dsyequb.f"
		i__2 = *n;
#line 250 "dsyequb.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 251 "dsyequb.f"
		    work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    j];
#line 252 "dsyequb.f"
		    work[j] += (d__1 = a[i__ + j * a_dim1], abs(d__1)) * s[
			    i__];
#line 253 "dsyequb.f"
		}
#line 254 "dsyequb.f"
	    }
#line 255 "dsyequb.f"
	}
/*        avg = s^T beta / n */
#line 258 "dsyequb.f"
	avg = 0.;
#line 259 "dsyequb.f"
	i__1 = *n;
#line 259 "dsyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 260 "dsyequb.f"
	    avg += s[i__] * work[i__];
#line 261 "dsyequb.f"
	}
#line 262 "dsyequb.f"
	avg /= *n;
#line 264 "dsyequb.f"
	std = 0.;
#line 265 "dsyequb.f"
	i__1 = *n << 1;
#line 265 "dsyequb.f"
	for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 266 "dsyequb.f"
	    work[i__] = s[i__ - *n] * work[i__ - *n] - avg;
#line 267 "dsyequb.f"
	}
#line 268 "dsyequb.f"
	dlassq_(n, &work[*n + 1], &c__1, &scale, &sumsq);
#line 269 "dsyequb.f"
	std = scale * sqrt(sumsq / *n);
#line 271 "dsyequb.f"
	if (std < tol * avg) {
#line 271 "dsyequb.f"
	    goto L999;
#line 271 "dsyequb.f"
	}
#line 273 "dsyequb.f"
	i__1 = *n;
#line 273 "dsyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "dsyequb.f"
	    t = (d__1 = a[i__ + i__ * a_dim1], abs(d__1));
#line 275 "dsyequb.f"
	    si = s[i__];
#line 276 "dsyequb.f"
	    c2 = (*n - 1) * t;
#line 277 "dsyequb.f"
	    c1 = (*n - 2) * (work[i__] - t * si);
#line 278 "dsyequb.f"
	    c0 = -(t * si) * si + work[i__] * 2 * si - *n * avg;
#line 279 "dsyequb.f"
	    d__ = c1 * c1 - c0 * 4 * c2;
#line 281 "dsyequb.f"
	    if (d__ <= 0.) {
#line 282 "dsyequb.f"
		*info = -1;
#line 283 "dsyequb.f"
		return 0;
#line 284 "dsyequb.f"
	    }
#line 285 "dsyequb.f"
	    si = c0 * -2 / (c1 + sqrt(d__));
#line 287 "dsyequb.f"
	    d__ = si - s[i__];
#line 288 "dsyequb.f"
	    u = 0.;
#line 289 "dsyequb.f"
	    if (up) {
#line 290 "dsyequb.f"
		i__2 = i__;
#line 290 "dsyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 291 "dsyequb.f"
		    t = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 292 "dsyequb.f"
		    u += s[j] * t;
#line 293 "dsyequb.f"
		    work[j] += d__ * t;
#line 294 "dsyequb.f"
		}
#line 295 "dsyequb.f"
		i__2 = *n;
#line 295 "dsyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 296 "dsyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 297 "dsyequb.f"
		    u += s[j] * t;
#line 298 "dsyequb.f"
		    work[j] += d__ * t;
#line 299 "dsyequb.f"
		}
#line 300 "dsyequb.f"
	    } else {
#line 301 "dsyequb.f"
		i__2 = i__;
#line 301 "dsyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 302 "dsyequb.f"
		    t = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 303 "dsyequb.f"
		    u += s[j] * t;
#line 304 "dsyequb.f"
		    work[j] += d__ * t;
#line 305 "dsyequb.f"
		}
#line 306 "dsyequb.f"
		i__2 = *n;
#line 306 "dsyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 307 "dsyequb.f"
		    t = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 308 "dsyequb.f"
		    u += s[j] * t;
#line 309 "dsyequb.f"
		    work[j] += d__ * t;
#line 310 "dsyequb.f"
		}
#line 311 "dsyequb.f"
	    }
#line 313 "dsyequb.f"
	    avg += (u + work[i__]) * d__ / *n;
#line 314 "dsyequb.f"
	    s[i__] = si;
#line 315 "dsyequb.f"
	}
#line 316 "dsyequb.f"
    }
#line 318 "dsyequb.f"
L999:
#line 320 "dsyequb.f"
    smlnum = dlamch_("SAFEMIN", (ftnlen)7);
#line 321 "dsyequb.f"
    bignum = 1. / smlnum;
#line 322 "dsyequb.f"
    smin = bignum;
#line 323 "dsyequb.f"
    smax = 0.;
#line 324 "dsyequb.f"
    t = 1. / sqrt(avg);
#line 325 "dsyequb.f"
    base = dlamch_("B", (ftnlen)1);
#line 326 "dsyequb.f"
    u = 1. / log(base);
#line 327 "dsyequb.f"
    i__1 = *n;
#line 327 "dsyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 328 "dsyequb.f"
	i__2 = (integer) (u * log(s[i__] * t));
#line 328 "dsyequb.f"
	s[i__] = pow_di(&base, &i__2);
/* Computing MIN */
#line 329 "dsyequb.f"
	d__1 = smin, d__2 = s[i__];
#line 329 "dsyequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 330 "dsyequb.f"
	d__1 = smax, d__2 = s[i__];
#line 330 "dsyequb.f"
	smax = max(d__1,d__2);
#line 331 "dsyequb.f"
    }
#line 332 "dsyequb.f"
    *scond = max(smin,smlnum) / min(smax,bignum);

#line 334 "dsyequb.f"
    return 0;
} /* dsyequb_ */

