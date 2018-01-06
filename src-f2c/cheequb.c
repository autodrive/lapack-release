#line 1 "cheequb.f"
/* cheequb.f -- translated by f2c (version 20100827).
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

#line 1 "cheequb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHEEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       REAL               AMAX, SCOND */
/*       CHARACTER          UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       REAL               S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEQUB computes row and column scalings intended to equilibrate a */
/* > Hermitian matrix A (with respect to the Euclidean norm) and reduce */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The N-by-N Hermitian matrix whose scaling factors are to be */
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
/* >          S is REAL array, dimension (N) */
/* >          If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* >          SCOND is REAL */
/* >          If INFO = 0, S contains the ratio of the smallest S(i) to */
/* >          the largest S(i). If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is REAL */
/* >          Largest absolute value of any matrix element. If AMAX is */
/* >          very close to overflow or very close to underflow, the */
/* >          matrix should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N) */
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

/* > \date April 2012 */

/* > \ingroup complexHEcomputational */

/* > \par References: */
/*  ================ */
/* > */
/* >  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* >  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* >  DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* >  Tech report version: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.1679 */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cheequb_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal), log(doublereal), pow_di(
	    doublereal *, integer *);

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
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 186 "cheequb.f"
    /* Parameter adjustments */
#line 186 "cheequb.f"
    a_dim1 = *lda;
#line 186 "cheequb.f"
    a_offset = 1 + a_dim1;
#line 186 "cheequb.f"
    a -= a_offset;
#line 186 "cheequb.f"
    --s;
#line 186 "cheequb.f"
    --work;
#line 186 "cheequb.f"

#line 186 "cheequb.f"
    /* Function Body */
#line 186 "cheequb.f"
    *info = 0;
#line 187 "cheequb.f"
    if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1))) {
#line 188 "cheequb.f"
	*info = -1;
#line 189 "cheequb.f"
    } else if (*n < 0) {
#line 190 "cheequb.f"
	*info = -2;
#line 191 "cheequb.f"
    } else if (*lda < max(1,*n)) {
#line 192 "cheequb.f"
	*info = -4;
#line 193 "cheequb.f"
    }
#line 194 "cheequb.f"
    if (*info != 0) {
#line 195 "cheequb.f"
	i__1 = -(*info);
#line 195 "cheequb.f"
	xerbla_("CHEEQUB", &i__1, (ftnlen)7);
#line 196 "cheequb.f"
	return 0;
#line 197 "cheequb.f"
    }
#line 199 "cheequb.f"
    up = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 200 "cheequb.f"
    *amax = 0.;

/*     Quick return if possible. */

#line 204 "cheequb.f"
    if (*n == 0) {
#line 205 "cheequb.f"
	*scond = 1.;
#line 206 "cheequb.f"
	return 0;
#line 207 "cheequb.f"
    }
#line 209 "cheequb.f"
    i__1 = *n;
#line 209 "cheequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 210 "cheequb.f"
	s[i__] = 0.;
#line 211 "cheequb.f"
    }
#line 213 "cheequb.f"
    *amax = 0.;
#line 214 "cheequb.f"
    if (up) {
#line 215 "cheequb.f"
	i__1 = *n;
#line 215 "cheequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 216 "cheequb.f"
	    i__2 = j - 1;
#line 216 "cheequb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 217 "cheequb.f"
		i__3 = i__ + j * a_dim1;
#line 217 "cheequb.f"
		d__3 = s[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 217 "cheequb.f"
		s[i__] = max(d__3,d__4);
/* Computing MAX */
#line 218 "cheequb.f"
		i__3 = i__ + j * a_dim1;
#line 218 "cheequb.f"
		d__3 = s[j], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 218 "cheequb.f"
		s[j] = max(d__3,d__4);
/* Computing MAX */
#line 219 "cheequb.f"
		i__3 = i__ + j * a_dim1;
#line 219 "cheequb.f"
		d__3 = *amax, d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 219 "cheequb.f"
		*amax = max(d__3,d__4);
#line 220 "cheequb.f"
	    }
/* Computing MAX */
#line 221 "cheequb.f"
	    i__2 = j + j * a_dim1;
#line 221 "cheequb.f"
	    d__3 = s[j], d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 221 "cheequb.f"
	    s[j] = max(d__3,d__4);
/* Computing MAX */
#line 222 "cheequb.f"
	    i__2 = j + j * a_dim1;
#line 222 "cheequb.f"
	    d__3 = *amax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 222 "cheequb.f"
	    *amax = max(d__3,d__4);
#line 223 "cheequb.f"
	}
#line 224 "cheequb.f"
    } else {
#line 225 "cheequb.f"
	i__1 = *n;
#line 225 "cheequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 226 "cheequb.f"
	    i__2 = j + j * a_dim1;
#line 226 "cheequb.f"
	    d__3 = s[j], d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 226 "cheequb.f"
	    s[j] = max(d__3,d__4);
/* Computing MAX */
#line 227 "cheequb.f"
	    i__2 = j + j * a_dim1;
#line 227 "cheequb.f"
	    d__3 = *amax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 227 "cheequb.f"
	    *amax = max(d__3,d__4);
#line 228 "cheequb.f"
	    i__2 = *n;
#line 228 "cheequb.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 229 "cheequb.f"
		i__3 = i__ + j * a_dim1;
#line 229 "cheequb.f"
		d__3 = s[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 229 "cheequb.f"
		s[i__] = max(d__3,d__4);
/* Computing MAX */
#line 230 "cheequb.f"
		i__3 = i__ + j * a_dim1;
#line 230 "cheequb.f"
		d__3 = s[j], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 230 "cheequb.f"
		s[j] = max(d__3,d__4);
/* Computing MAX */
#line 231 "cheequb.f"
		i__3 = i__ + j * a_dim1;
#line 231 "cheequb.f"
		d__3 = *amax, d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 231 "cheequb.f"
		*amax = max(d__3,d__4);
#line 232 "cheequb.f"
	    }
#line 233 "cheequb.f"
	}
#line 234 "cheequb.f"
    }
#line 235 "cheequb.f"
    i__1 = *n;
#line 235 "cheequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 236 "cheequb.f"
	s[j] = 1. / s[j];
#line 237 "cheequb.f"
    }
#line 239 "cheequb.f"
    tol = 1. / sqrt(*n * 2.);
#line 241 "cheequb.f"
    for (iter = 1; iter <= 100; ++iter) {
#line 242 "cheequb.f"
	scale = 0.;
#line 243 "cheequb.f"
	sumsq = 0.;
/*        beta = |A|s */
#line 245 "cheequb.f"
	i__1 = *n;
#line 245 "cheequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "cheequb.f"
	    i__2 = i__;
#line 246 "cheequb.f"
	    work[i__2].r = 0., work[i__2].i = 0.;
#line 247 "cheequb.f"
	}
#line 248 "cheequb.f"
	if (up) {
#line 249 "cheequb.f"
	    i__1 = *n;
#line 249 "cheequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 250 "cheequb.f"
		i__2 = j - 1;
#line 250 "cheequb.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 251 "cheequb.f"
		    i__3 = i__;
#line 251 "cheequb.f"
		    i__4 = i__;
#line 251 "cheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 251 "cheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[j];
#line 251 "cheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 251 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 252 "cheequb.f"
		    i__3 = j;
#line 252 "cheequb.f"
		    i__4 = j;
#line 252 "cheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 252 "cheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[i__];
#line 252 "cheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 252 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 253 "cheequb.f"
		}
#line 254 "cheequb.f"
		i__2 = j;
#line 254 "cheequb.f"
		i__3 = j;
#line 254 "cheequb.f"
		i__4 = j + j * a_dim1;
#line 254 "cheequb.f"
		d__3 = ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			j * a_dim1]), abs(d__2))) * s[j];
#line 254 "cheequb.f"
		z__1.r = work[i__3].r + d__3, z__1.i = work[i__3].i;
#line 254 "cheequb.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 255 "cheequb.f"
	    }
#line 256 "cheequb.f"
	} else {
#line 257 "cheequb.f"
	    i__1 = *n;
#line 257 "cheequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "cheequb.f"
		i__2 = j;
#line 258 "cheequb.f"
		i__3 = j;
#line 258 "cheequb.f"
		i__4 = j + j * a_dim1;
#line 258 "cheequb.f"
		d__3 = ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			j * a_dim1]), abs(d__2))) * s[j];
#line 258 "cheequb.f"
		z__1.r = work[i__3].r + d__3, z__1.i = work[i__3].i;
#line 258 "cheequb.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 259 "cheequb.f"
		i__2 = *n;
#line 259 "cheequb.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 260 "cheequb.f"
		    i__3 = i__;
#line 260 "cheequb.f"
		    i__4 = i__;
#line 260 "cheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 260 "cheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[j];
#line 260 "cheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 260 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 261 "cheequb.f"
		    i__3 = j;
#line 261 "cheequb.f"
		    i__4 = j;
#line 261 "cheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 261 "cheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[i__];
#line 261 "cheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 261 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 262 "cheequb.f"
		}
#line 263 "cheequb.f"
	    }
#line 264 "cheequb.f"
	}
/*        avg = s^T beta / n */
#line 267 "cheequb.f"
	avg = 0.;
#line 268 "cheequb.f"
	i__1 = *n;
#line 268 "cheequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "cheequb.f"
	    i__2 = i__;
#line 269 "cheequb.f"
	    i__3 = i__;
#line 269 "cheequb.f"
	    z__2.r = s[i__2] * work[i__3].r, z__2.i = s[i__2] * work[i__3].i;
#line 269 "cheequb.f"
	    z__1.r = avg + z__2.r, z__1.i = z__2.i;
#line 269 "cheequb.f"
	    avg = z__1.r;
#line 270 "cheequb.f"
	}
#line 271 "cheequb.f"
	avg /= *n;
#line 273 "cheequb.f"
	std = 0.;
#line 274 "cheequb.f"
	i__1 = *n << 1;
#line 274 "cheequb.f"
	for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 275 "cheequb.f"
	    i__2 = i__;
#line 275 "cheequb.f"
	    i__3 = i__ - *n;
#line 275 "cheequb.f"
	    i__4 = i__ - *n;
#line 275 "cheequb.f"
	    z__2.r = s[i__3] * work[i__4].r, z__2.i = s[i__3] * work[i__4].i;
#line 275 "cheequb.f"
	    z__1.r = z__2.r - avg, z__1.i = z__2.i;
#line 275 "cheequb.f"
	    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 276 "cheequb.f"
	}
#line 277 "cheequb.f"
	classq_(n, &work[*n + 1], &c__1, &scale, &sumsq);
#line 278 "cheequb.f"
	std = scale * sqrt(sumsq / *n);
#line 280 "cheequb.f"
	if (std < tol * avg) {
#line 280 "cheequb.f"
	    goto L999;
#line 280 "cheequb.f"
	}
#line 282 "cheequb.f"
	i__1 = *n;
#line 282 "cheequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 283 "cheequb.f"
	    i__2 = i__ + i__ * a_dim1;
#line 283 "cheequb.f"
	    t = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + i__ * 
		    a_dim1]), abs(d__2));
#line 284 "cheequb.f"
	    si = s[i__];
#line 285 "cheequb.f"
	    c2 = (*n - 1) * t;
#line 286 "cheequb.f"
	    i__2 = *n - 2;
#line 286 "cheequb.f"
	    i__3 = i__;
#line 286 "cheequb.f"
	    d__1 = t * si;
#line 286 "cheequb.f"
	    z__2.r = work[i__3].r - d__1, z__2.i = work[i__3].i;
#line 286 "cheequb.f"
	    d__2 = (doublereal) i__2;
#line 286 "cheequb.f"
	    z__1.r = d__2 * z__2.r, z__1.i = d__2 * z__2.i;
#line 286 "cheequb.f"
	    c1 = z__1.r;
#line 287 "cheequb.f"
	    d__1 = -(t * si) * si;
#line 287 "cheequb.f"
	    i__2 = i__;
#line 287 "cheequb.f"
	    d__2 = 2.;
#line 287 "cheequb.f"
	    z__4.r = d__2 * work[i__2].r, z__4.i = d__2 * work[i__2].i;
#line 287 "cheequb.f"
	    z__3.r = si * z__4.r, z__3.i = si * z__4.i;
#line 287 "cheequb.f"
	    z__2.r = d__1 + z__3.r, z__2.i = z__3.i;
#line 287 "cheequb.f"
	    d__3 = *n * avg;
#line 287 "cheequb.f"
	    z__1.r = z__2.r - d__3, z__1.i = z__2.i;
#line 287 "cheequb.f"
	    c0 = z__1.r;
#line 288 "cheequb.f"
	    d__ = c1 * c1 - c0 * 4 * c2;
#line 290 "cheequb.f"
	    if (d__ <= 0.) {
#line 291 "cheequb.f"
		*info = -1;
#line 292 "cheequb.f"
		return 0;
#line 293 "cheequb.f"
	    }
#line 294 "cheequb.f"
	    si = c0 * -2 / (c1 + sqrt(d__));
#line 296 "cheequb.f"
	    d__ = si - s[i__];
#line 297 "cheequb.f"
	    u = 0.;
#line 298 "cheequb.f"
	    if (up) {
#line 299 "cheequb.f"
		i__2 = i__;
#line 299 "cheequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 300 "cheequb.f"
		    i__3 = j + i__ * a_dim1;
#line 300 "cheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			    i__ * a_dim1]), abs(d__2));
#line 301 "cheequb.f"
		    u += s[j] * t;
#line 302 "cheequb.f"
		    i__3 = j;
#line 302 "cheequb.f"
		    i__4 = j;
#line 302 "cheequb.f"
		    d__1 = d__ * t;
#line 302 "cheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 302 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 303 "cheequb.f"
		}
#line 304 "cheequb.f"
		i__2 = *n;
#line 304 "cheequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 305 "cheequb.f"
		    i__3 = i__ + j * a_dim1;
#line 305 "cheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			    + j * a_dim1]), abs(d__2));
#line 306 "cheequb.f"
		    u += s[j] * t;
#line 307 "cheequb.f"
		    i__3 = j;
#line 307 "cheequb.f"
		    i__4 = j;
#line 307 "cheequb.f"
		    d__1 = d__ * t;
#line 307 "cheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 307 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 308 "cheequb.f"
		}
#line 309 "cheequb.f"
	    } else {
#line 310 "cheequb.f"
		i__2 = i__;
#line 310 "cheequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 311 "cheequb.f"
		    i__3 = i__ + j * a_dim1;
#line 311 "cheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			    + j * a_dim1]), abs(d__2));
#line 312 "cheequb.f"
		    u += s[j] * t;
#line 313 "cheequb.f"
		    i__3 = j;
#line 313 "cheequb.f"
		    i__4 = j;
#line 313 "cheequb.f"
		    d__1 = d__ * t;
#line 313 "cheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 313 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 314 "cheequb.f"
		}
#line 315 "cheequb.f"
		i__2 = *n;
#line 315 "cheequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 316 "cheequb.f"
		    i__3 = j + i__ * a_dim1;
#line 316 "cheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			    i__ * a_dim1]), abs(d__2));
#line 317 "cheequb.f"
		    u += s[j] * t;
#line 318 "cheequb.f"
		    i__3 = j;
#line 318 "cheequb.f"
		    i__4 = j;
#line 318 "cheequb.f"
		    d__1 = d__ * t;
#line 318 "cheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 318 "cheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 319 "cheequb.f"
		}
#line 320 "cheequb.f"
	    }
#line 322 "cheequb.f"
	    i__2 = i__;
#line 322 "cheequb.f"
	    z__4.r = u + work[i__2].r, z__4.i = work[i__2].i;
#line 322 "cheequb.f"
	    z__3.r = d__ * z__4.r, z__3.i = d__ * z__4.i;
#line 322 "cheequb.f"
	    d__1 = (doublereal) (*n);
#line 322 "cheequb.f"
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
#line 322 "cheequb.f"
	    z__1.r = avg + z__2.r, z__1.i = z__2.i;
#line 322 "cheequb.f"
	    avg = z__1.r;
#line 323 "cheequb.f"
	    s[i__] = si;
#line 324 "cheequb.f"
	}
#line 325 "cheequb.f"
    }
#line 327 "cheequb.f"
L999:
#line 329 "cheequb.f"
    smlnum = slamch_("SAFEMIN", (ftnlen)7);
#line 330 "cheequb.f"
    bignum = 1. / smlnum;
#line 331 "cheequb.f"
    smin = bignum;
#line 332 "cheequb.f"
    smax = 0.;
#line 333 "cheequb.f"
    t = 1. / sqrt(avg);
#line 334 "cheequb.f"
    base = slamch_("B", (ftnlen)1);
#line 335 "cheequb.f"
    u = 1. / log(base);
#line 336 "cheequb.f"
    i__1 = *n;
#line 336 "cheequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 337 "cheequb.f"
	i__2 = (integer) (u * log(s[i__] * t));
#line 337 "cheequb.f"
	s[i__] = pow_di(&base, &i__2);
/* Computing MIN */
#line 338 "cheequb.f"
	d__1 = smin, d__2 = s[i__];
#line 338 "cheequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 339 "cheequb.f"
	d__1 = smax, d__2 = s[i__];
#line 339 "cheequb.f"
	smax = max(d__1,d__2);
#line 340 "cheequb.f"
    }
#line 341 "cheequb.f"
    *scond = max(smin,smlnum) / min(smax,bignum);

#line 343 "cheequb.f"
    return 0;
} /* cheequb_ */

