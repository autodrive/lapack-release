#line 1 "csyequb.f"
/* csyequb.f -- translated by f2c (version 20100827).
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

#line 1 "csyequb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CSYEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */

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
/* > CSYEQUB computes row and column scalings intended to equilibrate a */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \date November 2017 */

/* > \ingroup complexSYcomputational */

/* > \par References: */
/*  ================ */
/* > */
/* >  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* >  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* >  DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* >  Tech report version: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.1679 */
/* > */
/*  ===================================================================== */
/* Subroutine */ int csyequb_(char *uplo, integer *n, doublecomplex *a, 
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 186 "csyequb.f"
    /* Parameter adjustments */
#line 186 "csyequb.f"
    a_dim1 = *lda;
#line 186 "csyequb.f"
    a_offset = 1 + a_dim1;
#line 186 "csyequb.f"
    a -= a_offset;
#line 186 "csyequb.f"
    --s;
#line 186 "csyequb.f"
    --work;
#line 186 "csyequb.f"

#line 186 "csyequb.f"
    /* Function Body */
#line 186 "csyequb.f"
    *info = 0;
#line 187 "csyequb.f"
    if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1))) {
#line 188 "csyequb.f"
	*info = -1;
#line 189 "csyequb.f"
    } else if (*n < 0) {
#line 190 "csyequb.f"
	*info = -2;
#line 191 "csyequb.f"
    } else if (*lda < max(1,*n)) {
#line 192 "csyequb.f"
	*info = -4;
#line 193 "csyequb.f"
    }
#line 194 "csyequb.f"
    if (*info != 0) {
#line 195 "csyequb.f"
	i__1 = -(*info);
#line 195 "csyequb.f"
	xerbla_("CSYEQUB", &i__1, (ftnlen)7);
#line 196 "csyequb.f"
	return 0;
#line 197 "csyequb.f"
    }
#line 199 "csyequb.f"
    up = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 200 "csyequb.f"
    *amax = 0.;

/*     Quick return if possible. */

#line 204 "csyequb.f"
    if (*n == 0) {
#line 205 "csyequb.f"
	*scond = 1.;
#line 206 "csyequb.f"
	return 0;
#line 207 "csyequb.f"
    }
#line 209 "csyequb.f"
    i__1 = *n;
#line 209 "csyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 210 "csyequb.f"
	s[i__] = 0.;
#line 211 "csyequb.f"
    }
#line 213 "csyequb.f"
    *amax = 0.;
#line 214 "csyequb.f"
    if (up) {
#line 215 "csyequb.f"
	i__1 = *n;
#line 215 "csyequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 216 "csyequb.f"
	    i__2 = j - 1;
#line 216 "csyequb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 217 "csyequb.f"
		i__3 = i__ + j * a_dim1;
#line 217 "csyequb.f"
		d__3 = s[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 217 "csyequb.f"
		s[i__] = max(d__3,d__4);
/* Computing MAX */
#line 218 "csyequb.f"
		i__3 = i__ + j * a_dim1;
#line 218 "csyequb.f"
		d__3 = s[j], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 218 "csyequb.f"
		s[j] = max(d__3,d__4);
/* Computing MAX */
#line 219 "csyequb.f"
		i__3 = i__ + j * a_dim1;
#line 219 "csyequb.f"
		d__3 = *amax, d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 219 "csyequb.f"
		*amax = max(d__3,d__4);
#line 220 "csyequb.f"
	    }
/* Computing MAX */
#line 221 "csyequb.f"
	    i__2 = j + j * a_dim1;
#line 221 "csyequb.f"
	    d__3 = s[j], d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 221 "csyequb.f"
	    s[j] = max(d__3,d__4);
/* Computing MAX */
#line 222 "csyequb.f"
	    i__2 = j + j * a_dim1;
#line 222 "csyequb.f"
	    d__3 = *amax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 222 "csyequb.f"
	    *amax = max(d__3,d__4);
#line 223 "csyequb.f"
	}
#line 224 "csyequb.f"
    } else {
#line 225 "csyequb.f"
	i__1 = *n;
#line 225 "csyequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 226 "csyequb.f"
	    i__2 = j + j * a_dim1;
#line 226 "csyequb.f"
	    d__3 = s[j], d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 226 "csyequb.f"
	    s[j] = max(d__3,d__4);
/* Computing MAX */
#line 227 "csyequb.f"
	    i__2 = j + j * a_dim1;
#line 227 "csyequb.f"
	    d__3 = *amax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 227 "csyequb.f"
	    *amax = max(d__3,d__4);
#line 228 "csyequb.f"
	    i__2 = *n;
#line 228 "csyequb.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 229 "csyequb.f"
		i__3 = i__ + j * a_dim1;
#line 229 "csyequb.f"
		d__3 = s[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 229 "csyequb.f"
		s[i__] = max(d__3,d__4);
/* Computing MAX */
#line 230 "csyequb.f"
		i__3 = i__ + j * a_dim1;
#line 230 "csyequb.f"
		d__3 = s[j], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 230 "csyequb.f"
		s[j] = max(d__3,d__4);
/* Computing MAX */
#line 231 "csyequb.f"
		i__3 = i__ + j * a_dim1;
#line 231 "csyequb.f"
		d__3 = *amax, d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 231 "csyequb.f"
		*amax = max(d__3,d__4);
#line 232 "csyequb.f"
	    }
#line 233 "csyequb.f"
	}
#line 234 "csyequb.f"
    }
#line 235 "csyequb.f"
    i__1 = *n;
#line 235 "csyequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 236 "csyequb.f"
	s[j] = 1. / s[j];
#line 237 "csyequb.f"
    }
#line 239 "csyequb.f"
    tol = 1. / sqrt(*n * 2.);
#line 241 "csyequb.f"
    for (iter = 1; iter <= 100; ++iter) {
#line 242 "csyequb.f"
	scale = 0.;
#line 243 "csyequb.f"
	sumsq = 0.;
/*        beta = |A|s */
#line 245 "csyequb.f"
	i__1 = *n;
#line 245 "csyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "csyequb.f"
	    i__2 = i__;
#line 246 "csyequb.f"
	    work[i__2].r = 0., work[i__2].i = 0.;
#line 247 "csyequb.f"
	}
#line 248 "csyequb.f"
	if (up) {
#line 249 "csyequb.f"
	    i__1 = *n;
#line 249 "csyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 250 "csyequb.f"
		i__2 = j - 1;
#line 250 "csyequb.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 251 "csyequb.f"
		    i__3 = i__;
#line 251 "csyequb.f"
		    i__4 = i__;
#line 251 "csyequb.f"
		    i__5 = i__ + j * a_dim1;
#line 251 "csyequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[j];
#line 251 "csyequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 251 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 252 "csyequb.f"
		    i__3 = j;
#line 252 "csyequb.f"
		    i__4 = j;
#line 252 "csyequb.f"
		    i__5 = i__ + j * a_dim1;
#line 252 "csyequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[i__];
#line 252 "csyequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 252 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 253 "csyequb.f"
		}
#line 254 "csyequb.f"
		i__2 = j;
#line 254 "csyequb.f"
		i__3 = j;
#line 254 "csyequb.f"
		i__4 = j + j * a_dim1;
#line 254 "csyequb.f"
		d__3 = ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			j * a_dim1]), abs(d__2))) * s[j];
#line 254 "csyequb.f"
		z__1.r = work[i__3].r + d__3, z__1.i = work[i__3].i;
#line 254 "csyequb.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 255 "csyequb.f"
	    }
#line 256 "csyequb.f"
	} else {
#line 257 "csyequb.f"
	    i__1 = *n;
#line 257 "csyequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "csyequb.f"
		i__2 = j;
#line 258 "csyequb.f"
		i__3 = j;
#line 258 "csyequb.f"
		i__4 = j + j * a_dim1;
#line 258 "csyequb.f"
		d__3 = ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			j * a_dim1]), abs(d__2))) * s[j];
#line 258 "csyequb.f"
		z__1.r = work[i__3].r + d__3, z__1.i = work[i__3].i;
#line 258 "csyequb.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 259 "csyequb.f"
		i__2 = *n;
#line 259 "csyequb.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 260 "csyequb.f"
		    i__3 = i__;
#line 260 "csyequb.f"
		    i__4 = i__;
#line 260 "csyequb.f"
		    i__5 = i__ + j * a_dim1;
#line 260 "csyequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[j];
#line 260 "csyequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 260 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 261 "csyequb.f"
		    i__3 = j;
#line 261 "csyequb.f"
		    i__4 = j;
#line 261 "csyequb.f"
		    i__5 = i__ + j * a_dim1;
#line 261 "csyequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[i__];
#line 261 "csyequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 261 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 262 "csyequb.f"
		}
#line 263 "csyequb.f"
	    }
#line 264 "csyequb.f"
	}
/*        avg = s^T beta / n */
#line 267 "csyequb.f"
	avg = 0.;
#line 268 "csyequb.f"
	i__1 = *n;
#line 268 "csyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "csyequb.f"
	    i__2 = i__;
#line 269 "csyequb.f"
	    i__3 = i__;
#line 269 "csyequb.f"
	    z__2.r = s[i__2] * work[i__3].r, z__2.i = s[i__2] * work[i__3].i;
#line 269 "csyequb.f"
	    z__1.r = avg + z__2.r, z__1.i = z__2.i;
#line 269 "csyequb.f"
	    avg = z__1.r;
#line 270 "csyequb.f"
	}
#line 271 "csyequb.f"
	avg /= *n;
#line 273 "csyequb.f"
	std = 0.;
#line 274 "csyequb.f"
	i__1 = *n << 1;
#line 274 "csyequb.f"
	for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 275 "csyequb.f"
	    i__2 = i__;
#line 275 "csyequb.f"
	    i__3 = i__ - *n;
#line 275 "csyequb.f"
	    i__4 = i__ - *n;
#line 275 "csyequb.f"
	    z__2.r = s[i__3] * work[i__4].r, z__2.i = s[i__3] * work[i__4].i;
#line 275 "csyequb.f"
	    z__1.r = z__2.r - avg, z__1.i = z__2.i;
#line 275 "csyequb.f"
	    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 276 "csyequb.f"
	}
#line 277 "csyequb.f"
	classq_(n, &work[*n + 1], &c__1, &scale, &sumsq);
#line 278 "csyequb.f"
	std = scale * sqrt(sumsq / *n);
#line 280 "csyequb.f"
	if (std < tol * avg) {
#line 280 "csyequb.f"
	    goto L999;
#line 280 "csyequb.f"
	}
#line 282 "csyequb.f"
	i__1 = *n;
#line 282 "csyequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 283 "csyequb.f"
	    i__2 = i__ + i__ * a_dim1;
#line 283 "csyequb.f"
	    t = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + i__ * 
		    a_dim1]), abs(d__2));
#line 284 "csyequb.f"
	    si = s[i__];
#line 285 "csyequb.f"
	    c2 = (*n - 1) * t;
#line 286 "csyequb.f"
	    i__2 = *n - 2;
#line 286 "csyequb.f"
	    i__3 = i__;
#line 286 "csyequb.f"
	    d__1 = t * si;
#line 286 "csyequb.f"
	    z__2.r = work[i__3].r - d__1, z__2.i = work[i__3].i;
#line 286 "csyequb.f"
	    d__2 = (doublereal) i__2;
#line 286 "csyequb.f"
	    z__1.r = d__2 * z__2.r, z__1.i = d__2 * z__2.i;
#line 286 "csyequb.f"
	    c1 = z__1.r;
#line 287 "csyequb.f"
	    d__1 = -(t * si) * si;
#line 287 "csyequb.f"
	    i__2 = i__;
#line 287 "csyequb.f"
	    d__2 = 2.;
#line 287 "csyequb.f"
	    z__4.r = d__2 * work[i__2].r, z__4.i = d__2 * work[i__2].i;
#line 287 "csyequb.f"
	    z__3.r = si * z__4.r, z__3.i = si * z__4.i;
#line 287 "csyequb.f"
	    z__2.r = d__1 + z__3.r, z__2.i = z__3.i;
#line 287 "csyequb.f"
	    d__3 = *n * avg;
#line 287 "csyequb.f"
	    z__1.r = z__2.r - d__3, z__1.i = z__2.i;
#line 287 "csyequb.f"
	    c0 = z__1.r;
#line 288 "csyequb.f"
	    d__ = c1 * c1 - c0 * 4 * c2;
#line 290 "csyequb.f"
	    if (d__ <= 0.) {
#line 291 "csyequb.f"
		*info = -1;
#line 292 "csyequb.f"
		return 0;
#line 293 "csyequb.f"
	    }
#line 294 "csyequb.f"
	    si = c0 * -2 / (c1 + sqrt(d__));
#line 296 "csyequb.f"
	    d__ = si - s[i__];
#line 297 "csyequb.f"
	    u = 0.;
#line 298 "csyequb.f"
	    if (up) {
#line 299 "csyequb.f"
		i__2 = i__;
#line 299 "csyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 300 "csyequb.f"
		    i__3 = j + i__ * a_dim1;
#line 300 "csyequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			    i__ * a_dim1]), abs(d__2));
#line 301 "csyequb.f"
		    u += s[j] * t;
#line 302 "csyequb.f"
		    i__3 = j;
#line 302 "csyequb.f"
		    i__4 = j;
#line 302 "csyequb.f"
		    d__1 = d__ * t;
#line 302 "csyequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 302 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 303 "csyequb.f"
		}
#line 304 "csyequb.f"
		i__2 = *n;
#line 304 "csyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 305 "csyequb.f"
		    i__3 = i__ + j * a_dim1;
#line 305 "csyequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			    + j * a_dim1]), abs(d__2));
#line 306 "csyequb.f"
		    u += s[j] * t;
#line 307 "csyequb.f"
		    i__3 = j;
#line 307 "csyequb.f"
		    i__4 = j;
#line 307 "csyequb.f"
		    d__1 = d__ * t;
#line 307 "csyequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 307 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 308 "csyequb.f"
		}
#line 309 "csyequb.f"
	    } else {
#line 310 "csyequb.f"
		i__2 = i__;
#line 310 "csyequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 311 "csyequb.f"
		    i__3 = i__ + j * a_dim1;
#line 311 "csyequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			    + j * a_dim1]), abs(d__2));
#line 312 "csyequb.f"
		    u += s[j] * t;
#line 313 "csyequb.f"
		    i__3 = j;
#line 313 "csyequb.f"
		    i__4 = j;
#line 313 "csyequb.f"
		    d__1 = d__ * t;
#line 313 "csyequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 313 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 314 "csyequb.f"
		}
#line 315 "csyequb.f"
		i__2 = *n;
#line 315 "csyequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 316 "csyequb.f"
		    i__3 = j + i__ * a_dim1;
#line 316 "csyequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			    i__ * a_dim1]), abs(d__2));
#line 317 "csyequb.f"
		    u += s[j] * t;
#line 318 "csyequb.f"
		    i__3 = j;
#line 318 "csyequb.f"
		    i__4 = j;
#line 318 "csyequb.f"
		    d__1 = d__ * t;
#line 318 "csyequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 318 "csyequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 319 "csyequb.f"
		}
#line 320 "csyequb.f"
	    }
#line 322 "csyequb.f"
	    i__2 = i__;
#line 322 "csyequb.f"
	    z__4.r = u + work[i__2].r, z__4.i = work[i__2].i;
#line 322 "csyequb.f"
	    z__3.r = d__ * z__4.r, z__3.i = d__ * z__4.i;
#line 322 "csyequb.f"
	    d__1 = (doublereal) (*n);
#line 322 "csyequb.f"
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
#line 322 "csyequb.f"
	    z__1.r = avg + z__2.r, z__1.i = z__2.i;
#line 322 "csyequb.f"
	    avg = z__1.r;
#line 323 "csyequb.f"
	    s[i__] = si;
#line 324 "csyequb.f"
	}
#line 325 "csyequb.f"
    }
#line 327 "csyequb.f"
L999:
#line 329 "csyequb.f"
    smlnum = slamch_("SAFEMIN", (ftnlen)7);
#line 330 "csyequb.f"
    bignum = 1. / smlnum;
#line 331 "csyequb.f"
    smin = bignum;
#line 332 "csyequb.f"
    smax = 0.;
#line 333 "csyequb.f"
    t = 1. / sqrt(avg);
#line 334 "csyequb.f"
    base = slamch_("B", (ftnlen)1);
#line 335 "csyequb.f"
    u = 1. / log(base);
#line 336 "csyequb.f"
    i__1 = *n;
#line 336 "csyequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 337 "csyequb.f"
	i__2 = (integer) (u * log(s[i__] * t));
#line 337 "csyequb.f"
	s[i__] = pow_di(&base, &i__2);
/* Computing MIN */
#line 338 "csyequb.f"
	d__1 = smin, d__2 = s[i__];
#line 338 "csyequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 339 "csyequb.f"
	d__1 = smax, d__2 = s[i__];
#line 339 "csyequb.f"
	smax = max(d__1,d__2);
#line 340 "csyequb.f"
    }
#line 341 "csyequb.f"
    *scond = max(smin,smlnum) / min(smax,bignum);

#line 343 "csyequb.f"
    return 0;
} /* csyequb_ */

