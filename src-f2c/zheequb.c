#line 1 "zheequb.f"
/* zheequb.f -- translated by f2c (version 20100827).
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

#line 1 "zheequb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHEEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       CHARACTER          UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       DOUBLE PRECISION   S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEQUB computes row and column scalings intended to equilibrate a */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
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

/* > \ingroup complex16HEcomputational */

/* > \par References: */
/*  ================ */
/* > */
/* >  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* >  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* >  DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* >  Tech report version: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.1679 */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zheequb_(char *uplo, integer *n, doublecomplex *a, 
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
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);


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

#line 186 "zheequb.f"
    /* Parameter adjustments */
#line 186 "zheequb.f"
    a_dim1 = *lda;
#line 186 "zheequb.f"
    a_offset = 1 + a_dim1;
#line 186 "zheequb.f"
    a -= a_offset;
#line 186 "zheequb.f"
    --s;
#line 186 "zheequb.f"
    --work;
#line 186 "zheequb.f"

#line 186 "zheequb.f"
    /* Function Body */
#line 186 "zheequb.f"
    *info = 0;
#line 187 "zheequb.f"
    if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1))) {
#line 188 "zheequb.f"
	*info = -1;
#line 189 "zheequb.f"
    } else if (*n < 0) {
#line 190 "zheequb.f"
	*info = -2;
#line 191 "zheequb.f"
    } else if (*lda < max(1,*n)) {
#line 192 "zheequb.f"
	*info = -4;
#line 193 "zheequb.f"
    }
#line 194 "zheequb.f"
    if (*info != 0) {
#line 195 "zheequb.f"
	i__1 = -(*info);
#line 195 "zheequb.f"
	xerbla_("ZHEEQUB", &i__1, (ftnlen)7);
#line 196 "zheequb.f"
	return 0;
#line 197 "zheequb.f"
    }
#line 199 "zheequb.f"
    up = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 200 "zheequb.f"
    *amax = 0.;

/*     Quick return if possible. */

#line 204 "zheequb.f"
    if (*n == 0) {
#line 205 "zheequb.f"
	*scond = 1.;
#line 206 "zheequb.f"
	return 0;
#line 207 "zheequb.f"
    }
#line 209 "zheequb.f"
    i__1 = *n;
#line 209 "zheequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 210 "zheequb.f"
	s[i__] = 0.;
#line 211 "zheequb.f"
    }
#line 213 "zheequb.f"
    *amax = 0.;
#line 214 "zheequb.f"
    if (up) {
#line 215 "zheequb.f"
	i__1 = *n;
#line 215 "zheequb.f"
	for (j = 1; j <= i__1; ++j) {
#line 216 "zheequb.f"
	    i__2 = j - 1;
#line 216 "zheequb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 217 "zheequb.f"
		i__3 = i__ + j * a_dim1;
#line 217 "zheequb.f"
		d__3 = s[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 217 "zheequb.f"
		s[i__] = max(d__3,d__4);
/* Computing MAX */
#line 218 "zheequb.f"
		i__3 = i__ + j * a_dim1;
#line 218 "zheequb.f"
		d__3 = s[j], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 218 "zheequb.f"
		s[j] = max(d__3,d__4);
/* Computing MAX */
#line 219 "zheequb.f"
		i__3 = i__ + j * a_dim1;
#line 219 "zheequb.f"
		d__3 = *amax, d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 219 "zheequb.f"
		*amax = max(d__3,d__4);
#line 220 "zheequb.f"
	    }
/* Computing MAX */
#line 221 "zheequb.f"
	    i__2 = j + j * a_dim1;
#line 221 "zheequb.f"
	    d__3 = s[j], d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 221 "zheequb.f"
	    s[j] = max(d__3,d__4);
/* Computing MAX */
#line 222 "zheequb.f"
	    i__2 = j + j * a_dim1;
#line 222 "zheequb.f"
	    d__3 = *amax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 222 "zheequb.f"
	    *amax = max(d__3,d__4);
#line 223 "zheequb.f"
	}
#line 224 "zheequb.f"
    } else {
#line 225 "zheequb.f"
	i__1 = *n;
#line 225 "zheequb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 226 "zheequb.f"
	    i__2 = j + j * a_dim1;
#line 226 "zheequb.f"
	    d__3 = s[j], d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 226 "zheequb.f"
	    s[j] = max(d__3,d__4);
/* Computing MAX */
#line 227 "zheequb.f"
	    i__2 = j + j * a_dim1;
#line 227 "zheequb.f"
	    d__3 = *amax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[j + j * a_dim1]), abs(d__2));
#line 227 "zheequb.f"
	    *amax = max(d__3,d__4);
#line 228 "zheequb.f"
	    i__2 = *n;
#line 228 "zheequb.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 229 "zheequb.f"
		i__3 = i__ + j * a_dim1;
#line 229 "zheequb.f"
		d__3 = s[i__], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 229 "zheequb.f"
		s[i__] = max(d__3,d__4);
/* Computing MAX */
#line 230 "zheequb.f"
		i__3 = i__ + j * a_dim1;
#line 230 "zheequb.f"
		d__3 = s[j], d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 230 "zheequb.f"
		s[j] = max(d__3,d__4);
/* Computing MAX */
#line 231 "zheequb.f"
		i__3 = i__ + j * a_dim1;
#line 231 "zheequb.f"
		d__3 = *amax, d__4 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&a[i__ + j * a_dim1]), abs(d__2));
#line 231 "zheequb.f"
		*amax = max(d__3,d__4);
#line 232 "zheequb.f"
	    }
#line 233 "zheequb.f"
	}
#line 234 "zheequb.f"
    }
#line 235 "zheequb.f"
    i__1 = *n;
#line 235 "zheequb.f"
    for (j = 1; j <= i__1; ++j) {
#line 236 "zheequb.f"
	s[j] = 1. / s[j];
#line 237 "zheequb.f"
    }
#line 239 "zheequb.f"
    tol = 1. / sqrt(*n * 2.);
#line 241 "zheequb.f"
    for (iter = 1; iter <= 100; ++iter) {
#line 242 "zheequb.f"
	scale = 0.;
#line 243 "zheequb.f"
	sumsq = 0.;
/*        beta = |A|s */
#line 245 "zheequb.f"
	i__1 = *n;
#line 245 "zheequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "zheequb.f"
	    i__2 = i__;
#line 246 "zheequb.f"
	    work[i__2].r = 0., work[i__2].i = 0.;
#line 247 "zheequb.f"
	}
#line 248 "zheequb.f"
	if (up) {
#line 249 "zheequb.f"
	    i__1 = *n;
#line 249 "zheequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 250 "zheequb.f"
		i__2 = j - 1;
#line 250 "zheequb.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 251 "zheequb.f"
		    i__3 = i__;
#line 251 "zheequb.f"
		    i__4 = i__;
#line 251 "zheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 251 "zheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[j];
#line 251 "zheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 251 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 252 "zheequb.f"
		    i__3 = j;
#line 252 "zheequb.f"
		    i__4 = j;
#line 252 "zheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 252 "zheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[i__];
#line 252 "zheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 252 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 253 "zheequb.f"
		}
#line 254 "zheequb.f"
		i__2 = j;
#line 254 "zheequb.f"
		i__3 = j;
#line 254 "zheequb.f"
		i__4 = j + j * a_dim1;
#line 254 "zheequb.f"
		d__3 = ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			j * a_dim1]), abs(d__2))) * s[j];
#line 254 "zheequb.f"
		z__1.r = work[i__3].r + d__3, z__1.i = work[i__3].i;
#line 254 "zheequb.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 255 "zheequb.f"
	    }
#line 256 "zheequb.f"
	} else {
#line 257 "zheequb.f"
	    i__1 = *n;
#line 257 "zheequb.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "zheequb.f"
		i__2 = j;
#line 258 "zheequb.f"
		i__3 = j;
#line 258 "zheequb.f"
		i__4 = j + j * a_dim1;
#line 258 "zheequb.f"
		d__3 = ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			j * a_dim1]), abs(d__2))) * s[j];
#line 258 "zheequb.f"
		z__1.r = work[i__3].r + d__3, z__1.i = work[i__3].i;
#line 258 "zheequb.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 259 "zheequb.f"
		i__2 = *n;
#line 259 "zheequb.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 260 "zheequb.f"
		    i__3 = i__;
#line 260 "zheequb.f"
		    i__4 = i__;
#line 260 "zheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 260 "zheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[j];
#line 260 "zheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 260 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 261 "zheequb.f"
		    i__3 = j;
#line 261 "zheequb.f"
		    i__4 = j;
#line 261 "zheequb.f"
		    i__5 = i__ + j * a_dim1;
#line 261 "zheequb.f"
		    d__3 = ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + j * a_dim1]), abs(d__2))) * s[i__];
#line 261 "zheequb.f"
		    z__1.r = work[i__4].r + d__3, z__1.i = work[i__4].i;
#line 261 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 262 "zheequb.f"
		}
#line 263 "zheequb.f"
	    }
#line 264 "zheequb.f"
	}
/*        avg = s^T beta / n */
#line 267 "zheequb.f"
	avg = 0.;
#line 268 "zheequb.f"
	i__1 = *n;
#line 268 "zheequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "zheequb.f"
	    i__2 = i__;
#line 269 "zheequb.f"
	    i__3 = i__;
#line 269 "zheequb.f"
	    z__2.r = s[i__2] * work[i__3].r, z__2.i = s[i__2] * work[i__3].i;
#line 269 "zheequb.f"
	    z__1.r = avg + z__2.r, z__1.i = z__2.i;
#line 269 "zheequb.f"
	    avg = z__1.r;
#line 270 "zheequb.f"
	}
#line 271 "zheequb.f"
	avg /= *n;
#line 273 "zheequb.f"
	std = 0.;
#line 274 "zheequb.f"
	i__1 = *n;
#line 274 "zheequb.f"
	for (i__ = *n + 1; i__ <= i__1; ++i__) {
#line 275 "zheequb.f"
	    i__2 = i__;
#line 275 "zheequb.f"
	    i__3 = i__ - *n;
#line 275 "zheequb.f"
	    i__4 = i__ - *n;
#line 275 "zheequb.f"
	    z__2.r = s[i__3] * work[i__4].r, z__2.i = s[i__3] * work[i__4].i;
#line 275 "zheequb.f"
	    z__1.r = z__2.r - avg, z__1.i = z__2.i;
#line 275 "zheequb.f"
	    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 276 "zheequb.f"
	}
#line 277 "zheequb.f"
	zlassq_(n, &work[*n + 1], &c__1, &scale, &sumsq);
#line 278 "zheequb.f"
	std = scale * sqrt(sumsq / *n);
#line 280 "zheequb.f"
	if (std < tol * avg) {
#line 280 "zheequb.f"
	    goto L999;
#line 280 "zheequb.f"
	}
#line 282 "zheequb.f"
	i__1 = *n;
#line 282 "zheequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 283 "zheequb.f"
	    i__2 = i__ + i__ * a_dim1;
#line 283 "zheequb.f"
	    t = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + i__ * 
		    a_dim1]), abs(d__2));
#line 284 "zheequb.f"
	    si = s[i__];
#line 285 "zheequb.f"
	    c2 = (*n - 1) * t;
#line 286 "zheequb.f"
	    i__2 = *n - 2;
#line 286 "zheequb.f"
	    i__3 = i__;
#line 286 "zheequb.f"
	    d__1 = t * si;
#line 286 "zheequb.f"
	    z__2.r = work[i__3].r - d__1, z__2.i = work[i__3].i;
#line 286 "zheequb.f"
	    d__2 = (doublereal) i__2;
#line 286 "zheequb.f"
	    z__1.r = d__2 * z__2.r, z__1.i = d__2 * z__2.i;
#line 286 "zheequb.f"
	    c1 = z__1.r;
#line 287 "zheequb.f"
	    d__1 = -(t * si) * si;
#line 287 "zheequb.f"
	    i__2 = i__;
#line 287 "zheequb.f"
	    d__2 = 2.;
#line 287 "zheequb.f"
	    z__4.r = d__2 * work[i__2].r, z__4.i = d__2 * work[i__2].i;
#line 287 "zheequb.f"
	    z__3.r = si * z__4.r, z__3.i = si * z__4.i;
#line 287 "zheequb.f"
	    z__2.r = d__1 + z__3.r, z__2.i = z__3.i;
#line 287 "zheequb.f"
	    d__3 = *n * avg;
#line 287 "zheequb.f"
	    z__1.r = z__2.r - d__3, z__1.i = z__2.i;
#line 287 "zheequb.f"
	    c0 = z__1.r;
#line 288 "zheequb.f"
	    d__ = c1 * c1 - c0 * 4 * c2;
#line 290 "zheequb.f"
	    if (d__ <= 0.) {
#line 291 "zheequb.f"
		*info = -1;
#line 292 "zheequb.f"
		return 0;
#line 293 "zheequb.f"
	    }
#line 294 "zheequb.f"
	    si = c0 * -2 / (c1 + sqrt(d__));
#line 296 "zheequb.f"
	    d__ = si - s[i__];
#line 297 "zheequb.f"
	    u = 0.;
#line 298 "zheequb.f"
	    if (up) {
#line 299 "zheequb.f"
		i__2 = i__;
#line 299 "zheequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 300 "zheequb.f"
		    i__3 = j + i__ * a_dim1;
#line 300 "zheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			    i__ * a_dim1]), abs(d__2));
#line 301 "zheequb.f"
		    u += s[j] * t;
#line 302 "zheequb.f"
		    i__3 = j;
#line 302 "zheequb.f"
		    i__4 = j;
#line 302 "zheequb.f"
		    d__1 = d__ * t;
#line 302 "zheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 302 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 303 "zheequb.f"
		}
#line 304 "zheequb.f"
		i__2 = *n;
#line 304 "zheequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 305 "zheequb.f"
		    i__3 = i__ + j * a_dim1;
#line 305 "zheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			    + j * a_dim1]), abs(d__2));
#line 306 "zheequb.f"
		    u += s[j] * t;
#line 307 "zheequb.f"
		    i__3 = j;
#line 307 "zheequb.f"
		    i__4 = j;
#line 307 "zheequb.f"
		    d__1 = d__ * t;
#line 307 "zheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 307 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 308 "zheequb.f"
		}
#line 309 "zheequb.f"
	    } else {
#line 310 "zheequb.f"
		i__2 = i__;
#line 310 "zheequb.f"
		for (j = 1; j <= i__2; ++j) {
#line 311 "zheequb.f"
		    i__3 = i__ + j * a_dim1;
#line 311 "zheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ 
			    + j * a_dim1]), abs(d__2));
#line 312 "zheequb.f"
		    u += s[j] * t;
#line 313 "zheequb.f"
		    i__3 = j;
#line 313 "zheequb.f"
		    i__4 = j;
#line 313 "zheequb.f"
		    d__1 = d__ * t;
#line 313 "zheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 313 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 314 "zheequb.f"
		}
#line 315 "zheequb.f"
		i__2 = *n;
#line 315 "zheequb.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 316 "zheequb.f"
		    i__3 = j + i__ * a_dim1;
#line 316 "zheequb.f"
		    t = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + 
			    i__ * a_dim1]), abs(d__2));
#line 317 "zheequb.f"
		    u += s[j] * t;
#line 318 "zheequb.f"
		    i__3 = j;
#line 318 "zheequb.f"
		    i__4 = j;
#line 318 "zheequb.f"
		    d__1 = d__ * t;
#line 318 "zheequb.f"
		    z__1.r = work[i__4].r + d__1, z__1.i = work[i__4].i;
#line 318 "zheequb.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 319 "zheequb.f"
		}
#line 320 "zheequb.f"
	    }
#line 322 "zheequb.f"
	    i__2 = i__;
#line 322 "zheequb.f"
	    z__4.r = u + work[i__2].r, z__4.i = work[i__2].i;
#line 322 "zheequb.f"
	    z__3.r = d__ * z__4.r, z__3.i = d__ * z__4.i;
#line 322 "zheequb.f"
	    d__1 = (doublereal) (*n);
#line 322 "zheequb.f"
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
#line 322 "zheequb.f"
	    z__1.r = avg + z__2.r, z__1.i = z__2.i;
#line 322 "zheequb.f"
	    avg = z__1.r;
#line 323 "zheequb.f"
	    s[i__] = si;
#line 324 "zheequb.f"
	}
#line 325 "zheequb.f"
    }
#line 327 "zheequb.f"
L999:
#line 329 "zheequb.f"
    smlnum = dlamch_("SAFEMIN", (ftnlen)7);
#line 330 "zheequb.f"
    bignum = 1. / smlnum;
#line 331 "zheequb.f"
    smin = bignum;
#line 332 "zheequb.f"
    smax = 0.;
#line 333 "zheequb.f"
    t = 1. / sqrt(avg);
#line 334 "zheequb.f"
    base = dlamch_("B", (ftnlen)1);
#line 335 "zheequb.f"
    u = 1. / log(base);
#line 336 "zheequb.f"
    i__1 = *n;
#line 336 "zheequb.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 337 "zheequb.f"
	i__2 = (integer) (u * log(s[i__] * t));
#line 337 "zheequb.f"
	s[i__] = pow_di(&base, &i__2);
/* Computing MIN */
#line 338 "zheequb.f"
	d__1 = smin, d__2 = s[i__];
#line 338 "zheequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 339 "zheequb.f"
	d__1 = smax, d__2 = s[i__];
#line 339 "zheequb.f"
	smax = max(d__1,d__2);
#line 340 "zheequb.f"
    }
#line 341 "zheequb.f"
    *scond = max(smin,smlnum) / min(smax,bignum);

#line 343 "zheequb.f"
    return 0;
} /* zheequb_ */

