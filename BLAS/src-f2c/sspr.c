#line 1 "sspr.f"
/* sspr.f -- translated by f2c (version 20100827).
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

#line 1 "sspr.f"
/* > \brief \b SSPR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER INCX,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPR    performs the symmetric rank 1 operation */
/* > */
/* >    A := alpha*x*x**T + A, */
/* > */
/* > where alpha is a real scalar, x is an n element vector and A is an */
/* > n by n symmetric matrix, supplied in packed form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the matrix A is supplied in the packed */
/* >           array AP as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   The upper triangular part of A is */
/* >                                  supplied in AP. */
/* > */
/* >              UPLO = 'L' or 'l'   The lower triangular part of A is */
/* >                                  supplied in AP. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the upper triangular part of the */
/* >           updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the lower triangular part of the */
/* >           updated matrix. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sspr_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
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

/*     Test the input parameters. */

#line 164 "sspr.f"
    /* Parameter adjustments */
#line 164 "sspr.f"
    --ap;
#line 164 "sspr.f"
    --x;
#line 164 "sspr.f"

#line 164 "sspr.f"
    /* Function Body */
#line 164 "sspr.f"
    info = 0;
#line 165 "sspr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 166 "sspr.f"
	info = 1;
#line 167 "sspr.f"
    } else if (*n < 0) {
#line 168 "sspr.f"
	info = 2;
#line 169 "sspr.f"
    } else if (*incx == 0) {
#line 170 "sspr.f"
	info = 5;
#line 171 "sspr.f"
    }
#line 172 "sspr.f"
    if (info != 0) {
#line 173 "sspr.f"
	xerbla_("SSPR  ", &info, (ftnlen)6);
#line 174 "sspr.f"
	return 0;
#line 175 "sspr.f"
    }

/*     Quick return if possible. */

#line 179 "sspr.f"
    if (*n == 0 || *alpha == 0.) {
#line 179 "sspr.f"
	return 0;
#line 179 "sspr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 183 "sspr.f"
    if (*incx <= 0) {
#line 184 "sspr.f"
	kx = 1 - (*n - 1) * *incx;
#line 185 "sspr.f"
    } else if (*incx != 1) {
#line 186 "sspr.f"
	kx = 1;
#line 187 "sspr.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 192 "sspr.f"
    kk = 1;
#line 193 "sspr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 197 "sspr.f"
	if (*incx == 1) {
#line 198 "sspr.f"
	    i__1 = *n;
#line 198 "sspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 199 "sspr.f"
		if (x[j] != 0.) {
#line 200 "sspr.f"
		    temp = *alpha * x[j];
#line 201 "sspr.f"
		    k = kk;
#line 202 "sspr.f"
		    i__2 = j;
#line 202 "sspr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 203 "sspr.f"
			ap[k] += x[i__] * temp;
#line 204 "sspr.f"
			++k;
#line 205 "sspr.f"
/* L10: */
#line 205 "sspr.f"
		    }
#line 206 "sspr.f"
		}
#line 207 "sspr.f"
		kk += j;
#line 208 "sspr.f"
/* L20: */
#line 208 "sspr.f"
	    }
#line 209 "sspr.f"
	} else {
#line 210 "sspr.f"
	    jx = kx;
#line 211 "sspr.f"
	    i__1 = *n;
#line 211 "sspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 212 "sspr.f"
		if (x[jx] != 0.) {
#line 213 "sspr.f"
		    temp = *alpha * x[jx];
#line 214 "sspr.f"
		    ix = kx;
#line 215 "sspr.f"
		    i__2 = kk + j - 1;
#line 215 "sspr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 216 "sspr.f"
			ap[k] += x[ix] * temp;
#line 217 "sspr.f"
			ix += *incx;
#line 218 "sspr.f"
/* L30: */
#line 218 "sspr.f"
		    }
#line 219 "sspr.f"
		}
#line 220 "sspr.f"
		jx += *incx;
#line 221 "sspr.f"
		kk += j;
#line 222 "sspr.f"
/* L40: */
#line 222 "sspr.f"
	    }
#line 223 "sspr.f"
	}
#line 224 "sspr.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 228 "sspr.f"
	if (*incx == 1) {
#line 229 "sspr.f"
	    i__1 = *n;
#line 229 "sspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 230 "sspr.f"
		if (x[j] != 0.) {
#line 231 "sspr.f"
		    temp = *alpha * x[j];
#line 232 "sspr.f"
		    k = kk;
#line 233 "sspr.f"
		    i__2 = *n;
#line 233 "sspr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 234 "sspr.f"
			ap[k] += x[i__] * temp;
#line 235 "sspr.f"
			++k;
#line 236 "sspr.f"
/* L50: */
#line 236 "sspr.f"
		    }
#line 237 "sspr.f"
		}
#line 238 "sspr.f"
		kk = kk + *n - j + 1;
#line 239 "sspr.f"
/* L60: */
#line 239 "sspr.f"
	    }
#line 240 "sspr.f"
	} else {
#line 241 "sspr.f"
	    jx = kx;
#line 242 "sspr.f"
	    i__1 = *n;
#line 242 "sspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 243 "sspr.f"
		if (x[jx] != 0.) {
#line 244 "sspr.f"
		    temp = *alpha * x[jx];
#line 245 "sspr.f"
		    ix = jx;
#line 246 "sspr.f"
		    i__2 = kk + *n - j;
#line 246 "sspr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 247 "sspr.f"
			ap[k] += x[ix] * temp;
#line 248 "sspr.f"
			ix += *incx;
#line 249 "sspr.f"
/* L70: */
#line 249 "sspr.f"
		    }
#line 250 "sspr.f"
		}
#line 251 "sspr.f"
		jx += *incx;
#line 252 "sspr.f"
		kk = kk + *n - j + 1;
#line 253 "sspr.f"
/* L80: */
#line 253 "sspr.f"
	    }
#line 254 "sspr.f"
	}
#line 255 "sspr.f"
    }

#line 257 "sspr.f"
    return 0;

/*     End of SSPR  . */

} /* sspr_ */

