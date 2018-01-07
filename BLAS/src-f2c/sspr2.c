#line 1 "sspr2.f"
/* sspr2.f -- translated by f2c (version 20100827).
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

#line 1 "sspr2.f"
/* > \brief \b SSPR2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPR2  performs the symmetric rank 2 operation */
/* > */
/* >    A := alpha*x*y**T + alpha*y*x**T + A, */
/* > */
/* > where alpha is a scalar, x and y are n element vectors and A is an */
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
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
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
/* Subroutine */ int sspr2_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    static doublereal temp1, temp2;
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

#line 179 "sspr2.f"
    /* Parameter adjustments */
#line 179 "sspr2.f"
    --ap;
#line 179 "sspr2.f"
    --y;
#line 179 "sspr2.f"
    --x;
#line 179 "sspr2.f"

#line 179 "sspr2.f"
    /* Function Body */
#line 179 "sspr2.f"
    info = 0;
#line 180 "sspr2.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 181 "sspr2.f"
	info = 1;
#line 182 "sspr2.f"
    } else if (*n < 0) {
#line 183 "sspr2.f"
	info = 2;
#line 184 "sspr2.f"
    } else if (*incx == 0) {
#line 185 "sspr2.f"
	info = 5;
#line 186 "sspr2.f"
    } else if (*incy == 0) {
#line 187 "sspr2.f"
	info = 7;
#line 188 "sspr2.f"
    }
#line 189 "sspr2.f"
    if (info != 0) {
#line 190 "sspr2.f"
	xerbla_("SSPR2 ", &info, (ftnlen)6);
#line 191 "sspr2.f"
	return 0;
#line 192 "sspr2.f"
    }

/*     Quick return if possible. */

#line 196 "sspr2.f"
    if (*n == 0 || *alpha == 0.) {
#line 196 "sspr2.f"
	return 0;
#line 196 "sspr2.f"
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

#line 201 "sspr2.f"
    if (*incx != 1 || *incy != 1) {
#line 202 "sspr2.f"
	if (*incx > 0) {
#line 203 "sspr2.f"
	    kx = 1;
#line 204 "sspr2.f"
	} else {
#line 205 "sspr2.f"
	    kx = 1 - (*n - 1) * *incx;
#line 206 "sspr2.f"
	}
#line 207 "sspr2.f"
	if (*incy > 0) {
#line 208 "sspr2.f"
	    ky = 1;
#line 209 "sspr2.f"
	} else {
#line 210 "sspr2.f"
	    ky = 1 - (*n - 1) * *incy;
#line 211 "sspr2.f"
	}
#line 212 "sspr2.f"
	jx = kx;
#line 213 "sspr2.f"
	jy = ky;
#line 214 "sspr2.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 219 "sspr2.f"
    kk = 1;
#line 220 "sspr2.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 224 "sspr2.f"
	if (*incx == 1 && *incy == 1) {
#line 225 "sspr2.f"
	    i__1 = *n;
#line 225 "sspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 226 "sspr2.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 227 "sspr2.f"
		    temp1 = *alpha * y[j];
#line 228 "sspr2.f"
		    temp2 = *alpha * x[j];
#line 229 "sspr2.f"
		    k = kk;
#line 230 "sspr2.f"
		    i__2 = j;
#line 230 "sspr2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 231 "sspr2.f"
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
#line 232 "sspr2.f"
			++k;
#line 233 "sspr2.f"
/* L10: */
#line 233 "sspr2.f"
		    }
#line 234 "sspr2.f"
		}
#line 235 "sspr2.f"
		kk += j;
#line 236 "sspr2.f"
/* L20: */
#line 236 "sspr2.f"
	    }
#line 237 "sspr2.f"
	} else {
#line 238 "sspr2.f"
	    i__1 = *n;
#line 238 "sspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "sspr2.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 240 "sspr2.f"
		    temp1 = *alpha * y[jy];
#line 241 "sspr2.f"
		    temp2 = *alpha * x[jx];
#line 242 "sspr2.f"
		    ix = kx;
#line 243 "sspr2.f"
		    iy = ky;
#line 244 "sspr2.f"
		    i__2 = kk + j - 1;
#line 244 "sspr2.f"
		    for (k = kk; k <= i__2; ++k) {
#line 245 "sspr2.f"
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
#line 246 "sspr2.f"
			ix += *incx;
#line 247 "sspr2.f"
			iy += *incy;
#line 248 "sspr2.f"
/* L30: */
#line 248 "sspr2.f"
		    }
#line 249 "sspr2.f"
		}
#line 250 "sspr2.f"
		jx += *incx;
#line 251 "sspr2.f"
		jy += *incy;
#line 252 "sspr2.f"
		kk += j;
#line 253 "sspr2.f"
/* L40: */
#line 253 "sspr2.f"
	    }
#line 254 "sspr2.f"
	}
#line 255 "sspr2.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 259 "sspr2.f"
	if (*incx == 1 && *incy == 1) {
#line 260 "sspr2.f"
	    i__1 = *n;
#line 260 "sspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 261 "sspr2.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 262 "sspr2.f"
		    temp1 = *alpha * y[j];
#line 263 "sspr2.f"
		    temp2 = *alpha * x[j];
#line 264 "sspr2.f"
		    k = kk;
#line 265 "sspr2.f"
		    i__2 = *n;
#line 265 "sspr2.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 266 "sspr2.f"
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
#line 267 "sspr2.f"
			++k;
#line 268 "sspr2.f"
/* L50: */
#line 268 "sspr2.f"
		    }
#line 269 "sspr2.f"
		}
#line 270 "sspr2.f"
		kk = kk + *n - j + 1;
#line 271 "sspr2.f"
/* L60: */
#line 271 "sspr2.f"
	    }
#line 272 "sspr2.f"
	} else {
#line 273 "sspr2.f"
	    i__1 = *n;
#line 273 "sspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 274 "sspr2.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 275 "sspr2.f"
		    temp1 = *alpha * y[jy];
#line 276 "sspr2.f"
		    temp2 = *alpha * x[jx];
#line 277 "sspr2.f"
		    ix = jx;
#line 278 "sspr2.f"
		    iy = jy;
#line 279 "sspr2.f"
		    i__2 = kk + *n - j;
#line 279 "sspr2.f"
		    for (k = kk; k <= i__2; ++k) {
#line 280 "sspr2.f"
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
#line 281 "sspr2.f"
			ix += *incx;
#line 282 "sspr2.f"
			iy += *incy;
#line 283 "sspr2.f"
/* L70: */
#line 283 "sspr2.f"
		    }
#line 284 "sspr2.f"
		}
#line 285 "sspr2.f"
		jx += *incx;
#line 286 "sspr2.f"
		jy += *incy;
#line 287 "sspr2.f"
		kk = kk + *n - j + 1;
#line 288 "sspr2.f"
/* L80: */
#line 288 "sspr2.f"
	    }
#line 289 "sspr2.f"
	}
#line 290 "sspr2.f"
    }

#line 292 "sspr2.f"
    return 0;

/*     End of SSPR2 . */

} /* sspr2_ */

