#line 1 "dspr2.f"
/* dspr2.f -- translated by f2c (version 20100827).
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

#line 1 "dspr2.f"
/* > \brief \b DSPR2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPR2  performs the symmetric rank 2 operation */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension at least */
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
/* >          Y is DOUBLE PRECISION array, dimension at least */
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
/* >          AP is DOUBLE PRECISION array, dimension at least */
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

/* > \date December 2016 */

/* > \ingroup double_blas_level2 */

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
/* Subroutine */ int dspr2_(char *uplo, integer *n, doublereal *alpha, 
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


/*  -- Reference BLAS level2 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
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

/*     Test the input parameters. */

#line 179 "dspr2.f"
    /* Parameter adjustments */
#line 179 "dspr2.f"
    --ap;
#line 179 "dspr2.f"
    --y;
#line 179 "dspr2.f"
    --x;
#line 179 "dspr2.f"

#line 179 "dspr2.f"
    /* Function Body */
#line 179 "dspr2.f"
    info = 0;
#line 180 "dspr2.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 181 "dspr2.f"
	info = 1;
#line 182 "dspr2.f"
    } else if (*n < 0) {
#line 183 "dspr2.f"
	info = 2;
#line 184 "dspr2.f"
    } else if (*incx == 0) {
#line 185 "dspr2.f"
	info = 5;
#line 186 "dspr2.f"
    } else if (*incy == 0) {
#line 187 "dspr2.f"
	info = 7;
#line 188 "dspr2.f"
    }
#line 189 "dspr2.f"
    if (info != 0) {
#line 190 "dspr2.f"
	xerbla_("DSPR2 ", &info, (ftnlen)6);
#line 191 "dspr2.f"
	return 0;
#line 192 "dspr2.f"
    }

/*     Quick return if possible. */

#line 196 "dspr2.f"
    if (*n == 0 || *alpha == 0.) {
#line 196 "dspr2.f"
	return 0;
#line 196 "dspr2.f"
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

#line 201 "dspr2.f"
    if (*incx != 1 || *incy != 1) {
#line 202 "dspr2.f"
	if (*incx > 0) {
#line 203 "dspr2.f"
	    kx = 1;
#line 204 "dspr2.f"
	} else {
#line 205 "dspr2.f"
	    kx = 1 - (*n - 1) * *incx;
#line 206 "dspr2.f"
	}
#line 207 "dspr2.f"
	if (*incy > 0) {
#line 208 "dspr2.f"
	    ky = 1;
#line 209 "dspr2.f"
	} else {
#line 210 "dspr2.f"
	    ky = 1 - (*n - 1) * *incy;
#line 211 "dspr2.f"
	}
#line 212 "dspr2.f"
	jx = kx;
#line 213 "dspr2.f"
	jy = ky;
#line 214 "dspr2.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 219 "dspr2.f"
    kk = 1;
#line 220 "dspr2.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 224 "dspr2.f"
	if (*incx == 1 && *incy == 1) {
#line 225 "dspr2.f"
	    i__1 = *n;
#line 225 "dspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 226 "dspr2.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 227 "dspr2.f"
		    temp1 = *alpha * y[j];
#line 228 "dspr2.f"
		    temp2 = *alpha * x[j];
#line 229 "dspr2.f"
		    k = kk;
#line 230 "dspr2.f"
		    i__2 = j;
#line 230 "dspr2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 231 "dspr2.f"
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
#line 232 "dspr2.f"
			++k;
#line 233 "dspr2.f"
/* L10: */
#line 233 "dspr2.f"
		    }
#line 234 "dspr2.f"
		}
#line 235 "dspr2.f"
		kk += j;
#line 236 "dspr2.f"
/* L20: */
#line 236 "dspr2.f"
	    }
#line 237 "dspr2.f"
	} else {
#line 238 "dspr2.f"
	    i__1 = *n;
#line 238 "dspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "dspr2.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 240 "dspr2.f"
		    temp1 = *alpha * y[jy];
#line 241 "dspr2.f"
		    temp2 = *alpha * x[jx];
#line 242 "dspr2.f"
		    ix = kx;
#line 243 "dspr2.f"
		    iy = ky;
#line 244 "dspr2.f"
		    i__2 = kk + j - 1;
#line 244 "dspr2.f"
		    for (k = kk; k <= i__2; ++k) {
#line 245 "dspr2.f"
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
#line 246 "dspr2.f"
			ix += *incx;
#line 247 "dspr2.f"
			iy += *incy;
#line 248 "dspr2.f"
/* L30: */
#line 248 "dspr2.f"
		    }
#line 249 "dspr2.f"
		}
#line 250 "dspr2.f"
		jx += *incx;
#line 251 "dspr2.f"
		jy += *incy;
#line 252 "dspr2.f"
		kk += j;
#line 253 "dspr2.f"
/* L40: */
#line 253 "dspr2.f"
	    }
#line 254 "dspr2.f"
	}
#line 255 "dspr2.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 259 "dspr2.f"
	if (*incx == 1 && *incy == 1) {
#line 260 "dspr2.f"
	    i__1 = *n;
#line 260 "dspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 261 "dspr2.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 262 "dspr2.f"
		    temp1 = *alpha * y[j];
#line 263 "dspr2.f"
		    temp2 = *alpha * x[j];
#line 264 "dspr2.f"
		    k = kk;
#line 265 "dspr2.f"
		    i__2 = *n;
#line 265 "dspr2.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 266 "dspr2.f"
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
#line 267 "dspr2.f"
			++k;
#line 268 "dspr2.f"
/* L50: */
#line 268 "dspr2.f"
		    }
#line 269 "dspr2.f"
		}
#line 270 "dspr2.f"
		kk = kk + *n - j + 1;
#line 271 "dspr2.f"
/* L60: */
#line 271 "dspr2.f"
	    }
#line 272 "dspr2.f"
	} else {
#line 273 "dspr2.f"
	    i__1 = *n;
#line 273 "dspr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 274 "dspr2.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 275 "dspr2.f"
		    temp1 = *alpha * y[jy];
#line 276 "dspr2.f"
		    temp2 = *alpha * x[jx];
#line 277 "dspr2.f"
		    ix = jx;
#line 278 "dspr2.f"
		    iy = jy;
#line 279 "dspr2.f"
		    i__2 = kk + *n - j;
#line 279 "dspr2.f"
		    for (k = kk; k <= i__2; ++k) {
#line 280 "dspr2.f"
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
#line 281 "dspr2.f"
			ix += *incx;
#line 282 "dspr2.f"
			iy += *incy;
#line 283 "dspr2.f"
/* L70: */
#line 283 "dspr2.f"
		    }
#line 284 "dspr2.f"
		}
#line 285 "dspr2.f"
		jx += *incx;
#line 286 "dspr2.f"
		jy += *incy;
#line 287 "dspr2.f"
		kk = kk + *n - j + 1;
#line 288 "dspr2.f"
/* L80: */
#line 288 "dspr2.f"
	    }
#line 289 "dspr2.f"
	}
#line 290 "dspr2.f"
    }

#line 292 "dspr2.f"
    return 0;

/*     End of DSPR2 . */

} /* dspr2_ */

