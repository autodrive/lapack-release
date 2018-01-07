#line 1 "dsymv.f"
/* dsymv.f -- translated by f2c (version 20100827).
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

#line 1 "dsymv.f"
/* > \brief \b DSYMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYMV  performs the matrix-vector  operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array A is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/* >                                  is to be referenced. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( LDA, N ) */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           upper triangular part of A is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
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
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. On exit, Y is overwritten by the updated */
/* >           vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
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
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsymv_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal 
	*beta, doublereal *y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
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
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

#line 192 "dsymv.f"
    /* Parameter adjustments */
#line 192 "dsymv.f"
    a_dim1 = *lda;
#line 192 "dsymv.f"
    a_offset = 1 + a_dim1;
#line 192 "dsymv.f"
    a -= a_offset;
#line 192 "dsymv.f"
    --x;
#line 192 "dsymv.f"
    --y;
#line 192 "dsymv.f"

#line 192 "dsymv.f"
    /* Function Body */
#line 192 "dsymv.f"
    info = 0;
#line 193 "dsymv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 194 "dsymv.f"
	info = 1;
#line 195 "dsymv.f"
    } else if (*n < 0) {
#line 196 "dsymv.f"
	info = 2;
#line 197 "dsymv.f"
    } else if (*lda < max(1,*n)) {
#line 198 "dsymv.f"
	info = 5;
#line 199 "dsymv.f"
    } else if (*incx == 0) {
#line 200 "dsymv.f"
	info = 7;
#line 201 "dsymv.f"
    } else if (*incy == 0) {
#line 202 "dsymv.f"
	info = 10;
#line 203 "dsymv.f"
    }
#line 204 "dsymv.f"
    if (info != 0) {
#line 205 "dsymv.f"
	xerbla_("DSYMV ", &info, (ftnlen)6);
#line 206 "dsymv.f"
	return 0;
#line 207 "dsymv.f"
    }

/*     Quick return if possible. */

#line 211 "dsymv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 211 "dsymv.f"
	return 0;
#line 211 "dsymv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 215 "dsymv.f"
    if (*incx > 0) {
#line 216 "dsymv.f"
	kx = 1;
#line 217 "dsymv.f"
    } else {
#line 218 "dsymv.f"
	kx = 1 - (*n - 1) * *incx;
#line 219 "dsymv.f"
    }
#line 220 "dsymv.f"
    if (*incy > 0) {
#line 221 "dsymv.f"
	ky = 1;
#line 222 "dsymv.f"
    } else {
#line 223 "dsymv.f"
	ky = 1 - (*n - 1) * *incy;
#line 224 "dsymv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

#line 232 "dsymv.f"
    if (*beta != 1.) {
#line 233 "dsymv.f"
	if (*incy == 1) {
#line 234 "dsymv.f"
	    if (*beta == 0.) {
#line 235 "dsymv.f"
		i__1 = *n;
#line 235 "dsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "dsymv.f"
		    y[i__] = 0.;
#line 237 "dsymv.f"
/* L10: */
#line 237 "dsymv.f"
		}
#line 238 "dsymv.f"
	    } else {
#line 239 "dsymv.f"
		i__1 = *n;
#line 239 "dsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "dsymv.f"
		    y[i__] = *beta * y[i__];
#line 241 "dsymv.f"
/* L20: */
#line 241 "dsymv.f"
		}
#line 242 "dsymv.f"
	    }
#line 243 "dsymv.f"
	} else {
#line 244 "dsymv.f"
	    iy = ky;
#line 245 "dsymv.f"
	    if (*beta == 0.) {
#line 246 "dsymv.f"
		i__1 = *n;
#line 246 "dsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 247 "dsymv.f"
		    y[iy] = 0.;
#line 248 "dsymv.f"
		    iy += *incy;
#line 249 "dsymv.f"
/* L30: */
#line 249 "dsymv.f"
		}
#line 250 "dsymv.f"
	    } else {
#line 251 "dsymv.f"
		i__1 = *n;
#line 251 "dsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 252 "dsymv.f"
		    y[iy] = *beta * y[iy];
#line 253 "dsymv.f"
		    iy += *incy;
#line 254 "dsymv.f"
/* L40: */
#line 254 "dsymv.f"
		}
#line 255 "dsymv.f"
	    }
#line 256 "dsymv.f"
	}
#line 257 "dsymv.f"
    }
#line 258 "dsymv.f"
    if (*alpha == 0.) {
#line 258 "dsymv.f"
	return 0;
#line 258 "dsymv.f"
    }
#line 259 "dsymv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when A is stored in upper triangle. */

#line 263 "dsymv.f"
	if (*incx == 1 && *incy == 1) {
#line 264 "dsymv.f"
	    i__1 = *n;
#line 264 "dsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 265 "dsymv.f"
		temp1 = *alpha * x[j];
#line 266 "dsymv.f"
		temp2 = 0.;
#line 267 "dsymv.f"
		i__2 = j - 1;
#line 267 "dsymv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 268 "dsymv.f"
		    y[i__] += temp1 * a[i__ + j * a_dim1];
#line 269 "dsymv.f"
		    temp2 += a[i__ + j * a_dim1] * x[i__];
#line 270 "dsymv.f"
/* L50: */
#line 270 "dsymv.f"
		}
#line 271 "dsymv.f"
		y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
#line 272 "dsymv.f"
/* L60: */
#line 272 "dsymv.f"
	    }
#line 273 "dsymv.f"
	} else {
#line 274 "dsymv.f"
	    jx = kx;
#line 275 "dsymv.f"
	    jy = ky;
#line 276 "dsymv.f"
	    i__1 = *n;
#line 276 "dsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 277 "dsymv.f"
		temp1 = *alpha * x[jx];
#line 278 "dsymv.f"
		temp2 = 0.;
#line 279 "dsymv.f"
		ix = kx;
#line 280 "dsymv.f"
		iy = ky;
#line 281 "dsymv.f"
		i__2 = j - 1;
#line 281 "dsymv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 282 "dsymv.f"
		    y[iy] += temp1 * a[i__ + j * a_dim1];
#line 283 "dsymv.f"
		    temp2 += a[i__ + j * a_dim1] * x[ix];
#line 284 "dsymv.f"
		    ix += *incx;
#line 285 "dsymv.f"
		    iy += *incy;
#line 286 "dsymv.f"
/* L70: */
#line 286 "dsymv.f"
		}
#line 287 "dsymv.f"
		y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
#line 288 "dsymv.f"
		jx += *incx;
#line 289 "dsymv.f"
		jy += *incy;
#line 290 "dsymv.f"
/* L80: */
#line 290 "dsymv.f"
	    }
#line 291 "dsymv.f"
	}
#line 292 "dsymv.f"
    } else {

/*        Form  y  when A is stored in lower triangle. */

#line 296 "dsymv.f"
	if (*incx == 1 && *incy == 1) {
#line 297 "dsymv.f"
	    i__1 = *n;
#line 297 "dsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 298 "dsymv.f"
		temp1 = *alpha * x[j];
#line 299 "dsymv.f"
		temp2 = 0.;
#line 300 "dsymv.f"
		y[j] += temp1 * a[j + j * a_dim1];
#line 301 "dsymv.f"
		i__2 = *n;
#line 301 "dsymv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 302 "dsymv.f"
		    y[i__] += temp1 * a[i__ + j * a_dim1];
#line 303 "dsymv.f"
		    temp2 += a[i__ + j * a_dim1] * x[i__];
#line 304 "dsymv.f"
/* L90: */
#line 304 "dsymv.f"
		}
#line 305 "dsymv.f"
		y[j] += *alpha * temp2;
#line 306 "dsymv.f"
/* L100: */
#line 306 "dsymv.f"
	    }
#line 307 "dsymv.f"
	} else {
#line 308 "dsymv.f"
	    jx = kx;
#line 309 "dsymv.f"
	    jy = ky;
#line 310 "dsymv.f"
	    i__1 = *n;
#line 310 "dsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 311 "dsymv.f"
		temp1 = *alpha * x[jx];
#line 312 "dsymv.f"
		temp2 = 0.;
#line 313 "dsymv.f"
		y[jy] += temp1 * a[j + j * a_dim1];
#line 314 "dsymv.f"
		ix = jx;
#line 315 "dsymv.f"
		iy = jy;
#line 316 "dsymv.f"
		i__2 = *n;
#line 316 "dsymv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 317 "dsymv.f"
		    ix += *incx;
#line 318 "dsymv.f"
		    iy += *incy;
#line 319 "dsymv.f"
		    y[iy] += temp1 * a[i__ + j * a_dim1];
#line 320 "dsymv.f"
		    temp2 += a[i__ + j * a_dim1] * x[ix];
#line 321 "dsymv.f"
/* L110: */
#line 321 "dsymv.f"
		}
#line 322 "dsymv.f"
		y[jy] += *alpha * temp2;
#line 323 "dsymv.f"
		jx += *incx;
#line 324 "dsymv.f"
		jy += *incy;
#line 325 "dsymv.f"
/* L120: */
#line 325 "dsymv.f"
	    }
#line 326 "dsymv.f"
	}
#line 327 "dsymv.f"
    }

#line 329 "dsymv.f"
    return 0;

/*     End of DSYMV . */

} /* dsymv_ */

