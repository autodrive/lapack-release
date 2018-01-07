#line 1 "dsyr2.f"
/* dsyr2.f -- translated by f2c (version 20100827).
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

#line 1 "dsyr2.f"
/* > \brief \b DSYR2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
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
/* > DSYR2  performs the symmetric rank 2 operation */
/* > */
/* >    A := alpha*x*y**T + alpha*y*x**T + A, */
/* > */
/* > where alpha is a scalar, x and y are n element vectors and A is an n */
/* > by n symmetric matrix. */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array of dimension at least */
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
/* >          Y is DOUBLE PRECISION array of dimension at least */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. On exit, the */
/* >           upper triangular part of the array A is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           upper triangular part of A is not referenced. On exit, the */
/* >           lower triangular part of the array A is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* Subroutine */ int dsyr2_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
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
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

#line 187 "dsyr2.f"
    /* Parameter adjustments */
#line 187 "dsyr2.f"
    --x;
#line 187 "dsyr2.f"
    --y;
#line 187 "dsyr2.f"
    a_dim1 = *lda;
#line 187 "dsyr2.f"
    a_offset = 1 + a_dim1;
#line 187 "dsyr2.f"
    a -= a_offset;
#line 187 "dsyr2.f"

#line 187 "dsyr2.f"
    /* Function Body */
#line 187 "dsyr2.f"
    info = 0;
#line 188 "dsyr2.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 189 "dsyr2.f"
	info = 1;
#line 190 "dsyr2.f"
    } else if (*n < 0) {
#line 191 "dsyr2.f"
	info = 2;
#line 192 "dsyr2.f"
    } else if (*incx == 0) {
#line 193 "dsyr2.f"
	info = 5;
#line 194 "dsyr2.f"
    } else if (*incy == 0) {
#line 195 "dsyr2.f"
	info = 7;
#line 196 "dsyr2.f"
    } else if (*lda < max(1,*n)) {
#line 197 "dsyr2.f"
	info = 9;
#line 198 "dsyr2.f"
    }
#line 199 "dsyr2.f"
    if (info != 0) {
#line 200 "dsyr2.f"
	xerbla_("DSYR2 ", &info, (ftnlen)6);
#line 201 "dsyr2.f"
	return 0;
#line 202 "dsyr2.f"
    }

/*     Quick return if possible. */

#line 206 "dsyr2.f"
    if (*n == 0 || *alpha == 0.) {
#line 206 "dsyr2.f"
	return 0;
#line 206 "dsyr2.f"
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

#line 211 "dsyr2.f"
    if (*incx != 1 || *incy != 1) {
#line 212 "dsyr2.f"
	if (*incx > 0) {
#line 213 "dsyr2.f"
	    kx = 1;
#line 214 "dsyr2.f"
	} else {
#line 215 "dsyr2.f"
	    kx = 1 - (*n - 1) * *incx;
#line 216 "dsyr2.f"
	}
#line 217 "dsyr2.f"
	if (*incy > 0) {
#line 218 "dsyr2.f"
	    ky = 1;
#line 219 "dsyr2.f"
	} else {
#line 220 "dsyr2.f"
	    ky = 1 - (*n - 1) * *incy;
#line 221 "dsyr2.f"
	}
#line 222 "dsyr2.f"
	jx = kx;
#line 223 "dsyr2.f"
	jy = ky;
#line 224 "dsyr2.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 230 "dsyr2.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in the upper triangle. */

#line 234 "dsyr2.f"
	if (*incx == 1 && *incy == 1) {
#line 235 "dsyr2.f"
	    i__1 = *n;
#line 235 "dsyr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 236 "dsyr2.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 237 "dsyr2.f"
		    temp1 = *alpha * y[j];
#line 238 "dsyr2.f"
		    temp2 = *alpha * x[j];
#line 239 "dsyr2.f"
		    i__2 = j;
#line 239 "dsyr2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 240 "dsyr2.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
#line 241 "dsyr2.f"
/* L10: */
#line 241 "dsyr2.f"
		    }
#line 242 "dsyr2.f"
		}
#line 243 "dsyr2.f"
/* L20: */
#line 243 "dsyr2.f"
	    }
#line 244 "dsyr2.f"
	} else {
#line 245 "dsyr2.f"
	    i__1 = *n;
#line 245 "dsyr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 246 "dsyr2.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 247 "dsyr2.f"
		    temp1 = *alpha * y[jy];
#line 248 "dsyr2.f"
		    temp2 = *alpha * x[jx];
#line 249 "dsyr2.f"
		    ix = kx;
#line 250 "dsyr2.f"
		    iy = ky;
#line 251 "dsyr2.f"
		    i__2 = j;
#line 251 "dsyr2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 252 "dsyr2.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
#line 253 "dsyr2.f"
			ix += *incx;
#line 254 "dsyr2.f"
			iy += *incy;
#line 255 "dsyr2.f"
/* L30: */
#line 255 "dsyr2.f"
		    }
#line 256 "dsyr2.f"
		}
#line 257 "dsyr2.f"
		jx += *incx;
#line 258 "dsyr2.f"
		jy += *incy;
#line 259 "dsyr2.f"
/* L40: */
#line 259 "dsyr2.f"
	    }
#line 260 "dsyr2.f"
	}
#line 261 "dsyr2.f"
    } else {

/*        Form  A  when A is stored in the lower triangle. */

#line 265 "dsyr2.f"
	if (*incx == 1 && *incy == 1) {
#line 266 "dsyr2.f"
	    i__1 = *n;
#line 266 "dsyr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 267 "dsyr2.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 268 "dsyr2.f"
		    temp1 = *alpha * y[j];
#line 269 "dsyr2.f"
		    temp2 = *alpha * x[j];
#line 270 "dsyr2.f"
		    i__2 = *n;
#line 270 "dsyr2.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 271 "dsyr2.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
#line 272 "dsyr2.f"
/* L50: */
#line 272 "dsyr2.f"
		    }
#line 273 "dsyr2.f"
		}
#line 274 "dsyr2.f"
/* L60: */
#line 274 "dsyr2.f"
	    }
#line 275 "dsyr2.f"
	} else {
#line 276 "dsyr2.f"
	    i__1 = *n;
#line 276 "dsyr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 277 "dsyr2.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 278 "dsyr2.f"
		    temp1 = *alpha * y[jy];
#line 279 "dsyr2.f"
		    temp2 = *alpha * x[jx];
#line 280 "dsyr2.f"
		    ix = jx;
#line 281 "dsyr2.f"
		    iy = jy;
#line 282 "dsyr2.f"
		    i__2 = *n;
#line 282 "dsyr2.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 283 "dsyr2.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
#line 284 "dsyr2.f"
			ix += *incx;
#line 285 "dsyr2.f"
			iy += *incy;
#line 286 "dsyr2.f"
/* L70: */
#line 286 "dsyr2.f"
		    }
#line 287 "dsyr2.f"
		}
#line 288 "dsyr2.f"
		jx += *incx;
#line 289 "dsyr2.f"
		jy += *incy;
#line 290 "dsyr2.f"
/* L80: */
#line 290 "dsyr2.f"
	    }
#line 291 "dsyr2.f"
	}
#line 292 "dsyr2.f"
    }

#line 294 "dsyr2.f"
    return 0;

/*     End of DSYR2 . */

} /* dsyr2_ */

