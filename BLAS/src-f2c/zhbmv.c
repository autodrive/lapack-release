#line 1 "zhbmv.f"
/* zhbmv.f -- translated by f2c (version 20100827).
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

#line 1 "zhbmv.f"
/* > \brief \b ZHBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER INCX,INCY,K,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBMV  performs the matrix-vector  operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n hermitian band matrix, with k super-diagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the band matrix A is being supplied as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   The upper triangular part of A is */
/* >                                  being supplied. */
/* > */
/* >              UPLO = 'L' or 'l'   The lower triangular part of A is */
/* >                                  being supplied. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry, K specifies the number of super-diagonals of the */
/* >           matrix A. K must satisfy  0 .le. K. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
/* >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the upper triangular */
/* >           band part of the hermitian matrix, supplied column by */
/* >           column, with the leading diagonal of the matrix in row */
/* >           ( k + 1 ) of the array, the first super-diagonal starting at */
/* >           position 2 in row k, and so on. The top left k by k triangle */
/* >           of the array A is not referenced. */
/* >           The following program segment will transfer the upper */
/* >           triangular part of a hermitian band matrix from conventional */
/* >           full matrix storage to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = K + 1 - J */
/* >                    DO 10, I = MAX( 1, J - K ), J */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the lower triangular */
/* >           band part of the hermitian matrix, supplied column by */
/* >           column, with the leading diagonal of the matrix in row 1 of */
/* >           the array, the first sub-diagonal starting at position 1 in */
/* >           row 2, and so on. The bottom right k by k triangle of the */
/* >           array A is not referenced. */
/* >           The following program segment will transfer the lower */
/* >           triangular part of a hermitian band matrix from conventional */
/* >           full matrix storage to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = 1 - J */
/* >                    DO 10, I = J, MIN( N, J + K ) */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           ( k + 1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
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
/* >          BETA is COMPLEX*16 */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the */
/* >           vector y. On exit, Y is overwritten by the updated vector y. */
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

/* > \ingroup complex16_blas_level2 */

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
/* Subroutine */ int zhbmv_(char *uplo, integer *n, integer *k, doublecomplex 
	*alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *
	incx, doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, ix, iy, jx, jy, kx, ky, info;
    static doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
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

#line 229 "zhbmv.f"
    /* Parameter adjustments */
#line 229 "zhbmv.f"
    a_dim1 = *lda;
#line 229 "zhbmv.f"
    a_offset = 1 + a_dim1;
#line 229 "zhbmv.f"
    a -= a_offset;
#line 229 "zhbmv.f"
    --x;
#line 229 "zhbmv.f"
    --y;
#line 229 "zhbmv.f"

#line 229 "zhbmv.f"
    /* Function Body */
#line 229 "zhbmv.f"
    info = 0;
#line 230 "zhbmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 231 "zhbmv.f"
	info = 1;
#line 232 "zhbmv.f"
    } else if (*n < 0) {
#line 233 "zhbmv.f"
	info = 2;
#line 234 "zhbmv.f"
    } else if (*k < 0) {
#line 235 "zhbmv.f"
	info = 3;
#line 236 "zhbmv.f"
    } else if (*lda < *k + 1) {
#line 237 "zhbmv.f"
	info = 6;
#line 238 "zhbmv.f"
    } else if (*incx == 0) {
#line 239 "zhbmv.f"
	info = 8;
#line 240 "zhbmv.f"
    } else if (*incy == 0) {
#line 241 "zhbmv.f"
	info = 11;
#line 242 "zhbmv.f"
    }
#line 243 "zhbmv.f"
    if (info != 0) {
#line 244 "zhbmv.f"
	xerbla_("ZHBMV ", &info, (ftnlen)6);
#line 245 "zhbmv.f"
	return 0;
#line 246 "zhbmv.f"
    }

/*     Quick return if possible. */

#line 250 "zhbmv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 250 "zhbmv.f"
	return 0;
#line 250 "zhbmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 254 "zhbmv.f"
    if (*incx > 0) {
#line 255 "zhbmv.f"
	kx = 1;
#line 256 "zhbmv.f"
    } else {
#line 257 "zhbmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 258 "zhbmv.f"
    }
#line 259 "zhbmv.f"
    if (*incy > 0) {
#line 260 "zhbmv.f"
	ky = 1;
#line 261 "zhbmv.f"
    } else {
#line 262 "zhbmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 263 "zhbmv.f"
    }

/*     Start the operations. In this version the elements of the array A */
/*     are accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 270 "zhbmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 271 "zhbmv.f"
	if (*incy == 1) {
#line 272 "zhbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 273 "zhbmv.f"
		i__1 = *n;
#line 273 "zhbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "zhbmv.f"
		    i__2 = i__;
#line 274 "zhbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 275 "zhbmv.f"
/* L10: */
#line 275 "zhbmv.f"
		}
#line 276 "zhbmv.f"
	    } else {
#line 277 "zhbmv.f"
		i__1 = *n;
#line 277 "zhbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "zhbmv.f"
		    i__2 = i__;
#line 278 "zhbmv.f"
		    i__3 = i__;
#line 278 "zhbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 278 "zhbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 279 "zhbmv.f"
/* L20: */
#line 279 "zhbmv.f"
		}
#line 280 "zhbmv.f"
	    }
#line 281 "zhbmv.f"
	} else {
#line 282 "zhbmv.f"
	    iy = ky;
#line 283 "zhbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 284 "zhbmv.f"
		i__1 = *n;
#line 284 "zhbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "zhbmv.f"
		    i__2 = iy;
#line 285 "zhbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 286 "zhbmv.f"
		    iy += *incy;
#line 287 "zhbmv.f"
/* L30: */
#line 287 "zhbmv.f"
		}
#line 288 "zhbmv.f"
	    } else {
#line 289 "zhbmv.f"
		i__1 = *n;
#line 289 "zhbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "zhbmv.f"
		    i__2 = iy;
#line 290 "zhbmv.f"
		    i__3 = iy;
#line 290 "zhbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 290 "zhbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 291 "zhbmv.f"
		    iy += *incy;
#line 292 "zhbmv.f"
/* L40: */
#line 292 "zhbmv.f"
		}
#line 293 "zhbmv.f"
	    }
#line 294 "zhbmv.f"
	}
#line 295 "zhbmv.f"
    }
#line 296 "zhbmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 296 "zhbmv.f"
	return 0;
#line 296 "zhbmv.f"
    }
#line 297 "zhbmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when upper triangle of A is stored. */

#line 301 "zhbmv.f"
	kplus1 = *k + 1;
#line 302 "zhbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 303 "zhbmv.f"
	    i__1 = *n;
#line 303 "zhbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 304 "zhbmv.f"
		i__2 = j;
#line 304 "zhbmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 304 "zhbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 305 "zhbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 306 "zhbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 307 "zhbmv.f"
		i__2 = 1, i__3 = j - *k;
#line 307 "zhbmv.f"
		i__4 = j - 1;
#line 307 "zhbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 308 "zhbmv.f"
		    i__2 = i__;
#line 308 "zhbmv.f"
		    i__3 = i__;
#line 308 "zhbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 308 "zhbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 308 "zhbmv.f"
		    z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 308 "zhbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 309 "zhbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 309 "zhbmv.f"
		    i__2 = i__;
#line 309 "zhbmv.f"
		    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, z__2.i =
			     z__3.r * x[i__2].i + z__3.i * x[i__2].r;
#line 309 "zhbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 309 "zhbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 310 "zhbmv.f"
/* L50: */
#line 310 "zhbmv.f"
		}
#line 311 "zhbmv.f"
		i__4 = j;
#line 311 "zhbmv.f"
		i__2 = j;
#line 311 "zhbmv.f"
		i__3 = kplus1 + j * a_dim1;
#line 311 "zhbmv.f"
		d__1 = a[i__3].r;
#line 311 "zhbmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 311 "zhbmv.f"
		z__2.r = y[i__2].r + z__3.r, z__2.i = y[i__2].i + z__3.i;
#line 311 "zhbmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 311 "zhbmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 311 "zhbmv.f"
		y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 312 "zhbmv.f"
/* L60: */
#line 312 "zhbmv.f"
	    }
#line 313 "zhbmv.f"
	} else {
#line 314 "zhbmv.f"
	    jx = kx;
#line 315 "zhbmv.f"
	    jy = ky;
#line 316 "zhbmv.f"
	    i__1 = *n;
#line 316 "zhbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 317 "zhbmv.f"
		i__4 = jx;
#line 317 "zhbmv.f"
		z__1.r = alpha->r * x[i__4].r - alpha->i * x[i__4].i, z__1.i =
			 alpha->r * x[i__4].i + alpha->i * x[i__4].r;
#line 317 "zhbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 318 "zhbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 319 "zhbmv.f"
		ix = kx;
#line 320 "zhbmv.f"
		iy = ky;
#line 321 "zhbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 322 "zhbmv.f"
		i__4 = 1, i__2 = j - *k;
#line 322 "zhbmv.f"
		i__3 = j - 1;
#line 322 "zhbmv.f"
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 323 "zhbmv.f"
		    i__4 = iy;
#line 323 "zhbmv.f"
		    i__2 = iy;
#line 323 "zhbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 323 "zhbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 323 "zhbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 323 "zhbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 324 "zhbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 324 "zhbmv.f"
		    i__4 = ix;
#line 324 "zhbmv.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 324 "zhbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 324 "zhbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 325 "zhbmv.f"
		    ix += *incx;
#line 326 "zhbmv.f"
		    iy += *incy;
#line 327 "zhbmv.f"
/* L70: */
#line 327 "zhbmv.f"
		}
#line 328 "zhbmv.f"
		i__3 = jy;
#line 328 "zhbmv.f"
		i__4 = jy;
#line 328 "zhbmv.f"
		i__2 = kplus1 + j * a_dim1;
#line 328 "zhbmv.f"
		d__1 = a[i__2].r;
#line 328 "zhbmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 328 "zhbmv.f"
		z__2.r = y[i__4].r + z__3.r, z__2.i = y[i__4].i + z__3.i;
#line 328 "zhbmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 328 "zhbmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 328 "zhbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 329 "zhbmv.f"
		jx += *incx;
#line 330 "zhbmv.f"
		jy += *incy;
#line 331 "zhbmv.f"
		if (j > *k) {
#line 332 "zhbmv.f"
		    kx += *incx;
#line 333 "zhbmv.f"
		    ky += *incy;
#line 334 "zhbmv.f"
		}
#line 335 "zhbmv.f"
/* L80: */
#line 335 "zhbmv.f"
	    }
#line 336 "zhbmv.f"
	}
#line 337 "zhbmv.f"
    } else {

/*        Form  y  when lower triangle of A is stored. */

#line 341 "zhbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 342 "zhbmv.f"
	    i__1 = *n;
#line 342 "zhbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 343 "zhbmv.f"
		i__3 = j;
#line 343 "zhbmv.f"
		z__1.r = alpha->r * x[i__3].r - alpha->i * x[i__3].i, z__1.i =
			 alpha->r * x[i__3].i + alpha->i * x[i__3].r;
#line 343 "zhbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 344 "zhbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 345 "zhbmv.f"
		i__3 = j;
#line 345 "zhbmv.f"
		i__4 = j;
#line 345 "zhbmv.f"
		i__2 = j * a_dim1 + 1;
#line 345 "zhbmv.f"
		d__1 = a[i__2].r;
#line 345 "zhbmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 345 "zhbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 345 "zhbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 346 "zhbmv.f"
		l = 1 - j;
/* Computing MIN */
#line 347 "zhbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 347 "zhbmv.f"
		i__3 = min(i__4,i__2);
#line 347 "zhbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 348 "zhbmv.f"
		    i__4 = i__;
#line 348 "zhbmv.f"
		    i__2 = i__;
#line 348 "zhbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 348 "zhbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 348 "zhbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 348 "zhbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 349 "zhbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 349 "zhbmv.f"
		    i__4 = i__;
#line 349 "zhbmv.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 349 "zhbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 349 "zhbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 350 "zhbmv.f"
/* L90: */
#line 350 "zhbmv.f"
		}
#line 351 "zhbmv.f"
		i__3 = j;
#line 351 "zhbmv.f"
		i__4 = j;
#line 351 "zhbmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 351 "zhbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 351 "zhbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 352 "zhbmv.f"
/* L100: */
#line 352 "zhbmv.f"
	    }
#line 353 "zhbmv.f"
	} else {
#line 354 "zhbmv.f"
	    jx = kx;
#line 355 "zhbmv.f"
	    jy = ky;
#line 356 "zhbmv.f"
	    i__1 = *n;
#line 356 "zhbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 357 "zhbmv.f"
		i__3 = jx;
#line 357 "zhbmv.f"
		z__1.r = alpha->r * x[i__3].r - alpha->i * x[i__3].i, z__1.i =
			 alpha->r * x[i__3].i + alpha->i * x[i__3].r;
#line 357 "zhbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 358 "zhbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 359 "zhbmv.f"
		i__3 = jy;
#line 359 "zhbmv.f"
		i__4 = jy;
#line 359 "zhbmv.f"
		i__2 = j * a_dim1 + 1;
#line 359 "zhbmv.f"
		d__1 = a[i__2].r;
#line 359 "zhbmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 359 "zhbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 359 "zhbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 360 "zhbmv.f"
		l = 1 - j;
#line 361 "zhbmv.f"
		ix = jx;
#line 362 "zhbmv.f"
		iy = jy;
/* Computing MIN */
#line 363 "zhbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 363 "zhbmv.f"
		i__3 = min(i__4,i__2);
#line 363 "zhbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 364 "zhbmv.f"
		    ix += *incx;
#line 365 "zhbmv.f"
		    iy += *incy;
#line 366 "zhbmv.f"
		    i__4 = iy;
#line 366 "zhbmv.f"
		    i__2 = iy;
#line 366 "zhbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 366 "zhbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 366 "zhbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 366 "zhbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 367 "zhbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 367 "zhbmv.f"
		    i__4 = ix;
#line 367 "zhbmv.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 367 "zhbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 367 "zhbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 368 "zhbmv.f"
/* L110: */
#line 368 "zhbmv.f"
		}
#line 369 "zhbmv.f"
		i__3 = jy;
#line 369 "zhbmv.f"
		i__4 = jy;
#line 369 "zhbmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 369 "zhbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 369 "zhbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 370 "zhbmv.f"
		jx += *incx;
#line 371 "zhbmv.f"
		jy += *incy;
#line 372 "zhbmv.f"
/* L120: */
#line 372 "zhbmv.f"
	    }
#line 373 "zhbmv.f"
	}
#line 374 "zhbmv.f"
    }

#line 376 "zhbmv.f"
    return 0;

/*     End of ZHBMV . */

} /* zhbmv_ */

