#line 1 "chbmv.f"
/* chbmv.f -- translated by f2c (version 20100827).
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

#line 1 "chbmv.f"
/* > \brief \b CHBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER INCX,INCY,K,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBMV  performs the matrix-vector  operation */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension ( LDA, N ) */
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
/* >          X is COMPLEX array, dimension at least */
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
/* >          BETA is COMPLEX */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension at least */
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

/* > \ingroup complex_blas_level2 */

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
/* Subroutine */ int chbmv_(char *uplo, integer *n, integer *k, doublecomplex 
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

#line 229 "chbmv.f"
    /* Parameter adjustments */
#line 229 "chbmv.f"
    a_dim1 = *lda;
#line 229 "chbmv.f"
    a_offset = 1 + a_dim1;
#line 229 "chbmv.f"
    a -= a_offset;
#line 229 "chbmv.f"
    --x;
#line 229 "chbmv.f"
    --y;
#line 229 "chbmv.f"

#line 229 "chbmv.f"
    /* Function Body */
#line 229 "chbmv.f"
    info = 0;
#line 230 "chbmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 231 "chbmv.f"
	info = 1;
#line 232 "chbmv.f"
    } else if (*n < 0) {
#line 233 "chbmv.f"
	info = 2;
#line 234 "chbmv.f"
    } else if (*k < 0) {
#line 235 "chbmv.f"
	info = 3;
#line 236 "chbmv.f"
    } else if (*lda < *k + 1) {
#line 237 "chbmv.f"
	info = 6;
#line 238 "chbmv.f"
    } else if (*incx == 0) {
#line 239 "chbmv.f"
	info = 8;
#line 240 "chbmv.f"
    } else if (*incy == 0) {
#line 241 "chbmv.f"
	info = 11;
#line 242 "chbmv.f"
    }
#line 243 "chbmv.f"
    if (info != 0) {
#line 244 "chbmv.f"
	xerbla_("CHBMV ", &info, (ftnlen)6);
#line 245 "chbmv.f"
	return 0;
#line 246 "chbmv.f"
    }

/*     Quick return if possible. */

#line 250 "chbmv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 250 "chbmv.f"
	return 0;
#line 250 "chbmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 254 "chbmv.f"
    if (*incx > 0) {
#line 255 "chbmv.f"
	kx = 1;
#line 256 "chbmv.f"
    } else {
#line 257 "chbmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 258 "chbmv.f"
    }
#line 259 "chbmv.f"
    if (*incy > 0) {
#line 260 "chbmv.f"
	ky = 1;
#line 261 "chbmv.f"
    } else {
#line 262 "chbmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 263 "chbmv.f"
    }

/*     Start the operations. In this version the elements of the array A */
/*     are accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 270 "chbmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 271 "chbmv.f"
	if (*incy == 1) {
#line 272 "chbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 273 "chbmv.f"
		i__1 = *n;
#line 273 "chbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "chbmv.f"
		    i__2 = i__;
#line 274 "chbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 275 "chbmv.f"
/* L10: */
#line 275 "chbmv.f"
		}
#line 276 "chbmv.f"
	    } else {
#line 277 "chbmv.f"
		i__1 = *n;
#line 277 "chbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "chbmv.f"
		    i__2 = i__;
#line 278 "chbmv.f"
		    i__3 = i__;
#line 278 "chbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 278 "chbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 279 "chbmv.f"
/* L20: */
#line 279 "chbmv.f"
		}
#line 280 "chbmv.f"
	    }
#line 281 "chbmv.f"
	} else {
#line 282 "chbmv.f"
	    iy = ky;
#line 283 "chbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 284 "chbmv.f"
		i__1 = *n;
#line 284 "chbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "chbmv.f"
		    i__2 = iy;
#line 285 "chbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 286 "chbmv.f"
		    iy += *incy;
#line 287 "chbmv.f"
/* L30: */
#line 287 "chbmv.f"
		}
#line 288 "chbmv.f"
	    } else {
#line 289 "chbmv.f"
		i__1 = *n;
#line 289 "chbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "chbmv.f"
		    i__2 = iy;
#line 290 "chbmv.f"
		    i__3 = iy;
#line 290 "chbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 290 "chbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 291 "chbmv.f"
		    iy += *incy;
#line 292 "chbmv.f"
/* L40: */
#line 292 "chbmv.f"
		}
#line 293 "chbmv.f"
	    }
#line 294 "chbmv.f"
	}
#line 295 "chbmv.f"
    }
#line 296 "chbmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 296 "chbmv.f"
	return 0;
#line 296 "chbmv.f"
    }
#line 297 "chbmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when upper triangle of A is stored. */

#line 301 "chbmv.f"
	kplus1 = *k + 1;
#line 302 "chbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 303 "chbmv.f"
	    i__1 = *n;
#line 303 "chbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 304 "chbmv.f"
		i__2 = j;
#line 304 "chbmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 304 "chbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 305 "chbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 306 "chbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 307 "chbmv.f"
		i__2 = 1, i__3 = j - *k;
#line 307 "chbmv.f"
		i__4 = j - 1;
#line 307 "chbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 308 "chbmv.f"
		    i__2 = i__;
#line 308 "chbmv.f"
		    i__3 = i__;
#line 308 "chbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 308 "chbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 308 "chbmv.f"
		    z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 308 "chbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 309 "chbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 309 "chbmv.f"
		    i__2 = i__;
#line 309 "chbmv.f"
		    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, z__2.i =
			     z__3.r * x[i__2].i + z__3.i * x[i__2].r;
#line 309 "chbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 309 "chbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 310 "chbmv.f"
/* L50: */
#line 310 "chbmv.f"
		}
#line 311 "chbmv.f"
		i__4 = j;
#line 311 "chbmv.f"
		i__2 = j;
#line 311 "chbmv.f"
		i__3 = kplus1 + j * a_dim1;
#line 311 "chbmv.f"
		d__1 = a[i__3].r;
#line 311 "chbmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 311 "chbmv.f"
		z__2.r = y[i__2].r + z__3.r, z__2.i = y[i__2].i + z__3.i;
#line 311 "chbmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 311 "chbmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 311 "chbmv.f"
		y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 312 "chbmv.f"
/* L60: */
#line 312 "chbmv.f"
	    }
#line 313 "chbmv.f"
	} else {
#line 314 "chbmv.f"
	    jx = kx;
#line 315 "chbmv.f"
	    jy = ky;
#line 316 "chbmv.f"
	    i__1 = *n;
#line 316 "chbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 317 "chbmv.f"
		i__4 = jx;
#line 317 "chbmv.f"
		z__1.r = alpha->r * x[i__4].r - alpha->i * x[i__4].i, z__1.i =
			 alpha->r * x[i__4].i + alpha->i * x[i__4].r;
#line 317 "chbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 318 "chbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 319 "chbmv.f"
		ix = kx;
#line 320 "chbmv.f"
		iy = ky;
#line 321 "chbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 322 "chbmv.f"
		i__4 = 1, i__2 = j - *k;
#line 322 "chbmv.f"
		i__3 = j - 1;
#line 322 "chbmv.f"
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 323 "chbmv.f"
		    i__4 = iy;
#line 323 "chbmv.f"
		    i__2 = iy;
#line 323 "chbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 323 "chbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 323 "chbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 323 "chbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 324 "chbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 324 "chbmv.f"
		    i__4 = ix;
#line 324 "chbmv.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 324 "chbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 324 "chbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 325 "chbmv.f"
		    ix += *incx;
#line 326 "chbmv.f"
		    iy += *incy;
#line 327 "chbmv.f"
/* L70: */
#line 327 "chbmv.f"
		}
#line 328 "chbmv.f"
		i__3 = jy;
#line 328 "chbmv.f"
		i__4 = jy;
#line 328 "chbmv.f"
		i__2 = kplus1 + j * a_dim1;
#line 328 "chbmv.f"
		d__1 = a[i__2].r;
#line 328 "chbmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 328 "chbmv.f"
		z__2.r = y[i__4].r + z__3.r, z__2.i = y[i__4].i + z__3.i;
#line 328 "chbmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 328 "chbmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 328 "chbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 329 "chbmv.f"
		jx += *incx;
#line 330 "chbmv.f"
		jy += *incy;
#line 331 "chbmv.f"
		if (j > *k) {
#line 332 "chbmv.f"
		    kx += *incx;
#line 333 "chbmv.f"
		    ky += *incy;
#line 334 "chbmv.f"
		}
#line 335 "chbmv.f"
/* L80: */
#line 335 "chbmv.f"
	    }
#line 336 "chbmv.f"
	}
#line 337 "chbmv.f"
    } else {

/*        Form  y  when lower triangle of A is stored. */

#line 341 "chbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 342 "chbmv.f"
	    i__1 = *n;
#line 342 "chbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 343 "chbmv.f"
		i__3 = j;
#line 343 "chbmv.f"
		z__1.r = alpha->r * x[i__3].r - alpha->i * x[i__3].i, z__1.i =
			 alpha->r * x[i__3].i + alpha->i * x[i__3].r;
#line 343 "chbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 344 "chbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 345 "chbmv.f"
		i__3 = j;
#line 345 "chbmv.f"
		i__4 = j;
#line 345 "chbmv.f"
		i__2 = j * a_dim1 + 1;
#line 345 "chbmv.f"
		d__1 = a[i__2].r;
#line 345 "chbmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 345 "chbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 345 "chbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 346 "chbmv.f"
		l = 1 - j;
/* Computing MIN */
#line 347 "chbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 347 "chbmv.f"
		i__3 = min(i__4,i__2);
#line 347 "chbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 348 "chbmv.f"
		    i__4 = i__;
#line 348 "chbmv.f"
		    i__2 = i__;
#line 348 "chbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 348 "chbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 348 "chbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 348 "chbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 349 "chbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 349 "chbmv.f"
		    i__4 = i__;
#line 349 "chbmv.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 349 "chbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 349 "chbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 350 "chbmv.f"
/* L90: */
#line 350 "chbmv.f"
		}
#line 351 "chbmv.f"
		i__3 = j;
#line 351 "chbmv.f"
		i__4 = j;
#line 351 "chbmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 351 "chbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 351 "chbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 352 "chbmv.f"
/* L100: */
#line 352 "chbmv.f"
	    }
#line 353 "chbmv.f"
	} else {
#line 354 "chbmv.f"
	    jx = kx;
#line 355 "chbmv.f"
	    jy = ky;
#line 356 "chbmv.f"
	    i__1 = *n;
#line 356 "chbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 357 "chbmv.f"
		i__3 = jx;
#line 357 "chbmv.f"
		z__1.r = alpha->r * x[i__3].r - alpha->i * x[i__3].i, z__1.i =
			 alpha->r * x[i__3].i + alpha->i * x[i__3].r;
#line 357 "chbmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 358 "chbmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 359 "chbmv.f"
		i__3 = jy;
#line 359 "chbmv.f"
		i__4 = jy;
#line 359 "chbmv.f"
		i__2 = j * a_dim1 + 1;
#line 359 "chbmv.f"
		d__1 = a[i__2].r;
#line 359 "chbmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 359 "chbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 359 "chbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 360 "chbmv.f"
		l = 1 - j;
#line 361 "chbmv.f"
		ix = jx;
#line 362 "chbmv.f"
		iy = jy;
/* Computing MIN */
#line 363 "chbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 363 "chbmv.f"
		i__3 = min(i__4,i__2);
#line 363 "chbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 364 "chbmv.f"
		    ix += *incx;
#line 365 "chbmv.f"
		    iy += *incy;
#line 366 "chbmv.f"
		    i__4 = iy;
#line 366 "chbmv.f"
		    i__2 = iy;
#line 366 "chbmv.f"
		    i__5 = l + i__ + j * a_dim1;
#line 366 "chbmv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 366 "chbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 366 "chbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 367 "chbmv.f"
		    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 367 "chbmv.f"
		    i__4 = ix;
#line 367 "chbmv.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 367 "chbmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 367 "chbmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 368 "chbmv.f"
/* L110: */
#line 368 "chbmv.f"
		}
#line 369 "chbmv.f"
		i__3 = jy;
#line 369 "chbmv.f"
		i__4 = jy;
#line 369 "chbmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 369 "chbmv.f"
		z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 369 "chbmv.f"
		y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 370 "chbmv.f"
		jx += *incx;
#line 371 "chbmv.f"
		jy += *incy;
#line 372 "chbmv.f"
/* L120: */
#line 372 "chbmv.f"
	    }
#line 373 "chbmv.f"
	}
#line 374 "chbmv.f"
    }

#line 376 "chbmv.f"
    return 0;

/*     End of CHBMV . */

} /* chbmv_ */

