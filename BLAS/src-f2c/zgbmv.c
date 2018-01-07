#line 1 "zgbmv.f"
/* zgbmv.f -- translated by f2c (version 20100827).
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

#line 1 "zgbmv.f"
/* > \brief \b ZGBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER INCX,INCY,KL,KU,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or */
/* > */
/* >    y := alpha*A**H*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n band matrix, with kl sub-diagonals and ku super-diagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */
/* > */
/* >              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */
/* > */
/* >              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >           On entry, KL specifies the number of sub-diagonals of the */
/* >           matrix A. KL must satisfy  0 .le. KL. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >           On entry, KU specifies the number of super-diagonals of the */
/* >           matrix A. KU must satisfy  0 .le. KU. */
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
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ). */
/* >           Before entry, the leading ( kl + ku + 1 ) by n part of the */
/* >           array A must contain the matrix of coefficients, supplied */
/* >           column by column, with the leading diagonal of the matrix in */
/* >           row ( ku + 1 ) of the array, the first super-diagonal */
/* >           starting at position 2 in row ku, the first sub-diagonal */
/* >           starting at position 1 in row ( ku + 2 ), and so on. */
/* >           Elements in the array A that do not correspond to elements */
/* >           in the band matrix (such as the top left ku by ku triangle) */
/* >           are not referenced. */
/* >           The following program segment will transfer a band matrix */
/* >           from conventional full matrix storage to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    K = KU + 1 - J */
/* >                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL ) */
/* >                       A( K + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           ( kl + ku + 1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of DIMENSION at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
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
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array of DIMENSION at least */
/* >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
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

/* > \date November 2011 */

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
/* Subroutine */ int zgbmv_(char *trans, integer *m, integer *n, integer *kl, 
	integer *ku, doublecomplex *alpha, doublecomplex *a, integer *lda, 
	doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *
	y, integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
    static doublecomplex temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconj;


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

#line 230 "zgbmv.f"
    /* Parameter adjustments */
#line 230 "zgbmv.f"
    a_dim1 = *lda;
#line 230 "zgbmv.f"
    a_offset = 1 + a_dim1;
#line 230 "zgbmv.f"
    a -= a_offset;
#line 230 "zgbmv.f"
    --x;
#line 230 "zgbmv.f"
    --y;
#line 230 "zgbmv.f"

#line 230 "zgbmv.f"
    /* Function Body */
#line 230 "zgbmv.f"
    info = 0;
#line 231 "zgbmv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 233 "zgbmv.f"
	info = 1;
#line 234 "zgbmv.f"
    } else if (*m < 0) {
#line 235 "zgbmv.f"
	info = 2;
#line 236 "zgbmv.f"
    } else if (*n < 0) {
#line 237 "zgbmv.f"
	info = 3;
#line 238 "zgbmv.f"
    } else if (*kl < 0) {
#line 239 "zgbmv.f"
	info = 4;
#line 240 "zgbmv.f"
    } else if (*ku < 0) {
#line 241 "zgbmv.f"
	info = 5;
#line 242 "zgbmv.f"
    } else if (*lda < *kl + *ku + 1) {
#line 243 "zgbmv.f"
	info = 8;
#line 244 "zgbmv.f"
    } else if (*incx == 0) {
#line 245 "zgbmv.f"
	info = 10;
#line 246 "zgbmv.f"
    } else if (*incy == 0) {
#line 247 "zgbmv.f"
	info = 13;
#line 248 "zgbmv.f"
    }
#line 249 "zgbmv.f"
    if (info != 0) {
#line 250 "zgbmv.f"
	xerbla_("ZGBMV ", &info, (ftnlen)6);
#line 251 "zgbmv.f"
	return 0;
#line 252 "zgbmv.f"
    }

/*     Quick return if possible. */

#line 256 "zgbmv.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
#line 256 "zgbmv.f"
	return 0;
#line 256 "zgbmv.f"
    }

#line 259 "zgbmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 264 "zgbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 265 "zgbmv.f"
	lenx = *n;
#line 266 "zgbmv.f"
	leny = *m;
#line 267 "zgbmv.f"
    } else {
#line 268 "zgbmv.f"
	lenx = *m;
#line 269 "zgbmv.f"
	leny = *n;
#line 270 "zgbmv.f"
    }
#line 271 "zgbmv.f"
    if (*incx > 0) {
#line 272 "zgbmv.f"
	kx = 1;
#line 273 "zgbmv.f"
    } else {
#line 274 "zgbmv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 275 "zgbmv.f"
    }
#line 276 "zgbmv.f"
    if (*incy > 0) {
#line 277 "zgbmv.f"
	ky = 1;
#line 278 "zgbmv.f"
    } else {
#line 279 "zgbmv.f"
	ky = 1 - (leny - 1) * *incy;
#line 280 "zgbmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the band part of A. */

/*     First form  y := beta*y. */

#line 287 "zgbmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 288 "zgbmv.f"
	if (*incy == 1) {
#line 289 "zgbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 290 "zgbmv.f"
		i__1 = leny;
#line 290 "zgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "zgbmv.f"
		    i__2 = i__;
#line 291 "zgbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 292 "zgbmv.f"
/* L10: */
#line 292 "zgbmv.f"
		}
#line 293 "zgbmv.f"
	    } else {
#line 294 "zgbmv.f"
		i__1 = leny;
#line 294 "zgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 295 "zgbmv.f"
		    i__2 = i__;
#line 295 "zgbmv.f"
		    i__3 = i__;
#line 295 "zgbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 295 "zgbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 296 "zgbmv.f"
/* L20: */
#line 296 "zgbmv.f"
		}
#line 297 "zgbmv.f"
	    }
#line 298 "zgbmv.f"
	} else {
#line 299 "zgbmv.f"
	    iy = ky;
#line 300 "zgbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 301 "zgbmv.f"
		i__1 = leny;
#line 301 "zgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "zgbmv.f"
		    i__2 = iy;
#line 302 "zgbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 303 "zgbmv.f"
		    iy += *incy;
#line 304 "zgbmv.f"
/* L30: */
#line 304 "zgbmv.f"
		}
#line 305 "zgbmv.f"
	    } else {
#line 306 "zgbmv.f"
		i__1 = leny;
#line 306 "zgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 307 "zgbmv.f"
		    i__2 = iy;
#line 307 "zgbmv.f"
		    i__3 = iy;
#line 307 "zgbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 307 "zgbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 308 "zgbmv.f"
		    iy += *incy;
#line 309 "zgbmv.f"
/* L40: */
#line 309 "zgbmv.f"
		}
#line 310 "zgbmv.f"
	    }
#line 311 "zgbmv.f"
	}
#line 312 "zgbmv.f"
    }
#line 313 "zgbmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 313 "zgbmv.f"
	return 0;
#line 313 "zgbmv.f"
    }
#line 314 "zgbmv.f"
    kup1 = *ku + 1;
#line 315 "zgbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 319 "zgbmv.f"
	jx = kx;
#line 320 "zgbmv.f"
	if (*incy == 1) {
#line 321 "zgbmv.f"
	    i__1 = *n;
#line 321 "zgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 322 "zgbmv.f"
		i__2 = jx;
#line 322 "zgbmv.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 323 "zgbmv.f"
		    i__2 = jx;
#line 323 "zgbmv.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 323 "zgbmv.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 324 "zgbmv.f"
		    k = kup1 - j;
/* Computing MAX */
#line 325 "zgbmv.f"
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 325 "zgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 325 "zgbmv.f"
		    i__4 = min(i__5,i__6);
#line 325 "zgbmv.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 326 "zgbmv.f"
			i__2 = i__;
#line 326 "zgbmv.f"
			i__3 = i__;
#line 326 "zgbmv.f"
			i__5 = k + i__ + j * a_dim1;
#line 326 "zgbmv.f"
			z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				z__2.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
#line 326 "zgbmv.f"
			z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + 
				z__2.i;
#line 326 "zgbmv.f"
			y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 327 "zgbmv.f"
/* L50: */
#line 327 "zgbmv.f"
		    }
#line 328 "zgbmv.f"
		}
#line 329 "zgbmv.f"
		jx += *incx;
#line 330 "zgbmv.f"
/* L60: */
#line 330 "zgbmv.f"
	    }
#line 331 "zgbmv.f"
	} else {
#line 332 "zgbmv.f"
	    i__1 = *n;
#line 332 "zgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 333 "zgbmv.f"
		i__4 = jx;
#line 333 "zgbmv.f"
		if (x[i__4].r != 0. || x[i__4].i != 0.) {
#line 334 "zgbmv.f"
		    i__4 = jx;
#line 334 "zgbmv.f"
		    z__1.r = alpha->r * x[i__4].r - alpha->i * x[i__4].i, 
			    z__1.i = alpha->r * x[i__4].i + alpha->i * x[i__4]
			    .r;
#line 334 "zgbmv.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 335 "zgbmv.f"
		    iy = ky;
#line 336 "zgbmv.f"
		    k = kup1 - j;
/* Computing MAX */
#line 337 "zgbmv.f"
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 337 "zgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 337 "zgbmv.f"
		    i__3 = min(i__5,i__6);
#line 337 "zgbmv.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 338 "zgbmv.f"
			i__4 = iy;
#line 338 "zgbmv.f"
			i__2 = iy;
#line 338 "zgbmv.f"
			i__5 = k + i__ + j * a_dim1;
#line 338 "zgbmv.f"
			z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				z__2.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
#line 338 "zgbmv.f"
			z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + 
				z__2.i;
#line 338 "zgbmv.f"
			y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 339 "zgbmv.f"
			iy += *incy;
#line 340 "zgbmv.f"
/* L70: */
#line 340 "zgbmv.f"
		    }
#line 341 "zgbmv.f"
		}
#line 342 "zgbmv.f"
		jx += *incx;
#line 343 "zgbmv.f"
		if (j > *ku) {
#line 343 "zgbmv.f"
		    ky += *incy;
#line 343 "zgbmv.f"
		}
#line 344 "zgbmv.f"
/* L80: */
#line 344 "zgbmv.f"
	    }
#line 345 "zgbmv.f"
	}
#line 346 "zgbmv.f"
    } else {

/*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y. */

#line 350 "zgbmv.f"
	jy = ky;
#line 351 "zgbmv.f"
	if (*incx == 1) {
#line 352 "zgbmv.f"
	    i__1 = *n;
#line 352 "zgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 353 "zgbmv.f"
		temp.r = 0., temp.i = 0.;
#line 354 "zgbmv.f"
		k = kup1 - j;
#line 355 "zgbmv.f"
		if (noconj) {
/* Computing MAX */
#line 356 "zgbmv.f"
		    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 356 "zgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 356 "zgbmv.f"
		    i__2 = min(i__5,i__6);
#line 356 "zgbmv.f"
		    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 357 "zgbmv.f"
			i__3 = k + i__ + j * a_dim1;
#line 357 "zgbmv.f"
			i__4 = i__;
#line 357 "zgbmv.f"
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
#line 357 "zgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 357 "zgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 358 "zgbmv.f"
/* L90: */
#line 358 "zgbmv.f"
		    }
#line 359 "zgbmv.f"
		} else {
/* Computing MAX */
#line 360 "zgbmv.f"
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 360 "zgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 360 "zgbmv.f"
		    i__4 = min(i__5,i__6);
#line 360 "zgbmv.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 361 "zgbmv.f"
			d_cnjg(&z__3, &a[k + i__ + j * a_dim1]);
#line 361 "zgbmv.f"
			i__2 = i__;
#line 361 "zgbmv.f"
			z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				z__2.i = z__3.r * x[i__2].i + z__3.i * x[i__2]
				.r;
#line 361 "zgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 361 "zgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 362 "zgbmv.f"
/* L100: */
#line 362 "zgbmv.f"
		    }
#line 363 "zgbmv.f"
		}
#line 364 "zgbmv.f"
		i__4 = jy;
#line 364 "zgbmv.f"
		i__2 = jy;
#line 364 "zgbmv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 364 "zgbmv.f"
		z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 364 "zgbmv.f"
		y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 365 "zgbmv.f"
		jy += *incy;
#line 366 "zgbmv.f"
/* L110: */
#line 366 "zgbmv.f"
	    }
#line 367 "zgbmv.f"
	} else {
#line 368 "zgbmv.f"
	    i__1 = *n;
#line 368 "zgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 369 "zgbmv.f"
		temp.r = 0., temp.i = 0.;
#line 370 "zgbmv.f"
		ix = kx;
#line 371 "zgbmv.f"
		k = kup1 - j;
#line 372 "zgbmv.f"
		if (noconj) {
/* Computing MAX */
#line 373 "zgbmv.f"
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 373 "zgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 373 "zgbmv.f"
		    i__3 = min(i__5,i__6);
#line 373 "zgbmv.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 374 "zgbmv.f"
			i__4 = k + i__ + j * a_dim1;
#line 374 "zgbmv.f"
			i__2 = ix;
#line 374 "zgbmv.f"
			z__2.r = a[i__4].r * x[i__2].r - a[i__4].i * x[i__2]
				.i, z__2.i = a[i__4].r * x[i__2].i + a[i__4]
				.i * x[i__2].r;
#line 374 "zgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 374 "zgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 375 "zgbmv.f"
			ix += *incx;
#line 376 "zgbmv.f"
/* L120: */
#line 376 "zgbmv.f"
		    }
#line 377 "zgbmv.f"
		} else {
/* Computing MAX */
#line 378 "zgbmv.f"
		    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 378 "zgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 378 "zgbmv.f"
		    i__2 = min(i__5,i__6);
#line 378 "zgbmv.f"
		    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 379 "zgbmv.f"
			d_cnjg(&z__3, &a[k + i__ + j * a_dim1]);
#line 379 "zgbmv.f"
			i__3 = ix;
#line 379 "zgbmv.f"
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
#line 379 "zgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 379 "zgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 380 "zgbmv.f"
			ix += *incx;
#line 381 "zgbmv.f"
/* L130: */
#line 381 "zgbmv.f"
		    }
#line 382 "zgbmv.f"
		}
#line 383 "zgbmv.f"
		i__2 = jy;
#line 383 "zgbmv.f"
		i__3 = jy;
#line 383 "zgbmv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 383 "zgbmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 383 "zgbmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 384 "zgbmv.f"
		jy += *incy;
#line 385 "zgbmv.f"
		if (j > *ku) {
#line 385 "zgbmv.f"
		    kx += *incx;
#line 385 "zgbmv.f"
		}
#line 386 "zgbmv.f"
/* L140: */
#line 386 "zgbmv.f"
	    }
#line 387 "zgbmv.f"
	}
#line 388 "zgbmv.f"
    }

#line 390 "zgbmv.f"
    return 0;

/*     End of ZGBMV . */

} /* zgbmv_ */

