#line 1 "cgbmv.f"
/* cgbmv.f -- translated by f2c (version 20100827).
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

#line 1 "cgbmv.f"
/* > \brief \b CGBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER INCX,INCY,KL,KU,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBMV  performs one of the matrix-vector operations */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
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
/* >          X is COMPLEX array of DIMENSION at least */
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
/* >          BETA is COMPLEX */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array of DIMENSION at least */
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
/* Subroutine */ int cgbmv_(char *trans, integer *m, integer *n, integer *kl, 
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

#line 230 "cgbmv.f"
    /* Parameter adjustments */
#line 230 "cgbmv.f"
    a_dim1 = *lda;
#line 230 "cgbmv.f"
    a_offset = 1 + a_dim1;
#line 230 "cgbmv.f"
    a -= a_offset;
#line 230 "cgbmv.f"
    --x;
#line 230 "cgbmv.f"
    --y;
#line 230 "cgbmv.f"

#line 230 "cgbmv.f"
    /* Function Body */
#line 230 "cgbmv.f"
    info = 0;
#line 231 "cgbmv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 233 "cgbmv.f"
	info = 1;
#line 234 "cgbmv.f"
    } else if (*m < 0) {
#line 235 "cgbmv.f"
	info = 2;
#line 236 "cgbmv.f"
    } else if (*n < 0) {
#line 237 "cgbmv.f"
	info = 3;
#line 238 "cgbmv.f"
    } else if (*kl < 0) {
#line 239 "cgbmv.f"
	info = 4;
#line 240 "cgbmv.f"
    } else if (*ku < 0) {
#line 241 "cgbmv.f"
	info = 5;
#line 242 "cgbmv.f"
    } else if (*lda < *kl + *ku + 1) {
#line 243 "cgbmv.f"
	info = 8;
#line 244 "cgbmv.f"
    } else if (*incx == 0) {
#line 245 "cgbmv.f"
	info = 10;
#line 246 "cgbmv.f"
    } else if (*incy == 0) {
#line 247 "cgbmv.f"
	info = 13;
#line 248 "cgbmv.f"
    }
#line 249 "cgbmv.f"
    if (info != 0) {
#line 250 "cgbmv.f"
	xerbla_("CGBMV ", &info, (ftnlen)6);
#line 251 "cgbmv.f"
	return 0;
#line 252 "cgbmv.f"
    }

/*     Quick return if possible. */

#line 256 "cgbmv.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
#line 256 "cgbmv.f"
	return 0;
#line 256 "cgbmv.f"
    }

#line 259 "cgbmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 264 "cgbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 265 "cgbmv.f"
	lenx = *n;
#line 266 "cgbmv.f"
	leny = *m;
#line 267 "cgbmv.f"
    } else {
#line 268 "cgbmv.f"
	lenx = *m;
#line 269 "cgbmv.f"
	leny = *n;
#line 270 "cgbmv.f"
    }
#line 271 "cgbmv.f"
    if (*incx > 0) {
#line 272 "cgbmv.f"
	kx = 1;
#line 273 "cgbmv.f"
    } else {
#line 274 "cgbmv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 275 "cgbmv.f"
    }
#line 276 "cgbmv.f"
    if (*incy > 0) {
#line 277 "cgbmv.f"
	ky = 1;
#line 278 "cgbmv.f"
    } else {
#line 279 "cgbmv.f"
	ky = 1 - (leny - 1) * *incy;
#line 280 "cgbmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the band part of A. */

/*     First form  y := beta*y. */

#line 287 "cgbmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 288 "cgbmv.f"
	if (*incy == 1) {
#line 289 "cgbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 290 "cgbmv.f"
		i__1 = leny;
#line 290 "cgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "cgbmv.f"
		    i__2 = i__;
#line 291 "cgbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 292 "cgbmv.f"
/* L10: */
#line 292 "cgbmv.f"
		}
#line 293 "cgbmv.f"
	    } else {
#line 294 "cgbmv.f"
		i__1 = leny;
#line 294 "cgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 295 "cgbmv.f"
		    i__2 = i__;
#line 295 "cgbmv.f"
		    i__3 = i__;
#line 295 "cgbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 295 "cgbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 296 "cgbmv.f"
/* L20: */
#line 296 "cgbmv.f"
		}
#line 297 "cgbmv.f"
	    }
#line 298 "cgbmv.f"
	} else {
#line 299 "cgbmv.f"
	    iy = ky;
#line 300 "cgbmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 301 "cgbmv.f"
		i__1 = leny;
#line 301 "cgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "cgbmv.f"
		    i__2 = iy;
#line 302 "cgbmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 303 "cgbmv.f"
		    iy += *incy;
#line 304 "cgbmv.f"
/* L30: */
#line 304 "cgbmv.f"
		}
#line 305 "cgbmv.f"
	    } else {
#line 306 "cgbmv.f"
		i__1 = leny;
#line 306 "cgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 307 "cgbmv.f"
		    i__2 = iy;
#line 307 "cgbmv.f"
		    i__3 = iy;
#line 307 "cgbmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 307 "cgbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 308 "cgbmv.f"
		    iy += *incy;
#line 309 "cgbmv.f"
/* L40: */
#line 309 "cgbmv.f"
		}
#line 310 "cgbmv.f"
	    }
#line 311 "cgbmv.f"
	}
#line 312 "cgbmv.f"
    }
#line 313 "cgbmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 313 "cgbmv.f"
	return 0;
#line 313 "cgbmv.f"
    }
#line 314 "cgbmv.f"
    kup1 = *ku + 1;
#line 315 "cgbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 319 "cgbmv.f"
	jx = kx;
#line 320 "cgbmv.f"
	if (*incy == 1) {
#line 321 "cgbmv.f"
	    i__1 = *n;
#line 321 "cgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 322 "cgbmv.f"
		i__2 = jx;
#line 322 "cgbmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 322 "cgbmv.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 323 "cgbmv.f"
		k = kup1 - j;
/* Computing MAX */
#line 324 "cgbmv.f"
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 324 "cgbmv.f"
		i__5 = *m, i__6 = j + *kl;
#line 324 "cgbmv.f"
		i__4 = min(i__5,i__6);
#line 324 "cgbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 325 "cgbmv.f"
		    i__2 = i__;
#line 325 "cgbmv.f"
		    i__3 = i__;
#line 325 "cgbmv.f"
		    i__5 = k + i__ + j * a_dim1;
#line 325 "cgbmv.f"
		    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, z__2.i =
			     temp.r * a[i__5].i + temp.i * a[i__5].r;
#line 325 "cgbmv.f"
		    z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 325 "cgbmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 326 "cgbmv.f"
/* L50: */
#line 326 "cgbmv.f"
		}
#line 327 "cgbmv.f"
		jx += *incx;
#line 328 "cgbmv.f"
/* L60: */
#line 328 "cgbmv.f"
	    }
#line 329 "cgbmv.f"
	} else {
#line 330 "cgbmv.f"
	    i__1 = *n;
#line 330 "cgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 331 "cgbmv.f"
		i__4 = jx;
#line 331 "cgbmv.f"
		z__1.r = alpha->r * x[i__4].r - alpha->i * x[i__4].i, z__1.i =
			 alpha->r * x[i__4].i + alpha->i * x[i__4].r;
#line 331 "cgbmv.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 332 "cgbmv.f"
		iy = ky;
#line 333 "cgbmv.f"
		k = kup1 - j;
/* Computing MAX */
#line 334 "cgbmv.f"
		i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 334 "cgbmv.f"
		i__5 = *m, i__6 = j + *kl;
#line 334 "cgbmv.f"
		i__3 = min(i__5,i__6);
#line 334 "cgbmv.f"
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 335 "cgbmv.f"
		    i__4 = iy;
#line 335 "cgbmv.f"
		    i__2 = iy;
#line 335 "cgbmv.f"
		    i__5 = k + i__ + j * a_dim1;
#line 335 "cgbmv.f"
		    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, z__2.i =
			     temp.r * a[i__5].i + temp.i * a[i__5].r;
#line 335 "cgbmv.f"
		    z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 335 "cgbmv.f"
		    y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 336 "cgbmv.f"
		    iy += *incy;
#line 337 "cgbmv.f"
/* L70: */
#line 337 "cgbmv.f"
		}
#line 338 "cgbmv.f"
		jx += *incx;
#line 339 "cgbmv.f"
		if (j > *ku) {
#line 339 "cgbmv.f"
		    ky += *incy;
#line 339 "cgbmv.f"
		}
#line 340 "cgbmv.f"
/* L80: */
#line 340 "cgbmv.f"
	    }
#line 341 "cgbmv.f"
	}
#line 342 "cgbmv.f"
    } else {

/*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y. */

#line 346 "cgbmv.f"
	jy = ky;
#line 347 "cgbmv.f"
	if (*incx == 1) {
#line 348 "cgbmv.f"
	    i__1 = *n;
#line 348 "cgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 349 "cgbmv.f"
		temp.r = 0., temp.i = 0.;
#line 350 "cgbmv.f"
		k = kup1 - j;
#line 351 "cgbmv.f"
		if (noconj) {
/* Computing MAX */
#line 352 "cgbmv.f"
		    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 352 "cgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 352 "cgbmv.f"
		    i__2 = min(i__5,i__6);
#line 352 "cgbmv.f"
		    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 353 "cgbmv.f"
			i__3 = k + i__ + j * a_dim1;
#line 353 "cgbmv.f"
			i__4 = i__;
#line 353 "cgbmv.f"
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
#line 353 "cgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 353 "cgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 354 "cgbmv.f"
/* L90: */
#line 354 "cgbmv.f"
		    }
#line 355 "cgbmv.f"
		} else {
/* Computing MAX */
#line 356 "cgbmv.f"
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 356 "cgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 356 "cgbmv.f"
		    i__4 = min(i__5,i__6);
#line 356 "cgbmv.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 357 "cgbmv.f"
			d_cnjg(&z__3, &a[k + i__ + j * a_dim1]);
#line 357 "cgbmv.f"
			i__2 = i__;
#line 357 "cgbmv.f"
			z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				z__2.i = z__3.r * x[i__2].i + z__3.i * x[i__2]
				.r;
#line 357 "cgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 357 "cgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 358 "cgbmv.f"
/* L100: */
#line 358 "cgbmv.f"
		    }
#line 359 "cgbmv.f"
		}
#line 360 "cgbmv.f"
		i__4 = jy;
#line 360 "cgbmv.f"
		i__2 = jy;
#line 360 "cgbmv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 360 "cgbmv.f"
		z__1.r = y[i__2].r + z__2.r, z__1.i = y[i__2].i + z__2.i;
#line 360 "cgbmv.f"
		y[i__4].r = z__1.r, y[i__4].i = z__1.i;
#line 361 "cgbmv.f"
		jy += *incy;
#line 362 "cgbmv.f"
/* L110: */
#line 362 "cgbmv.f"
	    }
#line 363 "cgbmv.f"
	} else {
#line 364 "cgbmv.f"
	    i__1 = *n;
#line 364 "cgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 365 "cgbmv.f"
		temp.r = 0., temp.i = 0.;
#line 366 "cgbmv.f"
		ix = kx;
#line 367 "cgbmv.f"
		k = kup1 - j;
#line 368 "cgbmv.f"
		if (noconj) {
/* Computing MAX */
#line 369 "cgbmv.f"
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 369 "cgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 369 "cgbmv.f"
		    i__3 = min(i__5,i__6);
#line 369 "cgbmv.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 370 "cgbmv.f"
			i__4 = k + i__ + j * a_dim1;
#line 370 "cgbmv.f"
			i__2 = ix;
#line 370 "cgbmv.f"
			z__2.r = a[i__4].r * x[i__2].r - a[i__4].i * x[i__2]
				.i, z__2.i = a[i__4].r * x[i__2].i + a[i__4]
				.i * x[i__2].r;
#line 370 "cgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 370 "cgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 371 "cgbmv.f"
			ix += *incx;
#line 372 "cgbmv.f"
/* L120: */
#line 372 "cgbmv.f"
		    }
#line 373 "cgbmv.f"
		} else {
/* Computing MAX */
#line 374 "cgbmv.f"
		    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 374 "cgbmv.f"
		    i__5 = *m, i__6 = j + *kl;
#line 374 "cgbmv.f"
		    i__2 = min(i__5,i__6);
#line 374 "cgbmv.f"
		    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 375 "cgbmv.f"
			d_cnjg(&z__3, &a[k + i__ + j * a_dim1]);
#line 375 "cgbmv.f"
			i__3 = ix;
#line 375 "cgbmv.f"
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
#line 375 "cgbmv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 375 "cgbmv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 376 "cgbmv.f"
			ix += *incx;
#line 377 "cgbmv.f"
/* L130: */
#line 377 "cgbmv.f"
		    }
#line 378 "cgbmv.f"
		}
#line 379 "cgbmv.f"
		i__2 = jy;
#line 379 "cgbmv.f"
		i__3 = jy;
#line 379 "cgbmv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 379 "cgbmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 379 "cgbmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 380 "cgbmv.f"
		jy += *incy;
#line 381 "cgbmv.f"
		if (j > *ku) {
#line 381 "cgbmv.f"
		    kx += *incx;
#line 381 "cgbmv.f"
		}
#line 382 "cgbmv.f"
/* L140: */
#line 382 "cgbmv.f"
	    }
#line 383 "cgbmv.f"
	}
#line 384 "cgbmv.f"
    }

#line 386 "cgbmv.f"
    return 0;

/*     End of CGBMV . */

} /* cgbmv_ */

