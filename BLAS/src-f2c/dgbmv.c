#line 1 "dgbmv.f"
/* dgbmv.f -- translated by f2c (version 20100827).
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

#line 1 "dgbmv.f"
/* > \brief \b DGBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER INCX,INCY,KL,KU,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */
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
/* >              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
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
/* >          X is DOUBLE PRECISION array of DIMENSION at least */
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
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array of DIMENSION at least */
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

/* > \date November 2015 */

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
/* Subroutine */ int dgbmv_(char *trans, integer *m, integer *n, integer *kl, 
	integer *ku, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *x, integer *incx, doublereal *beta, doublereal *y, 
	integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
    static doublereal temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.6.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 225 "dgbmv.f"
    /* Parameter adjustments */
#line 225 "dgbmv.f"
    a_dim1 = *lda;
#line 225 "dgbmv.f"
    a_offset = 1 + a_dim1;
#line 225 "dgbmv.f"
    a -= a_offset;
#line 225 "dgbmv.f"
    --x;
#line 225 "dgbmv.f"
    --y;
#line 225 "dgbmv.f"

#line 225 "dgbmv.f"
    /* Function Body */
#line 225 "dgbmv.f"
    info = 0;
#line 226 "dgbmv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 228 "dgbmv.f"
	info = 1;
#line 229 "dgbmv.f"
    } else if (*m < 0) {
#line 230 "dgbmv.f"
	info = 2;
#line 231 "dgbmv.f"
    } else if (*n < 0) {
#line 232 "dgbmv.f"
	info = 3;
#line 233 "dgbmv.f"
    } else if (*kl < 0) {
#line 234 "dgbmv.f"
	info = 4;
#line 235 "dgbmv.f"
    } else if (*ku < 0) {
#line 236 "dgbmv.f"
	info = 5;
#line 237 "dgbmv.f"
    } else if (*lda < *kl + *ku + 1) {
#line 238 "dgbmv.f"
	info = 8;
#line 239 "dgbmv.f"
    } else if (*incx == 0) {
#line 240 "dgbmv.f"
	info = 10;
#line 241 "dgbmv.f"
    } else if (*incy == 0) {
#line 242 "dgbmv.f"
	info = 13;
#line 243 "dgbmv.f"
    }
#line 244 "dgbmv.f"
    if (info != 0) {
#line 245 "dgbmv.f"
	xerbla_("DGBMV ", &info, (ftnlen)6);
#line 246 "dgbmv.f"
	return 0;
#line 247 "dgbmv.f"
    }

/*     Quick return if possible. */

#line 251 "dgbmv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 251 "dgbmv.f"
	return 0;
#line 251 "dgbmv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 257 "dgbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 258 "dgbmv.f"
	lenx = *n;
#line 259 "dgbmv.f"
	leny = *m;
#line 260 "dgbmv.f"
    } else {
#line 261 "dgbmv.f"
	lenx = *m;
#line 262 "dgbmv.f"
	leny = *n;
#line 263 "dgbmv.f"
    }
#line 264 "dgbmv.f"
    if (*incx > 0) {
#line 265 "dgbmv.f"
	kx = 1;
#line 266 "dgbmv.f"
    } else {
#line 267 "dgbmv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 268 "dgbmv.f"
    }
#line 269 "dgbmv.f"
    if (*incy > 0) {
#line 270 "dgbmv.f"
	ky = 1;
#line 271 "dgbmv.f"
    } else {
#line 272 "dgbmv.f"
	ky = 1 - (leny - 1) * *incy;
#line 273 "dgbmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the band part of A. */

/*     First form  y := beta*y. */

#line 280 "dgbmv.f"
    if (*beta != 1.) {
#line 281 "dgbmv.f"
	if (*incy == 1) {
#line 282 "dgbmv.f"
	    if (*beta == 0.) {
#line 283 "dgbmv.f"
		i__1 = leny;
#line 283 "dgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 284 "dgbmv.f"
		    y[i__] = 0.;
#line 285 "dgbmv.f"
/* L10: */
#line 285 "dgbmv.f"
		}
#line 286 "dgbmv.f"
	    } else {
#line 287 "dgbmv.f"
		i__1 = leny;
#line 287 "dgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 288 "dgbmv.f"
		    y[i__] = *beta * y[i__];
#line 289 "dgbmv.f"
/* L20: */
#line 289 "dgbmv.f"
		}
#line 290 "dgbmv.f"
	    }
#line 291 "dgbmv.f"
	} else {
#line 292 "dgbmv.f"
	    iy = ky;
#line 293 "dgbmv.f"
	    if (*beta == 0.) {
#line 294 "dgbmv.f"
		i__1 = leny;
#line 294 "dgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 295 "dgbmv.f"
		    y[iy] = 0.;
#line 296 "dgbmv.f"
		    iy += *incy;
#line 297 "dgbmv.f"
/* L30: */
#line 297 "dgbmv.f"
		}
#line 298 "dgbmv.f"
	    } else {
#line 299 "dgbmv.f"
		i__1 = leny;
#line 299 "dgbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 300 "dgbmv.f"
		    y[iy] = *beta * y[iy];
#line 301 "dgbmv.f"
		    iy += *incy;
#line 302 "dgbmv.f"
/* L40: */
#line 302 "dgbmv.f"
		}
#line 303 "dgbmv.f"
	    }
#line 304 "dgbmv.f"
	}
#line 305 "dgbmv.f"
    }
#line 306 "dgbmv.f"
    if (*alpha == 0.) {
#line 306 "dgbmv.f"
	return 0;
#line 306 "dgbmv.f"
    }
#line 307 "dgbmv.f"
    kup1 = *ku + 1;
#line 308 "dgbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 312 "dgbmv.f"
	jx = kx;
#line 313 "dgbmv.f"
	if (*incy == 1) {
#line 314 "dgbmv.f"
	    i__1 = *n;
#line 314 "dgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 315 "dgbmv.f"
		temp = *alpha * x[jx];
#line 316 "dgbmv.f"
		k = kup1 - j;
/* Computing MAX */
#line 317 "dgbmv.f"
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 317 "dgbmv.f"
		i__5 = *m, i__6 = j + *kl;
#line 317 "dgbmv.f"
		i__4 = min(i__5,i__6);
#line 317 "dgbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 318 "dgbmv.f"
		    y[i__] += temp * a[k + i__ + j * a_dim1];
#line 319 "dgbmv.f"
/* L50: */
#line 319 "dgbmv.f"
		}
#line 320 "dgbmv.f"
		jx += *incx;
#line 321 "dgbmv.f"
/* L60: */
#line 321 "dgbmv.f"
	    }
#line 322 "dgbmv.f"
	} else {
#line 323 "dgbmv.f"
	    i__1 = *n;
#line 323 "dgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 324 "dgbmv.f"
		temp = *alpha * x[jx];
#line 325 "dgbmv.f"
		iy = ky;
#line 326 "dgbmv.f"
		k = kup1 - j;
/* Computing MAX */
#line 327 "dgbmv.f"
		i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
#line 327 "dgbmv.f"
		i__5 = *m, i__6 = j + *kl;
#line 327 "dgbmv.f"
		i__3 = min(i__5,i__6);
#line 327 "dgbmv.f"
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 328 "dgbmv.f"
		    y[iy] += temp * a[k + i__ + j * a_dim1];
#line 329 "dgbmv.f"
		    iy += *incy;
#line 330 "dgbmv.f"
/* L70: */
#line 330 "dgbmv.f"
		}
#line 331 "dgbmv.f"
		jx += *incx;
#line 332 "dgbmv.f"
		if (j > *ku) {
#line 332 "dgbmv.f"
		    ky += *incy;
#line 332 "dgbmv.f"
		}
#line 333 "dgbmv.f"
/* L80: */
#line 333 "dgbmv.f"
	    }
#line 334 "dgbmv.f"
	}
#line 335 "dgbmv.f"
    } else {

/*        Form  y := alpha*A**T*x + y. */

#line 339 "dgbmv.f"
	jy = ky;
#line 340 "dgbmv.f"
	if (*incx == 1) {
#line 341 "dgbmv.f"
	    i__1 = *n;
#line 341 "dgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 342 "dgbmv.f"
		temp = 0.;
#line 343 "dgbmv.f"
		k = kup1 - j;
/* Computing MAX */
#line 344 "dgbmv.f"
		i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
#line 344 "dgbmv.f"
		i__5 = *m, i__6 = j + *kl;
#line 344 "dgbmv.f"
		i__2 = min(i__5,i__6);
#line 344 "dgbmv.f"
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 345 "dgbmv.f"
		    temp += a[k + i__ + j * a_dim1] * x[i__];
#line 346 "dgbmv.f"
/* L90: */
#line 346 "dgbmv.f"
		}
#line 347 "dgbmv.f"
		y[jy] += *alpha * temp;
#line 348 "dgbmv.f"
		jy += *incy;
#line 349 "dgbmv.f"
/* L100: */
#line 349 "dgbmv.f"
	    }
#line 350 "dgbmv.f"
	} else {
#line 351 "dgbmv.f"
	    i__1 = *n;
#line 351 "dgbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 352 "dgbmv.f"
		temp = 0.;
#line 353 "dgbmv.f"
		ix = kx;
#line 354 "dgbmv.f"
		k = kup1 - j;
/* Computing MAX */
#line 355 "dgbmv.f"
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 355 "dgbmv.f"
		i__5 = *m, i__6 = j + *kl;
#line 355 "dgbmv.f"
		i__4 = min(i__5,i__6);
#line 355 "dgbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 356 "dgbmv.f"
		    temp += a[k + i__ + j * a_dim1] * x[ix];
#line 357 "dgbmv.f"
		    ix += *incx;
#line 358 "dgbmv.f"
/* L110: */
#line 358 "dgbmv.f"
		}
#line 359 "dgbmv.f"
		y[jy] += *alpha * temp;
#line 360 "dgbmv.f"
		jy += *incy;
#line 361 "dgbmv.f"
		if (j > *ku) {
#line 361 "dgbmv.f"
		    kx += *incx;
#line 361 "dgbmv.f"
		}
#line 362 "dgbmv.f"
/* L120: */
#line 362 "dgbmv.f"
	    }
#line 363 "dgbmv.f"
	}
#line 364 "dgbmv.f"
    }

#line 366 "dgbmv.f"
    return 0;

/*     End of DGBMV . */

} /* dgbmv_ */

