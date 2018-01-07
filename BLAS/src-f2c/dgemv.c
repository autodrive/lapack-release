#line 1 "dgemv.f"
/* dgemv.f -- translated by f2c (version 20100827).
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

#line 1 "dgemv.f"
/* > \brief \b DGEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,M,N */
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
/* > DGEMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
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
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, m ). */
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
/* >           Before entry with BETA non-zero, the incremented array Y */
/* >           must contain the vector y. On exit, Y is overwritten by the */
/* >           updated vector y. */
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
/* Subroutine */ int dgemv_(char *trans, integer *m, integer *n, doublereal *
	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static doublereal temp;
    static integer lenx, leny;
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

#line 196 "dgemv.f"
    /* Parameter adjustments */
#line 196 "dgemv.f"
    a_dim1 = *lda;
#line 196 "dgemv.f"
    a_offset = 1 + a_dim1;
#line 196 "dgemv.f"
    a -= a_offset;
#line 196 "dgemv.f"
    --x;
#line 196 "dgemv.f"
    --y;
#line 196 "dgemv.f"

#line 196 "dgemv.f"
    /* Function Body */
#line 196 "dgemv.f"
    info = 0;
#line 197 "dgemv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 199 "dgemv.f"
	info = 1;
#line 200 "dgemv.f"
    } else if (*m < 0) {
#line 201 "dgemv.f"
	info = 2;
#line 202 "dgemv.f"
    } else if (*n < 0) {
#line 203 "dgemv.f"
	info = 3;
#line 204 "dgemv.f"
    } else if (*lda < max(1,*m)) {
#line 205 "dgemv.f"
	info = 6;
#line 206 "dgemv.f"
    } else if (*incx == 0) {
#line 207 "dgemv.f"
	info = 8;
#line 208 "dgemv.f"
    } else if (*incy == 0) {
#line 209 "dgemv.f"
	info = 11;
#line 210 "dgemv.f"
    }
#line 211 "dgemv.f"
    if (info != 0) {
#line 212 "dgemv.f"
	xerbla_("DGEMV ", &info, (ftnlen)6);
#line 213 "dgemv.f"
	return 0;
#line 214 "dgemv.f"
    }

/*     Quick return if possible. */

#line 218 "dgemv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 218 "dgemv.f"
	return 0;
#line 218 "dgemv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 224 "dgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 225 "dgemv.f"
	lenx = *n;
#line 226 "dgemv.f"
	leny = *m;
#line 227 "dgemv.f"
    } else {
#line 228 "dgemv.f"
	lenx = *m;
#line 229 "dgemv.f"
	leny = *n;
#line 230 "dgemv.f"
    }
#line 231 "dgemv.f"
    if (*incx > 0) {
#line 232 "dgemv.f"
	kx = 1;
#line 233 "dgemv.f"
    } else {
#line 234 "dgemv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 235 "dgemv.f"
    }
#line 236 "dgemv.f"
    if (*incy > 0) {
#line 237 "dgemv.f"
	ky = 1;
#line 238 "dgemv.f"
    } else {
#line 239 "dgemv.f"
	ky = 1 - (leny - 1) * *incy;
#line 240 "dgemv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 247 "dgemv.f"
    if (*beta != 1.) {
#line 248 "dgemv.f"
	if (*incy == 1) {
#line 249 "dgemv.f"
	    if (*beta == 0.) {
#line 250 "dgemv.f"
		i__1 = leny;
#line 250 "dgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "dgemv.f"
		    y[i__] = 0.;
#line 252 "dgemv.f"
/* L10: */
#line 252 "dgemv.f"
		}
#line 253 "dgemv.f"
	    } else {
#line 254 "dgemv.f"
		i__1 = leny;
#line 254 "dgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 255 "dgemv.f"
		    y[i__] = *beta * y[i__];
#line 256 "dgemv.f"
/* L20: */
#line 256 "dgemv.f"
		}
#line 257 "dgemv.f"
	    }
#line 258 "dgemv.f"
	} else {
#line 259 "dgemv.f"
	    iy = ky;
#line 260 "dgemv.f"
	    if (*beta == 0.) {
#line 261 "dgemv.f"
		i__1 = leny;
#line 261 "dgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "dgemv.f"
		    y[iy] = 0.;
#line 263 "dgemv.f"
		    iy += *incy;
#line 264 "dgemv.f"
/* L30: */
#line 264 "dgemv.f"
		}
#line 265 "dgemv.f"
	    } else {
#line 266 "dgemv.f"
		i__1 = leny;
#line 266 "dgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "dgemv.f"
		    y[iy] = *beta * y[iy];
#line 268 "dgemv.f"
		    iy += *incy;
#line 269 "dgemv.f"
/* L40: */
#line 269 "dgemv.f"
		}
#line 270 "dgemv.f"
	    }
#line 271 "dgemv.f"
	}
#line 272 "dgemv.f"
    }
#line 273 "dgemv.f"
    if (*alpha == 0.) {
#line 273 "dgemv.f"
	return 0;
#line 273 "dgemv.f"
    }
#line 274 "dgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 278 "dgemv.f"
	jx = kx;
#line 279 "dgemv.f"
	if (*incy == 1) {
#line 280 "dgemv.f"
	    i__1 = *n;
#line 280 "dgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 281 "dgemv.f"
		if (x[jx] != 0.) {
#line 282 "dgemv.f"
		    temp = *alpha * x[jx];
#line 283 "dgemv.f"
		    i__2 = *m;
#line 283 "dgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 284 "dgemv.f"
			y[i__] += temp * a[i__ + j * a_dim1];
#line 285 "dgemv.f"
/* L50: */
#line 285 "dgemv.f"
		    }
#line 286 "dgemv.f"
		}
#line 287 "dgemv.f"
		jx += *incx;
#line 288 "dgemv.f"
/* L60: */
#line 288 "dgemv.f"
	    }
#line 289 "dgemv.f"
	} else {
#line 290 "dgemv.f"
	    i__1 = *n;
#line 290 "dgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 291 "dgemv.f"
		if (x[jx] != 0.) {
#line 292 "dgemv.f"
		    temp = *alpha * x[jx];
#line 293 "dgemv.f"
		    iy = ky;
#line 294 "dgemv.f"
		    i__2 = *m;
#line 294 "dgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "dgemv.f"
			y[iy] += temp * a[i__ + j * a_dim1];
#line 296 "dgemv.f"
			iy += *incy;
#line 297 "dgemv.f"
/* L70: */
#line 297 "dgemv.f"
		    }
#line 298 "dgemv.f"
		}
#line 299 "dgemv.f"
		jx += *incx;
#line 300 "dgemv.f"
/* L80: */
#line 300 "dgemv.f"
	    }
#line 301 "dgemv.f"
	}
#line 302 "dgemv.f"
    } else {

/*        Form  y := alpha*A**T*x + y. */

#line 306 "dgemv.f"
	jy = ky;
#line 307 "dgemv.f"
	if (*incx == 1) {
#line 308 "dgemv.f"
	    i__1 = *n;
#line 308 "dgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 309 "dgemv.f"
		temp = 0.;
#line 310 "dgemv.f"
		i__2 = *m;
#line 310 "dgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 311 "dgemv.f"
		    temp += a[i__ + j * a_dim1] * x[i__];
#line 312 "dgemv.f"
/* L90: */
#line 312 "dgemv.f"
		}
#line 313 "dgemv.f"
		y[jy] += *alpha * temp;
#line 314 "dgemv.f"
		jy += *incy;
#line 315 "dgemv.f"
/* L100: */
#line 315 "dgemv.f"
	    }
#line 316 "dgemv.f"
	} else {
#line 317 "dgemv.f"
	    i__1 = *n;
#line 317 "dgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 318 "dgemv.f"
		temp = 0.;
#line 319 "dgemv.f"
		ix = kx;
#line 320 "dgemv.f"
		i__2 = *m;
#line 320 "dgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 321 "dgemv.f"
		    temp += a[i__ + j * a_dim1] * x[ix];
#line 322 "dgemv.f"
		    ix += *incx;
#line 323 "dgemv.f"
/* L110: */
#line 323 "dgemv.f"
		}
#line 324 "dgemv.f"
		y[jy] += *alpha * temp;
#line 325 "dgemv.f"
		jy += *incy;
#line 326 "dgemv.f"
/* L120: */
#line 326 "dgemv.f"
	    }
#line 327 "dgemv.f"
	}
#line 328 "dgemv.f"
    }

#line 330 "dgemv.f"
    return 0;

/*     End of DGEMV . */

} /* dgemv_ */

