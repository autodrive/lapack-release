#line 1 "sgemv.f"
/* sgemv.f -- translated by f2c (version 20100827).
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

#line 1 "sgemv.f"
/* > \brief \b SGEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMV  performs one of the matrix-vector operations */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
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
/* >          X is REAL array of DIMENSION at least */
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
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array of DIMENSION at least */
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

/* > \date December 2016 */

/* > \ingroup single_blas_level2 */

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
/* Subroutine */ int sgemv_(char *trans, integer *m, integer *n, doublereal *
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

#line 196 "sgemv.f"
    /* Parameter adjustments */
#line 196 "sgemv.f"
    a_dim1 = *lda;
#line 196 "sgemv.f"
    a_offset = 1 + a_dim1;
#line 196 "sgemv.f"
    a -= a_offset;
#line 196 "sgemv.f"
    --x;
#line 196 "sgemv.f"
    --y;
#line 196 "sgemv.f"

#line 196 "sgemv.f"
    /* Function Body */
#line 196 "sgemv.f"
    info = 0;
#line 197 "sgemv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 199 "sgemv.f"
	info = 1;
#line 200 "sgemv.f"
    } else if (*m < 0) {
#line 201 "sgemv.f"
	info = 2;
#line 202 "sgemv.f"
    } else if (*n < 0) {
#line 203 "sgemv.f"
	info = 3;
#line 204 "sgemv.f"
    } else if (*lda < max(1,*m)) {
#line 205 "sgemv.f"
	info = 6;
#line 206 "sgemv.f"
    } else if (*incx == 0) {
#line 207 "sgemv.f"
	info = 8;
#line 208 "sgemv.f"
    } else if (*incy == 0) {
#line 209 "sgemv.f"
	info = 11;
#line 210 "sgemv.f"
    }
#line 211 "sgemv.f"
    if (info != 0) {
#line 212 "sgemv.f"
	xerbla_("SGEMV ", &info, (ftnlen)6);
#line 213 "sgemv.f"
	return 0;
#line 214 "sgemv.f"
    }

/*     Quick return if possible. */

#line 218 "sgemv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 218 "sgemv.f"
	return 0;
#line 218 "sgemv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 224 "sgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 225 "sgemv.f"
	lenx = *n;
#line 226 "sgemv.f"
	leny = *m;
#line 227 "sgemv.f"
    } else {
#line 228 "sgemv.f"
	lenx = *m;
#line 229 "sgemv.f"
	leny = *n;
#line 230 "sgemv.f"
    }
#line 231 "sgemv.f"
    if (*incx > 0) {
#line 232 "sgemv.f"
	kx = 1;
#line 233 "sgemv.f"
    } else {
#line 234 "sgemv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 235 "sgemv.f"
    }
#line 236 "sgemv.f"
    if (*incy > 0) {
#line 237 "sgemv.f"
	ky = 1;
#line 238 "sgemv.f"
    } else {
#line 239 "sgemv.f"
	ky = 1 - (leny - 1) * *incy;
#line 240 "sgemv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 247 "sgemv.f"
    if (*beta != 1.) {
#line 248 "sgemv.f"
	if (*incy == 1) {
#line 249 "sgemv.f"
	    if (*beta == 0.) {
#line 250 "sgemv.f"
		i__1 = leny;
#line 250 "sgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "sgemv.f"
		    y[i__] = 0.;
#line 252 "sgemv.f"
/* L10: */
#line 252 "sgemv.f"
		}
#line 253 "sgemv.f"
	    } else {
#line 254 "sgemv.f"
		i__1 = leny;
#line 254 "sgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 255 "sgemv.f"
		    y[i__] = *beta * y[i__];
#line 256 "sgemv.f"
/* L20: */
#line 256 "sgemv.f"
		}
#line 257 "sgemv.f"
	    }
#line 258 "sgemv.f"
	} else {
#line 259 "sgemv.f"
	    iy = ky;
#line 260 "sgemv.f"
	    if (*beta == 0.) {
#line 261 "sgemv.f"
		i__1 = leny;
#line 261 "sgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "sgemv.f"
		    y[iy] = 0.;
#line 263 "sgemv.f"
		    iy += *incy;
#line 264 "sgemv.f"
/* L30: */
#line 264 "sgemv.f"
		}
#line 265 "sgemv.f"
	    } else {
#line 266 "sgemv.f"
		i__1 = leny;
#line 266 "sgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "sgemv.f"
		    y[iy] = *beta * y[iy];
#line 268 "sgemv.f"
		    iy += *incy;
#line 269 "sgemv.f"
/* L40: */
#line 269 "sgemv.f"
		}
#line 270 "sgemv.f"
	    }
#line 271 "sgemv.f"
	}
#line 272 "sgemv.f"
    }
#line 273 "sgemv.f"
    if (*alpha == 0.) {
#line 273 "sgemv.f"
	return 0;
#line 273 "sgemv.f"
    }
#line 274 "sgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 278 "sgemv.f"
	jx = kx;
#line 279 "sgemv.f"
	if (*incy == 1) {
#line 280 "sgemv.f"
	    i__1 = *n;
#line 280 "sgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 281 "sgemv.f"
		temp = *alpha * x[jx];
#line 282 "sgemv.f"
		i__2 = *m;
#line 282 "sgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 283 "sgemv.f"
		    y[i__] += temp * a[i__ + j * a_dim1];
#line 284 "sgemv.f"
/* L50: */
#line 284 "sgemv.f"
		}
#line 285 "sgemv.f"
		jx += *incx;
#line 286 "sgemv.f"
/* L60: */
#line 286 "sgemv.f"
	    }
#line 287 "sgemv.f"
	} else {
#line 288 "sgemv.f"
	    i__1 = *n;
#line 288 "sgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 289 "sgemv.f"
		temp = *alpha * x[jx];
#line 290 "sgemv.f"
		iy = ky;
#line 291 "sgemv.f"
		i__2 = *m;
#line 291 "sgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 292 "sgemv.f"
		    y[iy] += temp * a[i__ + j * a_dim1];
#line 293 "sgemv.f"
		    iy += *incy;
#line 294 "sgemv.f"
/* L70: */
#line 294 "sgemv.f"
		}
#line 295 "sgemv.f"
		jx += *incx;
#line 296 "sgemv.f"
/* L80: */
#line 296 "sgemv.f"
	    }
#line 297 "sgemv.f"
	}
#line 298 "sgemv.f"
    } else {

/*        Form  y := alpha*A**T*x + y. */

#line 302 "sgemv.f"
	jy = ky;
#line 303 "sgemv.f"
	if (*incx == 1) {
#line 304 "sgemv.f"
	    i__1 = *n;
#line 304 "sgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 305 "sgemv.f"
		temp = 0.;
#line 306 "sgemv.f"
		i__2 = *m;
#line 306 "sgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 307 "sgemv.f"
		    temp += a[i__ + j * a_dim1] * x[i__];
#line 308 "sgemv.f"
/* L90: */
#line 308 "sgemv.f"
		}
#line 309 "sgemv.f"
		y[jy] += *alpha * temp;
#line 310 "sgemv.f"
		jy += *incy;
#line 311 "sgemv.f"
/* L100: */
#line 311 "sgemv.f"
	    }
#line 312 "sgemv.f"
	} else {
#line 313 "sgemv.f"
	    i__1 = *n;
#line 313 "sgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 314 "sgemv.f"
		temp = 0.;
#line 315 "sgemv.f"
		ix = kx;
#line 316 "sgemv.f"
		i__2 = *m;
#line 316 "sgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 317 "sgemv.f"
		    temp += a[i__ + j * a_dim1] * x[ix];
#line 318 "sgemv.f"
		    ix += *incx;
#line 319 "sgemv.f"
/* L110: */
#line 319 "sgemv.f"
		}
#line 320 "sgemv.f"
		y[jy] += *alpha * temp;
#line 321 "sgemv.f"
		jy += *incy;
#line 322 "sgemv.f"
/* L120: */
#line 322 "sgemv.f"
	    }
#line 323 "sgemv.f"
	}
#line 324 "sgemv.f"
    }

#line 326 "sgemv.f"
    return 0;

/*     End of SGEMV . */

} /* sgemv_ */

