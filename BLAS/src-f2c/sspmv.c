#line 1 "sspmv.f"
/* sspmv.f -- translated by f2c (version 20100827).
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

#line 1 "sspmv.f"
/* > \brief \b SSPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPMV  performs the matrix-vector operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix, supplied in packed form. */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array of dimension at least */
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
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array of dimension at least */
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
/* Subroutine */ int sspmv_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *ap, doublereal *x, integer *incx, doublereal *beta, 
	doublereal *y, integer *incy, ftnlen uplo_len)
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

#line 184 "sspmv.f"
    /* Parameter adjustments */
#line 184 "sspmv.f"
    --y;
#line 184 "sspmv.f"
    --x;
#line 184 "sspmv.f"
    --ap;
#line 184 "sspmv.f"

#line 184 "sspmv.f"
    /* Function Body */
#line 184 "sspmv.f"
    info = 0;
#line 185 "sspmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 186 "sspmv.f"
	info = 1;
#line 187 "sspmv.f"
    } else if (*n < 0) {
#line 188 "sspmv.f"
	info = 2;
#line 189 "sspmv.f"
    } else if (*incx == 0) {
#line 190 "sspmv.f"
	info = 6;
#line 191 "sspmv.f"
    } else if (*incy == 0) {
#line 192 "sspmv.f"
	info = 9;
#line 193 "sspmv.f"
    }
#line 194 "sspmv.f"
    if (info != 0) {
#line 195 "sspmv.f"
	xerbla_("SSPMV ", &info, (ftnlen)6);
#line 196 "sspmv.f"
	return 0;
#line 197 "sspmv.f"
    }

/*     Quick return if possible. */

#line 201 "sspmv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 201 "sspmv.f"
	return 0;
#line 201 "sspmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 205 "sspmv.f"
    if (*incx > 0) {
#line 206 "sspmv.f"
	kx = 1;
#line 207 "sspmv.f"
    } else {
#line 208 "sspmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 209 "sspmv.f"
    }
#line 210 "sspmv.f"
    if (*incy > 0) {
#line 211 "sspmv.f"
	ky = 1;
#line 212 "sspmv.f"
    } else {
#line 213 "sspmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 214 "sspmv.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

#line 221 "sspmv.f"
    if (*beta != 1.) {
#line 222 "sspmv.f"
	if (*incy == 1) {
#line 223 "sspmv.f"
	    if (*beta == 0.) {
#line 224 "sspmv.f"
		i__1 = *n;
#line 224 "sspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 225 "sspmv.f"
		    y[i__] = 0.;
#line 226 "sspmv.f"
/* L10: */
#line 226 "sspmv.f"
		}
#line 227 "sspmv.f"
	    } else {
#line 228 "sspmv.f"
		i__1 = *n;
#line 228 "sspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "sspmv.f"
		    y[i__] = *beta * y[i__];
#line 230 "sspmv.f"
/* L20: */
#line 230 "sspmv.f"
		}
#line 231 "sspmv.f"
	    }
#line 232 "sspmv.f"
	} else {
#line 233 "sspmv.f"
	    iy = ky;
#line 234 "sspmv.f"
	    if (*beta == 0.) {
#line 235 "sspmv.f"
		i__1 = *n;
#line 235 "sspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "sspmv.f"
		    y[iy] = 0.;
#line 237 "sspmv.f"
		    iy += *incy;
#line 238 "sspmv.f"
/* L30: */
#line 238 "sspmv.f"
		}
#line 239 "sspmv.f"
	    } else {
#line 240 "sspmv.f"
		i__1 = *n;
#line 240 "sspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "sspmv.f"
		    y[iy] = *beta * y[iy];
#line 242 "sspmv.f"
		    iy += *incy;
#line 243 "sspmv.f"
/* L40: */
#line 243 "sspmv.f"
		}
#line 244 "sspmv.f"
	    }
#line 245 "sspmv.f"
	}
#line 246 "sspmv.f"
    }
#line 247 "sspmv.f"
    if (*alpha == 0.) {
#line 247 "sspmv.f"
	return 0;
#line 247 "sspmv.f"
    }
#line 248 "sspmv.f"
    kk = 1;
#line 249 "sspmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

#line 253 "sspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 254 "sspmv.f"
	    i__1 = *n;
#line 254 "sspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 255 "sspmv.f"
		temp1 = *alpha * x[j];
#line 256 "sspmv.f"
		temp2 = 0.;
#line 257 "sspmv.f"
		k = kk;
#line 258 "sspmv.f"
		i__2 = j - 1;
#line 258 "sspmv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 259 "sspmv.f"
		    y[i__] += temp1 * ap[k];
#line 260 "sspmv.f"
		    temp2 += ap[k] * x[i__];
#line 261 "sspmv.f"
		    ++k;
#line 262 "sspmv.f"
/* L50: */
#line 262 "sspmv.f"
		}
#line 263 "sspmv.f"
		y[j] = y[j] + temp1 * ap[kk + j - 1] + *alpha * temp2;
#line 264 "sspmv.f"
		kk += j;
#line 265 "sspmv.f"
/* L60: */
#line 265 "sspmv.f"
	    }
#line 266 "sspmv.f"
	} else {
#line 267 "sspmv.f"
	    jx = kx;
#line 268 "sspmv.f"
	    jy = ky;
#line 269 "sspmv.f"
	    i__1 = *n;
#line 269 "sspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 270 "sspmv.f"
		temp1 = *alpha * x[jx];
#line 271 "sspmv.f"
		temp2 = 0.;
#line 272 "sspmv.f"
		ix = kx;
#line 273 "sspmv.f"
		iy = ky;
#line 274 "sspmv.f"
		i__2 = kk + j - 2;
#line 274 "sspmv.f"
		for (k = kk; k <= i__2; ++k) {
#line 275 "sspmv.f"
		    y[iy] += temp1 * ap[k];
#line 276 "sspmv.f"
		    temp2 += ap[k] * x[ix];
#line 277 "sspmv.f"
		    ix += *incx;
#line 278 "sspmv.f"
		    iy += *incy;
#line 279 "sspmv.f"
/* L70: */
#line 279 "sspmv.f"
		}
#line 280 "sspmv.f"
		y[jy] = y[jy] + temp1 * ap[kk + j - 1] + *alpha * temp2;
#line 281 "sspmv.f"
		jx += *incx;
#line 282 "sspmv.f"
		jy += *incy;
#line 283 "sspmv.f"
		kk += j;
#line 284 "sspmv.f"
/* L80: */
#line 284 "sspmv.f"
	    }
#line 285 "sspmv.f"
	}
#line 286 "sspmv.f"
    } else {

/*        Form  y  when AP contains the lower triangle. */

#line 290 "sspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 291 "sspmv.f"
	    i__1 = *n;
#line 291 "sspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "sspmv.f"
		temp1 = *alpha * x[j];
#line 293 "sspmv.f"
		temp2 = 0.;
#line 294 "sspmv.f"
		y[j] += temp1 * ap[kk];
#line 295 "sspmv.f"
		k = kk + 1;
#line 296 "sspmv.f"
		i__2 = *n;
#line 296 "sspmv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 297 "sspmv.f"
		    y[i__] += temp1 * ap[k];
#line 298 "sspmv.f"
		    temp2 += ap[k] * x[i__];
#line 299 "sspmv.f"
		    ++k;
#line 300 "sspmv.f"
/* L90: */
#line 300 "sspmv.f"
		}
#line 301 "sspmv.f"
		y[j] += *alpha * temp2;
#line 302 "sspmv.f"
		kk += *n - j + 1;
#line 303 "sspmv.f"
/* L100: */
#line 303 "sspmv.f"
	    }
#line 304 "sspmv.f"
	} else {
#line 305 "sspmv.f"
	    jx = kx;
#line 306 "sspmv.f"
	    jy = ky;
#line 307 "sspmv.f"
	    i__1 = *n;
#line 307 "sspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 308 "sspmv.f"
		temp1 = *alpha * x[jx];
#line 309 "sspmv.f"
		temp2 = 0.;
#line 310 "sspmv.f"
		y[jy] += temp1 * ap[kk];
#line 311 "sspmv.f"
		ix = jx;
#line 312 "sspmv.f"
		iy = jy;
#line 313 "sspmv.f"
		i__2 = kk + *n - j;
#line 313 "sspmv.f"
		for (k = kk + 1; k <= i__2; ++k) {
#line 314 "sspmv.f"
		    ix += *incx;
#line 315 "sspmv.f"
		    iy += *incy;
#line 316 "sspmv.f"
		    y[iy] += temp1 * ap[k];
#line 317 "sspmv.f"
		    temp2 += ap[k] * x[ix];
#line 318 "sspmv.f"
/* L110: */
#line 318 "sspmv.f"
		}
#line 319 "sspmv.f"
		y[jy] += *alpha * temp2;
#line 320 "sspmv.f"
		jx += *incx;
#line 321 "sspmv.f"
		jy += *incy;
#line 322 "sspmv.f"
		kk += *n - j + 1;
#line 323 "sspmv.f"
/* L120: */
#line 323 "sspmv.f"
	    }
#line 324 "sspmv.f"
	}
#line 325 "sspmv.f"
    }

#line 327 "sspmv.f"
    return 0;

/*     End of SSPMV . */

} /* sspmv_ */

