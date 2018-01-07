#line 1 "dspmv.f"
/* dspmv.f -- translated by f2c (version 20100827).
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

#line 1 "dspmv.f"
/* > \brief \b DSPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPMV  performs the matrix-vector operation */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array of DIMENSION at least */
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
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array of dimension at least */
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
/* Subroutine */ int dspmv_(char *uplo, integer *n, doublereal *alpha, 
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

/*     Test the input parameters. */

#line 184 "dspmv.f"
    /* Parameter adjustments */
#line 184 "dspmv.f"
    --y;
#line 184 "dspmv.f"
    --x;
#line 184 "dspmv.f"
    --ap;
#line 184 "dspmv.f"

#line 184 "dspmv.f"
    /* Function Body */
#line 184 "dspmv.f"
    info = 0;
#line 185 "dspmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 186 "dspmv.f"
	info = 1;
#line 187 "dspmv.f"
    } else if (*n < 0) {
#line 188 "dspmv.f"
	info = 2;
#line 189 "dspmv.f"
    } else if (*incx == 0) {
#line 190 "dspmv.f"
	info = 6;
#line 191 "dspmv.f"
    } else if (*incy == 0) {
#line 192 "dspmv.f"
	info = 9;
#line 193 "dspmv.f"
    }
#line 194 "dspmv.f"
    if (info != 0) {
#line 195 "dspmv.f"
	xerbla_("DSPMV ", &info, (ftnlen)6);
#line 196 "dspmv.f"
	return 0;
#line 197 "dspmv.f"
    }

/*     Quick return if possible. */

#line 201 "dspmv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 201 "dspmv.f"
	return 0;
#line 201 "dspmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 205 "dspmv.f"
    if (*incx > 0) {
#line 206 "dspmv.f"
	kx = 1;
#line 207 "dspmv.f"
    } else {
#line 208 "dspmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 209 "dspmv.f"
    }
#line 210 "dspmv.f"
    if (*incy > 0) {
#line 211 "dspmv.f"
	ky = 1;
#line 212 "dspmv.f"
    } else {
#line 213 "dspmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 214 "dspmv.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

#line 221 "dspmv.f"
    if (*beta != 1.) {
#line 222 "dspmv.f"
	if (*incy == 1) {
#line 223 "dspmv.f"
	    if (*beta == 0.) {
#line 224 "dspmv.f"
		i__1 = *n;
#line 224 "dspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 225 "dspmv.f"
		    y[i__] = 0.;
#line 226 "dspmv.f"
/* L10: */
#line 226 "dspmv.f"
		}
#line 227 "dspmv.f"
	    } else {
#line 228 "dspmv.f"
		i__1 = *n;
#line 228 "dspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "dspmv.f"
		    y[i__] = *beta * y[i__];
#line 230 "dspmv.f"
/* L20: */
#line 230 "dspmv.f"
		}
#line 231 "dspmv.f"
	    }
#line 232 "dspmv.f"
	} else {
#line 233 "dspmv.f"
	    iy = ky;
#line 234 "dspmv.f"
	    if (*beta == 0.) {
#line 235 "dspmv.f"
		i__1 = *n;
#line 235 "dspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "dspmv.f"
		    y[iy] = 0.;
#line 237 "dspmv.f"
		    iy += *incy;
#line 238 "dspmv.f"
/* L30: */
#line 238 "dspmv.f"
		}
#line 239 "dspmv.f"
	    } else {
#line 240 "dspmv.f"
		i__1 = *n;
#line 240 "dspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "dspmv.f"
		    y[iy] = *beta * y[iy];
#line 242 "dspmv.f"
		    iy += *incy;
#line 243 "dspmv.f"
/* L40: */
#line 243 "dspmv.f"
		}
#line 244 "dspmv.f"
	    }
#line 245 "dspmv.f"
	}
#line 246 "dspmv.f"
    }
#line 247 "dspmv.f"
    if (*alpha == 0.) {
#line 247 "dspmv.f"
	return 0;
#line 247 "dspmv.f"
    }
#line 248 "dspmv.f"
    kk = 1;
#line 249 "dspmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

#line 253 "dspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 254 "dspmv.f"
	    i__1 = *n;
#line 254 "dspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 255 "dspmv.f"
		temp1 = *alpha * x[j];
#line 256 "dspmv.f"
		temp2 = 0.;
#line 257 "dspmv.f"
		k = kk;
#line 258 "dspmv.f"
		i__2 = j - 1;
#line 258 "dspmv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 259 "dspmv.f"
		    y[i__] += temp1 * ap[k];
#line 260 "dspmv.f"
		    temp2 += ap[k] * x[i__];
#line 261 "dspmv.f"
		    ++k;
#line 262 "dspmv.f"
/* L50: */
#line 262 "dspmv.f"
		}
#line 263 "dspmv.f"
		y[j] = y[j] + temp1 * ap[kk + j - 1] + *alpha * temp2;
#line 264 "dspmv.f"
		kk += j;
#line 265 "dspmv.f"
/* L60: */
#line 265 "dspmv.f"
	    }
#line 266 "dspmv.f"
	} else {
#line 267 "dspmv.f"
	    jx = kx;
#line 268 "dspmv.f"
	    jy = ky;
#line 269 "dspmv.f"
	    i__1 = *n;
#line 269 "dspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 270 "dspmv.f"
		temp1 = *alpha * x[jx];
#line 271 "dspmv.f"
		temp2 = 0.;
#line 272 "dspmv.f"
		ix = kx;
#line 273 "dspmv.f"
		iy = ky;
#line 274 "dspmv.f"
		i__2 = kk + j - 2;
#line 274 "dspmv.f"
		for (k = kk; k <= i__2; ++k) {
#line 275 "dspmv.f"
		    y[iy] += temp1 * ap[k];
#line 276 "dspmv.f"
		    temp2 += ap[k] * x[ix];
#line 277 "dspmv.f"
		    ix += *incx;
#line 278 "dspmv.f"
		    iy += *incy;
#line 279 "dspmv.f"
/* L70: */
#line 279 "dspmv.f"
		}
#line 280 "dspmv.f"
		y[jy] = y[jy] + temp1 * ap[kk + j - 1] + *alpha * temp2;
#line 281 "dspmv.f"
		jx += *incx;
#line 282 "dspmv.f"
		jy += *incy;
#line 283 "dspmv.f"
		kk += j;
#line 284 "dspmv.f"
/* L80: */
#line 284 "dspmv.f"
	    }
#line 285 "dspmv.f"
	}
#line 286 "dspmv.f"
    } else {

/*        Form  y  when AP contains the lower triangle. */

#line 290 "dspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 291 "dspmv.f"
	    i__1 = *n;
#line 291 "dspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "dspmv.f"
		temp1 = *alpha * x[j];
#line 293 "dspmv.f"
		temp2 = 0.;
#line 294 "dspmv.f"
		y[j] += temp1 * ap[kk];
#line 295 "dspmv.f"
		k = kk + 1;
#line 296 "dspmv.f"
		i__2 = *n;
#line 296 "dspmv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 297 "dspmv.f"
		    y[i__] += temp1 * ap[k];
#line 298 "dspmv.f"
		    temp2 += ap[k] * x[i__];
#line 299 "dspmv.f"
		    ++k;
#line 300 "dspmv.f"
/* L90: */
#line 300 "dspmv.f"
		}
#line 301 "dspmv.f"
		y[j] += *alpha * temp2;
#line 302 "dspmv.f"
		kk += *n - j + 1;
#line 303 "dspmv.f"
/* L100: */
#line 303 "dspmv.f"
	    }
#line 304 "dspmv.f"
	} else {
#line 305 "dspmv.f"
	    jx = kx;
#line 306 "dspmv.f"
	    jy = ky;
#line 307 "dspmv.f"
	    i__1 = *n;
#line 307 "dspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 308 "dspmv.f"
		temp1 = *alpha * x[jx];
#line 309 "dspmv.f"
		temp2 = 0.;
#line 310 "dspmv.f"
		y[jy] += temp1 * ap[kk];
#line 311 "dspmv.f"
		ix = jx;
#line 312 "dspmv.f"
		iy = jy;
#line 313 "dspmv.f"
		i__2 = kk + *n - j;
#line 313 "dspmv.f"
		for (k = kk + 1; k <= i__2; ++k) {
#line 314 "dspmv.f"
		    ix += *incx;
#line 315 "dspmv.f"
		    iy += *incy;
#line 316 "dspmv.f"
		    y[iy] += temp1 * ap[k];
#line 317 "dspmv.f"
		    temp2 += ap[k] * x[ix];
#line 318 "dspmv.f"
/* L110: */
#line 318 "dspmv.f"
		}
#line 319 "dspmv.f"
		y[jy] += *alpha * temp2;
#line 320 "dspmv.f"
		jx += *incx;
#line 321 "dspmv.f"
		jy += *incy;
#line 322 "dspmv.f"
		kk += *n - j + 1;
#line 323 "dspmv.f"
/* L120: */
#line 323 "dspmv.f"
	    }
#line 324 "dspmv.f"
	}
#line 325 "dspmv.f"
    }

#line 327 "dspmv.f"
    return 0;

/*     End of DSPMV . */

} /* dspmv_ */

