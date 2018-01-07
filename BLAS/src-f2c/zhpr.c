#line 1 "zhpr.f"
/* zhpr.f -- translated by f2c (version 20100827).
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

#line 1 "zhpr.f"
/* > \brief \b ZHPR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER INCX,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPR    performs the hermitian rank 1 operation */
/* > */
/* >    A := alpha*x*x**H + A, */
/* > */
/* > where alpha is a real scalar, x is an n element vector and A is an */
/* > n by n hermitian matrix, supplied in packed form. */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of dimension at least */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the hermitian matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the upper triangular part of the */
/* >           updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the hermitian matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the lower triangular part of the */
/* >           updated matrix. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set, they are assumed to be zero, and on exit they */
/* >           are set to zero. */
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
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhpr_(char *uplo, integer *n, doublereal *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static doublecomplex temp;
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

#line 170 "zhpr.f"
    /* Parameter adjustments */
#line 170 "zhpr.f"
    --ap;
#line 170 "zhpr.f"
    --x;
#line 170 "zhpr.f"

#line 170 "zhpr.f"
    /* Function Body */
#line 170 "zhpr.f"
    info = 0;
#line 171 "zhpr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 172 "zhpr.f"
	info = 1;
#line 173 "zhpr.f"
    } else if (*n < 0) {
#line 174 "zhpr.f"
	info = 2;
#line 175 "zhpr.f"
    } else if (*incx == 0) {
#line 176 "zhpr.f"
	info = 5;
#line 177 "zhpr.f"
    }
#line 178 "zhpr.f"
    if (info != 0) {
#line 179 "zhpr.f"
	xerbla_("ZHPR  ", &info, (ftnlen)6);
#line 180 "zhpr.f"
	return 0;
#line 181 "zhpr.f"
    }

/*     Quick return if possible. */

#line 185 "zhpr.f"
    if (*n == 0 || *alpha == 0.) {
#line 185 "zhpr.f"
	return 0;
#line 185 "zhpr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 189 "zhpr.f"
    if (*incx <= 0) {
#line 190 "zhpr.f"
	kx = 1 - (*n - 1) * *incx;
#line 191 "zhpr.f"
    } else if (*incx != 1) {
#line 192 "zhpr.f"
	kx = 1;
#line 193 "zhpr.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 198 "zhpr.f"
    kk = 1;
#line 199 "zhpr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 203 "zhpr.f"
	if (*incx == 1) {
#line 204 "zhpr.f"
	    i__1 = *n;
#line 204 "zhpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "zhpr.f"
		i__2 = j;
#line 205 "zhpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 206 "zhpr.f"
		    d_cnjg(&z__2, &x[j]);
#line 206 "zhpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 206 "zhpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 207 "zhpr.f"
		    k = kk;
#line 208 "zhpr.f"
		    i__2 = j - 1;
#line 208 "zhpr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 209 "zhpr.f"
			i__3 = k;
#line 209 "zhpr.f"
			i__4 = k;
#line 209 "zhpr.f"
			i__5 = i__;
#line 209 "zhpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 209 "zhpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 209 "zhpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 210 "zhpr.f"
			++k;
#line 211 "zhpr.f"
/* L10: */
#line 211 "zhpr.f"
		    }
#line 212 "zhpr.f"
		    i__2 = kk + j - 1;
#line 212 "zhpr.f"
		    i__3 = kk + j - 1;
#line 212 "zhpr.f"
		    i__4 = j;
#line 212 "zhpr.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 212 "zhpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 212 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 213 "zhpr.f"
		} else {
#line 214 "zhpr.f"
		    i__2 = kk + j - 1;
#line 214 "zhpr.f"
		    i__3 = kk + j - 1;
#line 214 "zhpr.f"
		    d__1 = ap[i__3].r;
#line 214 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 215 "zhpr.f"
		}
#line 216 "zhpr.f"
		kk += j;
#line 217 "zhpr.f"
/* L20: */
#line 217 "zhpr.f"
	    }
#line 218 "zhpr.f"
	} else {
#line 219 "zhpr.f"
	    jx = kx;
#line 220 "zhpr.f"
	    i__1 = *n;
#line 220 "zhpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 221 "zhpr.f"
		i__2 = jx;
#line 221 "zhpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 222 "zhpr.f"
		    d_cnjg(&z__2, &x[jx]);
#line 222 "zhpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 222 "zhpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 223 "zhpr.f"
		    ix = kx;
#line 224 "zhpr.f"
		    i__2 = kk + j - 2;
#line 224 "zhpr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 225 "zhpr.f"
			i__3 = k;
#line 225 "zhpr.f"
			i__4 = k;
#line 225 "zhpr.f"
			i__5 = ix;
#line 225 "zhpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 225 "zhpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 225 "zhpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 226 "zhpr.f"
			ix += *incx;
#line 227 "zhpr.f"
/* L30: */
#line 227 "zhpr.f"
		    }
#line 228 "zhpr.f"
		    i__2 = kk + j - 1;
#line 228 "zhpr.f"
		    i__3 = kk + j - 1;
#line 228 "zhpr.f"
		    i__4 = jx;
#line 228 "zhpr.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 228 "zhpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 228 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 229 "zhpr.f"
		} else {
#line 230 "zhpr.f"
		    i__2 = kk + j - 1;
#line 230 "zhpr.f"
		    i__3 = kk + j - 1;
#line 230 "zhpr.f"
		    d__1 = ap[i__3].r;
#line 230 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 231 "zhpr.f"
		}
#line 232 "zhpr.f"
		jx += *incx;
#line 233 "zhpr.f"
		kk += j;
#line 234 "zhpr.f"
/* L40: */
#line 234 "zhpr.f"
	    }
#line 235 "zhpr.f"
	}
#line 236 "zhpr.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 240 "zhpr.f"
	if (*incx == 1) {
#line 241 "zhpr.f"
	    i__1 = *n;
#line 241 "zhpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 242 "zhpr.f"
		i__2 = j;
#line 242 "zhpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 243 "zhpr.f"
		    d_cnjg(&z__2, &x[j]);
#line 243 "zhpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 243 "zhpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 244 "zhpr.f"
		    i__2 = kk;
#line 244 "zhpr.f"
		    i__3 = kk;
#line 244 "zhpr.f"
		    i__4 = j;
#line 244 "zhpr.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 244 "zhpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 244 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 245 "zhpr.f"
		    k = kk + 1;
#line 246 "zhpr.f"
		    i__2 = *n;
#line 246 "zhpr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 247 "zhpr.f"
			i__3 = k;
#line 247 "zhpr.f"
			i__4 = k;
#line 247 "zhpr.f"
			i__5 = i__;
#line 247 "zhpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 247 "zhpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 247 "zhpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 248 "zhpr.f"
			++k;
#line 249 "zhpr.f"
/* L50: */
#line 249 "zhpr.f"
		    }
#line 250 "zhpr.f"
		} else {
#line 251 "zhpr.f"
		    i__2 = kk;
#line 251 "zhpr.f"
		    i__3 = kk;
#line 251 "zhpr.f"
		    d__1 = ap[i__3].r;
#line 251 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 252 "zhpr.f"
		}
#line 253 "zhpr.f"
		kk = kk + *n - j + 1;
#line 254 "zhpr.f"
/* L60: */
#line 254 "zhpr.f"
	    }
#line 255 "zhpr.f"
	} else {
#line 256 "zhpr.f"
	    jx = kx;
#line 257 "zhpr.f"
	    i__1 = *n;
#line 257 "zhpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "zhpr.f"
		i__2 = jx;
#line 258 "zhpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 259 "zhpr.f"
		    d_cnjg(&z__2, &x[jx]);
#line 259 "zhpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 259 "zhpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 260 "zhpr.f"
		    i__2 = kk;
#line 260 "zhpr.f"
		    i__3 = kk;
#line 260 "zhpr.f"
		    i__4 = jx;
#line 260 "zhpr.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 260 "zhpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 260 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 261 "zhpr.f"
		    ix = jx;
#line 262 "zhpr.f"
		    i__2 = kk + *n - j;
#line 262 "zhpr.f"
		    for (k = kk + 1; k <= i__2; ++k) {
#line 263 "zhpr.f"
			ix += *incx;
#line 264 "zhpr.f"
			i__3 = k;
#line 264 "zhpr.f"
			i__4 = k;
#line 264 "zhpr.f"
			i__5 = ix;
#line 264 "zhpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 264 "zhpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 264 "zhpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 265 "zhpr.f"
/* L70: */
#line 265 "zhpr.f"
		    }
#line 266 "zhpr.f"
		} else {
#line 267 "zhpr.f"
		    i__2 = kk;
#line 267 "zhpr.f"
		    i__3 = kk;
#line 267 "zhpr.f"
		    d__1 = ap[i__3].r;
#line 267 "zhpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 268 "zhpr.f"
		}
#line 269 "zhpr.f"
		jx += *incx;
#line 270 "zhpr.f"
		kk = kk + *n - j + 1;
#line 271 "zhpr.f"
/* L80: */
#line 271 "zhpr.f"
	    }
#line 272 "zhpr.f"
	}
#line 273 "zhpr.f"
    }

#line 275 "zhpr.f"
    return 0;

/*     End of ZHPR  . */

} /* zhpr_ */

