#line 1 "chpr.f"
/* chpr.f -- translated by f2c (version 20100827).
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

#line 1 "chpr.f"
/* > \brief \b CHPR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER INCX,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPR    performs the hermitian rank 1 operation */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array of dimension at least */
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
/* >          AP is COMPLEX array of DIMENSION at least */
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

/* > \ingroup complex_blas_level2 */

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
/* Subroutine */ int chpr_(char *uplo, integer *n, doublereal *alpha, 
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

#line 170 "chpr.f"
    /* Parameter adjustments */
#line 170 "chpr.f"
    --ap;
#line 170 "chpr.f"
    --x;
#line 170 "chpr.f"

#line 170 "chpr.f"
    /* Function Body */
#line 170 "chpr.f"
    info = 0;
#line 171 "chpr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 172 "chpr.f"
	info = 1;
#line 173 "chpr.f"
    } else if (*n < 0) {
#line 174 "chpr.f"
	info = 2;
#line 175 "chpr.f"
    } else if (*incx == 0) {
#line 176 "chpr.f"
	info = 5;
#line 177 "chpr.f"
    }
#line 178 "chpr.f"
    if (info != 0) {
#line 179 "chpr.f"
	xerbla_("CHPR  ", &info, (ftnlen)6);
#line 180 "chpr.f"
	return 0;
#line 181 "chpr.f"
    }

/*     Quick return if possible. */

#line 185 "chpr.f"
    if (*n == 0 || *alpha == 0.) {
#line 185 "chpr.f"
	return 0;
#line 185 "chpr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 189 "chpr.f"
    if (*incx <= 0) {
#line 190 "chpr.f"
	kx = 1 - (*n - 1) * *incx;
#line 191 "chpr.f"
    } else if (*incx != 1) {
#line 192 "chpr.f"
	kx = 1;
#line 193 "chpr.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 198 "chpr.f"
    kk = 1;
#line 199 "chpr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 203 "chpr.f"
	if (*incx == 1) {
#line 204 "chpr.f"
	    i__1 = *n;
#line 204 "chpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "chpr.f"
		i__2 = j;
#line 205 "chpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 206 "chpr.f"
		    d_cnjg(&z__2, &x[j]);
#line 206 "chpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 206 "chpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 207 "chpr.f"
		    k = kk;
#line 208 "chpr.f"
		    i__2 = j - 1;
#line 208 "chpr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 209 "chpr.f"
			i__3 = k;
#line 209 "chpr.f"
			i__4 = k;
#line 209 "chpr.f"
			i__5 = i__;
#line 209 "chpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 209 "chpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 209 "chpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 210 "chpr.f"
			++k;
#line 211 "chpr.f"
/* L10: */
#line 211 "chpr.f"
		    }
#line 212 "chpr.f"
		    i__2 = kk + j - 1;
#line 212 "chpr.f"
		    i__3 = kk + j - 1;
#line 212 "chpr.f"
		    i__4 = j;
#line 212 "chpr.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 212 "chpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 212 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 213 "chpr.f"
		} else {
#line 214 "chpr.f"
		    i__2 = kk + j - 1;
#line 214 "chpr.f"
		    i__3 = kk + j - 1;
#line 214 "chpr.f"
		    d__1 = ap[i__3].r;
#line 214 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 215 "chpr.f"
		}
#line 216 "chpr.f"
		kk += j;
#line 217 "chpr.f"
/* L20: */
#line 217 "chpr.f"
	    }
#line 218 "chpr.f"
	} else {
#line 219 "chpr.f"
	    jx = kx;
#line 220 "chpr.f"
	    i__1 = *n;
#line 220 "chpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 221 "chpr.f"
		i__2 = jx;
#line 221 "chpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 222 "chpr.f"
		    d_cnjg(&z__2, &x[jx]);
#line 222 "chpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 222 "chpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 223 "chpr.f"
		    ix = kx;
#line 224 "chpr.f"
		    i__2 = kk + j - 2;
#line 224 "chpr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 225 "chpr.f"
			i__3 = k;
#line 225 "chpr.f"
			i__4 = k;
#line 225 "chpr.f"
			i__5 = ix;
#line 225 "chpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 225 "chpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 225 "chpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 226 "chpr.f"
			ix += *incx;
#line 227 "chpr.f"
/* L30: */
#line 227 "chpr.f"
		    }
#line 228 "chpr.f"
		    i__2 = kk + j - 1;
#line 228 "chpr.f"
		    i__3 = kk + j - 1;
#line 228 "chpr.f"
		    i__4 = jx;
#line 228 "chpr.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 228 "chpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 228 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 229 "chpr.f"
		} else {
#line 230 "chpr.f"
		    i__2 = kk + j - 1;
#line 230 "chpr.f"
		    i__3 = kk + j - 1;
#line 230 "chpr.f"
		    d__1 = ap[i__3].r;
#line 230 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 231 "chpr.f"
		}
#line 232 "chpr.f"
		jx += *incx;
#line 233 "chpr.f"
		kk += j;
#line 234 "chpr.f"
/* L40: */
#line 234 "chpr.f"
	    }
#line 235 "chpr.f"
	}
#line 236 "chpr.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 240 "chpr.f"
	if (*incx == 1) {
#line 241 "chpr.f"
	    i__1 = *n;
#line 241 "chpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 242 "chpr.f"
		i__2 = j;
#line 242 "chpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 243 "chpr.f"
		    d_cnjg(&z__2, &x[j]);
#line 243 "chpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 243 "chpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 244 "chpr.f"
		    i__2 = kk;
#line 244 "chpr.f"
		    i__3 = kk;
#line 244 "chpr.f"
		    i__4 = j;
#line 244 "chpr.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 244 "chpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 244 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 245 "chpr.f"
		    k = kk + 1;
#line 246 "chpr.f"
		    i__2 = *n;
#line 246 "chpr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 247 "chpr.f"
			i__3 = k;
#line 247 "chpr.f"
			i__4 = k;
#line 247 "chpr.f"
			i__5 = i__;
#line 247 "chpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 247 "chpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 247 "chpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 248 "chpr.f"
			++k;
#line 249 "chpr.f"
/* L50: */
#line 249 "chpr.f"
		    }
#line 250 "chpr.f"
		} else {
#line 251 "chpr.f"
		    i__2 = kk;
#line 251 "chpr.f"
		    i__3 = kk;
#line 251 "chpr.f"
		    d__1 = ap[i__3].r;
#line 251 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 252 "chpr.f"
		}
#line 253 "chpr.f"
		kk = kk + *n - j + 1;
#line 254 "chpr.f"
/* L60: */
#line 254 "chpr.f"
	    }
#line 255 "chpr.f"
	} else {
#line 256 "chpr.f"
	    jx = kx;
#line 257 "chpr.f"
	    i__1 = *n;
#line 257 "chpr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "chpr.f"
		i__2 = jx;
#line 258 "chpr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 259 "chpr.f"
		    d_cnjg(&z__2, &x[jx]);
#line 259 "chpr.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 259 "chpr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 260 "chpr.f"
		    i__2 = kk;
#line 260 "chpr.f"
		    i__3 = kk;
#line 260 "chpr.f"
		    i__4 = jx;
#line 260 "chpr.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 260 "chpr.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 260 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 261 "chpr.f"
		    ix = jx;
#line 262 "chpr.f"
		    i__2 = kk + *n - j;
#line 262 "chpr.f"
		    for (k = kk + 1; k <= i__2; ++k) {
#line 263 "chpr.f"
			ix += *incx;
#line 264 "chpr.f"
			i__3 = k;
#line 264 "chpr.f"
			i__4 = k;
#line 264 "chpr.f"
			i__5 = ix;
#line 264 "chpr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 264 "chpr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 264 "chpr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 265 "chpr.f"
/* L70: */
#line 265 "chpr.f"
		    }
#line 266 "chpr.f"
		} else {
#line 267 "chpr.f"
		    i__2 = kk;
#line 267 "chpr.f"
		    i__3 = kk;
#line 267 "chpr.f"
		    d__1 = ap[i__3].r;
#line 267 "chpr.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 268 "chpr.f"
		}
#line 269 "chpr.f"
		jx += *incx;
#line 270 "chpr.f"
		kk = kk + *n - j + 1;
#line 271 "chpr.f"
/* L80: */
#line 271 "chpr.f"
	    }
#line 272 "chpr.f"
	}
#line 273 "chpr.f"
    }

#line 275 "chpr.f"
    return 0;

/*     End of CHPR  . */

} /* chpr_ */

