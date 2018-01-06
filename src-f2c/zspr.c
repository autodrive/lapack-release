#line 1 "zspr.f"
/* zspr.f -- translated by f2c (version 20100827).
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

#line 1 "zspr.f"
/* > \brief \b ZSPR performs the symmetrical rank-1 update of a complex symmetric packed matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSPR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspr.f"
> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspr.f"
> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspr.f"
> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSPR( UPLO, N, ALPHA, X, INCX, AP ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCX, N */
/*       COMPLEX*16         ALPHA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSPR    performs the symmetric rank 1 operation */
/* > */
/* >    A := alpha*x*x**H + A, */
/* > */
/* > where alpha is a complex scalar, x is an n element vector and A is an */
/* > n by n symmetric matrix, supplied in packed form. */
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
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the N- */
/* >           element vector x. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension at least */
/* >           ( ( N*( N + 1 ) )/2 ). */
/* >           Before entry, with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the upper triangular part of the */
/* >           updated matrix. */
/* >           Before entry, with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
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

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zspr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 170 "zspr.f"
    /* Parameter adjustments */
#line 170 "zspr.f"
    --ap;
#line 170 "zspr.f"
    --x;
#line 170 "zspr.f"

#line 170 "zspr.f"
    /* Function Body */
#line 170 "zspr.f"
    info = 0;
#line 171 "zspr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 172 "zspr.f"
	info = 1;
#line 173 "zspr.f"
    } else if (*n < 0) {
#line 174 "zspr.f"
	info = 2;
#line 175 "zspr.f"
    } else if (*incx == 0) {
#line 176 "zspr.f"
	info = 5;
#line 177 "zspr.f"
    }
#line 178 "zspr.f"
    if (info != 0) {
#line 179 "zspr.f"
	xerbla_("ZSPR  ", &info, (ftnlen)6);
#line 180 "zspr.f"
	return 0;
#line 181 "zspr.f"
    }

/*     Quick return if possible. */

#line 185 "zspr.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
#line 185 "zspr.f"
	return 0;
#line 185 "zspr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 190 "zspr.f"
    if (*incx <= 0) {
#line 191 "zspr.f"
	kx = 1 - (*n - 1) * *incx;
#line 192 "zspr.f"
    } else if (*incx != 1) {
#line 193 "zspr.f"
	kx = 1;
#line 194 "zspr.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 199 "zspr.f"
    kk = 1;
#line 200 "zspr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 204 "zspr.f"
	if (*incx == 1) {
#line 205 "zspr.f"
	    i__1 = *n;
#line 205 "zspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 206 "zspr.f"
		i__2 = j;
#line 206 "zspr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 207 "zspr.f"
		    i__2 = j;
#line 207 "zspr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 207 "zspr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 208 "zspr.f"
		    k = kk;
#line 209 "zspr.f"
		    i__2 = j - 1;
#line 209 "zspr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 210 "zspr.f"
			i__3 = k;
#line 210 "zspr.f"
			i__4 = k;
#line 210 "zspr.f"
			i__5 = i__;
#line 210 "zspr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 210 "zspr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 210 "zspr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 211 "zspr.f"
			++k;
#line 212 "zspr.f"
/* L10: */
#line 212 "zspr.f"
		    }
#line 213 "zspr.f"
		    i__2 = kk + j - 1;
#line 213 "zspr.f"
		    i__3 = kk + j - 1;
#line 213 "zspr.f"
		    i__4 = j;
#line 213 "zspr.f"
		    z__2.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__2.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 213 "zspr.f"
		    z__1.r = ap[i__3].r + z__2.r, z__1.i = ap[i__3].i + 
			    z__2.i;
#line 213 "zspr.f"
		    ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 214 "zspr.f"
		} else {
#line 215 "zspr.f"
		    i__2 = kk + j - 1;
#line 215 "zspr.f"
		    i__3 = kk + j - 1;
#line 215 "zspr.f"
		    ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 216 "zspr.f"
		}
#line 217 "zspr.f"
		kk += j;
#line 218 "zspr.f"
/* L20: */
#line 218 "zspr.f"
	    }
#line 219 "zspr.f"
	} else {
#line 220 "zspr.f"
	    jx = kx;
#line 221 "zspr.f"
	    i__1 = *n;
#line 221 "zspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 222 "zspr.f"
		i__2 = jx;
#line 222 "zspr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 223 "zspr.f"
		    i__2 = jx;
#line 223 "zspr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 223 "zspr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 224 "zspr.f"
		    ix = kx;
#line 225 "zspr.f"
		    i__2 = kk + j - 2;
#line 225 "zspr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 226 "zspr.f"
			i__3 = k;
#line 226 "zspr.f"
			i__4 = k;
#line 226 "zspr.f"
			i__5 = ix;
#line 226 "zspr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 226 "zspr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 226 "zspr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 227 "zspr.f"
			ix += *incx;
#line 228 "zspr.f"
/* L30: */
#line 228 "zspr.f"
		    }
#line 229 "zspr.f"
		    i__2 = kk + j - 1;
#line 229 "zspr.f"
		    i__3 = kk + j - 1;
#line 229 "zspr.f"
		    i__4 = jx;
#line 229 "zspr.f"
		    z__2.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__2.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 229 "zspr.f"
		    z__1.r = ap[i__3].r + z__2.r, z__1.i = ap[i__3].i + 
			    z__2.i;
#line 229 "zspr.f"
		    ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 230 "zspr.f"
		} else {
#line 231 "zspr.f"
		    i__2 = kk + j - 1;
#line 231 "zspr.f"
		    i__3 = kk + j - 1;
#line 231 "zspr.f"
		    ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 232 "zspr.f"
		}
#line 233 "zspr.f"
		jx += *incx;
#line 234 "zspr.f"
		kk += j;
#line 235 "zspr.f"
/* L40: */
#line 235 "zspr.f"
	    }
#line 236 "zspr.f"
	}
#line 237 "zspr.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 241 "zspr.f"
	if (*incx == 1) {
#line 242 "zspr.f"
	    i__1 = *n;
#line 242 "zspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 243 "zspr.f"
		i__2 = j;
#line 243 "zspr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 244 "zspr.f"
		    i__2 = j;
#line 244 "zspr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 244 "zspr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 245 "zspr.f"
		    i__2 = kk;
#line 245 "zspr.f"
		    i__3 = kk;
#line 245 "zspr.f"
		    i__4 = j;
#line 245 "zspr.f"
		    z__2.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__2.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 245 "zspr.f"
		    z__1.r = ap[i__3].r + z__2.r, z__1.i = ap[i__3].i + 
			    z__2.i;
#line 245 "zspr.f"
		    ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 246 "zspr.f"
		    k = kk + 1;
#line 247 "zspr.f"
		    i__2 = *n;
#line 247 "zspr.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 248 "zspr.f"
			i__3 = k;
#line 248 "zspr.f"
			i__4 = k;
#line 248 "zspr.f"
			i__5 = i__;
#line 248 "zspr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 248 "zspr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 248 "zspr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 249 "zspr.f"
			++k;
#line 250 "zspr.f"
/* L50: */
#line 250 "zspr.f"
		    }
#line 251 "zspr.f"
		} else {
#line 252 "zspr.f"
		    i__2 = kk;
#line 252 "zspr.f"
		    i__3 = kk;
#line 252 "zspr.f"
		    ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 253 "zspr.f"
		}
#line 254 "zspr.f"
		kk = kk + *n - j + 1;
#line 255 "zspr.f"
/* L60: */
#line 255 "zspr.f"
	    }
#line 256 "zspr.f"
	} else {
#line 257 "zspr.f"
	    jx = kx;
#line 258 "zspr.f"
	    i__1 = *n;
#line 258 "zspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 259 "zspr.f"
		i__2 = jx;
#line 259 "zspr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 260 "zspr.f"
		    i__2 = jx;
#line 260 "zspr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 260 "zspr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 261 "zspr.f"
		    i__2 = kk;
#line 261 "zspr.f"
		    i__3 = kk;
#line 261 "zspr.f"
		    i__4 = jx;
#line 261 "zspr.f"
		    z__2.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__2.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 261 "zspr.f"
		    z__1.r = ap[i__3].r + z__2.r, z__1.i = ap[i__3].i + 
			    z__2.i;
#line 261 "zspr.f"
		    ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 262 "zspr.f"
		    ix = jx;
#line 263 "zspr.f"
		    i__2 = kk + *n - j;
#line 263 "zspr.f"
		    for (k = kk + 1; k <= i__2; ++k) {
#line 264 "zspr.f"
			ix += *incx;
#line 265 "zspr.f"
			i__3 = k;
#line 265 "zspr.f"
			i__4 = k;
#line 265 "zspr.f"
			i__5 = ix;
#line 265 "zspr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 265 "zspr.f"
			z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + 
				z__2.i;
#line 265 "zspr.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 266 "zspr.f"
/* L70: */
#line 266 "zspr.f"
		    }
#line 267 "zspr.f"
		} else {
#line 268 "zspr.f"
		    i__2 = kk;
#line 268 "zspr.f"
		    i__3 = kk;
#line 268 "zspr.f"
		    ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 269 "zspr.f"
		}
#line 270 "zspr.f"
		jx += *incx;
#line 271 "zspr.f"
		kk = kk + *n - j + 1;
#line 272 "zspr.f"
/* L80: */
#line 272 "zspr.f"
	    }
#line 273 "zspr.f"
	}
#line 274 "zspr.f"
    }

#line 276 "zspr.f"
    return 0;

/*     End of ZSPR */

} /* zspr_ */

