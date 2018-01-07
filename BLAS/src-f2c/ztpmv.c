#line 1 "ztpmv.f"
/* ztpmv.f -- translated by f2c (version 20100827).
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

#line 1 "ztpmv.f"
/* > \brief \b ZTPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x,   or   x := A**H*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix, supplied in packed form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   x := A*x. */
/* > */
/* >              TRANS = 'T' or 't'   x := A**T*x. */
/* > */
/* >              TRANS = 'C' or 'c'   x := A**H*x. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit */
/* >           triangular as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular matrix packed sequentially, */
/* >           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/* >           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/* >           respectively, and so on. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular matrix packed sequentially, */
/* >           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/* >           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/* >           respectively, and so on. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced, but are assumed to be unity. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. On exit, X is overwritten with the */
/* >           transformed vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
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
/* Subroutine */ int ztpmv_(char *uplo, char *trans, char *diag, integer *n, 
	doublecomplex *ap, doublecomplex *x, integer *incx, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconj, nounit;


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

#line 182 "ztpmv.f"
    /* Parameter adjustments */
#line 182 "ztpmv.f"
    --x;
#line 182 "ztpmv.f"
    --ap;
#line 182 "ztpmv.f"

#line 182 "ztpmv.f"
    /* Function Body */
#line 182 "ztpmv.f"
    info = 0;
#line 183 "ztpmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 184 "ztpmv.f"
	info = 1;
#line 185 "ztpmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 187 "ztpmv.f"
	info = 2;
#line 188 "ztpmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 189 "ztpmv.f"
	info = 3;
#line 190 "ztpmv.f"
    } else if (*n < 0) {
#line 191 "ztpmv.f"
	info = 4;
#line 192 "ztpmv.f"
    } else if (*incx == 0) {
#line 193 "ztpmv.f"
	info = 7;
#line 194 "ztpmv.f"
    }
#line 195 "ztpmv.f"
    if (info != 0) {
#line 196 "ztpmv.f"
	xerbla_("ZTPMV ", &info, (ftnlen)6);
#line 197 "ztpmv.f"
	return 0;
#line 198 "ztpmv.f"
    }

/*     Quick return if possible. */

#line 202 "ztpmv.f"
    if (*n == 0) {
#line 202 "ztpmv.f"
	return 0;
#line 202 "ztpmv.f"
    }

#line 204 "ztpmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 205 "ztpmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 210 "ztpmv.f"
    if (*incx <= 0) {
#line 211 "ztpmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 212 "ztpmv.f"
    } else if (*incx != 1) {
#line 213 "ztpmv.f"
	kx = 1;
#line 214 "ztpmv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 219 "ztpmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x:= A*x. */

#line 223 "ztpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 224 "ztpmv.f"
	    kk = 1;
#line 225 "ztpmv.f"
	    if (*incx == 1) {
#line 226 "ztpmv.f"
		i__1 = *n;
#line 226 "ztpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 227 "ztpmv.f"
		    i__2 = j;
#line 227 "ztpmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 228 "ztpmv.f"
			i__2 = j;
#line 228 "ztpmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 229 "ztpmv.f"
			k = kk;
#line 230 "ztpmv.f"
			i__2 = j - 1;
#line 230 "ztpmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 231 "ztpmv.f"
			    i__3 = i__;
#line 231 "ztpmv.f"
			    i__4 = i__;
#line 231 "ztpmv.f"
			    i__5 = k;
#line 231 "ztpmv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 231 "ztpmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 231 "ztpmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 232 "ztpmv.f"
			    ++k;
#line 233 "ztpmv.f"
/* L10: */
#line 233 "ztpmv.f"
			}
#line 234 "ztpmv.f"
			if (nounit) {
#line 234 "ztpmv.f"
			    i__2 = j;
#line 234 "ztpmv.f"
			    i__3 = j;
#line 234 "ztpmv.f"
			    i__4 = kk + j - 1;
#line 234 "ztpmv.f"
			    z__1.r = x[i__3].r * ap[i__4].r - x[i__3].i * ap[
				    i__4].i, z__1.i = x[i__3].r * ap[i__4].i 
				    + x[i__3].i * ap[i__4].r;
#line 234 "ztpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 234 "ztpmv.f"
			}
#line 235 "ztpmv.f"
		    }
#line 236 "ztpmv.f"
		    kk += j;
#line 237 "ztpmv.f"
/* L20: */
#line 237 "ztpmv.f"
		}
#line 238 "ztpmv.f"
	    } else {
#line 239 "ztpmv.f"
		jx = kx;
#line 240 "ztpmv.f"
		i__1 = *n;
#line 240 "ztpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 241 "ztpmv.f"
		    i__2 = jx;
#line 241 "ztpmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 242 "ztpmv.f"
			i__2 = jx;
#line 242 "ztpmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 243 "ztpmv.f"
			ix = kx;
#line 244 "ztpmv.f"
			i__2 = kk + j - 2;
#line 244 "ztpmv.f"
			for (k = kk; k <= i__2; ++k) {
#line 245 "ztpmv.f"
			    i__3 = ix;
#line 245 "ztpmv.f"
			    i__4 = ix;
#line 245 "ztpmv.f"
			    i__5 = k;
#line 245 "ztpmv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 245 "ztpmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 245 "ztpmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 246 "ztpmv.f"
			    ix += *incx;
#line 247 "ztpmv.f"
/* L30: */
#line 247 "ztpmv.f"
			}
#line 248 "ztpmv.f"
			if (nounit) {
#line 248 "ztpmv.f"
			    i__2 = jx;
#line 248 "ztpmv.f"
			    i__3 = jx;
#line 248 "ztpmv.f"
			    i__4 = kk + j - 1;
#line 248 "ztpmv.f"
			    z__1.r = x[i__3].r * ap[i__4].r - x[i__3].i * ap[
				    i__4].i, z__1.i = x[i__3].r * ap[i__4].i 
				    + x[i__3].i * ap[i__4].r;
#line 248 "ztpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 248 "ztpmv.f"
			}
#line 249 "ztpmv.f"
		    }
#line 250 "ztpmv.f"
		    jx += *incx;
#line 251 "ztpmv.f"
		    kk += j;
#line 252 "ztpmv.f"
/* L40: */
#line 252 "ztpmv.f"
		}
#line 253 "ztpmv.f"
	    }
#line 254 "ztpmv.f"
	} else {
#line 255 "ztpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 256 "ztpmv.f"
	    if (*incx == 1) {
#line 257 "ztpmv.f"
		for (j = *n; j >= 1; --j) {
#line 258 "ztpmv.f"
		    i__1 = j;
#line 258 "ztpmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 259 "ztpmv.f"
			i__1 = j;
#line 259 "ztpmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 260 "ztpmv.f"
			k = kk;
#line 261 "ztpmv.f"
			i__1 = j + 1;
#line 261 "ztpmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 262 "ztpmv.f"
			    i__2 = i__;
#line 262 "ztpmv.f"
			    i__3 = i__;
#line 262 "ztpmv.f"
			    i__4 = k;
#line 262 "ztpmv.f"
			    z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4]
				    .i, z__2.i = temp.r * ap[i__4].i + temp.i 
				    * ap[i__4].r;
#line 262 "ztpmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 262 "ztpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 263 "ztpmv.f"
			    --k;
#line 264 "ztpmv.f"
/* L50: */
#line 264 "ztpmv.f"
			}
#line 265 "ztpmv.f"
			if (nounit) {
#line 265 "ztpmv.f"
			    i__1 = j;
#line 265 "ztpmv.f"
			    i__2 = j;
#line 265 "ztpmv.f"
			    i__3 = kk - *n + j;
#line 265 "ztpmv.f"
			    z__1.r = x[i__2].r * ap[i__3].r - x[i__2].i * ap[
				    i__3].i, z__1.i = x[i__2].r * ap[i__3].i 
				    + x[i__2].i * ap[i__3].r;
#line 265 "ztpmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 265 "ztpmv.f"
			}
#line 266 "ztpmv.f"
		    }
#line 267 "ztpmv.f"
		    kk -= *n - j + 1;
#line 268 "ztpmv.f"
/* L60: */
#line 268 "ztpmv.f"
		}
#line 269 "ztpmv.f"
	    } else {
#line 270 "ztpmv.f"
		kx += (*n - 1) * *incx;
#line 271 "ztpmv.f"
		jx = kx;
#line 272 "ztpmv.f"
		for (j = *n; j >= 1; --j) {
#line 273 "ztpmv.f"
		    i__1 = jx;
#line 273 "ztpmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 274 "ztpmv.f"
			i__1 = jx;
#line 274 "ztpmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 275 "ztpmv.f"
			ix = kx;
#line 276 "ztpmv.f"
			i__1 = kk - (*n - (j + 1));
#line 276 "ztpmv.f"
			for (k = kk; k >= i__1; --k) {
#line 277 "ztpmv.f"
			    i__2 = ix;
#line 277 "ztpmv.f"
			    i__3 = ix;
#line 277 "ztpmv.f"
			    i__4 = k;
#line 277 "ztpmv.f"
			    z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4]
				    .i, z__2.i = temp.r * ap[i__4].i + temp.i 
				    * ap[i__4].r;
#line 277 "ztpmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 277 "ztpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 278 "ztpmv.f"
			    ix -= *incx;
#line 279 "ztpmv.f"
/* L70: */
#line 279 "ztpmv.f"
			}
#line 280 "ztpmv.f"
			if (nounit) {
#line 280 "ztpmv.f"
			    i__1 = jx;
#line 280 "ztpmv.f"
			    i__2 = jx;
#line 280 "ztpmv.f"
			    i__3 = kk - *n + j;
#line 280 "ztpmv.f"
			    z__1.r = x[i__2].r * ap[i__3].r - x[i__2].i * ap[
				    i__3].i, z__1.i = x[i__2].r * ap[i__3].i 
				    + x[i__2].i * ap[i__3].r;
#line 280 "ztpmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 280 "ztpmv.f"
			}
#line 281 "ztpmv.f"
		    }
#line 282 "ztpmv.f"
		    jx -= *incx;
#line 283 "ztpmv.f"
		    kk -= *n - j + 1;
#line 284 "ztpmv.f"
/* L80: */
#line 284 "ztpmv.f"
		}
#line 285 "ztpmv.f"
	    }
#line 286 "ztpmv.f"
	}
#line 287 "ztpmv.f"
    } else {

/*        Form  x := A**T*x  or  x := A**H*x. */

#line 291 "ztpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 292 "ztpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 293 "ztpmv.f"
	    if (*incx == 1) {
#line 294 "ztpmv.f"
		for (j = *n; j >= 1; --j) {
#line 295 "ztpmv.f"
		    i__1 = j;
#line 295 "ztpmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 296 "ztpmv.f"
		    k = kk - 1;
#line 297 "ztpmv.f"
		    if (noconj) {
#line 298 "ztpmv.f"
			if (nounit) {
#line 298 "ztpmv.f"
			    i__1 = kk;
#line 298 "ztpmv.f"
			    z__1.r = temp.r * ap[i__1].r - temp.i * ap[i__1]
				    .i, z__1.i = temp.r * ap[i__1].i + temp.i 
				    * ap[i__1].r;
#line 298 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 298 "ztpmv.f"
			}
#line 299 "ztpmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 300 "ztpmv.f"
			    i__1 = k;
#line 300 "ztpmv.f"
			    i__2 = i__;
#line 300 "ztpmv.f"
			    z__2.r = ap[i__1].r * x[i__2].r - ap[i__1].i * x[
				    i__2].i, z__2.i = ap[i__1].r * x[i__2].i 
				    + ap[i__1].i * x[i__2].r;
#line 300 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 300 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ztpmv.f"
			    --k;
#line 302 "ztpmv.f"
/* L90: */
#line 302 "ztpmv.f"
			}
#line 303 "ztpmv.f"
		    } else {
#line 304 "ztpmv.f"
			if (nounit) {
#line 304 "ztpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 304 "ztpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 304 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 304 "ztpmv.f"
			}
#line 305 "ztpmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 306 "ztpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 306 "ztpmv.f"
			    i__1 = i__;
#line 306 "ztpmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 306 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 306 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 307 "ztpmv.f"
			    --k;
#line 308 "ztpmv.f"
/* L100: */
#line 308 "ztpmv.f"
			}
#line 309 "ztpmv.f"
		    }
#line 310 "ztpmv.f"
		    i__1 = j;
#line 310 "ztpmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 311 "ztpmv.f"
		    kk -= j;
#line 312 "ztpmv.f"
/* L110: */
#line 312 "ztpmv.f"
		}
#line 313 "ztpmv.f"
	    } else {
#line 314 "ztpmv.f"
		jx = kx + (*n - 1) * *incx;
#line 315 "ztpmv.f"
		for (j = *n; j >= 1; --j) {
#line 316 "ztpmv.f"
		    i__1 = jx;
#line 316 "ztpmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 317 "ztpmv.f"
		    ix = jx;
#line 318 "ztpmv.f"
		    if (noconj) {
#line 319 "ztpmv.f"
			if (nounit) {
#line 319 "ztpmv.f"
			    i__1 = kk;
#line 319 "ztpmv.f"
			    z__1.r = temp.r * ap[i__1].r - temp.i * ap[i__1]
				    .i, z__1.i = temp.r * ap[i__1].i + temp.i 
				    * ap[i__1].r;
#line 319 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 319 "ztpmv.f"
			}
#line 320 "ztpmv.f"
			i__1 = kk - j + 1;
#line 320 "ztpmv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 321 "ztpmv.f"
			    ix -= *incx;
#line 322 "ztpmv.f"
			    i__2 = k;
#line 322 "ztpmv.f"
			    i__3 = ix;
#line 322 "ztpmv.f"
			    z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[
				    i__3].i, z__2.i = ap[i__2].r * x[i__3].i 
				    + ap[i__2].i * x[i__3].r;
#line 322 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 322 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 323 "ztpmv.f"
/* L120: */
#line 323 "ztpmv.f"
			}
#line 324 "ztpmv.f"
		    } else {
#line 325 "ztpmv.f"
			if (nounit) {
#line 325 "ztpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 325 "ztpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 325 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 325 "ztpmv.f"
			}
#line 326 "ztpmv.f"
			i__1 = kk - j + 1;
#line 326 "ztpmv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 327 "ztpmv.f"
			    ix -= *incx;
#line 328 "ztpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 328 "ztpmv.f"
			    i__2 = ix;
#line 328 "ztpmv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 328 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 328 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 329 "ztpmv.f"
/* L130: */
#line 329 "ztpmv.f"
			}
#line 330 "ztpmv.f"
		    }
#line 331 "ztpmv.f"
		    i__1 = jx;
#line 331 "ztpmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 332 "ztpmv.f"
		    jx -= *incx;
#line 333 "ztpmv.f"
		    kk -= j;
#line 334 "ztpmv.f"
/* L140: */
#line 334 "ztpmv.f"
		}
#line 335 "ztpmv.f"
	    }
#line 336 "ztpmv.f"
	} else {
#line 337 "ztpmv.f"
	    kk = 1;
#line 338 "ztpmv.f"
	    if (*incx == 1) {
#line 339 "ztpmv.f"
		i__1 = *n;
#line 339 "ztpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 340 "ztpmv.f"
		    i__2 = j;
#line 340 "ztpmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 341 "ztpmv.f"
		    k = kk + 1;
#line 342 "ztpmv.f"
		    if (noconj) {
#line 343 "ztpmv.f"
			if (nounit) {
#line 343 "ztpmv.f"
			    i__2 = kk;
#line 343 "ztpmv.f"
			    z__1.r = temp.r * ap[i__2].r - temp.i * ap[i__2]
				    .i, z__1.i = temp.r * ap[i__2].i + temp.i 
				    * ap[i__2].r;
#line 343 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 343 "ztpmv.f"
			}
#line 344 "ztpmv.f"
			i__2 = *n;
#line 344 "ztpmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 345 "ztpmv.f"
			    i__3 = k;
#line 345 "ztpmv.f"
			    i__4 = i__;
#line 345 "ztpmv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 345 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 345 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 346 "ztpmv.f"
			    ++k;
#line 347 "ztpmv.f"
/* L150: */
#line 347 "ztpmv.f"
			}
#line 348 "ztpmv.f"
		    } else {
#line 349 "ztpmv.f"
			if (nounit) {
#line 349 "ztpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 349 "ztpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 349 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 349 "ztpmv.f"
			}
#line 350 "ztpmv.f"
			i__2 = *n;
#line 350 "ztpmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 351 "ztpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 351 "ztpmv.f"
			    i__3 = i__;
#line 351 "ztpmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 351 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 351 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 352 "ztpmv.f"
			    ++k;
#line 353 "ztpmv.f"
/* L160: */
#line 353 "ztpmv.f"
			}
#line 354 "ztpmv.f"
		    }
#line 355 "ztpmv.f"
		    i__2 = j;
#line 355 "ztpmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 356 "ztpmv.f"
		    kk += *n - j + 1;
#line 357 "ztpmv.f"
/* L170: */
#line 357 "ztpmv.f"
		}
#line 358 "ztpmv.f"
	    } else {
#line 359 "ztpmv.f"
		jx = kx;
#line 360 "ztpmv.f"
		i__1 = *n;
#line 360 "ztpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 361 "ztpmv.f"
		    i__2 = jx;
#line 361 "ztpmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 362 "ztpmv.f"
		    ix = jx;
#line 363 "ztpmv.f"
		    if (noconj) {
#line 364 "ztpmv.f"
			if (nounit) {
#line 364 "ztpmv.f"
			    i__2 = kk;
#line 364 "ztpmv.f"
			    z__1.r = temp.r * ap[i__2].r - temp.i * ap[i__2]
				    .i, z__1.i = temp.r * ap[i__2].i + temp.i 
				    * ap[i__2].r;
#line 364 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 364 "ztpmv.f"
			}
#line 365 "ztpmv.f"
			i__2 = kk + *n - j;
#line 365 "ztpmv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 366 "ztpmv.f"
			    ix += *incx;
#line 367 "ztpmv.f"
			    i__3 = k;
#line 367 "ztpmv.f"
			    i__4 = ix;
#line 367 "ztpmv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 367 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 367 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 368 "ztpmv.f"
/* L180: */
#line 368 "ztpmv.f"
			}
#line 369 "ztpmv.f"
		    } else {
#line 370 "ztpmv.f"
			if (nounit) {
#line 370 "ztpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 370 "ztpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 370 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ztpmv.f"
			}
#line 371 "ztpmv.f"
			i__2 = kk + *n - j;
#line 371 "ztpmv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 372 "ztpmv.f"
			    ix += *incx;
#line 373 "ztpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 373 "ztpmv.f"
			    i__3 = ix;
#line 373 "ztpmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 373 "ztpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 373 "ztpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 374 "ztpmv.f"
/* L190: */
#line 374 "ztpmv.f"
			}
#line 375 "ztpmv.f"
		    }
#line 376 "ztpmv.f"
		    i__2 = jx;
#line 376 "ztpmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 377 "ztpmv.f"
		    jx += *incx;
#line 378 "ztpmv.f"
		    kk += *n - j + 1;
#line 379 "ztpmv.f"
/* L200: */
#line 379 "ztpmv.f"
		}
#line 380 "ztpmv.f"
	    }
#line 381 "ztpmv.f"
	}
#line 382 "ztpmv.f"
    }

#line 384 "ztpmv.f"
    return 0;

/*     End of ZTPMV . */

} /* ztpmv_ */

