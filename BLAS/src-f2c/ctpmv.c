#line 1 "ctpmv.f"
/* ctpmv.f -- translated by f2c (version 20100827).
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

#line 1 "ctpmv.f"
/* > \brief \b CTPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPMV  performs one of the matrix-vector operations */
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
/* >          AP is COMPLEX array of DIMENSION at least */
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
/* >          X is COMPLEX array of dimension at least */
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
/* Subroutine */ int ctpmv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 182 "ctpmv.f"
    /* Parameter adjustments */
#line 182 "ctpmv.f"
    --x;
#line 182 "ctpmv.f"
    --ap;
#line 182 "ctpmv.f"

#line 182 "ctpmv.f"
    /* Function Body */
#line 182 "ctpmv.f"
    info = 0;
#line 183 "ctpmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 184 "ctpmv.f"
	info = 1;
#line 185 "ctpmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 187 "ctpmv.f"
	info = 2;
#line 188 "ctpmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 189 "ctpmv.f"
	info = 3;
#line 190 "ctpmv.f"
    } else if (*n < 0) {
#line 191 "ctpmv.f"
	info = 4;
#line 192 "ctpmv.f"
    } else if (*incx == 0) {
#line 193 "ctpmv.f"
	info = 7;
#line 194 "ctpmv.f"
    }
#line 195 "ctpmv.f"
    if (info != 0) {
#line 196 "ctpmv.f"
	xerbla_("CTPMV ", &info, (ftnlen)6);
#line 197 "ctpmv.f"
	return 0;
#line 198 "ctpmv.f"
    }

/*     Quick return if possible. */

#line 202 "ctpmv.f"
    if (*n == 0) {
#line 202 "ctpmv.f"
	return 0;
#line 202 "ctpmv.f"
    }

#line 204 "ctpmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 205 "ctpmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 210 "ctpmv.f"
    if (*incx <= 0) {
#line 211 "ctpmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 212 "ctpmv.f"
    } else if (*incx != 1) {
#line 213 "ctpmv.f"
	kx = 1;
#line 214 "ctpmv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 219 "ctpmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x:= A*x. */

#line 223 "ctpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 224 "ctpmv.f"
	    kk = 1;
#line 225 "ctpmv.f"
	    if (*incx == 1) {
#line 226 "ctpmv.f"
		i__1 = *n;
#line 226 "ctpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 227 "ctpmv.f"
		    i__2 = j;
#line 227 "ctpmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 228 "ctpmv.f"
			i__2 = j;
#line 228 "ctpmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 229 "ctpmv.f"
			k = kk;
#line 230 "ctpmv.f"
			i__2 = j - 1;
#line 230 "ctpmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 231 "ctpmv.f"
			    i__3 = i__;
#line 231 "ctpmv.f"
			    i__4 = i__;
#line 231 "ctpmv.f"
			    i__5 = k;
#line 231 "ctpmv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 231 "ctpmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 231 "ctpmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 232 "ctpmv.f"
			    ++k;
#line 233 "ctpmv.f"
/* L10: */
#line 233 "ctpmv.f"
			}
#line 234 "ctpmv.f"
			if (nounit) {
#line 234 "ctpmv.f"
			    i__2 = j;
#line 234 "ctpmv.f"
			    i__3 = j;
#line 234 "ctpmv.f"
			    i__4 = kk + j - 1;
#line 234 "ctpmv.f"
			    z__1.r = x[i__3].r * ap[i__4].r - x[i__3].i * ap[
				    i__4].i, z__1.i = x[i__3].r * ap[i__4].i 
				    + x[i__3].i * ap[i__4].r;
#line 234 "ctpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 234 "ctpmv.f"
			}
#line 235 "ctpmv.f"
		    }
#line 236 "ctpmv.f"
		    kk += j;
#line 237 "ctpmv.f"
/* L20: */
#line 237 "ctpmv.f"
		}
#line 238 "ctpmv.f"
	    } else {
#line 239 "ctpmv.f"
		jx = kx;
#line 240 "ctpmv.f"
		i__1 = *n;
#line 240 "ctpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 241 "ctpmv.f"
		    i__2 = jx;
#line 241 "ctpmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 242 "ctpmv.f"
			i__2 = jx;
#line 242 "ctpmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 243 "ctpmv.f"
			ix = kx;
#line 244 "ctpmv.f"
			i__2 = kk + j - 2;
#line 244 "ctpmv.f"
			for (k = kk; k <= i__2; ++k) {
#line 245 "ctpmv.f"
			    i__3 = ix;
#line 245 "ctpmv.f"
			    i__4 = ix;
#line 245 "ctpmv.f"
			    i__5 = k;
#line 245 "ctpmv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 245 "ctpmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 245 "ctpmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 246 "ctpmv.f"
			    ix += *incx;
#line 247 "ctpmv.f"
/* L30: */
#line 247 "ctpmv.f"
			}
#line 248 "ctpmv.f"
			if (nounit) {
#line 248 "ctpmv.f"
			    i__2 = jx;
#line 248 "ctpmv.f"
			    i__3 = jx;
#line 248 "ctpmv.f"
			    i__4 = kk + j - 1;
#line 248 "ctpmv.f"
			    z__1.r = x[i__3].r * ap[i__4].r - x[i__3].i * ap[
				    i__4].i, z__1.i = x[i__3].r * ap[i__4].i 
				    + x[i__3].i * ap[i__4].r;
#line 248 "ctpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 248 "ctpmv.f"
			}
#line 249 "ctpmv.f"
		    }
#line 250 "ctpmv.f"
		    jx += *incx;
#line 251 "ctpmv.f"
		    kk += j;
#line 252 "ctpmv.f"
/* L40: */
#line 252 "ctpmv.f"
		}
#line 253 "ctpmv.f"
	    }
#line 254 "ctpmv.f"
	} else {
#line 255 "ctpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 256 "ctpmv.f"
	    if (*incx == 1) {
#line 257 "ctpmv.f"
		for (j = *n; j >= 1; --j) {
#line 258 "ctpmv.f"
		    i__1 = j;
#line 258 "ctpmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 259 "ctpmv.f"
			i__1 = j;
#line 259 "ctpmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 260 "ctpmv.f"
			k = kk;
#line 261 "ctpmv.f"
			i__1 = j + 1;
#line 261 "ctpmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 262 "ctpmv.f"
			    i__2 = i__;
#line 262 "ctpmv.f"
			    i__3 = i__;
#line 262 "ctpmv.f"
			    i__4 = k;
#line 262 "ctpmv.f"
			    z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4]
				    .i, z__2.i = temp.r * ap[i__4].i + temp.i 
				    * ap[i__4].r;
#line 262 "ctpmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 262 "ctpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 263 "ctpmv.f"
			    --k;
#line 264 "ctpmv.f"
/* L50: */
#line 264 "ctpmv.f"
			}
#line 265 "ctpmv.f"
			if (nounit) {
#line 265 "ctpmv.f"
			    i__1 = j;
#line 265 "ctpmv.f"
			    i__2 = j;
#line 265 "ctpmv.f"
			    i__3 = kk - *n + j;
#line 265 "ctpmv.f"
			    z__1.r = x[i__2].r * ap[i__3].r - x[i__2].i * ap[
				    i__3].i, z__1.i = x[i__2].r * ap[i__3].i 
				    + x[i__2].i * ap[i__3].r;
#line 265 "ctpmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 265 "ctpmv.f"
			}
#line 266 "ctpmv.f"
		    }
#line 267 "ctpmv.f"
		    kk -= *n - j + 1;
#line 268 "ctpmv.f"
/* L60: */
#line 268 "ctpmv.f"
		}
#line 269 "ctpmv.f"
	    } else {
#line 270 "ctpmv.f"
		kx += (*n - 1) * *incx;
#line 271 "ctpmv.f"
		jx = kx;
#line 272 "ctpmv.f"
		for (j = *n; j >= 1; --j) {
#line 273 "ctpmv.f"
		    i__1 = jx;
#line 273 "ctpmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 274 "ctpmv.f"
			i__1 = jx;
#line 274 "ctpmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 275 "ctpmv.f"
			ix = kx;
#line 276 "ctpmv.f"
			i__1 = kk - (*n - (j + 1));
#line 276 "ctpmv.f"
			for (k = kk; k >= i__1; --k) {
#line 277 "ctpmv.f"
			    i__2 = ix;
#line 277 "ctpmv.f"
			    i__3 = ix;
#line 277 "ctpmv.f"
			    i__4 = k;
#line 277 "ctpmv.f"
			    z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4]
				    .i, z__2.i = temp.r * ap[i__4].i + temp.i 
				    * ap[i__4].r;
#line 277 "ctpmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 277 "ctpmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 278 "ctpmv.f"
			    ix -= *incx;
#line 279 "ctpmv.f"
/* L70: */
#line 279 "ctpmv.f"
			}
#line 280 "ctpmv.f"
			if (nounit) {
#line 280 "ctpmv.f"
			    i__1 = jx;
#line 280 "ctpmv.f"
			    i__2 = jx;
#line 280 "ctpmv.f"
			    i__3 = kk - *n + j;
#line 280 "ctpmv.f"
			    z__1.r = x[i__2].r * ap[i__3].r - x[i__2].i * ap[
				    i__3].i, z__1.i = x[i__2].r * ap[i__3].i 
				    + x[i__2].i * ap[i__3].r;
#line 280 "ctpmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 280 "ctpmv.f"
			}
#line 281 "ctpmv.f"
		    }
#line 282 "ctpmv.f"
		    jx -= *incx;
#line 283 "ctpmv.f"
		    kk -= *n - j + 1;
#line 284 "ctpmv.f"
/* L80: */
#line 284 "ctpmv.f"
		}
#line 285 "ctpmv.f"
	    }
#line 286 "ctpmv.f"
	}
#line 287 "ctpmv.f"
    } else {

/*        Form  x := A**T*x  or  x := A**H*x. */

#line 291 "ctpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 292 "ctpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 293 "ctpmv.f"
	    if (*incx == 1) {
#line 294 "ctpmv.f"
		for (j = *n; j >= 1; --j) {
#line 295 "ctpmv.f"
		    i__1 = j;
#line 295 "ctpmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 296 "ctpmv.f"
		    k = kk - 1;
#line 297 "ctpmv.f"
		    if (noconj) {
#line 298 "ctpmv.f"
			if (nounit) {
#line 298 "ctpmv.f"
			    i__1 = kk;
#line 298 "ctpmv.f"
			    z__1.r = temp.r * ap[i__1].r - temp.i * ap[i__1]
				    .i, z__1.i = temp.r * ap[i__1].i + temp.i 
				    * ap[i__1].r;
#line 298 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 298 "ctpmv.f"
			}
#line 299 "ctpmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 300 "ctpmv.f"
			    i__1 = k;
#line 300 "ctpmv.f"
			    i__2 = i__;
#line 300 "ctpmv.f"
			    z__2.r = ap[i__1].r * x[i__2].r - ap[i__1].i * x[
				    i__2].i, z__2.i = ap[i__1].r * x[i__2].i 
				    + ap[i__1].i * x[i__2].r;
#line 300 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 300 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ctpmv.f"
			    --k;
#line 302 "ctpmv.f"
/* L90: */
#line 302 "ctpmv.f"
			}
#line 303 "ctpmv.f"
		    } else {
#line 304 "ctpmv.f"
			if (nounit) {
#line 304 "ctpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 304 "ctpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 304 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 304 "ctpmv.f"
			}
#line 305 "ctpmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 306 "ctpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 306 "ctpmv.f"
			    i__1 = i__;
#line 306 "ctpmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 306 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 306 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 307 "ctpmv.f"
			    --k;
#line 308 "ctpmv.f"
/* L100: */
#line 308 "ctpmv.f"
			}
#line 309 "ctpmv.f"
		    }
#line 310 "ctpmv.f"
		    i__1 = j;
#line 310 "ctpmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 311 "ctpmv.f"
		    kk -= j;
#line 312 "ctpmv.f"
/* L110: */
#line 312 "ctpmv.f"
		}
#line 313 "ctpmv.f"
	    } else {
#line 314 "ctpmv.f"
		jx = kx + (*n - 1) * *incx;
#line 315 "ctpmv.f"
		for (j = *n; j >= 1; --j) {
#line 316 "ctpmv.f"
		    i__1 = jx;
#line 316 "ctpmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 317 "ctpmv.f"
		    ix = jx;
#line 318 "ctpmv.f"
		    if (noconj) {
#line 319 "ctpmv.f"
			if (nounit) {
#line 319 "ctpmv.f"
			    i__1 = kk;
#line 319 "ctpmv.f"
			    z__1.r = temp.r * ap[i__1].r - temp.i * ap[i__1]
				    .i, z__1.i = temp.r * ap[i__1].i + temp.i 
				    * ap[i__1].r;
#line 319 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 319 "ctpmv.f"
			}
#line 320 "ctpmv.f"
			i__1 = kk - j + 1;
#line 320 "ctpmv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 321 "ctpmv.f"
			    ix -= *incx;
#line 322 "ctpmv.f"
			    i__2 = k;
#line 322 "ctpmv.f"
			    i__3 = ix;
#line 322 "ctpmv.f"
			    z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[
				    i__3].i, z__2.i = ap[i__2].r * x[i__3].i 
				    + ap[i__2].i * x[i__3].r;
#line 322 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 322 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 323 "ctpmv.f"
/* L120: */
#line 323 "ctpmv.f"
			}
#line 324 "ctpmv.f"
		    } else {
#line 325 "ctpmv.f"
			if (nounit) {
#line 325 "ctpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 325 "ctpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 325 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 325 "ctpmv.f"
			}
#line 326 "ctpmv.f"
			i__1 = kk - j + 1;
#line 326 "ctpmv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 327 "ctpmv.f"
			    ix -= *incx;
#line 328 "ctpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 328 "ctpmv.f"
			    i__2 = ix;
#line 328 "ctpmv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 328 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 328 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 329 "ctpmv.f"
/* L130: */
#line 329 "ctpmv.f"
			}
#line 330 "ctpmv.f"
		    }
#line 331 "ctpmv.f"
		    i__1 = jx;
#line 331 "ctpmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 332 "ctpmv.f"
		    jx -= *incx;
#line 333 "ctpmv.f"
		    kk -= j;
#line 334 "ctpmv.f"
/* L140: */
#line 334 "ctpmv.f"
		}
#line 335 "ctpmv.f"
	    }
#line 336 "ctpmv.f"
	} else {
#line 337 "ctpmv.f"
	    kk = 1;
#line 338 "ctpmv.f"
	    if (*incx == 1) {
#line 339 "ctpmv.f"
		i__1 = *n;
#line 339 "ctpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 340 "ctpmv.f"
		    i__2 = j;
#line 340 "ctpmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 341 "ctpmv.f"
		    k = kk + 1;
#line 342 "ctpmv.f"
		    if (noconj) {
#line 343 "ctpmv.f"
			if (nounit) {
#line 343 "ctpmv.f"
			    i__2 = kk;
#line 343 "ctpmv.f"
			    z__1.r = temp.r * ap[i__2].r - temp.i * ap[i__2]
				    .i, z__1.i = temp.r * ap[i__2].i + temp.i 
				    * ap[i__2].r;
#line 343 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 343 "ctpmv.f"
			}
#line 344 "ctpmv.f"
			i__2 = *n;
#line 344 "ctpmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 345 "ctpmv.f"
			    i__3 = k;
#line 345 "ctpmv.f"
			    i__4 = i__;
#line 345 "ctpmv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 345 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 345 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 346 "ctpmv.f"
			    ++k;
#line 347 "ctpmv.f"
/* L150: */
#line 347 "ctpmv.f"
			}
#line 348 "ctpmv.f"
		    } else {
#line 349 "ctpmv.f"
			if (nounit) {
#line 349 "ctpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 349 "ctpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 349 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 349 "ctpmv.f"
			}
#line 350 "ctpmv.f"
			i__2 = *n;
#line 350 "ctpmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 351 "ctpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 351 "ctpmv.f"
			    i__3 = i__;
#line 351 "ctpmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 351 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 351 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 352 "ctpmv.f"
			    ++k;
#line 353 "ctpmv.f"
/* L160: */
#line 353 "ctpmv.f"
			}
#line 354 "ctpmv.f"
		    }
#line 355 "ctpmv.f"
		    i__2 = j;
#line 355 "ctpmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 356 "ctpmv.f"
		    kk += *n - j + 1;
#line 357 "ctpmv.f"
/* L170: */
#line 357 "ctpmv.f"
		}
#line 358 "ctpmv.f"
	    } else {
#line 359 "ctpmv.f"
		jx = kx;
#line 360 "ctpmv.f"
		i__1 = *n;
#line 360 "ctpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 361 "ctpmv.f"
		    i__2 = jx;
#line 361 "ctpmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 362 "ctpmv.f"
		    ix = jx;
#line 363 "ctpmv.f"
		    if (noconj) {
#line 364 "ctpmv.f"
			if (nounit) {
#line 364 "ctpmv.f"
			    i__2 = kk;
#line 364 "ctpmv.f"
			    z__1.r = temp.r * ap[i__2].r - temp.i * ap[i__2]
				    .i, z__1.i = temp.r * ap[i__2].i + temp.i 
				    * ap[i__2].r;
#line 364 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 364 "ctpmv.f"
			}
#line 365 "ctpmv.f"
			i__2 = kk + *n - j;
#line 365 "ctpmv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 366 "ctpmv.f"
			    ix += *incx;
#line 367 "ctpmv.f"
			    i__3 = k;
#line 367 "ctpmv.f"
			    i__4 = ix;
#line 367 "ctpmv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 367 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 367 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 368 "ctpmv.f"
/* L180: */
#line 368 "ctpmv.f"
			}
#line 369 "ctpmv.f"
		    } else {
#line 370 "ctpmv.f"
			if (nounit) {
#line 370 "ctpmv.f"
			    d_cnjg(&z__2, &ap[kk]);
#line 370 "ctpmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 370 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ctpmv.f"
			}
#line 371 "ctpmv.f"
			i__2 = kk + *n - j;
#line 371 "ctpmv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 372 "ctpmv.f"
			    ix += *incx;
#line 373 "ctpmv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 373 "ctpmv.f"
			    i__3 = ix;
#line 373 "ctpmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 373 "ctpmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 373 "ctpmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 374 "ctpmv.f"
/* L190: */
#line 374 "ctpmv.f"
			}
#line 375 "ctpmv.f"
		    }
#line 376 "ctpmv.f"
		    i__2 = jx;
#line 376 "ctpmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 377 "ctpmv.f"
		    jx += *incx;
#line 378 "ctpmv.f"
		    kk += *n - j + 1;
#line 379 "ctpmv.f"
/* L200: */
#line 379 "ctpmv.f"
		}
#line 380 "ctpmv.f"
	    }
#line 381 "ctpmv.f"
	}
#line 382 "ctpmv.f"
    }

#line 384 "ctpmv.f"
    return 0;

/*     End of CTPMV . */

} /* ctpmv_ */

