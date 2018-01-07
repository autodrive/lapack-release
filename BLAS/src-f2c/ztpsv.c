#line 1 "ztpsv.f"
/* ztpsv.f -- translated by f2c (version 20100827).
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

#line 1 "ztpsv.f"
/* > \brief \b ZTPSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

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
/* > ZTPSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b,   or   A**H*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular matrix, supplied in packed form. */
/* > */
/* > No test for singularity or near-singularity is included in this */
/* > routine. Such tests must be performed before calling this routine. */
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
/* >           On entry, TRANS specifies the equations to be solved as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   A*x = b. */
/* > */
/* >              TRANS = 'T' or 't'   A**T*x = b. */
/* > */
/* >              TRANS = 'C' or 'c'   A**H*x = b. */
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
/* >           element right-hand side vector b. On exit, X is overwritten */
/* >           with the solution vector x. */
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
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztpsv_(char *uplo, char *trans, char *diag, integer *n, 
	doublecomplex *ap, doublecomplex *x, integer *incx, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

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

#line 184 "ztpsv.f"
    /* Parameter adjustments */
#line 184 "ztpsv.f"
    --x;
#line 184 "ztpsv.f"
    --ap;
#line 184 "ztpsv.f"

#line 184 "ztpsv.f"
    /* Function Body */
#line 184 "ztpsv.f"
    info = 0;
#line 185 "ztpsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 186 "ztpsv.f"
	info = 1;
#line 187 "ztpsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 189 "ztpsv.f"
	info = 2;
#line 190 "ztpsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 191 "ztpsv.f"
	info = 3;
#line 192 "ztpsv.f"
    } else if (*n < 0) {
#line 193 "ztpsv.f"
	info = 4;
#line 194 "ztpsv.f"
    } else if (*incx == 0) {
#line 195 "ztpsv.f"
	info = 7;
#line 196 "ztpsv.f"
    }
#line 197 "ztpsv.f"
    if (info != 0) {
#line 198 "ztpsv.f"
	xerbla_("ZTPSV ", &info, (ftnlen)6);
#line 199 "ztpsv.f"
	return 0;
#line 200 "ztpsv.f"
    }

/*     Quick return if possible. */

#line 204 "ztpsv.f"
    if (*n == 0) {
#line 204 "ztpsv.f"
	return 0;
#line 204 "ztpsv.f"
    }

#line 206 "ztpsv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 207 "ztpsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 212 "ztpsv.f"
    if (*incx <= 0) {
#line 213 "ztpsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 214 "ztpsv.f"
    } else if (*incx != 1) {
#line 215 "ztpsv.f"
	kx = 1;
#line 216 "ztpsv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 221 "ztpsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 225 "ztpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 226 "ztpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 227 "ztpsv.f"
	    if (*incx == 1) {
#line 228 "ztpsv.f"
		for (j = *n; j >= 1; --j) {
#line 229 "ztpsv.f"
		    i__1 = j;
#line 229 "ztpsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 230 "ztpsv.f"
			if (nounit) {
#line 230 "ztpsv.f"
			    i__1 = j;
#line 230 "ztpsv.f"
			    z_div(&z__1, &x[j], &ap[kk]);
#line 230 "ztpsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 230 "ztpsv.f"
			}
#line 231 "ztpsv.f"
			i__1 = j;
#line 231 "ztpsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 232 "ztpsv.f"
			k = kk - 1;
#line 233 "ztpsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 234 "ztpsv.f"
			    i__1 = i__;
#line 234 "ztpsv.f"
			    i__2 = i__;
#line 234 "ztpsv.f"
			    i__3 = k;
#line 234 "ztpsv.f"
			    z__2.r = temp.r * ap[i__3].r - temp.i * ap[i__3]
				    .i, z__2.i = temp.r * ap[i__3].i + temp.i 
				    * ap[i__3].r;
#line 234 "ztpsv.f"
			    z__1.r = x[i__2].r - z__2.r, z__1.i = x[i__2].i - 
				    z__2.i;
#line 234 "ztpsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 235 "ztpsv.f"
			    --k;
#line 236 "ztpsv.f"
/* L10: */
#line 236 "ztpsv.f"
			}
#line 237 "ztpsv.f"
		    }
#line 238 "ztpsv.f"
		    kk -= j;
#line 239 "ztpsv.f"
/* L20: */
#line 239 "ztpsv.f"
		}
#line 240 "ztpsv.f"
	    } else {
#line 241 "ztpsv.f"
		jx = kx + (*n - 1) * *incx;
#line 242 "ztpsv.f"
		for (j = *n; j >= 1; --j) {
#line 243 "ztpsv.f"
		    i__1 = jx;
#line 243 "ztpsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 244 "ztpsv.f"
			if (nounit) {
#line 244 "ztpsv.f"
			    i__1 = jx;
#line 244 "ztpsv.f"
			    z_div(&z__1, &x[jx], &ap[kk]);
#line 244 "ztpsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 244 "ztpsv.f"
			}
#line 245 "ztpsv.f"
			i__1 = jx;
#line 245 "ztpsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 246 "ztpsv.f"
			ix = jx;
#line 247 "ztpsv.f"
			i__1 = kk - j + 1;
#line 247 "ztpsv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 248 "ztpsv.f"
			    ix -= *incx;
#line 249 "ztpsv.f"
			    i__2 = ix;
#line 249 "ztpsv.f"
			    i__3 = ix;
#line 249 "ztpsv.f"
			    i__4 = k;
#line 249 "ztpsv.f"
			    z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4]
				    .i, z__2.i = temp.r * ap[i__4].i + temp.i 
				    * ap[i__4].r;
#line 249 "ztpsv.f"
			    z__1.r = x[i__3].r - z__2.r, z__1.i = x[i__3].i - 
				    z__2.i;
#line 249 "ztpsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 250 "ztpsv.f"
/* L30: */
#line 250 "ztpsv.f"
			}
#line 251 "ztpsv.f"
		    }
#line 252 "ztpsv.f"
		    jx -= *incx;
#line 253 "ztpsv.f"
		    kk -= j;
#line 254 "ztpsv.f"
/* L40: */
#line 254 "ztpsv.f"
		}
#line 255 "ztpsv.f"
	    }
#line 256 "ztpsv.f"
	} else {
#line 257 "ztpsv.f"
	    kk = 1;
#line 258 "ztpsv.f"
	    if (*incx == 1) {
#line 259 "ztpsv.f"
		i__1 = *n;
#line 259 "ztpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "ztpsv.f"
		    i__2 = j;
#line 260 "ztpsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 261 "ztpsv.f"
			if (nounit) {
#line 261 "ztpsv.f"
			    i__2 = j;
#line 261 "ztpsv.f"
			    z_div(&z__1, &x[j], &ap[kk]);
#line 261 "ztpsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 261 "ztpsv.f"
			}
#line 262 "ztpsv.f"
			i__2 = j;
#line 262 "ztpsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 263 "ztpsv.f"
			k = kk + 1;
#line 264 "ztpsv.f"
			i__2 = *n;
#line 264 "ztpsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 265 "ztpsv.f"
			    i__3 = i__;
#line 265 "ztpsv.f"
			    i__4 = i__;
#line 265 "ztpsv.f"
			    i__5 = k;
#line 265 "ztpsv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 265 "ztpsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 265 "ztpsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 266 "ztpsv.f"
			    ++k;
#line 267 "ztpsv.f"
/* L50: */
#line 267 "ztpsv.f"
			}
#line 268 "ztpsv.f"
		    }
#line 269 "ztpsv.f"
		    kk += *n - j + 1;
#line 270 "ztpsv.f"
/* L60: */
#line 270 "ztpsv.f"
		}
#line 271 "ztpsv.f"
	    } else {
#line 272 "ztpsv.f"
		jx = kx;
#line 273 "ztpsv.f"
		i__1 = *n;
#line 273 "ztpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "ztpsv.f"
		    i__2 = jx;
#line 274 "ztpsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 275 "ztpsv.f"
			if (nounit) {
#line 275 "ztpsv.f"
			    i__2 = jx;
#line 275 "ztpsv.f"
			    z_div(&z__1, &x[jx], &ap[kk]);
#line 275 "ztpsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 275 "ztpsv.f"
			}
#line 276 "ztpsv.f"
			i__2 = jx;
#line 276 "ztpsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 277 "ztpsv.f"
			ix = jx;
#line 278 "ztpsv.f"
			i__2 = kk + *n - j;
#line 278 "ztpsv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 279 "ztpsv.f"
			    ix += *incx;
#line 280 "ztpsv.f"
			    i__3 = ix;
#line 280 "ztpsv.f"
			    i__4 = ix;
#line 280 "ztpsv.f"
			    i__5 = k;
#line 280 "ztpsv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 280 "ztpsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 280 "ztpsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 281 "ztpsv.f"
/* L70: */
#line 281 "ztpsv.f"
			}
#line 282 "ztpsv.f"
		    }
#line 283 "ztpsv.f"
		    jx += *incx;
#line 284 "ztpsv.f"
		    kk += *n - j + 1;
#line 285 "ztpsv.f"
/* L80: */
#line 285 "ztpsv.f"
		}
#line 286 "ztpsv.f"
	    }
#line 287 "ztpsv.f"
	}
#line 288 "ztpsv.f"
    } else {

/*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x. */

#line 292 "ztpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 293 "ztpsv.f"
	    kk = 1;
#line 294 "ztpsv.f"
	    if (*incx == 1) {
#line 295 "ztpsv.f"
		i__1 = *n;
#line 295 "ztpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "ztpsv.f"
		    i__2 = j;
#line 296 "ztpsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 297 "ztpsv.f"
		    k = kk;
#line 298 "ztpsv.f"
		    if (noconj) {
#line 299 "ztpsv.f"
			i__2 = j - 1;
#line 299 "ztpsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 300 "ztpsv.f"
			    i__3 = k;
#line 300 "ztpsv.f"
			    i__4 = i__;
#line 300 "ztpsv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 300 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 300 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ztpsv.f"
			    ++k;
#line 302 "ztpsv.f"
/* L90: */
#line 302 "ztpsv.f"
			}
#line 303 "ztpsv.f"
			if (nounit) {
#line 303 "ztpsv.f"
			    z_div(&z__1, &temp, &ap[kk + j - 1]);
#line 303 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 303 "ztpsv.f"
			}
#line 304 "ztpsv.f"
		    } else {
#line 305 "ztpsv.f"
			i__2 = j - 1;
#line 305 "ztpsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 306 "ztpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 306 "ztpsv.f"
			    i__3 = i__;
#line 306 "ztpsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 306 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 306 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 307 "ztpsv.f"
			    ++k;
#line 308 "ztpsv.f"
/* L100: */
#line 308 "ztpsv.f"
			}
#line 309 "ztpsv.f"
			if (nounit) {
#line 309 "ztpsv.f"
			    d_cnjg(&z__2, &ap[kk + j - 1]);
#line 309 "ztpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 309 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 309 "ztpsv.f"
			}
#line 310 "ztpsv.f"
		    }
#line 311 "ztpsv.f"
		    i__2 = j;
#line 311 "ztpsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 312 "ztpsv.f"
		    kk += j;
#line 313 "ztpsv.f"
/* L110: */
#line 313 "ztpsv.f"
		}
#line 314 "ztpsv.f"
	    } else {
#line 315 "ztpsv.f"
		jx = kx;
#line 316 "ztpsv.f"
		i__1 = *n;
#line 316 "ztpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 317 "ztpsv.f"
		    i__2 = jx;
#line 317 "ztpsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 318 "ztpsv.f"
		    ix = kx;
#line 319 "ztpsv.f"
		    if (noconj) {
#line 320 "ztpsv.f"
			i__2 = kk + j - 2;
#line 320 "ztpsv.f"
			for (k = kk; k <= i__2; ++k) {
#line 321 "ztpsv.f"
			    i__3 = k;
#line 321 "ztpsv.f"
			    i__4 = ix;
#line 321 "ztpsv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 321 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 321 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 322 "ztpsv.f"
			    ix += *incx;
#line 323 "ztpsv.f"
/* L120: */
#line 323 "ztpsv.f"
			}
#line 324 "ztpsv.f"
			if (nounit) {
#line 324 "ztpsv.f"
			    z_div(&z__1, &temp, &ap[kk + j - 1]);
#line 324 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 324 "ztpsv.f"
			}
#line 325 "ztpsv.f"
		    } else {
#line 326 "ztpsv.f"
			i__2 = kk + j - 2;
#line 326 "ztpsv.f"
			for (k = kk; k <= i__2; ++k) {
#line 327 "ztpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 327 "ztpsv.f"
			    i__3 = ix;
#line 327 "ztpsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 327 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 327 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 328 "ztpsv.f"
			    ix += *incx;
#line 329 "ztpsv.f"
/* L130: */
#line 329 "ztpsv.f"
			}
#line 330 "ztpsv.f"
			if (nounit) {
#line 330 "ztpsv.f"
			    d_cnjg(&z__2, &ap[kk + j - 1]);
#line 330 "ztpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 330 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 330 "ztpsv.f"
			}
#line 331 "ztpsv.f"
		    }
#line 332 "ztpsv.f"
		    i__2 = jx;
#line 332 "ztpsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 333 "ztpsv.f"
		    jx += *incx;
#line 334 "ztpsv.f"
		    kk += j;
#line 335 "ztpsv.f"
/* L140: */
#line 335 "ztpsv.f"
		}
#line 336 "ztpsv.f"
	    }
#line 337 "ztpsv.f"
	} else {
#line 338 "ztpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 339 "ztpsv.f"
	    if (*incx == 1) {
#line 340 "ztpsv.f"
		for (j = *n; j >= 1; --j) {
#line 341 "ztpsv.f"
		    i__1 = j;
#line 341 "ztpsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 342 "ztpsv.f"
		    k = kk;
#line 343 "ztpsv.f"
		    if (noconj) {
#line 344 "ztpsv.f"
			i__1 = j + 1;
#line 344 "ztpsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 345 "ztpsv.f"
			    i__2 = k;
#line 345 "ztpsv.f"
			    i__3 = i__;
#line 345 "ztpsv.f"
			    z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[
				    i__3].i, z__2.i = ap[i__2].r * x[i__3].i 
				    + ap[i__2].i * x[i__3].r;
#line 345 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 345 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 346 "ztpsv.f"
			    --k;
#line 347 "ztpsv.f"
/* L150: */
#line 347 "ztpsv.f"
			}
#line 348 "ztpsv.f"
			if (nounit) {
#line 348 "ztpsv.f"
			    z_div(&z__1, &temp, &ap[kk - *n + j]);
#line 348 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 348 "ztpsv.f"
			}
#line 349 "ztpsv.f"
		    } else {
#line 350 "ztpsv.f"
			i__1 = j + 1;
#line 350 "ztpsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 351 "ztpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 351 "ztpsv.f"
			    i__2 = i__;
#line 351 "ztpsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 351 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 351 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 352 "ztpsv.f"
			    --k;
#line 353 "ztpsv.f"
/* L160: */
#line 353 "ztpsv.f"
			}
#line 354 "ztpsv.f"
			if (nounit) {
#line 354 "ztpsv.f"
			    d_cnjg(&z__2, &ap[kk - *n + j]);
#line 354 "ztpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 354 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ztpsv.f"
			}
#line 355 "ztpsv.f"
		    }
#line 356 "ztpsv.f"
		    i__1 = j;
#line 356 "ztpsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 357 "ztpsv.f"
		    kk -= *n - j + 1;
#line 358 "ztpsv.f"
/* L170: */
#line 358 "ztpsv.f"
		}
#line 359 "ztpsv.f"
	    } else {
#line 360 "ztpsv.f"
		kx += (*n - 1) * *incx;
#line 361 "ztpsv.f"
		jx = kx;
#line 362 "ztpsv.f"
		for (j = *n; j >= 1; --j) {
#line 363 "ztpsv.f"
		    i__1 = jx;
#line 363 "ztpsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 364 "ztpsv.f"
		    ix = kx;
#line 365 "ztpsv.f"
		    if (noconj) {
#line 366 "ztpsv.f"
			i__1 = kk - (*n - (j + 1));
#line 366 "ztpsv.f"
			for (k = kk; k >= i__1; --k) {
#line 367 "ztpsv.f"
			    i__2 = k;
#line 367 "ztpsv.f"
			    i__3 = ix;
#line 367 "ztpsv.f"
			    z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[
				    i__3].i, z__2.i = ap[i__2].r * x[i__3].i 
				    + ap[i__2].i * x[i__3].r;
#line 367 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 367 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 368 "ztpsv.f"
			    ix -= *incx;
#line 369 "ztpsv.f"
/* L180: */
#line 369 "ztpsv.f"
			}
#line 370 "ztpsv.f"
			if (nounit) {
#line 370 "ztpsv.f"
			    z_div(&z__1, &temp, &ap[kk - *n + j]);
#line 370 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ztpsv.f"
			}
#line 371 "ztpsv.f"
		    } else {
#line 372 "ztpsv.f"
			i__1 = kk - (*n - (j + 1));
#line 372 "ztpsv.f"
			for (k = kk; k >= i__1; --k) {
#line 373 "ztpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 373 "ztpsv.f"
			    i__2 = ix;
#line 373 "ztpsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 373 "ztpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 373 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 374 "ztpsv.f"
			    ix -= *incx;
#line 375 "ztpsv.f"
/* L190: */
#line 375 "ztpsv.f"
			}
#line 376 "ztpsv.f"
			if (nounit) {
#line 376 "ztpsv.f"
			    d_cnjg(&z__2, &ap[kk - *n + j]);
#line 376 "ztpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 376 "ztpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 376 "ztpsv.f"
			}
#line 377 "ztpsv.f"
		    }
#line 378 "ztpsv.f"
		    i__1 = jx;
#line 378 "ztpsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 379 "ztpsv.f"
		    jx -= *incx;
#line 380 "ztpsv.f"
		    kk -= *n - j + 1;
#line 381 "ztpsv.f"
/* L200: */
#line 381 "ztpsv.f"
		}
#line 382 "ztpsv.f"
	    }
#line 383 "ztpsv.f"
	}
#line 384 "ztpsv.f"
    }

#line 386 "ztpsv.f"
    return 0;

/*     End of ZTPSV . */

} /* ztpsv_ */

