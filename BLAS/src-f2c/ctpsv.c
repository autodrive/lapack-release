#line 1 "ctpsv.f"
/* ctpsv.f -- translated by f2c (version 20100827).
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

#line 1 "ctpsv.f"
/* > \brief \b CTPSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

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
/* > CTPSV  solves one of the systems of equations */
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
/* >          AP is COMPLEX array, dimension at least */
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
/* >          X is COMPLEX array, dimension at least */
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
/* Subroutine */ int ctpsv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 184 "ctpsv.f"
    /* Parameter adjustments */
#line 184 "ctpsv.f"
    --x;
#line 184 "ctpsv.f"
    --ap;
#line 184 "ctpsv.f"

#line 184 "ctpsv.f"
    /* Function Body */
#line 184 "ctpsv.f"
    info = 0;
#line 185 "ctpsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 186 "ctpsv.f"
	info = 1;
#line 187 "ctpsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 189 "ctpsv.f"
	info = 2;
#line 190 "ctpsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 191 "ctpsv.f"
	info = 3;
#line 192 "ctpsv.f"
    } else if (*n < 0) {
#line 193 "ctpsv.f"
	info = 4;
#line 194 "ctpsv.f"
    } else if (*incx == 0) {
#line 195 "ctpsv.f"
	info = 7;
#line 196 "ctpsv.f"
    }
#line 197 "ctpsv.f"
    if (info != 0) {
#line 198 "ctpsv.f"
	xerbla_("CTPSV ", &info, (ftnlen)6);
#line 199 "ctpsv.f"
	return 0;
#line 200 "ctpsv.f"
    }

/*     Quick return if possible. */

#line 204 "ctpsv.f"
    if (*n == 0) {
#line 204 "ctpsv.f"
	return 0;
#line 204 "ctpsv.f"
    }

#line 206 "ctpsv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 207 "ctpsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 212 "ctpsv.f"
    if (*incx <= 0) {
#line 213 "ctpsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 214 "ctpsv.f"
    } else if (*incx != 1) {
#line 215 "ctpsv.f"
	kx = 1;
#line 216 "ctpsv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 221 "ctpsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 225 "ctpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 226 "ctpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 227 "ctpsv.f"
	    if (*incx == 1) {
#line 228 "ctpsv.f"
		for (j = *n; j >= 1; --j) {
#line 229 "ctpsv.f"
		    i__1 = j;
#line 229 "ctpsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 230 "ctpsv.f"
			if (nounit) {
#line 230 "ctpsv.f"
			    i__1 = j;
#line 230 "ctpsv.f"
			    z_div(&z__1, &x[j], &ap[kk]);
#line 230 "ctpsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 230 "ctpsv.f"
			}
#line 231 "ctpsv.f"
			i__1 = j;
#line 231 "ctpsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 232 "ctpsv.f"
			k = kk - 1;
#line 233 "ctpsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 234 "ctpsv.f"
			    i__1 = i__;
#line 234 "ctpsv.f"
			    i__2 = i__;
#line 234 "ctpsv.f"
			    i__3 = k;
#line 234 "ctpsv.f"
			    z__2.r = temp.r * ap[i__3].r - temp.i * ap[i__3]
				    .i, z__2.i = temp.r * ap[i__3].i + temp.i 
				    * ap[i__3].r;
#line 234 "ctpsv.f"
			    z__1.r = x[i__2].r - z__2.r, z__1.i = x[i__2].i - 
				    z__2.i;
#line 234 "ctpsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 235 "ctpsv.f"
			    --k;
#line 236 "ctpsv.f"
/* L10: */
#line 236 "ctpsv.f"
			}
#line 237 "ctpsv.f"
		    }
#line 238 "ctpsv.f"
		    kk -= j;
#line 239 "ctpsv.f"
/* L20: */
#line 239 "ctpsv.f"
		}
#line 240 "ctpsv.f"
	    } else {
#line 241 "ctpsv.f"
		jx = kx + (*n - 1) * *incx;
#line 242 "ctpsv.f"
		for (j = *n; j >= 1; --j) {
#line 243 "ctpsv.f"
		    i__1 = jx;
#line 243 "ctpsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 244 "ctpsv.f"
			if (nounit) {
#line 244 "ctpsv.f"
			    i__1 = jx;
#line 244 "ctpsv.f"
			    z_div(&z__1, &x[jx], &ap[kk]);
#line 244 "ctpsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 244 "ctpsv.f"
			}
#line 245 "ctpsv.f"
			i__1 = jx;
#line 245 "ctpsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 246 "ctpsv.f"
			ix = jx;
#line 247 "ctpsv.f"
			i__1 = kk - j + 1;
#line 247 "ctpsv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 248 "ctpsv.f"
			    ix -= *incx;
#line 249 "ctpsv.f"
			    i__2 = ix;
#line 249 "ctpsv.f"
			    i__3 = ix;
#line 249 "ctpsv.f"
			    i__4 = k;
#line 249 "ctpsv.f"
			    z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4]
				    .i, z__2.i = temp.r * ap[i__4].i + temp.i 
				    * ap[i__4].r;
#line 249 "ctpsv.f"
			    z__1.r = x[i__3].r - z__2.r, z__1.i = x[i__3].i - 
				    z__2.i;
#line 249 "ctpsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 250 "ctpsv.f"
/* L30: */
#line 250 "ctpsv.f"
			}
#line 251 "ctpsv.f"
		    }
#line 252 "ctpsv.f"
		    jx -= *incx;
#line 253 "ctpsv.f"
		    kk -= j;
#line 254 "ctpsv.f"
/* L40: */
#line 254 "ctpsv.f"
		}
#line 255 "ctpsv.f"
	    }
#line 256 "ctpsv.f"
	} else {
#line 257 "ctpsv.f"
	    kk = 1;
#line 258 "ctpsv.f"
	    if (*incx == 1) {
#line 259 "ctpsv.f"
		i__1 = *n;
#line 259 "ctpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "ctpsv.f"
		    i__2 = j;
#line 260 "ctpsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 261 "ctpsv.f"
			if (nounit) {
#line 261 "ctpsv.f"
			    i__2 = j;
#line 261 "ctpsv.f"
			    z_div(&z__1, &x[j], &ap[kk]);
#line 261 "ctpsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 261 "ctpsv.f"
			}
#line 262 "ctpsv.f"
			i__2 = j;
#line 262 "ctpsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 263 "ctpsv.f"
			k = kk + 1;
#line 264 "ctpsv.f"
			i__2 = *n;
#line 264 "ctpsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 265 "ctpsv.f"
			    i__3 = i__;
#line 265 "ctpsv.f"
			    i__4 = i__;
#line 265 "ctpsv.f"
			    i__5 = k;
#line 265 "ctpsv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 265 "ctpsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 265 "ctpsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 266 "ctpsv.f"
			    ++k;
#line 267 "ctpsv.f"
/* L50: */
#line 267 "ctpsv.f"
			}
#line 268 "ctpsv.f"
		    }
#line 269 "ctpsv.f"
		    kk += *n - j + 1;
#line 270 "ctpsv.f"
/* L60: */
#line 270 "ctpsv.f"
		}
#line 271 "ctpsv.f"
	    } else {
#line 272 "ctpsv.f"
		jx = kx;
#line 273 "ctpsv.f"
		i__1 = *n;
#line 273 "ctpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "ctpsv.f"
		    i__2 = jx;
#line 274 "ctpsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 275 "ctpsv.f"
			if (nounit) {
#line 275 "ctpsv.f"
			    i__2 = jx;
#line 275 "ctpsv.f"
			    z_div(&z__1, &x[jx], &ap[kk]);
#line 275 "ctpsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 275 "ctpsv.f"
			}
#line 276 "ctpsv.f"
			i__2 = jx;
#line 276 "ctpsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 277 "ctpsv.f"
			ix = jx;
#line 278 "ctpsv.f"
			i__2 = kk + *n - j;
#line 278 "ctpsv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 279 "ctpsv.f"
			    ix += *incx;
#line 280 "ctpsv.f"
			    i__3 = ix;
#line 280 "ctpsv.f"
			    i__4 = ix;
#line 280 "ctpsv.f"
			    i__5 = k;
#line 280 "ctpsv.f"
			    z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5]
				    .i, z__2.i = temp.r * ap[i__5].i + temp.i 
				    * ap[i__5].r;
#line 280 "ctpsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 280 "ctpsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 281 "ctpsv.f"
/* L70: */
#line 281 "ctpsv.f"
			}
#line 282 "ctpsv.f"
		    }
#line 283 "ctpsv.f"
		    jx += *incx;
#line 284 "ctpsv.f"
		    kk += *n - j + 1;
#line 285 "ctpsv.f"
/* L80: */
#line 285 "ctpsv.f"
		}
#line 286 "ctpsv.f"
	    }
#line 287 "ctpsv.f"
	}
#line 288 "ctpsv.f"
    } else {

/*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x. */

#line 292 "ctpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 293 "ctpsv.f"
	    kk = 1;
#line 294 "ctpsv.f"
	    if (*incx == 1) {
#line 295 "ctpsv.f"
		i__1 = *n;
#line 295 "ctpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "ctpsv.f"
		    i__2 = j;
#line 296 "ctpsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 297 "ctpsv.f"
		    k = kk;
#line 298 "ctpsv.f"
		    if (noconj) {
#line 299 "ctpsv.f"
			i__2 = j - 1;
#line 299 "ctpsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 300 "ctpsv.f"
			    i__3 = k;
#line 300 "ctpsv.f"
			    i__4 = i__;
#line 300 "ctpsv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 300 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 300 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ctpsv.f"
			    ++k;
#line 302 "ctpsv.f"
/* L90: */
#line 302 "ctpsv.f"
			}
#line 303 "ctpsv.f"
			if (nounit) {
#line 303 "ctpsv.f"
			    z_div(&z__1, &temp, &ap[kk + j - 1]);
#line 303 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 303 "ctpsv.f"
			}
#line 304 "ctpsv.f"
		    } else {
#line 305 "ctpsv.f"
			i__2 = j - 1;
#line 305 "ctpsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 306 "ctpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 306 "ctpsv.f"
			    i__3 = i__;
#line 306 "ctpsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 306 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 306 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 307 "ctpsv.f"
			    ++k;
#line 308 "ctpsv.f"
/* L100: */
#line 308 "ctpsv.f"
			}
#line 309 "ctpsv.f"
			if (nounit) {
#line 309 "ctpsv.f"
			    d_cnjg(&z__2, &ap[kk + j - 1]);
#line 309 "ctpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 309 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 309 "ctpsv.f"
			}
#line 310 "ctpsv.f"
		    }
#line 311 "ctpsv.f"
		    i__2 = j;
#line 311 "ctpsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 312 "ctpsv.f"
		    kk += j;
#line 313 "ctpsv.f"
/* L110: */
#line 313 "ctpsv.f"
		}
#line 314 "ctpsv.f"
	    } else {
#line 315 "ctpsv.f"
		jx = kx;
#line 316 "ctpsv.f"
		i__1 = *n;
#line 316 "ctpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 317 "ctpsv.f"
		    i__2 = jx;
#line 317 "ctpsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 318 "ctpsv.f"
		    ix = kx;
#line 319 "ctpsv.f"
		    if (noconj) {
#line 320 "ctpsv.f"
			i__2 = kk + j - 2;
#line 320 "ctpsv.f"
			for (k = kk; k <= i__2; ++k) {
#line 321 "ctpsv.f"
			    i__3 = k;
#line 321 "ctpsv.f"
			    i__4 = ix;
#line 321 "ctpsv.f"
			    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[
				    i__4].i, z__2.i = ap[i__3].r * x[i__4].i 
				    + ap[i__3].i * x[i__4].r;
#line 321 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 321 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 322 "ctpsv.f"
			    ix += *incx;
#line 323 "ctpsv.f"
/* L120: */
#line 323 "ctpsv.f"
			}
#line 324 "ctpsv.f"
			if (nounit) {
#line 324 "ctpsv.f"
			    z_div(&z__1, &temp, &ap[kk + j - 1]);
#line 324 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 324 "ctpsv.f"
			}
#line 325 "ctpsv.f"
		    } else {
#line 326 "ctpsv.f"
			i__2 = kk + j - 2;
#line 326 "ctpsv.f"
			for (k = kk; k <= i__2; ++k) {
#line 327 "ctpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 327 "ctpsv.f"
			    i__3 = ix;
#line 327 "ctpsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 327 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 327 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 328 "ctpsv.f"
			    ix += *incx;
#line 329 "ctpsv.f"
/* L130: */
#line 329 "ctpsv.f"
			}
#line 330 "ctpsv.f"
			if (nounit) {
#line 330 "ctpsv.f"
			    d_cnjg(&z__2, &ap[kk + j - 1]);
#line 330 "ctpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 330 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 330 "ctpsv.f"
			}
#line 331 "ctpsv.f"
		    }
#line 332 "ctpsv.f"
		    i__2 = jx;
#line 332 "ctpsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 333 "ctpsv.f"
		    jx += *incx;
#line 334 "ctpsv.f"
		    kk += j;
#line 335 "ctpsv.f"
/* L140: */
#line 335 "ctpsv.f"
		}
#line 336 "ctpsv.f"
	    }
#line 337 "ctpsv.f"
	} else {
#line 338 "ctpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 339 "ctpsv.f"
	    if (*incx == 1) {
#line 340 "ctpsv.f"
		for (j = *n; j >= 1; --j) {
#line 341 "ctpsv.f"
		    i__1 = j;
#line 341 "ctpsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 342 "ctpsv.f"
		    k = kk;
#line 343 "ctpsv.f"
		    if (noconj) {
#line 344 "ctpsv.f"
			i__1 = j + 1;
#line 344 "ctpsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 345 "ctpsv.f"
			    i__2 = k;
#line 345 "ctpsv.f"
			    i__3 = i__;
#line 345 "ctpsv.f"
			    z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[
				    i__3].i, z__2.i = ap[i__2].r * x[i__3].i 
				    + ap[i__2].i * x[i__3].r;
#line 345 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 345 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 346 "ctpsv.f"
			    --k;
#line 347 "ctpsv.f"
/* L150: */
#line 347 "ctpsv.f"
			}
#line 348 "ctpsv.f"
			if (nounit) {
#line 348 "ctpsv.f"
			    z_div(&z__1, &temp, &ap[kk - *n + j]);
#line 348 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 348 "ctpsv.f"
			}
#line 349 "ctpsv.f"
		    } else {
#line 350 "ctpsv.f"
			i__1 = j + 1;
#line 350 "ctpsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 351 "ctpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 351 "ctpsv.f"
			    i__2 = i__;
#line 351 "ctpsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 351 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 351 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 352 "ctpsv.f"
			    --k;
#line 353 "ctpsv.f"
/* L160: */
#line 353 "ctpsv.f"
			}
#line 354 "ctpsv.f"
			if (nounit) {
#line 354 "ctpsv.f"
			    d_cnjg(&z__2, &ap[kk - *n + j]);
#line 354 "ctpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 354 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ctpsv.f"
			}
#line 355 "ctpsv.f"
		    }
#line 356 "ctpsv.f"
		    i__1 = j;
#line 356 "ctpsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 357 "ctpsv.f"
		    kk -= *n - j + 1;
#line 358 "ctpsv.f"
/* L170: */
#line 358 "ctpsv.f"
		}
#line 359 "ctpsv.f"
	    } else {
#line 360 "ctpsv.f"
		kx += (*n - 1) * *incx;
#line 361 "ctpsv.f"
		jx = kx;
#line 362 "ctpsv.f"
		for (j = *n; j >= 1; --j) {
#line 363 "ctpsv.f"
		    i__1 = jx;
#line 363 "ctpsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 364 "ctpsv.f"
		    ix = kx;
#line 365 "ctpsv.f"
		    if (noconj) {
#line 366 "ctpsv.f"
			i__1 = kk - (*n - (j + 1));
#line 366 "ctpsv.f"
			for (k = kk; k >= i__1; --k) {
#line 367 "ctpsv.f"
			    i__2 = k;
#line 367 "ctpsv.f"
			    i__3 = ix;
#line 367 "ctpsv.f"
			    z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[
				    i__3].i, z__2.i = ap[i__2].r * x[i__3].i 
				    + ap[i__2].i * x[i__3].r;
#line 367 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 367 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 368 "ctpsv.f"
			    ix -= *incx;
#line 369 "ctpsv.f"
/* L180: */
#line 369 "ctpsv.f"
			}
#line 370 "ctpsv.f"
			if (nounit) {
#line 370 "ctpsv.f"
			    z_div(&z__1, &temp, &ap[kk - *n + j]);
#line 370 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ctpsv.f"
			}
#line 371 "ctpsv.f"
		    } else {
#line 372 "ctpsv.f"
			i__1 = kk - (*n - (j + 1));
#line 372 "ctpsv.f"
			for (k = kk; k >= i__1; --k) {
#line 373 "ctpsv.f"
			    d_cnjg(&z__3, &ap[k]);
#line 373 "ctpsv.f"
			    i__2 = ix;
#line 373 "ctpsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 373 "ctpsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 373 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 374 "ctpsv.f"
			    ix -= *incx;
#line 375 "ctpsv.f"
/* L190: */
#line 375 "ctpsv.f"
			}
#line 376 "ctpsv.f"
			if (nounit) {
#line 376 "ctpsv.f"
			    d_cnjg(&z__2, &ap[kk - *n + j]);
#line 376 "ctpsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 376 "ctpsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 376 "ctpsv.f"
			}
#line 377 "ctpsv.f"
		    }
#line 378 "ctpsv.f"
		    i__1 = jx;
#line 378 "ctpsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 379 "ctpsv.f"
		    jx -= *incx;
#line 380 "ctpsv.f"
		    kk -= *n - j + 1;
#line 381 "ctpsv.f"
/* L200: */
#line 381 "ctpsv.f"
		}
#line 382 "ctpsv.f"
	    }
#line 383 "ctpsv.f"
	}
#line 384 "ctpsv.f"
    }

#line 386 "ctpsv.f"
    return 0;

/*     End of CTPSV . */

} /* ctpsv_ */

