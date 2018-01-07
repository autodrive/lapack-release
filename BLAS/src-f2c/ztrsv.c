#line 1 "ztrsv.f"
/* ztrsv.f -- translated by f2c (version 20100827).
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

#line 1 "ztrsv.f"
/* > \brief \b ZTRSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b,   or   A**H*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular matrix. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular matrix and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular matrix and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced either, but are assumed to be unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of dimension at least */
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
/* Subroutine */ int ztrsv_(char *uplo, char *trans, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
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

#line 189 "ztrsv.f"
    /* Parameter adjustments */
#line 189 "ztrsv.f"
    a_dim1 = *lda;
#line 189 "ztrsv.f"
    a_offset = 1 + a_dim1;
#line 189 "ztrsv.f"
    a -= a_offset;
#line 189 "ztrsv.f"
    --x;
#line 189 "ztrsv.f"

#line 189 "ztrsv.f"
    /* Function Body */
#line 189 "ztrsv.f"
    info = 0;
#line 190 "ztrsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 191 "ztrsv.f"
	info = 1;
#line 192 "ztrsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 194 "ztrsv.f"
	info = 2;
#line 195 "ztrsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 196 "ztrsv.f"
	info = 3;
#line 197 "ztrsv.f"
    } else if (*n < 0) {
#line 198 "ztrsv.f"
	info = 4;
#line 199 "ztrsv.f"
    } else if (*lda < max(1,*n)) {
#line 200 "ztrsv.f"
	info = 6;
#line 201 "ztrsv.f"
    } else if (*incx == 0) {
#line 202 "ztrsv.f"
	info = 8;
#line 203 "ztrsv.f"
    }
#line 204 "ztrsv.f"
    if (info != 0) {
#line 205 "ztrsv.f"
	xerbla_("ZTRSV ", &info, (ftnlen)6);
#line 206 "ztrsv.f"
	return 0;
#line 207 "ztrsv.f"
    }

/*     Quick return if possible. */

#line 211 "ztrsv.f"
    if (*n == 0) {
#line 211 "ztrsv.f"
	return 0;
#line 211 "ztrsv.f"
    }

#line 213 "ztrsv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 214 "ztrsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 219 "ztrsv.f"
    if (*incx <= 0) {
#line 220 "ztrsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 221 "ztrsv.f"
    } else if (*incx != 1) {
#line 222 "ztrsv.f"
	kx = 1;
#line 223 "ztrsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 228 "ztrsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 232 "ztrsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 233 "ztrsv.f"
	    if (*incx == 1) {
#line 234 "ztrsv.f"
		for (j = *n; j >= 1; --j) {
#line 235 "ztrsv.f"
		    i__1 = j;
#line 235 "ztrsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 236 "ztrsv.f"
			if (nounit) {
#line 236 "ztrsv.f"
			    i__1 = j;
#line 236 "ztrsv.f"
			    z_div(&z__1, &x[j], &a[j + j * a_dim1]);
#line 236 "ztrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 236 "ztrsv.f"
			}
#line 237 "ztrsv.f"
			i__1 = j;
#line 237 "ztrsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 238 "ztrsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 239 "ztrsv.f"
			    i__1 = i__;
#line 239 "ztrsv.f"
			    i__2 = i__;
#line 239 "ztrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 239 "ztrsv.f"
			    z__2.r = temp.r * a[i__3].r - temp.i * a[i__3].i, 
				    z__2.i = temp.r * a[i__3].i + temp.i * a[
				    i__3].r;
#line 239 "ztrsv.f"
			    z__1.r = x[i__2].r - z__2.r, z__1.i = x[i__2].i - 
				    z__2.i;
#line 239 "ztrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 240 "ztrsv.f"
/* L10: */
#line 240 "ztrsv.f"
			}
#line 241 "ztrsv.f"
		    }
#line 242 "ztrsv.f"
/* L20: */
#line 242 "ztrsv.f"
		}
#line 243 "ztrsv.f"
	    } else {
#line 244 "ztrsv.f"
		jx = kx + (*n - 1) * *incx;
#line 245 "ztrsv.f"
		for (j = *n; j >= 1; --j) {
#line 246 "ztrsv.f"
		    i__1 = jx;
#line 246 "ztrsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 247 "ztrsv.f"
			if (nounit) {
#line 247 "ztrsv.f"
			    i__1 = jx;
#line 247 "ztrsv.f"
			    z_div(&z__1, &x[jx], &a[j + j * a_dim1]);
#line 247 "ztrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 247 "ztrsv.f"
			}
#line 248 "ztrsv.f"
			i__1 = jx;
#line 248 "ztrsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 249 "ztrsv.f"
			ix = jx;
#line 250 "ztrsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 251 "ztrsv.f"
			    ix -= *incx;
#line 252 "ztrsv.f"
			    i__1 = ix;
#line 252 "ztrsv.f"
			    i__2 = ix;
#line 252 "ztrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 252 "ztrsv.f"
			    z__2.r = temp.r * a[i__3].r - temp.i * a[i__3].i, 
				    z__2.i = temp.r * a[i__3].i + temp.i * a[
				    i__3].r;
#line 252 "ztrsv.f"
			    z__1.r = x[i__2].r - z__2.r, z__1.i = x[i__2].i - 
				    z__2.i;
#line 252 "ztrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 253 "ztrsv.f"
/* L30: */
#line 253 "ztrsv.f"
			}
#line 254 "ztrsv.f"
		    }
#line 255 "ztrsv.f"
		    jx -= *incx;
#line 256 "ztrsv.f"
/* L40: */
#line 256 "ztrsv.f"
		}
#line 257 "ztrsv.f"
	    }
#line 258 "ztrsv.f"
	} else {
#line 259 "ztrsv.f"
	    if (*incx == 1) {
#line 260 "ztrsv.f"
		i__1 = *n;
#line 260 "ztrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 261 "ztrsv.f"
		    i__2 = j;
#line 261 "ztrsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 262 "ztrsv.f"
			if (nounit) {
#line 262 "ztrsv.f"
			    i__2 = j;
#line 262 "ztrsv.f"
			    z_div(&z__1, &x[j], &a[j + j * a_dim1]);
#line 262 "ztrsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 262 "ztrsv.f"
			}
#line 263 "ztrsv.f"
			i__2 = j;
#line 263 "ztrsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 264 "ztrsv.f"
			i__2 = *n;
#line 264 "ztrsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 265 "ztrsv.f"
			    i__3 = i__;
#line 265 "ztrsv.f"
			    i__4 = i__;
#line 265 "ztrsv.f"
			    i__5 = i__ + j * a_dim1;
#line 265 "ztrsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 265 "ztrsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 265 "ztrsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 266 "ztrsv.f"
/* L50: */
#line 266 "ztrsv.f"
			}
#line 267 "ztrsv.f"
		    }
#line 268 "ztrsv.f"
/* L60: */
#line 268 "ztrsv.f"
		}
#line 269 "ztrsv.f"
	    } else {
#line 270 "ztrsv.f"
		jx = kx;
#line 271 "ztrsv.f"
		i__1 = *n;
#line 271 "ztrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 272 "ztrsv.f"
		    i__2 = jx;
#line 272 "ztrsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 273 "ztrsv.f"
			if (nounit) {
#line 273 "ztrsv.f"
			    i__2 = jx;
#line 273 "ztrsv.f"
			    z_div(&z__1, &x[jx], &a[j + j * a_dim1]);
#line 273 "ztrsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 273 "ztrsv.f"
			}
#line 274 "ztrsv.f"
			i__2 = jx;
#line 274 "ztrsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 275 "ztrsv.f"
			ix = jx;
#line 276 "ztrsv.f"
			i__2 = *n;
#line 276 "ztrsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 277 "ztrsv.f"
			    ix += *incx;
#line 278 "ztrsv.f"
			    i__3 = ix;
#line 278 "ztrsv.f"
			    i__4 = ix;
#line 278 "ztrsv.f"
			    i__5 = i__ + j * a_dim1;
#line 278 "ztrsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 278 "ztrsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 278 "ztrsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 279 "ztrsv.f"
/* L70: */
#line 279 "ztrsv.f"
			}
#line 280 "ztrsv.f"
		    }
#line 281 "ztrsv.f"
		    jx += *incx;
#line 282 "ztrsv.f"
/* L80: */
#line 282 "ztrsv.f"
		}
#line 283 "ztrsv.f"
	    }
#line 284 "ztrsv.f"
	}
#line 285 "ztrsv.f"
    } else {

/*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x. */

#line 289 "ztrsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 290 "ztrsv.f"
	    if (*incx == 1) {
#line 291 "ztrsv.f"
		i__1 = *n;
#line 291 "ztrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "ztrsv.f"
		    i__2 = j;
#line 292 "ztrsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 293 "ztrsv.f"
		    if (noconj) {
#line 294 "ztrsv.f"
			i__2 = j - 1;
#line 294 "ztrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "ztrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 295 "ztrsv.f"
			    i__4 = i__;
#line 295 "ztrsv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 295 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 295 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 296 "ztrsv.f"
/* L90: */
#line 296 "ztrsv.f"
			}
#line 297 "ztrsv.f"
			if (nounit) {
#line 297 "ztrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 297 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 297 "ztrsv.f"
			}
#line 298 "ztrsv.f"
		    } else {
#line 299 "ztrsv.f"
			i__2 = j - 1;
#line 299 "ztrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 300 "ztrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 300 "ztrsv.f"
			    i__3 = i__;
#line 300 "ztrsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 300 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 300 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ztrsv.f"
/* L100: */
#line 301 "ztrsv.f"
			}
#line 302 "ztrsv.f"
			if (nounit) {
#line 302 "ztrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 302 "ztrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 302 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 302 "ztrsv.f"
			}
#line 303 "ztrsv.f"
		    }
#line 304 "ztrsv.f"
		    i__2 = j;
#line 304 "ztrsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 305 "ztrsv.f"
/* L110: */
#line 305 "ztrsv.f"
		}
#line 306 "ztrsv.f"
	    } else {
#line 307 "ztrsv.f"
		jx = kx;
#line 308 "ztrsv.f"
		i__1 = *n;
#line 308 "ztrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 309 "ztrsv.f"
		    ix = kx;
#line 310 "ztrsv.f"
		    i__2 = jx;
#line 310 "ztrsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 311 "ztrsv.f"
		    if (noconj) {
#line 312 "ztrsv.f"
			i__2 = j - 1;
#line 312 "ztrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 313 "ztrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 313 "ztrsv.f"
			    i__4 = ix;
#line 313 "ztrsv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 313 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 313 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 314 "ztrsv.f"
			    ix += *incx;
#line 315 "ztrsv.f"
/* L120: */
#line 315 "ztrsv.f"
			}
#line 316 "ztrsv.f"
			if (nounit) {
#line 316 "ztrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 316 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 316 "ztrsv.f"
			}
#line 317 "ztrsv.f"
		    } else {
#line 318 "ztrsv.f"
			i__2 = j - 1;
#line 318 "ztrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "ztrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 319 "ztrsv.f"
			    i__3 = ix;
#line 319 "ztrsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 319 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 319 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 320 "ztrsv.f"
			    ix += *incx;
#line 321 "ztrsv.f"
/* L130: */
#line 321 "ztrsv.f"
			}
#line 322 "ztrsv.f"
			if (nounit) {
#line 322 "ztrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 322 "ztrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 322 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 322 "ztrsv.f"
			}
#line 323 "ztrsv.f"
		    }
#line 324 "ztrsv.f"
		    i__2 = jx;
#line 324 "ztrsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 325 "ztrsv.f"
		    jx += *incx;
#line 326 "ztrsv.f"
/* L140: */
#line 326 "ztrsv.f"
		}
#line 327 "ztrsv.f"
	    }
#line 328 "ztrsv.f"
	} else {
#line 329 "ztrsv.f"
	    if (*incx == 1) {
#line 330 "ztrsv.f"
		for (j = *n; j >= 1; --j) {
#line 331 "ztrsv.f"
		    i__1 = j;
#line 331 "ztrsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 332 "ztrsv.f"
		    if (noconj) {
#line 333 "ztrsv.f"
			i__1 = j + 1;
#line 333 "ztrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 334 "ztrsv.f"
			    i__2 = i__ + j * a_dim1;
#line 334 "ztrsv.f"
			    i__3 = i__;
#line 334 "ztrsv.f"
			    z__2.r = a[i__2].r * x[i__3].r - a[i__2].i * x[
				    i__3].i, z__2.i = a[i__2].r * x[i__3].i + 
				    a[i__2].i * x[i__3].r;
#line 334 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 334 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 335 "ztrsv.f"
/* L150: */
#line 335 "ztrsv.f"
			}
#line 336 "ztrsv.f"
			if (nounit) {
#line 336 "ztrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 336 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 336 "ztrsv.f"
			}
#line 337 "ztrsv.f"
		    } else {
#line 338 "ztrsv.f"
			i__1 = j + 1;
#line 338 "ztrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 339 "ztrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 339 "ztrsv.f"
			    i__2 = i__;
#line 339 "ztrsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 339 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 339 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 340 "ztrsv.f"
/* L160: */
#line 340 "ztrsv.f"
			}
#line 341 "ztrsv.f"
			if (nounit) {
#line 341 "ztrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 341 "ztrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 341 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 341 "ztrsv.f"
			}
#line 342 "ztrsv.f"
		    }
#line 343 "ztrsv.f"
		    i__1 = j;
#line 343 "ztrsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 344 "ztrsv.f"
/* L170: */
#line 344 "ztrsv.f"
		}
#line 345 "ztrsv.f"
	    } else {
#line 346 "ztrsv.f"
		kx += (*n - 1) * *incx;
#line 347 "ztrsv.f"
		jx = kx;
#line 348 "ztrsv.f"
		for (j = *n; j >= 1; --j) {
#line 349 "ztrsv.f"
		    ix = kx;
#line 350 "ztrsv.f"
		    i__1 = jx;
#line 350 "ztrsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 351 "ztrsv.f"
		    if (noconj) {
#line 352 "ztrsv.f"
			i__1 = j + 1;
#line 352 "ztrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 353 "ztrsv.f"
			    i__2 = i__ + j * a_dim1;
#line 353 "ztrsv.f"
			    i__3 = ix;
#line 353 "ztrsv.f"
			    z__2.r = a[i__2].r * x[i__3].r - a[i__2].i * x[
				    i__3].i, z__2.i = a[i__2].r * x[i__3].i + 
				    a[i__2].i * x[i__3].r;
#line 353 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 353 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ztrsv.f"
			    ix -= *incx;
#line 355 "ztrsv.f"
/* L180: */
#line 355 "ztrsv.f"
			}
#line 356 "ztrsv.f"
			if (nounit) {
#line 356 "ztrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 356 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 356 "ztrsv.f"
			}
#line 357 "ztrsv.f"
		    } else {
#line 358 "ztrsv.f"
			i__1 = j + 1;
#line 358 "ztrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 359 "ztrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 359 "ztrsv.f"
			    i__2 = ix;
#line 359 "ztrsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 359 "ztrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 359 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 360 "ztrsv.f"
			    ix -= *incx;
#line 361 "ztrsv.f"
/* L190: */
#line 361 "ztrsv.f"
			}
#line 362 "ztrsv.f"
			if (nounit) {
#line 362 "ztrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 362 "ztrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 362 "ztrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 362 "ztrsv.f"
			}
#line 363 "ztrsv.f"
		    }
#line 364 "ztrsv.f"
		    i__1 = jx;
#line 364 "ztrsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 365 "ztrsv.f"
		    jx -= *incx;
#line 366 "ztrsv.f"
/* L200: */
#line 366 "ztrsv.f"
		}
#line 367 "ztrsv.f"
	    }
#line 368 "ztrsv.f"
	}
#line 369 "ztrsv.f"
    }

#line 371 "ztrsv.f"
    return 0;

/*     End of ZTRSV . */

} /* ztrsv_ */

