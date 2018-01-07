#line 1 "ztrmv.f"
/* ztrmv.f -- translated by f2c (version 20100827).
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

#line 1 "ztrmv.f"
/* > \brief \b ZTRMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

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
/* > ZTRMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x,   or   x := A**H*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix. */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is (input/output) COMPLEX*16 array of dimension at least */
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
/* Subroutine */ int ztrmv_(char *uplo, char *trans, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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

#line 187 "ztrmv.f"
    /* Parameter adjustments */
#line 187 "ztrmv.f"
    a_dim1 = *lda;
#line 187 "ztrmv.f"
    a_offset = 1 + a_dim1;
#line 187 "ztrmv.f"
    a -= a_offset;
#line 187 "ztrmv.f"
    --x;
#line 187 "ztrmv.f"

#line 187 "ztrmv.f"
    /* Function Body */
#line 187 "ztrmv.f"
    info = 0;
#line 188 "ztrmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 189 "ztrmv.f"
	info = 1;
#line 190 "ztrmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 192 "ztrmv.f"
	info = 2;
#line 193 "ztrmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 194 "ztrmv.f"
	info = 3;
#line 195 "ztrmv.f"
    } else if (*n < 0) {
#line 196 "ztrmv.f"
	info = 4;
#line 197 "ztrmv.f"
    } else if (*lda < max(1,*n)) {
#line 198 "ztrmv.f"
	info = 6;
#line 199 "ztrmv.f"
    } else if (*incx == 0) {
#line 200 "ztrmv.f"
	info = 8;
#line 201 "ztrmv.f"
    }
#line 202 "ztrmv.f"
    if (info != 0) {
#line 203 "ztrmv.f"
	xerbla_("ZTRMV ", &info, (ftnlen)6);
#line 204 "ztrmv.f"
	return 0;
#line 205 "ztrmv.f"
    }

/*     Quick return if possible. */

#line 209 "ztrmv.f"
    if (*n == 0) {
#line 209 "ztrmv.f"
	return 0;
#line 209 "ztrmv.f"
    }

#line 211 "ztrmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 212 "ztrmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 217 "ztrmv.f"
    if (*incx <= 0) {
#line 218 "ztrmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 219 "ztrmv.f"
    } else if (*incx != 1) {
#line 220 "ztrmv.f"
	kx = 1;
#line 221 "ztrmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 226 "ztrmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := A*x. */

#line 230 "ztrmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 231 "ztrmv.f"
	    if (*incx == 1) {
#line 232 "ztrmv.f"
		i__1 = *n;
#line 232 "ztrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 233 "ztrmv.f"
		    i__2 = j;
#line 233 "ztrmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 234 "ztrmv.f"
			i__2 = j;
#line 234 "ztrmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 235 "ztrmv.f"
			i__2 = j - 1;
#line 235 "ztrmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 236 "ztrmv.f"
			    i__3 = i__;
#line 236 "ztrmv.f"
			    i__4 = i__;
#line 236 "ztrmv.f"
			    i__5 = i__ + j * a_dim1;
#line 236 "ztrmv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 236 "ztrmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 236 "ztrmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 237 "ztrmv.f"
/* L10: */
#line 237 "ztrmv.f"
			}
#line 238 "ztrmv.f"
			if (nounit) {
#line 238 "ztrmv.f"
			    i__2 = j;
#line 238 "ztrmv.f"
			    i__3 = j;
#line 238 "ztrmv.f"
			    i__4 = j + j * a_dim1;
#line 238 "ztrmv.f"
			    z__1.r = x[i__3].r * a[i__4].r - x[i__3].i * a[
				    i__4].i, z__1.i = x[i__3].r * a[i__4].i + 
				    x[i__3].i * a[i__4].r;
#line 238 "ztrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 238 "ztrmv.f"
			}
#line 239 "ztrmv.f"
		    }
#line 240 "ztrmv.f"
/* L20: */
#line 240 "ztrmv.f"
		}
#line 241 "ztrmv.f"
	    } else {
#line 242 "ztrmv.f"
		jx = kx;
#line 243 "ztrmv.f"
		i__1 = *n;
#line 243 "ztrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 244 "ztrmv.f"
		    i__2 = jx;
#line 244 "ztrmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 245 "ztrmv.f"
			i__2 = jx;
#line 245 "ztrmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 246 "ztrmv.f"
			ix = kx;
#line 247 "ztrmv.f"
			i__2 = j - 1;
#line 247 "ztrmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 248 "ztrmv.f"
			    i__3 = ix;
#line 248 "ztrmv.f"
			    i__4 = ix;
#line 248 "ztrmv.f"
			    i__5 = i__ + j * a_dim1;
#line 248 "ztrmv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 248 "ztrmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 248 "ztrmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 249 "ztrmv.f"
			    ix += *incx;
#line 250 "ztrmv.f"
/* L30: */
#line 250 "ztrmv.f"
			}
#line 251 "ztrmv.f"
			if (nounit) {
#line 251 "ztrmv.f"
			    i__2 = jx;
#line 251 "ztrmv.f"
			    i__3 = jx;
#line 251 "ztrmv.f"
			    i__4 = j + j * a_dim1;
#line 251 "ztrmv.f"
			    z__1.r = x[i__3].r * a[i__4].r - x[i__3].i * a[
				    i__4].i, z__1.i = x[i__3].r * a[i__4].i + 
				    x[i__3].i * a[i__4].r;
#line 251 "ztrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 251 "ztrmv.f"
			}
#line 252 "ztrmv.f"
		    }
#line 253 "ztrmv.f"
		    jx += *incx;
#line 254 "ztrmv.f"
/* L40: */
#line 254 "ztrmv.f"
		}
#line 255 "ztrmv.f"
	    }
#line 256 "ztrmv.f"
	} else {
#line 257 "ztrmv.f"
	    if (*incx == 1) {
#line 258 "ztrmv.f"
		for (j = *n; j >= 1; --j) {
#line 259 "ztrmv.f"
		    i__1 = j;
#line 259 "ztrmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 260 "ztrmv.f"
			i__1 = j;
#line 260 "ztrmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 261 "ztrmv.f"
			i__1 = j + 1;
#line 261 "ztrmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 262 "ztrmv.f"
			    i__2 = i__;
#line 262 "ztrmv.f"
			    i__3 = i__;
#line 262 "ztrmv.f"
			    i__4 = i__ + j * a_dim1;
#line 262 "ztrmv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 262 "ztrmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 262 "ztrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 263 "ztrmv.f"
/* L50: */
#line 263 "ztrmv.f"
			}
#line 264 "ztrmv.f"
			if (nounit) {
#line 264 "ztrmv.f"
			    i__1 = j;
#line 264 "ztrmv.f"
			    i__2 = j;
#line 264 "ztrmv.f"
			    i__3 = j + j * a_dim1;
#line 264 "ztrmv.f"
			    z__1.r = x[i__2].r * a[i__3].r - x[i__2].i * a[
				    i__3].i, z__1.i = x[i__2].r * a[i__3].i + 
				    x[i__2].i * a[i__3].r;
#line 264 "ztrmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 264 "ztrmv.f"
			}
#line 265 "ztrmv.f"
		    }
#line 266 "ztrmv.f"
/* L60: */
#line 266 "ztrmv.f"
		}
#line 267 "ztrmv.f"
	    } else {
#line 268 "ztrmv.f"
		kx += (*n - 1) * *incx;
#line 269 "ztrmv.f"
		jx = kx;
#line 270 "ztrmv.f"
		for (j = *n; j >= 1; --j) {
#line 271 "ztrmv.f"
		    i__1 = jx;
#line 271 "ztrmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 272 "ztrmv.f"
			i__1 = jx;
#line 272 "ztrmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 273 "ztrmv.f"
			ix = kx;
#line 274 "ztrmv.f"
			i__1 = j + 1;
#line 274 "ztrmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 275 "ztrmv.f"
			    i__2 = ix;
#line 275 "ztrmv.f"
			    i__3 = ix;
#line 275 "ztrmv.f"
			    i__4 = i__ + j * a_dim1;
#line 275 "ztrmv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 275 "ztrmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 275 "ztrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 276 "ztrmv.f"
			    ix -= *incx;
#line 277 "ztrmv.f"
/* L70: */
#line 277 "ztrmv.f"
			}
#line 278 "ztrmv.f"
			if (nounit) {
#line 278 "ztrmv.f"
			    i__1 = jx;
#line 278 "ztrmv.f"
			    i__2 = jx;
#line 278 "ztrmv.f"
			    i__3 = j + j * a_dim1;
#line 278 "ztrmv.f"
			    z__1.r = x[i__2].r * a[i__3].r - x[i__2].i * a[
				    i__3].i, z__1.i = x[i__2].r * a[i__3].i + 
				    x[i__2].i * a[i__3].r;
#line 278 "ztrmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 278 "ztrmv.f"
			}
#line 279 "ztrmv.f"
		    }
#line 280 "ztrmv.f"
		    jx -= *incx;
#line 281 "ztrmv.f"
/* L80: */
#line 281 "ztrmv.f"
		}
#line 282 "ztrmv.f"
	    }
#line 283 "ztrmv.f"
	}
#line 284 "ztrmv.f"
    } else {

/*        Form  x := A**T*x  or  x := A**H*x. */

#line 288 "ztrmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 289 "ztrmv.f"
	    if (*incx == 1) {
#line 290 "ztrmv.f"
		for (j = *n; j >= 1; --j) {
#line 291 "ztrmv.f"
		    i__1 = j;
#line 291 "ztrmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 292 "ztrmv.f"
		    if (noconj) {
#line 293 "ztrmv.f"
			if (nounit) {
#line 293 "ztrmv.f"
			    i__1 = j + j * a_dim1;
#line 293 "ztrmv.f"
			    z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				    z__1.i = temp.r * a[i__1].i + temp.i * a[
				    i__1].r;
#line 293 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 293 "ztrmv.f"
			}
#line 294 "ztrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 295 "ztrmv.f"
			    i__1 = i__ + j * a_dim1;
#line 295 "ztrmv.f"
			    i__2 = i__;
#line 295 "ztrmv.f"
			    z__2.r = a[i__1].r * x[i__2].r - a[i__1].i * x[
				    i__2].i, z__2.i = a[i__1].r * x[i__2].i + 
				    a[i__1].i * x[i__2].r;
#line 295 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 295 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 296 "ztrmv.f"
/* L90: */
#line 296 "ztrmv.f"
			}
#line 297 "ztrmv.f"
		    } else {
#line 298 "ztrmv.f"
			if (nounit) {
#line 298 "ztrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 298 "ztrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 298 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 298 "ztrmv.f"
			}
#line 299 "ztrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 300 "ztrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 300 "ztrmv.f"
			    i__1 = i__;
#line 300 "ztrmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 300 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 300 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ztrmv.f"
/* L100: */
#line 301 "ztrmv.f"
			}
#line 302 "ztrmv.f"
		    }
#line 303 "ztrmv.f"
		    i__1 = j;
#line 303 "ztrmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 304 "ztrmv.f"
/* L110: */
#line 304 "ztrmv.f"
		}
#line 305 "ztrmv.f"
	    } else {
#line 306 "ztrmv.f"
		jx = kx + (*n - 1) * *incx;
#line 307 "ztrmv.f"
		for (j = *n; j >= 1; --j) {
#line 308 "ztrmv.f"
		    i__1 = jx;
#line 308 "ztrmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 309 "ztrmv.f"
		    ix = jx;
#line 310 "ztrmv.f"
		    if (noconj) {
#line 311 "ztrmv.f"
			if (nounit) {
#line 311 "ztrmv.f"
			    i__1 = j + j * a_dim1;
#line 311 "ztrmv.f"
			    z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				    z__1.i = temp.r * a[i__1].i + temp.i * a[
				    i__1].r;
#line 311 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 311 "ztrmv.f"
			}
#line 312 "ztrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 313 "ztrmv.f"
			    ix -= *incx;
#line 314 "ztrmv.f"
			    i__1 = i__ + j * a_dim1;
#line 314 "ztrmv.f"
			    i__2 = ix;
#line 314 "ztrmv.f"
			    z__2.r = a[i__1].r * x[i__2].r - a[i__1].i * x[
				    i__2].i, z__2.i = a[i__1].r * x[i__2].i + 
				    a[i__1].i * x[i__2].r;
#line 314 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 314 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 315 "ztrmv.f"
/* L120: */
#line 315 "ztrmv.f"
			}
#line 316 "ztrmv.f"
		    } else {
#line 317 "ztrmv.f"
			if (nounit) {
#line 317 "ztrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 317 "ztrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 317 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 317 "ztrmv.f"
			}
#line 318 "ztrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 319 "ztrmv.f"
			    ix -= *incx;
#line 320 "ztrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 320 "ztrmv.f"
			    i__1 = ix;
#line 320 "ztrmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 320 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 320 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 321 "ztrmv.f"
/* L130: */
#line 321 "ztrmv.f"
			}
#line 322 "ztrmv.f"
		    }
#line 323 "ztrmv.f"
		    i__1 = jx;
#line 323 "ztrmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 324 "ztrmv.f"
		    jx -= *incx;
#line 325 "ztrmv.f"
/* L140: */
#line 325 "ztrmv.f"
		}
#line 326 "ztrmv.f"
	    }
#line 327 "ztrmv.f"
	} else {
#line 328 "ztrmv.f"
	    if (*incx == 1) {
#line 329 "ztrmv.f"
		i__1 = *n;
#line 329 "ztrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 330 "ztrmv.f"
		    i__2 = j;
#line 330 "ztrmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 331 "ztrmv.f"
		    if (noconj) {
#line 332 "ztrmv.f"
			if (nounit) {
#line 332 "ztrmv.f"
			    i__2 = j + j * a_dim1;
#line 332 "ztrmv.f"
			    z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__1.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 332 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 332 "ztrmv.f"
			}
#line 333 "ztrmv.f"
			i__2 = *n;
#line 333 "ztrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 334 "ztrmv.f"
			    i__3 = i__ + j * a_dim1;
#line 334 "ztrmv.f"
			    i__4 = i__;
#line 334 "ztrmv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 334 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 334 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 335 "ztrmv.f"
/* L150: */
#line 335 "ztrmv.f"
			}
#line 336 "ztrmv.f"
		    } else {
#line 337 "ztrmv.f"
			if (nounit) {
#line 337 "ztrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 337 "ztrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 337 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 337 "ztrmv.f"
			}
#line 338 "ztrmv.f"
			i__2 = *n;
#line 338 "ztrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 339 "ztrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 339 "ztrmv.f"
			    i__3 = i__;
#line 339 "ztrmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 339 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 339 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 340 "ztrmv.f"
/* L160: */
#line 340 "ztrmv.f"
			}
#line 341 "ztrmv.f"
		    }
#line 342 "ztrmv.f"
		    i__2 = j;
#line 342 "ztrmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 343 "ztrmv.f"
/* L170: */
#line 343 "ztrmv.f"
		}
#line 344 "ztrmv.f"
	    } else {
#line 345 "ztrmv.f"
		jx = kx;
#line 346 "ztrmv.f"
		i__1 = *n;
#line 346 "ztrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 347 "ztrmv.f"
		    i__2 = jx;
#line 347 "ztrmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 348 "ztrmv.f"
		    ix = jx;
#line 349 "ztrmv.f"
		    if (noconj) {
#line 350 "ztrmv.f"
			if (nounit) {
#line 350 "ztrmv.f"
			    i__2 = j + j * a_dim1;
#line 350 "ztrmv.f"
			    z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__1.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 350 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 350 "ztrmv.f"
			}
#line 351 "ztrmv.f"
			i__2 = *n;
#line 351 "ztrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 352 "ztrmv.f"
			    ix += *incx;
#line 353 "ztrmv.f"
			    i__3 = i__ + j * a_dim1;
#line 353 "ztrmv.f"
			    i__4 = ix;
#line 353 "ztrmv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 353 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 353 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ztrmv.f"
/* L180: */
#line 354 "ztrmv.f"
			}
#line 355 "ztrmv.f"
		    } else {
#line 356 "ztrmv.f"
			if (nounit) {
#line 356 "ztrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 356 "ztrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 356 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 356 "ztrmv.f"
			}
#line 357 "ztrmv.f"
			i__2 = *n;
#line 357 "ztrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 358 "ztrmv.f"
			    ix += *incx;
#line 359 "ztrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 359 "ztrmv.f"
			    i__3 = ix;
#line 359 "ztrmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 359 "ztrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 359 "ztrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 360 "ztrmv.f"
/* L190: */
#line 360 "ztrmv.f"
			}
#line 361 "ztrmv.f"
		    }
#line 362 "ztrmv.f"
		    i__2 = jx;
#line 362 "ztrmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 363 "ztrmv.f"
		    jx += *incx;
#line 364 "ztrmv.f"
/* L200: */
#line 364 "ztrmv.f"
		}
#line 365 "ztrmv.f"
	    }
#line 366 "ztrmv.f"
	}
#line 367 "ztrmv.f"
    }

#line 369 "ztrmv.f"
    return 0;

/*     End of ZTRMV . */

} /* ztrmv_ */

