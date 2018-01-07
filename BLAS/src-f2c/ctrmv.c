#line 1 "ctrmv.f"
/* ctrmv.f -- translated by f2c (version 20100827).
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

#line 1 "ctrmv.f"
/* > \brief \b CTRMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRMV  performs one of the matrix-vector operations */
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
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
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
/* Subroutine */ int ctrmv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 187 "ctrmv.f"
    /* Parameter adjustments */
#line 187 "ctrmv.f"
    a_dim1 = *lda;
#line 187 "ctrmv.f"
    a_offset = 1 + a_dim1;
#line 187 "ctrmv.f"
    a -= a_offset;
#line 187 "ctrmv.f"
    --x;
#line 187 "ctrmv.f"

#line 187 "ctrmv.f"
    /* Function Body */
#line 187 "ctrmv.f"
    info = 0;
#line 188 "ctrmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 189 "ctrmv.f"
	info = 1;
#line 190 "ctrmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 192 "ctrmv.f"
	info = 2;
#line 193 "ctrmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 194 "ctrmv.f"
	info = 3;
#line 195 "ctrmv.f"
    } else if (*n < 0) {
#line 196 "ctrmv.f"
	info = 4;
#line 197 "ctrmv.f"
    } else if (*lda < max(1,*n)) {
#line 198 "ctrmv.f"
	info = 6;
#line 199 "ctrmv.f"
    } else if (*incx == 0) {
#line 200 "ctrmv.f"
	info = 8;
#line 201 "ctrmv.f"
    }
#line 202 "ctrmv.f"
    if (info != 0) {
#line 203 "ctrmv.f"
	xerbla_("CTRMV ", &info, (ftnlen)6);
#line 204 "ctrmv.f"
	return 0;
#line 205 "ctrmv.f"
    }

/*     Quick return if possible. */

#line 209 "ctrmv.f"
    if (*n == 0) {
#line 209 "ctrmv.f"
	return 0;
#line 209 "ctrmv.f"
    }

#line 211 "ctrmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 212 "ctrmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 217 "ctrmv.f"
    if (*incx <= 0) {
#line 218 "ctrmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 219 "ctrmv.f"
    } else if (*incx != 1) {
#line 220 "ctrmv.f"
	kx = 1;
#line 221 "ctrmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 226 "ctrmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := A*x. */

#line 230 "ctrmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 231 "ctrmv.f"
	    if (*incx == 1) {
#line 232 "ctrmv.f"
		i__1 = *n;
#line 232 "ctrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 233 "ctrmv.f"
		    i__2 = j;
#line 233 "ctrmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 234 "ctrmv.f"
			i__2 = j;
#line 234 "ctrmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 235 "ctrmv.f"
			i__2 = j - 1;
#line 235 "ctrmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 236 "ctrmv.f"
			    i__3 = i__;
#line 236 "ctrmv.f"
			    i__4 = i__;
#line 236 "ctrmv.f"
			    i__5 = i__ + j * a_dim1;
#line 236 "ctrmv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 236 "ctrmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 236 "ctrmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 237 "ctrmv.f"
/* L10: */
#line 237 "ctrmv.f"
			}
#line 238 "ctrmv.f"
			if (nounit) {
#line 238 "ctrmv.f"
			    i__2 = j;
#line 238 "ctrmv.f"
			    i__3 = j;
#line 238 "ctrmv.f"
			    i__4 = j + j * a_dim1;
#line 238 "ctrmv.f"
			    z__1.r = x[i__3].r * a[i__4].r - x[i__3].i * a[
				    i__4].i, z__1.i = x[i__3].r * a[i__4].i + 
				    x[i__3].i * a[i__4].r;
#line 238 "ctrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 238 "ctrmv.f"
			}
#line 239 "ctrmv.f"
		    }
#line 240 "ctrmv.f"
/* L20: */
#line 240 "ctrmv.f"
		}
#line 241 "ctrmv.f"
	    } else {
#line 242 "ctrmv.f"
		jx = kx;
#line 243 "ctrmv.f"
		i__1 = *n;
#line 243 "ctrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 244 "ctrmv.f"
		    i__2 = jx;
#line 244 "ctrmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 245 "ctrmv.f"
			i__2 = jx;
#line 245 "ctrmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 246 "ctrmv.f"
			ix = kx;
#line 247 "ctrmv.f"
			i__2 = j - 1;
#line 247 "ctrmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 248 "ctrmv.f"
			    i__3 = ix;
#line 248 "ctrmv.f"
			    i__4 = ix;
#line 248 "ctrmv.f"
			    i__5 = i__ + j * a_dim1;
#line 248 "ctrmv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 248 "ctrmv.f"
			    z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + 
				    z__2.i;
#line 248 "ctrmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 249 "ctrmv.f"
			    ix += *incx;
#line 250 "ctrmv.f"
/* L30: */
#line 250 "ctrmv.f"
			}
#line 251 "ctrmv.f"
			if (nounit) {
#line 251 "ctrmv.f"
			    i__2 = jx;
#line 251 "ctrmv.f"
			    i__3 = jx;
#line 251 "ctrmv.f"
			    i__4 = j + j * a_dim1;
#line 251 "ctrmv.f"
			    z__1.r = x[i__3].r * a[i__4].r - x[i__3].i * a[
				    i__4].i, z__1.i = x[i__3].r * a[i__4].i + 
				    x[i__3].i * a[i__4].r;
#line 251 "ctrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 251 "ctrmv.f"
			}
#line 252 "ctrmv.f"
		    }
#line 253 "ctrmv.f"
		    jx += *incx;
#line 254 "ctrmv.f"
/* L40: */
#line 254 "ctrmv.f"
		}
#line 255 "ctrmv.f"
	    }
#line 256 "ctrmv.f"
	} else {
#line 257 "ctrmv.f"
	    if (*incx == 1) {
#line 258 "ctrmv.f"
		for (j = *n; j >= 1; --j) {
#line 259 "ctrmv.f"
		    i__1 = j;
#line 259 "ctrmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 260 "ctrmv.f"
			i__1 = j;
#line 260 "ctrmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 261 "ctrmv.f"
			i__1 = j + 1;
#line 261 "ctrmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 262 "ctrmv.f"
			    i__2 = i__;
#line 262 "ctrmv.f"
			    i__3 = i__;
#line 262 "ctrmv.f"
			    i__4 = i__ + j * a_dim1;
#line 262 "ctrmv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 262 "ctrmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 262 "ctrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 263 "ctrmv.f"
/* L50: */
#line 263 "ctrmv.f"
			}
#line 264 "ctrmv.f"
			if (nounit) {
#line 264 "ctrmv.f"
			    i__1 = j;
#line 264 "ctrmv.f"
			    i__2 = j;
#line 264 "ctrmv.f"
			    i__3 = j + j * a_dim1;
#line 264 "ctrmv.f"
			    z__1.r = x[i__2].r * a[i__3].r - x[i__2].i * a[
				    i__3].i, z__1.i = x[i__2].r * a[i__3].i + 
				    x[i__2].i * a[i__3].r;
#line 264 "ctrmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 264 "ctrmv.f"
			}
#line 265 "ctrmv.f"
		    }
#line 266 "ctrmv.f"
/* L60: */
#line 266 "ctrmv.f"
		}
#line 267 "ctrmv.f"
	    } else {
#line 268 "ctrmv.f"
		kx += (*n - 1) * *incx;
#line 269 "ctrmv.f"
		jx = kx;
#line 270 "ctrmv.f"
		for (j = *n; j >= 1; --j) {
#line 271 "ctrmv.f"
		    i__1 = jx;
#line 271 "ctrmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 272 "ctrmv.f"
			i__1 = jx;
#line 272 "ctrmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 273 "ctrmv.f"
			ix = kx;
#line 274 "ctrmv.f"
			i__1 = j + 1;
#line 274 "ctrmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 275 "ctrmv.f"
			    i__2 = ix;
#line 275 "ctrmv.f"
			    i__3 = ix;
#line 275 "ctrmv.f"
			    i__4 = i__ + j * a_dim1;
#line 275 "ctrmv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 275 "ctrmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 275 "ctrmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 276 "ctrmv.f"
			    ix -= *incx;
#line 277 "ctrmv.f"
/* L70: */
#line 277 "ctrmv.f"
			}
#line 278 "ctrmv.f"
			if (nounit) {
#line 278 "ctrmv.f"
			    i__1 = jx;
#line 278 "ctrmv.f"
			    i__2 = jx;
#line 278 "ctrmv.f"
			    i__3 = j + j * a_dim1;
#line 278 "ctrmv.f"
			    z__1.r = x[i__2].r * a[i__3].r - x[i__2].i * a[
				    i__3].i, z__1.i = x[i__2].r * a[i__3].i + 
				    x[i__2].i * a[i__3].r;
#line 278 "ctrmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 278 "ctrmv.f"
			}
#line 279 "ctrmv.f"
		    }
#line 280 "ctrmv.f"
		    jx -= *incx;
#line 281 "ctrmv.f"
/* L80: */
#line 281 "ctrmv.f"
		}
#line 282 "ctrmv.f"
	    }
#line 283 "ctrmv.f"
	}
#line 284 "ctrmv.f"
    } else {

/*        Form  x := A**T*x  or  x := A**H*x. */

#line 288 "ctrmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 289 "ctrmv.f"
	    if (*incx == 1) {
#line 290 "ctrmv.f"
		for (j = *n; j >= 1; --j) {
#line 291 "ctrmv.f"
		    i__1 = j;
#line 291 "ctrmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 292 "ctrmv.f"
		    if (noconj) {
#line 293 "ctrmv.f"
			if (nounit) {
#line 293 "ctrmv.f"
			    i__1 = j + j * a_dim1;
#line 293 "ctrmv.f"
			    z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				    z__1.i = temp.r * a[i__1].i + temp.i * a[
				    i__1].r;
#line 293 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 293 "ctrmv.f"
			}
#line 294 "ctrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 295 "ctrmv.f"
			    i__1 = i__ + j * a_dim1;
#line 295 "ctrmv.f"
			    i__2 = i__;
#line 295 "ctrmv.f"
			    z__2.r = a[i__1].r * x[i__2].r - a[i__1].i * x[
				    i__2].i, z__2.i = a[i__1].r * x[i__2].i + 
				    a[i__1].i * x[i__2].r;
#line 295 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 295 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 296 "ctrmv.f"
/* L90: */
#line 296 "ctrmv.f"
			}
#line 297 "ctrmv.f"
		    } else {
#line 298 "ctrmv.f"
			if (nounit) {
#line 298 "ctrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 298 "ctrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 298 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 298 "ctrmv.f"
			}
#line 299 "ctrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 300 "ctrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 300 "ctrmv.f"
			    i__1 = i__;
#line 300 "ctrmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 300 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 300 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ctrmv.f"
/* L100: */
#line 301 "ctrmv.f"
			}
#line 302 "ctrmv.f"
		    }
#line 303 "ctrmv.f"
		    i__1 = j;
#line 303 "ctrmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 304 "ctrmv.f"
/* L110: */
#line 304 "ctrmv.f"
		}
#line 305 "ctrmv.f"
	    } else {
#line 306 "ctrmv.f"
		jx = kx + (*n - 1) * *incx;
#line 307 "ctrmv.f"
		for (j = *n; j >= 1; --j) {
#line 308 "ctrmv.f"
		    i__1 = jx;
#line 308 "ctrmv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 309 "ctrmv.f"
		    ix = jx;
#line 310 "ctrmv.f"
		    if (noconj) {
#line 311 "ctrmv.f"
			if (nounit) {
#line 311 "ctrmv.f"
			    i__1 = j + j * a_dim1;
#line 311 "ctrmv.f"
			    z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				    z__1.i = temp.r * a[i__1].i + temp.i * a[
				    i__1].r;
#line 311 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 311 "ctrmv.f"
			}
#line 312 "ctrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 313 "ctrmv.f"
			    ix -= *incx;
#line 314 "ctrmv.f"
			    i__1 = i__ + j * a_dim1;
#line 314 "ctrmv.f"
			    i__2 = ix;
#line 314 "ctrmv.f"
			    z__2.r = a[i__1].r * x[i__2].r - a[i__1].i * x[
				    i__2].i, z__2.i = a[i__1].r * x[i__2].i + 
				    a[i__1].i * x[i__2].r;
#line 314 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 314 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 315 "ctrmv.f"
/* L120: */
#line 315 "ctrmv.f"
			}
#line 316 "ctrmv.f"
		    } else {
#line 317 "ctrmv.f"
			if (nounit) {
#line 317 "ctrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 317 "ctrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 317 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 317 "ctrmv.f"
			}
#line 318 "ctrmv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 319 "ctrmv.f"
			    ix -= *incx;
#line 320 "ctrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 320 "ctrmv.f"
			    i__1 = ix;
#line 320 "ctrmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 320 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 320 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 321 "ctrmv.f"
/* L130: */
#line 321 "ctrmv.f"
			}
#line 322 "ctrmv.f"
		    }
#line 323 "ctrmv.f"
		    i__1 = jx;
#line 323 "ctrmv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 324 "ctrmv.f"
		    jx -= *incx;
#line 325 "ctrmv.f"
/* L140: */
#line 325 "ctrmv.f"
		}
#line 326 "ctrmv.f"
	    }
#line 327 "ctrmv.f"
	} else {
#line 328 "ctrmv.f"
	    if (*incx == 1) {
#line 329 "ctrmv.f"
		i__1 = *n;
#line 329 "ctrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 330 "ctrmv.f"
		    i__2 = j;
#line 330 "ctrmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 331 "ctrmv.f"
		    if (noconj) {
#line 332 "ctrmv.f"
			if (nounit) {
#line 332 "ctrmv.f"
			    i__2 = j + j * a_dim1;
#line 332 "ctrmv.f"
			    z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__1.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 332 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 332 "ctrmv.f"
			}
#line 333 "ctrmv.f"
			i__2 = *n;
#line 333 "ctrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 334 "ctrmv.f"
			    i__3 = i__ + j * a_dim1;
#line 334 "ctrmv.f"
			    i__4 = i__;
#line 334 "ctrmv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 334 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 334 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 335 "ctrmv.f"
/* L150: */
#line 335 "ctrmv.f"
			}
#line 336 "ctrmv.f"
		    } else {
#line 337 "ctrmv.f"
			if (nounit) {
#line 337 "ctrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 337 "ctrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 337 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 337 "ctrmv.f"
			}
#line 338 "ctrmv.f"
			i__2 = *n;
#line 338 "ctrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 339 "ctrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 339 "ctrmv.f"
			    i__3 = i__;
#line 339 "ctrmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 339 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 339 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 340 "ctrmv.f"
/* L160: */
#line 340 "ctrmv.f"
			}
#line 341 "ctrmv.f"
		    }
#line 342 "ctrmv.f"
		    i__2 = j;
#line 342 "ctrmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 343 "ctrmv.f"
/* L170: */
#line 343 "ctrmv.f"
		}
#line 344 "ctrmv.f"
	    } else {
#line 345 "ctrmv.f"
		jx = kx;
#line 346 "ctrmv.f"
		i__1 = *n;
#line 346 "ctrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 347 "ctrmv.f"
		    i__2 = jx;
#line 347 "ctrmv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 348 "ctrmv.f"
		    ix = jx;
#line 349 "ctrmv.f"
		    if (noconj) {
#line 350 "ctrmv.f"
			if (nounit) {
#line 350 "ctrmv.f"
			    i__2 = j + j * a_dim1;
#line 350 "ctrmv.f"
			    z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__1.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 350 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 350 "ctrmv.f"
			}
#line 351 "ctrmv.f"
			i__2 = *n;
#line 351 "ctrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 352 "ctrmv.f"
			    ix += *incx;
#line 353 "ctrmv.f"
			    i__3 = i__ + j * a_dim1;
#line 353 "ctrmv.f"
			    i__4 = ix;
#line 353 "ctrmv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 353 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 353 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ctrmv.f"
/* L180: */
#line 354 "ctrmv.f"
			}
#line 355 "ctrmv.f"
		    } else {
#line 356 "ctrmv.f"
			if (nounit) {
#line 356 "ctrmv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 356 "ctrmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 356 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 356 "ctrmv.f"
			}
#line 357 "ctrmv.f"
			i__2 = *n;
#line 357 "ctrmv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 358 "ctrmv.f"
			    ix += *incx;
#line 359 "ctrmv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 359 "ctrmv.f"
			    i__3 = ix;
#line 359 "ctrmv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 359 "ctrmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 359 "ctrmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 360 "ctrmv.f"
/* L190: */
#line 360 "ctrmv.f"
			}
#line 361 "ctrmv.f"
		    }
#line 362 "ctrmv.f"
		    i__2 = jx;
#line 362 "ctrmv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 363 "ctrmv.f"
		    jx += *incx;
#line 364 "ctrmv.f"
/* L200: */
#line 364 "ctrmv.f"
		}
#line 365 "ctrmv.f"
	    }
#line 366 "ctrmv.f"
	}
#line 367 "ctrmv.f"
    }

#line 369 "ctrmv.f"
    return 0;

/*     End of CTRMV . */

} /* ctrmv_ */

