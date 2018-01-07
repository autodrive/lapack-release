#line 1 "ctrsv.f"
/* ctrsv.f -- translated by f2c (version 20100827).
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

#line 1 "ctrsv.f"
/* > \brief \b CTRSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

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
/* > CTRSV  solves one of the systems of equations */
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
/* >          A is COMPLEX array, dimension ( LDA, N ) */
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
/* Subroutine */ int ctrsv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 189 "ctrsv.f"
    /* Parameter adjustments */
#line 189 "ctrsv.f"
    a_dim1 = *lda;
#line 189 "ctrsv.f"
    a_offset = 1 + a_dim1;
#line 189 "ctrsv.f"
    a -= a_offset;
#line 189 "ctrsv.f"
    --x;
#line 189 "ctrsv.f"

#line 189 "ctrsv.f"
    /* Function Body */
#line 189 "ctrsv.f"
    info = 0;
#line 190 "ctrsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 191 "ctrsv.f"
	info = 1;
#line 192 "ctrsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 194 "ctrsv.f"
	info = 2;
#line 195 "ctrsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 196 "ctrsv.f"
	info = 3;
#line 197 "ctrsv.f"
    } else if (*n < 0) {
#line 198 "ctrsv.f"
	info = 4;
#line 199 "ctrsv.f"
    } else if (*lda < max(1,*n)) {
#line 200 "ctrsv.f"
	info = 6;
#line 201 "ctrsv.f"
    } else if (*incx == 0) {
#line 202 "ctrsv.f"
	info = 8;
#line 203 "ctrsv.f"
    }
#line 204 "ctrsv.f"
    if (info != 0) {
#line 205 "ctrsv.f"
	xerbla_("CTRSV ", &info, (ftnlen)6);
#line 206 "ctrsv.f"
	return 0;
#line 207 "ctrsv.f"
    }

/*     Quick return if possible. */

#line 211 "ctrsv.f"
    if (*n == 0) {
#line 211 "ctrsv.f"
	return 0;
#line 211 "ctrsv.f"
    }

#line 213 "ctrsv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 214 "ctrsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 219 "ctrsv.f"
    if (*incx <= 0) {
#line 220 "ctrsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 221 "ctrsv.f"
    } else if (*incx != 1) {
#line 222 "ctrsv.f"
	kx = 1;
#line 223 "ctrsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 228 "ctrsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 232 "ctrsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 233 "ctrsv.f"
	    if (*incx == 1) {
#line 234 "ctrsv.f"
		for (j = *n; j >= 1; --j) {
#line 235 "ctrsv.f"
		    i__1 = j;
#line 235 "ctrsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 236 "ctrsv.f"
			if (nounit) {
#line 236 "ctrsv.f"
			    i__1 = j;
#line 236 "ctrsv.f"
			    z_div(&z__1, &x[j], &a[j + j * a_dim1]);
#line 236 "ctrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 236 "ctrsv.f"
			}
#line 237 "ctrsv.f"
			i__1 = j;
#line 237 "ctrsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 238 "ctrsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 239 "ctrsv.f"
			    i__1 = i__;
#line 239 "ctrsv.f"
			    i__2 = i__;
#line 239 "ctrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 239 "ctrsv.f"
			    z__2.r = temp.r * a[i__3].r - temp.i * a[i__3].i, 
				    z__2.i = temp.r * a[i__3].i + temp.i * a[
				    i__3].r;
#line 239 "ctrsv.f"
			    z__1.r = x[i__2].r - z__2.r, z__1.i = x[i__2].i - 
				    z__2.i;
#line 239 "ctrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 240 "ctrsv.f"
/* L10: */
#line 240 "ctrsv.f"
			}
#line 241 "ctrsv.f"
		    }
#line 242 "ctrsv.f"
/* L20: */
#line 242 "ctrsv.f"
		}
#line 243 "ctrsv.f"
	    } else {
#line 244 "ctrsv.f"
		jx = kx + (*n - 1) * *incx;
#line 245 "ctrsv.f"
		for (j = *n; j >= 1; --j) {
#line 246 "ctrsv.f"
		    i__1 = jx;
#line 246 "ctrsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 247 "ctrsv.f"
			if (nounit) {
#line 247 "ctrsv.f"
			    i__1 = jx;
#line 247 "ctrsv.f"
			    z_div(&z__1, &x[jx], &a[j + j * a_dim1]);
#line 247 "ctrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 247 "ctrsv.f"
			}
#line 248 "ctrsv.f"
			i__1 = jx;
#line 248 "ctrsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 249 "ctrsv.f"
			ix = jx;
#line 250 "ctrsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 251 "ctrsv.f"
			    ix -= *incx;
#line 252 "ctrsv.f"
			    i__1 = ix;
#line 252 "ctrsv.f"
			    i__2 = ix;
#line 252 "ctrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 252 "ctrsv.f"
			    z__2.r = temp.r * a[i__3].r - temp.i * a[i__3].i, 
				    z__2.i = temp.r * a[i__3].i + temp.i * a[
				    i__3].r;
#line 252 "ctrsv.f"
			    z__1.r = x[i__2].r - z__2.r, z__1.i = x[i__2].i - 
				    z__2.i;
#line 252 "ctrsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 253 "ctrsv.f"
/* L30: */
#line 253 "ctrsv.f"
			}
#line 254 "ctrsv.f"
		    }
#line 255 "ctrsv.f"
		    jx -= *incx;
#line 256 "ctrsv.f"
/* L40: */
#line 256 "ctrsv.f"
		}
#line 257 "ctrsv.f"
	    }
#line 258 "ctrsv.f"
	} else {
#line 259 "ctrsv.f"
	    if (*incx == 1) {
#line 260 "ctrsv.f"
		i__1 = *n;
#line 260 "ctrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 261 "ctrsv.f"
		    i__2 = j;
#line 261 "ctrsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 262 "ctrsv.f"
			if (nounit) {
#line 262 "ctrsv.f"
			    i__2 = j;
#line 262 "ctrsv.f"
			    z_div(&z__1, &x[j], &a[j + j * a_dim1]);
#line 262 "ctrsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 262 "ctrsv.f"
			}
#line 263 "ctrsv.f"
			i__2 = j;
#line 263 "ctrsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 264 "ctrsv.f"
			i__2 = *n;
#line 264 "ctrsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 265 "ctrsv.f"
			    i__3 = i__;
#line 265 "ctrsv.f"
			    i__4 = i__;
#line 265 "ctrsv.f"
			    i__5 = i__ + j * a_dim1;
#line 265 "ctrsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 265 "ctrsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 265 "ctrsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 266 "ctrsv.f"
/* L50: */
#line 266 "ctrsv.f"
			}
#line 267 "ctrsv.f"
		    }
#line 268 "ctrsv.f"
/* L60: */
#line 268 "ctrsv.f"
		}
#line 269 "ctrsv.f"
	    } else {
#line 270 "ctrsv.f"
		jx = kx;
#line 271 "ctrsv.f"
		i__1 = *n;
#line 271 "ctrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 272 "ctrsv.f"
		    i__2 = jx;
#line 272 "ctrsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 273 "ctrsv.f"
			if (nounit) {
#line 273 "ctrsv.f"
			    i__2 = jx;
#line 273 "ctrsv.f"
			    z_div(&z__1, &x[jx], &a[j + j * a_dim1]);
#line 273 "ctrsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 273 "ctrsv.f"
			}
#line 274 "ctrsv.f"
			i__2 = jx;
#line 274 "ctrsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 275 "ctrsv.f"
			ix = jx;
#line 276 "ctrsv.f"
			i__2 = *n;
#line 276 "ctrsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 277 "ctrsv.f"
			    ix += *incx;
#line 278 "ctrsv.f"
			    i__3 = ix;
#line 278 "ctrsv.f"
			    i__4 = ix;
#line 278 "ctrsv.f"
			    i__5 = i__ + j * a_dim1;
#line 278 "ctrsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 278 "ctrsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 278 "ctrsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 279 "ctrsv.f"
/* L70: */
#line 279 "ctrsv.f"
			}
#line 280 "ctrsv.f"
		    }
#line 281 "ctrsv.f"
		    jx += *incx;
#line 282 "ctrsv.f"
/* L80: */
#line 282 "ctrsv.f"
		}
#line 283 "ctrsv.f"
	    }
#line 284 "ctrsv.f"
	}
#line 285 "ctrsv.f"
    } else {

/*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x. */

#line 289 "ctrsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 290 "ctrsv.f"
	    if (*incx == 1) {
#line 291 "ctrsv.f"
		i__1 = *n;
#line 291 "ctrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "ctrsv.f"
		    i__2 = j;
#line 292 "ctrsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 293 "ctrsv.f"
		    if (noconj) {
#line 294 "ctrsv.f"
			i__2 = j - 1;
#line 294 "ctrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "ctrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 295 "ctrsv.f"
			    i__4 = i__;
#line 295 "ctrsv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 295 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 295 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 296 "ctrsv.f"
/* L90: */
#line 296 "ctrsv.f"
			}
#line 297 "ctrsv.f"
			if (nounit) {
#line 297 "ctrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 297 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 297 "ctrsv.f"
			}
#line 298 "ctrsv.f"
		    } else {
#line 299 "ctrsv.f"
			i__2 = j - 1;
#line 299 "ctrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 300 "ctrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 300 "ctrsv.f"
			    i__3 = i__;
#line 300 "ctrsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 300 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 300 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 301 "ctrsv.f"
/* L100: */
#line 301 "ctrsv.f"
			}
#line 302 "ctrsv.f"
			if (nounit) {
#line 302 "ctrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 302 "ctrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 302 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 302 "ctrsv.f"
			}
#line 303 "ctrsv.f"
		    }
#line 304 "ctrsv.f"
		    i__2 = j;
#line 304 "ctrsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 305 "ctrsv.f"
/* L110: */
#line 305 "ctrsv.f"
		}
#line 306 "ctrsv.f"
	    } else {
#line 307 "ctrsv.f"
		jx = kx;
#line 308 "ctrsv.f"
		i__1 = *n;
#line 308 "ctrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 309 "ctrsv.f"
		    ix = kx;
#line 310 "ctrsv.f"
		    i__2 = jx;
#line 310 "ctrsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 311 "ctrsv.f"
		    if (noconj) {
#line 312 "ctrsv.f"
			i__2 = j - 1;
#line 312 "ctrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 313 "ctrsv.f"
			    i__3 = i__ + j * a_dim1;
#line 313 "ctrsv.f"
			    i__4 = ix;
#line 313 "ctrsv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 313 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 313 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 314 "ctrsv.f"
			    ix += *incx;
#line 315 "ctrsv.f"
/* L120: */
#line 315 "ctrsv.f"
			}
#line 316 "ctrsv.f"
			if (nounit) {
#line 316 "ctrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 316 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 316 "ctrsv.f"
			}
#line 317 "ctrsv.f"
		    } else {
#line 318 "ctrsv.f"
			i__2 = j - 1;
#line 318 "ctrsv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "ctrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 319 "ctrsv.f"
			    i__3 = ix;
#line 319 "ctrsv.f"
			    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				    z__2.i = z__3.r * x[i__3].i + z__3.i * x[
				    i__3].r;
#line 319 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 319 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 320 "ctrsv.f"
			    ix += *incx;
#line 321 "ctrsv.f"
/* L130: */
#line 321 "ctrsv.f"
			}
#line 322 "ctrsv.f"
			if (nounit) {
#line 322 "ctrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 322 "ctrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 322 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 322 "ctrsv.f"
			}
#line 323 "ctrsv.f"
		    }
#line 324 "ctrsv.f"
		    i__2 = jx;
#line 324 "ctrsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 325 "ctrsv.f"
		    jx += *incx;
#line 326 "ctrsv.f"
/* L140: */
#line 326 "ctrsv.f"
		}
#line 327 "ctrsv.f"
	    }
#line 328 "ctrsv.f"
	} else {
#line 329 "ctrsv.f"
	    if (*incx == 1) {
#line 330 "ctrsv.f"
		for (j = *n; j >= 1; --j) {
#line 331 "ctrsv.f"
		    i__1 = j;
#line 331 "ctrsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 332 "ctrsv.f"
		    if (noconj) {
#line 333 "ctrsv.f"
			i__1 = j + 1;
#line 333 "ctrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 334 "ctrsv.f"
			    i__2 = i__ + j * a_dim1;
#line 334 "ctrsv.f"
			    i__3 = i__;
#line 334 "ctrsv.f"
			    z__2.r = a[i__2].r * x[i__3].r - a[i__2].i * x[
				    i__3].i, z__2.i = a[i__2].r * x[i__3].i + 
				    a[i__2].i * x[i__3].r;
#line 334 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 334 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 335 "ctrsv.f"
/* L150: */
#line 335 "ctrsv.f"
			}
#line 336 "ctrsv.f"
			if (nounit) {
#line 336 "ctrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 336 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 336 "ctrsv.f"
			}
#line 337 "ctrsv.f"
		    } else {
#line 338 "ctrsv.f"
			i__1 = j + 1;
#line 338 "ctrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 339 "ctrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 339 "ctrsv.f"
			    i__2 = i__;
#line 339 "ctrsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 339 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 339 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 340 "ctrsv.f"
/* L160: */
#line 340 "ctrsv.f"
			}
#line 341 "ctrsv.f"
			if (nounit) {
#line 341 "ctrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 341 "ctrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 341 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 341 "ctrsv.f"
			}
#line 342 "ctrsv.f"
		    }
#line 343 "ctrsv.f"
		    i__1 = j;
#line 343 "ctrsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 344 "ctrsv.f"
/* L170: */
#line 344 "ctrsv.f"
		}
#line 345 "ctrsv.f"
	    } else {
#line 346 "ctrsv.f"
		kx += (*n - 1) * *incx;
#line 347 "ctrsv.f"
		jx = kx;
#line 348 "ctrsv.f"
		for (j = *n; j >= 1; --j) {
#line 349 "ctrsv.f"
		    ix = kx;
#line 350 "ctrsv.f"
		    i__1 = jx;
#line 350 "ctrsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 351 "ctrsv.f"
		    if (noconj) {
#line 352 "ctrsv.f"
			i__1 = j + 1;
#line 352 "ctrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 353 "ctrsv.f"
			    i__2 = i__ + j * a_dim1;
#line 353 "ctrsv.f"
			    i__3 = ix;
#line 353 "ctrsv.f"
			    z__2.r = a[i__2].r * x[i__3].r - a[i__2].i * x[
				    i__3].i, z__2.i = a[i__2].r * x[i__3].i + 
				    a[i__2].i * x[i__3].r;
#line 353 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 353 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ctrsv.f"
			    ix -= *incx;
#line 355 "ctrsv.f"
/* L180: */
#line 355 "ctrsv.f"
			}
#line 356 "ctrsv.f"
			if (nounit) {
#line 356 "ctrsv.f"
			    z_div(&z__1, &temp, &a[j + j * a_dim1]);
#line 356 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 356 "ctrsv.f"
			}
#line 357 "ctrsv.f"
		    } else {
#line 358 "ctrsv.f"
			i__1 = j + 1;
#line 358 "ctrsv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 359 "ctrsv.f"
			    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 359 "ctrsv.f"
			    i__2 = ix;
#line 359 "ctrsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 359 "ctrsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 359 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 360 "ctrsv.f"
			    ix -= *incx;
#line 361 "ctrsv.f"
/* L190: */
#line 361 "ctrsv.f"
			}
#line 362 "ctrsv.f"
			if (nounit) {
#line 362 "ctrsv.f"
			    d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 362 "ctrsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 362 "ctrsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 362 "ctrsv.f"
			}
#line 363 "ctrsv.f"
		    }
#line 364 "ctrsv.f"
		    i__1 = jx;
#line 364 "ctrsv.f"
		    x[i__1].r = temp.r, x[i__1].i = temp.i;
#line 365 "ctrsv.f"
		    jx -= *incx;
#line 366 "ctrsv.f"
/* L200: */
#line 366 "ctrsv.f"
		}
#line 367 "ctrsv.f"
	    }
#line 368 "ctrsv.f"
	}
#line 369 "ctrsv.f"
    }

#line 371 "ctrsv.f"
    return 0;

/*     End of CTRSV . */

} /* ctrsv_ */

