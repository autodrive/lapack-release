#line 1 "ztbsv.f"
/* ztbsv.f -- translated by f2c (version 20100827).
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

#line 1 "ztbsv.f"
/* > \brief \b ZTBSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,K,LDA,N */
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
/* > ZTBSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b,   or   A**H*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular band matrix, with ( k + 1 ) */
/* > diagonals. */
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
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with UPLO = 'U' or 'u', K specifies the number of */
/* >           super-diagonals of the matrix A. */
/* >           On entry with UPLO = 'L' or 'l', K specifies the number of */
/* >           sub-diagonals of the matrix A. */
/* >           K must satisfy  0 .le. K. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
/* >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the upper triangular */
/* >           band part of the matrix of coefficients, supplied column by */
/* >           column, with the leading diagonal of the matrix in row */
/* >           ( k + 1 ) of the array, the first super-diagonal starting at */
/* >           position 2 in row k, and so on. The top left k by k triangle */
/* >           of the array A is not referenced. */
/* >           The following program segment will transfer an upper */
/* >           triangular band matrix from conventional full matrix storage */
/* >           to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = K + 1 - J */
/* >                    DO 10, I = MAX( 1, J - K ), J */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the lower triangular */
/* >           band part of the matrix of coefficients, supplied column by */
/* >           column, with the leading diagonal of the matrix in row 1 of */
/* >           the array, the first sub-diagonal starting at position 1 in */
/* >           row 2, and so on. The bottom right k by k triangle of the */
/* >           array A is not referenced. */
/* >           The following program segment will transfer a lower */
/* >           triangular band matrix from conventional full matrix storage */
/* >           to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = 1 - J */
/* >                    DO 10, I = J, MIN( N, J + K ) */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Note that when DIAG = 'U' or 'u' the elements of the array A */
/* >           corresponding to the diagonal elements of the matrix are not */
/* >           referenced, but are assumed to be unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           ( k + 1 ). */
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
/* Subroutine */ int ztbsv_(char *uplo, char *trans, char *diag, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *x, integer 
	*incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, ix, jx, kx, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
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

#line 229 "ztbsv.f"
    /* Parameter adjustments */
#line 229 "ztbsv.f"
    a_dim1 = *lda;
#line 229 "ztbsv.f"
    a_offset = 1 + a_dim1;
#line 229 "ztbsv.f"
    a -= a_offset;
#line 229 "ztbsv.f"
    --x;
#line 229 "ztbsv.f"

#line 229 "ztbsv.f"
    /* Function Body */
#line 229 "ztbsv.f"
    info = 0;
#line 230 "ztbsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 231 "ztbsv.f"
	info = 1;
#line 232 "ztbsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 234 "ztbsv.f"
	info = 2;
#line 235 "ztbsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 236 "ztbsv.f"
	info = 3;
#line 237 "ztbsv.f"
    } else if (*n < 0) {
#line 238 "ztbsv.f"
	info = 4;
#line 239 "ztbsv.f"
    } else if (*k < 0) {
#line 240 "ztbsv.f"
	info = 5;
#line 241 "ztbsv.f"
    } else if (*lda < *k + 1) {
#line 242 "ztbsv.f"
	info = 7;
#line 243 "ztbsv.f"
    } else if (*incx == 0) {
#line 244 "ztbsv.f"
	info = 9;
#line 245 "ztbsv.f"
    }
#line 246 "ztbsv.f"
    if (info != 0) {
#line 247 "ztbsv.f"
	xerbla_("ZTBSV ", &info, (ftnlen)6);
#line 248 "ztbsv.f"
	return 0;
#line 249 "ztbsv.f"
    }

/*     Quick return if possible. */

#line 253 "ztbsv.f"
    if (*n == 0) {
#line 253 "ztbsv.f"
	return 0;
#line 253 "ztbsv.f"
    }

#line 255 "ztbsv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 256 "ztbsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 261 "ztbsv.f"
    if (*incx <= 0) {
#line 262 "ztbsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 263 "ztbsv.f"
    } else if (*incx != 1) {
#line 264 "ztbsv.f"
	kx = 1;
#line 265 "ztbsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed by sequentially with one pass through A. */

#line 270 "ztbsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 274 "ztbsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 275 "ztbsv.f"
	    kplus1 = *k + 1;
#line 276 "ztbsv.f"
	    if (*incx == 1) {
#line 277 "ztbsv.f"
		for (j = *n; j >= 1; --j) {
#line 278 "ztbsv.f"
		    i__1 = j;
#line 278 "ztbsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 279 "ztbsv.f"
			l = kplus1 - j;
#line 280 "ztbsv.f"
			if (nounit) {
#line 280 "ztbsv.f"
			    i__1 = j;
#line 280 "ztbsv.f"
			    z_div(&z__1, &x[j], &a[kplus1 + j * a_dim1]);
#line 280 "ztbsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 280 "ztbsv.f"
			}
#line 281 "ztbsv.f"
			i__1 = j;
#line 281 "ztbsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
/* Computing MAX */
#line 282 "ztbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 282 "ztbsv.f"
			i__1 = max(i__2,i__3);
#line 282 "ztbsv.f"
			for (i__ = j - 1; i__ >= i__1; --i__) {
#line 283 "ztbsv.f"
			    i__2 = i__;
#line 283 "ztbsv.f"
			    i__3 = i__;
#line 283 "ztbsv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 283 "ztbsv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 283 "ztbsv.f"
			    z__1.r = x[i__3].r - z__2.r, z__1.i = x[i__3].i - 
				    z__2.i;
#line 283 "ztbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 284 "ztbsv.f"
/* L10: */
#line 284 "ztbsv.f"
			}
#line 285 "ztbsv.f"
		    }
#line 286 "ztbsv.f"
/* L20: */
#line 286 "ztbsv.f"
		}
#line 287 "ztbsv.f"
	    } else {
#line 288 "ztbsv.f"
		kx += (*n - 1) * *incx;
#line 289 "ztbsv.f"
		jx = kx;
#line 290 "ztbsv.f"
		for (j = *n; j >= 1; --j) {
#line 291 "ztbsv.f"
		    kx -= *incx;
#line 292 "ztbsv.f"
		    i__1 = jx;
#line 292 "ztbsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 293 "ztbsv.f"
			ix = kx;
#line 294 "ztbsv.f"
			l = kplus1 - j;
#line 295 "ztbsv.f"
			if (nounit) {
#line 295 "ztbsv.f"
			    i__1 = jx;
#line 295 "ztbsv.f"
			    z_div(&z__1, &x[jx], &a[kplus1 + j * a_dim1]);
#line 295 "ztbsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 295 "ztbsv.f"
			}
#line 296 "ztbsv.f"
			i__1 = jx;
#line 296 "ztbsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
/* Computing MAX */
#line 297 "ztbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 297 "ztbsv.f"
			i__1 = max(i__2,i__3);
#line 297 "ztbsv.f"
			for (i__ = j - 1; i__ >= i__1; --i__) {
#line 298 "ztbsv.f"
			    i__2 = ix;
#line 298 "ztbsv.f"
			    i__3 = ix;
#line 298 "ztbsv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 298 "ztbsv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 298 "ztbsv.f"
			    z__1.r = x[i__3].r - z__2.r, z__1.i = x[i__3].i - 
				    z__2.i;
#line 298 "ztbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 299 "ztbsv.f"
			    ix -= *incx;
#line 300 "ztbsv.f"
/* L30: */
#line 300 "ztbsv.f"
			}
#line 301 "ztbsv.f"
		    }
#line 302 "ztbsv.f"
		    jx -= *incx;
#line 303 "ztbsv.f"
/* L40: */
#line 303 "ztbsv.f"
		}
#line 304 "ztbsv.f"
	    }
#line 305 "ztbsv.f"
	} else {
#line 306 "ztbsv.f"
	    if (*incx == 1) {
#line 307 "ztbsv.f"
		i__1 = *n;
#line 307 "ztbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 308 "ztbsv.f"
		    i__2 = j;
#line 308 "ztbsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 309 "ztbsv.f"
			l = 1 - j;
#line 310 "ztbsv.f"
			if (nounit) {
#line 310 "ztbsv.f"
			    i__2 = j;
#line 310 "ztbsv.f"
			    z_div(&z__1, &x[j], &a[j * a_dim1 + 1]);
#line 310 "ztbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 310 "ztbsv.f"
			}
#line 311 "ztbsv.f"
			i__2 = j;
#line 311 "ztbsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
/* Computing MIN */
#line 312 "ztbsv.f"
			i__3 = *n, i__4 = j + *k;
#line 312 "ztbsv.f"
			i__2 = min(i__3,i__4);
#line 312 "ztbsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 313 "ztbsv.f"
			    i__3 = i__;
#line 313 "ztbsv.f"
			    i__4 = i__;
#line 313 "ztbsv.f"
			    i__5 = l + i__ + j * a_dim1;
#line 313 "ztbsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 313 "ztbsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 313 "ztbsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 314 "ztbsv.f"
/* L50: */
#line 314 "ztbsv.f"
			}
#line 315 "ztbsv.f"
		    }
#line 316 "ztbsv.f"
/* L60: */
#line 316 "ztbsv.f"
		}
#line 317 "ztbsv.f"
	    } else {
#line 318 "ztbsv.f"
		jx = kx;
#line 319 "ztbsv.f"
		i__1 = *n;
#line 319 "ztbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 320 "ztbsv.f"
		    kx += *incx;
#line 321 "ztbsv.f"
		    i__2 = jx;
#line 321 "ztbsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 322 "ztbsv.f"
			ix = kx;
#line 323 "ztbsv.f"
			l = 1 - j;
#line 324 "ztbsv.f"
			if (nounit) {
#line 324 "ztbsv.f"
			    i__2 = jx;
#line 324 "ztbsv.f"
			    z_div(&z__1, &x[jx], &a[j * a_dim1 + 1]);
#line 324 "ztbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 324 "ztbsv.f"
			}
#line 325 "ztbsv.f"
			i__2 = jx;
#line 325 "ztbsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
/* Computing MIN */
#line 326 "ztbsv.f"
			i__3 = *n, i__4 = j + *k;
#line 326 "ztbsv.f"
			i__2 = min(i__3,i__4);
#line 326 "ztbsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 327 "ztbsv.f"
			    i__3 = ix;
#line 327 "ztbsv.f"
			    i__4 = ix;
#line 327 "ztbsv.f"
			    i__5 = l + i__ + j * a_dim1;
#line 327 "ztbsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 327 "ztbsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 327 "ztbsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 328 "ztbsv.f"
			    ix += *incx;
#line 329 "ztbsv.f"
/* L70: */
#line 329 "ztbsv.f"
			}
#line 330 "ztbsv.f"
		    }
#line 331 "ztbsv.f"
		    jx += *incx;
#line 332 "ztbsv.f"
/* L80: */
#line 332 "ztbsv.f"
		}
#line 333 "ztbsv.f"
	    }
#line 334 "ztbsv.f"
	}
#line 335 "ztbsv.f"
    } else {

/*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x. */

#line 339 "ztbsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 340 "ztbsv.f"
	    kplus1 = *k + 1;
#line 341 "ztbsv.f"
	    if (*incx == 1) {
#line 342 "ztbsv.f"
		i__1 = *n;
#line 342 "ztbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 343 "ztbsv.f"
		    i__2 = j;
#line 343 "ztbsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 344 "ztbsv.f"
		    l = kplus1 - j;
#line 345 "ztbsv.f"
		    if (noconj) {
/* Computing MAX */
#line 346 "ztbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 346 "ztbsv.f"
			i__4 = j - 1;
#line 346 "ztbsv.f"
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 347 "ztbsv.f"
			    i__2 = l + i__ + j * a_dim1;
#line 347 "ztbsv.f"
			    i__3 = i__;
#line 347 "ztbsv.f"
			    z__2.r = a[i__2].r * x[i__3].r - a[i__2].i * x[
				    i__3].i, z__2.i = a[i__2].r * x[i__3].i + 
				    a[i__2].i * x[i__3].r;
#line 347 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 347 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 348 "ztbsv.f"
/* L90: */
#line 348 "ztbsv.f"
			}
#line 349 "ztbsv.f"
			if (nounit) {
#line 349 "ztbsv.f"
			    z_div(&z__1, &temp, &a[kplus1 + j * a_dim1]);
#line 349 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 349 "ztbsv.f"
			}
#line 350 "ztbsv.f"
		    } else {
/* Computing MAX */
#line 351 "ztbsv.f"
			i__4 = 1, i__2 = j - *k;
#line 351 "ztbsv.f"
			i__3 = j - 1;
#line 351 "ztbsv.f"
			for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 352 "ztbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 352 "ztbsv.f"
			    i__4 = i__;
#line 352 "ztbsv.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 352 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 352 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 353 "ztbsv.f"
/* L100: */
#line 353 "ztbsv.f"
			}
#line 354 "ztbsv.f"
			if (nounit) {
#line 354 "ztbsv.f"
			    d_cnjg(&z__2, &a[kplus1 + j * a_dim1]);
#line 354 "ztbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 354 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ztbsv.f"
			}
#line 355 "ztbsv.f"
		    }
#line 356 "ztbsv.f"
		    i__3 = j;
#line 356 "ztbsv.f"
		    x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 357 "ztbsv.f"
/* L110: */
#line 357 "ztbsv.f"
		}
#line 358 "ztbsv.f"
	    } else {
#line 359 "ztbsv.f"
		jx = kx;
#line 360 "ztbsv.f"
		i__1 = *n;
#line 360 "ztbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 361 "ztbsv.f"
		    i__3 = jx;
#line 361 "ztbsv.f"
		    temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 362 "ztbsv.f"
		    ix = kx;
#line 363 "ztbsv.f"
		    l = kplus1 - j;
#line 364 "ztbsv.f"
		    if (noconj) {
/* Computing MAX */
#line 365 "ztbsv.f"
			i__3 = 1, i__4 = j - *k;
#line 365 "ztbsv.f"
			i__2 = j - 1;
#line 365 "ztbsv.f"
			for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 366 "ztbsv.f"
			    i__3 = l + i__ + j * a_dim1;
#line 366 "ztbsv.f"
			    i__4 = ix;
#line 366 "ztbsv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 366 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 366 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 367 "ztbsv.f"
			    ix += *incx;
#line 368 "ztbsv.f"
/* L120: */
#line 368 "ztbsv.f"
			}
#line 369 "ztbsv.f"
			if (nounit) {
#line 369 "ztbsv.f"
			    z_div(&z__1, &temp, &a[kplus1 + j * a_dim1]);
#line 369 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 369 "ztbsv.f"
			}
#line 370 "ztbsv.f"
		    } else {
/* Computing MAX */
#line 371 "ztbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 371 "ztbsv.f"
			i__4 = j - 1;
#line 371 "ztbsv.f"
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 372 "ztbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 372 "ztbsv.f"
			    i__2 = ix;
#line 372 "ztbsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 372 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 372 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 373 "ztbsv.f"
			    ix += *incx;
#line 374 "ztbsv.f"
/* L130: */
#line 374 "ztbsv.f"
			}
#line 375 "ztbsv.f"
			if (nounit) {
#line 375 "ztbsv.f"
			    d_cnjg(&z__2, &a[kplus1 + j * a_dim1]);
#line 375 "ztbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 375 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 375 "ztbsv.f"
			}
#line 376 "ztbsv.f"
		    }
#line 377 "ztbsv.f"
		    i__4 = jx;
#line 377 "ztbsv.f"
		    x[i__4].r = temp.r, x[i__4].i = temp.i;
#line 378 "ztbsv.f"
		    jx += *incx;
#line 379 "ztbsv.f"
		    if (j > *k) {
#line 379 "ztbsv.f"
			kx += *incx;
#line 379 "ztbsv.f"
		    }
#line 380 "ztbsv.f"
/* L140: */
#line 380 "ztbsv.f"
		}
#line 381 "ztbsv.f"
	    }
#line 382 "ztbsv.f"
	} else {
#line 383 "ztbsv.f"
	    if (*incx == 1) {
#line 384 "ztbsv.f"
		for (j = *n; j >= 1; --j) {
#line 385 "ztbsv.f"
		    i__1 = j;
#line 385 "ztbsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 386 "ztbsv.f"
		    l = 1 - j;
#line 387 "ztbsv.f"
		    if (noconj) {
/* Computing MIN */
#line 388 "ztbsv.f"
			i__1 = *n, i__4 = j + *k;
#line 388 "ztbsv.f"
			i__2 = j + 1;
#line 388 "ztbsv.f"
			for (i__ = min(i__1,i__4); i__ >= i__2; --i__) {
#line 389 "ztbsv.f"
			    i__1 = l + i__ + j * a_dim1;
#line 389 "ztbsv.f"
			    i__4 = i__;
#line 389 "ztbsv.f"
			    z__2.r = a[i__1].r * x[i__4].r - a[i__1].i * x[
				    i__4].i, z__2.i = a[i__1].r * x[i__4].i + 
				    a[i__1].i * x[i__4].r;
#line 389 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 389 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 390 "ztbsv.f"
/* L150: */
#line 390 "ztbsv.f"
			}
#line 391 "ztbsv.f"
			if (nounit) {
#line 391 "ztbsv.f"
			    z_div(&z__1, &temp, &a[j * a_dim1 + 1]);
#line 391 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 391 "ztbsv.f"
			}
#line 392 "ztbsv.f"
		    } else {
/* Computing MIN */
#line 393 "ztbsv.f"
			i__2 = *n, i__1 = j + *k;
#line 393 "ztbsv.f"
			i__4 = j + 1;
#line 393 "ztbsv.f"
			for (i__ = min(i__2,i__1); i__ >= i__4; --i__) {
#line 394 "ztbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 394 "ztbsv.f"
			    i__2 = i__;
#line 394 "ztbsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 394 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 394 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 395 "ztbsv.f"
/* L160: */
#line 395 "ztbsv.f"
			}
#line 396 "ztbsv.f"
			if (nounit) {
#line 396 "ztbsv.f"
			    d_cnjg(&z__2, &a[j * a_dim1 + 1]);
#line 396 "ztbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 396 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 396 "ztbsv.f"
			}
#line 397 "ztbsv.f"
		    }
#line 398 "ztbsv.f"
		    i__4 = j;
#line 398 "ztbsv.f"
		    x[i__4].r = temp.r, x[i__4].i = temp.i;
#line 399 "ztbsv.f"
/* L170: */
#line 399 "ztbsv.f"
		}
#line 400 "ztbsv.f"
	    } else {
#line 401 "ztbsv.f"
		kx += (*n - 1) * *incx;
#line 402 "ztbsv.f"
		jx = kx;
#line 403 "ztbsv.f"
		for (j = *n; j >= 1; --j) {
#line 404 "ztbsv.f"
		    i__4 = jx;
#line 404 "ztbsv.f"
		    temp.r = x[i__4].r, temp.i = x[i__4].i;
#line 405 "ztbsv.f"
		    ix = kx;
#line 406 "ztbsv.f"
		    l = 1 - j;
#line 407 "ztbsv.f"
		    if (noconj) {
/* Computing MIN */
#line 408 "ztbsv.f"
			i__4 = *n, i__2 = j + *k;
#line 408 "ztbsv.f"
			i__1 = j + 1;
#line 408 "ztbsv.f"
			for (i__ = min(i__4,i__2); i__ >= i__1; --i__) {
#line 409 "ztbsv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 409 "ztbsv.f"
			    i__2 = ix;
#line 409 "ztbsv.f"
			    z__2.r = a[i__4].r * x[i__2].r - a[i__4].i * x[
				    i__2].i, z__2.i = a[i__4].r * x[i__2].i + 
				    a[i__4].i * x[i__2].r;
#line 409 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 409 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 410 "ztbsv.f"
			    ix -= *incx;
#line 411 "ztbsv.f"
/* L180: */
#line 411 "ztbsv.f"
			}
#line 412 "ztbsv.f"
			if (nounit) {
#line 412 "ztbsv.f"
			    z_div(&z__1, &temp, &a[j * a_dim1 + 1]);
#line 412 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 412 "ztbsv.f"
			}
#line 413 "ztbsv.f"
		    } else {
/* Computing MIN */
#line 414 "ztbsv.f"
			i__1 = *n, i__4 = j + *k;
#line 414 "ztbsv.f"
			i__2 = j + 1;
#line 414 "ztbsv.f"
			for (i__ = min(i__1,i__4); i__ >= i__2; --i__) {
#line 415 "ztbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 415 "ztbsv.f"
			    i__1 = ix;
#line 415 "ztbsv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 415 "ztbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 415 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 416 "ztbsv.f"
			    ix -= *incx;
#line 417 "ztbsv.f"
/* L190: */
#line 417 "ztbsv.f"
			}
#line 418 "ztbsv.f"
			if (nounit) {
#line 418 "ztbsv.f"
			    d_cnjg(&z__2, &a[j * a_dim1 + 1]);
#line 418 "ztbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 418 "ztbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 418 "ztbsv.f"
			}
#line 419 "ztbsv.f"
		    }
#line 420 "ztbsv.f"
		    i__2 = jx;
#line 420 "ztbsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 421 "ztbsv.f"
		    jx -= *incx;
#line 422 "ztbsv.f"
		    if (*n - j >= *k) {
#line 422 "ztbsv.f"
			kx -= *incx;
#line 422 "ztbsv.f"
		    }
#line 423 "ztbsv.f"
/* L200: */
#line 423 "ztbsv.f"
		}
#line 424 "ztbsv.f"
	    }
#line 425 "ztbsv.f"
	}
#line 426 "ztbsv.f"
    }

#line 428 "ztbsv.f"
    return 0;

/*     End of ZTBSV . */

} /* ztbsv_ */

