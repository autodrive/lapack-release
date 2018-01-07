#line 1 "ctbsv.f"
/* ctbsv.f -- translated by f2c (version 20100827).
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

#line 1 "ctbsv.f"
/* > \brief \b CTBSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,K,LDA,N */
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
/* > CTBSV  solves one of the systems of equations */
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
/* >          A is COMPLEX array, dimension ( LDA, N ) */
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
/* Subroutine */ int ctbsv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 229 "ctbsv.f"
    /* Parameter adjustments */
#line 229 "ctbsv.f"
    a_dim1 = *lda;
#line 229 "ctbsv.f"
    a_offset = 1 + a_dim1;
#line 229 "ctbsv.f"
    a -= a_offset;
#line 229 "ctbsv.f"
    --x;
#line 229 "ctbsv.f"

#line 229 "ctbsv.f"
    /* Function Body */
#line 229 "ctbsv.f"
    info = 0;
#line 230 "ctbsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 231 "ctbsv.f"
	info = 1;
#line 232 "ctbsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 234 "ctbsv.f"
	info = 2;
#line 235 "ctbsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 236 "ctbsv.f"
	info = 3;
#line 237 "ctbsv.f"
    } else if (*n < 0) {
#line 238 "ctbsv.f"
	info = 4;
#line 239 "ctbsv.f"
    } else if (*k < 0) {
#line 240 "ctbsv.f"
	info = 5;
#line 241 "ctbsv.f"
    } else if (*lda < *k + 1) {
#line 242 "ctbsv.f"
	info = 7;
#line 243 "ctbsv.f"
    } else if (*incx == 0) {
#line 244 "ctbsv.f"
	info = 9;
#line 245 "ctbsv.f"
    }
#line 246 "ctbsv.f"
    if (info != 0) {
#line 247 "ctbsv.f"
	xerbla_("CTBSV ", &info, (ftnlen)6);
#line 248 "ctbsv.f"
	return 0;
#line 249 "ctbsv.f"
    }

/*     Quick return if possible. */

#line 253 "ctbsv.f"
    if (*n == 0) {
#line 253 "ctbsv.f"
	return 0;
#line 253 "ctbsv.f"
    }

#line 255 "ctbsv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 256 "ctbsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 261 "ctbsv.f"
    if (*incx <= 0) {
#line 262 "ctbsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 263 "ctbsv.f"
    } else if (*incx != 1) {
#line 264 "ctbsv.f"
	kx = 1;
#line 265 "ctbsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed by sequentially with one pass through A. */

#line 270 "ctbsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 274 "ctbsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 275 "ctbsv.f"
	    kplus1 = *k + 1;
#line 276 "ctbsv.f"
	    if (*incx == 1) {
#line 277 "ctbsv.f"
		for (j = *n; j >= 1; --j) {
#line 278 "ctbsv.f"
		    i__1 = j;
#line 278 "ctbsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 279 "ctbsv.f"
			l = kplus1 - j;
#line 280 "ctbsv.f"
			if (nounit) {
#line 280 "ctbsv.f"
			    i__1 = j;
#line 280 "ctbsv.f"
			    z_div(&z__1, &x[j], &a[kplus1 + j * a_dim1]);
#line 280 "ctbsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 280 "ctbsv.f"
			}
#line 281 "ctbsv.f"
			i__1 = j;
#line 281 "ctbsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
/* Computing MAX */
#line 282 "ctbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 282 "ctbsv.f"
			i__1 = max(i__2,i__3);
#line 282 "ctbsv.f"
			for (i__ = j - 1; i__ >= i__1; --i__) {
#line 283 "ctbsv.f"
			    i__2 = i__;
#line 283 "ctbsv.f"
			    i__3 = i__;
#line 283 "ctbsv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 283 "ctbsv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 283 "ctbsv.f"
			    z__1.r = x[i__3].r - z__2.r, z__1.i = x[i__3].i - 
				    z__2.i;
#line 283 "ctbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 284 "ctbsv.f"
/* L10: */
#line 284 "ctbsv.f"
			}
#line 285 "ctbsv.f"
		    }
#line 286 "ctbsv.f"
/* L20: */
#line 286 "ctbsv.f"
		}
#line 287 "ctbsv.f"
	    } else {
#line 288 "ctbsv.f"
		kx += (*n - 1) * *incx;
#line 289 "ctbsv.f"
		jx = kx;
#line 290 "ctbsv.f"
		for (j = *n; j >= 1; --j) {
#line 291 "ctbsv.f"
		    kx -= *incx;
#line 292 "ctbsv.f"
		    i__1 = jx;
#line 292 "ctbsv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 293 "ctbsv.f"
			ix = kx;
#line 294 "ctbsv.f"
			l = kplus1 - j;
#line 295 "ctbsv.f"
			if (nounit) {
#line 295 "ctbsv.f"
			    i__1 = jx;
#line 295 "ctbsv.f"
			    z_div(&z__1, &x[jx], &a[kplus1 + j * a_dim1]);
#line 295 "ctbsv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 295 "ctbsv.f"
			}
#line 296 "ctbsv.f"
			i__1 = jx;
#line 296 "ctbsv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
/* Computing MAX */
#line 297 "ctbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 297 "ctbsv.f"
			i__1 = max(i__2,i__3);
#line 297 "ctbsv.f"
			for (i__ = j - 1; i__ >= i__1; --i__) {
#line 298 "ctbsv.f"
			    i__2 = ix;
#line 298 "ctbsv.f"
			    i__3 = ix;
#line 298 "ctbsv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 298 "ctbsv.f"
			    z__2.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__2.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 298 "ctbsv.f"
			    z__1.r = x[i__3].r - z__2.r, z__1.i = x[i__3].i - 
				    z__2.i;
#line 298 "ctbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 299 "ctbsv.f"
			    ix -= *incx;
#line 300 "ctbsv.f"
/* L30: */
#line 300 "ctbsv.f"
			}
#line 301 "ctbsv.f"
		    }
#line 302 "ctbsv.f"
		    jx -= *incx;
#line 303 "ctbsv.f"
/* L40: */
#line 303 "ctbsv.f"
		}
#line 304 "ctbsv.f"
	    }
#line 305 "ctbsv.f"
	} else {
#line 306 "ctbsv.f"
	    if (*incx == 1) {
#line 307 "ctbsv.f"
		i__1 = *n;
#line 307 "ctbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 308 "ctbsv.f"
		    i__2 = j;
#line 308 "ctbsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 309 "ctbsv.f"
			l = 1 - j;
#line 310 "ctbsv.f"
			if (nounit) {
#line 310 "ctbsv.f"
			    i__2 = j;
#line 310 "ctbsv.f"
			    z_div(&z__1, &x[j], &a[j * a_dim1 + 1]);
#line 310 "ctbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 310 "ctbsv.f"
			}
#line 311 "ctbsv.f"
			i__2 = j;
#line 311 "ctbsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
/* Computing MIN */
#line 312 "ctbsv.f"
			i__3 = *n, i__4 = j + *k;
#line 312 "ctbsv.f"
			i__2 = min(i__3,i__4);
#line 312 "ctbsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 313 "ctbsv.f"
			    i__3 = i__;
#line 313 "ctbsv.f"
			    i__4 = i__;
#line 313 "ctbsv.f"
			    i__5 = l + i__ + j * a_dim1;
#line 313 "ctbsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 313 "ctbsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 313 "ctbsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 314 "ctbsv.f"
/* L50: */
#line 314 "ctbsv.f"
			}
#line 315 "ctbsv.f"
		    }
#line 316 "ctbsv.f"
/* L60: */
#line 316 "ctbsv.f"
		}
#line 317 "ctbsv.f"
	    } else {
#line 318 "ctbsv.f"
		jx = kx;
#line 319 "ctbsv.f"
		i__1 = *n;
#line 319 "ctbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 320 "ctbsv.f"
		    kx += *incx;
#line 321 "ctbsv.f"
		    i__2 = jx;
#line 321 "ctbsv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 322 "ctbsv.f"
			ix = kx;
#line 323 "ctbsv.f"
			l = 1 - j;
#line 324 "ctbsv.f"
			if (nounit) {
#line 324 "ctbsv.f"
			    i__2 = jx;
#line 324 "ctbsv.f"
			    z_div(&z__1, &x[jx], &a[j * a_dim1 + 1]);
#line 324 "ctbsv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 324 "ctbsv.f"
			}
#line 325 "ctbsv.f"
			i__2 = jx;
#line 325 "ctbsv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
/* Computing MIN */
#line 326 "ctbsv.f"
			i__3 = *n, i__4 = j + *k;
#line 326 "ctbsv.f"
			i__2 = min(i__3,i__4);
#line 326 "ctbsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 327 "ctbsv.f"
			    i__3 = ix;
#line 327 "ctbsv.f"
			    i__4 = ix;
#line 327 "ctbsv.f"
			    i__5 = l + i__ + j * a_dim1;
#line 327 "ctbsv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 327 "ctbsv.f"
			    z__1.r = x[i__4].r - z__2.r, z__1.i = x[i__4].i - 
				    z__2.i;
#line 327 "ctbsv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 328 "ctbsv.f"
			    ix += *incx;
#line 329 "ctbsv.f"
/* L70: */
#line 329 "ctbsv.f"
			}
#line 330 "ctbsv.f"
		    }
#line 331 "ctbsv.f"
		    jx += *incx;
#line 332 "ctbsv.f"
/* L80: */
#line 332 "ctbsv.f"
		}
#line 333 "ctbsv.f"
	    }
#line 334 "ctbsv.f"
	}
#line 335 "ctbsv.f"
    } else {

/*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x. */

#line 339 "ctbsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 340 "ctbsv.f"
	    kplus1 = *k + 1;
#line 341 "ctbsv.f"
	    if (*incx == 1) {
#line 342 "ctbsv.f"
		i__1 = *n;
#line 342 "ctbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 343 "ctbsv.f"
		    i__2 = j;
#line 343 "ctbsv.f"
		    temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 344 "ctbsv.f"
		    l = kplus1 - j;
#line 345 "ctbsv.f"
		    if (noconj) {
/* Computing MAX */
#line 346 "ctbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 346 "ctbsv.f"
			i__4 = j - 1;
#line 346 "ctbsv.f"
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 347 "ctbsv.f"
			    i__2 = l + i__ + j * a_dim1;
#line 347 "ctbsv.f"
			    i__3 = i__;
#line 347 "ctbsv.f"
			    z__2.r = a[i__2].r * x[i__3].r - a[i__2].i * x[
				    i__3].i, z__2.i = a[i__2].r * x[i__3].i + 
				    a[i__2].i * x[i__3].r;
#line 347 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 347 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 348 "ctbsv.f"
/* L90: */
#line 348 "ctbsv.f"
			}
#line 349 "ctbsv.f"
			if (nounit) {
#line 349 "ctbsv.f"
			    z_div(&z__1, &temp, &a[kplus1 + j * a_dim1]);
#line 349 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 349 "ctbsv.f"
			}
#line 350 "ctbsv.f"
		    } else {
/* Computing MAX */
#line 351 "ctbsv.f"
			i__4 = 1, i__2 = j - *k;
#line 351 "ctbsv.f"
			i__3 = j - 1;
#line 351 "ctbsv.f"
			for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 352 "ctbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 352 "ctbsv.f"
			    i__4 = i__;
#line 352 "ctbsv.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 352 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 352 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 353 "ctbsv.f"
/* L100: */
#line 353 "ctbsv.f"
			}
#line 354 "ctbsv.f"
			if (nounit) {
#line 354 "ctbsv.f"
			    d_cnjg(&z__2, &a[kplus1 + j * a_dim1]);
#line 354 "ctbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 354 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ctbsv.f"
			}
#line 355 "ctbsv.f"
		    }
#line 356 "ctbsv.f"
		    i__3 = j;
#line 356 "ctbsv.f"
		    x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 357 "ctbsv.f"
/* L110: */
#line 357 "ctbsv.f"
		}
#line 358 "ctbsv.f"
	    } else {
#line 359 "ctbsv.f"
		jx = kx;
#line 360 "ctbsv.f"
		i__1 = *n;
#line 360 "ctbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 361 "ctbsv.f"
		    i__3 = jx;
#line 361 "ctbsv.f"
		    temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 362 "ctbsv.f"
		    ix = kx;
#line 363 "ctbsv.f"
		    l = kplus1 - j;
#line 364 "ctbsv.f"
		    if (noconj) {
/* Computing MAX */
#line 365 "ctbsv.f"
			i__3 = 1, i__4 = j - *k;
#line 365 "ctbsv.f"
			i__2 = j - 1;
#line 365 "ctbsv.f"
			for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
#line 366 "ctbsv.f"
			    i__3 = l + i__ + j * a_dim1;
#line 366 "ctbsv.f"
			    i__4 = ix;
#line 366 "ctbsv.f"
			    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[
				    i__4].i, z__2.i = a[i__3].r * x[i__4].i + 
				    a[i__3].i * x[i__4].r;
#line 366 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 366 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 367 "ctbsv.f"
			    ix += *incx;
#line 368 "ctbsv.f"
/* L120: */
#line 368 "ctbsv.f"
			}
#line 369 "ctbsv.f"
			if (nounit) {
#line 369 "ctbsv.f"
			    z_div(&z__1, &temp, &a[kplus1 + j * a_dim1]);
#line 369 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 369 "ctbsv.f"
			}
#line 370 "ctbsv.f"
		    } else {
/* Computing MAX */
#line 371 "ctbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 371 "ctbsv.f"
			i__4 = j - 1;
#line 371 "ctbsv.f"
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 372 "ctbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 372 "ctbsv.f"
			    i__2 = ix;
#line 372 "ctbsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 372 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 372 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 373 "ctbsv.f"
			    ix += *incx;
#line 374 "ctbsv.f"
/* L130: */
#line 374 "ctbsv.f"
			}
#line 375 "ctbsv.f"
			if (nounit) {
#line 375 "ctbsv.f"
			    d_cnjg(&z__2, &a[kplus1 + j * a_dim1]);
#line 375 "ctbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 375 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 375 "ctbsv.f"
			}
#line 376 "ctbsv.f"
		    }
#line 377 "ctbsv.f"
		    i__4 = jx;
#line 377 "ctbsv.f"
		    x[i__4].r = temp.r, x[i__4].i = temp.i;
#line 378 "ctbsv.f"
		    jx += *incx;
#line 379 "ctbsv.f"
		    if (j > *k) {
#line 379 "ctbsv.f"
			kx += *incx;
#line 379 "ctbsv.f"
		    }
#line 380 "ctbsv.f"
/* L140: */
#line 380 "ctbsv.f"
		}
#line 381 "ctbsv.f"
	    }
#line 382 "ctbsv.f"
	} else {
#line 383 "ctbsv.f"
	    if (*incx == 1) {
#line 384 "ctbsv.f"
		for (j = *n; j >= 1; --j) {
#line 385 "ctbsv.f"
		    i__1 = j;
#line 385 "ctbsv.f"
		    temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 386 "ctbsv.f"
		    l = 1 - j;
#line 387 "ctbsv.f"
		    if (noconj) {
/* Computing MIN */
#line 388 "ctbsv.f"
			i__1 = *n, i__4 = j + *k;
#line 388 "ctbsv.f"
			i__2 = j + 1;
#line 388 "ctbsv.f"
			for (i__ = min(i__1,i__4); i__ >= i__2; --i__) {
#line 389 "ctbsv.f"
			    i__1 = l + i__ + j * a_dim1;
#line 389 "ctbsv.f"
			    i__4 = i__;
#line 389 "ctbsv.f"
			    z__2.r = a[i__1].r * x[i__4].r - a[i__1].i * x[
				    i__4].i, z__2.i = a[i__1].r * x[i__4].i + 
				    a[i__1].i * x[i__4].r;
#line 389 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 389 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 390 "ctbsv.f"
/* L150: */
#line 390 "ctbsv.f"
			}
#line 391 "ctbsv.f"
			if (nounit) {
#line 391 "ctbsv.f"
			    z_div(&z__1, &temp, &a[j * a_dim1 + 1]);
#line 391 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 391 "ctbsv.f"
			}
#line 392 "ctbsv.f"
		    } else {
/* Computing MIN */
#line 393 "ctbsv.f"
			i__2 = *n, i__1 = j + *k;
#line 393 "ctbsv.f"
			i__4 = j + 1;
#line 393 "ctbsv.f"
			for (i__ = min(i__2,i__1); i__ >= i__4; --i__) {
#line 394 "ctbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 394 "ctbsv.f"
			    i__2 = i__;
#line 394 "ctbsv.f"
			    z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				    z__2.i = z__3.r * x[i__2].i + z__3.i * x[
				    i__2].r;
#line 394 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 394 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 395 "ctbsv.f"
/* L160: */
#line 395 "ctbsv.f"
			}
#line 396 "ctbsv.f"
			if (nounit) {
#line 396 "ctbsv.f"
			    d_cnjg(&z__2, &a[j * a_dim1 + 1]);
#line 396 "ctbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 396 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 396 "ctbsv.f"
			}
#line 397 "ctbsv.f"
		    }
#line 398 "ctbsv.f"
		    i__4 = j;
#line 398 "ctbsv.f"
		    x[i__4].r = temp.r, x[i__4].i = temp.i;
#line 399 "ctbsv.f"
/* L170: */
#line 399 "ctbsv.f"
		}
#line 400 "ctbsv.f"
	    } else {
#line 401 "ctbsv.f"
		kx += (*n - 1) * *incx;
#line 402 "ctbsv.f"
		jx = kx;
#line 403 "ctbsv.f"
		for (j = *n; j >= 1; --j) {
#line 404 "ctbsv.f"
		    i__4 = jx;
#line 404 "ctbsv.f"
		    temp.r = x[i__4].r, temp.i = x[i__4].i;
#line 405 "ctbsv.f"
		    ix = kx;
#line 406 "ctbsv.f"
		    l = 1 - j;
#line 407 "ctbsv.f"
		    if (noconj) {
/* Computing MIN */
#line 408 "ctbsv.f"
			i__4 = *n, i__2 = j + *k;
#line 408 "ctbsv.f"
			i__1 = j + 1;
#line 408 "ctbsv.f"
			for (i__ = min(i__4,i__2); i__ >= i__1; --i__) {
#line 409 "ctbsv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 409 "ctbsv.f"
			    i__2 = ix;
#line 409 "ctbsv.f"
			    z__2.r = a[i__4].r * x[i__2].r - a[i__4].i * x[
				    i__2].i, z__2.i = a[i__4].r * x[i__2].i + 
				    a[i__4].i * x[i__2].r;
#line 409 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 409 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 410 "ctbsv.f"
			    ix -= *incx;
#line 411 "ctbsv.f"
/* L180: */
#line 411 "ctbsv.f"
			}
#line 412 "ctbsv.f"
			if (nounit) {
#line 412 "ctbsv.f"
			    z_div(&z__1, &temp, &a[j * a_dim1 + 1]);
#line 412 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 412 "ctbsv.f"
			}
#line 413 "ctbsv.f"
		    } else {
/* Computing MIN */
#line 414 "ctbsv.f"
			i__1 = *n, i__4 = j + *k;
#line 414 "ctbsv.f"
			i__2 = j + 1;
#line 414 "ctbsv.f"
			for (i__ = min(i__1,i__4); i__ >= i__2; --i__) {
#line 415 "ctbsv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 415 "ctbsv.f"
			    i__1 = ix;
#line 415 "ctbsv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 415 "ctbsv.f"
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
#line 415 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 416 "ctbsv.f"
			    ix -= *incx;
#line 417 "ctbsv.f"
/* L190: */
#line 417 "ctbsv.f"
			}
#line 418 "ctbsv.f"
			if (nounit) {
#line 418 "ctbsv.f"
			    d_cnjg(&z__2, &a[j * a_dim1 + 1]);
#line 418 "ctbsv.f"
			    z_div(&z__1, &temp, &z__2);
#line 418 "ctbsv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 418 "ctbsv.f"
			}
#line 419 "ctbsv.f"
		    }
#line 420 "ctbsv.f"
		    i__2 = jx;
#line 420 "ctbsv.f"
		    x[i__2].r = temp.r, x[i__2].i = temp.i;
#line 421 "ctbsv.f"
		    jx -= *incx;
#line 422 "ctbsv.f"
		    if (*n - j >= *k) {
#line 422 "ctbsv.f"
			kx -= *incx;
#line 422 "ctbsv.f"
		    }
#line 423 "ctbsv.f"
/* L200: */
#line 423 "ctbsv.f"
		}
#line 424 "ctbsv.f"
	    }
#line 425 "ctbsv.f"
	}
#line 426 "ctbsv.f"
    }

#line 428 "ctbsv.f"
    return 0;

/*     End of CTBSV . */

} /* ctbsv_ */

