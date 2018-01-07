#line 1 "ztbmv.f"
/* ztbmv.f -- translated by f2c (version 20100827).
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

#line 1 "ztbmv.f"
/* > \brief \b ZTBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) */

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
/* > ZTBMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x,   or   x := A**H*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular band matrix, with ( k + 1 ) diagonals. */
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
/* >          A is COMPLEX*16 array, dimension ( LDA, N ). */
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
/* Subroutine */ int ztbmv_(char *uplo, char *trans, char *diag, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *x, integer 
	*incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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

#line 226 "ztbmv.f"
    /* Parameter adjustments */
#line 226 "ztbmv.f"
    a_dim1 = *lda;
#line 226 "ztbmv.f"
    a_offset = 1 + a_dim1;
#line 226 "ztbmv.f"
    a -= a_offset;
#line 226 "ztbmv.f"
    --x;
#line 226 "ztbmv.f"

#line 226 "ztbmv.f"
    /* Function Body */
#line 226 "ztbmv.f"
    info = 0;
#line 227 "ztbmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 228 "ztbmv.f"
	info = 1;
#line 229 "ztbmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 231 "ztbmv.f"
	info = 2;
#line 232 "ztbmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 233 "ztbmv.f"
	info = 3;
#line 234 "ztbmv.f"
    } else if (*n < 0) {
#line 235 "ztbmv.f"
	info = 4;
#line 236 "ztbmv.f"
    } else if (*k < 0) {
#line 237 "ztbmv.f"
	info = 5;
#line 238 "ztbmv.f"
    } else if (*lda < *k + 1) {
#line 239 "ztbmv.f"
	info = 7;
#line 240 "ztbmv.f"
    } else if (*incx == 0) {
#line 241 "ztbmv.f"
	info = 9;
#line 242 "ztbmv.f"
    }
#line 243 "ztbmv.f"
    if (info != 0) {
#line 244 "ztbmv.f"
	xerbla_("ZTBMV ", &info, (ftnlen)6);
#line 245 "ztbmv.f"
	return 0;
#line 246 "ztbmv.f"
    }

/*     Quick return if possible. */

#line 250 "ztbmv.f"
    if (*n == 0) {
#line 250 "ztbmv.f"
	return 0;
#line 250 "ztbmv.f"
    }

#line 252 "ztbmv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 253 "ztbmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX   too small for descending loops. */

#line 258 "ztbmv.f"
    if (*incx <= 0) {
#line 259 "ztbmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 260 "ztbmv.f"
    } else if (*incx != 1) {
#line 261 "ztbmv.f"
	kx = 1;
#line 262 "ztbmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 267 "ztbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*         Form  x := A*x. */

#line 271 "ztbmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 272 "ztbmv.f"
	    kplus1 = *k + 1;
#line 273 "ztbmv.f"
	    if (*incx == 1) {
#line 274 "ztbmv.f"
		i__1 = *n;
#line 274 "ztbmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 275 "ztbmv.f"
		    i__2 = j;
#line 275 "ztbmv.f"
		    if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 276 "ztbmv.f"
			i__2 = j;
#line 276 "ztbmv.f"
			temp.r = x[i__2].r, temp.i = x[i__2].i;
#line 277 "ztbmv.f"
			l = kplus1 - j;
/* Computing MAX */
#line 278 "ztbmv.f"
			i__2 = 1, i__3 = j - *k;
#line 278 "ztbmv.f"
			i__4 = j - 1;
#line 278 "ztbmv.f"
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 279 "ztbmv.f"
			    i__2 = i__;
#line 279 "ztbmv.f"
			    i__3 = i__;
#line 279 "ztbmv.f"
			    i__5 = l + i__ + j * a_dim1;
#line 279 "ztbmv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 279 "ztbmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 279 "ztbmv.f"
			    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 280 "ztbmv.f"
/* L10: */
#line 280 "ztbmv.f"
			}
#line 281 "ztbmv.f"
			if (nounit) {
#line 281 "ztbmv.f"
			    i__4 = j;
#line 281 "ztbmv.f"
			    i__2 = j;
#line 281 "ztbmv.f"
			    i__3 = kplus1 + j * a_dim1;
#line 281 "ztbmv.f"
			    z__1.r = x[i__2].r * a[i__3].r - x[i__2].i * a[
				    i__3].i, z__1.i = x[i__2].r * a[i__3].i + 
				    x[i__2].i * a[i__3].r;
#line 281 "ztbmv.f"
			    x[i__4].r = z__1.r, x[i__4].i = z__1.i;
#line 281 "ztbmv.f"
			}
#line 282 "ztbmv.f"
		    }
#line 283 "ztbmv.f"
/* L20: */
#line 283 "ztbmv.f"
		}
#line 284 "ztbmv.f"
	    } else {
#line 285 "ztbmv.f"
		jx = kx;
#line 286 "ztbmv.f"
		i__1 = *n;
#line 286 "ztbmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 287 "ztbmv.f"
		    i__4 = jx;
#line 287 "ztbmv.f"
		    if (x[i__4].r != 0. || x[i__4].i != 0.) {
#line 288 "ztbmv.f"
			i__4 = jx;
#line 288 "ztbmv.f"
			temp.r = x[i__4].r, temp.i = x[i__4].i;
#line 289 "ztbmv.f"
			ix = kx;
#line 290 "ztbmv.f"
			l = kplus1 - j;
/* Computing MAX */
#line 291 "ztbmv.f"
			i__4 = 1, i__2 = j - *k;
#line 291 "ztbmv.f"
			i__3 = j - 1;
#line 291 "ztbmv.f"
			for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 292 "ztbmv.f"
			    i__4 = ix;
#line 292 "ztbmv.f"
			    i__2 = ix;
#line 292 "ztbmv.f"
			    i__5 = l + i__ + j * a_dim1;
#line 292 "ztbmv.f"
			    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				    z__2.i = temp.r * a[i__5].i + temp.i * a[
				    i__5].r;
#line 292 "ztbmv.f"
			    z__1.r = x[i__2].r + z__2.r, z__1.i = x[i__2].i + 
				    z__2.i;
#line 292 "ztbmv.f"
			    x[i__4].r = z__1.r, x[i__4].i = z__1.i;
#line 293 "ztbmv.f"
			    ix += *incx;
#line 294 "ztbmv.f"
/* L30: */
#line 294 "ztbmv.f"
			}
#line 295 "ztbmv.f"
			if (nounit) {
#line 295 "ztbmv.f"
			    i__3 = jx;
#line 295 "ztbmv.f"
			    i__4 = jx;
#line 295 "ztbmv.f"
			    i__2 = kplus1 + j * a_dim1;
#line 295 "ztbmv.f"
			    z__1.r = x[i__4].r * a[i__2].r - x[i__4].i * a[
				    i__2].i, z__1.i = x[i__4].r * a[i__2].i + 
				    x[i__4].i * a[i__2].r;
#line 295 "ztbmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 295 "ztbmv.f"
			}
#line 296 "ztbmv.f"
		    }
#line 297 "ztbmv.f"
		    jx += *incx;
#line 298 "ztbmv.f"
		    if (j > *k) {
#line 298 "ztbmv.f"
			kx += *incx;
#line 298 "ztbmv.f"
		    }
#line 299 "ztbmv.f"
/* L40: */
#line 299 "ztbmv.f"
		}
#line 300 "ztbmv.f"
	    }
#line 301 "ztbmv.f"
	} else {
#line 302 "ztbmv.f"
	    if (*incx == 1) {
#line 303 "ztbmv.f"
		for (j = *n; j >= 1; --j) {
#line 304 "ztbmv.f"
		    i__1 = j;
#line 304 "ztbmv.f"
		    if (x[i__1].r != 0. || x[i__1].i != 0.) {
#line 305 "ztbmv.f"
			i__1 = j;
#line 305 "ztbmv.f"
			temp.r = x[i__1].r, temp.i = x[i__1].i;
#line 306 "ztbmv.f"
			l = 1 - j;
/* Computing MIN */
#line 307 "ztbmv.f"
			i__1 = *n, i__3 = j + *k;
#line 307 "ztbmv.f"
			i__4 = j + 1;
#line 307 "ztbmv.f"
			for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
#line 308 "ztbmv.f"
			    i__1 = i__;
#line 308 "ztbmv.f"
			    i__3 = i__;
#line 308 "ztbmv.f"
			    i__2 = l + i__ + j * a_dim1;
#line 308 "ztbmv.f"
			    z__2.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__2.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 308 "ztbmv.f"
			    z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + 
				    z__2.i;
#line 308 "ztbmv.f"
			    x[i__1].r = z__1.r, x[i__1].i = z__1.i;
#line 309 "ztbmv.f"
/* L50: */
#line 309 "ztbmv.f"
			}
#line 310 "ztbmv.f"
			if (nounit) {
#line 310 "ztbmv.f"
			    i__4 = j;
#line 310 "ztbmv.f"
			    i__1 = j;
#line 310 "ztbmv.f"
			    i__3 = j * a_dim1 + 1;
#line 310 "ztbmv.f"
			    z__1.r = x[i__1].r * a[i__3].r - x[i__1].i * a[
				    i__3].i, z__1.i = x[i__1].r * a[i__3].i + 
				    x[i__1].i * a[i__3].r;
#line 310 "ztbmv.f"
			    x[i__4].r = z__1.r, x[i__4].i = z__1.i;
#line 310 "ztbmv.f"
			}
#line 311 "ztbmv.f"
		    }
#line 312 "ztbmv.f"
/* L60: */
#line 312 "ztbmv.f"
		}
#line 313 "ztbmv.f"
	    } else {
#line 314 "ztbmv.f"
		kx += (*n - 1) * *incx;
#line 315 "ztbmv.f"
		jx = kx;
#line 316 "ztbmv.f"
		for (j = *n; j >= 1; --j) {
#line 317 "ztbmv.f"
		    i__4 = jx;
#line 317 "ztbmv.f"
		    if (x[i__4].r != 0. || x[i__4].i != 0.) {
#line 318 "ztbmv.f"
			i__4 = jx;
#line 318 "ztbmv.f"
			temp.r = x[i__4].r, temp.i = x[i__4].i;
#line 319 "ztbmv.f"
			ix = kx;
#line 320 "ztbmv.f"
			l = 1 - j;
/* Computing MIN */
#line 321 "ztbmv.f"
			i__4 = *n, i__1 = j + *k;
#line 321 "ztbmv.f"
			i__3 = j + 1;
#line 321 "ztbmv.f"
			for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
#line 322 "ztbmv.f"
			    i__4 = ix;
#line 322 "ztbmv.f"
			    i__1 = ix;
#line 322 "ztbmv.f"
			    i__2 = l + i__ + j * a_dim1;
#line 322 "ztbmv.f"
			    z__2.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__2.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 322 "ztbmv.f"
			    z__1.r = x[i__1].r + z__2.r, z__1.i = x[i__1].i + 
				    z__2.i;
#line 322 "ztbmv.f"
			    x[i__4].r = z__1.r, x[i__4].i = z__1.i;
#line 323 "ztbmv.f"
			    ix -= *incx;
#line 324 "ztbmv.f"
/* L70: */
#line 324 "ztbmv.f"
			}
#line 325 "ztbmv.f"
			if (nounit) {
#line 325 "ztbmv.f"
			    i__3 = jx;
#line 325 "ztbmv.f"
			    i__4 = jx;
#line 325 "ztbmv.f"
			    i__1 = j * a_dim1 + 1;
#line 325 "ztbmv.f"
			    z__1.r = x[i__4].r * a[i__1].r - x[i__4].i * a[
				    i__1].i, z__1.i = x[i__4].r * a[i__1].i + 
				    x[i__4].i * a[i__1].r;
#line 325 "ztbmv.f"
			    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 325 "ztbmv.f"
			}
#line 326 "ztbmv.f"
		    }
#line 327 "ztbmv.f"
		    jx -= *incx;
#line 328 "ztbmv.f"
		    if (*n - j >= *k) {
#line 328 "ztbmv.f"
			kx -= *incx;
#line 328 "ztbmv.f"
		    }
#line 329 "ztbmv.f"
/* L80: */
#line 329 "ztbmv.f"
		}
#line 330 "ztbmv.f"
	    }
#line 331 "ztbmv.f"
	}
#line 332 "ztbmv.f"
    } else {

/*        Form  x := A**T*x  or  x := A**H*x. */

#line 336 "ztbmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 337 "ztbmv.f"
	    kplus1 = *k + 1;
#line 338 "ztbmv.f"
	    if (*incx == 1) {
#line 339 "ztbmv.f"
		for (j = *n; j >= 1; --j) {
#line 340 "ztbmv.f"
		    i__3 = j;
#line 340 "ztbmv.f"
		    temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 341 "ztbmv.f"
		    l = kplus1 - j;
#line 342 "ztbmv.f"
		    if (noconj) {
#line 343 "ztbmv.f"
			if (nounit) {
#line 343 "ztbmv.f"
			    i__3 = kplus1 + j * a_dim1;
#line 343 "ztbmv.f"
			    z__1.r = temp.r * a[i__3].r - temp.i * a[i__3].i, 
				    z__1.i = temp.r * a[i__3].i + temp.i * a[
				    i__3].r;
#line 343 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 343 "ztbmv.f"
			}
/* Computing MAX */
#line 344 "ztbmv.f"
			i__4 = 1, i__1 = j - *k;
#line 344 "ztbmv.f"
			i__3 = max(i__4,i__1);
#line 344 "ztbmv.f"
			for (i__ = j - 1; i__ >= i__3; --i__) {
#line 345 "ztbmv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 345 "ztbmv.f"
			    i__1 = i__;
#line 345 "ztbmv.f"
			    z__2.r = a[i__4].r * x[i__1].r - a[i__4].i * x[
				    i__1].i, z__2.i = a[i__4].r * x[i__1].i + 
				    a[i__4].i * x[i__1].r;
#line 345 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 345 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 346 "ztbmv.f"
/* L90: */
#line 346 "ztbmv.f"
			}
#line 347 "ztbmv.f"
		    } else {
#line 348 "ztbmv.f"
			if (nounit) {
#line 348 "ztbmv.f"
			    d_cnjg(&z__2, &a[kplus1 + j * a_dim1]);
#line 348 "ztbmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 348 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 348 "ztbmv.f"
			}
/* Computing MAX */
#line 349 "ztbmv.f"
			i__4 = 1, i__1 = j - *k;
#line 349 "ztbmv.f"
			i__3 = max(i__4,i__1);
#line 349 "ztbmv.f"
			for (i__ = j - 1; i__ >= i__3; --i__) {
#line 350 "ztbmv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 350 "ztbmv.f"
			    i__4 = i__;
#line 350 "ztbmv.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 350 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 350 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 351 "ztbmv.f"
/* L100: */
#line 351 "ztbmv.f"
			}
#line 352 "ztbmv.f"
		    }
#line 353 "ztbmv.f"
		    i__3 = j;
#line 353 "ztbmv.f"
		    x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 354 "ztbmv.f"
/* L110: */
#line 354 "ztbmv.f"
		}
#line 355 "ztbmv.f"
	    } else {
#line 356 "ztbmv.f"
		kx += (*n - 1) * *incx;
#line 357 "ztbmv.f"
		jx = kx;
#line 358 "ztbmv.f"
		for (j = *n; j >= 1; --j) {
#line 359 "ztbmv.f"
		    i__3 = jx;
#line 359 "ztbmv.f"
		    temp.r = x[i__3].r, temp.i = x[i__3].i;
#line 360 "ztbmv.f"
		    kx -= *incx;
#line 361 "ztbmv.f"
		    ix = kx;
#line 362 "ztbmv.f"
		    l = kplus1 - j;
#line 363 "ztbmv.f"
		    if (noconj) {
#line 364 "ztbmv.f"
			if (nounit) {
#line 364 "ztbmv.f"
			    i__3 = kplus1 + j * a_dim1;
#line 364 "ztbmv.f"
			    z__1.r = temp.r * a[i__3].r - temp.i * a[i__3].i, 
				    z__1.i = temp.r * a[i__3].i + temp.i * a[
				    i__3].r;
#line 364 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 364 "ztbmv.f"
			}
/* Computing MAX */
#line 365 "ztbmv.f"
			i__4 = 1, i__1 = j - *k;
#line 365 "ztbmv.f"
			i__3 = max(i__4,i__1);
#line 365 "ztbmv.f"
			for (i__ = j - 1; i__ >= i__3; --i__) {
#line 366 "ztbmv.f"
			    i__4 = l + i__ + j * a_dim1;
#line 366 "ztbmv.f"
			    i__1 = ix;
#line 366 "ztbmv.f"
			    z__2.r = a[i__4].r * x[i__1].r - a[i__4].i * x[
				    i__1].i, z__2.i = a[i__4].r * x[i__1].i + 
				    a[i__4].i * x[i__1].r;
#line 366 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 366 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 367 "ztbmv.f"
			    ix -= *incx;
#line 368 "ztbmv.f"
/* L120: */
#line 368 "ztbmv.f"
			}
#line 369 "ztbmv.f"
		    } else {
#line 370 "ztbmv.f"
			if (nounit) {
#line 370 "ztbmv.f"
			    d_cnjg(&z__2, &a[kplus1 + j * a_dim1]);
#line 370 "ztbmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 370 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ztbmv.f"
			}
/* Computing MAX */
#line 371 "ztbmv.f"
			i__4 = 1, i__1 = j - *k;
#line 371 "ztbmv.f"
			i__3 = max(i__4,i__1);
#line 371 "ztbmv.f"
			for (i__ = j - 1; i__ >= i__3; --i__) {
#line 372 "ztbmv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 372 "ztbmv.f"
			    i__4 = ix;
#line 372 "ztbmv.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 372 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 372 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 373 "ztbmv.f"
			    ix -= *incx;
#line 374 "ztbmv.f"
/* L130: */
#line 374 "ztbmv.f"
			}
#line 375 "ztbmv.f"
		    }
#line 376 "ztbmv.f"
		    i__3 = jx;
#line 376 "ztbmv.f"
		    x[i__3].r = temp.r, x[i__3].i = temp.i;
#line 377 "ztbmv.f"
		    jx -= *incx;
#line 378 "ztbmv.f"
/* L140: */
#line 378 "ztbmv.f"
		}
#line 379 "ztbmv.f"
	    }
#line 380 "ztbmv.f"
	} else {
#line 381 "ztbmv.f"
	    if (*incx == 1) {
#line 382 "ztbmv.f"
		i__3 = *n;
#line 382 "ztbmv.f"
		for (j = 1; j <= i__3; ++j) {
#line 383 "ztbmv.f"
		    i__4 = j;
#line 383 "ztbmv.f"
		    temp.r = x[i__4].r, temp.i = x[i__4].i;
#line 384 "ztbmv.f"
		    l = 1 - j;
#line 385 "ztbmv.f"
		    if (noconj) {
#line 386 "ztbmv.f"
			if (nounit) {
#line 386 "ztbmv.f"
			    i__4 = j * a_dim1 + 1;
#line 386 "ztbmv.f"
			    z__1.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__1.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 386 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 386 "ztbmv.f"
			}
/* Computing MIN */
#line 387 "ztbmv.f"
			i__1 = *n, i__2 = j + *k;
#line 387 "ztbmv.f"
			i__4 = min(i__1,i__2);
#line 387 "ztbmv.f"
			for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 388 "ztbmv.f"
			    i__1 = l + i__ + j * a_dim1;
#line 388 "ztbmv.f"
			    i__2 = i__;
#line 388 "ztbmv.f"
			    z__2.r = a[i__1].r * x[i__2].r - a[i__1].i * x[
				    i__2].i, z__2.i = a[i__1].r * x[i__2].i + 
				    a[i__1].i * x[i__2].r;
#line 388 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 388 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 389 "ztbmv.f"
/* L150: */
#line 389 "ztbmv.f"
			}
#line 390 "ztbmv.f"
		    } else {
#line 391 "ztbmv.f"
			if (nounit) {
#line 391 "ztbmv.f"
			    d_cnjg(&z__2, &a[j * a_dim1 + 1]);
#line 391 "ztbmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 391 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 391 "ztbmv.f"
			}
/* Computing MIN */
#line 392 "ztbmv.f"
			i__1 = *n, i__2 = j + *k;
#line 392 "ztbmv.f"
			i__4 = min(i__1,i__2);
#line 392 "ztbmv.f"
			for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 393 "ztbmv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 393 "ztbmv.f"
			    i__1 = i__;
#line 393 "ztbmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 393 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 393 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 394 "ztbmv.f"
/* L160: */
#line 394 "ztbmv.f"
			}
#line 395 "ztbmv.f"
		    }
#line 396 "ztbmv.f"
		    i__4 = j;
#line 396 "ztbmv.f"
		    x[i__4].r = temp.r, x[i__4].i = temp.i;
#line 397 "ztbmv.f"
/* L170: */
#line 397 "ztbmv.f"
		}
#line 398 "ztbmv.f"
	    } else {
#line 399 "ztbmv.f"
		jx = kx;
#line 400 "ztbmv.f"
		i__3 = *n;
#line 400 "ztbmv.f"
		for (j = 1; j <= i__3; ++j) {
#line 401 "ztbmv.f"
		    i__4 = jx;
#line 401 "ztbmv.f"
		    temp.r = x[i__4].r, temp.i = x[i__4].i;
#line 402 "ztbmv.f"
		    kx += *incx;
#line 403 "ztbmv.f"
		    ix = kx;
#line 404 "ztbmv.f"
		    l = 1 - j;
#line 405 "ztbmv.f"
		    if (noconj) {
#line 406 "ztbmv.f"
			if (nounit) {
#line 406 "ztbmv.f"
			    i__4 = j * a_dim1 + 1;
#line 406 "ztbmv.f"
			    z__1.r = temp.r * a[i__4].r - temp.i * a[i__4].i, 
				    z__1.i = temp.r * a[i__4].i + temp.i * a[
				    i__4].r;
#line 406 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 406 "ztbmv.f"
			}
/* Computing MIN */
#line 407 "ztbmv.f"
			i__1 = *n, i__2 = j + *k;
#line 407 "ztbmv.f"
			i__4 = min(i__1,i__2);
#line 407 "ztbmv.f"
			for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 408 "ztbmv.f"
			    i__1 = l + i__ + j * a_dim1;
#line 408 "ztbmv.f"
			    i__2 = ix;
#line 408 "ztbmv.f"
			    z__2.r = a[i__1].r * x[i__2].r - a[i__1].i * x[
				    i__2].i, z__2.i = a[i__1].r * x[i__2].i + 
				    a[i__1].i * x[i__2].r;
#line 408 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 408 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 409 "ztbmv.f"
			    ix += *incx;
#line 410 "ztbmv.f"
/* L180: */
#line 410 "ztbmv.f"
			}
#line 411 "ztbmv.f"
		    } else {
#line 412 "ztbmv.f"
			if (nounit) {
#line 412 "ztbmv.f"
			    d_cnjg(&z__2, &a[j * a_dim1 + 1]);
#line 412 "ztbmv.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 412 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 412 "ztbmv.f"
			}
/* Computing MIN */
#line 413 "ztbmv.f"
			i__1 = *n, i__2 = j + *k;
#line 413 "ztbmv.f"
			i__4 = min(i__1,i__2);
#line 413 "ztbmv.f"
			for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 414 "ztbmv.f"
			    d_cnjg(&z__3, &a[l + i__ + j * a_dim1]);
#line 414 "ztbmv.f"
			    i__1 = ix;
#line 414 "ztbmv.f"
			    z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i, 
				    z__2.i = z__3.r * x[i__1].i + z__3.i * x[
				    i__1].r;
#line 414 "ztbmv.f"
			    z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
				    z__2.i;
#line 414 "ztbmv.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 415 "ztbmv.f"
			    ix += *incx;
#line 416 "ztbmv.f"
/* L190: */
#line 416 "ztbmv.f"
			}
#line 417 "ztbmv.f"
		    }
#line 418 "ztbmv.f"
		    i__4 = jx;
#line 418 "ztbmv.f"
		    x[i__4].r = temp.r, x[i__4].i = temp.i;
#line 419 "ztbmv.f"
		    jx += *incx;
#line 420 "ztbmv.f"
/* L200: */
#line 420 "ztbmv.f"
		}
#line 421 "ztbmv.f"
	    }
#line 422 "ztbmv.f"
	}
#line 423 "ztbmv.f"
    }

#line 425 "ztbmv.f"
    return 0;

/*     End of ZTBMV . */

} /* ztbmv_ */

