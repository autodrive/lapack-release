#line 1 "ctrsm.f"
/* ctrsm.f -- translated by f2c (version 20100827).
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

#line 1 "ctrsm.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b CTRSM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA */
/*       INTEGER LDA,LDB,M,N */
/*       CHARACTER DIAG,SIDE,TRANSA,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),B(LDB,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRSM  solves one of the matrix equations */
/* > */
/* >    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, */
/* > */
/* > where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H. */
/* > */
/* > The matrix X is overwritten on B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry, SIDE specifies whether op( A ) appears on the left */
/* >           or right of X as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */
/* > */
/* >              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix A is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n'   op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't'   op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c'   op( A ) = A**H. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit triangular */
/* >           as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of B. M must be at */
/* >           least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of B.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, k ), */
/* >           where k is m when SIDE = 'L' or 'l' */
/* >             and k is n when SIDE = 'R' or 'r'. */
/* >           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
/* >           upper triangular part of the array  A must contain the upper */
/* >           triangular matrix  and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
/* >           lower triangular part of the array  A must contain the lower */
/* >           triangular matrix  and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
/* >           A  are not referenced either,  but are assumed to be  unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/* >           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' */
/* >           then LDA must be at least max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array of DIMENSION ( LDB, n ). */
/* >           Before entry,  the leading  m by n part of the array  B must */
/* >           contain  the  right-hand  side  matrix  B,  and  on exit  is */
/* >           overwritten by the solution matrix  X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           max( 1, m ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex_blas_level3 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 3 Blas routine. */
/* > */
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctrsm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, ftnlen side_len, ftnlen 
	uplo_len, ftnlen transa_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lside;
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconj, nounit;


/*  -- Reference BLAS level3 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

#line 223 "ctrsm.f"
    /* Parameter adjustments */
#line 223 "ctrsm.f"
    a_dim1 = *lda;
#line 223 "ctrsm.f"
    a_offset = 1 + a_dim1;
#line 223 "ctrsm.f"
    a -= a_offset;
#line 223 "ctrsm.f"
    b_dim1 = *ldb;
#line 223 "ctrsm.f"
    b_offset = 1 + b_dim1;
#line 223 "ctrsm.f"
    b -= b_offset;
#line 223 "ctrsm.f"

#line 223 "ctrsm.f"
    /* Function Body */
#line 223 "ctrsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 224 "ctrsm.f"
    if (lside) {
#line 225 "ctrsm.f"
	nrowa = *m;
#line 226 "ctrsm.f"
    } else {
#line 227 "ctrsm.f"
	nrowa = *n;
#line 228 "ctrsm.f"
    }
#line 229 "ctrsm.f"
    noconj = lsame_(transa, "T", (ftnlen)1, (ftnlen)1);
#line 230 "ctrsm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 231 "ctrsm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 233 "ctrsm.f"
    info = 0;
#line 234 "ctrsm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 235 "ctrsm.f"
	info = 1;
#line 236 "ctrsm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 237 "ctrsm.f"
	info = 2;
#line 238 "ctrsm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 241 "ctrsm.f"
	info = 3;
#line 242 "ctrsm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 243 "ctrsm.f"
	info = 4;
#line 244 "ctrsm.f"
    } else if (*m < 0) {
#line 245 "ctrsm.f"
	info = 5;
#line 246 "ctrsm.f"
    } else if (*n < 0) {
#line 247 "ctrsm.f"
	info = 6;
#line 248 "ctrsm.f"
    } else if (*lda < max(1,nrowa)) {
#line 249 "ctrsm.f"
	info = 9;
#line 250 "ctrsm.f"
    } else if (*ldb < max(1,*m)) {
#line 251 "ctrsm.f"
	info = 11;
#line 252 "ctrsm.f"
    }
#line 253 "ctrsm.f"
    if (info != 0) {
#line 254 "ctrsm.f"
	xerbla_("CTRSM ", &info, (ftnlen)6);
#line 255 "ctrsm.f"
	return 0;
#line 256 "ctrsm.f"
    }

/*     Quick return if possible. */

#line 260 "ctrsm.f"
    if (*m == 0 || *n == 0) {
#line 260 "ctrsm.f"
	return 0;
#line 260 "ctrsm.f"
    }

/*     And when  alpha.eq.zero. */

#line 264 "ctrsm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 265 "ctrsm.f"
	i__1 = *n;
#line 265 "ctrsm.f"
	for (j = 1; j <= i__1; ++j) {
#line 266 "ctrsm.f"
	    i__2 = *m;
#line 266 "ctrsm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 267 "ctrsm.f"
		i__3 = i__ + j * b_dim1;
#line 267 "ctrsm.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 268 "ctrsm.f"
/* L10: */
#line 268 "ctrsm.f"
	    }
#line 269 "ctrsm.f"
/* L20: */
#line 269 "ctrsm.f"
	}
#line 270 "ctrsm.f"
	return 0;
#line 271 "ctrsm.f"
    }

/*     Start the operations. */

#line 275 "ctrsm.f"
    if (lside) {
#line 276 "ctrsm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*inv( A )*B. */

#line 280 "ctrsm.f"
	    if (upper) {
#line 281 "ctrsm.f"
		i__1 = *n;
#line 281 "ctrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 282 "ctrsm.f"
		    if (alpha->r != 1. || alpha->i != 0.) {
#line 283 "ctrsm.f"
			i__2 = *m;
#line 283 "ctrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 284 "ctrsm.f"
			    i__3 = i__ + j * b_dim1;
#line 284 "ctrsm.f"
			    i__4 = i__ + j * b_dim1;
#line 284 "ctrsm.f"
			    z__1.r = alpha->r * b[i__4].r - alpha->i * b[i__4]
				    .i, z__1.i = alpha->r * b[i__4].i + 
				    alpha->i * b[i__4].r;
#line 284 "ctrsm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 285 "ctrsm.f"
/* L30: */
#line 285 "ctrsm.f"
			}
#line 286 "ctrsm.f"
		    }
#line 287 "ctrsm.f"
		    for (k = *m; k >= 1; --k) {
#line 288 "ctrsm.f"
			i__2 = k + j * b_dim1;
#line 288 "ctrsm.f"
			if (b[i__2].r != 0. || b[i__2].i != 0.) {
#line 289 "ctrsm.f"
			    if (nounit) {
#line 289 "ctrsm.f"
				i__2 = k + j * b_dim1;
#line 289 "ctrsm.f"
				z_div(&z__1, &b[k + j * b_dim1], &a[k + k * 
					a_dim1]);
#line 289 "ctrsm.f"
				b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 289 "ctrsm.f"
			    }
#line 290 "ctrsm.f"
			    i__2 = k - 1;
#line 290 "ctrsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 291 "ctrsm.f"
				i__3 = i__ + j * b_dim1;
#line 291 "ctrsm.f"
				i__4 = i__ + j * b_dim1;
#line 291 "ctrsm.f"
				i__5 = k + j * b_dim1;
#line 291 "ctrsm.f"
				i__6 = i__ + k * a_dim1;
#line 291 "ctrsm.f"
				z__2.r = b[i__5].r * a[i__6].r - b[i__5].i * 
					a[i__6].i, z__2.i = b[i__5].r * a[
					i__6].i + b[i__5].i * a[i__6].r;
#line 291 "ctrsm.f"
				z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4]
					.i - z__2.i;
#line 291 "ctrsm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 292 "ctrsm.f"
/* L40: */
#line 292 "ctrsm.f"
			    }
#line 293 "ctrsm.f"
			}
#line 294 "ctrsm.f"
/* L50: */
#line 294 "ctrsm.f"
		    }
#line 295 "ctrsm.f"
/* L60: */
#line 295 "ctrsm.f"
		}
#line 296 "ctrsm.f"
	    } else {
#line 297 "ctrsm.f"
		i__1 = *n;
#line 297 "ctrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 298 "ctrsm.f"
		    if (alpha->r != 1. || alpha->i != 0.) {
#line 299 "ctrsm.f"
			i__2 = *m;
#line 299 "ctrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 300 "ctrsm.f"
			    i__3 = i__ + j * b_dim1;
#line 300 "ctrsm.f"
			    i__4 = i__ + j * b_dim1;
#line 300 "ctrsm.f"
			    z__1.r = alpha->r * b[i__4].r - alpha->i * b[i__4]
				    .i, z__1.i = alpha->r * b[i__4].i + 
				    alpha->i * b[i__4].r;
#line 300 "ctrsm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 301 "ctrsm.f"
/* L70: */
#line 301 "ctrsm.f"
			}
#line 302 "ctrsm.f"
		    }
#line 303 "ctrsm.f"
		    i__2 = *m;
#line 303 "ctrsm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 304 "ctrsm.f"
			i__3 = k + j * b_dim1;
#line 304 "ctrsm.f"
			if (b[i__3].r != 0. || b[i__3].i != 0.) {
#line 305 "ctrsm.f"
			    if (nounit) {
#line 305 "ctrsm.f"
				i__3 = k + j * b_dim1;
#line 305 "ctrsm.f"
				z_div(&z__1, &b[k + j * b_dim1], &a[k + k * 
					a_dim1]);
#line 305 "ctrsm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 305 "ctrsm.f"
			    }
#line 306 "ctrsm.f"
			    i__3 = *m;
#line 306 "ctrsm.f"
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 307 "ctrsm.f"
				i__4 = i__ + j * b_dim1;
#line 307 "ctrsm.f"
				i__5 = i__ + j * b_dim1;
#line 307 "ctrsm.f"
				i__6 = k + j * b_dim1;
#line 307 "ctrsm.f"
				i__7 = i__ + k * a_dim1;
#line 307 "ctrsm.f"
				z__2.r = b[i__6].r * a[i__7].r - b[i__6].i * 
					a[i__7].i, z__2.i = b[i__6].r * a[
					i__7].i + b[i__6].i * a[i__7].r;
#line 307 "ctrsm.f"
				z__1.r = b[i__5].r - z__2.r, z__1.i = b[i__5]
					.i - z__2.i;
#line 307 "ctrsm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 308 "ctrsm.f"
/* L80: */
#line 308 "ctrsm.f"
			    }
#line 309 "ctrsm.f"
			}
#line 310 "ctrsm.f"
/* L90: */
#line 310 "ctrsm.f"
		    }
#line 311 "ctrsm.f"
/* L100: */
#line 311 "ctrsm.f"
		}
#line 312 "ctrsm.f"
	    }
#line 313 "ctrsm.f"
	} else {

/*           Form  B := alpha*inv( A**T )*B */
/*           or    B := alpha*inv( A**H )*B. */

#line 318 "ctrsm.f"
	    if (upper) {
#line 319 "ctrsm.f"
		i__1 = *n;
#line 319 "ctrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 320 "ctrsm.f"
		    i__2 = *m;
#line 320 "ctrsm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 321 "ctrsm.f"
			i__3 = i__ + j * b_dim1;
#line 321 "ctrsm.f"
			z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				z__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
#line 321 "ctrsm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 322 "ctrsm.f"
			if (noconj) {
#line 323 "ctrsm.f"
			    i__3 = i__ - 1;
#line 323 "ctrsm.f"
			    for (k = 1; k <= i__3; ++k) {
#line 324 "ctrsm.f"
				i__4 = k + i__ * a_dim1;
#line 324 "ctrsm.f"
				i__5 = k + j * b_dim1;
#line 324 "ctrsm.f"
				z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * 
					b[i__5].i, z__2.i = a[i__4].r * b[
					i__5].i + a[i__4].i * b[i__5].r;
#line 324 "ctrsm.f"
				z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
					z__2.i;
#line 324 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 325 "ctrsm.f"
/* L110: */
#line 325 "ctrsm.f"
			    }
#line 326 "ctrsm.f"
			    if (nounit) {
#line 326 "ctrsm.f"
				z_div(&z__1, &temp, &a[i__ + i__ * a_dim1]);
#line 326 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 326 "ctrsm.f"
			    }
#line 327 "ctrsm.f"
			} else {
#line 328 "ctrsm.f"
			    i__3 = i__ - 1;
#line 328 "ctrsm.f"
			    for (k = 1; k <= i__3; ++k) {
#line 329 "ctrsm.f"
				d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 329 "ctrsm.f"
				i__4 = k + j * b_dim1;
#line 329 "ctrsm.f"
				z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4]
					.i, z__2.i = z__3.r * b[i__4].i + 
					z__3.i * b[i__4].r;
#line 329 "ctrsm.f"
				z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
					z__2.i;
#line 329 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 330 "ctrsm.f"
/* L120: */
#line 330 "ctrsm.f"
			    }
#line 331 "ctrsm.f"
			    if (nounit) {
#line 331 "ctrsm.f"
				d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
#line 331 "ctrsm.f"
				z_div(&z__1, &temp, &z__2);
#line 331 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 331 "ctrsm.f"
			    }
#line 332 "ctrsm.f"
			}
#line 333 "ctrsm.f"
			i__3 = i__ + j * b_dim1;
#line 333 "ctrsm.f"
			b[i__3].r = temp.r, b[i__3].i = temp.i;
#line 334 "ctrsm.f"
/* L130: */
#line 334 "ctrsm.f"
		    }
#line 335 "ctrsm.f"
/* L140: */
#line 335 "ctrsm.f"
		}
#line 336 "ctrsm.f"
	    } else {
#line 337 "ctrsm.f"
		i__1 = *n;
#line 337 "ctrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 338 "ctrsm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 339 "ctrsm.f"
			i__2 = i__ + j * b_dim1;
#line 339 "ctrsm.f"
			z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2].i, 
				z__1.i = alpha->r * b[i__2].i + alpha->i * b[
				i__2].r;
#line 339 "ctrsm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 340 "ctrsm.f"
			if (noconj) {
#line 341 "ctrsm.f"
			    i__2 = *m;
#line 341 "ctrsm.f"
			    for (k = i__ + 1; k <= i__2; ++k) {
#line 342 "ctrsm.f"
				i__3 = k + i__ * a_dim1;
#line 342 "ctrsm.f"
				i__4 = k + j * b_dim1;
#line 342 "ctrsm.f"
				z__2.r = a[i__3].r * b[i__4].r - a[i__3].i * 
					b[i__4].i, z__2.i = a[i__3].r * b[
					i__4].i + a[i__3].i * b[i__4].r;
#line 342 "ctrsm.f"
				z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
					z__2.i;
#line 342 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 343 "ctrsm.f"
/* L150: */
#line 343 "ctrsm.f"
			    }
#line 344 "ctrsm.f"
			    if (nounit) {
#line 344 "ctrsm.f"
				z_div(&z__1, &temp, &a[i__ + i__ * a_dim1]);
#line 344 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 344 "ctrsm.f"
			    }
#line 345 "ctrsm.f"
			} else {
#line 346 "ctrsm.f"
			    i__2 = *m;
#line 346 "ctrsm.f"
			    for (k = i__ + 1; k <= i__2; ++k) {
#line 347 "ctrsm.f"
				d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 347 "ctrsm.f"
				i__3 = k + j * b_dim1;
#line 347 "ctrsm.f"
				z__2.r = z__3.r * b[i__3].r - z__3.i * b[i__3]
					.i, z__2.i = z__3.r * b[i__3].i + 
					z__3.i * b[i__3].r;
#line 347 "ctrsm.f"
				z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
					z__2.i;
#line 347 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 348 "ctrsm.f"
/* L160: */
#line 348 "ctrsm.f"
			    }
#line 349 "ctrsm.f"
			    if (nounit) {
#line 349 "ctrsm.f"
				d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
#line 349 "ctrsm.f"
				z_div(&z__1, &temp, &z__2);
#line 349 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 349 "ctrsm.f"
			    }
#line 350 "ctrsm.f"
			}
#line 351 "ctrsm.f"
			i__2 = i__ + j * b_dim1;
#line 351 "ctrsm.f"
			b[i__2].r = temp.r, b[i__2].i = temp.i;
#line 352 "ctrsm.f"
/* L170: */
#line 352 "ctrsm.f"
		    }
#line 353 "ctrsm.f"
/* L180: */
#line 353 "ctrsm.f"
		}
#line 354 "ctrsm.f"
	    }
#line 355 "ctrsm.f"
	}
#line 356 "ctrsm.f"
    } else {
#line 357 "ctrsm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*inv( A ). */

#line 361 "ctrsm.f"
	    if (upper) {
#line 362 "ctrsm.f"
		i__1 = *n;
#line 362 "ctrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 363 "ctrsm.f"
		    if (alpha->r != 1. || alpha->i != 0.) {
#line 364 "ctrsm.f"
			i__2 = *m;
#line 364 "ctrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 365 "ctrsm.f"
			    i__3 = i__ + j * b_dim1;
#line 365 "ctrsm.f"
			    i__4 = i__ + j * b_dim1;
#line 365 "ctrsm.f"
			    z__1.r = alpha->r * b[i__4].r - alpha->i * b[i__4]
				    .i, z__1.i = alpha->r * b[i__4].i + 
				    alpha->i * b[i__4].r;
#line 365 "ctrsm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 366 "ctrsm.f"
/* L190: */
#line 366 "ctrsm.f"
			}
#line 367 "ctrsm.f"
		    }
#line 368 "ctrsm.f"
		    i__2 = j - 1;
#line 368 "ctrsm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 369 "ctrsm.f"
			i__3 = k + j * a_dim1;
#line 369 "ctrsm.f"
			if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 370 "ctrsm.f"
			    i__3 = *m;
#line 370 "ctrsm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 371 "ctrsm.f"
				i__4 = i__ + j * b_dim1;
#line 371 "ctrsm.f"
				i__5 = i__ + j * b_dim1;
#line 371 "ctrsm.f"
				i__6 = k + j * a_dim1;
#line 371 "ctrsm.f"
				i__7 = i__ + k * b_dim1;
#line 371 "ctrsm.f"
				z__2.r = a[i__6].r * b[i__7].r - a[i__6].i * 
					b[i__7].i, z__2.i = a[i__6].r * b[
					i__7].i + a[i__6].i * b[i__7].r;
#line 371 "ctrsm.f"
				z__1.r = b[i__5].r - z__2.r, z__1.i = b[i__5]
					.i - z__2.i;
#line 371 "ctrsm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 372 "ctrsm.f"
/* L200: */
#line 372 "ctrsm.f"
			    }
#line 373 "ctrsm.f"
			}
#line 374 "ctrsm.f"
/* L210: */
#line 374 "ctrsm.f"
		    }
#line 375 "ctrsm.f"
		    if (nounit) {
#line 376 "ctrsm.f"
			z_div(&z__1, &c_b1, &a[j + j * a_dim1]);
#line 376 "ctrsm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 377 "ctrsm.f"
			i__2 = *m;
#line 377 "ctrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 378 "ctrsm.f"
			    i__3 = i__ + j * b_dim1;
#line 378 "ctrsm.f"
			    i__4 = i__ + j * b_dim1;
#line 378 "ctrsm.f"
			    z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i, 
				    z__1.i = temp.r * b[i__4].i + temp.i * b[
				    i__4].r;
#line 378 "ctrsm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 379 "ctrsm.f"
/* L220: */
#line 379 "ctrsm.f"
			}
#line 380 "ctrsm.f"
		    }
#line 381 "ctrsm.f"
/* L230: */
#line 381 "ctrsm.f"
		}
#line 382 "ctrsm.f"
	    } else {
#line 383 "ctrsm.f"
		for (j = *n; j >= 1; --j) {
#line 384 "ctrsm.f"
		    if (alpha->r != 1. || alpha->i != 0.) {
#line 385 "ctrsm.f"
			i__1 = *m;
#line 385 "ctrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "ctrsm.f"
			    i__2 = i__ + j * b_dim1;
#line 386 "ctrsm.f"
			    i__3 = i__ + j * b_dim1;
#line 386 "ctrsm.f"
			    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3]
				    .i, z__1.i = alpha->r * b[i__3].i + 
				    alpha->i * b[i__3].r;
#line 386 "ctrsm.f"
			    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 387 "ctrsm.f"
/* L240: */
#line 387 "ctrsm.f"
			}
#line 388 "ctrsm.f"
		    }
#line 389 "ctrsm.f"
		    i__1 = *n;
#line 389 "ctrsm.f"
		    for (k = j + 1; k <= i__1; ++k) {
#line 390 "ctrsm.f"
			i__2 = k + j * a_dim1;
#line 390 "ctrsm.f"
			if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 391 "ctrsm.f"
			    i__2 = *m;
#line 391 "ctrsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 392 "ctrsm.f"
				i__3 = i__ + j * b_dim1;
#line 392 "ctrsm.f"
				i__4 = i__ + j * b_dim1;
#line 392 "ctrsm.f"
				i__5 = k + j * a_dim1;
#line 392 "ctrsm.f"
				i__6 = i__ + k * b_dim1;
#line 392 "ctrsm.f"
				z__2.r = a[i__5].r * b[i__6].r - a[i__5].i * 
					b[i__6].i, z__2.i = a[i__5].r * b[
					i__6].i + a[i__5].i * b[i__6].r;
#line 392 "ctrsm.f"
				z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4]
					.i - z__2.i;
#line 392 "ctrsm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 393 "ctrsm.f"
/* L250: */
#line 393 "ctrsm.f"
			    }
#line 394 "ctrsm.f"
			}
#line 395 "ctrsm.f"
/* L260: */
#line 395 "ctrsm.f"
		    }
#line 396 "ctrsm.f"
		    if (nounit) {
#line 397 "ctrsm.f"
			z_div(&z__1, &c_b1, &a[j + j * a_dim1]);
#line 397 "ctrsm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 398 "ctrsm.f"
			i__1 = *m;
#line 398 "ctrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 399 "ctrsm.f"
			    i__2 = i__ + j * b_dim1;
#line 399 "ctrsm.f"
			    i__3 = i__ + j * b_dim1;
#line 399 "ctrsm.f"
			    z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i, 
				    z__1.i = temp.r * b[i__3].i + temp.i * b[
				    i__3].r;
#line 399 "ctrsm.f"
			    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 400 "ctrsm.f"
/* L270: */
#line 400 "ctrsm.f"
			}
#line 401 "ctrsm.f"
		    }
#line 402 "ctrsm.f"
/* L280: */
#line 402 "ctrsm.f"
		}
#line 403 "ctrsm.f"
	    }
#line 404 "ctrsm.f"
	} else {

/*           Form  B := alpha*B*inv( A**T ) */
/*           or    B := alpha*B*inv( A**H ). */

#line 409 "ctrsm.f"
	    if (upper) {
#line 410 "ctrsm.f"
		for (k = *n; k >= 1; --k) {
#line 411 "ctrsm.f"
		    if (nounit) {
#line 412 "ctrsm.f"
			if (noconj) {
#line 413 "ctrsm.f"
			    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 413 "ctrsm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 414 "ctrsm.f"
			} else {
#line 415 "ctrsm.f"
			    d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 415 "ctrsm.f"
			    z_div(&z__1, &c_b1, &z__2);
#line 415 "ctrsm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 416 "ctrsm.f"
			}
#line 417 "ctrsm.f"
			i__1 = *m;
#line 417 "ctrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 418 "ctrsm.f"
			    i__2 = i__ + k * b_dim1;
#line 418 "ctrsm.f"
			    i__3 = i__ + k * b_dim1;
#line 418 "ctrsm.f"
			    z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i, 
				    z__1.i = temp.r * b[i__3].i + temp.i * b[
				    i__3].r;
#line 418 "ctrsm.f"
			    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 419 "ctrsm.f"
/* L290: */
#line 419 "ctrsm.f"
			}
#line 420 "ctrsm.f"
		    }
#line 421 "ctrsm.f"
		    i__1 = k - 1;
#line 421 "ctrsm.f"
		    for (j = 1; j <= i__1; ++j) {
#line 422 "ctrsm.f"
			i__2 = j + k * a_dim1;
#line 422 "ctrsm.f"
			if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 423 "ctrsm.f"
			    if (noconj) {
#line 424 "ctrsm.f"
				i__2 = j + k * a_dim1;
#line 424 "ctrsm.f"
				temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 425 "ctrsm.f"
			    } else {
#line 426 "ctrsm.f"
				d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 426 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 427 "ctrsm.f"
			    }
#line 428 "ctrsm.f"
			    i__2 = *m;
#line 428 "ctrsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 429 "ctrsm.f"
				i__3 = i__ + j * b_dim1;
#line 429 "ctrsm.f"
				i__4 = i__ + j * b_dim1;
#line 429 "ctrsm.f"
				i__5 = i__ + k * b_dim1;
#line 429 "ctrsm.f"
				z__2.r = temp.r * b[i__5].r - temp.i * b[i__5]
					.i, z__2.i = temp.r * b[i__5].i + 
					temp.i * b[i__5].r;
#line 429 "ctrsm.f"
				z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4]
					.i - z__2.i;
#line 429 "ctrsm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 430 "ctrsm.f"
/* L300: */
#line 430 "ctrsm.f"
			    }
#line 431 "ctrsm.f"
			}
#line 432 "ctrsm.f"
/* L310: */
#line 432 "ctrsm.f"
		    }
#line 433 "ctrsm.f"
		    if (alpha->r != 1. || alpha->i != 0.) {
#line 434 "ctrsm.f"
			i__1 = *m;
#line 434 "ctrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "ctrsm.f"
			    i__2 = i__ + k * b_dim1;
#line 435 "ctrsm.f"
			    i__3 = i__ + k * b_dim1;
#line 435 "ctrsm.f"
			    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3]
				    .i, z__1.i = alpha->r * b[i__3].i + 
				    alpha->i * b[i__3].r;
#line 435 "ctrsm.f"
			    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 436 "ctrsm.f"
/* L320: */
#line 436 "ctrsm.f"
			}
#line 437 "ctrsm.f"
		    }
#line 438 "ctrsm.f"
/* L330: */
#line 438 "ctrsm.f"
		}
#line 439 "ctrsm.f"
	    } else {
#line 440 "ctrsm.f"
		i__1 = *n;
#line 440 "ctrsm.f"
		for (k = 1; k <= i__1; ++k) {
#line 441 "ctrsm.f"
		    if (nounit) {
#line 442 "ctrsm.f"
			if (noconj) {
#line 443 "ctrsm.f"
			    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 443 "ctrsm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 444 "ctrsm.f"
			} else {
#line 445 "ctrsm.f"
			    d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 445 "ctrsm.f"
			    z_div(&z__1, &c_b1, &z__2);
#line 445 "ctrsm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 446 "ctrsm.f"
			}
#line 447 "ctrsm.f"
			i__2 = *m;
#line 447 "ctrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "ctrsm.f"
			    i__3 = i__ + k * b_dim1;
#line 448 "ctrsm.f"
			    i__4 = i__ + k * b_dim1;
#line 448 "ctrsm.f"
			    z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i, 
				    z__1.i = temp.r * b[i__4].i + temp.i * b[
				    i__4].r;
#line 448 "ctrsm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 449 "ctrsm.f"
/* L340: */
#line 449 "ctrsm.f"
			}
#line 450 "ctrsm.f"
		    }
#line 451 "ctrsm.f"
		    i__2 = *n;
#line 451 "ctrsm.f"
		    for (j = k + 1; j <= i__2; ++j) {
#line 452 "ctrsm.f"
			i__3 = j + k * a_dim1;
#line 452 "ctrsm.f"
			if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 453 "ctrsm.f"
			    if (noconj) {
#line 454 "ctrsm.f"
				i__3 = j + k * a_dim1;
#line 454 "ctrsm.f"
				temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 455 "ctrsm.f"
			    } else {
#line 456 "ctrsm.f"
				d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 456 "ctrsm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 457 "ctrsm.f"
			    }
#line 458 "ctrsm.f"
			    i__3 = *m;
#line 458 "ctrsm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 459 "ctrsm.f"
				i__4 = i__ + j * b_dim1;
#line 459 "ctrsm.f"
				i__5 = i__ + j * b_dim1;
#line 459 "ctrsm.f"
				i__6 = i__ + k * b_dim1;
#line 459 "ctrsm.f"
				z__2.r = temp.r * b[i__6].r - temp.i * b[i__6]
					.i, z__2.i = temp.r * b[i__6].i + 
					temp.i * b[i__6].r;
#line 459 "ctrsm.f"
				z__1.r = b[i__5].r - z__2.r, z__1.i = b[i__5]
					.i - z__2.i;
#line 459 "ctrsm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 460 "ctrsm.f"
/* L350: */
#line 460 "ctrsm.f"
			    }
#line 461 "ctrsm.f"
			}
#line 462 "ctrsm.f"
/* L360: */
#line 462 "ctrsm.f"
		    }
#line 463 "ctrsm.f"
		    if (alpha->r != 1. || alpha->i != 0.) {
#line 464 "ctrsm.f"
			i__2 = *m;
#line 464 "ctrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 465 "ctrsm.f"
			    i__3 = i__ + k * b_dim1;
#line 465 "ctrsm.f"
			    i__4 = i__ + k * b_dim1;
#line 465 "ctrsm.f"
			    z__1.r = alpha->r * b[i__4].r - alpha->i * b[i__4]
				    .i, z__1.i = alpha->r * b[i__4].i + 
				    alpha->i * b[i__4].r;
#line 465 "ctrsm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 466 "ctrsm.f"
/* L370: */
#line 466 "ctrsm.f"
			}
#line 467 "ctrsm.f"
		    }
#line 468 "ctrsm.f"
/* L380: */
#line 468 "ctrsm.f"
		}
#line 469 "ctrsm.f"
	    }
#line 470 "ctrsm.f"
	}
#line 471 "ctrsm.f"
    }

#line 473 "ctrsm.f"
    return 0;

/*     End of CTRSM . */

} /* ctrsm_ */

