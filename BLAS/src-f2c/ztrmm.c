#line 1 "ztrmm.f"
/* ztrmm.f -- translated by f2c (version 20100827).
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

#line 1 "ztrmm.f"
/* > \brief \b ZTRMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA */
/*       INTEGER LDA,LDB,M,N */
/*       CHARACTER DIAG,SIDE,TRANSA,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),B(LDB,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRMM  performs one of the matrix-matrix operations */
/* > */
/* >    B := alpha*op( A )*B,   or   B := alpha*B*op( A ) */
/* > */
/* > where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE specifies whether  op( A ) multiplies B from */
/* >           the left or right as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */
/* > */
/* >              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */
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
/* >          ALPHA is COMPLEX*16 */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, k ), where k is m */
/* >           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
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
/* >          B is COMPLEX*16 array, dimension ( LDB, N ). */
/* >           Before entry,  the leading  m by n part of the array  B must */
/* >           contain the matrix  B,  and  on exit  is overwritten  by the */
/* >           transformed matrix. */
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

/* > \date December 2016 */

/* > \ingroup complex16_blas_level3 */

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
/* Subroutine */ int ztrmm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, ftnlen side_len, ftnlen 
	uplo_len, ftnlen transa_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, info;
    static doublecomplex temp;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconj, nounit;


/*  -- Reference BLAS level3 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 220 "ztrmm.f"
    /* Parameter adjustments */
#line 220 "ztrmm.f"
    a_dim1 = *lda;
#line 220 "ztrmm.f"
    a_offset = 1 + a_dim1;
#line 220 "ztrmm.f"
    a -= a_offset;
#line 220 "ztrmm.f"
    b_dim1 = *ldb;
#line 220 "ztrmm.f"
    b_offset = 1 + b_dim1;
#line 220 "ztrmm.f"
    b -= b_offset;
#line 220 "ztrmm.f"

#line 220 "ztrmm.f"
    /* Function Body */
#line 220 "ztrmm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 221 "ztrmm.f"
    if (lside) {
#line 222 "ztrmm.f"
	nrowa = *m;
#line 223 "ztrmm.f"
    } else {
#line 224 "ztrmm.f"
	nrowa = *n;
#line 225 "ztrmm.f"
    }
#line 226 "ztrmm.f"
    noconj = lsame_(transa, "T", (ftnlen)1, (ftnlen)1);
#line 227 "ztrmm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 228 "ztrmm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 230 "ztrmm.f"
    info = 0;
#line 231 "ztrmm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 232 "ztrmm.f"
	info = 1;
#line 233 "ztrmm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 234 "ztrmm.f"
	info = 2;
#line 235 "ztrmm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 238 "ztrmm.f"
	info = 3;
#line 239 "ztrmm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 240 "ztrmm.f"
	info = 4;
#line 241 "ztrmm.f"
    } else if (*m < 0) {
#line 242 "ztrmm.f"
	info = 5;
#line 243 "ztrmm.f"
    } else if (*n < 0) {
#line 244 "ztrmm.f"
	info = 6;
#line 245 "ztrmm.f"
    } else if (*lda < max(1,nrowa)) {
#line 246 "ztrmm.f"
	info = 9;
#line 247 "ztrmm.f"
    } else if (*ldb < max(1,*m)) {
#line 248 "ztrmm.f"
	info = 11;
#line 249 "ztrmm.f"
    }
#line 250 "ztrmm.f"
    if (info != 0) {
#line 251 "ztrmm.f"
	xerbla_("ZTRMM ", &info, (ftnlen)6);
#line 252 "ztrmm.f"
	return 0;
#line 253 "ztrmm.f"
    }

/*     Quick return if possible. */

#line 257 "ztrmm.f"
    if (*m == 0 || *n == 0) {
#line 257 "ztrmm.f"
	return 0;
#line 257 "ztrmm.f"
    }

/*     And when  alpha.eq.zero. */

#line 261 "ztrmm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 262 "ztrmm.f"
	i__1 = *n;
#line 262 "ztrmm.f"
	for (j = 1; j <= i__1; ++j) {
#line 263 "ztrmm.f"
	    i__2 = *m;
#line 263 "ztrmm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 264 "ztrmm.f"
		i__3 = i__ + j * b_dim1;
#line 264 "ztrmm.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 265 "ztrmm.f"
/* L10: */
#line 265 "ztrmm.f"
	    }
#line 266 "ztrmm.f"
/* L20: */
#line 266 "ztrmm.f"
	}
#line 267 "ztrmm.f"
	return 0;
#line 268 "ztrmm.f"
    }

/*     Start the operations. */

#line 272 "ztrmm.f"
    if (lside) {
#line 273 "ztrmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*A*B. */

#line 277 "ztrmm.f"
	    if (upper) {
#line 278 "ztrmm.f"
		i__1 = *n;
#line 278 "ztrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 279 "ztrmm.f"
		    i__2 = *m;
#line 279 "ztrmm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 280 "ztrmm.f"
			i__3 = k + j * b_dim1;
#line 280 "ztrmm.f"
			if (b[i__3].r != 0. || b[i__3].i != 0.) {
#line 281 "ztrmm.f"
			    i__3 = k + j * b_dim1;
#line 281 "ztrmm.f"
			    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3]
				    .i, z__1.i = alpha->r * b[i__3].i + 
				    alpha->i * b[i__3].r;
#line 281 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 282 "ztrmm.f"
			    i__3 = k - 1;
#line 282 "ztrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 283 "ztrmm.f"
				i__4 = i__ + j * b_dim1;
#line 283 "ztrmm.f"
				i__5 = i__ + j * b_dim1;
#line 283 "ztrmm.f"
				i__6 = i__ + k * a_dim1;
#line 283 "ztrmm.f"
				z__2.r = temp.r * a[i__6].r - temp.i * a[i__6]
					.i, z__2.i = temp.r * a[i__6].i + 
					temp.i * a[i__6].r;
#line 283 "ztrmm.f"
				z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5]
					.i + z__2.i;
#line 283 "ztrmm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 284 "ztrmm.f"
/* L30: */
#line 284 "ztrmm.f"
			    }
#line 285 "ztrmm.f"
			    if (nounit) {
#line 285 "ztrmm.f"
				i__3 = k + k * a_dim1;
#line 285 "ztrmm.f"
				z__1.r = temp.r * a[i__3].r - temp.i * a[i__3]
					.i, z__1.i = temp.r * a[i__3].i + 
					temp.i * a[i__3].r;
#line 285 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 285 "ztrmm.f"
			    }
#line 286 "ztrmm.f"
			    i__3 = k + j * b_dim1;
#line 286 "ztrmm.f"
			    b[i__3].r = temp.r, b[i__3].i = temp.i;
#line 287 "ztrmm.f"
			}
#line 288 "ztrmm.f"
/* L40: */
#line 288 "ztrmm.f"
		    }
#line 289 "ztrmm.f"
/* L50: */
#line 289 "ztrmm.f"
		}
#line 290 "ztrmm.f"
	    } else {
#line 291 "ztrmm.f"
		i__1 = *n;
#line 291 "ztrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "ztrmm.f"
		    for (k = *m; k >= 1; --k) {
#line 293 "ztrmm.f"
			i__2 = k + j * b_dim1;
#line 293 "ztrmm.f"
			if (b[i__2].r != 0. || b[i__2].i != 0.) {
#line 294 "ztrmm.f"
			    i__2 = k + j * b_dim1;
#line 294 "ztrmm.f"
			    z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2]
				    .i, z__1.i = alpha->r * b[i__2].i + 
				    alpha->i * b[i__2].r;
#line 294 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 295 "ztrmm.f"
			    i__2 = k + j * b_dim1;
#line 295 "ztrmm.f"
			    b[i__2].r = temp.r, b[i__2].i = temp.i;
#line 296 "ztrmm.f"
			    if (nounit) {
#line 296 "ztrmm.f"
				i__2 = k + j * b_dim1;
#line 296 "ztrmm.f"
				i__3 = k + j * b_dim1;
#line 296 "ztrmm.f"
				i__4 = k + k * a_dim1;
#line 296 "ztrmm.f"
				z__1.r = b[i__3].r * a[i__4].r - b[i__3].i * 
					a[i__4].i, z__1.i = b[i__3].r * a[
					i__4].i + b[i__3].i * a[i__4].r;
#line 296 "ztrmm.f"
				b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 296 "ztrmm.f"
			    }
#line 297 "ztrmm.f"
			    i__2 = *m;
#line 297 "ztrmm.f"
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 298 "ztrmm.f"
				i__3 = i__ + j * b_dim1;
#line 298 "ztrmm.f"
				i__4 = i__ + j * b_dim1;
#line 298 "ztrmm.f"
				i__5 = i__ + k * a_dim1;
#line 298 "ztrmm.f"
				z__2.r = temp.r * a[i__5].r - temp.i * a[i__5]
					.i, z__2.i = temp.r * a[i__5].i + 
					temp.i * a[i__5].r;
#line 298 "ztrmm.f"
				z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4]
					.i + z__2.i;
#line 298 "ztrmm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 299 "ztrmm.f"
/* L60: */
#line 299 "ztrmm.f"
			    }
#line 300 "ztrmm.f"
			}
#line 301 "ztrmm.f"
/* L70: */
#line 301 "ztrmm.f"
		    }
#line 302 "ztrmm.f"
/* L80: */
#line 302 "ztrmm.f"
		}
#line 303 "ztrmm.f"
	    }
#line 304 "ztrmm.f"
	} else {

/*           Form  B := alpha*A**T*B   or   B := alpha*A**H*B. */

#line 308 "ztrmm.f"
	    if (upper) {
#line 309 "ztrmm.f"
		i__1 = *n;
#line 309 "ztrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 310 "ztrmm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 311 "ztrmm.f"
			i__2 = i__ + j * b_dim1;
#line 311 "ztrmm.f"
			temp.r = b[i__2].r, temp.i = b[i__2].i;
#line 312 "ztrmm.f"
			if (noconj) {
#line 313 "ztrmm.f"
			    if (nounit) {
#line 313 "ztrmm.f"
				i__2 = i__ + i__ * a_dim1;
#line 313 "ztrmm.f"
				z__1.r = temp.r * a[i__2].r - temp.i * a[i__2]
					.i, z__1.i = temp.r * a[i__2].i + 
					temp.i * a[i__2].r;
#line 313 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 313 "ztrmm.f"
			    }
#line 314 "ztrmm.f"
			    i__2 = i__ - 1;
#line 314 "ztrmm.f"
			    for (k = 1; k <= i__2; ++k) {
#line 315 "ztrmm.f"
				i__3 = k + i__ * a_dim1;
#line 315 "ztrmm.f"
				i__4 = k + j * b_dim1;
#line 315 "ztrmm.f"
				z__2.r = a[i__3].r * b[i__4].r - a[i__3].i * 
					b[i__4].i, z__2.i = a[i__3].r * b[
					i__4].i + a[i__3].i * b[i__4].r;
#line 315 "ztrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 315 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 316 "ztrmm.f"
/* L90: */
#line 316 "ztrmm.f"
			    }
#line 317 "ztrmm.f"
			} else {
#line 318 "ztrmm.f"
			    if (nounit) {
#line 318 "ztrmm.f"
				d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
#line 318 "ztrmm.f"
				z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
					z__1.i = temp.r * z__2.i + temp.i * 
					z__2.r;
#line 318 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 318 "ztrmm.f"
			    }
#line 319 "ztrmm.f"
			    i__2 = i__ - 1;
#line 319 "ztrmm.f"
			    for (k = 1; k <= i__2; ++k) {
#line 320 "ztrmm.f"
				d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 320 "ztrmm.f"
				i__3 = k + j * b_dim1;
#line 320 "ztrmm.f"
				z__2.r = z__3.r * b[i__3].r - z__3.i * b[i__3]
					.i, z__2.i = z__3.r * b[i__3].i + 
					z__3.i * b[i__3].r;
#line 320 "ztrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 320 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 321 "ztrmm.f"
/* L100: */
#line 321 "ztrmm.f"
			    }
#line 322 "ztrmm.f"
			}
#line 323 "ztrmm.f"
			i__2 = i__ + j * b_dim1;
#line 323 "ztrmm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 323 "ztrmm.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 324 "ztrmm.f"
/* L110: */
#line 324 "ztrmm.f"
		    }
#line 325 "ztrmm.f"
/* L120: */
#line 325 "ztrmm.f"
		}
#line 326 "ztrmm.f"
	    } else {
#line 327 "ztrmm.f"
		i__1 = *n;
#line 327 "ztrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 328 "ztrmm.f"
		    i__2 = *m;
#line 328 "ztrmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 329 "ztrmm.f"
			i__3 = i__ + j * b_dim1;
#line 329 "ztrmm.f"
			temp.r = b[i__3].r, temp.i = b[i__3].i;
#line 330 "ztrmm.f"
			if (noconj) {
#line 331 "ztrmm.f"
			    if (nounit) {
#line 331 "ztrmm.f"
				i__3 = i__ + i__ * a_dim1;
#line 331 "ztrmm.f"
				z__1.r = temp.r * a[i__3].r - temp.i * a[i__3]
					.i, z__1.i = temp.r * a[i__3].i + 
					temp.i * a[i__3].r;
#line 331 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 331 "ztrmm.f"
			    }
#line 332 "ztrmm.f"
			    i__3 = *m;
#line 332 "ztrmm.f"
			    for (k = i__ + 1; k <= i__3; ++k) {
#line 333 "ztrmm.f"
				i__4 = k + i__ * a_dim1;
#line 333 "ztrmm.f"
				i__5 = k + j * b_dim1;
#line 333 "ztrmm.f"
				z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * 
					b[i__5].i, z__2.i = a[i__4].r * b[
					i__5].i + a[i__4].i * b[i__5].r;
#line 333 "ztrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 333 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 334 "ztrmm.f"
/* L130: */
#line 334 "ztrmm.f"
			    }
#line 335 "ztrmm.f"
			} else {
#line 336 "ztrmm.f"
			    if (nounit) {
#line 336 "ztrmm.f"
				d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
#line 336 "ztrmm.f"
				z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
					z__1.i = temp.r * z__2.i + temp.i * 
					z__2.r;
#line 336 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 336 "ztrmm.f"
			    }
#line 337 "ztrmm.f"
			    i__3 = *m;
#line 337 "ztrmm.f"
			    for (k = i__ + 1; k <= i__3; ++k) {
#line 338 "ztrmm.f"
				d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 338 "ztrmm.f"
				i__4 = k + j * b_dim1;
#line 338 "ztrmm.f"
				z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4]
					.i, z__2.i = z__3.r * b[i__4].i + 
					z__3.i * b[i__4].r;
#line 338 "ztrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 338 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 339 "ztrmm.f"
/* L140: */
#line 339 "ztrmm.f"
			    }
#line 340 "ztrmm.f"
			}
#line 341 "ztrmm.f"
			i__3 = i__ + j * b_dim1;
#line 341 "ztrmm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 341 "ztrmm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 342 "ztrmm.f"
/* L150: */
#line 342 "ztrmm.f"
		    }
#line 343 "ztrmm.f"
/* L160: */
#line 343 "ztrmm.f"
		}
#line 344 "ztrmm.f"
	    }
#line 345 "ztrmm.f"
	}
#line 346 "ztrmm.f"
    } else {
#line 347 "ztrmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*A. */

#line 351 "ztrmm.f"
	    if (upper) {
#line 352 "ztrmm.f"
		for (j = *n; j >= 1; --j) {
#line 353 "ztrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 354 "ztrmm.f"
		    if (nounit) {
#line 354 "ztrmm.f"
			i__1 = j + j * a_dim1;
#line 354 "ztrmm.f"
			z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				z__1.i = temp.r * a[i__1].i + temp.i * a[i__1]
				.r;
#line 354 "ztrmm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ztrmm.f"
		    }
#line 355 "ztrmm.f"
		    i__1 = *m;
#line 355 "ztrmm.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 356 "ztrmm.f"
			i__2 = i__ + j * b_dim1;
#line 356 "ztrmm.f"
			i__3 = i__ + j * b_dim1;
#line 356 "ztrmm.f"
			z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i, 
				z__1.i = temp.r * b[i__3].i + temp.i * b[i__3]
				.r;
#line 356 "ztrmm.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 357 "ztrmm.f"
/* L170: */
#line 357 "ztrmm.f"
		    }
#line 358 "ztrmm.f"
		    i__1 = j - 1;
#line 358 "ztrmm.f"
		    for (k = 1; k <= i__1; ++k) {
#line 359 "ztrmm.f"
			i__2 = k + j * a_dim1;
#line 359 "ztrmm.f"
			if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 360 "ztrmm.f"
			    i__2 = k + j * a_dim1;
#line 360 "ztrmm.f"
			    z__1.r = alpha->r * a[i__2].r - alpha->i * a[i__2]
				    .i, z__1.i = alpha->r * a[i__2].i + 
				    alpha->i * a[i__2].r;
#line 360 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 361 "ztrmm.f"
			    i__2 = *m;
#line 361 "ztrmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 362 "ztrmm.f"
				i__3 = i__ + j * b_dim1;
#line 362 "ztrmm.f"
				i__4 = i__ + j * b_dim1;
#line 362 "ztrmm.f"
				i__5 = i__ + k * b_dim1;
#line 362 "ztrmm.f"
				z__2.r = temp.r * b[i__5].r - temp.i * b[i__5]
					.i, z__2.i = temp.r * b[i__5].i + 
					temp.i * b[i__5].r;
#line 362 "ztrmm.f"
				z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4]
					.i + z__2.i;
#line 362 "ztrmm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 363 "ztrmm.f"
/* L180: */
#line 363 "ztrmm.f"
			    }
#line 364 "ztrmm.f"
			}
#line 365 "ztrmm.f"
/* L190: */
#line 365 "ztrmm.f"
		    }
#line 366 "ztrmm.f"
/* L200: */
#line 366 "ztrmm.f"
		}
#line 367 "ztrmm.f"
	    } else {
#line 368 "ztrmm.f"
		i__1 = *n;
#line 368 "ztrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 369 "ztrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 370 "ztrmm.f"
		    if (nounit) {
#line 370 "ztrmm.f"
			i__2 = j + j * a_dim1;
#line 370 "ztrmm.f"
			z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				z__1.i = temp.r * a[i__2].i + temp.i * a[i__2]
				.r;
#line 370 "ztrmm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ztrmm.f"
		    }
#line 371 "ztrmm.f"
		    i__2 = *m;
#line 371 "ztrmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 372 "ztrmm.f"
			i__3 = i__ + j * b_dim1;
#line 372 "ztrmm.f"
			i__4 = i__ + j * b_dim1;
#line 372 "ztrmm.f"
			z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i, 
				z__1.i = temp.r * b[i__4].i + temp.i * b[i__4]
				.r;
#line 372 "ztrmm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 373 "ztrmm.f"
/* L210: */
#line 373 "ztrmm.f"
		    }
#line 374 "ztrmm.f"
		    i__2 = *n;
#line 374 "ztrmm.f"
		    for (k = j + 1; k <= i__2; ++k) {
#line 375 "ztrmm.f"
			i__3 = k + j * a_dim1;
#line 375 "ztrmm.f"
			if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 376 "ztrmm.f"
			    i__3 = k + j * a_dim1;
#line 376 "ztrmm.f"
			    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3]
				    .i, z__1.i = alpha->r * a[i__3].i + 
				    alpha->i * a[i__3].r;
#line 376 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 377 "ztrmm.f"
			    i__3 = *m;
#line 377 "ztrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 378 "ztrmm.f"
				i__4 = i__ + j * b_dim1;
#line 378 "ztrmm.f"
				i__5 = i__ + j * b_dim1;
#line 378 "ztrmm.f"
				i__6 = i__ + k * b_dim1;
#line 378 "ztrmm.f"
				z__2.r = temp.r * b[i__6].r - temp.i * b[i__6]
					.i, z__2.i = temp.r * b[i__6].i + 
					temp.i * b[i__6].r;
#line 378 "ztrmm.f"
				z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5]
					.i + z__2.i;
#line 378 "ztrmm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 379 "ztrmm.f"
/* L220: */
#line 379 "ztrmm.f"
			    }
#line 380 "ztrmm.f"
			}
#line 381 "ztrmm.f"
/* L230: */
#line 381 "ztrmm.f"
		    }
#line 382 "ztrmm.f"
/* L240: */
#line 382 "ztrmm.f"
		}
#line 383 "ztrmm.f"
	    }
#line 384 "ztrmm.f"
	} else {

/*           Form  B := alpha*B*A**T   or   B := alpha*B*A**H. */

#line 388 "ztrmm.f"
	    if (upper) {
#line 389 "ztrmm.f"
		i__1 = *n;
#line 389 "ztrmm.f"
		for (k = 1; k <= i__1; ++k) {
#line 390 "ztrmm.f"
		    i__2 = k - 1;
#line 390 "ztrmm.f"
		    for (j = 1; j <= i__2; ++j) {
#line 391 "ztrmm.f"
			i__3 = j + k * a_dim1;
#line 391 "ztrmm.f"
			if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 392 "ztrmm.f"
			    if (noconj) {
#line 393 "ztrmm.f"
				i__3 = j + k * a_dim1;
#line 393 "ztrmm.f"
				z__1.r = alpha->r * a[i__3].r - alpha->i * a[
					i__3].i, z__1.i = alpha->r * a[i__3]
					.i + alpha->i * a[i__3].r;
#line 393 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 394 "ztrmm.f"
			    } else {
#line 395 "ztrmm.f"
				d_cnjg(&z__2, &a[j + k * a_dim1]);
#line 395 "ztrmm.f"
				z__1.r = alpha->r * z__2.r - alpha->i * 
					z__2.i, z__1.i = alpha->r * z__2.i + 
					alpha->i * z__2.r;
#line 395 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 396 "ztrmm.f"
			    }
#line 397 "ztrmm.f"
			    i__3 = *m;
#line 397 "ztrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 398 "ztrmm.f"
				i__4 = i__ + j * b_dim1;
#line 398 "ztrmm.f"
				i__5 = i__ + j * b_dim1;
#line 398 "ztrmm.f"
				i__6 = i__ + k * b_dim1;
#line 398 "ztrmm.f"
				z__2.r = temp.r * b[i__6].r - temp.i * b[i__6]
					.i, z__2.i = temp.r * b[i__6].i + 
					temp.i * b[i__6].r;
#line 398 "ztrmm.f"
				z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5]
					.i + z__2.i;
#line 398 "ztrmm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 399 "ztrmm.f"
/* L250: */
#line 399 "ztrmm.f"
			    }
#line 400 "ztrmm.f"
			}
#line 401 "ztrmm.f"
/* L260: */
#line 401 "ztrmm.f"
		    }
#line 402 "ztrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 403 "ztrmm.f"
		    if (nounit) {
#line 404 "ztrmm.f"
			if (noconj) {
#line 405 "ztrmm.f"
			    i__2 = k + k * a_dim1;
#line 405 "ztrmm.f"
			    z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__1.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 405 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 406 "ztrmm.f"
			} else {
#line 407 "ztrmm.f"
			    d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 407 "ztrmm.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 407 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 408 "ztrmm.f"
			}
#line 409 "ztrmm.f"
		    }
#line 410 "ztrmm.f"
		    if (temp.r != 1. || temp.i != 0.) {
#line 411 "ztrmm.f"
			i__2 = *m;
#line 411 "ztrmm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 412 "ztrmm.f"
			    i__3 = i__ + k * b_dim1;
#line 412 "ztrmm.f"
			    i__4 = i__ + k * b_dim1;
#line 412 "ztrmm.f"
			    z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i, 
				    z__1.i = temp.r * b[i__4].i + temp.i * b[
				    i__4].r;
#line 412 "ztrmm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 413 "ztrmm.f"
/* L270: */
#line 413 "ztrmm.f"
			}
#line 414 "ztrmm.f"
		    }
#line 415 "ztrmm.f"
/* L280: */
#line 415 "ztrmm.f"
		}
#line 416 "ztrmm.f"
	    } else {
#line 417 "ztrmm.f"
		for (k = *n; k >= 1; --k) {
#line 418 "ztrmm.f"
		    i__1 = *n;
#line 418 "ztrmm.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 419 "ztrmm.f"
			i__2 = j + k * a_dim1;
#line 419 "ztrmm.f"
			if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 420 "ztrmm.f"
			    if (noconj) {
#line 421 "ztrmm.f"
				i__2 = j + k * a_dim1;
#line 421 "ztrmm.f"
				z__1.r = alpha->r * a[i__2].r - alpha->i * a[
					i__2].i, z__1.i = alpha->r * a[i__2]
					.i + alpha->i * a[i__2].r;
#line 421 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 422 "ztrmm.f"
			    } else {
#line 423 "ztrmm.f"
				d_cnjg(&z__2, &a[j + k * a_dim1]);
#line 423 "ztrmm.f"
				z__1.r = alpha->r * z__2.r - alpha->i * 
					z__2.i, z__1.i = alpha->r * z__2.i + 
					alpha->i * z__2.r;
#line 423 "ztrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 424 "ztrmm.f"
			    }
#line 425 "ztrmm.f"
			    i__2 = *m;
#line 425 "ztrmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 426 "ztrmm.f"
				i__3 = i__ + j * b_dim1;
#line 426 "ztrmm.f"
				i__4 = i__ + j * b_dim1;
#line 426 "ztrmm.f"
				i__5 = i__ + k * b_dim1;
#line 426 "ztrmm.f"
				z__2.r = temp.r * b[i__5].r - temp.i * b[i__5]
					.i, z__2.i = temp.r * b[i__5].i + 
					temp.i * b[i__5].r;
#line 426 "ztrmm.f"
				z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4]
					.i + z__2.i;
#line 426 "ztrmm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 427 "ztrmm.f"
/* L290: */
#line 427 "ztrmm.f"
			    }
#line 428 "ztrmm.f"
			}
#line 429 "ztrmm.f"
/* L300: */
#line 429 "ztrmm.f"
		    }
#line 430 "ztrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 431 "ztrmm.f"
		    if (nounit) {
#line 432 "ztrmm.f"
			if (noconj) {
#line 433 "ztrmm.f"
			    i__1 = k + k * a_dim1;
#line 433 "ztrmm.f"
			    z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				    z__1.i = temp.r * a[i__1].i + temp.i * a[
				    i__1].r;
#line 433 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 434 "ztrmm.f"
			} else {
#line 435 "ztrmm.f"
			    d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 435 "ztrmm.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 435 "ztrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 436 "ztrmm.f"
			}
#line 437 "ztrmm.f"
		    }
#line 438 "ztrmm.f"
		    if (temp.r != 1. || temp.i != 0.) {
#line 439 "ztrmm.f"
			i__1 = *m;
#line 439 "ztrmm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 440 "ztrmm.f"
			    i__2 = i__ + k * b_dim1;
#line 440 "ztrmm.f"
			    i__3 = i__ + k * b_dim1;
#line 440 "ztrmm.f"
			    z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i, 
				    z__1.i = temp.r * b[i__3].i + temp.i * b[
				    i__3].r;
#line 440 "ztrmm.f"
			    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 441 "ztrmm.f"
/* L310: */
#line 441 "ztrmm.f"
			}
#line 442 "ztrmm.f"
		    }
#line 443 "ztrmm.f"
/* L320: */
#line 443 "ztrmm.f"
		}
#line 444 "ztrmm.f"
	    }
#line 445 "ztrmm.f"
	}
#line 446 "ztrmm.f"
    }

#line 448 "ztrmm.f"
    return 0;

/*     End of ZTRMM . */

} /* ztrmm_ */

