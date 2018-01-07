#line 1 "ctrmm.f"
/* ctrmm.f -- translated by f2c (version 20100827).
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

#line 1 "ctrmm.f"
/* > \brief \b CTRMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

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
/* > CTRMM  performs one of the matrix-matrix operations */
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
/* >          ALPHA is COMPLEX */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, k ), where k is m */
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
/* >          B is COMPLEX array of DIMENSION ( LDB, n ). */
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
/* Subroutine */ int ctrmm_(char *side, char *uplo, char *transa, char *diag, 
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lside;
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

#line 220 "ctrmm.f"
    /* Parameter adjustments */
#line 220 "ctrmm.f"
    a_dim1 = *lda;
#line 220 "ctrmm.f"
    a_offset = 1 + a_dim1;
#line 220 "ctrmm.f"
    a -= a_offset;
#line 220 "ctrmm.f"
    b_dim1 = *ldb;
#line 220 "ctrmm.f"
    b_offset = 1 + b_dim1;
#line 220 "ctrmm.f"
    b -= b_offset;
#line 220 "ctrmm.f"

#line 220 "ctrmm.f"
    /* Function Body */
#line 220 "ctrmm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 221 "ctrmm.f"
    if (lside) {
#line 222 "ctrmm.f"
	nrowa = *m;
#line 223 "ctrmm.f"
    } else {
#line 224 "ctrmm.f"
	nrowa = *n;
#line 225 "ctrmm.f"
    }
#line 226 "ctrmm.f"
    noconj = lsame_(transa, "T", (ftnlen)1, (ftnlen)1);
#line 227 "ctrmm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 228 "ctrmm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 230 "ctrmm.f"
    info = 0;
#line 231 "ctrmm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 232 "ctrmm.f"
	info = 1;
#line 233 "ctrmm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 234 "ctrmm.f"
	info = 2;
#line 235 "ctrmm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 238 "ctrmm.f"
	info = 3;
#line 239 "ctrmm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 240 "ctrmm.f"
	info = 4;
#line 241 "ctrmm.f"
    } else if (*m < 0) {
#line 242 "ctrmm.f"
	info = 5;
#line 243 "ctrmm.f"
    } else if (*n < 0) {
#line 244 "ctrmm.f"
	info = 6;
#line 245 "ctrmm.f"
    } else if (*lda < max(1,nrowa)) {
#line 246 "ctrmm.f"
	info = 9;
#line 247 "ctrmm.f"
    } else if (*ldb < max(1,*m)) {
#line 248 "ctrmm.f"
	info = 11;
#line 249 "ctrmm.f"
    }
#line 250 "ctrmm.f"
    if (info != 0) {
#line 251 "ctrmm.f"
	xerbla_("CTRMM ", &info, (ftnlen)6);
#line 252 "ctrmm.f"
	return 0;
#line 253 "ctrmm.f"
    }

/*     Quick return if possible. */

#line 257 "ctrmm.f"
    if (*m == 0 || *n == 0) {
#line 257 "ctrmm.f"
	return 0;
#line 257 "ctrmm.f"
    }

/*     And when  alpha.eq.zero. */

#line 261 "ctrmm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 262 "ctrmm.f"
	i__1 = *n;
#line 262 "ctrmm.f"
	for (j = 1; j <= i__1; ++j) {
#line 263 "ctrmm.f"
	    i__2 = *m;
#line 263 "ctrmm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 264 "ctrmm.f"
		i__3 = i__ + j * b_dim1;
#line 264 "ctrmm.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 265 "ctrmm.f"
/* L10: */
#line 265 "ctrmm.f"
	    }
#line 266 "ctrmm.f"
/* L20: */
#line 266 "ctrmm.f"
	}
#line 267 "ctrmm.f"
	return 0;
#line 268 "ctrmm.f"
    }

/*     Start the operations. */

#line 272 "ctrmm.f"
    if (lside) {
#line 273 "ctrmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*A*B. */

#line 277 "ctrmm.f"
	    if (upper) {
#line 278 "ctrmm.f"
		i__1 = *n;
#line 278 "ctrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 279 "ctrmm.f"
		    i__2 = *m;
#line 279 "ctrmm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 280 "ctrmm.f"
			i__3 = k + j * b_dim1;
#line 280 "ctrmm.f"
			if (b[i__3].r != 0. || b[i__3].i != 0.) {
#line 281 "ctrmm.f"
			    i__3 = k + j * b_dim1;
#line 281 "ctrmm.f"
			    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3]
				    .i, z__1.i = alpha->r * b[i__3].i + 
				    alpha->i * b[i__3].r;
#line 281 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 282 "ctrmm.f"
			    i__3 = k - 1;
#line 282 "ctrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 283 "ctrmm.f"
				i__4 = i__ + j * b_dim1;
#line 283 "ctrmm.f"
				i__5 = i__ + j * b_dim1;
#line 283 "ctrmm.f"
				i__6 = i__ + k * a_dim1;
#line 283 "ctrmm.f"
				z__2.r = temp.r * a[i__6].r - temp.i * a[i__6]
					.i, z__2.i = temp.r * a[i__6].i + 
					temp.i * a[i__6].r;
#line 283 "ctrmm.f"
				z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5]
					.i + z__2.i;
#line 283 "ctrmm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 284 "ctrmm.f"
/* L30: */
#line 284 "ctrmm.f"
			    }
#line 285 "ctrmm.f"
			    if (nounit) {
#line 285 "ctrmm.f"
				i__3 = k + k * a_dim1;
#line 285 "ctrmm.f"
				z__1.r = temp.r * a[i__3].r - temp.i * a[i__3]
					.i, z__1.i = temp.r * a[i__3].i + 
					temp.i * a[i__3].r;
#line 285 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 285 "ctrmm.f"
			    }
#line 286 "ctrmm.f"
			    i__3 = k + j * b_dim1;
#line 286 "ctrmm.f"
			    b[i__3].r = temp.r, b[i__3].i = temp.i;
#line 287 "ctrmm.f"
			}
#line 288 "ctrmm.f"
/* L40: */
#line 288 "ctrmm.f"
		    }
#line 289 "ctrmm.f"
/* L50: */
#line 289 "ctrmm.f"
		}
#line 290 "ctrmm.f"
	    } else {
#line 291 "ctrmm.f"
		i__1 = *n;
#line 291 "ctrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "ctrmm.f"
		    for (k = *m; k >= 1; --k) {
#line 293 "ctrmm.f"
			i__2 = k + j * b_dim1;
#line 293 "ctrmm.f"
			if (b[i__2].r != 0. || b[i__2].i != 0.) {
#line 294 "ctrmm.f"
			    i__2 = k + j * b_dim1;
#line 294 "ctrmm.f"
			    z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2]
				    .i, z__1.i = alpha->r * b[i__2].i + 
				    alpha->i * b[i__2].r;
#line 294 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 295 "ctrmm.f"
			    i__2 = k + j * b_dim1;
#line 295 "ctrmm.f"
			    b[i__2].r = temp.r, b[i__2].i = temp.i;
#line 296 "ctrmm.f"
			    if (nounit) {
#line 296 "ctrmm.f"
				i__2 = k + j * b_dim1;
#line 296 "ctrmm.f"
				i__3 = k + j * b_dim1;
#line 296 "ctrmm.f"
				i__4 = k + k * a_dim1;
#line 296 "ctrmm.f"
				z__1.r = b[i__3].r * a[i__4].r - b[i__3].i * 
					a[i__4].i, z__1.i = b[i__3].r * a[
					i__4].i + b[i__3].i * a[i__4].r;
#line 296 "ctrmm.f"
				b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 296 "ctrmm.f"
			    }
#line 297 "ctrmm.f"
			    i__2 = *m;
#line 297 "ctrmm.f"
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 298 "ctrmm.f"
				i__3 = i__ + j * b_dim1;
#line 298 "ctrmm.f"
				i__4 = i__ + j * b_dim1;
#line 298 "ctrmm.f"
				i__5 = i__ + k * a_dim1;
#line 298 "ctrmm.f"
				z__2.r = temp.r * a[i__5].r - temp.i * a[i__5]
					.i, z__2.i = temp.r * a[i__5].i + 
					temp.i * a[i__5].r;
#line 298 "ctrmm.f"
				z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4]
					.i + z__2.i;
#line 298 "ctrmm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 299 "ctrmm.f"
/* L60: */
#line 299 "ctrmm.f"
			    }
#line 300 "ctrmm.f"
			}
#line 301 "ctrmm.f"
/* L70: */
#line 301 "ctrmm.f"
		    }
#line 302 "ctrmm.f"
/* L80: */
#line 302 "ctrmm.f"
		}
#line 303 "ctrmm.f"
	    }
#line 304 "ctrmm.f"
	} else {

/*           Form  B := alpha*A**T*B   or   B := alpha*A**H*B. */

#line 308 "ctrmm.f"
	    if (upper) {
#line 309 "ctrmm.f"
		i__1 = *n;
#line 309 "ctrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 310 "ctrmm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 311 "ctrmm.f"
			i__2 = i__ + j * b_dim1;
#line 311 "ctrmm.f"
			temp.r = b[i__2].r, temp.i = b[i__2].i;
#line 312 "ctrmm.f"
			if (noconj) {
#line 313 "ctrmm.f"
			    if (nounit) {
#line 313 "ctrmm.f"
				i__2 = i__ + i__ * a_dim1;
#line 313 "ctrmm.f"
				z__1.r = temp.r * a[i__2].r - temp.i * a[i__2]
					.i, z__1.i = temp.r * a[i__2].i + 
					temp.i * a[i__2].r;
#line 313 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 313 "ctrmm.f"
			    }
#line 314 "ctrmm.f"
			    i__2 = i__ - 1;
#line 314 "ctrmm.f"
			    for (k = 1; k <= i__2; ++k) {
#line 315 "ctrmm.f"
				i__3 = k + i__ * a_dim1;
#line 315 "ctrmm.f"
				i__4 = k + j * b_dim1;
#line 315 "ctrmm.f"
				z__2.r = a[i__3].r * b[i__4].r - a[i__3].i * 
					b[i__4].i, z__2.i = a[i__3].r * b[
					i__4].i + a[i__3].i * b[i__4].r;
#line 315 "ctrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 315 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 316 "ctrmm.f"
/* L90: */
#line 316 "ctrmm.f"
			    }
#line 317 "ctrmm.f"
			} else {
#line 318 "ctrmm.f"
			    if (nounit) {
#line 318 "ctrmm.f"
				d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
#line 318 "ctrmm.f"
				z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
					z__1.i = temp.r * z__2.i + temp.i * 
					z__2.r;
#line 318 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 318 "ctrmm.f"
			    }
#line 319 "ctrmm.f"
			    i__2 = i__ - 1;
#line 319 "ctrmm.f"
			    for (k = 1; k <= i__2; ++k) {
#line 320 "ctrmm.f"
				d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 320 "ctrmm.f"
				i__3 = k + j * b_dim1;
#line 320 "ctrmm.f"
				z__2.r = z__3.r * b[i__3].r - z__3.i * b[i__3]
					.i, z__2.i = z__3.r * b[i__3].i + 
					z__3.i * b[i__3].r;
#line 320 "ctrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 320 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 321 "ctrmm.f"
/* L100: */
#line 321 "ctrmm.f"
			    }
#line 322 "ctrmm.f"
			}
#line 323 "ctrmm.f"
			i__2 = i__ + j * b_dim1;
#line 323 "ctrmm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 323 "ctrmm.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 324 "ctrmm.f"
/* L110: */
#line 324 "ctrmm.f"
		    }
#line 325 "ctrmm.f"
/* L120: */
#line 325 "ctrmm.f"
		}
#line 326 "ctrmm.f"
	    } else {
#line 327 "ctrmm.f"
		i__1 = *n;
#line 327 "ctrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 328 "ctrmm.f"
		    i__2 = *m;
#line 328 "ctrmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 329 "ctrmm.f"
			i__3 = i__ + j * b_dim1;
#line 329 "ctrmm.f"
			temp.r = b[i__3].r, temp.i = b[i__3].i;
#line 330 "ctrmm.f"
			if (noconj) {
#line 331 "ctrmm.f"
			    if (nounit) {
#line 331 "ctrmm.f"
				i__3 = i__ + i__ * a_dim1;
#line 331 "ctrmm.f"
				z__1.r = temp.r * a[i__3].r - temp.i * a[i__3]
					.i, z__1.i = temp.r * a[i__3].i + 
					temp.i * a[i__3].r;
#line 331 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 331 "ctrmm.f"
			    }
#line 332 "ctrmm.f"
			    i__3 = *m;
#line 332 "ctrmm.f"
			    for (k = i__ + 1; k <= i__3; ++k) {
#line 333 "ctrmm.f"
				i__4 = k + i__ * a_dim1;
#line 333 "ctrmm.f"
				i__5 = k + j * b_dim1;
#line 333 "ctrmm.f"
				z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * 
					b[i__5].i, z__2.i = a[i__4].r * b[
					i__5].i + a[i__4].i * b[i__5].r;
#line 333 "ctrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 333 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 334 "ctrmm.f"
/* L130: */
#line 334 "ctrmm.f"
			    }
#line 335 "ctrmm.f"
			} else {
#line 336 "ctrmm.f"
			    if (nounit) {
#line 336 "ctrmm.f"
				d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
#line 336 "ctrmm.f"
				z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
					z__1.i = temp.r * z__2.i + temp.i * 
					z__2.r;
#line 336 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 336 "ctrmm.f"
			    }
#line 337 "ctrmm.f"
			    i__3 = *m;
#line 337 "ctrmm.f"
			    for (k = i__ + 1; k <= i__3; ++k) {
#line 338 "ctrmm.f"
				d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 338 "ctrmm.f"
				i__4 = k + j * b_dim1;
#line 338 "ctrmm.f"
				z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4]
					.i, z__2.i = z__3.r * b[i__4].i + 
					z__3.i * b[i__4].r;
#line 338 "ctrmm.f"
				z__1.r = temp.r + z__2.r, z__1.i = temp.i + 
					z__2.i;
#line 338 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 339 "ctrmm.f"
/* L140: */
#line 339 "ctrmm.f"
			    }
#line 340 "ctrmm.f"
			}
#line 341 "ctrmm.f"
			i__3 = i__ + j * b_dim1;
#line 341 "ctrmm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 341 "ctrmm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 342 "ctrmm.f"
/* L150: */
#line 342 "ctrmm.f"
		    }
#line 343 "ctrmm.f"
/* L160: */
#line 343 "ctrmm.f"
		}
#line 344 "ctrmm.f"
	    }
#line 345 "ctrmm.f"
	}
#line 346 "ctrmm.f"
    } else {
#line 347 "ctrmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*A. */

#line 351 "ctrmm.f"
	    if (upper) {
#line 352 "ctrmm.f"
		for (j = *n; j >= 1; --j) {
#line 353 "ctrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 354 "ctrmm.f"
		    if (nounit) {
#line 354 "ctrmm.f"
			i__1 = j + j * a_dim1;
#line 354 "ctrmm.f"
			z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				z__1.i = temp.r * a[i__1].i + temp.i * a[i__1]
				.r;
#line 354 "ctrmm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 354 "ctrmm.f"
		    }
#line 355 "ctrmm.f"
		    i__1 = *m;
#line 355 "ctrmm.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 356 "ctrmm.f"
			i__2 = i__ + j * b_dim1;
#line 356 "ctrmm.f"
			i__3 = i__ + j * b_dim1;
#line 356 "ctrmm.f"
			z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i, 
				z__1.i = temp.r * b[i__3].i + temp.i * b[i__3]
				.r;
#line 356 "ctrmm.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 357 "ctrmm.f"
/* L170: */
#line 357 "ctrmm.f"
		    }
#line 358 "ctrmm.f"
		    i__1 = j - 1;
#line 358 "ctrmm.f"
		    for (k = 1; k <= i__1; ++k) {
#line 359 "ctrmm.f"
			i__2 = k + j * a_dim1;
#line 359 "ctrmm.f"
			if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 360 "ctrmm.f"
			    i__2 = k + j * a_dim1;
#line 360 "ctrmm.f"
			    z__1.r = alpha->r * a[i__2].r - alpha->i * a[i__2]
				    .i, z__1.i = alpha->r * a[i__2].i + 
				    alpha->i * a[i__2].r;
#line 360 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 361 "ctrmm.f"
			    i__2 = *m;
#line 361 "ctrmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 362 "ctrmm.f"
				i__3 = i__ + j * b_dim1;
#line 362 "ctrmm.f"
				i__4 = i__ + j * b_dim1;
#line 362 "ctrmm.f"
				i__5 = i__ + k * b_dim1;
#line 362 "ctrmm.f"
				z__2.r = temp.r * b[i__5].r - temp.i * b[i__5]
					.i, z__2.i = temp.r * b[i__5].i + 
					temp.i * b[i__5].r;
#line 362 "ctrmm.f"
				z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4]
					.i + z__2.i;
#line 362 "ctrmm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 363 "ctrmm.f"
/* L180: */
#line 363 "ctrmm.f"
			    }
#line 364 "ctrmm.f"
			}
#line 365 "ctrmm.f"
/* L190: */
#line 365 "ctrmm.f"
		    }
#line 366 "ctrmm.f"
/* L200: */
#line 366 "ctrmm.f"
		}
#line 367 "ctrmm.f"
	    } else {
#line 368 "ctrmm.f"
		i__1 = *n;
#line 368 "ctrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 369 "ctrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 370 "ctrmm.f"
		    if (nounit) {
#line 370 "ctrmm.f"
			i__2 = j + j * a_dim1;
#line 370 "ctrmm.f"
			z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				z__1.i = temp.r * a[i__2].i + temp.i * a[i__2]
				.r;
#line 370 "ctrmm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 370 "ctrmm.f"
		    }
#line 371 "ctrmm.f"
		    i__2 = *m;
#line 371 "ctrmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 372 "ctrmm.f"
			i__3 = i__ + j * b_dim1;
#line 372 "ctrmm.f"
			i__4 = i__ + j * b_dim1;
#line 372 "ctrmm.f"
			z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i, 
				z__1.i = temp.r * b[i__4].i + temp.i * b[i__4]
				.r;
#line 372 "ctrmm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 373 "ctrmm.f"
/* L210: */
#line 373 "ctrmm.f"
		    }
#line 374 "ctrmm.f"
		    i__2 = *n;
#line 374 "ctrmm.f"
		    for (k = j + 1; k <= i__2; ++k) {
#line 375 "ctrmm.f"
			i__3 = k + j * a_dim1;
#line 375 "ctrmm.f"
			if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 376 "ctrmm.f"
			    i__3 = k + j * a_dim1;
#line 376 "ctrmm.f"
			    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3]
				    .i, z__1.i = alpha->r * a[i__3].i + 
				    alpha->i * a[i__3].r;
#line 376 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 377 "ctrmm.f"
			    i__3 = *m;
#line 377 "ctrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 378 "ctrmm.f"
				i__4 = i__ + j * b_dim1;
#line 378 "ctrmm.f"
				i__5 = i__ + j * b_dim1;
#line 378 "ctrmm.f"
				i__6 = i__ + k * b_dim1;
#line 378 "ctrmm.f"
				z__2.r = temp.r * b[i__6].r - temp.i * b[i__6]
					.i, z__2.i = temp.r * b[i__6].i + 
					temp.i * b[i__6].r;
#line 378 "ctrmm.f"
				z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5]
					.i + z__2.i;
#line 378 "ctrmm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 379 "ctrmm.f"
/* L220: */
#line 379 "ctrmm.f"
			    }
#line 380 "ctrmm.f"
			}
#line 381 "ctrmm.f"
/* L230: */
#line 381 "ctrmm.f"
		    }
#line 382 "ctrmm.f"
/* L240: */
#line 382 "ctrmm.f"
		}
#line 383 "ctrmm.f"
	    }
#line 384 "ctrmm.f"
	} else {

/*           Form  B := alpha*B*A**T   or   B := alpha*B*A**H. */

#line 388 "ctrmm.f"
	    if (upper) {
#line 389 "ctrmm.f"
		i__1 = *n;
#line 389 "ctrmm.f"
		for (k = 1; k <= i__1; ++k) {
#line 390 "ctrmm.f"
		    i__2 = k - 1;
#line 390 "ctrmm.f"
		    for (j = 1; j <= i__2; ++j) {
#line 391 "ctrmm.f"
			i__3 = j + k * a_dim1;
#line 391 "ctrmm.f"
			if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 392 "ctrmm.f"
			    if (noconj) {
#line 393 "ctrmm.f"
				i__3 = j + k * a_dim1;
#line 393 "ctrmm.f"
				z__1.r = alpha->r * a[i__3].r - alpha->i * a[
					i__3].i, z__1.i = alpha->r * a[i__3]
					.i + alpha->i * a[i__3].r;
#line 393 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 394 "ctrmm.f"
			    } else {
#line 395 "ctrmm.f"
				d_cnjg(&z__2, &a[j + k * a_dim1]);
#line 395 "ctrmm.f"
				z__1.r = alpha->r * z__2.r - alpha->i * 
					z__2.i, z__1.i = alpha->r * z__2.i + 
					alpha->i * z__2.r;
#line 395 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 396 "ctrmm.f"
			    }
#line 397 "ctrmm.f"
			    i__3 = *m;
#line 397 "ctrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 398 "ctrmm.f"
				i__4 = i__ + j * b_dim1;
#line 398 "ctrmm.f"
				i__5 = i__ + j * b_dim1;
#line 398 "ctrmm.f"
				i__6 = i__ + k * b_dim1;
#line 398 "ctrmm.f"
				z__2.r = temp.r * b[i__6].r - temp.i * b[i__6]
					.i, z__2.i = temp.r * b[i__6].i + 
					temp.i * b[i__6].r;
#line 398 "ctrmm.f"
				z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5]
					.i + z__2.i;
#line 398 "ctrmm.f"
				b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 399 "ctrmm.f"
/* L250: */
#line 399 "ctrmm.f"
			    }
#line 400 "ctrmm.f"
			}
#line 401 "ctrmm.f"
/* L260: */
#line 401 "ctrmm.f"
		    }
#line 402 "ctrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 403 "ctrmm.f"
		    if (nounit) {
#line 404 "ctrmm.f"
			if (noconj) {
#line 405 "ctrmm.f"
			    i__2 = k + k * a_dim1;
#line 405 "ctrmm.f"
			    z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i, 
				    z__1.i = temp.r * a[i__2].i + temp.i * a[
				    i__2].r;
#line 405 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 406 "ctrmm.f"
			} else {
#line 407 "ctrmm.f"
			    d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 407 "ctrmm.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 407 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 408 "ctrmm.f"
			}
#line 409 "ctrmm.f"
		    }
#line 410 "ctrmm.f"
		    if (temp.r != 1. || temp.i != 0.) {
#line 411 "ctrmm.f"
			i__2 = *m;
#line 411 "ctrmm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 412 "ctrmm.f"
			    i__3 = i__ + k * b_dim1;
#line 412 "ctrmm.f"
			    i__4 = i__ + k * b_dim1;
#line 412 "ctrmm.f"
			    z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i, 
				    z__1.i = temp.r * b[i__4].i + temp.i * b[
				    i__4].r;
#line 412 "ctrmm.f"
			    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 413 "ctrmm.f"
/* L270: */
#line 413 "ctrmm.f"
			}
#line 414 "ctrmm.f"
		    }
#line 415 "ctrmm.f"
/* L280: */
#line 415 "ctrmm.f"
		}
#line 416 "ctrmm.f"
	    } else {
#line 417 "ctrmm.f"
		for (k = *n; k >= 1; --k) {
#line 418 "ctrmm.f"
		    i__1 = *n;
#line 418 "ctrmm.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 419 "ctrmm.f"
			i__2 = j + k * a_dim1;
#line 419 "ctrmm.f"
			if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 420 "ctrmm.f"
			    if (noconj) {
#line 421 "ctrmm.f"
				i__2 = j + k * a_dim1;
#line 421 "ctrmm.f"
				z__1.r = alpha->r * a[i__2].r - alpha->i * a[
					i__2].i, z__1.i = alpha->r * a[i__2]
					.i + alpha->i * a[i__2].r;
#line 421 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 422 "ctrmm.f"
			    } else {
#line 423 "ctrmm.f"
				d_cnjg(&z__2, &a[j + k * a_dim1]);
#line 423 "ctrmm.f"
				z__1.r = alpha->r * z__2.r - alpha->i * 
					z__2.i, z__1.i = alpha->r * z__2.i + 
					alpha->i * z__2.r;
#line 423 "ctrmm.f"
				temp.r = z__1.r, temp.i = z__1.i;
#line 424 "ctrmm.f"
			    }
#line 425 "ctrmm.f"
			    i__2 = *m;
#line 425 "ctrmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 426 "ctrmm.f"
				i__3 = i__ + j * b_dim1;
#line 426 "ctrmm.f"
				i__4 = i__ + j * b_dim1;
#line 426 "ctrmm.f"
				i__5 = i__ + k * b_dim1;
#line 426 "ctrmm.f"
				z__2.r = temp.r * b[i__5].r - temp.i * b[i__5]
					.i, z__2.i = temp.r * b[i__5].i + 
					temp.i * b[i__5].r;
#line 426 "ctrmm.f"
				z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4]
					.i + z__2.i;
#line 426 "ctrmm.f"
				b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 427 "ctrmm.f"
/* L290: */
#line 427 "ctrmm.f"
			    }
#line 428 "ctrmm.f"
			}
#line 429 "ctrmm.f"
/* L300: */
#line 429 "ctrmm.f"
		    }
#line 430 "ctrmm.f"
		    temp.r = alpha->r, temp.i = alpha->i;
#line 431 "ctrmm.f"
		    if (nounit) {
#line 432 "ctrmm.f"
			if (noconj) {
#line 433 "ctrmm.f"
			    i__1 = k + k * a_dim1;
#line 433 "ctrmm.f"
			    z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i, 
				    z__1.i = temp.r * a[i__1].i + temp.i * a[
				    i__1].r;
#line 433 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 434 "ctrmm.f"
			} else {
#line 435 "ctrmm.f"
			    d_cnjg(&z__2, &a[k + k * a_dim1]);
#line 435 "ctrmm.f"
			    z__1.r = temp.r * z__2.r - temp.i * z__2.i, 
				    z__1.i = temp.r * z__2.i + temp.i * 
				    z__2.r;
#line 435 "ctrmm.f"
			    temp.r = z__1.r, temp.i = z__1.i;
#line 436 "ctrmm.f"
			}
#line 437 "ctrmm.f"
		    }
#line 438 "ctrmm.f"
		    if (temp.r != 1. || temp.i != 0.) {
#line 439 "ctrmm.f"
			i__1 = *m;
#line 439 "ctrmm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 440 "ctrmm.f"
			    i__2 = i__ + k * b_dim1;
#line 440 "ctrmm.f"
			    i__3 = i__ + k * b_dim1;
#line 440 "ctrmm.f"
			    z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i, 
				    z__1.i = temp.r * b[i__3].i + temp.i * b[
				    i__3].r;
#line 440 "ctrmm.f"
			    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 441 "ctrmm.f"
/* L310: */
#line 441 "ctrmm.f"
			}
#line 442 "ctrmm.f"
		    }
#line 443 "ctrmm.f"
/* L320: */
#line 443 "ctrmm.f"
		}
#line 444 "ctrmm.f"
	    }
#line 445 "ctrmm.f"
	}
#line 446 "ctrmm.f"
    }

#line 448 "ctrmm.f"
    return 0;

/*     End of CTRMM . */

} /* ctrmm_ */

