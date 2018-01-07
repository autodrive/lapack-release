#line 1 "zhemm.f"
/* zhemm.f -- translated by f2c (version 20100827).
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

#line 1 "zhemm.f"
/* > \brief \b ZHEMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER LDA,LDB,LDC,M,N */
/*       CHARACTER SIDE,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*A*B + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*B*A + beta*C, */
/* > */
/* > where alpha and beta are scalars, A is an hermitian matrix and  B and */
/* > C are m by n matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE  specifies whether  the  hermitian matrix  A */
/* >           appears on the  left or right  in the  operation as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C, */
/* > */
/* >              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C, */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of  the  hermitian  matrix   A  is  to  be */
/* >           referenced as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of the */
/* >                                  hermitian matrix is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of the */
/* >                                  hermitian matrix is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry,  M  specifies the number of rows of the matrix  C. */
/* >           M  must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix C. */
/* >           N  must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is */
/* >           m  when  SIDE = 'L' or 'l'  and is n  otherwise. */
/* >           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of */
/* >           the array  A  must contain the  hermitian matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading m by m upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  hermitian matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  m by m  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  hermitian */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* >           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of */
/* >           the array  A  must contain the  hermitian matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading n by n upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  hermitian matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  n by n  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  hermitian */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* >           Note that the imaginary parts  of the diagonal elements need */
/* >           not be set, they are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then */
/* >           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/* >           least max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension ( LDB, N ) */
/* >           Before entry, the leading  m by n part of the array  B  must */
/* >           contain the matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           max( 1, m ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension ( LDC, N ) */
/* >           Before entry, the leading  m by n  part of the array  C must */
/* >           contain the matrix  C,  except when  beta  is zero, in which */
/* >           case C need not be set on entry. */
/* >           On exit, the array  C  is overwritten by the  m by n updated */
/* >           matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
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
/* Subroutine */ int zhemm_(char *side, char *uplo, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *
	ldc, ftnlen side_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, info;
    static doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     Set NROWA as the number of rows of A. */

#line 234 "zhemm.f"
    /* Parameter adjustments */
#line 234 "zhemm.f"
    a_dim1 = *lda;
#line 234 "zhemm.f"
    a_offset = 1 + a_dim1;
#line 234 "zhemm.f"
    a -= a_offset;
#line 234 "zhemm.f"
    b_dim1 = *ldb;
#line 234 "zhemm.f"
    b_offset = 1 + b_dim1;
#line 234 "zhemm.f"
    b -= b_offset;
#line 234 "zhemm.f"
    c_dim1 = *ldc;
#line 234 "zhemm.f"
    c_offset = 1 + c_dim1;
#line 234 "zhemm.f"
    c__ -= c_offset;
#line 234 "zhemm.f"

#line 234 "zhemm.f"
    /* Function Body */
#line 234 "zhemm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 235 "zhemm.f"
	nrowa = *m;
#line 236 "zhemm.f"
    } else {
#line 237 "zhemm.f"
	nrowa = *n;
#line 238 "zhemm.f"
    }
#line 239 "zhemm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 243 "zhemm.f"
    info = 0;
#line 244 "zhemm.f"
    if (! lsame_(side, "L", (ftnlen)1, (ftnlen)1) && ! lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1)) {
#line 245 "zhemm.f"
	info = 1;
#line 246 "zhemm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 247 "zhemm.f"
	info = 2;
#line 248 "zhemm.f"
    } else if (*m < 0) {
#line 249 "zhemm.f"
	info = 3;
#line 250 "zhemm.f"
    } else if (*n < 0) {
#line 251 "zhemm.f"
	info = 4;
#line 252 "zhemm.f"
    } else if (*lda < max(1,nrowa)) {
#line 253 "zhemm.f"
	info = 7;
#line 254 "zhemm.f"
    } else if (*ldb < max(1,*m)) {
#line 255 "zhemm.f"
	info = 9;
#line 256 "zhemm.f"
    } else if (*ldc < max(1,*m)) {
#line 257 "zhemm.f"
	info = 12;
#line 258 "zhemm.f"
    }
#line 259 "zhemm.f"
    if (info != 0) {
#line 260 "zhemm.f"
	xerbla_("ZHEMM ", &info, (ftnlen)6);
#line 261 "zhemm.f"
	return 0;
#line 262 "zhemm.f"
    }

/*     Quick return if possible. */

#line 266 "zhemm.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
#line 266 "zhemm.f"
	return 0;
#line 266 "zhemm.f"
    }

/*     And when  alpha.eq.zero. */

#line 271 "zhemm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 272 "zhemm.f"
	if (beta->r == 0. && beta->i == 0.) {
#line 273 "zhemm.f"
	    i__1 = *n;
#line 273 "zhemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 274 "zhemm.f"
		i__2 = *m;
#line 274 "zhemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 275 "zhemm.f"
		    i__3 = i__ + j * c_dim1;
#line 275 "zhemm.f"
		    c__[i__3].r = 0., c__[i__3].i = 0.;
#line 276 "zhemm.f"
/* L10: */
#line 276 "zhemm.f"
		}
#line 277 "zhemm.f"
/* L20: */
#line 277 "zhemm.f"
	    }
#line 278 "zhemm.f"
	} else {
#line 279 "zhemm.f"
	    i__1 = *n;
#line 279 "zhemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 280 "zhemm.f"
		i__2 = *m;
#line 280 "zhemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 281 "zhemm.f"
		    i__3 = i__ + j * c_dim1;
#line 281 "zhemm.f"
		    i__4 = i__ + j * c_dim1;
#line 281 "zhemm.f"
		    z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    z__1.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
#line 281 "zhemm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 282 "zhemm.f"
/* L30: */
#line 282 "zhemm.f"
		}
#line 283 "zhemm.f"
/* L40: */
#line 283 "zhemm.f"
	    }
#line 284 "zhemm.f"
	}
#line 285 "zhemm.f"
	return 0;
#line 286 "zhemm.f"
    }

/*     Start the operations. */

#line 290 "zhemm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B + beta*C. */

#line 294 "zhemm.f"
	if (upper) {
#line 295 "zhemm.f"
	    i__1 = *n;
#line 295 "zhemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 296 "zhemm.f"
		i__2 = *m;
#line 296 "zhemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 297 "zhemm.f"
		    i__3 = i__ + j * b_dim1;
#line 297 "zhemm.f"
		    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
			    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3]
			    .r;
#line 297 "zhemm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 298 "zhemm.f"
		    temp2.r = 0., temp2.i = 0.;
#line 299 "zhemm.f"
		    i__3 = i__ - 1;
#line 299 "zhemm.f"
		    for (k = 1; k <= i__3; ++k) {
#line 300 "zhemm.f"
			i__4 = k + j * c_dim1;
#line 300 "zhemm.f"
			i__5 = k + j * c_dim1;
#line 300 "zhemm.f"
			i__6 = k + i__ * a_dim1;
#line 300 "zhemm.f"
			z__2.r = temp1.r * a[i__6].r - temp1.i * a[i__6].i, 
				z__2.i = temp1.r * a[i__6].i + temp1.i * a[
				i__6].r;
#line 300 "zhemm.f"
			z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
				z__2.i;
#line 300 "zhemm.f"
			c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 301 "zhemm.f"
			i__4 = k + j * b_dim1;
#line 301 "zhemm.f"
			d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 301 "zhemm.f"
			z__2.r = b[i__4].r * z__3.r - b[i__4].i * z__3.i, 
				z__2.i = b[i__4].r * z__3.i + b[i__4].i * 
				z__3.r;
#line 301 "zhemm.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 301 "zhemm.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 302 "zhemm.f"
/* L50: */
#line 302 "zhemm.f"
		    }
#line 303 "zhemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 304 "zhemm.f"
			i__3 = i__ + j * c_dim1;
#line 304 "zhemm.f"
			i__4 = i__ + i__ * a_dim1;
#line 304 "zhemm.f"
			d__1 = a[i__4].r;
#line 304 "zhemm.f"
			z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 304 "zhemm.f"
			z__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 304 "zhemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 304 "zhemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 305 "zhemm.f"
		    } else {
#line 306 "zhemm.f"
			i__3 = i__ + j * c_dim1;
#line 306 "zhemm.f"
			i__4 = i__ + j * c_dim1;
#line 306 "zhemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 306 "zhemm.f"
			i__5 = i__ + i__ * a_dim1;
#line 306 "zhemm.f"
			d__1 = a[i__5].r;
#line 306 "zhemm.f"
			z__4.r = d__1 * temp1.r, z__4.i = d__1 * temp1.i;
#line 306 "zhemm.f"
			z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 306 "zhemm.f"
			z__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 306 "zhemm.f"
			z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 306 "zhemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 308 "zhemm.f"
		    }
#line 309 "zhemm.f"
/* L60: */
#line 309 "zhemm.f"
		}
#line 310 "zhemm.f"
/* L70: */
#line 310 "zhemm.f"
	    }
#line 311 "zhemm.f"
	} else {
#line 312 "zhemm.f"
	    i__1 = *n;
#line 312 "zhemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 313 "zhemm.f"
		for (i__ = *m; i__ >= 1; --i__) {
#line 314 "zhemm.f"
		    i__2 = i__ + j * b_dim1;
#line 314 "zhemm.f"
		    z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2].i, 
			    z__1.i = alpha->r * b[i__2].i + alpha->i * b[i__2]
			    .r;
#line 314 "zhemm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 315 "zhemm.f"
		    temp2.r = 0., temp2.i = 0.;
#line 316 "zhemm.f"
		    i__2 = *m;
#line 316 "zhemm.f"
		    for (k = i__ + 1; k <= i__2; ++k) {
#line 317 "zhemm.f"
			i__3 = k + j * c_dim1;
#line 317 "zhemm.f"
			i__4 = k + j * c_dim1;
#line 317 "zhemm.f"
			i__5 = k + i__ * a_dim1;
#line 317 "zhemm.f"
			z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
				z__2.i = temp1.r * a[i__5].i + temp1.i * a[
				i__5].r;
#line 317 "zhemm.f"
			z__1.r = c__[i__4].r + z__2.r, z__1.i = c__[i__4].i + 
				z__2.i;
#line 317 "zhemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 318 "zhemm.f"
			i__3 = k + j * b_dim1;
#line 318 "zhemm.f"
			d_cnjg(&z__3, &a[k + i__ * a_dim1]);
#line 318 "zhemm.f"
			z__2.r = b[i__3].r * z__3.r - b[i__3].i * z__3.i, 
				z__2.i = b[i__3].r * z__3.i + b[i__3].i * 
				z__3.r;
#line 318 "zhemm.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 318 "zhemm.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 319 "zhemm.f"
/* L80: */
#line 319 "zhemm.f"
		    }
#line 320 "zhemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 321 "zhemm.f"
			i__2 = i__ + j * c_dim1;
#line 321 "zhemm.f"
			i__3 = i__ + i__ * a_dim1;
#line 321 "zhemm.f"
			d__1 = a[i__3].r;
#line 321 "zhemm.f"
			z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 321 "zhemm.f"
			z__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 321 "zhemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 321 "zhemm.f"
			c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 322 "zhemm.f"
		    } else {
#line 323 "zhemm.f"
			i__2 = i__ + j * c_dim1;
#line 323 "zhemm.f"
			i__3 = i__ + j * c_dim1;
#line 323 "zhemm.f"
			z__3.r = beta->r * c__[i__3].r - beta->i * c__[i__3]
				.i, z__3.i = beta->r * c__[i__3].i + beta->i *
				 c__[i__3].r;
#line 323 "zhemm.f"
			i__4 = i__ + i__ * a_dim1;
#line 323 "zhemm.f"
			d__1 = a[i__4].r;
#line 323 "zhemm.f"
			z__4.r = d__1 * temp1.r, z__4.i = d__1 * temp1.i;
#line 323 "zhemm.f"
			z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 323 "zhemm.f"
			z__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 323 "zhemm.f"
			z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 323 "zhemm.f"
			c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 325 "zhemm.f"
		    }
#line 326 "zhemm.f"
/* L90: */
#line 326 "zhemm.f"
		}
#line 327 "zhemm.f"
/* L100: */
#line 327 "zhemm.f"
	    }
#line 328 "zhemm.f"
	}
#line 329 "zhemm.f"
    } else {

/*        Form  C := alpha*B*A + beta*C. */

#line 333 "zhemm.f"
	i__1 = *n;
#line 333 "zhemm.f"
	for (j = 1; j <= i__1; ++j) {
#line 334 "zhemm.f"
	    i__2 = j + j * a_dim1;
#line 334 "zhemm.f"
	    d__1 = a[i__2].r;
#line 334 "zhemm.f"
	    z__1.r = d__1 * alpha->r, z__1.i = d__1 * alpha->i;
#line 334 "zhemm.f"
	    temp1.r = z__1.r, temp1.i = z__1.i;
#line 335 "zhemm.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 336 "zhemm.f"
		i__2 = *m;
#line 336 "zhemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 337 "zhemm.f"
		    i__3 = i__ + j * c_dim1;
#line 337 "zhemm.f"
		    i__4 = i__ + j * b_dim1;
#line 337 "zhemm.f"
		    z__1.r = temp1.r * b[i__4].r - temp1.i * b[i__4].i, 
			    z__1.i = temp1.r * b[i__4].i + temp1.i * b[i__4]
			    .r;
#line 337 "zhemm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 338 "zhemm.f"
/* L110: */
#line 338 "zhemm.f"
		}
#line 339 "zhemm.f"
	    } else {
#line 340 "zhemm.f"
		i__2 = *m;
#line 340 "zhemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 341 "zhemm.f"
		    i__3 = i__ + j * c_dim1;
#line 341 "zhemm.f"
		    i__4 = i__ + j * c_dim1;
#line 341 "zhemm.f"
		    z__2.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    z__2.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
#line 341 "zhemm.f"
		    i__5 = i__ + j * b_dim1;
#line 341 "zhemm.f"
		    z__3.r = temp1.r * b[i__5].r - temp1.i * b[i__5].i, 
			    z__3.i = temp1.r * b[i__5].i + temp1.i * b[i__5]
			    .r;
#line 341 "zhemm.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 341 "zhemm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 342 "zhemm.f"
/* L120: */
#line 342 "zhemm.f"
		}
#line 343 "zhemm.f"
	    }
#line 344 "zhemm.f"
	    i__2 = j - 1;
#line 344 "zhemm.f"
	    for (k = 1; k <= i__2; ++k) {
#line 345 "zhemm.f"
		if (upper) {
#line 346 "zhemm.f"
		    i__3 = k + j * a_dim1;
#line 346 "zhemm.f"
		    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
			    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3]
			    .r;
#line 346 "zhemm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 347 "zhemm.f"
		} else {
#line 348 "zhemm.f"
		    d_cnjg(&z__2, &a[j + k * a_dim1]);
#line 348 "zhemm.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 348 "zhemm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 349 "zhemm.f"
		}
#line 350 "zhemm.f"
		i__3 = *m;
#line 350 "zhemm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 351 "zhemm.f"
		    i__4 = i__ + j * c_dim1;
#line 351 "zhemm.f"
		    i__5 = i__ + j * c_dim1;
#line 351 "zhemm.f"
		    i__6 = i__ + k * b_dim1;
#line 351 "zhemm.f"
		    z__2.r = temp1.r * b[i__6].r - temp1.i * b[i__6].i, 
			    z__2.i = temp1.r * b[i__6].i + temp1.i * b[i__6]
			    .r;
#line 351 "zhemm.f"
		    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
			    z__2.i;
#line 351 "zhemm.f"
		    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 352 "zhemm.f"
/* L130: */
#line 352 "zhemm.f"
		}
#line 353 "zhemm.f"
/* L140: */
#line 353 "zhemm.f"
	    }
#line 354 "zhemm.f"
	    i__2 = *n;
#line 354 "zhemm.f"
	    for (k = j + 1; k <= i__2; ++k) {
#line 355 "zhemm.f"
		if (upper) {
#line 356 "zhemm.f"
		    d_cnjg(&z__2, &a[j + k * a_dim1]);
#line 356 "zhemm.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 356 "zhemm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 357 "zhemm.f"
		} else {
#line 358 "zhemm.f"
		    i__3 = k + j * a_dim1;
#line 358 "zhemm.f"
		    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
			    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3]
			    .r;
#line 358 "zhemm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 359 "zhemm.f"
		}
#line 360 "zhemm.f"
		i__3 = *m;
#line 360 "zhemm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 361 "zhemm.f"
		    i__4 = i__ + j * c_dim1;
#line 361 "zhemm.f"
		    i__5 = i__ + j * c_dim1;
#line 361 "zhemm.f"
		    i__6 = i__ + k * b_dim1;
#line 361 "zhemm.f"
		    z__2.r = temp1.r * b[i__6].r - temp1.i * b[i__6].i, 
			    z__2.i = temp1.r * b[i__6].i + temp1.i * b[i__6]
			    .r;
#line 361 "zhemm.f"
		    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
			    z__2.i;
#line 361 "zhemm.f"
		    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 362 "zhemm.f"
/* L150: */
#line 362 "zhemm.f"
		}
#line 363 "zhemm.f"
/* L160: */
#line 363 "zhemm.f"
	    }
#line 364 "zhemm.f"
/* L170: */
#line 364 "zhemm.f"
	}
#line 365 "zhemm.f"
    }

#line 367 "zhemm.f"
    return 0;

/*     End of ZHEMM . */

} /* zhemm_ */

