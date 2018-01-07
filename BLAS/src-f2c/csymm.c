#line 1 "csymm.f"
/* csymm.f -- translated by f2c (version 20100827).
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

#line 1 "csymm.f"
/* > \brief \b CSYMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER LDA,LDB,LDC,M,N */
/*       CHARACTER SIDE,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*A*B + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*B*A + beta*C, */
/* > */
/* > where  alpha and beta are scalars, A is a symmetric matrix and  B and */
/* > C are m by n matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE  specifies whether  the  symmetric matrix  A */
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
/* >           triangular  part  of  the  symmetric  matrix   A  is  to  be */
/* >           referenced as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of the */
/* >                                  symmetric matrix is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of the */
/* >                                  symmetric matrix is to be referenced. */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is */
/* >           m  when  SIDE = 'L' or 'l'  and is n  otherwise. */
/* >           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of */
/* >           the array  A  must contain the  symmetric matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading m by m upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  symmetric matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  m by m  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  symmetric */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* >           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of */
/* >           the array  A  must contain the  symmetric matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading n by n upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  symmetric matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  n by n  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  symmetric */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
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
/* >          B is COMPLEX array of DIMENSION ( LDB, n ). */
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
/* >          BETA is COMPLEX */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array of DIMENSION ( LDC, n ). */
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
/* Subroutine */ int csymm_(char *side, char *uplo, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *
	ldc, ftnlen side_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

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

#line 232 "csymm.f"
    /* Parameter adjustments */
#line 232 "csymm.f"
    a_dim1 = *lda;
#line 232 "csymm.f"
    a_offset = 1 + a_dim1;
#line 232 "csymm.f"
    a -= a_offset;
#line 232 "csymm.f"
    b_dim1 = *ldb;
#line 232 "csymm.f"
    b_offset = 1 + b_dim1;
#line 232 "csymm.f"
    b -= b_offset;
#line 232 "csymm.f"
    c_dim1 = *ldc;
#line 232 "csymm.f"
    c_offset = 1 + c_dim1;
#line 232 "csymm.f"
    c__ -= c_offset;
#line 232 "csymm.f"

#line 232 "csymm.f"
    /* Function Body */
#line 232 "csymm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 233 "csymm.f"
	nrowa = *m;
#line 234 "csymm.f"
    } else {
#line 235 "csymm.f"
	nrowa = *n;
#line 236 "csymm.f"
    }
#line 237 "csymm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 241 "csymm.f"
    info = 0;
#line 242 "csymm.f"
    if (! lsame_(side, "L", (ftnlen)1, (ftnlen)1) && ! lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1)) {
#line 243 "csymm.f"
	info = 1;
#line 244 "csymm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 245 "csymm.f"
	info = 2;
#line 246 "csymm.f"
    } else if (*m < 0) {
#line 247 "csymm.f"
	info = 3;
#line 248 "csymm.f"
    } else if (*n < 0) {
#line 249 "csymm.f"
	info = 4;
#line 250 "csymm.f"
    } else if (*lda < max(1,nrowa)) {
#line 251 "csymm.f"
	info = 7;
#line 252 "csymm.f"
    } else if (*ldb < max(1,*m)) {
#line 253 "csymm.f"
	info = 9;
#line 254 "csymm.f"
    } else if (*ldc < max(1,*m)) {
#line 255 "csymm.f"
	info = 12;
#line 256 "csymm.f"
    }
#line 257 "csymm.f"
    if (info != 0) {
#line 258 "csymm.f"
	xerbla_("CSYMM ", &info, (ftnlen)6);
#line 259 "csymm.f"
	return 0;
#line 260 "csymm.f"
    }

/*     Quick return if possible. */

#line 264 "csymm.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
#line 264 "csymm.f"
	return 0;
#line 264 "csymm.f"
    }

/*     And when  alpha.eq.zero. */

#line 269 "csymm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 270 "csymm.f"
	if (beta->r == 0. && beta->i == 0.) {
#line 271 "csymm.f"
	    i__1 = *n;
#line 271 "csymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 272 "csymm.f"
		i__2 = *m;
#line 272 "csymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 273 "csymm.f"
		    i__3 = i__ + j * c_dim1;
#line 273 "csymm.f"
		    c__[i__3].r = 0., c__[i__3].i = 0.;
#line 274 "csymm.f"
/* L10: */
#line 274 "csymm.f"
		}
#line 275 "csymm.f"
/* L20: */
#line 275 "csymm.f"
	    }
#line 276 "csymm.f"
	} else {
#line 277 "csymm.f"
	    i__1 = *n;
#line 277 "csymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 278 "csymm.f"
		i__2 = *m;
#line 278 "csymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 279 "csymm.f"
		    i__3 = i__ + j * c_dim1;
#line 279 "csymm.f"
		    i__4 = i__ + j * c_dim1;
#line 279 "csymm.f"
		    z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    z__1.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
#line 279 "csymm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 280 "csymm.f"
/* L30: */
#line 280 "csymm.f"
		}
#line 281 "csymm.f"
/* L40: */
#line 281 "csymm.f"
	    }
#line 282 "csymm.f"
	}
#line 283 "csymm.f"
	return 0;
#line 284 "csymm.f"
    }

/*     Start the operations. */

#line 288 "csymm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B + beta*C. */

#line 292 "csymm.f"
	if (upper) {
#line 293 "csymm.f"
	    i__1 = *n;
#line 293 "csymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 294 "csymm.f"
		i__2 = *m;
#line 294 "csymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "csymm.f"
		    i__3 = i__ + j * b_dim1;
#line 295 "csymm.f"
		    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
			    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3]
			    .r;
#line 295 "csymm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 296 "csymm.f"
		    temp2.r = 0., temp2.i = 0.;
#line 297 "csymm.f"
		    i__3 = i__ - 1;
#line 297 "csymm.f"
		    for (k = 1; k <= i__3; ++k) {
#line 298 "csymm.f"
			i__4 = k + j * c_dim1;
#line 298 "csymm.f"
			i__5 = k + j * c_dim1;
#line 298 "csymm.f"
			i__6 = k + i__ * a_dim1;
#line 298 "csymm.f"
			z__2.r = temp1.r * a[i__6].r - temp1.i * a[i__6].i, 
				z__2.i = temp1.r * a[i__6].i + temp1.i * a[
				i__6].r;
#line 298 "csymm.f"
			z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
				z__2.i;
#line 298 "csymm.f"
			c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 299 "csymm.f"
			i__4 = k + j * b_dim1;
#line 299 "csymm.f"
			i__5 = k + i__ * a_dim1;
#line 299 "csymm.f"
			z__2.r = b[i__4].r * a[i__5].r - b[i__4].i * a[i__5]
				.i, z__2.i = b[i__4].r * a[i__5].i + b[i__4]
				.i * a[i__5].r;
#line 299 "csymm.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 299 "csymm.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 300 "csymm.f"
/* L50: */
#line 300 "csymm.f"
		    }
#line 301 "csymm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 302 "csymm.f"
			i__3 = i__ + j * c_dim1;
#line 302 "csymm.f"
			i__4 = i__ + i__ * a_dim1;
#line 302 "csymm.f"
			z__2.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i, 
				z__2.i = temp1.r * a[i__4].i + temp1.i * a[
				i__4].r;
#line 302 "csymm.f"
			z__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 302 "csymm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 302 "csymm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 303 "csymm.f"
		    } else {
#line 304 "csymm.f"
			i__3 = i__ + j * c_dim1;
#line 304 "csymm.f"
			i__4 = i__ + j * c_dim1;
#line 304 "csymm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 304 "csymm.f"
			i__5 = i__ + i__ * a_dim1;
#line 304 "csymm.f"
			z__4.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
				z__4.i = temp1.r * a[i__5].i + temp1.i * a[
				i__5].r;
#line 304 "csymm.f"
			z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 304 "csymm.f"
			z__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 304 "csymm.f"
			z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 304 "csymm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 306 "csymm.f"
		    }
#line 307 "csymm.f"
/* L60: */
#line 307 "csymm.f"
		}
#line 308 "csymm.f"
/* L70: */
#line 308 "csymm.f"
	    }
#line 309 "csymm.f"
	} else {
#line 310 "csymm.f"
	    i__1 = *n;
#line 310 "csymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 311 "csymm.f"
		for (i__ = *m; i__ >= 1; --i__) {
#line 312 "csymm.f"
		    i__2 = i__ + j * b_dim1;
#line 312 "csymm.f"
		    z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2].i, 
			    z__1.i = alpha->r * b[i__2].i + alpha->i * b[i__2]
			    .r;
#line 312 "csymm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 313 "csymm.f"
		    temp2.r = 0., temp2.i = 0.;
#line 314 "csymm.f"
		    i__2 = *m;
#line 314 "csymm.f"
		    for (k = i__ + 1; k <= i__2; ++k) {
#line 315 "csymm.f"
			i__3 = k + j * c_dim1;
#line 315 "csymm.f"
			i__4 = k + j * c_dim1;
#line 315 "csymm.f"
			i__5 = k + i__ * a_dim1;
#line 315 "csymm.f"
			z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
				z__2.i = temp1.r * a[i__5].i + temp1.i * a[
				i__5].r;
#line 315 "csymm.f"
			z__1.r = c__[i__4].r + z__2.r, z__1.i = c__[i__4].i + 
				z__2.i;
#line 315 "csymm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 316 "csymm.f"
			i__3 = k + j * b_dim1;
#line 316 "csymm.f"
			i__4 = k + i__ * a_dim1;
#line 316 "csymm.f"
			z__2.r = b[i__3].r * a[i__4].r - b[i__3].i * a[i__4]
				.i, z__2.i = b[i__3].r * a[i__4].i + b[i__3]
				.i * a[i__4].r;
#line 316 "csymm.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 316 "csymm.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 317 "csymm.f"
/* L80: */
#line 317 "csymm.f"
		    }
#line 318 "csymm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 319 "csymm.f"
			i__2 = i__ + j * c_dim1;
#line 319 "csymm.f"
			i__3 = i__ + i__ * a_dim1;
#line 319 "csymm.f"
			z__2.r = temp1.r * a[i__3].r - temp1.i * a[i__3].i, 
				z__2.i = temp1.r * a[i__3].i + temp1.i * a[
				i__3].r;
#line 319 "csymm.f"
			z__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 319 "csymm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 319 "csymm.f"
			c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 320 "csymm.f"
		    } else {
#line 321 "csymm.f"
			i__2 = i__ + j * c_dim1;
#line 321 "csymm.f"
			i__3 = i__ + j * c_dim1;
#line 321 "csymm.f"
			z__3.r = beta->r * c__[i__3].r - beta->i * c__[i__3]
				.i, z__3.i = beta->r * c__[i__3].i + beta->i *
				 c__[i__3].r;
#line 321 "csymm.f"
			i__4 = i__ + i__ * a_dim1;
#line 321 "csymm.f"
			z__4.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i, 
				z__4.i = temp1.r * a[i__4].i + temp1.i * a[
				i__4].r;
#line 321 "csymm.f"
			z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 321 "csymm.f"
			z__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 321 "csymm.f"
			z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 321 "csymm.f"
			c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 323 "csymm.f"
		    }
#line 324 "csymm.f"
/* L90: */
#line 324 "csymm.f"
		}
#line 325 "csymm.f"
/* L100: */
#line 325 "csymm.f"
	    }
#line 326 "csymm.f"
	}
#line 327 "csymm.f"
    } else {

/*        Form  C := alpha*B*A + beta*C. */

#line 331 "csymm.f"
	i__1 = *n;
#line 331 "csymm.f"
	for (j = 1; j <= i__1; ++j) {
#line 332 "csymm.f"
	    i__2 = j + j * a_dim1;
#line 332 "csymm.f"
	    z__1.r = alpha->r * a[i__2].r - alpha->i * a[i__2].i, z__1.i = 
		    alpha->r * a[i__2].i + alpha->i * a[i__2].r;
#line 332 "csymm.f"
	    temp1.r = z__1.r, temp1.i = z__1.i;
#line 333 "csymm.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 334 "csymm.f"
		i__2 = *m;
#line 334 "csymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 335 "csymm.f"
		    i__3 = i__ + j * c_dim1;
#line 335 "csymm.f"
		    i__4 = i__ + j * b_dim1;
#line 335 "csymm.f"
		    z__1.r = temp1.r * b[i__4].r - temp1.i * b[i__4].i, 
			    z__1.i = temp1.r * b[i__4].i + temp1.i * b[i__4]
			    .r;
#line 335 "csymm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 336 "csymm.f"
/* L110: */
#line 336 "csymm.f"
		}
#line 337 "csymm.f"
	    } else {
#line 338 "csymm.f"
		i__2 = *m;
#line 338 "csymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 339 "csymm.f"
		    i__3 = i__ + j * c_dim1;
#line 339 "csymm.f"
		    i__4 = i__ + j * c_dim1;
#line 339 "csymm.f"
		    z__2.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    z__2.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
#line 339 "csymm.f"
		    i__5 = i__ + j * b_dim1;
#line 339 "csymm.f"
		    z__3.r = temp1.r * b[i__5].r - temp1.i * b[i__5].i, 
			    z__3.i = temp1.r * b[i__5].i + temp1.i * b[i__5]
			    .r;
#line 339 "csymm.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 339 "csymm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 340 "csymm.f"
/* L120: */
#line 340 "csymm.f"
		}
#line 341 "csymm.f"
	    }
#line 342 "csymm.f"
	    i__2 = j - 1;
#line 342 "csymm.f"
	    for (k = 1; k <= i__2; ++k) {
#line 343 "csymm.f"
		if (upper) {
#line 344 "csymm.f"
		    i__3 = k + j * a_dim1;
#line 344 "csymm.f"
		    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
			    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3]
			    .r;
#line 344 "csymm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 345 "csymm.f"
		} else {
#line 346 "csymm.f"
		    i__3 = j + k * a_dim1;
#line 346 "csymm.f"
		    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
			    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3]
			    .r;
#line 346 "csymm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 347 "csymm.f"
		}
#line 348 "csymm.f"
		i__3 = *m;
#line 348 "csymm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 349 "csymm.f"
		    i__4 = i__ + j * c_dim1;
#line 349 "csymm.f"
		    i__5 = i__ + j * c_dim1;
#line 349 "csymm.f"
		    i__6 = i__ + k * b_dim1;
#line 349 "csymm.f"
		    z__2.r = temp1.r * b[i__6].r - temp1.i * b[i__6].i, 
			    z__2.i = temp1.r * b[i__6].i + temp1.i * b[i__6]
			    .r;
#line 349 "csymm.f"
		    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
			    z__2.i;
#line 349 "csymm.f"
		    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 350 "csymm.f"
/* L130: */
#line 350 "csymm.f"
		}
#line 351 "csymm.f"
/* L140: */
#line 351 "csymm.f"
	    }
#line 352 "csymm.f"
	    i__2 = *n;
#line 352 "csymm.f"
	    for (k = j + 1; k <= i__2; ++k) {
#line 353 "csymm.f"
		if (upper) {
#line 354 "csymm.f"
		    i__3 = j + k * a_dim1;
#line 354 "csymm.f"
		    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
			    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3]
			    .r;
#line 354 "csymm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 355 "csymm.f"
		} else {
#line 356 "csymm.f"
		    i__3 = k + j * a_dim1;
#line 356 "csymm.f"
		    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
			    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3]
			    .r;
#line 356 "csymm.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 357 "csymm.f"
		}
#line 358 "csymm.f"
		i__3 = *m;
#line 358 "csymm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 359 "csymm.f"
		    i__4 = i__ + j * c_dim1;
#line 359 "csymm.f"
		    i__5 = i__ + j * c_dim1;
#line 359 "csymm.f"
		    i__6 = i__ + k * b_dim1;
#line 359 "csymm.f"
		    z__2.r = temp1.r * b[i__6].r - temp1.i * b[i__6].i, 
			    z__2.i = temp1.r * b[i__6].i + temp1.i * b[i__6]
			    .r;
#line 359 "csymm.f"
		    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
			    z__2.i;
#line 359 "csymm.f"
		    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 360 "csymm.f"
/* L150: */
#line 360 "csymm.f"
		}
#line 361 "csymm.f"
/* L160: */
#line 361 "csymm.f"
	    }
#line 362 "csymm.f"
/* L170: */
#line 362 "csymm.f"
	}
#line 363 "csymm.f"
    }

#line 365 "csymm.f"
    return 0;

/*     End of CSYMM . */

} /* csymm_ */

