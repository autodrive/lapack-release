#line 1 "zgemm.f"
/* zgemm.f -- translated by f2c (version 20100827).
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

#line 1 "zgemm.f"
/* > \brief \b ZGEMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,M,N */
/*       CHARACTER TRANSA,TRANSB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*op( A )*op( B ) + beta*C, */
/* > */
/* > where  op( X ) is one of */
/* > */
/* >    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H, */
/* > */
/* > alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/* > an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n',  op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't',  op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c',  op( A ) = A**H. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSB */
/* > \verbatim */
/* >          TRANSB is CHARACTER*1 */
/* >           On entry, TRANSB specifies the form of op( B ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSB = 'N' or 'n',  op( B ) = B. */
/* > */
/* >              TRANSB = 'T' or 't',  op( B ) = B**T. */
/* > */
/* >              TRANSB = 'C' or 'c',  op( B ) = B**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry,  M  specifies  the number  of rows  of the  matrix */
/* >           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N  specifies the number  of columns of the matrix */
/* >           op( B ) and the number of columns of the matrix C. N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry,  K  specifies  the number of columns of the matrix */
/* >           op( A ) and the number of rows of the matrix op( B ). K must */
/* >           be at least  zero. */
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
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is */
/* >           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/* >           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by m  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/* >           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/* >           least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is */
/* >           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/* >           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  n by k  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/* >           LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/* >           least  max( 1, n ). */
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
/* >          C is COMPLEX*16 array of DIMENSION ( LDC, n ). */
/* >           Before entry, the leading  m by n  part of the array  C must */
/* >           contain the matrix  C,  except when  beta  is zero, in which */
/* >           case C need not be set on entry. */
/* >           On exit, the array  C  is overwritten by the  m by n  matrix */
/* >           ( alpha*op( A )*op( B ) + beta*C ). */
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

/* > \date November 2011 */

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
/* Subroutine */ int zgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *
	c__, integer *ldc, ftnlen transa_len, ftnlen transb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, info;
    static logical nota, notb;
    static doublecomplex temp;
    static logical conja, conjb;
    static integer ncola;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and */
/*     B  respectively are to be  transposed but  not conjugated  and set */
/*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A */
/*     and the number of rows of  B  respectively. */

#line 234 "zgemm.f"
    /* Parameter adjustments */
#line 234 "zgemm.f"
    a_dim1 = *lda;
#line 234 "zgemm.f"
    a_offset = 1 + a_dim1;
#line 234 "zgemm.f"
    a -= a_offset;
#line 234 "zgemm.f"
    b_dim1 = *ldb;
#line 234 "zgemm.f"
    b_offset = 1 + b_dim1;
#line 234 "zgemm.f"
    b -= b_offset;
#line 234 "zgemm.f"
    c_dim1 = *ldc;
#line 234 "zgemm.f"
    c_offset = 1 + c_dim1;
#line 234 "zgemm.f"
    c__ -= c_offset;
#line 234 "zgemm.f"

#line 234 "zgemm.f"
    /* Function Body */
#line 234 "zgemm.f"
    nota = lsame_(transa, "N", (ftnlen)1, (ftnlen)1);
#line 235 "zgemm.f"
    notb = lsame_(transb, "N", (ftnlen)1, (ftnlen)1);
#line 236 "zgemm.f"
    conja = lsame_(transa, "C", (ftnlen)1, (ftnlen)1);
#line 237 "zgemm.f"
    conjb = lsame_(transb, "C", (ftnlen)1, (ftnlen)1);
#line 238 "zgemm.f"
    if (nota) {
#line 239 "zgemm.f"
	nrowa = *m;
#line 240 "zgemm.f"
	ncola = *k;
#line 241 "zgemm.f"
    } else {
#line 242 "zgemm.f"
	nrowa = *k;
#line 243 "zgemm.f"
	ncola = *m;
#line 244 "zgemm.f"
    }
#line 245 "zgemm.f"
    if (notb) {
#line 246 "zgemm.f"
	nrowb = *k;
#line 247 "zgemm.f"
    } else {
#line 248 "zgemm.f"
	nrowb = *n;
#line 249 "zgemm.f"
    }

/*     Test the input parameters. */

#line 253 "zgemm.f"
    info = 0;
#line 254 "zgemm.f"
    if (! nota && ! conja && ! lsame_(transa, "T", (ftnlen)1, (ftnlen)1)) {
#line 256 "zgemm.f"
	info = 1;
#line 257 "zgemm.f"
    } else if (! notb && ! conjb && ! lsame_(transb, "T", (ftnlen)1, (ftnlen)
	    1)) {
#line 259 "zgemm.f"
	info = 2;
#line 260 "zgemm.f"
    } else if (*m < 0) {
#line 261 "zgemm.f"
	info = 3;
#line 262 "zgemm.f"
    } else if (*n < 0) {
#line 263 "zgemm.f"
	info = 4;
#line 264 "zgemm.f"
    } else if (*k < 0) {
#line 265 "zgemm.f"
	info = 5;
#line 266 "zgemm.f"
    } else if (*lda < max(1,nrowa)) {
#line 267 "zgemm.f"
	info = 8;
#line 268 "zgemm.f"
    } else if (*ldb < max(1,nrowb)) {
#line 269 "zgemm.f"
	info = 10;
#line 270 "zgemm.f"
    } else if (*ldc < max(1,*m)) {
#line 271 "zgemm.f"
	info = 13;
#line 272 "zgemm.f"
    }
#line 273 "zgemm.f"
    if (info != 0) {
#line 274 "zgemm.f"
	xerbla_("ZGEMM ", &info, (ftnlen)6);
#line 275 "zgemm.f"
	return 0;
#line 276 "zgemm.f"
    }

/*     Quick return if possible. */

#line 280 "zgemm.f"
    if (*m == 0 || *n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) &&
	     (beta->r == 1. && beta->i == 0.)) {
#line 280 "zgemm.f"
	return 0;
#line 280 "zgemm.f"
    }

/*     And when  alpha.eq.zero. */

#line 285 "zgemm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 286 "zgemm.f"
	if (beta->r == 0. && beta->i == 0.) {
#line 287 "zgemm.f"
	    i__1 = *n;
#line 287 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 288 "zgemm.f"
		i__2 = *m;
#line 288 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "zgemm.f"
		    i__3 = i__ + j * c_dim1;
#line 289 "zgemm.f"
		    c__[i__3].r = 0., c__[i__3].i = 0.;
#line 290 "zgemm.f"
/* L10: */
#line 290 "zgemm.f"
		}
#line 291 "zgemm.f"
/* L20: */
#line 291 "zgemm.f"
	    }
#line 292 "zgemm.f"
	} else {
#line 293 "zgemm.f"
	    i__1 = *n;
#line 293 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 294 "zgemm.f"
		i__2 = *m;
#line 294 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "zgemm.f"
		    i__3 = i__ + j * c_dim1;
#line 295 "zgemm.f"
		    i__4 = i__ + j * c_dim1;
#line 295 "zgemm.f"
		    z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    z__1.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
#line 295 "zgemm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 296 "zgemm.f"
/* L30: */
#line 296 "zgemm.f"
		}
#line 297 "zgemm.f"
/* L40: */
#line 297 "zgemm.f"
	    }
#line 298 "zgemm.f"
	}
#line 299 "zgemm.f"
	return 0;
#line 300 "zgemm.f"
    }

/*     Start the operations. */

#line 304 "zgemm.f"
    if (notb) {
#line 305 "zgemm.f"
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

#line 309 "zgemm.f"
	    i__1 = *n;
#line 309 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 310 "zgemm.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 311 "zgemm.f"
		    i__2 = *m;
#line 311 "zgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 312 "zgemm.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 313 "zgemm.f"
/* L50: */
#line 313 "zgemm.f"
		    }
#line 314 "zgemm.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 315 "zgemm.f"
		    i__2 = *m;
#line 315 "zgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 316 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 316 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 316 "zgemm.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 316 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 317 "zgemm.f"
/* L60: */
#line 317 "zgemm.f"
		    }
#line 318 "zgemm.f"
		}
#line 319 "zgemm.f"
		i__2 = *k;
#line 319 "zgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 320 "zgemm.f"
		    i__3 = l + j * b_dim1;
#line 320 "zgemm.f"
		    if (b[i__3].r != 0. || b[i__3].i != 0.) {
#line 321 "zgemm.f"
			i__3 = l + j * b_dim1;
#line 321 "zgemm.f"
			z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				z__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
#line 321 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 322 "zgemm.f"
			i__3 = *m;
#line 322 "zgemm.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 323 "zgemm.f"
			    i__4 = i__ + j * c_dim1;
#line 323 "zgemm.f"
			    i__5 = i__ + j * c_dim1;
#line 323 "zgemm.f"
			    i__6 = i__ + l * a_dim1;
#line 323 "zgemm.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 323 "zgemm.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 323 "zgemm.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 324 "zgemm.f"
/* L70: */
#line 324 "zgemm.f"
			}
#line 325 "zgemm.f"
		    }
#line 326 "zgemm.f"
/* L80: */
#line 326 "zgemm.f"
		}
#line 327 "zgemm.f"
/* L90: */
#line 327 "zgemm.f"
	    }
#line 328 "zgemm.f"
	} else if (conja) {

/*           Form  C := alpha*A**H*B + beta*C. */

#line 332 "zgemm.f"
	    i__1 = *n;
#line 332 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 333 "zgemm.f"
		i__2 = *m;
#line 333 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 334 "zgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 335 "zgemm.f"
		    i__3 = *k;
#line 335 "zgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 336 "zgemm.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 336 "zgemm.f"
			i__4 = l + j * b_dim1;
#line 336 "zgemm.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 336 "zgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 336 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 337 "zgemm.f"
/* L100: */
#line 337 "zgemm.f"
		    }
#line 338 "zgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 339 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 339 "zgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 339 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 340 "zgemm.f"
		    } else {
#line 341 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 341 "zgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 341 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 341 "zgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 341 "zgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 341 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 342 "zgemm.f"
		    }
#line 343 "zgemm.f"
/* L110: */
#line 343 "zgemm.f"
		}
#line 344 "zgemm.f"
/* L120: */
#line 344 "zgemm.f"
	    }
#line 345 "zgemm.f"
	} else {

/*           Form  C := alpha*A**T*B + beta*C */

#line 349 "zgemm.f"
	    i__1 = *n;
#line 349 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 350 "zgemm.f"
		i__2 = *m;
#line 350 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 351 "zgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 352 "zgemm.f"
		    i__3 = *k;
#line 352 "zgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 353 "zgemm.f"
			i__4 = l + i__ * a_dim1;
#line 353 "zgemm.f"
			i__5 = l + j * b_dim1;
#line 353 "zgemm.f"
			z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, z__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
#line 353 "zgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 353 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 354 "zgemm.f"
/* L130: */
#line 354 "zgemm.f"
		    }
#line 355 "zgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 356 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 356 "zgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 356 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 357 "zgemm.f"
		    } else {
#line 358 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 358 "zgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 358 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 358 "zgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 358 "zgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 358 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 359 "zgemm.f"
		    }
#line 360 "zgemm.f"
/* L140: */
#line 360 "zgemm.f"
		}
#line 361 "zgemm.f"
/* L150: */
#line 361 "zgemm.f"
	    }
#line 362 "zgemm.f"
	}
#line 363 "zgemm.f"
    } else if (nota) {
#line 364 "zgemm.f"
	if (conjb) {

/*           Form  C := alpha*A*B**H + beta*C. */

#line 368 "zgemm.f"
	    i__1 = *n;
#line 368 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 369 "zgemm.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 370 "zgemm.f"
		    i__2 = *m;
#line 370 "zgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 371 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 371 "zgemm.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 372 "zgemm.f"
/* L160: */
#line 372 "zgemm.f"
		    }
#line 373 "zgemm.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 374 "zgemm.f"
		    i__2 = *m;
#line 374 "zgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 375 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 375 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 375 "zgemm.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 375 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 376 "zgemm.f"
/* L170: */
#line 376 "zgemm.f"
		    }
#line 377 "zgemm.f"
		}
#line 378 "zgemm.f"
		i__2 = *k;
#line 378 "zgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 379 "zgemm.f"
		    i__3 = j + l * b_dim1;
#line 379 "zgemm.f"
		    if (b[i__3].r != 0. || b[i__3].i != 0.) {
#line 380 "zgemm.f"
			d_cnjg(&z__2, &b[j + l * b_dim1]);
#line 380 "zgemm.f"
			z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, 
				z__1.i = alpha->r * z__2.i + alpha->i * 
				z__2.r;
#line 380 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 381 "zgemm.f"
			i__3 = *m;
#line 381 "zgemm.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 382 "zgemm.f"
			    i__4 = i__ + j * c_dim1;
#line 382 "zgemm.f"
			    i__5 = i__ + j * c_dim1;
#line 382 "zgemm.f"
			    i__6 = i__ + l * a_dim1;
#line 382 "zgemm.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 382 "zgemm.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 382 "zgemm.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 383 "zgemm.f"
/* L180: */
#line 383 "zgemm.f"
			}
#line 384 "zgemm.f"
		    }
#line 385 "zgemm.f"
/* L190: */
#line 385 "zgemm.f"
		}
#line 386 "zgemm.f"
/* L200: */
#line 386 "zgemm.f"
	    }
#line 387 "zgemm.f"
	} else {

/*           Form  C := alpha*A*B**T          + beta*C */

#line 391 "zgemm.f"
	    i__1 = *n;
#line 391 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 392 "zgemm.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 393 "zgemm.f"
		    i__2 = *m;
#line 393 "zgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 394 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 394 "zgemm.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 395 "zgemm.f"
/* L210: */
#line 395 "zgemm.f"
		    }
#line 396 "zgemm.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 397 "zgemm.f"
		    i__2 = *m;
#line 397 "zgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 398 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 398 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 398 "zgemm.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 398 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 399 "zgemm.f"
/* L220: */
#line 399 "zgemm.f"
		    }
#line 400 "zgemm.f"
		}
#line 401 "zgemm.f"
		i__2 = *k;
#line 401 "zgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 402 "zgemm.f"
		    i__3 = j + l * b_dim1;
#line 402 "zgemm.f"
		    if (b[i__3].r != 0. || b[i__3].i != 0.) {
#line 403 "zgemm.f"
			i__3 = j + l * b_dim1;
#line 403 "zgemm.f"
			z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				z__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
#line 403 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 404 "zgemm.f"
			i__3 = *m;
#line 404 "zgemm.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 405 "zgemm.f"
			    i__4 = i__ + j * c_dim1;
#line 405 "zgemm.f"
			    i__5 = i__ + j * c_dim1;
#line 405 "zgemm.f"
			    i__6 = i__ + l * a_dim1;
#line 405 "zgemm.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 405 "zgemm.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 405 "zgemm.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 406 "zgemm.f"
/* L230: */
#line 406 "zgemm.f"
			}
#line 407 "zgemm.f"
		    }
#line 408 "zgemm.f"
/* L240: */
#line 408 "zgemm.f"
		}
#line 409 "zgemm.f"
/* L250: */
#line 409 "zgemm.f"
	    }
#line 410 "zgemm.f"
	}
#line 411 "zgemm.f"
    } else if (conja) {
#line 412 "zgemm.f"
	if (conjb) {

/*           Form  C := alpha*A**H*B**H + beta*C. */

#line 416 "zgemm.f"
	    i__1 = *n;
#line 416 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 417 "zgemm.f"
		i__2 = *m;
#line 417 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "zgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 419 "zgemm.f"
		    i__3 = *k;
#line 419 "zgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 420 "zgemm.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 420 "zgemm.f"
			d_cnjg(&z__4, &b[j + l * b_dim1]);
#line 420 "zgemm.f"
			z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = 
				z__3.r * z__4.i + z__3.i * z__4.r;
#line 420 "zgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 420 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 421 "zgemm.f"
/* L260: */
#line 421 "zgemm.f"
		    }
#line 422 "zgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 423 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 423 "zgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 423 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 424 "zgemm.f"
		    } else {
#line 425 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 425 "zgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 425 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 425 "zgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 425 "zgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 425 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 426 "zgemm.f"
		    }
#line 427 "zgemm.f"
/* L270: */
#line 427 "zgemm.f"
		}
#line 428 "zgemm.f"
/* L280: */
#line 428 "zgemm.f"
	    }
#line 429 "zgemm.f"
	} else {

/*           Form  C := alpha*A**H*B**T + beta*C */

#line 433 "zgemm.f"
	    i__1 = *n;
#line 433 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 434 "zgemm.f"
		i__2 = *m;
#line 434 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 435 "zgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 436 "zgemm.f"
		    i__3 = *k;
#line 436 "zgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 437 "zgemm.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 437 "zgemm.f"
			i__4 = j + l * b_dim1;
#line 437 "zgemm.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 437 "zgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 437 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 438 "zgemm.f"
/* L290: */
#line 438 "zgemm.f"
		    }
#line 439 "zgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 440 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 440 "zgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 440 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 441 "zgemm.f"
		    } else {
#line 442 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 442 "zgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 442 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 442 "zgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 442 "zgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 442 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 443 "zgemm.f"
		    }
#line 444 "zgemm.f"
/* L300: */
#line 444 "zgemm.f"
		}
#line 445 "zgemm.f"
/* L310: */
#line 445 "zgemm.f"
	    }
#line 446 "zgemm.f"
	}
#line 447 "zgemm.f"
    } else {
#line 448 "zgemm.f"
	if (conjb) {

/*           Form  C := alpha*A**T*B**H + beta*C */

#line 452 "zgemm.f"
	    i__1 = *n;
#line 452 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 453 "zgemm.f"
		i__2 = *m;
#line 453 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 454 "zgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 455 "zgemm.f"
		    i__3 = *k;
#line 455 "zgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 456 "zgemm.f"
			i__4 = l + i__ * a_dim1;
#line 456 "zgemm.f"
			d_cnjg(&z__3, &b[j + l * b_dim1]);
#line 456 "zgemm.f"
			z__2.r = a[i__4].r * z__3.r - a[i__4].i * z__3.i, 
				z__2.i = a[i__4].r * z__3.i + a[i__4].i * 
				z__3.r;
#line 456 "zgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 456 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 457 "zgemm.f"
/* L320: */
#line 457 "zgemm.f"
		    }
#line 458 "zgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 459 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 459 "zgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 459 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 460 "zgemm.f"
		    } else {
#line 461 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 461 "zgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 461 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 461 "zgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 461 "zgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 461 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 462 "zgemm.f"
		    }
#line 463 "zgemm.f"
/* L330: */
#line 463 "zgemm.f"
		}
#line 464 "zgemm.f"
/* L340: */
#line 464 "zgemm.f"
	    }
#line 465 "zgemm.f"
	} else {

/*           Form  C := alpha*A**T*B**T + beta*C */

#line 469 "zgemm.f"
	    i__1 = *n;
#line 469 "zgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 470 "zgemm.f"
		i__2 = *m;
#line 470 "zgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 471 "zgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 472 "zgemm.f"
		    i__3 = *k;
#line 472 "zgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 473 "zgemm.f"
			i__4 = l + i__ * a_dim1;
#line 473 "zgemm.f"
			i__5 = j + l * b_dim1;
#line 473 "zgemm.f"
			z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, z__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
#line 473 "zgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 473 "zgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 474 "zgemm.f"
/* L350: */
#line 474 "zgemm.f"
		    }
#line 475 "zgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 476 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 476 "zgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 476 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 477 "zgemm.f"
		    } else {
#line 478 "zgemm.f"
			i__3 = i__ + j * c_dim1;
#line 478 "zgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 478 "zgemm.f"
			i__4 = i__ + j * c_dim1;
#line 478 "zgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 478 "zgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 478 "zgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 479 "zgemm.f"
		    }
#line 480 "zgemm.f"
/* L360: */
#line 480 "zgemm.f"
		}
#line 481 "zgemm.f"
/* L370: */
#line 481 "zgemm.f"
	    }
#line 482 "zgemm.f"
	}
#line 483 "zgemm.f"
    }

#line 485 "zgemm.f"
    return 0;

/*     End of ZGEMM . */

} /* zgemm_ */

