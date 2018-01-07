#line 1 "cgemm.f"
/* cgemm.f -- translated by f2c (version 20100827).
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

#line 1 "cgemm.f"
/* > \brief \b CGEMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,M,N */
/*       CHARACTER TRANSA,TRANSB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEMM  performs one of the matrix-matrix operations */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension ( LDA, ka ), where ka is */
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
/* >          B is COMPLEX array, dimension ( LDB, kb ), where kb is */
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
/* >          BETA is COMPLEX */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension ( LDC, N ) */
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
/* Subroutine */ int cgemm_(char *transa, char *transb, integer *m, integer *
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

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and */
/*     B  respectively are to be  transposed but  not conjugated  and set */
/*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A */
/*     and the number of rows of  B  respectively. */

#line 234 "cgemm.f"
    /* Parameter adjustments */
#line 234 "cgemm.f"
    a_dim1 = *lda;
#line 234 "cgemm.f"
    a_offset = 1 + a_dim1;
#line 234 "cgemm.f"
    a -= a_offset;
#line 234 "cgemm.f"
    b_dim1 = *ldb;
#line 234 "cgemm.f"
    b_offset = 1 + b_dim1;
#line 234 "cgemm.f"
    b -= b_offset;
#line 234 "cgemm.f"
    c_dim1 = *ldc;
#line 234 "cgemm.f"
    c_offset = 1 + c_dim1;
#line 234 "cgemm.f"
    c__ -= c_offset;
#line 234 "cgemm.f"

#line 234 "cgemm.f"
    /* Function Body */
#line 234 "cgemm.f"
    nota = lsame_(transa, "N", (ftnlen)1, (ftnlen)1);
#line 235 "cgemm.f"
    notb = lsame_(transb, "N", (ftnlen)1, (ftnlen)1);
#line 236 "cgemm.f"
    conja = lsame_(transa, "C", (ftnlen)1, (ftnlen)1);
#line 237 "cgemm.f"
    conjb = lsame_(transb, "C", (ftnlen)1, (ftnlen)1);
#line 238 "cgemm.f"
    if (nota) {
#line 239 "cgemm.f"
	nrowa = *m;
#line 240 "cgemm.f"
	ncola = *k;
#line 241 "cgemm.f"
    } else {
#line 242 "cgemm.f"
	nrowa = *k;
#line 243 "cgemm.f"
	ncola = *m;
#line 244 "cgemm.f"
    }
#line 245 "cgemm.f"
    if (notb) {
#line 246 "cgemm.f"
	nrowb = *k;
#line 247 "cgemm.f"
    } else {
#line 248 "cgemm.f"
	nrowb = *n;
#line 249 "cgemm.f"
    }

/*     Test the input parameters. */

#line 253 "cgemm.f"
    info = 0;
#line 254 "cgemm.f"
    if (! nota && ! conja && ! lsame_(transa, "T", (ftnlen)1, (ftnlen)1)) {
#line 256 "cgemm.f"
	info = 1;
#line 257 "cgemm.f"
    } else if (! notb && ! conjb && ! lsame_(transb, "T", (ftnlen)1, (ftnlen)
	    1)) {
#line 259 "cgemm.f"
	info = 2;
#line 260 "cgemm.f"
    } else if (*m < 0) {
#line 261 "cgemm.f"
	info = 3;
#line 262 "cgemm.f"
    } else if (*n < 0) {
#line 263 "cgemm.f"
	info = 4;
#line 264 "cgemm.f"
    } else if (*k < 0) {
#line 265 "cgemm.f"
	info = 5;
#line 266 "cgemm.f"
    } else if (*lda < max(1,nrowa)) {
#line 267 "cgemm.f"
	info = 8;
#line 268 "cgemm.f"
    } else if (*ldb < max(1,nrowb)) {
#line 269 "cgemm.f"
	info = 10;
#line 270 "cgemm.f"
    } else if (*ldc < max(1,*m)) {
#line 271 "cgemm.f"
	info = 13;
#line 272 "cgemm.f"
    }
#line 273 "cgemm.f"
    if (info != 0) {
#line 274 "cgemm.f"
	xerbla_("CGEMM ", &info, (ftnlen)6);
#line 275 "cgemm.f"
	return 0;
#line 276 "cgemm.f"
    }

/*     Quick return if possible. */

#line 280 "cgemm.f"
    if (*m == 0 || *n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) &&
	     (beta->r == 1. && beta->i == 0.)) {
#line 280 "cgemm.f"
	return 0;
#line 280 "cgemm.f"
    }

/*     And when  alpha.eq.zero. */

#line 285 "cgemm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 286 "cgemm.f"
	if (beta->r == 0. && beta->i == 0.) {
#line 287 "cgemm.f"
	    i__1 = *n;
#line 287 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 288 "cgemm.f"
		i__2 = *m;
#line 288 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "cgemm.f"
		    i__3 = i__ + j * c_dim1;
#line 289 "cgemm.f"
		    c__[i__3].r = 0., c__[i__3].i = 0.;
#line 290 "cgemm.f"
/* L10: */
#line 290 "cgemm.f"
		}
#line 291 "cgemm.f"
/* L20: */
#line 291 "cgemm.f"
	    }
#line 292 "cgemm.f"
	} else {
#line 293 "cgemm.f"
	    i__1 = *n;
#line 293 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 294 "cgemm.f"
		i__2 = *m;
#line 294 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "cgemm.f"
		    i__3 = i__ + j * c_dim1;
#line 295 "cgemm.f"
		    i__4 = i__ + j * c_dim1;
#line 295 "cgemm.f"
		    z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i, 
			    z__1.i = beta->r * c__[i__4].i + beta->i * c__[
			    i__4].r;
#line 295 "cgemm.f"
		    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 296 "cgemm.f"
/* L30: */
#line 296 "cgemm.f"
		}
#line 297 "cgemm.f"
/* L40: */
#line 297 "cgemm.f"
	    }
#line 298 "cgemm.f"
	}
#line 299 "cgemm.f"
	return 0;
#line 300 "cgemm.f"
    }

/*     Start the operations. */

#line 304 "cgemm.f"
    if (notb) {
#line 305 "cgemm.f"
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

#line 309 "cgemm.f"
	    i__1 = *n;
#line 309 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 310 "cgemm.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 311 "cgemm.f"
		    i__2 = *m;
#line 311 "cgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 312 "cgemm.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 313 "cgemm.f"
/* L50: */
#line 313 "cgemm.f"
		    }
#line 314 "cgemm.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 315 "cgemm.f"
		    i__2 = *m;
#line 315 "cgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 316 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 316 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 316 "cgemm.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 316 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 317 "cgemm.f"
/* L60: */
#line 317 "cgemm.f"
		    }
#line 318 "cgemm.f"
		}
#line 319 "cgemm.f"
		i__2 = *k;
#line 319 "cgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 320 "cgemm.f"
		    i__3 = l + j * b_dim1;
#line 320 "cgemm.f"
		    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
			    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3]
			    .r;
#line 320 "cgemm.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 321 "cgemm.f"
		    i__3 = *m;
#line 321 "cgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 322 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 322 "cgemm.f"
			i__5 = i__ + j * c_dim1;
#line 322 "cgemm.f"
			i__6 = i__ + l * a_dim1;
#line 322 "cgemm.f"
			z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				z__2.i = temp.r * a[i__6].i + temp.i * a[i__6]
				.r;
#line 322 "cgemm.f"
			z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
				z__2.i;
#line 322 "cgemm.f"
			c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 323 "cgemm.f"
/* L70: */
#line 323 "cgemm.f"
		    }
#line 324 "cgemm.f"
/* L80: */
#line 324 "cgemm.f"
		}
#line 325 "cgemm.f"
/* L90: */
#line 325 "cgemm.f"
	    }
#line 326 "cgemm.f"
	} else if (conja) {

/*           Form  C := alpha*A**H*B + beta*C. */

#line 330 "cgemm.f"
	    i__1 = *n;
#line 330 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 331 "cgemm.f"
		i__2 = *m;
#line 331 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 332 "cgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 333 "cgemm.f"
		    i__3 = *k;
#line 333 "cgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 334 "cgemm.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 334 "cgemm.f"
			i__4 = l + j * b_dim1;
#line 334 "cgemm.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 334 "cgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 334 "cgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 335 "cgemm.f"
/* L100: */
#line 335 "cgemm.f"
		    }
#line 336 "cgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 337 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 337 "cgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 337 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 338 "cgemm.f"
		    } else {
#line 339 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 339 "cgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 339 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 339 "cgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 339 "cgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 339 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 340 "cgemm.f"
		    }
#line 341 "cgemm.f"
/* L110: */
#line 341 "cgemm.f"
		}
#line 342 "cgemm.f"
/* L120: */
#line 342 "cgemm.f"
	    }
#line 343 "cgemm.f"
	} else {

/*           Form  C := alpha*A**T*B + beta*C */

#line 347 "cgemm.f"
	    i__1 = *n;
#line 347 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 348 "cgemm.f"
		i__2 = *m;
#line 348 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 349 "cgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 350 "cgemm.f"
		    i__3 = *k;
#line 350 "cgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 351 "cgemm.f"
			i__4 = l + i__ * a_dim1;
#line 351 "cgemm.f"
			i__5 = l + j * b_dim1;
#line 351 "cgemm.f"
			z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, z__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
#line 351 "cgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 351 "cgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 352 "cgemm.f"
/* L130: */
#line 352 "cgemm.f"
		    }
#line 353 "cgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 354 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 354 "cgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 354 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 355 "cgemm.f"
		    } else {
#line 356 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 356 "cgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 356 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 356 "cgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 356 "cgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 356 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 357 "cgemm.f"
		    }
#line 358 "cgemm.f"
/* L140: */
#line 358 "cgemm.f"
		}
#line 359 "cgemm.f"
/* L150: */
#line 359 "cgemm.f"
	    }
#line 360 "cgemm.f"
	}
#line 361 "cgemm.f"
    } else if (nota) {
#line 362 "cgemm.f"
	if (conjb) {

/*           Form  C := alpha*A*B**H + beta*C. */

#line 366 "cgemm.f"
	    i__1 = *n;
#line 366 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 367 "cgemm.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 368 "cgemm.f"
		    i__2 = *m;
#line 368 "cgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 369 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 369 "cgemm.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 370 "cgemm.f"
/* L160: */
#line 370 "cgemm.f"
		    }
#line 371 "cgemm.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 372 "cgemm.f"
		    i__2 = *m;
#line 372 "cgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 373 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 373 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 373 "cgemm.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 373 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 374 "cgemm.f"
/* L170: */
#line 374 "cgemm.f"
		    }
#line 375 "cgemm.f"
		}
#line 376 "cgemm.f"
		i__2 = *k;
#line 376 "cgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 377 "cgemm.f"
		    d_cnjg(&z__2, &b[j + l * b_dim1]);
#line 377 "cgemm.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 377 "cgemm.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 378 "cgemm.f"
		    i__3 = *m;
#line 378 "cgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 379 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 379 "cgemm.f"
			i__5 = i__ + j * c_dim1;
#line 379 "cgemm.f"
			i__6 = i__ + l * a_dim1;
#line 379 "cgemm.f"
			z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				z__2.i = temp.r * a[i__6].i + temp.i * a[i__6]
				.r;
#line 379 "cgemm.f"
			z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
				z__2.i;
#line 379 "cgemm.f"
			c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 380 "cgemm.f"
/* L180: */
#line 380 "cgemm.f"
		    }
#line 381 "cgemm.f"
/* L190: */
#line 381 "cgemm.f"
		}
#line 382 "cgemm.f"
/* L200: */
#line 382 "cgemm.f"
	    }
#line 383 "cgemm.f"
	} else {

/*           Form  C := alpha*A*B**T + beta*C */

#line 387 "cgemm.f"
	    i__1 = *n;
#line 387 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 388 "cgemm.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 389 "cgemm.f"
		    i__2 = *m;
#line 389 "cgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 390 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 390 "cgemm.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 391 "cgemm.f"
/* L210: */
#line 391 "cgemm.f"
		    }
#line 392 "cgemm.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 393 "cgemm.f"
		    i__2 = *m;
#line 393 "cgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 394 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 394 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 394 "cgemm.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 394 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 395 "cgemm.f"
/* L220: */
#line 395 "cgemm.f"
		    }
#line 396 "cgemm.f"
		}
#line 397 "cgemm.f"
		i__2 = *k;
#line 397 "cgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 398 "cgemm.f"
		    i__3 = j + l * b_dim1;
#line 398 "cgemm.f"
		    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
			    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3]
			    .r;
#line 398 "cgemm.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 399 "cgemm.f"
		    i__3 = *m;
#line 399 "cgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 400 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 400 "cgemm.f"
			i__5 = i__ + j * c_dim1;
#line 400 "cgemm.f"
			i__6 = i__ + l * a_dim1;
#line 400 "cgemm.f"
			z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				z__2.i = temp.r * a[i__6].i + temp.i * a[i__6]
				.r;
#line 400 "cgemm.f"
			z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + 
				z__2.i;
#line 400 "cgemm.f"
			c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 401 "cgemm.f"
/* L230: */
#line 401 "cgemm.f"
		    }
#line 402 "cgemm.f"
/* L240: */
#line 402 "cgemm.f"
		}
#line 403 "cgemm.f"
/* L250: */
#line 403 "cgemm.f"
	    }
#line 404 "cgemm.f"
	}
#line 405 "cgemm.f"
    } else if (conja) {
#line 406 "cgemm.f"
	if (conjb) {

/*           Form  C := alpha*A**H*B**H + beta*C. */

#line 410 "cgemm.f"
	    i__1 = *n;
#line 410 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 411 "cgemm.f"
		i__2 = *m;
#line 411 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 412 "cgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 413 "cgemm.f"
		    i__3 = *k;
#line 413 "cgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 414 "cgemm.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 414 "cgemm.f"
			d_cnjg(&z__4, &b[j + l * b_dim1]);
#line 414 "cgemm.f"
			z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = 
				z__3.r * z__4.i + z__3.i * z__4.r;
#line 414 "cgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 414 "cgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 415 "cgemm.f"
/* L260: */
#line 415 "cgemm.f"
		    }
#line 416 "cgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 417 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 417 "cgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 417 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 418 "cgemm.f"
		    } else {
#line 419 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 419 "cgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 419 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 419 "cgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 419 "cgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 419 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 420 "cgemm.f"
		    }
#line 421 "cgemm.f"
/* L270: */
#line 421 "cgemm.f"
		}
#line 422 "cgemm.f"
/* L280: */
#line 422 "cgemm.f"
	    }
#line 423 "cgemm.f"
	} else {

/*           Form  C := alpha*A**H*B**T + beta*C */

#line 427 "cgemm.f"
	    i__1 = *n;
#line 427 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 428 "cgemm.f"
		i__2 = *m;
#line 428 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 429 "cgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 430 "cgemm.f"
		    i__3 = *k;
#line 430 "cgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 431 "cgemm.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 431 "cgemm.f"
			i__4 = j + l * b_dim1;
#line 431 "cgemm.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 431 "cgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 431 "cgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 432 "cgemm.f"
/* L290: */
#line 432 "cgemm.f"
		    }
#line 433 "cgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 434 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 434 "cgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 434 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 435 "cgemm.f"
		    } else {
#line 436 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 436 "cgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 436 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 436 "cgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 436 "cgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 436 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 437 "cgemm.f"
		    }
#line 438 "cgemm.f"
/* L300: */
#line 438 "cgemm.f"
		}
#line 439 "cgemm.f"
/* L310: */
#line 439 "cgemm.f"
	    }
#line 440 "cgemm.f"
	}
#line 441 "cgemm.f"
    } else {
#line 442 "cgemm.f"
	if (conjb) {

/*           Form  C := alpha*A**T*B**H + beta*C */

#line 446 "cgemm.f"
	    i__1 = *n;
#line 446 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 447 "cgemm.f"
		i__2 = *m;
#line 447 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "cgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 449 "cgemm.f"
		    i__3 = *k;
#line 449 "cgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 450 "cgemm.f"
			i__4 = l + i__ * a_dim1;
#line 450 "cgemm.f"
			d_cnjg(&z__3, &b[j + l * b_dim1]);
#line 450 "cgemm.f"
			z__2.r = a[i__4].r * z__3.r - a[i__4].i * z__3.i, 
				z__2.i = a[i__4].r * z__3.i + a[i__4].i * 
				z__3.r;
#line 450 "cgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 450 "cgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 451 "cgemm.f"
/* L320: */
#line 451 "cgemm.f"
		    }
#line 452 "cgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 453 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 453 "cgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 453 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 454 "cgemm.f"
		    } else {
#line 455 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 455 "cgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 455 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 455 "cgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 455 "cgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 455 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 456 "cgemm.f"
		    }
#line 457 "cgemm.f"
/* L330: */
#line 457 "cgemm.f"
		}
#line 458 "cgemm.f"
/* L340: */
#line 458 "cgemm.f"
	    }
#line 459 "cgemm.f"
	} else {

/*           Form  C := alpha*A**T*B**T + beta*C */

#line 463 "cgemm.f"
	    i__1 = *n;
#line 463 "cgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 464 "cgemm.f"
		i__2 = *m;
#line 464 "cgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 465 "cgemm.f"
		    temp.r = 0., temp.i = 0.;
#line 466 "cgemm.f"
		    i__3 = *k;
#line 466 "cgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 467 "cgemm.f"
			i__4 = l + i__ * a_dim1;
#line 467 "cgemm.f"
			i__5 = j + l * b_dim1;
#line 467 "cgemm.f"
			z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, z__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
#line 467 "cgemm.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 467 "cgemm.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 468 "cgemm.f"
/* L350: */
#line 468 "cgemm.f"
		    }
#line 469 "cgemm.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 470 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 470 "cgemm.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 470 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 471 "cgemm.f"
		    } else {
#line 472 "cgemm.f"
			i__3 = i__ + j * c_dim1;
#line 472 "cgemm.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 472 "cgemm.f"
			i__4 = i__ + j * c_dim1;
#line 472 "cgemm.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 472 "cgemm.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 472 "cgemm.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 473 "cgemm.f"
		    }
#line 474 "cgemm.f"
/* L360: */
#line 474 "cgemm.f"
		}
#line 475 "cgemm.f"
/* L370: */
#line 475 "cgemm.f"
	    }
#line 476 "cgemm.f"
	}
#line 477 "cgemm.f"
    }

#line 479 "cgemm.f"
    return 0;

/*     End of CGEMM . */

} /* cgemm_ */

