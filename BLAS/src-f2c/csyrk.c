#line 1 "csyrk.f"
/* csyrk.f -- translated by f2c (version 20100827).
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

#line 1 "csyrk.f"
/* > \brief \b CSYRK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER K,LDA,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYRK  performs one of the symmetric rank k operations */
/* > */
/* >    C := alpha*A*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars,  C is an  n by n symmetric matrix */
/* > and  A  is an  n by k  matrix in the first case and a  k by n  matrix */
/* > in the second case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of the  array  C  is to be  referenced  as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/* >                                  is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry,  TRANS  specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C. */
/* > */
/* >              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N specifies the order of the matrix C.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/* >           of  columns   of  the   matrix   A,   and  on   entry   with */
/* >           TRANS = 'T' or 't',  K  specifies  the number of rows of the */
/* >           matrix A.  K must be at least zero. */
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
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by n  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
/* >           be at least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array of DIMENSION ( LDC, n ). */
/* >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/* >           upper triangular part of the array C must contain the upper */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           lower triangular part of C is not referenced.  On exit, the */
/* >           upper triangular part of the array  C is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/* >           lower triangular part of the array C must contain the lower */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           upper triangular part of C is not referenced.  On exit, the */
/* >           lower triangular part of the array  C is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/* >           max( 1, n ). */
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
/* Subroutine */ int csyrk_(char *uplo, char *trans, integer *n, integer *k, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, j, l, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
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

/*     Test the input parameters. */

#line 210 "csyrk.f"
    /* Parameter adjustments */
#line 210 "csyrk.f"
    a_dim1 = *lda;
#line 210 "csyrk.f"
    a_offset = 1 + a_dim1;
#line 210 "csyrk.f"
    a -= a_offset;
#line 210 "csyrk.f"
    c_dim1 = *ldc;
#line 210 "csyrk.f"
    c_offset = 1 + c_dim1;
#line 210 "csyrk.f"
    c__ -= c_offset;
#line 210 "csyrk.f"

#line 210 "csyrk.f"
    /* Function Body */
#line 210 "csyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 211 "csyrk.f"
	nrowa = *n;
#line 212 "csyrk.f"
    } else {
#line 213 "csyrk.f"
	nrowa = *k;
#line 214 "csyrk.f"
    }
#line 215 "csyrk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 217 "csyrk.f"
    info = 0;
#line 218 "csyrk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 219 "csyrk.f"
	info = 1;
#line 220 "csyrk.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1)) {
#line 222 "csyrk.f"
	info = 2;
#line 223 "csyrk.f"
    } else if (*n < 0) {
#line 224 "csyrk.f"
	info = 3;
#line 225 "csyrk.f"
    } else if (*k < 0) {
#line 226 "csyrk.f"
	info = 4;
#line 227 "csyrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 228 "csyrk.f"
	info = 7;
#line 229 "csyrk.f"
    } else if (*ldc < max(1,*n)) {
#line 230 "csyrk.f"
	info = 10;
#line 231 "csyrk.f"
    }
#line 232 "csyrk.f"
    if (info != 0) {
#line 233 "csyrk.f"
	xerbla_("CSYRK ", &info, (ftnlen)6);
#line 234 "csyrk.f"
	return 0;
#line 235 "csyrk.f"
    }

/*     Quick return if possible. */

#line 239 "csyrk.f"
    if (*n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) && (beta->r 
	    == 1. && beta->i == 0.)) {
#line 239 "csyrk.f"
	return 0;
#line 239 "csyrk.f"
    }

/*     And when  alpha.eq.zero. */

#line 244 "csyrk.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 245 "csyrk.f"
	if (upper) {
#line 246 "csyrk.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 247 "csyrk.f"
		i__1 = *n;
#line 247 "csyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 248 "csyrk.f"
		    i__2 = j;
#line 248 "csyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 249 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 249 "csyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 250 "csyrk.f"
/* L10: */
#line 250 "csyrk.f"
		    }
#line 251 "csyrk.f"
/* L20: */
#line 251 "csyrk.f"
		}
#line 252 "csyrk.f"
	    } else {
#line 253 "csyrk.f"
		i__1 = *n;
#line 253 "csyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 254 "csyrk.f"
		    i__2 = j;
#line 254 "csyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 255 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 255 "csyrk.f"
			i__4 = i__ + j * c_dim1;
#line 255 "csyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 255 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 256 "csyrk.f"
/* L30: */
#line 256 "csyrk.f"
		    }
#line 257 "csyrk.f"
/* L40: */
#line 257 "csyrk.f"
		}
#line 258 "csyrk.f"
	    }
#line 259 "csyrk.f"
	} else {
#line 260 "csyrk.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 261 "csyrk.f"
		i__1 = *n;
#line 261 "csyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 262 "csyrk.f"
		    i__2 = *n;
#line 262 "csyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 263 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 263 "csyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 264 "csyrk.f"
/* L50: */
#line 264 "csyrk.f"
		    }
#line 265 "csyrk.f"
/* L60: */
#line 265 "csyrk.f"
		}
#line 266 "csyrk.f"
	    } else {
#line 267 "csyrk.f"
		i__1 = *n;
#line 267 "csyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "csyrk.f"
		    i__2 = *n;
#line 268 "csyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 269 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 269 "csyrk.f"
			i__4 = i__ + j * c_dim1;
#line 269 "csyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 269 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 270 "csyrk.f"
/* L70: */
#line 270 "csyrk.f"
		    }
#line 271 "csyrk.f"
/* L80: */
#line 271 "csyrk.f"
		}
#line 272 "csyrk.f"
	    }
#line 273 "csyrk.f"
	}
#line 274 "csyrk.f"
	return 0;
#line 275 "csyrk.f"
    }

/*     Start the operations. */

#line 279 "csyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*A**T + beta*C. */

#line 283 "csyrk.f"
	if (upper) {
#line 284 "csyrk.f"
	    i__1 = *n;
#line 284 "csyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 285 "csyrk.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 286 "csyrk.f"
		    i__2 = j;
#line 286 "csyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 287 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 287 "csyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 288 "csyrk.f"
/* L90: */
#line 288 "csyrk.f"
		    }
#line 289 "csyrk.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 290 "csyrk.f"
		    i__2 = j;
#line 290 "csyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 291 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 291 "csyrk.f"
			i__4 = i__ + j * c_dim1;
#line 291 "csyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 291 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 292 "csyrk.f"
/* L100: */
#line 292 "csyrk.f"
		    }
#line 293 "csyrk.f"
		}
#line 294 "csyrk.f"
		i__2 = *k;
#line 294 "csyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 295 "csyrk.f"
		    i__3 = j + l * a_dim1;
#line 295 "csyrk.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 296 "csyrk.f"
			i__3 = j + l * a_dim1;
#line 296 "csyrk.f"
			z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__1.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 296 "csyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 297 "csyrk.f"
			i__3 = j;
#line 297 "csyrk.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 298 "csyrk.f"
			    i__4 = i__ + j * c_dim1;
#line 298 "csyrk.f"
			    i__5 = i__ + j * c_dim1;
#line 298 "csyrk.f"
			    i__6 = i__ + l * a_dim1;
#line 298 "csyrk.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 298 "csyrk.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 298 "csyrk.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 299 "csyrk.f"
/* L110: */
#line 299 "csyrk.f"
			}
#line 300 "csyrk.f"
		    }
#line 301 "csyrk.f"
/* L120: */
#line 301 "csyrk.f"
		}
#line 302 "csyrk.f"
/* L130: */
#line 302 "csyrk.f"
	    }
#line 303 "csyrk.f"
	} else {
#line 304 "csyrk.f"
	    i__1 = *n;
#line 304 "csyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 305 "csyrk.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 306 "csyrk.f"
		    i__2 = *n;
#line 306 "csyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 307 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 307 "csyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 308 "csyrk.f"
/* L140: */
#line 308 "csyrk.f"
		    }
#line 309 "csyrk.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 310 "csyrk.f"
		    i__2 = *n;
#line 310 "csyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 311 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 311 "csyrk.f"
			i__4 = i__ + j * c_dim1;
#line 311 "csyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 311 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 312 "csyrk.f"
/* L150: */
#line 312 "csyrk.f"
		    }
#line 313 "csyrk.f"
		}
#line 314 "csyrk.f"
		i__2 = *k;
#line 314 "csyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 315 "csyrk.f"
		    i__3 = j + l * a_dim1;
#line 315 "csyrk.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 316 "csyrk.f"
			i__3 = j + l * a_dim1;
#line 316 "csyrk.f"
			z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__1.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 316 "csyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 317 "csyrk.f"
			i__3 = *n;
#line 317 "csyrk.f"
			for (i__ = j; i__ <= i__3; ++i__) {
#line 318 "csyrk.f"
			    i__4 = i__ + j * c_dim1;
#line 318 "csyrk.f"
			    i__5 = i__ + j * c_dim1;
#line 318 "csyrk.f"
			    i__6 = i__ + l * a_dim1;
#line 318 "csyrk.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 318 "csyrk.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 318 "csyrk.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 319 "csyrk.f"
/* L160: */
#line 319 "csyrk.f"
			}
#line 320 "csyrk.f"
		    }
#line 321 "csyrk.f"
/* L170: */
#line 321 "csyrk.f"
		}
#line 322 "csyrk.f"
/* L180: */
#line 322 "csyrk.f"
	    }
#line 323 "csyrk.f"
	}
#line 324 "csyrk.f"
    } else {

/*        Form  C := alpha*A**T*A + beta*C. */

#line 328 "csyrk.f"
	if (upper) {
#line 329 "csyrk.f"
	    i__1 = *n;
#line 329 "csyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 330 "csyrk.f"
		i__2 = j;
#line 330 "csyrk.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 331 "csyrk.f"
		    temp.r = 0., temp.i = 0.;
#line 332 "csyrk.f"
		    i__3 = *k;
#line 332 "csyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 333 "csyrk.f"
			i__4 = l + i__ * a_dim1;
#line 333 "csyrk.f"
			i__5 = l + j * a_dim1;
#line 333 "csyrk.f"
			z__2.r = a[i__4].r * a[i__5].r - a[i__4].i * a[i__5]
				.i, z__2.i = a[i__4].r * a[i__5].i + a[i__4]
				.i * a[i__5].r;
#line 333 "csyrk.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 333 "csyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 334 "csyrk.f"
/* L190: */
#line 334 "csyrk.f"
		    }
#line 335 "csyrk.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 336 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 336 "csyrk.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 336 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 337 "csyrk.f"
		    } else {
#line 338 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 338 "csyrk.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 338 "csyrk.f"
			i__4 = i__ + j * c_dim1;
#line 338 "csyrk.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 338 "csyrk.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 338 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 339 "csyrk.f"
		    }
#line 340 "csyrk.f"
/* L200: */
#line 340 "csyrk.f"
		}
#line 341 "csyrk.f"
/* L210: */
#line 341 "csyrk.f"
	    }
#line 342 "csyrk.f"
	} else {
#line 343 "csyrk.f"
	    i__1 = *n;
#line 343 "csyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 344 "csyrk.f"
		i__2 = *n;
#line 344 "csyrk.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 345 "csyrk.f"
		    temp.r = 0., temp.i = 0.;
#line 346 "csyrk.f"
		    i__3 = *k;
#line 346 "csyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 347 "csyrk.f"
			i__4 = l + i__ * a_dim1;
#line 347 "csyrk.f"
			i__5 = l + j * a_dim1;
#line 347 "csyrk.f"
			z__2.r = a[i__4].r * a[i__5].r - a[i__4].i * a[i__5]
				.i, z__2.i = a[i__4].r * a[i__5].i + a[i__4]
				.i * a[i__5].r;
#line 347 "csyrk.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 347 "csyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 348 "csyrk.f"
/* L220: */
#line 348 "csyrk.f"
		    }
#line 349 "csyrk.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 350 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 350 "csyrk.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 350 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 351 "csyrk.f"
		    } else {
#line 352 "csyrk.f"
			i__3 = i__ + j * c_dim1;
#line 352 "csyrk.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 352 "csyrk.f"
			i__4 = i__ + j * c_dim1;
#line 352 "csyrk.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 352 "csyrk.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 352 "csyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 353 "csyrk.f"
		    }
#line 354 "csyrk.f"
/* L230: */
#line 354 "csyrk.f"
		}
#line 355 "csyrk.f"
/* L240: */
#line 355 "csyrk.f"
	    }
#line 356 "csyrk.f"
	}
#line 357 "csyrk.f"
    }

#line 359 "csyrk.f"
    return 0;

/*     End of CSYRK . */

} /* csyrk_ */

