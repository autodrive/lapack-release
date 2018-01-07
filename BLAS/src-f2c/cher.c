#line 1 "cher.f"
/* cher.f -- translated by f2c (version 20100827).
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

#line 1 "cher.f"
/* > \brief \b CHER */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHER   performs the hermitian rank 1 operation */
/* > */
/* >    A := alpha*x*x**H + A, */
/* > */
/* > where alpha is a real scalar, x is an n element vector and A is an */
/* > n by n hermitian matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array A is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/* >                                  is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the hermitian matrix and the strictly */
/* >           lower triangular part of A is not referenced. On exit, the */
/* >           upper triangular part of the array A is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the hermitian matrix and the strictly */
/* >           upper triangular part of A is not referenced. On exit, the */
/* >           lower triangular part of the array A is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set, they are assumed to be zero, and on exit they */
/* >           are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* Subroutine */ int cher_(char *uplo, integer *n, doublereal *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *a, integer *lda, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 175 "cher.f"
    /* Parameter adjustments */
#line 175 "cher.f"
    --x;
#line 175 "cher.f"
    a_dim1 = *lda;
#line 175 "cher.f"
    a_offset = 1 + a_dim1;
#line 175 "cher.f"
    a -= a_offset;
#line 175 "cher.f"

#line 175 "cher.f"
    /* Function Body */
#line 175 "cher.f"
    info = 0;
#line 176 "cher.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 177 "cher.f"
	info = 1;
#line 178 "cher.f"
    } else if (*n < 0) {
#line 179 "cher.f"
	info = 2;
#line 180 "cher.f"
    } else if (*incx == 0) {
#line 181 "cher.f"
	info = 5;
#line 182 "cher.f"
    } else if (*lda < max(1,*n)) {
#line 183 "cher.f"
	info = 7;
#line 184 "cher.f"
    }
#line 185 "cher.f"
    if (info != 0) {
#line 186 "cher.f"
	xerbla_("CHER  ", &info, (ftnlen)6);
#line 187 "cher.f"
	return 0;
#line 188 "cher.f"
    }

/*     Quick return if possible. */

#line 192 "cher.f"
    if (*n == 0 || *alpha == 0.) {
#line 192 "cher.f"
	return 0;
#line 192 "cher.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 196 "cher.f"
    if (*incx <= 0) {
#line 197 "cher.f"
	kx = 1 - (*n - 1) * *incx;
#line 198 "cher.f"
    } else if (*incx != 1) {
#line 199 "cher.f"
	kx = 1;
#line 200 "cher.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 206 "cher.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in upper triangle. */

#line 210 "cher.f"
	if (*incx == 1) {
#line 211 "cher.f"
	    i__1 = *n;
#line 211 "cher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 212 "cher.f"
		i__2 = j;
#line 212 "cher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 213 "cher.f"
		    d_cnjg(&z__2, &x[j]);
#line 213 "cher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 213 "cher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 214 "cher.f"
		    i__2 = j - 1;
#line 214 "cher.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "cher.f"
			i__3 = i__ + j * a_dim1;
#line 215 "cher.f"
			i__4 = i__ + j * a_dim1;
#line 215 "cher.f"
			i__5 = i__;
#line 215 "cher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 215 "cher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 215 "cher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 216 "cher.f"
/* L10: */
#line 216 "cher.f"
		    }
#line 217 "cher.f"
		    i__2 = j + j * a_dim1;
#line 217 "cher.f"
		    i__3 = j + j * a_dim1;
#line 217 "cher.f"
		    i__4 = j;
#line 217 "cher.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 217 "cher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 217 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 218 "cher.f"
		} else {
#line 219 "cher.f"
		    i__2 = j + j * a_dim1;
#line 219 "cher.f"
		    i__3 = j + j * a_dim1;
#line 219 "cher.f"
		    d__1 = a[i__3].r;
#line 219 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 220 "cher.f"
		}
#line 221 "cher.f"
/* L20: */
#line 221 "cher.f"
	    }
#line 222 "cher.f"
	} else {
#line 223 "cher.f"
	    jx = kx;
#line 224 "cher.f"
	    i__1 = *n;
#line 224 "cher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 225 "cher.f"
		i__2 = jx;
#line 225 "cher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 226 "cher.f"
		    d_cnjg(&z__2, &x[jx]);
#line 226 "cher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 226 "cher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 227 "cher.f"
		    ix = kx;
#line 228 "cher.f"
		    i__2 = j - 1;
#line 228 "cher.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 229 "cher.f"
			i__3 = i__ + j * a_dim1;
#line 229 "cher.f"
			i__4 = i__ + j * a_dim1;
#line 229 "cher.f"
			i__5 = ix;
#line 229 "cher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 229 "cher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 229 "cher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 230 "cher.f"
			ix += *incx;
#line 231 "cher.f"
/* L30: */
#line 231 "cher.f"
		    }
#line 232 "cher.f"
		    i__2 = j + j * a_dim1;
#line 232 "cher.f"
		    i__3 = j + j * a_dim1;
#line 232 "cher.f"
		    i__4 = jx;
#line 232 "cher.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 232 "cher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 232 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 233 "cher.f"
		} else {
#line 234 "cher.f"
		    i__2 = j + j * a_dim1;
#line 234 "cher.f"
		    i__3 = j + j * a_dim1;
#line 234 "cher.f"
		    d__1 = a[i__3].r;
#line 234 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 235 "cher.f"
		}
#line 236 "cher.f"
		jx += *incx;
#line 237 "cher.f"
/* L40: */
#line 237 "cher.f"
	    }
#line 238 "cher.f"
	}
#line 239 "cher.f"
    } else {

/*        Form  A  when A is stored in lower triangle. */

#line 243 "cher.f"
	if (*incx == 1) {
#line 244 "cher.f"
	    i__1 = *n;
#line 244 "cher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "cher.f"
		i__2 = j;
#line 245 "cher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 246 "cher.f"
		    d_cnjg(&z__2, &x[j]);
#line 246 "cher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 246 "cher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 247 "cher.f"
		    i__2 = j + j * a_dim1;
#line 247 "cher.f"
		    i__3 = j + j * a_dim1;
#line 247 "cher.f"
		    i__4 = j;
#line 247 "cher.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 247 "cher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 247 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 248 "cher.f"
		    i__2 = *n;
#line 248 "cher.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 249 "cher.f"
			i__3 = i__ + j * a_dim1;
#line 249 "cher.f"
			i__4 = i__ + j * a_dim1;
#line 249 "cher.f"
			i__5 = i__;
#line 249 "cher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 249 "cher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 249 "cher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 250 "cher.f"
/* L50: */
#line 250 "cher.f"
		    }
#line 251 "cher.f"
		} else {
#line 252 "cher.f"
		    i__2 = j + j * a_dim1;
#line 252 "cher.f"
		    i__3 = j + j * a_dim1;
#line 252 "cher.f"
		    d__1 = a[i__3].r;
#line 252 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 253 "cher.f"
		}
#line 254 "cher.f"
/* L60: */
#line 254 "cher.f"
	    }
#line 255 "cher.f"
	} else {
#line 256 "cher.f"
	    jx = kx;
#line 257 "cher.f"
	    i__1 = *n;
#line 257 "cher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "cher.f"
		i__2 = jx;
#line 258 "cher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 259 "cher.f"
		    d_cnjg(&z__2, &x[jx]);
#line 259 "cher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 259 "cher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 260 "cher.f"
		    i__2 = j + j * a_dim1;
#line 260 "cher.f"
		    i__3 = j + j * a_dim1;
#line 260 "cher.f"
		    i__4 = jx;
#line 260 "cher.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 260 "cher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 260 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 261 "cher.f"
		    ix = jx;
#line 262 "cher.f"
		    i__2 = *n;
#line 262 "cher.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 263 "cher.f"
			ix += *incx;
#line 264 "cher.f"
			i__3 = i__ + j * a_dim1;
#line 264 "cher.f"
			i__4 = i__ + j * a_dim1;
#line 264 "cher.f"
			i__5 = ix;
#line 264 "cher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 264 "cher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 264 "cher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 265 "cher.f"
/* L70: */
#line 265 "cher.f"
		    }
#line 266 "cher.f"
		} else {
#line 267 "cher.f"
		    i__2 = j + j * a_dim1;
#line 267 "cher.f"
		    i__3 = j + j * a_dim1;
#line 267 "cher.f"
		    d__1 = a[i__3].r;
#line 267 "cher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 268 "cher.f"
		}
#line 269 "cher.f"
		jx += *incx;
#line 270 "cher.f"
/* L80: */
#line 270 "cher.f"
	    }
#line 271 "cher.f"
	}
#line 272 "cher.f"
    }

#line 274 "cher.f"
    return 0;

/*     End of CHER  . */

} /* cher_ */

