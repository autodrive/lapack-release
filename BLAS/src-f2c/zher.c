#line 1 "zher.f"
/* zher.f -- translated by f2c (version 20100827).
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

#line 1 "zher.f"
/* > \brief \b ZHER */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHER   performs the hermitian rank 1 operation */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of dimension at least */
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
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ). */
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

/* > \date December 2016 */

/* > \ingroup complex16_blas_level2 */

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
/* Subroutine */ int zher_(char *uplo, integer *n, doublereal *alpha, 
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

#line 175 "zher.f"
    /* Parameter adjustments */
#line 175 "zher.f"
    --x;
#line 175 "zher.f"
    a_dim1 = *lda;
#line 175 "zher.f"
    a_offset = 1 + a_dim1;
#line 175 "zher.f"
    a -= a_offset;
#line 175 "zher.f"

#line 175 "zher.f"
    /* Function Body */
#line 175 "zher.f"
    info = 0;
#line 176 "zher.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 177 "zher.f"
	info = 1;
#line 178 "zher.f"
    } else if (*n < 0) {
#line 179 "zher.f"
	info = 2;
#line 180 "zher.f"
    } else if (*incx == 0) {
#line 181 "zher.f"
	info = 5;
#line 182 "zher.f"
    } else if (*lda < max(1,*n)) {
#line 183 "zher.f"
	info = 7;
#line 184 "zher.f"
    }
#line 185 "zher.f"
    if (info != 0) {
#line 186 "zher.f"
	xerbla_("ZHER  ", &info, (ftnlen)6);
#line 187 "zher.f"
	return 0;
#line 188 "zher.f"
    }

/*     Quick return if possible. */

#line 192 "zher.f"
    if (*n == 0 || *alpha == 0.) {
#line 192 "zher.f"
	return 0;
#line 192 "zher.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 196 "zher.f"
    if (*incx <= 0) {
#line 197 "zher.f"
	kx = 1 - (*n - 1) * *incx;
#line 198 "zher.f"
    } else if (*incx != 1) {
#line 199 "zher.f"
	kx = 1;
#line 200 "zher.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 206 "zher.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in upper triangle. */

#line 210 "zher.f"
	if (*incx == 1) {
#line 211 "zher.f"
	    i__1 = *n;
#line 211 "zher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 212 "zher.f"
		i__2 = j;
#line 212 "zher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 213 "zher.f"
		    d_cnjg(&z__2, &x[j]);
#line 213 "zher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 213 "zher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 214 "zher.f"
		    i__2 = j - 1;
#line 214 "zher.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "zher.f"
			i__3 = i__ + j * a_dim1;
#line 215 "zher.f"
			i__4 = i__ + j * a_dim1;
#line 215 "zher.f"
			i__5 = i__;
#line 215 "zher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 215 "zher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 215 "zher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 216 "zher.f"
/* L10: */
#line 216 "zher.f"
		    }
#line 217 "zher.f"
		    i__2 = j + j * a_dim1;
#line 217 "zher.f"
		    i__3 = j + j * a_dim1;
#line 217 "zher.f"
		    i__4 = j;
#line 217 "zher.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 217 "zher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 217 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 218 "zher.f"
		} else {
#line 219 "zher.f"
		    i__2 = j + j * a_dim1;
#line 219 "zher.f"
		    i__3 = j + j * a_dim1;
#line 219 "zher.f"
		    d__1 = a[i__3].r;
#line 219 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 220 "zher.f"
		}
#line 221 "zher.f"
/* L20: */
#line 221 "zher.f"
	    }
#line 222 "zher.f"
	} else {
#line 223 "zher.f"
	    jx = kx;
#line 224 "zher.f"
	    i__1 = *n;
#line 224 "zher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 225 "zher.f"
		i__2 = jx;
#line 225 "zher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 226 "zher.f"
		    d_cnjg(&z__2, &x[jx]);
#line 226 "zher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 226 "zher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 227 "zher.f"
		    ix = kx;
#line 228 "zher.f"
		    i__2 = j - 1;
#line 228 "zher.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 229 "zher.f"
			i__3 = i__ + j * a_dim1;
#line 229 "zher.f"
			i__4 = i__ + j * a_dim1;
#line 229 "zher.f"
			i__5 = ix;
#line 229 "zher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 229 "zher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 229 "zher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 230 "zher.f"
			ix += *incx;
#line 231 "zher.f"
/* L30: */
#line 231 "zher.f"
		    }
#line 232 "zher.f"
		    i__2 = j + j * a_dim1;
#line 232 "zher.f"
		    i__3 = j + j * a_dim1;
#line 232 "zher.f"
		    i__4 = jx;
#line 232 "zher.f"
		    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i, z__1.i =
			     x[i__4].r * temp.i + x[i__4].i * temp.r;
#line 232 "zher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 232 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 233 "zher.f"
		} else {
#line 234 "zher.f"
		    i__2 = j + j * a_dim1;
#line 234 "zher.f"
		    i__3 = j + j * a_dim1;
#line 234 "zher.f"
		    d__1 = a[i__3].r;
#line 234 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 235 "zher.f"
		}
#line 236 "zher.f"
		jx += *incx;
#line 237 "zher.f"
/* L40: */
#line 237 "zher.f"
	    }
#line 238 "zher.f"
	}
#line 239 "zher.f"
    } else {

/*        Form  A  when A is stored in lower triangle. */

#line 243 "zher.f"
	if (*incx == 1) {
#line 244 "zher.f"
	    i__1 = *n;
#line 244 "zher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "zher.f"
		i__2 = j;
#line 245 "zher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 246 "zher.f"
		    d_cnjg(&z__2, &x[j]);
#line 246 "zher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 246 "zher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 247 "zher.f"
		    i__2 = j + j * a_dim1;
#line 247 "zher.f"
		    i__3 = j + j * a_dim1;
#line 247 "zher.f"
		    i__4 = j;
#line 247 "zher.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 247 "zher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 247 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 248 "zher.f"
		    i__2 = *n;
#line 248 "zher.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 249 "zher.f"
			i__3 = i__ + j * a_dim1;
#line 249 "zher.f"
			i__4 = i__ + j * a_dim1;
#line 249 "zher.f"
			i__5 = i__;
#line 249 "zher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 249 "zher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 249 "zher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 250 "zher.f"
/* L50: */
#line 250 "zher.f"
		    }
#line 251 "zher.f"
		} else {
#line 252 "zher.f"
		    i__2 = j + j * a_dim1;
#line 252 "zher.f"
		    i__3 = j + j * a_dim1;
#line 252 "zher.f"
		    d__1 = a[i__3].r;
#line 252 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 253 "zher.f"
		}
#line 254 "zher.f"
/* L60: */
#line 254 "zher.f"
	    }
#line 255 "zher.f"
	} else {
#line 256 "zher.f"
	    jx = kx;
#line 257 "zher.f"
	    i__1 = *n;
#line 257 "zher.f"
	    for (j = 1; j <= i__1; ++j) {
#line 258 "zher.f"
		i__2 = jx;
#line 258 "zher.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 259 "zher.f"
		    d_cnjg(&z__2, &x[jx]);
#line 259 "zher.f"
		    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 259 "zher.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 260 "zher.f"
		    i__2 = j + j * a_dim1;
#line 260 "zher.f"
		    i__3 = j + j * a_dim1;
#line 260 "zher.f"
		    i__4 = jx;
#line 260 "zher.f"
		    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i, z__1.i =
			     temp.r * x[i__4].i + temp.i * x[i__4].r;
#line 260 "zher.f"
		    d__1 = a[i__3].r + z__1.r;
#line 260 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 261 "zher.f"
		    ix = jx;
#line 262 "zher.f"
		    i__2 = *n;
#line 262 "zher.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 263 "zher.f"
			ix += *incx;
#line 264 "zher.f"
			i__3 = i__ + j * a_dim1;
#line 264 "zher.f"
			i__4 = i__ + j * a_dim1;
#line 264 "zher.f"
			i__5 = ix;
#line 264 "zher.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 264 "zher.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 264 "zher.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 265 "zher.f"
/* L70: */
#line 265 "zher.f"
		    }
#line 266 "zher.f"
		} else {
#line 267 "zher.f"
		    i__2 = j + j * a_dim1;
#line 267 "zher.f"
		    i__3 = j + j * a_dim1;
#line 267 "zher.f"
		    d__1 = a[i__3].r;
#line 267 "zher.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 268 "zher.f"
		}
#line 269 "zher.f"
		jx += *incx;
#line 270 "zher.f"
/* L80: */
#line 270 "zher.f"
	    }
#line 271 "zher.f"
	}
#line 272 "zher.f"
    }

#line 274 "zher.f"
    return 0;

/*     End of ZHER  . */

} /* zher_ */

