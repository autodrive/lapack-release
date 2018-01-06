#line 1 "zsyr.f"
/* zsyr.f -- translated by f2c (version 20100827).
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

#line 1 "zsyr.f"
/* > \brief \b ZSYR performs the symmetric rank-1 update of a complex symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyr.f"
> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyr.f"
> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyr.f"
> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCX, LDA, N */
/*       COMPLEX*16         ALPHA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYR   performs the symmetric rank 1 operation */
/* > */
/* >    A := alpha*x*x**H + A, */
/* > */
/* > where alpha is a complex scalar, x is an n element vector and A is an */
/* > n by n symmetric matrix. */
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
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the N- */
/* >           element vector x. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
/* >           Before entry, with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. On exit, the */
/* >           upper triangular part of the array A is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry, with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           upper triangular part of A is not referenced. On exit, the */
/* >           lower triangular part of the array A is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, N ). */
/* >           Unchanged on exit. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16SYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zsyr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *a, integer *lda, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 176 "zsyr.f"
    /* Parameter adjustments */
#line 176 "zsyr.f"
    --x;
#line 176 "zsyr.f"
    a_dim1 = *lda;
#line 176 "zsyr.f"
    a_offset = 1 + a_dim1;
#line 176 "zsyr.f"
    a -= a_offset;
#line 176 "zsyr.f"

#line 176 "zsyr.f"
    /* Function Body */
#line 176 "zsyr.f"
    info = 0;
#line 177 "zsyr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 178 "zsyr.f"
	info = 1;
#line 179 "zsyr.f"
    } else if (*n < 0) {
#line 180 "zsyr.f"
	info = 2;
#line 181 "zsyr.f"
    } else if (*incx == 0) {
#line 182 "zsyr.f"
	info = 5;
#line 183 "zsyr.f"
    } else if (*lda < max(1,*n)) {
#line 184 "zsyr.f"
	info = 7;
#line 185 "zsyr.f"
    }
#line 186 "zsyr.f"
    if (info != 0) {
#line 187 "zsyr.f"
	xerbla_("ZSYR  ", &info, (ftnlen)6);
#line 188 "zsyr.f"
	return 0;
#line 189 "zsyr.f"
    }

/*     Quick return if possible. */

#line 193 "zsyr.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
#line 193 "zsyr.f"
	return 0;
#line 193 "zsyr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 198 "zsyr.f"
    if (*incx <= 0) {
#line 199 "zsyr.f"
	kx = 1 - (*n - 1) * *incx;
#line 200 "zsyr.f"
    } else if (*incx != 1) {
#line 201 "zsyr.f"
	kx = 1;
#line 202 "zsyr.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 208 "zsyr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in upper triangle. */

#line 212 "zsyr.f"
	if (*incx == 1) {
#line 213 "zsyr.f"
	    i__1 = *n;
#line 213 "zsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 214 "zsyr.f"
		i__2 = j;
#line 214 "zsyr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 215 "zsyr.f"
		    i__2 = j;
#line 215 "zsyr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 215 "zsyr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 216 "zsyr.f"
		    i__2 = j;
#line 216 "zsyr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 217 "zsyr.f"
			i__3 = i__ + j * a_dim1;
#line 217 "zsyr.f"
			i__4 = i__ + j * a_dim1;
#line 217 "zsyr.f"
			i__5 = i__;
#line 217 "zsyr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 217 "zsyr.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 217 "zsyr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 218 "zsyr.f"
/* L10: */
#line 218 "zsyr.f"
		    }
#line 219 "zsyr.f"
		}
#line 220 "zsyr.f"
/* L20: */
#line 220 "zsyr.f"
	    }
#line 221 "zsyr.f"
	} else {
#line 222 "zsyr.f"
	    jx = kx;
#line 223 "zsyr.f"
	    i__1 = *n;
#line 223 "zsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 224 "zsyr.f"
		i__2 = jx;
#line 224 "zsyr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 225 "zsyr.f"
		    i__2 = jx;
#line 225 "zsyr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 225 "zsyr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 226 "zsyr.f"
		    ix = kx;
#line 227 "zsyr.f"
		    i__2 = j;
#line 227 "zsyr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 228 "zsyr.f"
			i__3 = i__ + j * a_dim1;
#line 228 "zsyr.f"
			i__4 = i__ + j * a_dim1;
#line 228 "zsyr.f"
			i__5 = ix;
#line 228 "zsyr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 228 "zsyr.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 228 "zsyr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 229 "zsyr.f"
			ix += *incx;
#line 230 "zsyr.f"
/* L30: */
#line 230 "zsyr.f"
		    }
#line 231 "zsyr.f"
		}
#line 232 "zsyr.f"
		jx += *incx;
#line 233 "zsyr.f"
/* L40: */
#line 233 "zsyr.f"
	    }
#line 234 "zsyr.f"
	}
#line 235 "zsyr.f"
    } else {

/*        Form  A  when A is stored in lower triangle. */

#line 239 "zsyr.f"
	if (*incx == 1) {
#line 240 "zsyr.f"
	    i__1 = *n;
#line 240 "zsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 241 "zsyr.f"
		i__2 = j;
#line 241 "zsyr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 242 "zsyr.f"
		    i__2 = j;
#line 242 "zsyr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 242 "zsyr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 243 "zsyr.f"
		    i__2 = *n;
#line 243 "zsyr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 244 "zsyr.f"
			i__3 = i__ + j * a_dim1;
#line 244 "zsyr.f"
			i__4 = i__ + j * a_dim1;
#line 244 "zsyr.f"
			i__5 = i__;
#line 244 "zsyr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 244 "zsyr.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 244 "zsyr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 245 "zsyr.f"
/* L50: */
#line 245 "zsyr.f"
		    }
#line 246 "zsyr.f"
		}
#line 247 "zsyr.f"
/* L60: */
#line 247 "zsyr.f"
	    }
#line 248 "zsyr.f"
	} else {
#line 249 "zsyr.f"
	    jx = kx;
#line 250 "zsyr.f"
	    i__1 = *n;
#line 250 "zsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 251 "zsyr.f"
		i__2 = jx;
#line 251 "zsyr.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 252 "zsyr.f"
		    i__2 = jx;
#line 252 "zsyr.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 252 "zsyr.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 253 "zsyr.f"
		    ix = jx;
#line 254 "zsyr.f"
		    i__2 = *n;
#line 254 "zsyr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 255 "zsyr.f"
			i__3 = i__ + j * a_dim1;
#line 255 "zsyr.f"
			i__4 = i__ + j * a_dim1;
#line 255 "zsyr.f"
			i__5 = ix;
#line 255 "zsyr.f"
			z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, 
				z__2.i = x[i__5].r * temp.i + x[i__5].i * 
				temp.r;
#line 255 "zsyr.f"
			z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + 
				z__2.i;
#line 255 "zsyr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 256 "zsyr.f"
			ix += *incx;
#line 257 "zsyr.f"
/* L70: */
#line 257 "zsyr.f"
		    }
#line 258 "zsyr.f"
		}
#line 259 "zsyr.f"
		jx += *incx;
#line 260 "zsyr.f"
/* L80: */
#line 260 "zsyr.f"
	    }
#line 261 "zsyr.f"
	}
#line 262 "zsyr.f"
    }

#line 264 "zsyr.f"
    return 0;

/*     End of ZSYR */

} /* zsyr_ */

