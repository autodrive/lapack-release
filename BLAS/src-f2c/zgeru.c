#line 1 "zgeru.f"
/* zgeru.f -- translated by f2c (version 20100827).
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

#line 1 "zgeru.f"
/* > \brief \b ZGERU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGERU  performs the rank 1 operation */
/* > */
/* >    A := alpha*x*y**T + A, */
/* > */
/* > where alpha is a scalar, x is an m element vector, y is an n element */
/* > vector and A is an m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of dimension at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the m */
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
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ). */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. On exit, A is */
/* >           overwritten by the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, m ). */
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
/* Subroutine */ int zgeru_(integer *m, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
	doublecomplex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, ix, jy, kx, info;
    static doublecomplex temp;
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

#line 165 "zgeru.f"
    /* Parameter adjustments */
#line 165 "zgeru.f"
    --x;
#line 165 "zgeru.f"
    --y;
#line 165 "zgeru.f"
    a_dim1 = *lda;
#line 165 "zgeru.f"
    a_offset = 1 + a_dim1;
#line 165 "zgeru.f"
    a -= a_offset;
#line 165 "zgeru.f"

#line 165 "zgeru.f"
    /* Function Body */
#line 165 "zgeru.f"
    info = 0;
#line 166 "zgeru.f"
    if (*m < 0) {
#line 167 "zgeru.f"
	info = 1;
#line 168 "zgeru.f"
    } else if (*n < 0) {
#line 169 "zgeru.f"
	info = 2;
#line 170 "zgeru.f"
    } else if (*incx == 0) {
#line 171 "zgeru.f"
	info = 5;
#line 172 "zgeru.f"
    } else if (*incy == 0) {
#line 173 "zgeru.f"
	info = 7;
#line 174 "zgeru.f"
    } else if (*lda < max(1,*m)) {
#line 175 "zgeru.f"
	info = 9;
#line 176 "zgeru.f"
    }
#line 177 "zgeru.f"
    if (info != 0) {
#line 178 "zgeru.f"
	xerbla_("ZGERU ", &info, (ftnlen)6);
#line 179 "zgeru.f"
	return 0;
#line 180 "zgeru.f"
    }

/*     Quick return if possible. */

#line 184 "zgeru.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0.) {
#line 184 "zgeru.f"
	return 0;
#line 184 "zgeru.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 189 "zgeru.f"
    if (*incy > 0) {
#line 190 "zgeru.f"
	jy = 1;
#line 191 "zgeru.f"
    } else {
#line 192 "zgeru.f"
	jy = 1 - (*n - 1) * *incy;
#line 193 "zgeru.f"
    }
#line 194 "zgeru.f"
    if (*incx == 1) {
#line 195 "zgeru.f"
	i__1 = *n;
#line 195 "zgeru.f"
	for (j = 1; j <= i__1; ++j) {
#line 196 "zgeru.f"
	    i__2 = jy;
#line 196 "zgeru.f"
	    if (y[i__2].r != 0. || y[i__2].i != 0.) {
#line 197 "zgeru.f"
		i__2 = jy;
#line 197 "zgeru.f"
		z__1.r = alpha->r * y[i__2].r - alpha->i * y[i__2].i, z__1.i =
			 alpha->r * y[i__2].i + alpha->i * y[i__2].r;
#line 197 "zgeru.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 198 "zgeru.f"
		i__2 = *m;
#line 198 "zgeru.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 199 "zgeru.f"
		    i__3 = i__ + j * a_dim1;
#line 199 "zgeru.f"
		    i__4 = i__ + j * a_dim1;
#line 199 "zgeru.f"
		    i__5 = i__;
#line 199 "zgeru.f"
		    z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, z__2.i =
			     x[i__5].r * temp.i + x[i__5].i * temp.r;
#line 199 "zgeru.f"
		    z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
#line 199 "zgeru.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 200 "zgeru.f"
/* L10: */
#line 200 "zgeru.f"
		}
#line 201 "zgeru.f"
	    }
#line 202 "zgeru.f"
	    jy += *incy;
#line 203 "zgeru.f"
/* L20: */
#line 203 "zgeru.f"
	}
#line 204 "zgeru.f"
    } else {
#line 205 "zgeru.f"
	if (*incx > 0) {
#line 206 "zgeru.f"
	    kx = 1;
#line 207 "zgeru.f"
	} else {
#line 208 "zgeru.f"
	    kx = 1 - (*m - 1) * *incx;
#line 209 "zgeru.f"
	}
#line 210 "zgeru.f"
	i__1 = *n;
#line 210 "zgeru.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "zgeru.f"
	    i__2 = jy;
#line 211 "zgeru.f"
	    if (y[i__2].r != 0. || y[i__2].i != 0.) {
#line 212 "zgeru.f"
		i__2 = jy;
#line 212 "zgeru.f"
		z__1.r = alpha->r * y[i__2].r - alpha->i * y[i__2].i, z__1.i =
			 alpha->r * y[i__2].i + alpha->i * y[i__2].r;
#line 212 "zgeru.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 213 "zgeru.f"
		ix = kx;
#line 214 "zgeru.f"
		i__2 = *m;
#line 214 "zgeru.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "zgeru.f"
		    i__3 = i__ + j * a_dim1;
#line 215 "zgeru.f"
		    i__4 = i__ + j * a_dim1;
#line 215 "zgeru.f"
		    i__5 = ix;
#line 215 "zgeru.f"
		    z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, z__2.i =
			     x[i__5].r * temp.i + x[i__5].i * temp.r;
#line 215 "zgeru.f"
		    z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
#line 215 "zgeru.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 216 "zgeru.f"
		    ix += *incx;
#line 217 "zgeru.f"
/* L30: */
#line 217 "zgeru.f"
		}
#line 218 "zgeru.f"
	    }
#line 219 "zgeru.f"
	    jy += *incy;
#line 220 "zgeru.f"
/* L40: */
#line 220 "zgeru.f"
	}
#line 221 "zgeru.f"
    }

#line 223 "zgeru.f"
    return 0;

/*     End of ZGERU . */

} /* zgeru_ */

