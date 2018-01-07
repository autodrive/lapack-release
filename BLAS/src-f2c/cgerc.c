#line 1 "cgerc.f"
/* cgerc.f -- translated by f2c (version 20100827).
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

#line 1 "cgerc.f"
/* > \brief \b CGERC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGERC  performs the rank 1 operation */
/* > */
/* >    A := alpha*x*y**H + A, */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array of dimension at least */
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
/* >          Y is COMPLEX array of dimension at least */
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
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
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
/* Subroutine */ int cgerc_(integer *m, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
	doublecomplex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ix, jy, kx, info;
    static doublecomplex temp;
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

#line 165 "cgerc.f"
    /* Parameter adjustments */
#line 165 "cgerc.f"
    --x;
#line 165 "cgerc.f"
    --y;
#line 165 "cgerc.f"
    a_dim1 = *lda;
#line 165 "cgerc.f"
    a_offset = 1 + a_dim1;
#line 165 "cgerc.f"
    a -= a_offset;
#line 165 "cgerc.f"

#line 165 "cgerc.f"
    /* Function Body */
#line 165 "cgerc.f"
    info = 0;
#line 166 "cgerc.f"
    if (*m < 0) {
#line 167 "cgerc.f"
	info = 1;
#line 168 "cgerc.f"
    } else if (*n < 0) {
#line 169 "cgerc.f"
	info = 2;
#line 170 "cgerc.f"
    } else if (*incx == 0) {
#line 171 "cgerc.f"
	info = 5;
#line 172 "cgerc.f"
    } else if (*incy == 0) {
#line 173 "cgerc.f"
	info = 7;
#line 174 "cgerc.f"
    } else if (*lda < max(1,*m)) {
#line 175 "cgerc.f"
	info = 9;
#line 176 "cgerc.f"
    }
#line 177 "cgerc.f"
    if (info != 0) {
#line 178 "cgerc.f"
	xerbla_("CGERC ", &info, (ftnlen)6);
#line 179 "cgerc.f"
	return 0;
#line 180 "cgerc.f"
    }

/*     Quick return if possible. */

#line 184 "cgerc.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0.) {
#line 184 "cgerc.f"
	return 0;
#line 184 "cgerc.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 189 "cgerc.f"
    if (*incy > 0) {
#line 190 "cgerc.f"
	jy = 1;
#line 191 "cgerc.f"
    } else {
#line 192 "cgerc.f"
	jy = 1 - (*n - 1) * *incy;
#line 193 "cgerc.f"
    }
#line 194 "cgerc.f"
    if (*incx == 1) {
#line 195 "cgerc.f"
	i__1 = *n;
#line 195 "cgerc.f"
	for (j = 1; j <= i__1; ++j) {
#line 196 "cgerc.f"
	    i__2 = jy;
#line 196 "cgerc.f"
	    if (y[i__2].r != 0. || y[i__2].i != 0.) {
#line 197 "cgerc.f"
		d_cnjg(&z__2, &y[jy]);
#line 197 "cgerc.f"
		z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			alpha->r * z__2.i + alpha->i * z__2.r;
#line 197 "cgerc.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 198 "cgerc.f"
		i__2 = *m;
#line 198 "cgerc.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 199 "cgerc.f"
		    i__3 = i__ + j * a_dim1;
#line 199 "cgerc.f"
		    i__4 = i__ + j * a_dim1;
#line 199 "cgerc.f"
		    i__5 = i__;
#line 199 "cgerc.f"
		    z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, z__2.i =
			     x[i__5].r * temp.i + x[i__5].i * temp.r;
#line 199 "cgerc.f"
		    z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
#line 199 "cgerc.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 200 "cgerc.f"
/* L10: */
#line 200 "cgerc.f"
		}
#line 201 "cgerc.f"
	    }
#line 202 "cgerc.f"
	    jy += *incy;
#line 203 "cgerc.f"
/* L20: */
#line 203 "cgerc.f"
	}
#line 204 "cgerc.f"
    } else {
#line 205 "cgerc.f"
	if (*incx > 0) {
#line 206 "cgerc.f"
	    kx = 1;
#line 207 "cgerc.f"
	} else {
#line 208 "cgerc.f"
	    kx = 1 - (*m - 1) * *incx;
#line 209 "cgerc.f"
	}
#line 210 "cgerc.f"
	i__1 = *n;
#line 210 "cgerc.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "cgerc.f"
	    i__2 = jy;
#line 211 "cgerc.f"
	    if (y[i__2].r != 0. || y[i__2].i != 0.) {
#line 212 "cgerc.f"
		d_cnjg(&z__2, &y[jy]);
#line 212 "cgerc.f"
		z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			alpha->r * z__2.i + alpha->i * z__2.r;
#line 212 "cgerc.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 213 "cgerc.f"
		ix = kx;
#line 214 "cgerc.f"
		i__2 = *m;
#line 214 "cgerc.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "cgerc.f"
		    i__3 = i__ + j * a_dim1;
#line 215 "cgerc.f"
		    i__4 = i__ + j * a_dim1;
#line 215 "cgerc.f"
		    i__5 = ix;
#line 215 "cgerc.f"
		    z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i, z__2.i =
			     x[i__5].r * temp.i + x[i__5].i * temp.r;
#line 215 "cgerc.f"
		    z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
#line 215 "cgerc.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 216 "cgerc.f"
		    ix += *incx;
#line 217 "cgerc.f"
/* L30: */
#line 217 "cgerc.f"
		}
#line 218 "cgerc.f"
	    }
#line 219 "cgerc.f"
	    jy += *incy;
#line 220 "cgerc.f"
/* L40: */
#line 220 "cgerc.f"
	}
#line 221 "cgerc.f"
    }

#line 223 "cgerc.f"
    return 0;

/*     End of CGERC . */

} /* cgerc_ */

