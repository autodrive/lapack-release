#line 1 "dger.f"
/* dger.f -- translated by f2c (version 20100827).
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

#line 1 "dger.f"
/* > \brief \b DGER */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGER   performs the rank 1 operation */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension at least */
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
/* >          Y is DOUBLE PRECISION array, dimension at least */
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
/* >          A is DOUBLE PRECISION array, dimension ( LDA, N ) */
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

/* > \ingroup double_blas_level2 */

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
/* Subroutine */ int dger_(integer *m, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jy, kx, info;
    static doublereal temp;
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

#line 165 "dger.f"
    /* Parameter adjustments */
#line 165 "dger.f"
    --x;
#line 165 "dger.f"
    --y;
#line 165 "dger.f"
    a_dim1 = *lda;
#line 165 "dger.f"
    a_offset = 1 + a_dim1;
#line 165 "dger.f"
    a -= a_offset;
#line 165 "dger.f"

#line 165 "dger.f"
    /* Function Body */
#line 165 "dger.f"
    info = 0;
#line 166 "dger.f"
    if (*m < 0) {
#line 167 "dger.f"
	info = 1;
#line 168 "dger.f"
    } else if (*n < 0) {
#line 169 "dger.f"
	info = 2;
#line 170 "dger.f"
    } else if (*incx == 0) {
#line 171 "dger.f"
	info = 5;
#line 172 "dger.f"
    } else if (*incy == 0) {
#line 173 "dger.f"
	info = 7;
#line 174 "dger.f"
    } else if (*lda < max(1,*m)) {
#line 175 "dger.f"
	info = 9;
#line 176 "dger.f"
    }
#line 177 "dger.f"
    if (info != 0) {
#line 178 "dger.f"
	xerbla_("DGER  ", &info, (ftnlen)6);
#line 179 "dger.f"
	return 0;
#line 180 "dger.f"
    }

/*     Quick return if possible. */

#line 184 "dger.f"
    if (*m == 0 || *n == 0 || *alpha == 0.) {
#line 184 "dger.f"
	return 0;
#line 184 "dger.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 189 "dger.f"
    if (*incy > 0) {
#line 190 "dger.f"
	jy = 1;
#line 191 "dger.f"
    } else {
#line 192 "dger.f"
	jy = 1 - (*n - 1) * *incy;
#line 193 "dger.f"
    }
#line 194 "dger.f"
    if (*incx == 1) {
#line 195 "dger.f"
	i__1 = *n;
#line 195 "dger.f"
	for (j = 1; j <= i__1; ++j) {
#line 196 "dger.f"
	    if (y[jy] != 0.) {
#line 197 "dger.f"
		temp = *alpha * y[jy];
#line 198 "dger.f"
		i__2 = *m;
#line 198 "dger.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 199 "dger.f"
		    a[i__ + j * a_dim1] += x[i__] * temp;
#line 200 "dger.f"
/* L10: */
#line 200 "dger.f"
		}
#line 201 "dger.f"
	    }
#line 202 "dger.f"
	    jy += *incy;
#line 203 "dger.f"
/* L20: */
#line 203 "dger.f"
	}
#line 204 "dger.f"
    } else {
#line 205 "dger.f"
	if (*incx > 0) {
#line 206 "dger.f"
	    kx = 1;
#line 207 "dger.f"
	} else {
#line 208 "dger.f"
	    kx = 1 - (*m - 1) * *incx;
#line 209 "dger.f"
	}
#line 210 "dger.f"
	i__1 = *n;
#line 210 "dger.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "dger.f"
	    if (y[jy] != 0.) {
#line 212 "dger.f"
		temp = *alpha * y[jy];
#line 213 "dger.f"
		ix = kx;
#line 214 "dger.f"
		i__2 = *m;
#line 214 "dger.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "dger.f"
		    a[i__ + j * a_dim1] += x[ix] * temp;
#line 216 "dger.f"
		    ix += *incx;
#line 217 "dger.f"
/* L30: */
#line 217 "dger.f"
		}
#line 218 "dger.f"
	    }
#line 219 "dger.f"
	    jy += *incy;
#line 220 "dger.f"
/* L40: */
#line 220 "dger.f"
	}
#line 221 "dger.f"
    }

#line 223 "dger.f"
    return 0;

/*     End of DGER  . */

} /* dger_ */

