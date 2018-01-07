#line 1 "dsyr.f"
/* dsyr.f -- translated by f2c (version 20100827).
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

#line 1 "dsyr.f"
/* > \brief \b DSYR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYR   performs the symmetric rank 1 operation */
/* > */
/* >    A := alpha*x*x**T + A, */
/* > */
/* > where alpha is a real scalar, x is an n element vector and A is an */
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
/* >          X is DOUBLE PRECISION array of dimension at least */
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
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. On exit, the */
/* >           upper triangular part of the array A is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
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
/* >           max( 1, n ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* Subroutine */ int dsyr_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *a, integer *lda, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static doublereal temp;
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

#line 172 "dsyr.f"
    /* Parameter adjustments */
#line 172 "dsyr.f"
    --x;
#line 172 "dsyr.f"
    a_dim1 = *lda;
#line 172 "dsyr.f"
    a_offset = 1 + a_dim1;
#line 172 "dsyr.f"
    a -= a_offset;
#line 172 "dsyr.f"

#line 172 "dsyr.f"
    /* Function Body */
#line 172 "dsyr.f"
    info = 0;
#line 173 "dsyr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 174 "dsyr.f"
	info = 1;
#line 175 "dsyr.f"
    } else if (*n < 0) {
#line 176 "dsyr.f"
	info = 2;
#line 177 "dsyr.f"
    } else if (*incx == 0) {
#line 178 "dsyr.f"
	info = 5;
#line 179 "dsyr.f"
    } else if (*lda < max(1,*n)) {
#line 180 "dsyr.f"
	info = 7;
#line 181 "dsyr.f"
    }
#line 182 "dsyr.f"
    if (info != 0) {
#line 183 "dsyr.f"
	xerbla_("DSYR  ", &info, (ftnlen)6);
#line 184 "dsyr.f"
	return 0;
#line 185 "dsyr.f"
    }

/*     Quick return if possible. */

#line 189 "dsyr.f"
    if (*n == 0 || *alpha == 0.) {
#line 189 "dsyr.f"
	return 0;
#line 189 "dsyr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 193 "dsyr.f"
    if (*incx <= 0) {
#line 194 "dsyr.f"
	kx = 1 - (*n - 1) * *incx;
#line 195 "dsyr.f"
    } else if (*incx != 1) {
#line 196 "dsyr.f"
	kx = 1;
#line 197 "dsyr.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 203 "dsyr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in upper triangle. */

#line 207 "dsyr.f"
	if (*incx == 1) {
#line 208 "dsyr.f"
	    i__1 = *n;
#line 208 "dsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 209 "dsyr.f"
		if (x[j] != 0.) {
#line 210 "dsyr.f"
		    temp = *alpha * x[j];
#line 211 "dsyr.f"
		    i__2 = j;
#line 211 "dsyr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 212 "dsyr.f"
			a[i__ + j * a_dim1] += x[i__] * temp;
#line 213 "dsyr.f"
/* L10: */
#line 213 "dsyr.f"
		    }
#line 214 "dsyr.f"
		}
#line 215 "dsyr.f"
/* L20: */
#line 215 "dsyr.f"
	    }
#line 216 "dsyr.f"
	} else {
#line 217 "dsyr.f"
	    jx = kx;
#line 218 "dsyr.f"
	    i__1 = *n;
#line 218 "dsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 219 "dsyr.f"
		if (x[jx] != 0.) {
#line 220 "dsyr.f"
		    temp = *alpha * x[jx];
#line 221 "dsyr.f"
		    ix = kx;
#line 222 "dsyr.f"
		    i__2 = j;
#line 222 "dsyr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 223 "dsyr.f"
			a[i__ + j * a_dim1] += x[ix] * temp;
#line 224 "dsyr.f"
			ix += *incx;
#line 225 "dsyr.f"
/* L30: */
#line 225 "dsyr.f"
		    }
#line 226 "dsyr.f"
		}
#line 227 "dsyr.f"
		jx += *incx;
#line 228 "dsyr.f"
/* L40: */
#line 228 "dsyr.f"
	    }
#line 229 "dsyr.f"
	}
#line 230 "dsyr.f"
    } else {

/*        Form  A  when A is stored in lower triangle. */

#line 234 "dsyr.f"
	if (*incx == 1) {
#line 235 "dsyr.f"
	    i__1 = *n;
#line 235 "dsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 236 "dsyr.f"
		if (x[j] != 0.) {
#line 237 "dsyr.f"
		    temp = *alpha * x[j];
#line 238 "dsyr.f"
		    i__2 = *n;
#line 238 "dsyr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 239 "dsyr.f"
			a[i__ + j * a_dim1] += x[i__] * temp;
#line 240 "dsyr.f"
/* L50: */
#line 240 "dsyr.f"
		    }
#line 241 "dsyr.f"
		}
#line 242 "dsyr.f"
/* L60: */
#line 242 "dsyr.f"
	    }
#line 243 "dsyr.f"
	} else {
#line 244 "dsyr.f"
	    jx = kx;
#line 245 "dsyr.f"
	    i__1 = *n;
#line 245 "dsyr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 246 "dsyr.f"
		if (x[jx] != 0.) {
#line 247 "dsyr.f"
		    temp = *alpha * x[jx];
#line 248 "dsyr.f"
		    ix = jx;
#line 249 "dsyr.f"
		    i__2 = *n;
#line 249 "dsyr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 250 "dsyr.f"
			a[i__ + j * a_dim1] += x[ix] * temp;
#line 251 "dsyr.f"
			ix += *incx;
#line 252 "dsyr.f"
/* L70: */
#line 252 "dsyr.f"
		    }
#line 253 "dsyr.f"
		}
#line 254 "dsyr.f"
		jx += *incx;
#line 255 "dsyr.f"
/* L80: */
#line 255 "dsyr.f"
	    }
#line 256 "dsyr.f"
	}
#line 257 "dsyr.f"
    }

#line 259 "dsyr.f"
    return 0;

/*     End of DSYR  . */

} /* dsyr_ */

