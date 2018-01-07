#line 1 "dspr.f"
/* dspr.f -- translated by f2c (version 20100827).
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

#line 1 "dspr.f"
/* > \brief \b DSPR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER INCX,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPR    performs the symmetric rank 1 operation */
/* > */
/* >    A := alpha*x*x**T + A, */
/* > */
/* > where alpha is a real scalar, x is an n element vector and A is an */
/* > n by n symmetric matrix, supplied in packed form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the matrix A is supplied in the packed */
/* >           array AP as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   The upper triangular part of A is */
/* >                                  supplied in AP. */
/* > */
/* >              UPLO = 'L' or 'l'   The lower triangular part of A is */
/* >                                  supplied in AP. */
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
/* >          X is DOUBLE PRECISION array, dimension at least */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the upper triangular part of the */
/* >           updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the lower triangular part of the */
/* >           updated matrix. */
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
/* Subroutine */ int dspr_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static doublereal temp;
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

/*     Test the input parameters. */

#line 164 "dspr.f"
    /* Parameter adjustments */
#line 164 "dspr.f"
    --ap;
#line 164 "dspr.f"
    --x;
#line 164 "dspr.f"

#line 164 "dspr.f"
    /* Function Body */
#line 164 "dspr.f"
    info = 0;
#line 165 "dspr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 166 "dspr.f"
	info = 1;
#line 167 "dspr.f"
    } else if (*n < 0) {
#line 168 "dspr.f"
	info = 2;
#line 169 "dspr.f"
    } else if (*incx == 0) {
#line 170 "dspr.f"
	info = 5;
#line 171 "dspr.f"
    }
#line 172 "dspr.f"
    if (info != 0) {
#line 173 "dspr.f"
	xerbla_("DSPR  ", &info, (ftnlen)6);
#line 174 "dspr.f"
	return 0;
#line 175 "dspr.f"
    }

/*     Quick return if possible. */

#line 179 "dspr.f"
    if (*n == 0 || *alpha == 0.) {
#line 179 "dspr.f"
	return 0;
#line 179 "dspr.f"
    }

/*     Set the start point in X if the increment is not unity. */

#line 183 "dspr.f"
    if (*incx <= 0) {
#line 184 "dspr.f"
	kx = 1 - (*n - 1) * *incx;
#line 185 "dspr.f"
    } else if (*incx != 1) {
#line 186 "dspr.f"
	kx = 1;
#line 187 "dspr.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 192 "dspr.f"
    kk = 1;
#line 193 "dspr.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 197 "dspr.f"
	if (*incx == 1) {
#line 198 "dspr.f"
	    i__1 = *n;
#line 198 "dspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 199 "dspr.f"
		if (x[j] != 0.) {
#line 200 "dspr.f"
		    temp = *alpha * x[j];
#line 201 "dspr.f"
		    k = kk;
#line 202 "dspr.f"
		    i__2 = j;
#line 202 "dspr.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 203 "dspr.f"
			ap[k] += x[i__] * temp;
#line 204 "dspr.f"
			++k;
#line 205 "dspr.f"
/* L10: */
#line 205 "dspr.f"
		    }
#line 206 "dspr.f"
		}
#line 207 "dspr.f"
		kk += j;
#line 208 "dspr.f"
/* L20: */
#line 208 "dspr.f"
	    }
#line 209 "dspr.f"
	} else {
#line 210 "dspr.f"
	    jx = kx;
#line 211 "dspr.f"
	    i__1 = *n;
#line 211 "dspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 212 "dspr.f"
		if (x[jx] != 0.) {
#line 213 "dspr.f"
		    temp = *alpha * x[jx];
#line 214 "dspr.f"
		    ix = kx;
#line 215 "dspr.f"
		    i__2 = kk + j - 1;
#line 215 "dspr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 216 "dspr.f"
			ap[k] += x[ix] * temp;
#line 217 "dspr.f"
			ix += *incx;
#line 218 "dspr.f"
/* L30: */
#line 218 "dspr.f"
		    }
#line 219 "dspr.f"
		}
#line 220 "dspr.f"
		jx += *incx;
#line 221 "dspr.f"
		kk += j;
#line 222 "dspr.f"
/* L40: */
#line 222 "dspr.f"
	    }
#line 223 "dspr.f"
	}
#line 224 "dspr.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 228 "dspr.f"
	if (*incx == 1) {
#line 229 "dspr.f"
	    i__1 = *n;
#line 229 "dspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 230 "dspr.f"
		if (x[j] != 0.) {
#line 231 "dspr.f"
		    temp = *alpha * x[j];
#line 232 "dspr.f"
		    k = kk;
#line 233 "dspr.f"
		    i__2 = *n;
#line 233 "dspr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 234 "dspr.f"
			ap[k] += x[i__] * temp;
#line 235 "dspr.f"
			++k;
#line 236 "dspr.f"
/* L50: */
#line 236 "dspr.f"
		    }
#line 237 "dspr.f"
		}
#line 238 "dspr.f"
		kk = kk + *n - j + 1;
#line 239 "dspr.f"
/* L60: */
#line 239 "dspr.f"
	    }
#line 240 "dspr.f"
	} else {
#line 241 "dspr.f"
	    jx = kx;
#line 242 "dspr.f"
	    i__1 = *n;
#line 242 "dspr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 243 "dspr.f"
		if (x[jx] != 0.) {
#line 244 "dspr.f"
		    temp = *alpha * x[jx];
#line 245 "dspr.f"
		    ix = jx;
#line 246 "dspr.f"
		    i__2 = kk + *n - j;
#line 246 "dspr.f"
		    for (k = kk; k <= i__2; ++k) {
#line 247 "dspr.f"
			ap[k] += x[ix] * temp;
#line 248 "dspr.f"
			ix += *incx;
#line 249 "dspr.f"
/* L70: */
#line 249 "dspr.f"
		    }
#line 250 "dspr.f"
		}
#line 251 "dspr.f"
		jx += *incx;
#line 252 "dspr.f"
		kk = kk + *n - j + 1;
#line 253 "dspr.f"
/* L80: */
#line 253 "dspr.f"
	    }
#line 254 "dspr.f"
	}
#line 255 "dspr.f"
    }

#line 257 "dspr.f"
    return 0;

/*     End of DSPR  . */

} /* dspr_ */

