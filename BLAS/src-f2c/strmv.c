#line 1 "strmv.f"
/* strmv.f -- translated by f2c (version 20100827).
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

#line 1 "strmv.f"
/* > \brief \b STRMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   x := A*x. */
/* > */
/* >              TRANS = 'T' or 't'   x := A**T*x. */
/* > */
/* >              TRANS = 'C' or 'c'   x := A**T*x. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit */
/* >           triangular as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension ( LDA, N ) */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular matrix and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular matrix and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced either, but are assumed to be unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. On exit, X is overwritten with the */
/* >           transformed vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup single_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int strmv_(char *uplo, char *trans, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, ftnlen 
	uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, jx, kx, info;
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

#line 187 "strmv.f"
    /* Parameter adjustments */
#line 187 "strmv.f"
    a_dim1 = *lda;
#line 187 "strmv.f"
    a_offset = 1 + a_dim1;
#line 187 "strmv.f"
    a -= a_offset;
#line 187 "strmv.f"
    --x;
#line 187 "strmv.f"

#line 187 "strmv.f"
    /* Function Body */
#line 187 "strmv.f"
    info = 0;
#line 188 "strmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 189 "strmv.f"
	info = 1;
#line 190 "strmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 192 "strmv.f"
	info = 2;
#line 193 "strmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 194 "strmv.f"
	info = 3;
#line 195 "strmv.f"
    } else if (*n < 0) {
#line 196 "strmv.f"
	info = 4;
#line 197 "strmv.f"
    } else if (*lda < max(1,*n)) {
#line 198 "strmv.f"
	info = 6;
#line 199 "strmv.f"
    } else if (*incx == 0) {
#line 200 "strmv.f"
	info = 8;
#line 201 "strmv.f"
    }
#line 202 "strmv.f"
    if (info != 0) {
#line 203 "strmv.f"
	xerbla_("STRMV ", &info, (ftnlen)6);
#line 204 "strmv.f"
	return 0;
#line 205 "strmv.f"
    }

/*     Quick return if possible. */

#line 209 "strmv.f"
    if (*n == 0) {
#line 209 "strmv.f"
	return 0;
#line 209 "strmv.f"
    }

#line 211 "strmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 216 "strmv.f"
    if (*incx <= 0) {
#line 217 "strmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 218 "strmv.f"
    } else if (*incx != 1) {
#line 219 "strmv.f"
	kx = 1;
#line 220 "strmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 225 "strmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := A*x. */

#line 229 "strmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "strmv.f"
	    if (*incx == 1) {
#line 231 "strmv.f"
		i__1 = *n;
#line 231 "strmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 232 "strmv.f"
		    if (x[j] != 0.) {
#line 233 "strmv.f"
			temp = x[j];
#line 234 "strmv.f"
			i__2 = j - 1;
#line 234 "strmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 235 "strmv.f"
			    x[i__] += temp * a[i__ + j * a_dim1];
#line 236 "strmv.f"
/* L10: */
#line 236 "strmv.f"
			}
#line 237 "strmv.f"
			if (nounit) {
#line 237 "strmv.f"
			    x[j] *= a[j + j * a_dim1];
#line 237 "strmv.f"
			}
#line 238 "strmv.f"
		    }
#line 239 "strmv.f"
/* L20: */
#line 239 "strmv.f"
		}
#line 240 "strmv.f"
	    } else {
#line 241 "strmv.f"
		jx = kx;
#line 242 "strmv.f"
		i__1 = *n;
#line 242 "strmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 243 "strmv.f"
		    if (x[jx] != 0.) {
#line 244 "strmv.f"
			temp = x[jx];
#line 245 "strmv.f"
			ix = kx;
#line 246 "strmv.f"
			i__2 = j - 1;
#line 246 "strmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 247 "strmv.f"
			    x[ix] += temp * a[i__ + j * a_dim1];
#line 248 "strmv.f"
			    ix += *incx;
#line 249 "strmv.f"
/* L30: */
#line 249 "strmv.f"
			}
#line 250 "strmv.f"
			if (nounit) {
#line 250 "strmv.f"
			    x[jx] *= a[j + j * a_dim1];
#line 250 "strmv.f"
			}
#line 251 "strmv.f"
		    }
#line 252 "strmv.f"
		    jx += *incx;
#line 253 "strmv.f"
/* L40: */
#line 253 "strmv.f"
		}
#line 254 "strmv.f"
	    }
#line 255 "strmv.f"
	} else {
#line 256 "strmv.f"
	    if (*incx == 1) {
#line 257 "strmv.f"
		for (j = *n; j >= 1; --j) {
#line 258 "strmv.f"
		    if (x[j] != 0.) {
#line 259 "strmv.f"
			temp = x[j];
#line 260 "strmv.f"
			i__1 = j + 1;
#line 260 "strmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 261 "strmv.f"
			    x[i__] += temp * a[i__ + j * a_dim1];
#line 262 "strmv.f"
/* L50: */
#line 262 "strmv.f"
			}
#line 263 "strmv.f"
			if (nounit) {
#line 263 "strmv.f"
			    x[j] *= a[j + j * a_dim1];
#line 263 "strmv.f"
			}
#line 264 "strmv.f"
		    }
#line 265 "strmv.f"
/* L60: */
#line 265 "strmv.f"
		}
#line 266 "strmv.f"
	    } else {
#line 267 "strmv.f"
		kx += (*n - 1) * *incx;
#line 268 "strmv.f"
		jx = kx;
#line 269 "strmv.f"
		for (j = *n; j >= 1; --j) {
#line 270 "strmv.f"
		    if (x[jx] != 0.) {
#line 271 "strmv.f"
			temp = x[jx];
#line 272 "strmv.f"
			ix = kx;
#line 273 "strmv.f"
			i__1 = j + 1;
#line 273 "strmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 274 "strmv.f"
			    x[ix] += temp * a[i__ + j * a_dim1];
#line 275 "strmv.f"
			    ix -= *incx;
#line 276 "strmv.f"
/* L70: */
#line 276 "strmv.f"
			}
#line 277 "strmv.f"
			if (nounit) {
#line 277 "strmv.f"
			    x[jx] *= a[j + j * a_dim1];
#line 277 "strmv.f"
			}
#line 278 "strmv.f"
		    }
#line 279 "strmv.f"
		    jx -= *incx;
#line 280 "strmv.f"
/* L80: */
#line 280 "strmv.f"
		}
#line 281 "strmv.f"
	    }
#line 282 "strmv.f"
	}
#line 283 "strmv.f"
    } else {

/*        Form  x := A**T*x. */

#line 287 "strmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 288 "strmv.f"
	    if (*incx == 1) {
#line 289 "strmv.f"
		for (j = *n; j >= 1; --j) {
#line 290 "strmv.f"
		    temp = x[j];
#line 291 "strmv.f"
		    if (nounit) {
#line 291 "strmv.f"
			temp *= a[j + j * a_dim1];
#line 291 "strmv.f"
		    }
#line 292 "strmv.f"
		    for (i__ = j - 1; i__ >= 1; --i__) {
#line 293 "strmv.f"
			temp += a[i__ + j * a_dim1] * x[i__];
#line 294 "strmv.f"
/* L90: */
#line 294 "strmv.f"
		    }
#line 295 "strmv.f"
		    x[j] = temp;
#line 296 "strmv.f"
/* L100: */
#line 296 "strmv.f"
		}
#line 297 "strmv.f"
	    } else {
#line 298 "strmv.f"
		jx = kx + (*n - 1) * *incx;
#line 299 "strmv.f"
		for (j = *n; j >= 1; --j) {
#line 300 "strmv.f"
		    temp = x[jx];
#line 301 "strmv.f"
		    ix = jx;
#line 302 "strmv.f"
		    if (nounit) {
#line 302 "strmv.f"
			temp *= a[j + j * a_dim1];
#line 302 "strmv.f"
		    }
#line 303 "strmv.f"
		    for (i__ = j - 1; i__ >= 1; --i__) {
#line 304 "strmv.f"
			ix -= *incx;
#line 305 "strmv.f"
			temp += a[i__ + j * a_dim1] * x[ix];
#line 306 "strmv.f"
/* L110: */
#line 306 "strmv.f"
		    }
#line 307 "strmv.f"
		    x[jx] = temp;
#line 308 "strmv.f"
		    jx -= *incx;
#line 309 "strmv.f"
/* L120: */
#line 309 "strmv.f"
		}
#line 310 "strmv.f"
	    }
#line 311 "strmv.f"
	} else {
#line 312 "strmv.f"
	    if (*incx == 1) {
#line 313 "strmv.f"
		i__1 = *n;
#line 313 "strmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 314 "strmv.f"
		    temp = x[j];
#line 315 "strmv.f"
		    if (nounit) {
#line 315 "strmv.f"
			temp *= a[j + j * a_dim1];
#line 315 "strmv.f"
		    }
#line 316 "strmv.f"
		    i__2 = *n;
#line 316 "strmv.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 317 "strmv.f"
			temp += a[i__ + j * a_dim1] * x[i__];
#line 318 "strmv.f"
/* L130: */
#line 318 "strmv.f"
		    }
#line 319 "strmv.f"
		    x[j] = temp;
#line 320 "strmv.f"
/* L140: */
#line 320 "strmv.f"
		}
#line 321 "strmv.f"
	    } else {
#line 322 "strmv.f"
		jx = kx;
#line 323 "strmv.f"
		i__1 = *n;
#line 323 "strmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 324 "strmv.f"
		    temp = x[jx];
#line 325 "strmv.f"
		    ix = jx;
#line 326 "strmv.f"
		    if (nounit) {
#line 326 "strmv.f"
			temp *= a[j + j * a_dim1];
#line 326 "strmv.f"
		    }
#line 327 "strmv.f"
		    i__2 = *n;
#line 327 "strmv.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 328 "strmv.f"
			ix += *incx;
#line 329 "strmv.f"
			temp += a[i__ + j * a_dim1] * x[ix];
#line 330 "strmv.f"
/* L150: */
#line 330 "strmv.f"
		    }
#line 331 "strmv.f"
		    x[jx] = temp;
#line 332 "strmv.f"
		    jx += *incx;
#line 333 "strmv.f"
/* L160: */
#line 333 "strmv.f"
		}
#line 334 "strmv.f"
	    }
#line 335 "strmv.f"
	}
#line 336 "strmv.f"
    }

#line 338 "strmv.f"
    return 0;

/*     End of STRMV . */

} /* strmv_ */

