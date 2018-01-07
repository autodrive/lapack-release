#line 1 "dtrmv.f"
/* dtrmv.f -- translated by f2c (version 20100827).
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

#line 1 "dtrmv.f"
/* > \brief \b DTRMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRMV  performs one of the matrix-vector operations */
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
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
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
/* >          X is DOUBLE PRECISION array of dimension at least */
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

/* > \ingroup double_blas_level2 */

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
/* Subroutine */ int dtrmv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 187 "dtrmv.f"
    /* Parameter adjustments */
#line 187 "dtrmv.f"
    a_dim1 = *lda;
#line 187 "dtrmv.f"
    a_offset = 1 + a_dim1;
#line 187 "dtrmv.f"
    a -= a_offset;
#line 187 "dtrmv.f"
    --x;
#line 187 "dtrmv.f"

#line 187 "dtrmv.f"
    /* Function Body */
#line 187 "dtrmv.f"
    info = 0;
#line 188 "dtrmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 189 "dtrmv.f"
	info = 1;
#line 190 "dtrmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 192 "dtrmv.f"
	info = 2;
#line 193 "dtrmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 194 "dtrmv.f"
	info = 3;
#line 195 "dtrmv.f"
    } else if (*n < 0) {
#line 196 "dtrmv.f"
	info = 4;
#line 197 "dtrmv.f"
    } else if (*lda < max(1,*n)) {
#line 198 "dtrmv.f"
	info = 6;
#line 199 "dtrmv.f"
    } else if (*incx == 0) {
#line 200 "dtrmv.f"
	info = 8;
#line 201 "dtrmv.f"
    }
#line 202 "dtrmv.f"
    if (info != 0) {
#line 203 "dtrmv.f"
	xerbla_("DTRMV ", &info, (ftnlen)6);
#line 204 "dtrmv.f"
	return 0;
#line 205 "dtrmv.f"
    }

/*     Quick return if possible. */

#line 209 "dtrmv.f"
    if (*n == 0) {
#line 209 "dtrmv.f"
	return 0;
#line 209 "dtrmv.f"
    }

#line 211 "dtrmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 216 "dtrmv.f"
    if (*incx <= 0) {
#line 217 "dtrmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 218 "dtrmv.f"
    } else if (*incx != 1) {
#line 219 "dtrmv.f"
	kx = 1;
#line 220 "dtrmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 225 "dtrmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := A*x. */

#line 229 "dtrmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "dtrmv.f"
	    if (*incx == 1) {
#line 231 "dtrmv.f"
		i__1 = *n;
#line 231 "dtrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 232 "dtrmv.f"
		    if (x[j] != 0.) {
#line 233 "dtrmv.f"
			temp = x[j];
#line 234 "dtrmv.f"
			i__2 = j - 1;
#line 234 "dtrmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 235 "dtrmv.f"
			    x[i__] += temp * a[i__ + j * a_dim1];
#line 236 "dtrmv.f"
/* L10: */
#line 236 "dtrmv.f"
			}
#line 237 "dtrmv.f"
			if (nounit) {
#line 237 "dtrmv.f"
			    x[j] *= a[j + j * a_dim1];
#line 237 "dtrmv.f"
			}
#line 238 "dtrmv.f"
		    }
#line 239 "dtrmv.f"
/* L20: */
#line 239 "dtrmv.f"
		}
#line 240 "dtrmv.f"
	    } else {
#line 241 "dtrmv.f"
		jx = kx;
#line 242 "dtrmv.f"
		i__1 = *n;
#line 242 "dtrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 243 "dtrmv.f"
		    if (x[jx] != 0.) {
#line 244 "dtrmv.f"
			temp = x[jx];
#line 245 "dtrmv.f"
			ix = kx;
#line 246 "dtrmv.f"
			i__2 = j - 1;
#line 246 "dtrmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 247 "dtrmv.f"
			    x[ix] += temp * a[i__ + j * a_dim1];
#line 248 "dtrmv.f"
			    ix += *incx;
#line 249 "dtrmv.f"
/* L30: */
#line 249 "dtrmv.f"
			}
#line 250 "dtrmv.f"
			if (nounit) {
#line 250 "dtrmv.f"
			    x[jx] *= a[j + j * a_dim1];
#line 250 "dtrmv.f"
			}
#line 251 "dtrmv.f"
		    }
#line 252 "dtrmv.f"
		    jx += *incx;
#line 253 "dtrmv.f"
/* L40: */
#line 253 "dtrmv.f"
		}
#line 254 "dtrmv.f"
	    }
#line 255 "dtrmv.f"
	} else {
#line 256 "dtrmv.f"
	    if (*incx == 1) {
#line 257 "dtrmv.f"
		for (j = *n; j >= 1; --j) {
#line 258 "dtrmv.f"
		    if (x[j] != 0.) {
#line 259 "dtrmv.f"
			temp = x[j];
#line 260 "dtrmv.f"
			i__1 = j + 1;
#line 260 "dtrmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 261 "dtrmv.f"
			    x[i__] += temp * a[i__ + j * a_dim1];
#line 262 "dtrmv.f"
/* L50: */
#line 262 "dtrmv.f"
			}
#line 263 "dtrmv.f"
			if (nounit) {
#line 263 "dtrmv.f"
			    x[j] *= a[j + j * a_dim1];
#line 263 "dtrmv.f"
			}
#line 264 "dtrmv.f"
		    }
#line 265 "dtrmv.f"
/* L60: */
#line 265 "dtrmv.f"
		}
#line 266 "dtrmv.f"
	    } else {
#line 267 "dtrmv.f"
		kx += (*n - 1) * *incx;
#line 268 "dtrmv.f"
		jx = kx;
#line 269 "dtrmv.f"
		for (j = *n; j >= 1; --j) {
#line 270 "dtrmv.f"
		    if (x[jx] != 0.) {
#line 271 "dtrmv.f"
			temp = x[jx];
#line 272 "dtrmv.f"
			ix = kx;
#line 273 "dtrmv.f"
			i__1 = j + 1;
#line 273 "dtrmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 274 "dtrmv.f"
			    x[ix] += temp * a[i__ + j * a_dim1];
#line 275 "dtrmv.f"
			    ix -= *incx;
#line 276 "dtrmv.f"
/* L70: */
#line 276 "dtrmv.f"
			}
#line 277 "dtrmv.f"
			if (nounit) {
#line 277 "dtrmv.f"
			    x[jx] *= a[j + j * a_dim1];
#line 277 "dtrmv.f"
			}
#line 278 "dtrmv.f"
		    }
#line 279 "dtrmv.f"
		    jx -= *incx;
#line 280 "dtrmv.f"
/* L80: */
#line 280 "dtrmv.f"
		}
#line 281 "dtrmv.f"
	    }
#line 282 "dtrmv.f"
	}
#line 283 "dtrmv.f"
    } else {

/*        Form  x := A**T*x. */

#line 287 "dtrmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 288 "dtrmv.f"
	    if (*incx == 1) {
#line 289 "dtrmv.f"
		for (j = *n; j >= 1; --j) {
#line 290 "dtrmv.f"
		    temp = x[j];
#line 291 "dtrmv.f"
		    if (nounit) {
#line 291 "dtrmv.f"
			temp *= a[j + j * a_dim1];
#line 291 "dtrmv.f"
		    }
#line 292 "dtrmv.f"
		    for (i__ = j - 1; i__ >= 1; --i__) {
#line 293 "dtrmv.f"
			temp += a[i__ + j * a_dim1] * x[i__];
#line 294 "dtrmv.f"
/* L90: */
#line 294 "dtrmv.f"
		    }
#line 295 "dtrmv.f"
		    x[j] = temp;
#line 296 "dtrmv.f"
/* L100: */
#line 296 "dtrmv.f"
		}
#line 297 "dtrmv.f"
	    } else {
#line 298 "dtrmv.f"
		jx = kx + (*n - 1) * *incx;
#line 299 "dtrmv.f"
		for (j = *n; j >= 1; --j) {
#line 300 "dtrmv.f"
		    temp = x[jx];
#line 301 "dtrmv.f"
		    ix = jx;
#line 302 "dtrmv.f"
		    if (nounit) {
#line 302 "dtrmv.f"
			temp *= a[j + j * a_dim1];
#line 302 "dtrmv.f"
		    }
#line 303 "dtrmv.f"
		    for (i__ = j - 1; i__ >= 1; --i__) {
#line 304 "dtrmv.f"
			ix -= *incx;
#line 305 "dtrmv.f"
			temp += a[i__ + j * a_dim1] * x[ix];
#line 306 "dtrmv.f"
/* L110: */
#line 306 "dtrmv.f"
		    }
#line 307 "dtrmv.f"
		    x[jx] = temp;
#line 308 "dtrmv.f"
		    jx -= *incx;
#line 309 "dtrmv.f"
/* L120: */
#line 309 "dtrmv.f"
		}
#line 310 "dtrmv.f"
	    }
#line 311 "dtrmv.f"
	} else {
#line 312 "dtrmv.f"
	    if (*incx == 1) {
#line 313 "dtrmv.f"
		i__1 = *n;
#line 313 "dtrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 314 "dtrmv.f"
		    temp = x[j];
#line 315 "dtrmv.f"
		    if (nounit) {
#line 315 "dtrmv.f"
			temp *= a[j + j * a_dim1];
#line 315 "dtrmv.f"
		    }
#line 316 "dtrmv.f"
		    i__2 = *n;
#line 316 "dtrmv.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 317 "dtrmv.f"
			temp += a[i__ + j * a_dim1] * x[i__];
#line 318 "dtrmv.f"
/* L130: */
#line 318 "dtrmv.f"
		    }
#line 319 "dtrmv.f"
		    x[j] = temp;
#line 320 "dtrmv.f"
/* L140: */
#line 320 "dtrmv.f"
		}
#line 321 "dtrmv.f"
	    } else {
#line 322 "dtrmv.f"
		jx = kx;
#line 323 "dtrmv.f"
		i__1 = *n;
#line 323 "dtrmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 324 "dtrmv.f"
		    temp = x[jx];
#line 325 "dtrmv.f"
		    ix = jx;
#line 326 "dtrmv.f"
		    if (nounit) {
#line 326 "dtrmv.f"
			temp *= a[j + j * a_dim1];
#line 326 "dtrmv.f"
		    }
#line 327 "dtrmv.f"
		    i__2 = *n;
#line 327 "dtrmv.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 328 "dtrmv.f"
			ix += *incx;
#line 329 "dtrmv.f"
			temp += a[i__ + j * a_dim1] * x[ix];
#line 330 "dtrmv.f"
/* L150: */
#line 330 "dtrmv.f"
		    }
#line 331 "dtrmv.f"
		    x[jx] = temp;
#line 332 "dtrmv.f"
		    jx += *incx;
#line 333 "dtrmv.f"
/* L160: */
#line 333 "dtrmv.f"
		}
#line 334 "dtrmv.f"
	    }
#line 335 "dtrmv.f"
	}
#line 336 "dtrmv.f"
    }

#line 338 "dtrmv.f"
    return 0;

/*     End of DTRMV . */

} /* dtrmv_ */

