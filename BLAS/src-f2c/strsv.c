#line 1 "strsv.f"
/* strsv.f -- translated by f2c (version 20100827).
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

#line 1 "strsv.f"
/* > \brief \b STRSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

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
/* > STRSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular matrix. */
/* > */
/* > No test for singularity or near-singularity is included in this */
/* > routine. Such tests must be performed before calling this routine. */
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
/* >           On entry, TRANS specifies the equations to be solved as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   A*x = b. */
/* > */
/* >              TRANS = 'T' or 't'   A**T*x = b. */
/* > */
/* >              TRANS = 'C' or 'c'   A**T*x = b. */
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
/* >          A is REAL array of DIMENSION ( LDA, n ). */
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
/* >          X is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element right-hand side vector b. On exit, X is overwritten */
/* >           with the solution vector x. */
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
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int strsv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 189 "strsv.f"
    /* Parameter adjustments */
#line 189 "strsv.f"
    a_dim1 = *lda;
#line 189 "strsv.f"
    a_offset = 1 + a_dim1;
#line 189 "strsv.f"
    a -= a_offset;
#line 189 "strsv.f"
    --x;
#line 189 "strsv.f"

#line 189 "strsv.f"
    /* Function Body */
#line 189 "strsv.f"
    info = 0;
#line 190 "strsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 191 "strsv.f"
	info = 1;
#line 192 "strsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 194 "strsv.f"
	info = 2;
#line 195 "strsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 196 "strsv.f"
	info = 3;
#line 197 "strsv.f"
    } else if (*n < 0) {
#line 198 "strsv.f"
	info = 4;
#line 199 "strsv.f"
    } else if (*lda < max(1,*n)) {
#line 200 "strsv.f"
	info = 6;
#line 201 "strsv.f"
    } else if (*incx == 0) {
#line 202 "strsv.f"
	info = 8;
#line 203 "strsv.f"
    }
#line 204 "strsv.f"
    if (info != 0) {
#line 205 "strsv.f"
	xerbla_("STRSV ", &info, (ftnlen)6);
#line 206 "strsv.f"
	return 0;
#line 207 "strsv.f"
    }

/*     Quick return if possible. */

#line 211 "strsv.f"
    if (*n == 0) {
#line 211 "strsv.f"
	return 0;
#line 211 "strsv.f"
    }

#line 213 "strsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 218 "strsv.f"
    if (*incx <= 0) {
#line 219 "strsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 220 "strsv.f"
    } else if (*incx != 1) {
#line 221 "strsv.f"
	kx = 1;
#line 222 "strsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 227 "strsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 231 "strsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 232 "strsv.f"
	    if (*incx == 1) {
#line 233 "strsv.f"
		for (j = *n; j >= 1; --j) {
#line 234 "strsv.f"
		    if (x[j] != 0.) {
#line 235 "strsv.f"
			if (nounit) {
#line 235 "strsv.f"
			    x[j] /= a[j + j * a_dim1];
#line 235 "strsv.f"
			}
#line 236 "strsv.f"
			temp = x[j];
#line 237 "strsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 238 "strsv.f"
			    x[i__] -= temp * a[i__ + j * a_dim1];
#line 239 "strsv.f"
/* L10: */
#line 239 "strsv.f"
			}
#line 240 "strsv.f"
		    }
#line 241 "strsv.f"
/* L20: */
#line 241 "strsv.f"
		}
#line 242 "strsv.f"
	    } else {
#line 243 "strsv.f"
		jx = kx + (*n - 1) * *incx;
#line 244 "strsv.f"
		for (j = *n; j >= 1; --j) {
#line 245 "strsv.f"
		    if (x[jx] != 0.) {
#line 246 "strsv.f"
			if (nounit) {
#line 246 "strsv.f"
			    x[jx] /= a[j + j * a_dim1];
#line 246 "strsv.f"
			}
#line 247 "strsv.f"
			temp = x[jx];
#line 248 "strsv.f"
			ix = jx;
#line 249 "strsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 250 "strsv.f"
			    ix -= *incx;
#line 251 "strsv.f"
			    x[ix] -= temp * a[i__ + j * a_dim1];
#line 252 "strsv.f"
/* L30: */
#line 252 "strsv.f"
			}
#line 253 "strsv.f"
		    }
#line 254 "strsv.f"
		    jx -= *incx;
#line 255 "strsv.f"
/* L40: */
#line 255 "strsv.f"
		}
#line 256 "strsv.f"
	    }
#line 257 "strsv.f"
	} else {
#line 258 "strsv.f"
	    if (*incx == 1) {
#line 259 "strsv.f"
		i__1 = *n;
#line 259 "strsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "strsv.f"
		    if (x[j] != 0.) {
#line 261 "strsv.f"
			if (nounit) {
#line 261 "strsv.f"
			    x[j] /= a[j + j * a_dim1];
#line 261 "strsv.f"
			}
#line 262 "strsv.f"
			temp = x[j];
#line 263 "strsv.f"
			i__2 = *n;
#line 263 "strsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 264 "strsv.f"
			    x[i__] -= temp * a[i__ + j * a_dim1];
#line 265 "strsv.f"
/* L50: */
#line 265 "strsv.f"
			}
#line 266 "strsv.f"
		    }
#line 267 "strsv.f"
/* L60: */
#line 267 "strsv.f"
		}
#line 268 "strsv.f"
	    } else {
#line 269 "strsv.f"
		jx = kx;
#line 270 "strsv.f"
		i__1 = *n;
#line 270 "strsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 271 "strsv.f"
		    if (x[jx] != 0.) {
#line 272 "strsv.f"
			if (nounit) {
#line 272 "strsv.f"
			    x[jx] /= a[j + j * a_dim1];
#line 272 "strsv.f"
			}
#line 273 "strsv.f"
			temp = x[jx];
#line 274 "strsv.f"
			ix = jx;
#line 275 "strsv.f"
			i__2 = *n;
#line 275 "strsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 276 "strsv.f"
			    ix += *incx;
#line 277 "strsv.f"
			    x[ix] -= temp * a[i__ + j * a_dim1];
#line 278 "strsv.f"
/* L70: */
#line 278 "strsv.f"
			}
#line 279 "strsv.f"
		    }
#line 280 "strsv.f"
		    jx += *incx;
#line 281 "strsv.f"
/* L80: */
#line 281 "strsv.f"
		}
#line 282 "strsv.f"
	    }
#line 283 "strsv.f"
	}
#line 284 "strsv.f"
    } else {

/*        Form  x := inv( A**T )*x. */

#line 288 "strsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 289 "strsv.f"
	    if (*incx == 1) {
#line 290 "strsv.f"
		i__1 = *n;
#line 290 "strsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 291 "strsv.f"
		    temp = x[j];
#line 292 "strsv.f"
		    i__2 = j - 1;
#line 292 "strsv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 293 "strsv.f"
			temp -= a[i__ + j * a_dim1] * x[i__];
#line 294 "strsv.f"
/* L90: */
#line 294 "strsv.f"
		    }
#line 295 "strsv.f"
		    if (nounit) {
#line 295 "strsv.f"
			temp /= a[j + j * a_dim1];
#line 295 "strsv.f"
		    }
#line 296 "strsv.f"
		    x[j] = temp;
#line 297 "strsv.f"
/* L100: */
#line 297 "strsv.f"
		}
#line 298 "strsv.f"
	    } else {
#line 299 "strsv.f"
		jx = kx;
#line 300 "strsv.f"
		i__1 = *n;
#line 300 "strsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 301 "strsv.f"
		    temp = x[jx];
#line 302 "strsv.f"
		    ix = kx;
#line 303 "strsv.f"
		    i__2 = j - 1;
#line 303 "strsv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "strsv.f"
			temp -= a[i__ + j * a_dim1] * x[ix];
#line 305 "strsv.f"
			ix += *incx;
#line 306 "strsv.f"
/* L110: */
#line 306 "strsv.f"
		    }
#line 307 "strsv.f"
		    if (nounit) {
#line 307 "strsv.f"
			temp /= a[j + j * a_dim1];
#line 307 "strsv.f"
		    }
#line 308 "strsv.f"
		    x[jx] = temp;
#line 309 "strsv.f"
		    jx += *incx;
#line 310 "strsv.f"
/* L120: */
#line 310 "strsv.f"
		}
#line 311 "strsv.f"
	    }
#line 312 "strsv.f"
	} else {
#line 313 "strsv.f"
	    if (*incx == 1) {
#line 314 "strsv.f"
		for (j = *n; j >= 1; --j) {
#line 315 "strsv.f"
		    temp = x[j];
#line 316 "strsv.f"
		    i__1 = j + 1;
#line 316 "strsv.f"
		    for (i__ = *n; i__ >= i__1; --i__) {
#line 317 "strsv.f"
			temp -= a[i__ + j * a_dim1] * x[i__];
#line 318 "strsv.f"
/* L130: */
#line 318 "strsv.f"
		    }
#line 319 "strsv.f"
		    if (nounit) {
#line 319 "strsv.f"
			temp /= a[j + j * a_dim1];
#line 319 "strsv.f"
		    }
#line 320 "strsv.f"
		    x[j] = temp;
#line 321 "strsv.f"
/* L140: */
#line 321 "strsv.f"
		}
#line 322 "strsv.f"
	    } else {
#line 323 "strsv.f"
		kx += (*n - 1) * *incx;
#line 324 "strsv.f"
		jx = kx;
#line 325 "strsv.f"
		for (j = *n; j >= 1; --j) {
#line 326 "strsv.f"
		    temp = x[jx];
#line 327 "strsv.f"
		    ix = kx;
#line 328 "strsv.f"
		    i__1 = j + 1;
#line 328 "strsv.f"
		    for (i__ = *n; i__ >= i__1; --i__) {
#line 329 "strsv.f"
			temp -= a[i__ + j * a_dim1] * x[ix];
#line 330 "strsv.f"
			ix -= *incx;
#line 331 "strsv.f"
/* L150: */
#line 331 "strsv.f"
		    }
#line 332 "strsv.f"
		    if (nounit) {
#line 332 "strsv.f"
			temp /= a[j + j * a_dim1];
#line 332 "strsv.f"
		    }
#line 333 "strsv.f"
		    x[jx] = temp;
#line 334 "strsv.f"
		    jx -= *incx;
#line 335 "strsv.f"
/* L160: */
#line 335 "strsv.f"
		}
#line 336 "strsv.f"
	    }
#line 337 "strsv.f"
	}
#line 338 "strsv.f"
    }

#line 340 "strsv.f"
    return 0;

/*     End of STRSV . */

} /* strsv_ */

