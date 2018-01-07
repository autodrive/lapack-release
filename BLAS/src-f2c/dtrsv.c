#line 1 "dtrsv.f"
/* dtrsv.f -- translated by f2c (version 20100827).
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

#line 1 "dtrsv.f"
/* > \brief \b DTRSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

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
/* > DTRSV  solves one of the systems of equations */
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
/* >           element right-hand side vector b. On exit, X is overwritten */
/* >           with the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > */
/* >  Level 2 Blas routine. */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int dtrsv_(char *uplo, char *trans, char *diag, integer *n, 
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


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
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

#line 183 "dtrsv.f"
    /* Parameter adjustments */
#line 183 "dtrsv.f"
    a_dim1 = *lda;
#line 183 "dtrsv.f"
    a_offset = 1 + a_dim1;
#line 183 "dtrsv.f"
    a -= a_offset;
#line 183 "dtrsv.f"
    --x;
#line 183 "dtrsv.f"

#line 183 "dtrsv.f"
    /* Function Body */
#line 183 "dtrsv.f"
    info = 0;
#line 184 "dtrsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 185 "dtrsv.f"
	info = 1;
#line 186 "dtrsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 188 "dtrsv.f"
	info = 2;
#line 189 "dtrsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 190 "dtrsv.f"
	info = 3;
#line 191 "dtrsv.f"
    } else if (*n < 0) {
#line 192 "dtrsv.f"
	info = 4;
#line 193 "dtrsv.f"
    } else if (*lda < max(1,*n)) {
#line 194 "dtrsv.f"
	info = 6;
#line 195 "dtrsv.f"
    } else if (*incx == 0) {
#line 196 "dtrsv.f"
	info = 8;
#line 197 "dtrsv.f"
    }
#line 198 "dtrsv.f"
    if (info != 0) {
#line 199 "dtrsv.f"
	xerbla_("DTRSV ", &info, (ftnlen)6);
#line 200 "dtrsv.f"
	return 0;
#line 201 "dtrsv.f"
    }

/*     Quick return if possible. */

#line 205 "dtrsv.f"
    if (*n == 0) {
#line 205 "dtrsv.f"
	return 0;
#line 205 "dtrsv.f"
    }

#line 207 "dtrsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 212 "dtrsv.f"
    if (*incx <= 0) {
#line 213 "dtrsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 214 "dtrsv.f"
    } else if (*incx != 1) {
#line 215 "dtrsv.f"
	kx = 1;
#line 216 "dtrsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 221 "dtrsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 225 "dtrsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 226 "dtrsv.f"
	    if (*incx == 1) {
#line 227 "dtrsv.f"
		for (j = *n; j >= 1; --j) {
#line 228 "dtrsv.f"
		    if (x[j] != 0.) {
#line 229 "dtrsv.f"
			if (nounit) {
#line 229 "dtrsv.f"
			    x[j] /= a[j + j * a_dim1];
#line 229 "dtrsv.f"
			}
#line 230 "dtrsv.f"
			temp = x[j];
#line 231 "dtrsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 232 "dtrsv.f"
			    x[i__] -= temp * a[i__ + j * a_dim1];
#line 233 "dtrsv.f"
/* L10: */
#line 233 "dtrsv.f"
			}
#line 234 "dtrsv.f"
		    }
#line 235 "dtrsv.f"
/* L20: */
#line 235 "dtrsv.f"
		}
#line 236 "dtrsv.f"
	    } else {
#line 237 "dtrsv.f"
		jx = kx + (*n - 1) * *incx;
#line 238 "dtrsv.f"
		for (j = *n; j >= 1; --j) {
#line 239 "dtrsv.f"
		    if (x[jx] != 0.) {
#line 240 "dtrsv.f"
			if (nounit) {
#line 240 "dtrsv.f"
			    x[jx] /= a[j + j * a_dim1];
#line 240 "dtrsv.f"
			}
#line 241 "dtrsv.f"
			temp = x[jx];
#line 242 "dtrsv.f"
			ix = jx;
#line 243 "dtrsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 244 "dtrsv.f"
			    ix -= *incx;
#line 245 "dtrsv.f"
			    x[ix] -= temp * a[i__ + j * a_dim1];
#line 246 "dtrsv.f"
/* L30: */
#line 246 "dtrsv.f"
			}
#line 247 "dtrsv.f"
		    }
#line 248 "dtrsv.f"
		    jx -= *incx;
#line 249 "dtrsv.f"
/* L40: */
#line 249 "dtrsv.f"
		}
#line 250 "dtrsv.f"
	    }
#line 251 "dtrsv.f"
	} else {
#line 252 "dtrsv.f"
	    if (*incx == 1) {
#line 253 "dtrsv.f"
		i__1 = *n;
#line 253 "dtrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 254 "dtrsv.f"
		    if (x[j] != 0.) {
#line 255 "dtrsv.f"
			if (nounit) {
#line 255 "dtrsv.f"
			    x[j] /= a[j + j * a_dim1];
#line 255 "dtrsv.f"
			}
#line 256 "dtrsv.f"
			temp = x[j];
#line 257 "dtrsv.f"
			i__2 = *n;
#line 257 "dtrsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 258 "dtrsv.f"
			    x[i__] -= temp * a[i__ + j * a_dim1];
#line 259 "dtrsv.f"
/* L50: */
#line 259 "dtrsv.f"
			}
#line 260 "dtrsv.f"
		    }
#line 261 "dtrsv.f"
/* L60: */
#line 261 "dtrsv.f"
		}
#line 262 "dtrsv.f"
	    } else {
#line 263 "dtrsv.f"
		jx = kx;
#line 264 "dtrsv.f"
		i__1 = *n;
#line 264 "dtrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 265 "dtrsv.f"
		    if (x[jx] != 0.) {
#line 266 "dtrsv.f"
			if (nounit) {
#line 266 "dtrsv.f"
			    x[jx] /= a[j + j * a_dim1];
#line 266 "dtrsv.f"
			}
#line 267 "dtrsv.f"
			temp = x[jx];
#line 268 "dtrsv.f"
			ix = jx;
#line 269 "dtrsv.f"
			i__2 = *n;
#line 269 "dtrsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 270 "dtrsv.f"
			    ix += *incx;
#line 271 "dtrsv.f"
			    x[ix] -= temp * a[i__ + j * a_dim1];
#line 272 "dtrsv.f"
/* L70: */
#line 272 "dtrsv.f"
			}
#line 273 "dtrsv.f"
		    }
#line 274 "dtrsv.f"
		    jx += *incx;
#line 275 "dtrsv.f"
/* L80: */
#line 275 "dtrsv.f"
		}
#line 276 "dtrsv.f"
	    }
#line 277 "dtrsv.f"
	}
#line 278 "dtrsv.f"
    } else {

/*        Form  x := inv( A**T )*x. */

#line 282 "dtrsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 283 "dtrsv.f"
	    if (*incx == 1) {
#line 284 "dtrsv.f"
		i__1 = *n;
#line 284 "dtrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 285 "dtrsv.f"
		    temp = x[j];
#line 286 "dtrsv.f"
		    i__2 = j - 1;
#line 286 "dtrsv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 287 "dtrsv.f"
			temp -= a[i__ + j * a_dim1] * x[i__];
#line 288 "dtrsv.f"
/* L90: */
#line 288 "dtrsv.f"
		    }
#line 289 "dtrsv.f"
		    if (nounit) {
#line 289 "dtrsv.f"
			temp /= a[j + j * a_dim1];
#line 289 "dtrsv.f"
		    }
#line 290 "dtrsv.f"
		    x[j] = temp;
#line 291 "dtrsv.f"
/* L100: */
#line 291 "dtrsv.f"
		}
#line 292 "dtrsv.f"
	    } else {
#line 293 "dtrsv.f"
		jx = kx;
#line 294 "dtrsv.f"
		i__1 = *n;
#line 294 "dtrsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 295 "dtrsv.f"
		    temp = x[jx];
#line 296 "dtrsv.f"
		    ix = kx;
#line 297 "dtrsv.f"
		    i__2 = j - 1;
#line 297 "dtrsv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 298 "dtrsv.f"
			temp -= a[i__ + j * a_dim1] * x[ix];
#line 299 "dtrsv.f"
			ix += *incx;
#line 300 "dtrsv.f"
/* L110: */
#line 300 "dtrsv.f"
		    }
#line 301 "dtrsv.f"
		    if (nounit) {
#line 301 "dtrsv.f"
			temp /= a[j + j * a_dim1];
#line 301 "dtrsv.f"
		    }
#line 302 "dtrsv.f"
		    x[jx] = temp;
#line 303 "dtrsv.f"
		    jx += *incx;
#line 304 "dtrsv.f"
/* L120: */
#line 304 "dtrsv.f"
		}
#line 305 "dtrsv.f"
	    }
#line 306 "dtrsv.f"
	} else {
#line 307 "dtrsv.f"
	    if (*incx == 1) {
#line 308 "dtrsv.f"
		for (j = *n; j >= 1; --j) {
#line 309 "dtrsv.f"
		    temp = x[j];
#line 310 "dtrsv.f"
		    i__1 = j + 1;
#line 310 "dtrsv.f"
		    for (i__ = *n; i__ >= i__1; --i__) {
#line 311 "dtrsv.f"
			temp -= a[i__ + j * a_dim1] * x[i__];
#line 312 "dtrsv.f"
/* L130: */
#line 312 "dtrsv.f"
		    }
#line 313 "dtrsv.f"
		    if (nounit) {
#line 313 "dtrsv.f"
			temp /= a[j + j * a_dim1];
#line 313 "dtrsv.f"
		    }
#line 314 "dtrsv.f"
		    x[j] = temp;
#line 315 "dtrsv.f"
/* L140: */
#line 315 "dtrsv.f"
		}
#line 316 "dtrsv.f"
	    } else {
#line 317 "dtrsv.f"
		kx += (*n - 1) * *incx;
#line 318 "dtrsv.f"
		jx = kx;
#line 319 "dtrsv.f"
		for (j = *n; j >= 1; --j) {
#line 320 "dtrsv.f"
		    temp = x[jx];
#line 321 "dtrsv.f"
		    ix = kx;
#line 322 "dtrsv.f"
		    i__1 = j + 1;
#line 322 "dtrsv.f"
		    for (i__ = *n; i__ >= i__1; --i__) {
#line 323 "dtrsv.f"
			temp -= a[i__ + j * a_dim1] * x[ix];
#line 324 "dtrsv.f"
			ix -= *incx;
#line 325 "dtrsv.f"
/* L150: */
#line 325 "dtrsv.f"
		    }
#line 326 "dtrsv.f"
		    if (nounit) {
#line 326 "dtrsv.f"
			temp /= a[j + j * a_dim1];
#line 326 "dtrsv.f"
		    }
#line 327 "dtrsv.f"
		    x[jx] = temp;
#line 328 "dtrsv.f"
		    jx -= *incx;
#line 329 "dtrsv.f"
/* L160: */
#line 329 "dtrsv.f"
		}
#line 330 "dtrsv.f"
	    }
#line 331 "dtrsv.f"
	}
#line 332 "dtrsv.f"
    }

#line 334 "dtrsv.f"
    return 0;

/*     End of DTRSV . */

} /* dtrsv_ */

