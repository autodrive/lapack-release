#line 1 "dtpsv.f"
/* dtpsv.f -- translated by f2c (version 20100827).
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

#line 1 "dtpsv.f"
/* > \brief \b DTPSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular matrix, supplied in packed form. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular matrix packed sequentially, */
/* >           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/* >           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/* >           respectively, and so on. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular matrix packed sequentially, */
/* >           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/* >           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/* >           respectively, and so on. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced, but are assumed to be unity. */
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
/* Subroutine */ int dtpsv_(char *uplo, char *trans, char *diag, integer *n, 
	doublereal *ap, doublereal *x, integer *incx, ftnlen uplo_len, ftnlen 
	trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kk, ix, jx, kx, info;
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

/*     Test the input parameters. */

#line 181 "dtpsv.f"
    /* Parameter adjustments */
#line 181 "dtpsv.f"
    --x;
#line 181 "dtpsv.f"
    --ap;
#line 181 "dtpsv.f"

#line 181 "dtpsv.f"
    /* Function Body */
#line 181 "dtpsv.f"
    info = 0;
#line 182 "dtpsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 183 "dtpsv.f"
	info = 1;
#line 184 "dtpsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 186 "dtpsv.f"
	info = 2;
#line 187 "dtpsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 188 "dtpsv.f"
	info = 3;
#line 189 "dtpsv.f"
    } else if (*n < 0) {
#line 190 "dtpsv.f"
	info = 4;
#line 191 "dtpsv.f"
    } else if (*incx == 0) {
#line 192 "dtpsv.f"
	info = 7;
#line 193 "dtpsv.f"
    }
#line 194 "dtpsv.f"
    if (info != 0) {
#line 195 "dtpsv.f"
	xerbla_("DTPSV ", &info, (ftnlen)6);
#line 196 "dtpsv.f"
	return 0;
#line 197 "dtpsv.f"
    }

/*     Quick return if possible. */

#line 201 "dtpsv.f"
    if (*n == 0) {
#line 201 "dtpsv.f"
	return 0;
#line 201 "dtpsv.f"
    }

#line 203 "dtpsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 208 "dtpsv.f"
    if (*incx <= 0) {
#line 209 "dtpsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 210 "dtpsv.f"
    } else if (*incx != 1) {
#line 211 "dtpsv.f"
	kx = 1;
#line 212 "dtpsv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 217 "dtpsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 221 "dtpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 222 "dtpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 223 "dtpsv.f"
	    if (*incx == 1) {
#line 224 "dtpsv.f"
		for (j = *n; j >= 1; --j) {
#line 225 "dtpsv.f"
		    if (x[j] != 0.) {
#line 226 "dtpsv.f"
			if (nounit) {
#line 226 "dtpsv.f"
			    x[j] /= ap[kk];
#line 226 "dtpsv.f"
			}
#line 227 "dtpsv.f"
			temp = x[j];
#line 228 "dtpsv.f"
			k = kk - 1;
#line 229 "dtpsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 230 "dtpsv.f"
			    x[i__] -= temp * ap[k];
#line 231 "dtpsv.f"
			    --k;
#line 232 "dtpsv.f"
/* L10: */
#line 232 "dtpsv.f"
			}
#line 233 "dtpsv.f"
		    }
#line 234 "dtpsv.f"
		    kk -= j;
#line 235 "dtpsv.f"
/* L20: */
#line 235 "dtpsv.f"
		}
#line 236 "dtpsv.f"
	    } else {
#line 237 "dtpsv.f"
		jx = kx + (*n - 1) * *incx;
#line 238 "dtpsv.f"
		for (j = *n; j >= 1; --j) {
#line 239 "dtpsv.f"
		    if (x[jx] != 0.) {
#line 240 "dtpsv.f"
			if (nounit) {
#line 240 "dtpsv.f"
			    x[jx] /= ap[kk];
#line 240 "dtpsv.f"
			}
#line 241 "dtpsv.f"
			temp = x[jx];
#line 242 "dtpsv.f"
			ix = jx;
#line 243 "dtpsv.f"
			i__1 = kk - j + 1;
#line 243 "dtpsv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 244 "dtpsv.f"
			    ix -= *incx;
#line 245 "dtpsv.f"
			    x[ix] -= temp * ap[k];
#line 246 "dtpsv.f"
/* L30: */
#line 246 "dtpsv.f"
			}
#line 247 "dtpsv.f"
		    }
#line 248 "dtpsv.f"
		    jx -= *incx;
#line 249 "dtpsv.f"
		    kk -= j;
#line 250 "dtpsv.f"
/* L40: */
#line 250 "dtpsv.f"
		}
#line 251 "dtpsv.f"
	    }
#line 252 "dtpsv.f"
	} else {
#line 253 "dtpsv.f"
	    kk = 1;
#line 254 "dtpsv.f"
	    if (*incx == 1) {
#line 255 "dtpsv.f"
		i__1 = *n;
#line 255 "dtpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 256 "dtpsv.f"
		    if (x[j] != 0.) {
#line 257 "dtpsv.f"
			if (nounit) {
#line 257 "dtpsv.f"
			    x[j] /= ap[kk];
#line 257 "dtpsv.f"
			}
#line 258 "dtpsv.f"
			temp = x[j];
#line 259 "dtpsv.f"
			k = kk + 1;
#line 260 "dtpsv.f"
			i__2 = *n;
#line 260 "dtpsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 261 "dtpsv.f"
			    x[i__] -= temp * ap[k];
#line 262 "dtpsv.f"
			    ++k;
#line 263 "dtpsv.f"
/* L50: */
#line 263 "dtpsv.f"
			}
#line 264 "dtpsv.f"
		    }
#line 265 "dtpsv.f"
		    kk += *n - j + 1;
#line 266 "dtpsv.f"
/* L60: */
#line 266 "dtpsv.f"
		}
#line 267 "dtpsv.f"
	    } else {
#line 268 "dtpsv.f"
		jx = kx;
#line 269 "dtpsv.f"
		i__1 = *n;
#line 269 "dtpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 270 "dtpsv.f"
		    if (x[jx] != 0.) {
#line 271 "dtpsv.f"
			if (nounit) {
#line 271 "dtpsv.f"
			    x[jx] /= ap[kk];
#line 271 "dtpsv.f"
			}
#line 272 "dtpsv.f"
			temp = x[jx];
#line 273 "dtpsv.f"
			ix = jx;
#line 274 "dtpsv.f"
			i__2 = kk + *n - j;
#line 274 "dtpsv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 275 "dtpsv.f"
			    ix += *incx;
#line 276 "dtpsv.f"
			    x[ix] -= temp * ap[k];
#line 277 "dtpsv.f"
/* L70: */
#line 277 "dtpsv.f"
			}
#line 278 "dtpsv.f"
		    }
#line 279 "dtpsv.f"
		    jx += *incx;
#line 280 "dtpsv.f"
		    kk += *n - j + 1;
#line 281 "dtpsv.f"
/* L80: */
#line 281 "dtpsv.f"
		}
#line 282 "dtpsv.f"
	    }
#line 283 "dtpsv.f"
	}
#line 284 "dtpsv.f"
    } else {

/*        Form  x := inv( A**T )*x. */

#line 288 "dtpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 289 "dtpsv.f"
	    kk = 1;
#line 290 "dtpsv.f"
	    if (*incx == 1) {
#line 291 "dtpsv.f"
		i__1 = *n;
#line 291 "dtpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "dtpsv.f"
		    temp = x[j];
#line 293 "dtpsv.f"
		    k = kk;
#line 294 "dtpsv.f"
		    i__2 = j - 1;
#line 294 "dtpsv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "dtpsv.f"
			temp -= ap[k] * x[i__];
#line 296 "dtpsv.f"
			++k;
#line 297 "dtpsv.f"
/* L90: */
#line 297 "dtpsv.f"
		    }
#line 298 "dtpsv.f"
		    if (nounit) {
#line 298 "dtpsv.f"
			temp /= ap[kk + j - 1];
#line 298 "dtpsv.f"
		    }
#line 299 "dtpsv.f"
		    x[j] = temp;
#line 300 "dtpsv.f"
		    kk += j;
#line 301 "dtpsv.f"
/* L100: */
#line 301 "dtpsv.f"
		}
#line 302 "dtpsv.f"
	    } else {
#line 303 "dtpsv.f"
		jx = kx;
#line 304 "dtpsv.f"
		i__1 = *n;
#line 304 "dtpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 305 "dtpsv.f"
		    temp = x[jx];
#line 306 "dtpsv.f"
		    ix = kx;
#line 307 "dtpsv.f"
		    i__2 = kk + j - 2;
#line 307 "dtpsv.f"
		    for (k = kk; k <= i__2; ++k) {
#line 308 "dtpsv.f"
			temp -= ap[k] * x[ix];
#line 309 "dtpsv.f"
			ix += *incx;
#line 310 "dtpsv.f"
/* L110: */
#line 310 "dtpsv.f"
		    }
#line 311 "dtpsv.f"
		    if (nounit) {
#line 311 "dtpsv.f"
			temp /= ap[kk + j - 1];
#line 311 "dtpsv.f"
		    }
#line 312 "dtpsv.f"
		    x[jx] = temp;
#line 313 "dtpsv.f"
		    jx += *incx;
#line 314 "dtpsv.f"
		    kk += j;
#line 315 "dtpsv.f"
/* L120: */
#line 315 "dtpsv.f"
		}
#line 316 "dtpsv.f"
	    }
#line 317 "dtpsv.f"
	} else {
#line 318 "dtpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 319 "dtpsv.f"
	    if (*incx == 1) {
#line 320 "dtpsv.f"
		for (j = *n; j >= 1; --j) {
#line 321 "dtpsv.f"
		    temp = x[j];
#line 322 "dtpsv.f"
		    k = kk;
#line 323 "dtpsv.f"
		    i__1 = j + 1;
#line 323 "dtpsv.f"
		    for (i__ = *n; i__ >= i__1; --i__) {
#line 324 "dtpsv.f"
			temp -= ap[k] * x[i__];
#line 325 "dtpsv.f"
			--k;
#line 326 "dtpsv.f"
/* L130: */
#line 326 "dtpsv.f"
		    }
#line 327 "dtpsv.f"
		    if (nounit) {
#line 327 "dtpsv.f"
			temp /= ap[kk - *n + j];
#line 327 "dtpsv.f"
		    }
#line 328 "dtpsv.f"
		    x[j] = temp;
#line 329 "dtpsv.f"
		    kk -= *n - j + 1;
#line 330 "dtpsv.f"
/* L140: */
#line 330 "dtpsv.f"
		}
#line 331 "dtpsv.f"
	    } else {
#line 332 "dtpsv.f"
		kx += (*n - 1) * *incx;
#line 333 "dtpsv.f"
		jx = kx;
#line 334 "dtpsv.f"
		for (j = *n; j >= 1; --j) {
#line 335 "dtpsv.f"
		    temp = x[jx];
#line 336 "dtpsv.f"
		    ix = kx;
#line 337 "dtpsv.f"
		    i__1 = kk - (*n - (j + 1));
#line 337 "dtpsv.f"
		    for (k = kk; k >= i__1; --k) {
#line 338 "dtpsv.f"
			temp -= ap[k] * x[ix];
#line 339 "dtpsv.f"
			ix -= *incx;
#line 340 "dtpsv.f"
/* L150: */
#line 340 "dtpsv.f"
		    }
#line 341 "dtpsv.f"
		    if (nounit) {
#line 341 "dtpsv.f"
			temp /= ap[kk - *n + j];
#line 341 "dtpsv.f"
		    }
#line 342 "dtpsv.f"
		    x[jx] = temp;
#line 343 "dtpsv.f"
		    jx -= *incx;
#line 344 "dtpsv.f"
		    kk -= *n - j + 1;
#line 345 "dtpsv.f"
/* L160: */
#line 345 "dtpsv.f"
		}
#line 346 "dtpsv.f"
	    }
#line 347 "dtpsv.f"
	}
#line 348 "dtpsv.f"
    }

#line 350 "dtpsv.f"
    return 0;

/*     End of DTPSV . */

} /* dtpsv_ */

