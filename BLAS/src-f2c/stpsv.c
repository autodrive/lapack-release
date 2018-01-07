#line 1 "stpsv.f"
/* stpsv.f -- translated by f2c (version 20100827).
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

#line 1 "stpsv.f"
/* > \brief \b STPSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPSV  solves one of the systems of equations */
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
/* >          AP is REAL array of DIMENSION at least */
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
/* Subroutine */ int stpsv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 181 "stpsv.f"
    /* Parameter adjustments */
#line 181 "stpsv.f"
    --x;
#line 181 "stpsv.f"
    --ap;
#line 181 "stpsv.f"

#line 181 "stpsv.f"
    /* Function Body */
#line 181 "stpsv.f"
    info = 0;
#line 182 "stpsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 183 "stpsv.f"
	info = 1;
#line 184 "stpsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 186 "stpsv.f"
	info = 2;
#line 187 "stpsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 188 "stpsv.f"
	info = 3;
#line 189 "stpsv.f"
    } else if (*n < 0) {
#line 190 "stpsv.f"
	info = 4;
#line 191 "stpsv.f"
    } else if (*incx == 0) {
#line 192 "stpsv.f"
	info = 7;
#line 193 "stpsv.f"
    }
#line 194 "stpsv.f"
    if (info != 0) {
#line 195 "stpsv.f"
	xerbla_("STPSV ", &info, (ftnlen)6);
#line 196 "stpsv.f"
	return 0;
#line 197 "stpsv.f"
    }

/*     Quick return if possible. */

#line 201 "stpsv.f"
    if (*n == 0) {
#line 201 "stpsv.f"
	return 0;
#line 201 "stpsv.f"
    }

#line 203 "stpsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 208 "stpsv.f"
    if (*incx <= 0) {
#line 209 "stpsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 210 "stpsv.f"
    } else if (*incx != 1) {
#line 211 "stpsv.f"
	kx = 1;
#line 212 "stpsv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 217 "stpsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 221 "stpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 222 "stpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 223 "stpsv.f"
	    if (*incx == 1) {
#line 224 "stpsv.f"
		for (j = *n; j >= 1; --j) {
#line 225 "stpsv.f"
		    if (x[j] != 0.) {
#line 226 "stpsv.f"
			if (nounit) {
#line 226 "stpsv.f"
			    x[j] /= ap[kk];
#line 226 "stpsv.f"
			}
#line 227 "stpsv.f"
			temp = x[j];
#line 228 "stpsv.f"
			k = kk - 1;
#line 229 "stpsv.f"
			for (i__ = j - 1; i__ >= 1; --i__) {
#line 230 "stpsv.f"
			    x[i__] -= temp * ap[k];
#line 231 "stpsv.f"
			    --k;
#line 232 "stpsv.f"
/* L10: */
#line 232 "stpsv.f"
			}
#line 233 "stpsv.f"
		    }
#line 234 "stpsv.f"
		    kk -= j;
#line 235 "stpsv.f"
/* L20: */
#line 235 "stpsv.f"
		}
#line 236 "stpsv.f"
	    } else {
#line 237 "stpsv.f"
		jx = kx + (*n - 1) * *incx;
#line 238 "stpsv.f"
		for (j = *n; j >= 1; --j) {
#line 239 "stpsv.f"
		    if (x[jx] != 0.) {
#line 240 "stpsv.f"
			if (nounit) {
#line 240 "stpsv.f"
			    x[jx] /= ap[kk];
#line 240 "stpsv.f"
			}
#line 241 "stpsv.f"
			temp = x[jx];
#line 242 "stpsv.f"
			ix = jx;
#line 243 "stpsv.f"
			i__1 = kk - j + 1;
#line 243 "stpsv.f"
			for (k = kk - 1; k >= i__1; --k) {
#line 244 "stpsv.f"
			    ix -= *incx;
#line 245 "stpsv.f"
			    x[ix] -= temp * ap[k];
#line 246 "stpsv.f"
/* L30: */
#line 246 "stpsv.f"
			}
#line 247 "stpsv.f"
		    }
#line 248 "stpsv.f"
		    jx -= *incx;
#line 249 "stpsv.f"
		    kk -= j;
#line 250 "stpsv.f"
/* L40: */
#line 250 "stpsv.f"
		}
#line 251 "stpsv.f"
	    }
#line 252 "stpsv.f"
	} else {
#line 253 "stpsv.f"
	    kk = 1;
#line 254 "stpsv.f"
	    if (*incx == 1) {
#line 255 "stpsv.f"
		i__1 = *n;
#line 255 "stpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 256 "stpsv.f"
		    if (x[j] != 0.) {
#line 257 "stpsv.f"
			if (nounit) {
#line 257 "stpsv.f"
			    x[j] /= ap[kk];
#line 257 "stpsv.f"
			}
#line 258 "stpsv.f"
			temp = x[j];
#line 259 "stpsv.f"
			k = kk + 1;
#line 260 "stpsv.f"
			i__2 = *n;
#line 260 "stpsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 261 "stpsv.f"
			    x[i__] -= temp * ap[k];
#line 262 "stpsv.f"
			    ++k;
#line 263 "stpsv.f"
/* L50: */
#line 263 "stpsv.f"
			}
#line 264 "stpsv.f"
		    }
#line 265 "stpsv.f"
		    kk += *n - j + 1;
#line 266 "stpsv.f"
/* L60: */
#line 266 "stpsv.f"
		}
#line 267 "stpsv.f"
	    } else {
#line 268 "stpsv.f"
		jx = kx;
#line 269 "stpsv.f"
		i__1 = *n;
#line 269 "stpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 270 "stpsv.f"
		    if (x[jx] != 0.) {
#line 271 "stpsv.f"
			if (nounit) {
#line 271 "stpsv.f"
			    x[jx] /= ap[kk];
#line 271 "stpsv.f"
			}
#line 272 "stpsv.f"
			temp = x[jx];
#line 273 "stpsv.f"
			ix = jx;
#line 274 "stpsv.f"
			i__2 = kk + *n - j;
#line 274 "stpsv.f"
			for (k = kk + 1; k <= i__2; ++k) {
#line 275 "stpsv.f"
			    ix += *incx;
#line 276 "stpsv.f"
			    x[ix] -= temp * ap[k];
#line 277 "stpsv.f"
/* L70: */
#line 277 "stpsv.f"
			}
#line 278 "stpsv.f"
		    }
#line 279 "stpsv.f"
		    jx += *incx;
#line 280 "stpsv.f"
		    kk += *n - j + 1;
#line 281 "stpsv.f"
/* L80: */
#line 281 "stpsv.f"
		}
#line 282 "stpsv.f"
	    }
#line 283 "stpsv.f"
	}
#line 284 "stpsv.f"
    } else {

/*        Form  x := inv( A**T )*x. */

#line 288 "stpsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 289 "stpsv.f"
	    kk = 1;
#line 290 "stpsv.f"
	    if (*incx == 1) {
#line 291 "stpsv.f"
		i__1 = *n;
#line 291 "stpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 292 "stpsv.f"
		    temp = x[j];
#line 293 "stpsv.f"
		    k = kk;
#line 294 "stpsv.f"
		    i__2 = j - 1;
#line 294 "stpsv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 295 "stpsv.f"
			temp -= ap[k] * x[i__];
#line 296 "stpsv.f"
			++k;
#line 297 "stpsv.f"
/* L90: */
#line 297 "stpsv.f"
		    }
#line 298 "stpsv.f"
		    if (nounit) {
#line 298 "stpsv.f"
			temp /= ap[kk + j - 1];
#line 298 "stpsv.f"
		    }
#line 299 "stpsv.f"
		    x[j] = temp;
#line 300 "stpsv.f"
		    kk += j;
#line 301 "stpsv.f"
/* L100: */
#line 301 "stpsv.f"
		}
#line 302 "stpsv.f"
	    } else {
#line 303 "stpsv.f"
		jx = kx;
#line 304 "stpsv.f"
		i__1 = *n;
#line 304 "stpsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 305 "stpsv.f"
		    temp = x[jx];
#line 306 "stpsv.f"
		    ix = kx;
#line 307 "stpsv.f"
		    i__2 = kk + j - 2;
#line 307 "stpsv.f"
		    for (k = kk; k <= i__2; ++k) {
#line 308 "stpsv.f"
			temp -= ap[k] * x[ix];
#line 309 "stpsv.f"
			ix += *incx;
#line 310 "stpsv.f"
/* L110: */
#line 310 "stpsv.f"
		    }
#line 311 "stpsv.f"
		    if (nounit) {
#line 311 "stpsv.f"
			temp /= ap[kk + j - 1];
#line 311 "stpsv.f"
		    }
#line 312 "stpsv.f"
		    x[jx] = temp;
#line 313 "stpsv.f"
		    jx += *incx;
#line 314 "stpsv.f"
		    kk += j;
#line 315 "stpsv.f"
/* L120: */
#line 315 "stpsv.f"
		}
#line 316 "stpsv.f"
	    }
#line 317 "stpsv.f"
	} else {
#line 318 "stpsv.f"
	    kk = *n * (*n + 1) / 2;
#line 319 "stpsv.f"
	    if (*incx == 1) {
#line 320 "stpsv.f"
		for (j = *n; j >= 1; --j) {
#line 321 "stpsv.f"
		    temp = x[j];
#line 322 "stpsv.f"
		    k = kk;
#line 323 "stpsv.f"
		    i__1 = j + 1;
#line 323 "stpsv.f"
		    for (i__ = *n; i__ >= i__1; --i__) {
#line 324 "stpsv.f"
			temp -= ap[k] * x[i__];
#line 325 "stpsv.f"
			--k;
#line 326 "stpsv.f"
/* L130: */
#line 326 "stpsv.f"
		    }
#line 327 "stpsv.f"
		    if (nounit) {
#line 327 "stpsv.f"
			temp /= ap[kk - *n + j];
#line 327 "stpsv.f"
		    }
#line 328 "stpsv.f"
		    x[j] = temp;
#line 329 "stpsv.f"
		    kk -= *n - j + 1;
#line 330 "stpsv.f"
/* L140: */
#line 330 "stpsv.f"
		}
#line 331 "stpsv.f"
	    } else {
#line 332 "stpsv.f"
		kx += (*n - 1) * *incx;
#line 333 "stpsv.f"
		jx = kx;
#line 334 "stpsv.f"
		for (j = *n; j >= 1; --j) {
#line 335 "stpsv.f"
		    temp = x[jx];
#line 336 "stpsv.f"
		    ix = kx;
#line 337 "stpsv.f"
		    i__1 = kk - (*n - (j + 1));
#line 337 "stpsv.f"
		    for (k = kk; k >= i__1; --k) {
#line 338 "stpsv.f"
			temp -= ap[k] * x[ix];
#line 339 "stpsv.f"
			ix -= *incx;
#line 340 "stpsv.f"
/* L150: */
#line 340 "stpsv.f"
		    }
#line 341 "stpsv.f"
		    if (nounit) {
#line 341 "stpsv.f"
			temp /= ap[kk - *n + j];
#line 341 "stpsv.f"
		    }
#line 342 "stpsv.f"
		    x[jx] = temp;
#line 343 "stpsv.f"
		    jx -= *incx;
#line 344 "stpsv.f"
		    kk -= *n - j + 1;
#line 345 "stpsv.f"
/* L160: */
#line 345 "stpsv.f"
		}
#line 346 "stpsv.f"
	    }
#line 347 "stpsv.f"
	}
#line 348 "stpsv.f"
    }

#line 350 "stpsv.f"
    return 0;

/*     End of STPSV . */

} /* stpsv_ */

