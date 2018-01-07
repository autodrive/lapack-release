#line 1 "stpmv.f"
/* stpmv.f -- translated by f2c (version 20100827).
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

#line 1 "stpmv.f"
/* > \brief \b STPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

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
/* > STPMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix, supplied in packed form. */
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
/* Subroutine */ int stpmv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 179 "stpmv.f"
    /* Parameter adjustments */
#line 179 "stpmv.f"
    --x;
#line 179 "stpmv.f"
    --ap;
#line 179 "stpmv.f"

#line 179 "stpmv.f"
    /* Function Body */
#line 179 "stpmv.f"
    info = 0;
#line 180 "stpmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 181 "stpmv.f"
	info = 1;
#line 182 "stpmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 184 "stpmv.f"
	info = 2;
#line 185 "stpmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 186 "stpmv.f"
	info = 3;
#line 187 "stpmv.f"
    } else if (*n < 0) {
#line 188 "stpmv.f"
	info = 4;
#line 189 "stpmv.f"
    } else if (*incx == 0) {
#line 190 "stpmv.f"
	info = 7;
#line 191 "stpmv.f"
    }
#line 192 "stpmv.f"
    if (info != 0) {
#line 193 "stpmv.f"
	xerbla_("STPMV ", &info, (ftnlen)6);
#line 194 "stpmv.f"
	return 0;
#line 195 "stpmv.f"
    }

/*     Quick return if possible. */

#line 199 "stpmv.f"
    if (*n == 0) {
#line 199 "stpmv.f"
	return 0;
#line 199 "stpmv.f"
    }

#line 201 "stpmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 206 "stpmv.f"
    if (*incx <= 0) {
#line 207 "stpmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 208 "stpmv.f"
    } else if (*incx != 1) {
#line 209 "stpmv.f"
	kx = 1;
#line 210 "stpmv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 215 "stpmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x:= A*x. */

#line 219 "stpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 220 "stpmv.f"
	    kk = 1;
#line 221 "stpmv.f"
	    if (*incx == 1) {
#line 222 "stpmv.f"
		i__1 = *n;
#line 222 "stpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 223 "stpmv.f"
		    if (x[j] != 0.) {
#line 224 "stpmv.f"
			temp = x[j];
#line 225 "stpmv.f"
			k = kk;
#line 226 "stpmv.f"
			i__2 = j - 1;
#line 226 "stpmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 227 "stpmv.f"
			    x[i__] += temp * ap[k];
#line 228 "stpmv.f"
			    ++k;
#line 229 "stpmv.f"
/* L10: */
#line 229 "stpmv.f"
			}
#line 230 "stpmv.f"
			if (nounit) {
#line 230 "stpmv.f"
			    x[j] *= ap[kk + j - 1];
#line 230 "stpmv.f"
			}
#line 231 "stpmv.f"
		    }
#line 232 "stpmv.f"
		    kk += j;
#line 233 "stpmv.f"
/* L20: */
#line 233 "stpmv.f"
		}
#line 234 "stpmv.f"
	    } else {
#line 235 "stpmv.f"
		jx = kx;
#line 236 "stpmv.f"
		i__1 = *n;
#line 236 "stpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 237 "stpmv.f"
		    if (x[jx] != 0.) {
#line 238 "stpmv.f"
			temp = x[jx];
#line 239 "stpmv.f"
			ix = kx;
#line 240 "stpmv.f"
			i__2 = kk + j - 2;
#line 240 "stpmv.f"
			for (k = kk; k <= i__2; ++k) {
#line 241 "stpmv.f"
			    x[ix] += temp * ap[k];
#line 242 "stpmv.f"
			    ix += *incx;
#line 243 "stpmv.f"
/* L30: */
#line 243 "stpmv.f"
			}
#line 244 "stpmv.f"
			if (nounit) {
#line 244 "stpmv.f"
			    x[jx] *= ap[kk + j - 1];
#line 244 "stpmv.f"
			}
#line 245 "stpmv.f"
		    }
#line 246 "stpmv.f"
		    jx += *incx;
#line 247 "stpmv.f"
		    kk += j;
#line 248 "stpmv.f"
/* L40: */
#line 248 "stpmv.f"
		}
#line 249 "stpmv.f"
	    }
#line 250 "stpmv.f"
	} else {
#line 251 "stpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 252 "stpmv.f"
	    if (*incx == 1) {
#line 253 "stpmv.f"
		for (j = *n; j >= 1; --j) {
#line 254 "stpmv.f"
		    if (x[j] != 0.) {
#line 255 "stpmv.f"
			temp = x[j];
#line 256 "stpmv.f"
			k = kk;
#line 257 "stpmv.f"
			i__1 = j + 1;
#line 257 "stpmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 258 "stpmv.f"
			    x[i__] += temp * ap[k];
#line 259 "stpmv.f"
			    --k;
#line 260 "stpmv.f"
/* L50: */
#line 260 "stpmv.f"
			}
#line 261 "stpmv.f"
			if (nounit) {
#line 261 "stpmv.f"
			    x[j] *= ap[kk - *n + j];
#line 261 "stpmv.f"
			}
#line 262 "stpmv.f"
		    }
#line 263 "stpmv.f"
		    kk -= *n - j + 1;
#line 264 "stpmv.f"
/* L60: */
#line 264 "stpmv.f"
		}
#line 265 "stpmv.f"
	    } else {
#line 266 "stpmv.f"
		kx += (*n - 1) * *incx;
#line 267 "stpmv.f"
		jx = kx;
#line 268 "stpmv.f"
		for (j = *n; j >= 1; --j) {
#line 269 "stpmv.f"
		    if (x[jx] != 0.) {
#line 270 "stpmv.f"
			temp = x[jx];
#line 271 "stpmv.f"
			ix = kx;
#line 272 "stpmv.f"
			i__1 = kk - (*n - (j + 1));
#line 272 "stpmv.f"
			for (k = kk; k >= i__1; --k) {
#line 273 "stpmv.f"
			    x[ix] += temp * ap[k];
#line 274 "stpmv.f"
			    ix -= *incx;
#line 275 "stpmv.f"
/* L70: */
#line 275 "stpmv.f"
			}
#line 276 "stpmv.f"
			if (nounit) {
#line 276 "stpmv.f"
			    x[jx] *= ap[kk - *n + j];
#line 276 "stpmv.f"
			}
#line 277 "stpmv.f"
		    }
#line 278 "stpmv.f"
		    jx -= *incx;
#line 279 "stpmv.f"
		    kk -= *n - j + 1;
#line 280 "stpmv.f"
/* L80: */
#line 280 "stpmv.f"
		}
#line 281 "stpmv.f"
	    }
#line 282 "stpmv.f"
	}
#line 283 "stpmv.f"
    } else {

/*        Form  x := A**T*x. */

#line 287 "stpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 288 "stpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 289 "stpmv.f"
	    if (*incx == 1) {
#line 290 "stpmv.f"
		for (j = *n; j >= 1; --j) {
#line 291 "stpmv.f"
		    temp = x[j];
#line 292 "stpmv.f"
		    if (nounit) {
#line 292 "stpmv.f"
			temp *= ap[kk];
#line 292 "stpmv.f"
		    }
#line 293 "stpmv.f"
		    k = kk - 1;
#line 294 "stpmv.f"
		    for (i__ = j - 1; i__ >= 1; --i__) {
#line 295 "stpmv.f"
			temp += ap[k] * x[i__];
#line 296 "stpmv.f"
			--k;
#line 297 "stpmv.f"
/* L90: */
#line 297 "stpmv.f"
		    }
#line 298 "stpmv.f"
		    x[j] = temp;
#line 299 "stpmv.f"
		    kk -= j;
#line 300 "stpmv.f"
/* L100: */
#line 300 "stpmv.f"
		}
#line 301 "stpmv.f"
	    } else {
#line 302 "stpmv.f"
		jx = kx + (*n - 1) * *incx;
#line 303 "stpmv.f"
		for (j = *n; j >= 1; --j) {
#line 304 "stpmv.f"
		    temp = x[jx];
#line 305 "stpmv.f"
		    ix = jx;
#line 306 "stpmv.f"
		    if (nounit) {
#line 306 "stpmv.f"
			temp *= ap[kk];
#line 306 "stpmv.f"
		    }
#line 307 "stpmv.f"
		    i__1 = kk - j + 1;
#line 307 "stpmv.f"
		    for (k = kk - 1; k >= i__1; --k) {
#line 308 "stpmv.f"
			ix -= *incx;
#line 309 "stpmv.f"
			temp += ap[k] * x[ix];
#line 310 "stpmv.f"
/* L110: */
#line 310 "stpmv.f"
		    }
#line 311 "stpmv.f"
		    x[jx] = temp;
#line 312 "stpmv.f"
		    jx -= *incx;
#line 313 "stpmv.f"
		    kk -= j;
#line 314 "stpmv.f"
/* L120: */
#line 314 "stpmv.f"
		}
#line 315 "stpmv.f"
	    }
#line 316 "stpmv.f"
	} else {
#line 317 "stpmv.f"
	    kk = 1;
#line 318 "stpmv.f"
	    if (*incx == 1) {
#line 319 "stpmv.f"
		i__1 = *n;
#line 319 "stpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 320 "stpmv.f"
		    temp = x[j];
#line 321 "stpmv.f"
		    if (nounit) {
#line 321 "stpmv.f"
			temp *= ap[kk];
#line 321 "stpmv.f"
		    }
#line 322 "stpmv.f"
		    k = kk + 1;
#line 323 "stpmv.f"
		    i__2 = *n;
#line 323 "stpmv.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 324 "stpmv.f"
			temp += ap[k] * x[i__];
#line 325 "stpmv.f"
			++k;
#line 326 "stpmv.f"
/* L130: */
#line 326 "stpmv.f"
		    }
#line 327 "stpmv.f"
		    x[j] = temp;
#line 328 "stpmv.f"
		    kk += *n - j + 1;
#line 329 "stpmv.f"
/* L140: */
#line 329 "stpmv.f"
		}
#line 330 "stpmv.f"
	    } else {
#line 331 "stpmv.f"
		jx = kx;
#line 332 "stpmv.f"
		i__1 = *n;
#line 332 "stpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 333 "stpmv.f"
		    temp = x[jx];
#line 334 "stpmv.f"
		    ix = jx;
#line 335 "stpmv.f"
		    if (nounit) {
#line 335 "stpmv.f"
			temp *= ap[kk];
#line 335 "stpmv.f"
		    }
#line 336 "stpmv.f"
		    i__2 = kk + *n - j;
#line 336 "stpmv.f"
		    for (k = kk + 1; k <= i__2; ++k) {
#line 337 "stpmv.f"
			ix += *incx;
#line 338 "stpmv.f"
			temp += ap[k] * x[ix];
#line 339 "stpmv.f"
/* L150: */
#line 339 "stpmv.f"
		    }
#line 340 "stpmv.f"
		    x[jx] = temp;
#line 341 "stpmv.f"
		    jx += *incx;
#line 342 "stpmv.f"
		    kk += *n - j + 1;
#line 343 "stpmv.f"
/* L160: */
#line 343 "stpmv.f"
		}
#line 344 "stpmv.f"
	    }
#line 345 "stpmv.f"
	}
#line 346 "stpmv.f"
    }

#line 348 "stpmv.f"
    return 0;

/*     End of STPMV . */

} /* stpmv_ */

