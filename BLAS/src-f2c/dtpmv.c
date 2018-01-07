#line 1 "dtpmv.f"
/* dtpmv.f -- translated by f2c (version 20100827).
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

#line 1 "dtpmv.f"
/* > \brief \b DTPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

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
/* > DTPMV  performs one of the matrix-vector operations */
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
/* >          AP is DOUBLE PRECISION array, dimension at least */
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
/* >          X is DOUBLE PRECISION array, dimension at least */
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
/* Subroutine */ int dtpmv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 179 "dtpmv.f"
    /* Parameter adjustments */
#line 179 "dtpmv.f"
    --x;
#line 179 "dtpmv.f"
    --ap;
#line 179 "dtpmv.f"

#line 179 "dtpmv.f"
    /* Function Body */
#line 179 "dtpmv.f"
    info = 0;
#line 180 "dtpmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 181 "dtpmv.f"
	info = 1;
#line 182 "dtpmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 184 "dtpmv.f"
	info = 2;
#line 185 "dtpmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 186 "dtpmv.f"
	info = 3;
#line 187 "dtpmv.f"
    } else if (*n < 0) {
#line 188 "dtpmv.f"
	info = 4;
#line 189 "dtpmv.f"
    } else if (*incx == 0) {
#line 190 "dtpmv.f"
	info = 7;
#line 191 "dtpmv.f"
    }
#line 192 "dtpmv.f"
    if (info != 0) {
#line 193 "dtpmv.f"
	xerbla_("DTPMV ", &info, (ftnlen)6);
#line 194 "dtpmv.f"
	return 0;
#line 195 "dtpmv.f"
    }

/*     Quick return if possible. */

#line 199 "dtpmv.f"
    if (*n == 0) {
#line 199 "dtpmv.f"
	return 0;
#line 199 "dtpmv.f"
    }

#line 201 "dtpmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 206 "dtpmv.f"
    if (*incx <= 0) {
#line 207 "dtpmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 208 "dtpmv.f"
    } else if (*incx != 1) {
#line 209 "dtpmv.f"
	kx = 1;
#line 210 "dtpmv.f"
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

#line 215 "dtpmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x:= A*x. */

#line 219 "dtpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 220 "dtpmv.f"
	    kk = 1;
#line 221 "dtpmv.f"
	    if (*incx == 1) {
#line 222 "dtpmv.f"
		i__1 = *n;
#line 222 "dtpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 223 "dtpmv.f"
		    if (x[j] != 0.) {
#line 224 "dtpmv.f"
			temp = x[j];
#line 225 "dtpmv.f"
			k = kk;
#line 226 "dtpmv.f"
			i__2 = j - 1;
#line 226 "dtpmv.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 227 "dtpmv.f"
			    x[i__] += temp * ap[k];
#line 228 "dtpmv.f"
			    ++k;
#line 229 "dtpmv.f"
/* L10: */
#line 229 "dtpmv.f"
			}
#line 230 "dtpmv.f"
			if (nounit) {
#line 230 "dtpmv.f"
			    x[j] *= ap[kk + j - 1];
#line 230 "dtpmv.f"
			}
#line 231 "dtpmv.f"
		    }
#line 232 "dtpmv.f"
		    kk += j;
#line 233 "dtpmv.f"
/* L20: */
#line 233 "dtpmv.f"
		}
#line 234 "dtpmv.f"
	    } else {
#line 235 "dtpmv.f"
		jx = kx;
#line 236 "dtpmv.f"
		i__1 = *n;
#line 236 "dtpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 237 "dtpmv.f"
		    if (x[jx] != 0.) {
#line 238 "dtpmv.f"
			temp = x[jx];
#line 239 "dtpmv.f"
			ix = kx;
#line 240 "dtpmv.f"
			i__2 = kk + j - 2;
#line 240 "dtpmv.f"
			for (k = kk; k <= i__2; ++k) {
#line 241 "dtpmv.f"
			    x[ix] += temp * ap[k];
#line 242 "dtpmv.f"
			    ix += *incx;
#line 243 "dtpmv.f"
/* L30: */
#line 243 "dtpmv.f"
			}
#line 244 "dtpmv.f"
			if (nounit) {
#line 244 "dtpmv.f"
			    x[jx] *= ap[kk + j - 1];
#line 244 "dtpmv.f"
			}
#line 245 "dtpmv.f"
		    }
#line 246 "dtpmv.f"
		    jx += *incx;
#line 247 "dtpmv.f"
		    kk += j;
#line 248 "dtpmv.f"
/* L40: */
#line 248 "dtpmv.f"
		}
#line 249 "dtpmv.f"
	    }
#line 250 "dtpmv.f"
	} else {
#line 251 "dtpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 252 "dtpmv.f"
	    if (*incx == 1) {
#line 253 "dtpmv.f"
		for (j = *n; j >= 1; --j) {
#line 254 "dtpmv.f"
		    if (x[j] != 0.) {
#line 255 "dtpmv.f"
			temp = x[j];
#line 256 "dtpmv.f"
			k = kk;
#line 257 "dtpmv.f"
			i__1 = j + 1;
#line 257 "dtpmv.f"
			for (i__ = *n; i__ >= i__1; --i__) {
#line 258 "dtpmv.f"
			    x[i__] += temp * ap[k];
#line 259 "dtpmv.f"
			    --k;
#line 260 "dtpmv.f"
/* L50: */
#line 260 "dtpmv.f"
			}
#line 261 "dtpmv.f"
			if (nounit) {
#line 261 "dtpmv.f"
			    x[j] *= ap[kk - *n + j];
#line 261 "dtpmv.f"
			}
#line 262 "dtpmv.f"
		    }
#line 263 "dtpmv.f"
		    kk -= *n - j + 1;
#line 264 "dtpmv.f"
/* L60: */
#line 264 "dtpmv.f"
		}
#line 265 "dtpmv.f"
	    } else {
#line 266 "dtpmv.f"
		kx += (*n - 1) * *incx;
#line 267 "dtpmv.f"
		jx = kx;
#line 268 "dtpmv.f"
		for (j = *n; j >= 1; --j) {
#line 269 "dtpmv.f"
		    if (x[jx] != 0.) {
#line 270 "dtpmv.f"
			temp = x[jx];
#line 271 "dtpmv.f"
			ix = kx;
#line 272 "dtpmv.f"
			i__1 = kk - (*n - (j + 1));
#line 272 "dtpmv.f"
			for (k = kk; k >= i__1; --k) {
#line 273 "dtpmv.f"
			    x[ix] += temp * ap[k];
#line 274 "dtpmv.f"
			    ix -= *incx;
#line 275 "dtpmv.f"
/* L70: */
#line 275 "dtpmv.f"
			}
#line 276 "dtpmv.f"
			if (nounit) {
#line 276 "dtpmv.f"
			    x[jx] *= ap[kk - *n + j];
#line 276 "dtpmv.f"
			}
#line 277 "dtpmv.f"
		    }
#line 278 "dtpmv.f"
		    jx -= *incx;
#line 279 "dtpmv.f"
		    kk -= *n - j + 1;
#line 280 "dtpmv.f"
/* L80: */
#line 280 "dtpmv.f"
		}
#line 281 "dtpmv.f"
	    }
#line 282 "dtpmv.f"
	}
#line 283 "dtpmv.f"
    } else {

/*        Form  x := A**T*x. */

#line 287 "dtpmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 288 "dtpmv.f"
	    kk = *n * (*n + 1) / 2;
#line 289 "dtpmv.f"
	    if (*incx == 1) {
#line 290 "dtpmv.f"
		for (j = *n; j >= 1; --j) {
#line 291 "dtpmv.f"
		    temp = x[j];
#line 292 "dtpmv.f"
		    if (nounit) {
#line 292 "dtpmv.f"
			temp *= ap[kk];
#line 292 "dtpmv.f"
		    }
#line 293 "dtpmv.f"
		    k = kk - 1;
#line 294 "dtpmv.f"
		    for (i__ = j - 1; i__ >= 1; --i__) {
#line 295 "dtpmv.f"
			temp += ap[k] * x[i__];
#line 296 "dtpmv.f"
			--k;
#line 297 "dtpmv.f"
/* L90: */
#line 297 "dtpmv.f"
		    }
#line 298 "dtpmv.f"
		    x[j] = temp;
#line 299 "dtpmv.f"
		    kk -= j;
#line 300 "dtpmv.f"
/* L100: */
#line 300 "dtpmv.f"
		}
#line 301 "dtpmv.f"
	    } else {
#line 302 "dtpmv.f"
		jx = kx + (*n - 1) * *incx;
#line 303 "dtpmv.f"
		for (j = *n; j >= 1; --j) {
#line 304 "dtpmv.f"
		    temp = x[jx];
#line 305 "dtpmv.f"
		    ix = jx;
#line 306 "dtpmv.f"
		    if (nounit) {
#line 306 "dtpmv.f"
			temp *= ap[kk];
#line 306 "dtpmv.f"
		    }
#line 307 "dtpmv.f"
		    i__1 = kk - j + 1;
#line 307 "dtpmv.f"
		    for (k = kk - 1; k >= i__1; --k) {
#line 308 "dtpmv.f"
			ix -= *incx;
#line 309 "dtpmv.f"
			temp += ap[k] * x[ix];
#line 310 "dtpmv.f"
/* L110: */
#line 310 "dtpmv.f"
		    }
#line 311 "dtpmv.f"
		    x[jx] = temp;
#line 312 "dtpmv.f"
		    jx -= *incx;
#line 313 "dtpmv.f"
		    kk -= j;
#line 314 "dtpmv.f"
/* L120: */
#line 314 "dtpmv.f"
		}
#line 315 "dtpmv.f"
	    }
#line 316 "dtpmv.f"
	} else {
#line 317 "dtpmv.f"
	    kk = 1;
#line 318 "dtpmv.f"
	    if (*incx == 1) {
#line 319 "dtpmv.f"
		i__1 = *n;
#line 319 "dtpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 320 "dtpmv.f"
		    temp = x[j];
#line 321 "dtpmv.f"
		    if (nounit) {
#line 321 "dtpmv.f"
			temp *= ap[kk];
#line 321 "dtpmv.f"
		    }
#line 322 "dtpmv.f"
		    k = kk + 1;
#line 323 "dtpmv.f"
		    i__2 = *n;
#line 323 "dtpmv.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 324 "dtpmv.f"
			temp += ap[k] * x[i__];
#line 325 "dtpmv.f"
			++k;
#line 326 "dtpmv.f"
/* L130: */
#line 326 "dtpmv.f"
		    }
#line 327 "dtpmv.f"
		    x[j] = temp;
#line 328 "dtpmv.f"
		    kk += *n - j + 1;
#line 329 "dtpmv.f"
/* L140: */
#line 329 "dtpmv.f"
		}
#line 330 "dtpmv.f"
	    } else {
#line 331 "dtpmv.f"
		jx = kx;
#line 332 "dtpmv.f"
		i__1 = *n;
#line 332 "dtpmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 333 "dtpmv.f"
		    temp = x[jx];
#line 334 "dtpmv.f"
		    ix = jx;
#line 335 "dtpmv.f"
		    if (nounit) {
#line 335 "dtpmv.f"
			temp *= ap[kk];
#line 335 "dtpmv.f"
		    }
#line 336 "dtpmv.f"
		    i__2 = kk + *n - j;
#line 336 "dtpmv.f"
		    for (k = kk + 1; k <= i__2; ++k) {
#line 337 "dtpmv.f"
			ix += *incx;
#line 338 "dtpmv.f"
			temp += ap[k] * x[ix];
#line 339 "dtpmv.f"
/* L150: */
#line 339 "dtpmv.f"
		    }
#line 340 "dtpmv.f"
		    x[jx] = temp;
#line 341 "dtpmv.f"
		    jx += *incx;
#line 342 "dtpmv.f"
		    kk += *n - j + 1;
#line 343 "dtpmv.f"
/* L160: */
#line 343 "dtpmv.f"
		}
#line 344 "dtpmv.f"
	    }
#line 345 "dtpmv.f"
	}
#line 346 "dtpmv.f"
    }

#line 348 "dtpmv.f"
    return 0;

/*     End of DTPMV . */

} /* dtpmv_ */

