#line 1 "dlasr.f"
/* dlasr.f -- translated by f2c (version 20100827).
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

#line 1 "dlasr.f"
/* > \brief \b DLASR applies a sequence of plane rotations to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasr.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasr.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasr.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, PIVOT, SIDE */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASR applies a sequence of plane rotations to a real matrix A, */
/* > from either the left or the right. */
/* > */
/* > When SIDE = 'L', the transformation takes the form */
/* > */
/* >    A := P*A */
/* > */
/* > and when SIDE = 'R', the transformation takes the form */
/* > */
/* >    A := A*P**T */
/* > */
/* > where P is an orthogonal matrix consisting of a sequence of z plane */
/* > rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R', */
/* > and P**T is the transpose of P. */
/* > */
/* > When DIRECT = 'F' (Forward sequence), then */
/* > */
/* >    P = P(z-1) * ... * P(2) * P(1) */
/* > */
/* > and when DIRECT = 'B' (Backward sequence), then */
/* > */
/* >    P = P(1) * P(2) * ... * P(z-1) */
/* > */
/* > where P(k) is a plane rotation matrix defined by the 2-by-2 rotation */
/* > */
/* >    R(k) = (  c(k)  s(k) ) */
/* >         = ( -s(k)  c(k) ). */
/* > */
/* > When PIVOT = 'V' (Variable pivot), the rotation is performed */
/* > for the plane (k,k+1), i.e., P(k) has the form */
/* > */
/* >    P(k) = (  1                                            ) */
/* >           (       ...                                     ) */
/* >           (              1                                ) */
/* >           (                   c(k)  s(k)                  ) */
/* >           (                  -s(k)  c(k)                  ) */
/* >           (                                1              ) */
/* >           (                                     ...       ) */
/* >           (                                            1  ) */
/* > */
/* > where R(k) appears as a rank-2 modification to the identity matrix in */
/* > rows and columns k and k+1. */
/* > */
/* > When PIVOT = 'T' (Top pivot), the rotation is performed for the */
/* > plane (1,k+1), so P(k) has the form */
/* > */
/* >    P(k) = (  c(k)                    s(k)                 ) */
/* >           (         1                                     ) */
/* >           (              ...                              ) */
/* >           (                     1                         ) */
/* >           ( -s(k)                    c(k)                 ) */
/* >           (                                 1             ) */
/* >           (                                      ...      ) */
/* >           (                                             1 ) */
/* > */
/* > where R(k) appears in rows and columns 1 and k+1. */
/* > */
/* > Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is */
/* > performed for the plane (k,z), giving P(k) the form */
/* > */
/* >    P(k) = ( 1                                             ) */
/* >           (      ...                                      ) */
/* >           (             1                                 ) */
/* >           (                  c(k)                    s(k) ) */
/* >           (                         1                     ) */
/* >           (                              ...              ) */
/* >           (                                     1         ) */
/* >           (                 -s(k)                    c(k) ) */
/* > */
/* > where R(k) appears in rows and columns k and z.  The rotations are */
/* > performed without ever forming P(k) explicitly. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          Specifies whether the plane rotation matrix P is applied to */
/* >          A on the left or the right. */
/* >          = 'L':  Left, compute A := P*A */
/* >          = 'R':  Right, compute A:= A*P**T */
/* > \endverbatim */
/* > */
/* > \param[in] PIVOT */
/* > \verbatim */
/* >          PIVOT is CHARACTER*1 */
/* >          Specifies the plane for which P(k) is a plane rotation */
/* >          matrix. */
/* >          = 'V':  Variable pivot, the plane (k,k+1) */
/* >          = 'T':  Top pivot, the plane (1,k+1) */
/* >          = 'B':  Bottom pivot, the plane (k,z) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* >          DIRECT is CHARACTER*1 */
/* >          Specifies whether P is a forward or backward sequence of */
/* >          plane rotations. */
/* >          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1) */
/* >          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1) */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  If m <= 1, an immediate */
/* >          return is effected. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  If n <= 1, an */
/* >          immediate return is effected. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension */
/* >                  (M-1) if SIDE = 'L' */
/* >                  (N-1) if SIDE = 'R' */
/* >          The cosines c(k) of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension */
/* >                  (M-1) if SIDE = 'L' */
/* >                  (N-1) if SIDE = 'R' */
/* >          The sines s(k) of the plane rotations.  The 2-by-2 plane */
/* >          rotation part of the matrix P(k), R(k), has the form */
/* >          R(k) = (  c(k)  s(k) ) */
/* >                 ( -s(k)  c(k) ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The M-by-N matrix A.  On exit, A is overwritten by P*A if */
/* >          SIDE = 'R' or by A*P**T if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *
	lda, ftnlen side_len, ftnlen pivot_len, ftnlen direct_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, info;
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal ctemp, stemp;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 239 "dlasr.f"
    /* Parameter adjustments */
#line 239 "dlasr.f"
    --c__;
#line 239 "dlasr.f"
    --s;
#line 239 "dlasr.f"
    a_dim1 = *lda;
#line 239 "dlasr.f"
    a_offset = 1 + a_dim1;
#line 239 "dlasr.f"
    a -= a_offset;
#line 239 "dlasr.f"

#line 239 "dlasr.f"
    /* Function Body */
#line 239 "dlasr.f"
    info = 0;
#line 240 "dlasr.f"
    if (! (lsame_(side, "L", (ftnlen)1, (ftnlen)1) || lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1))) {
#line 241 "dlasr.f"
	info = 1;
#line 242 "dlasr.f"
    } else if (! (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1) || lsame_(pivot, 
	    "T", (ftnlen)1, (ftnlen)1) || lsame_(pivot, "B", (ftnlen)1, (
	    ftnlen)1))) {
#line 244 "dlasr.f"
	info = 2;
#line 245 "dlasr.f"
    } else if (! (lsame_(direct, "F", (ftnlen)1, (ftnlen)1) || lsame_(direct, 
	    "B", (ftnlen)1, (ftnlen)1))) {
#line 247 "dlasr.f"
	info = 3;
#line 248 "dlasr.f"
    } else if (*m < 0) {
#line 249 "dlasr.f"
	info = 4;
#line 250 "dlasr.f"
    } else if (*n < 0) {
#line 251 "dlasr.f"
	info = 5;
#line 252 "dlasr.f"
    } else if (*lda < max(1,*m)) {
#line 253 "dlasr.f"
	info = 9;
#line 254 "dlasr.f"
    }
#line 255 "dlasr.f"
    if (info != 0) {
#line 256 "dlasr.f"
	xerbla_("DLASR ", &info, (ftnlen)6);
#line 257 "dlasr.f"
	return 0;
#line 258 "dlasr.f"
    }

/*     Quick return if possible */

#line 262 "dlasr.f"
    if (*m == 0 || *n == 0) {
#line 262 "dlasr.f"
	return 0;
#line 262 "dlasr.f"
    }
#line 264 "dlasr.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  P * A */

#line 268 "dlasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 269 "dlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 270 "dlasr.f"
		i__1 = *m - 1;
#line 270 "dlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 271 "dlasr.f"
		    ctemp = c__[j];
#line 272 "dlasr.f"
		    stemp = s[j];
#line 273 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 274 "dlasr.f"
			i__2 = *n;
#line 274 "dlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 275 "dlasr.f"
			    temp = a[j + 1 + i__ * a_dim1];
#line 276 "dlasr.f"
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
#line 277 "dlasr.f"
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
#line 278 "dlasr.f"
/* L10: */
#line 278 "dlasr.f"
			}
#line 279 "dlasr.f"
		    }
#line 280 "dlasr.f"
/* L20: */
#line 280 "dlasr.f"
		}
#line 281 "dlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 282 "dlasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 283 "dlasr.f"
		    ctemp = c__[j];
#line 284 "dlasr.f"
		    stemp = s[j];
#line 285 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 286 "dlasr.f"
			i__1 = *n;
#line 286 "dlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "dlasr.f"
			    temp = a[j + 1 + i__ * a_dim1];
#line 288 "dlasr.f"
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
#line 289 "dlasr.f"
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
#line 290 "dlasr.f"
/* L30: */
#line 290 "dlasr.f"
			}
#line 291 "dlasr.f"
		    }
#line 292 "dlasr.f"
/* L40: */
#line 292 "dlasr.f"
		}
#line 293 "dlasr.f"
	    }
#line 294 "dlasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 295 "dlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 296 "dlasr.f"
		i__1 = *m;
#line 296 "dlasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 297 "dlasr.f"
		    ctemp = c__[j - 1];
#line 298 "dlasr.f"
		    stemp = s[j - 1];
#line 299 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 300 "dlasr.f"
			i__2 = *n;
#line 300 "dlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 301 "dlasr.f"
			    temp = a[j + i__ * a_dim1];
#line 302 "dlasr.f"
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
#line 303 "dlasr.f"
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
#line 304 "dlasr.f"
/* L50: */
#line 304 "dlasr.f"
			}
#line 305 "dlasr.f"
		    }
#line 306 "dlasr.f"
/* L60: */
#line 306 "dlasr.f"
		}
#line 307 "dlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 308 "dlasr.f"
		for (j = *m; j >= 2; --j) {
#line 309 "dlasr.f"
		    ctemp = c__[j - 1];
#line 310 "dlasr.f"
		    stemp = s[j - 1];
#line 311 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 312 "dlasr.f"
			i__1 = *n;
#line 312 "dlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "dlasr.f"
			    temp = a[j + i__ * a_dim1];
#line 314 "dlasr.f"
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
#line 315 "dlasr.f"
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
#line 316 "dlasr.f"
/* L70: */
#line 316 "dlasr.f"
			}
#line 317 "dlasr.f"
		    }
#line 318 "dlasr.f"
/* L80: */
#line 318 "dlasr.f"
		}
#line 319 "dlasr.f"
	    }
#line 320 "dlasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 321 "dlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 322 "dlasr.f"
		i__1 = *m - 1;
#line 322 "dlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 323 "dlasr.f"
		    ctemp = c__[j];
#line 324 "dlasr.f"
		    stemp = s[j];
#line 325 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 326 "dlasr.f"
			i__2 = *n;
#line 326 "dlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 327 "dlasr.f"
			    temp = a[j + i__ * a_dim1];
#line 328 "dlasr.f"
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
#line 329 "dlasr.f"
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
#line 330 "dlasr.f"
/* L90: */
#line 330 "dlasr.f"
			}
#line 331 "dlasr.f"
		    }
#line 332 "dlasr.f"
/* L100: */
#line 332 "dlasr.f"
		}
#line 333 "dlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 334 "dlasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 335 "dlasr.f"
		    ctemp = c__[j];
#line 336 "dlasr.f"
		    stemp = s[j];
#line 337 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 338 "dlasr.f"
			i__1 = *n;
#line 338 "dlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 339 "dlasr.f"
			    temp = a[j + i__ * a_dim1];
#line 340 "dlasr.f"
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
#line 341 "dlasr.f"
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
#line 342 "dlasr.f"
/* L110: */
#line 342 "dlasr.f"
			}
#line 343 "dlasr.f"
		    }
#line 344 "dlasr.f"
/* L120: */
#line 344 "dlasr.f"
		}
#line 345 "dlasr.f"
	    }
#line 346 "dlasr.f"
	}
#line 347 "dlasr.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form A * P**T */

#line 351 "dlasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 352 "dlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 353 "dlasr.f"
		i__1 = *n - 1;
#line 353 "dlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 354 "dlasr.f"
		    ctemp = c__[j];
#line 355 "dlasr.f"
		    stemp = s[j];
#line 356 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 357 "dlasr.f"
			i__2 = *m;
#line 357 "dlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 358 "dlasr.f"
			    temp = a[i__ + (j + 1) * a_dim1];
#line 359 "dlasr.f"
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
#line 360 "dlasr.f"
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
#line 361 "dlasr.f"
/* L130: */
#line 361 "dlasr.f"
			}
#line 362 "dlasr.f"
		    }
#line 363 "dlasr.f"
/* L140: */
#line 363 "dlasr.f"
		}
#line 364 "dlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 365 "dlasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 366 "dlasr.f"
		    ctemp = c__[j];
#line 367 "dlasr.f"
		    stemp = s[j];
#line 368 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 369 "dlasr.f"
			i__1 = *m;
#line 369 "dlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 370 "dlasr.f"
			    temp = a[i__ + (j + 1) * a_dim1];
#line 371 "dlasr.f"
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
#line 372 "dlasr.f"
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
#line 373 "dlasr.f"
/* L150: */
#line 373 "dlasr.f"
			}
#line 374 "dlasr.f"
		    }
#line 375 "dlasr.f"
/* L160: */
#line 375 "dlasr.f"
		}
#line 376 "dlasr.f"
	    }
#line 377 "dlasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 378 "dlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 379 "dlasr.f"
		i__1 = *n;
#line 379 "dlasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 380 "dlasr.f"
		    ctemp = c__[j - 1];
#line 381 "dlasr.f"
		    stemp = s[j - 1];
#line 382 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 383 "dlasr.f"
			i__2 = *m;
#line 383 "dlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 384 "dlasr.f"
			    temp = a[i__ + j * a_dim1];
#line 385 "dlasr.f"
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
#line 386 "dlasr.f"
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
#line 387 "dlasr.f"
/* L170: */
#line 387 "dlasr.f"
			}
#line 388 "dlasr.f"
		    }
#line 389 "dlasr.f"
/* L180: */
#line 389 "dlasr.f"
		}
#line 390 "dlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 391 "dlasr.f"
		for (j = *n; j >= 2; --j) {
#line 392 "dlasr.f"
		    ctemp = c__[j - 1];
#line 393 "dlasr.f"
		    stemp = s[j - 1];
#line 394 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 395 "dlasr.f"
			i__1 = *m;
#line 395 "dlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 396 "dlasr.f"
			    temp = a[i__ + j * a_dim1];
#line 397 "dlasr.f"
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
#line 398 "dlasr.f"
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
#line 399 "dlasr.f"
/* L190: */
#line 399 "dlasr.f"
			}
#line 400 "dlasr.f"
		    }
#line 401 "dlasr.f"
/* L200: */
#line 401 "dlasr.f"
		}
#line 402 "dlasr.f"
	    }
#line 403 "dlasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 404 "dlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 405 "dlasr.f"
		i__1 = *n - 1;
#line 405 "dlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 406 "dlasr.f"
		    ctemp = c__[j];
#line 407 "dlasr.f"
		    stemp = s[j];
#line 408 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 409 "dlasr.f"
			i__2 = *m;
#line 409 "dlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 410 "dlasr.f"
			    temp = a[i__ + j * a_dim1];
#line 411 "dlasr.f"
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
#line 412 "dlasr.f"
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
#line 413 "dlasr.f"
/* L210: */
#line 413 "dlasr.f"
			}
#line 414 "dlasr.f"
		    }
#line 415 "dlasr.f"
/* L220: */
#line 415 "dlasr.f"
		}
#line 416 "dlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 417 "dlasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 418 "dlasr.f"
		    ctemp = c__[j];
#line 419 "dlasr.f"
		    stemp = s[j];
#line 420 "dlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 421 "dlasr.f"
			i__1 = *m;
#line 421 "dlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 422 "dlasr.f"
			    temp = a[i__ + j * a_dim1];
#line 423 "dlasr.f"
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
#line 424 "dlasr.f"
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
#line 425 "dlasr.f"
/* L230: */
#line 425 "dlasr.f"
			}
#line 426 "dlasr.f"
		    }
#line 427 "dlasr.f"
/* L240: */
#line 427 "dlasr.f"
		}
#line 428 "dlasr.f"
	    }
#line 429 "dlasr.f"
	}
#line 430 "dlasr.f"
    }

#line 432 "dlasr.f"
    return 0;

/*     End of DLASR */

} /* dlasr_ */

