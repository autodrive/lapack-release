#line 1 "slasr.f"
/* slasr.f -- translated by f2c (version 20100827).
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

#line 1 "slasr.f"
/* > \brief \b SLASR applies a sequence of plane rotations to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasr.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasr.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasr.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, PIVOT, SIDE */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASR applies a sequence of plane rotations to a real matrix A, */
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
/* >          C is REAL array, dimension */
/* >                  (M-1) if SIDE = 'L' */
/* >                  (N-1) if SIDE = 'R' */
/* >          The cosines c(k) of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL array, dimension */
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
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slasr_(char *side, char *pivot, char *direct, integer *m,
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


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
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
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 239 "slasr.f"
    /* Parameter adjustments */
#line 239 "slasr.f"
    --c__;
#line 239 "slasr.f"
    --s;
#line 239 "slasr.f"
    a_dim1 = *lda;
#line 239 "slasr.f"
    a_offset = 1 + a_dim1;
#line 239 "slasr.f"
    a -= a_offset;
#line 239 "slasr.f"

#line 239 "slasr.f"
    /* Function Body */
#line 239 "slasr.f"
    info = 0;
#line 240 "slasr.f"
    if (! (lsame_(side, "L", (ftnlen)1, (ftnlen)1) || lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1))) {
#line 241 "slasr.f"
	info = 1;
#line 242 "slasr.f"
    } else if (! (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1) || lsame_(pivot, 
	    "T", (ftnlen)1, (ftnlen)1) || lsame_(pivot, "B", (ftnlen)1, (
	    ftnlen)1))) {
#line 244 "slasr.f"
	info = 2;
#line 245 "slasr.f"
    } else if (! (lsame_(direct, "F", (ftnlen)1, (ftnlen)1) || lsame_(direct, 
	    "B", (ftnlen)1, (ftnlen)1))) {
#line 247 "slasr.f"
	info = 3;
#line 248 "slasr.f"
    } else if (*m < 0) {
#line 249 "slasr.f"
	info = 4;
#line 250 "slasr.f"
    } else if (*n < 0) {
#line 251 "slasr.f"
	info = 5;
#line 252 "slasr.f"
    } else if (*lda < max(1,*m)) {
#line 253 "slasr.f"
	info = 9;
#line 254 "slasr.f"
    }
#line 255 "slasr.f"
    if (info != 0) {
#line 256 "slasr.f"
	xerbla_("SLASR ", &info, (ftnlen)6);
#line 257 "slasr.f"
	return 0;
#line 258 "slasr.f"
    }

/*     Quick return if possible */

#line 262 "slasr.f"
    if (*m == 0 || *n == 0) {
#line 262 "slasr.f"
	return 0;
#line 262 "slasr.f"
    }
#line 264 "slasr.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  P * A */

#line 268 "slasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 269 "slasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 270 "slasr.f"
		i__1 = *m - 1;
#line 270 "slasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 271 "slasr.f"
		    ctemp = c__[j];
#line 272 "slasr.f"
		    stemp = s[j];
#line 273 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 274 "slasr.f"
			i__2 = *n;
#line 274 "slasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 275 "slasr.f"
			    temp = a[j + 1 + i__ * a_dim1];
#line 276 "slasr.f"
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
#line 277 "slasr.f"
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
#line 278 "slasr.f"
/* L10: */
#line 278 "slasr.f"
			}
#line 279 "slasr.f"
		    }
#line 280 "slasr.f"
/* L20: */
#line 280 "slasr.f"
		}
#line 281 "slasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 282 "slasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 283 "slasr.f"
		    ctemp = c__[j];
#line 284 "slasr.f"
		    stemp = s[j];
#line 285 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 286 "slasr.f"
			i__1 = *n;
#line 286 "slasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "slasr.f"
			    temp = a[j + 1 + i__ * a_dim1];
#line 288 "slasr.f"
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
#line 289 "slasr.f"
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
#line 290 "slasr.f"
/* L30: */
#line 290 "slasr.f"
			}
#line 291 "slasr.f"
		    }
#line 292 "slasr.f"
/* L40: */
#line 292 "slasr.f"
		}
#line 293 "slasr.f"
	    }
#line 294 "slasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 295 "slasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 296 "slasr.f"
		i__1 = *m;
#line 296 "slasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 297 "slasr.f"
		    ctemp = c__[j - 1];
#line 298 "slasr.f"
		    stemp = s[j - 1];
#line 299 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 300 "slasr.f"
			i__2 = *n;
#line 300 "slasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 301 "slasr.f"
			    temp = a[j + i__ * a_dim1];
#line 302 "slasr.f"
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
#line 303 "slasr.f"
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
#line 304 "slasr.f"
/* L50: */
#line 304 "slasr.f"
			}
#line 305 "slasr.f"
		    }
#line 306 "slasr.f"
/* L60: */
#line 306 "slasr.f"
		}
#line 307 "slasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 308 "slasr.f"
		for (j = *m; j >= 2; --j) {
#line 309 "slasr.f"
		    ctemp = c__[j - 1];
#line 310 "slasr.f"
		    stemp = s[j - 1];
#line 311 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 312 "slasr.f"
			i__1 = *n;
#line 312 "slasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "slasr.f"
			    temp = a[j + i__ * a_dim1];
#line 314 "slasr.f"
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
#line 315 "slasr.f"
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
#line 316 "slasr.f"
/* L70: */
#line 316 "slasr.f"
			}
#line 317 "slasr.f"
		    }
#line 318 "slasr.f"
/* L80: */
#line 318 "slasr.f"
		}
#line 319 "slasr.f"
	    }
#line 320 "slasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 321 "slasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 322 "slasr.f"
		i__1 = *m - 1;
#line 322 "slasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 323 "slasr.f"
		    ctemp = c__[j];
#line 324 "slasr.f"
		    stemp = s[j];
#line 325 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 326 "slasr.f"
			i__2 = *n;
#line 326 "slasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 327 "slasr.f"
			    temp = a[j + i__ * a_dim1];
#line 328 "slasr.f"
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
#line 329 "slasr.f"
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
#line 330 "slasr.f"
/* L90: */
#line 330 "slasr.f"
			}
#line 331 "slasr.f"
		    }
#line 332 "slasr.f"
/* L100: */
#line 332 "slasr.f"
		}
#line 333 "slasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 334 "slasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 335 "slasr.f"
		    ctemp = c__[j];
#line 336 "slasr.f"
		    stemp = s[j];
#line 337 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 338 "slasr.f"
			i__1 = *n;
#line 338 "slasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 339 "slasr.f"
			    temp = a[j + i__ * a_dim1];
#line 340 "slasr.f"
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
#line 341 "slasr.f"
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
#line 342 "slasr.f"
/* L110: */
#line 342 "slasr.f"
			}
#line 343 "slasr.f"
		    }
#line 344 "slasr.f"
/* L120: */
#line 344 "slasr.f"
		}
#line 345 "slasr.f"
	    }
#line 346 "slasr.f"
	}
#line 347 "slasr.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form A * P**T */

#line 351 "slasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 352 "slasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 353 "slasr.f"
		i__1 = *n - 1;
#line 353 "slasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 354 "slasr.f"
		    ctemp = c__[j];
#line 355 "slasr.f"
		    stemp = s[j];
#line 356 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 357 "slasr.f"
			i__2 = *m;
#line 357 "slasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 358 "slasr.f"
			    temp = a[i__ + (j + 1) * a_dim1];
#line 359 "slasr.f"
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
#line 360 "slasr.f"
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
#line 361 "slasr.f"
/* L130: */
#line 361 "slasr.f"
			}
#line 362 "slasr.f"
		    }
#line 363 "slasr.f"
/* L140: */
#line 363 "slasr.f"
		}
#line 364 "slasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 365 "slasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 366 "slasr.f"
		    ctemp = c__[j];
#line 367 "slasr.f"
		    stemp = s[j];
#line 368 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 369 "slasr.f"
			i__1 = *m;
#line 369 "slasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 370 "slasr.f"
			    temp = a[i__ + (j + 1) * a_dim1];
#line 371 "slasr.f"
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
#line 372 "slasr.f"
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
#line 373 "slasr.f"
/* L150: */
#line 373 "slasr.f"
			}
#line 374 "slasr.f"
		    }
#line 375 "slasr.f"
/* L160: */
#line 375 "slasr.f"
		}
#line 376 "slasr.f"
	    }
#line 377 "slasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 378 "slasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 379 "slasr.f"
		i__1 = *n;
#line 379 "slasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 380 "slasr.f"
		    ctemp = c__[j - 1];
#line 381 "slasr.f"
		    stemp = s[j - 1];
#line 382 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 383 "slasr.f"
			i__2 = *m;
#line 383 "slasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 384 "slasr.f"
			    temp = a[i__ + j * a_dim1];
#line 385 "slasr.f"
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
#line 386 "slasr.f"
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
#line 387 "slasr.f"
/* L170: */
#line 387 "slasr.f"
			}
#line 388 "slasr.f"
		    }
#line 389 "slasr.f"
/* L180: */
#line 389 "slasr.f"
		}
#line 390 "slasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 391 "slasr.f"
		for (j = *n; j >= 2; --j) {
#line 392 "slasr.f"
		    ctemp = c__[j - 1];
#line 393 "slasr.f"
		    stemp = s[j - 1];
#line 394 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 395 "slasr.f"
			i__1 = *m;
#line 395 "slasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 396 "slasr.f"
			    temp = a[i__ + j * a_dim1];
#line 397 "slasr.f"
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
#line 398 "slasr.f"
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
#line 399 "slasr.f"
/* L190: */
#line 399 "slasr.f"
			}
#line 400 "slasr.f"
		    }
#line 401 "slasr.f"
/* L200: */
#line 401 "slasr.f"
		}
#line 402 "slasr.f"
	    }
#line 403 "slasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 404 "slasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 405 "slasr.f"
		i__1 = *n - 1;
#line 405 "slasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 406 "slasr.f"
		    ctemp = c__[j];
#line 407 "slasr.f"
		    stemp = s[j];
#line 408 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 409 "slasr.f"
			i__2 = *m;
#line 409 "slasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 410 "slasr.f"
			    temp = a[i__ + j * a_dim1];
#line 411 "slasr.f"
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
#line 412 "slasr.f"
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
#line 413 "slasr.f"
/* L210: */
#line 413 "slasr.f"
			}
#line 414 "slasr.f"
		    }
#line 415 "slasr.f"
/* L220: */
#line 415 "slasr.f"
		}
#line 416 "slasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 417 "slasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 418 "slasr.f"
		    ctemp = c__[j];
#line 419 "slasr.f"
		    stemp = s[j];
#line 420 "slasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 421 "slasr.f"
			i__1 = *m;
#line 421 "slasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 422 "slasr.f"
			    temp = a[i__ + j * a_dim1];
#line 423 "slasr.f"
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
#line 424 "slasr.f"
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
#line 425 "slasr.f"
/* L230: */
#line 425 "slasr.f"
			}
#line 426 "slasr.f"
		    }
#line 427 "slasr.f"
/* L240: */
#line 427 "slasr.f"
		}
#line 428 "slasr.f"
	    }
#line 429 "slasr.f"
	}
#line 430 "slasr.f"
    }

#line 432 "slasr.f"
    return 0;

/*     End of SLASR */

} /* slasr_ */

