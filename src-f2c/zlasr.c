#line 1 "zlasr.f"
/* zlasr.f -- translated by f2c (version 20100827).
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

#line 1 "zlasr.f"
/* > \brief \b ZLASR applies a sequence of plane rotations to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLASR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlasr.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlasr.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlasr.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, PIVOT, SIDE */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLASR applies a sequence of real plane rotations to a complex matrix */
/* > A, from either the left or the right. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c__, doublereal *s, doublecomplex *a, 
	integer *lda, ftnlen side_len, ftnlen pivot_len, ftnlen direct_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, j, info;
    static doublecomplex temp;
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 242 "zlasr.f"
    /* Parameter adjustments */
#line 242 "zlasr.f"
    --c__;
#line 242 "zlasr.f"
    --s;
#line 242 "zlasr.f"
    a_dim1 = *lda;
#line 242 "zlasr.f"
    a_offset = 1 + a_dim1;
#line 242 "zlasr.f"
    a -= a_offset;
#line 242 "zlasr.f"

#line 242 "zlasr.f"
    /* Function Body */
#line 242 "zlasr.f"
    info = 0;
#line 243 "zlasr.f"
    if (! (lsame_(side, "L", (ftnlen)1, (ftnlen)1) || lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1))) {
#line 244 "zlasr.f"
	info = 1;
#line 245 "zlasr.f"
    } else if (! (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1) || lsame_(pivot, 
	    "T", (ftnlen)1, (ftnlen)1) || lsame_(pivot, "B", (ftnlen)1, (
	    ftnlen)1))) {
#line 247 "zlasr.f"
	info = 2;
#line 248 "zlasr.f"
    } else if (! (lsame_(direct, "F", (ftnlen)1, (ftnlen)1) || lsame_(direct, 
	    "B", (ftnlen)1, (ftnlen)1))) {
#line 250 "zlasr.f"
	info = 3;
#line 251 "zlasr.f"
    } else if (*m < 0) {
#line 252 "zlasr.f"
	info = 4;
#line 253 "zlasr.f"
    } else if (*n < 0) {
#line 254 "zlasr.f"
	info = 5;
#line 255 "zlasr.f"
    } else if (*lda < max(1,*m)) {
#line 256 "zlasr.f"
	info = 9;
#line 257 "zlasr.f"
    }
#line 258 "zlasr.f"
    if (info != 0) {
#line 259 "zlasr.f"
	xerbla_("ZLASR ", &info, (ftnlen)6);
#line 260 "zlasr.f"
	return 0;
#line 261 "zlasr.f"
    }

/*     Quick return if possible */

#line 265 "zlasr.f"
    if (*m == 0 || *n == 0) {
#line 265 "zlasr.f"
	return 0;
#line 265 "zlasr.f"
    }
#line 267 "zlasr.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  P * A */

#line 271 "zlasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 272 "zlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 273 "zlasr.f"
		i__1 = *m - 1;
#line 273 "zlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "zlasr.f"
		    ctemp = c__[j];
#line 275 "zlasr.f"
		    stemp = s[j];
#line 276 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 277 "zlasr.f"
			i__2 = *n;
#line 277 "zlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "zlasr.f"
			    i__3 = j + 1 + i__ * a_dim1;
#line 278 "zlasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 279 "zlasr.f"
			    i__3 = j + 1 + i__ * a_dim1;
#line 279 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 279 "zlasr.f"
			    i__4 = j + i__ * a_dim1;
#line 279 "zlasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 279 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 279 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 280 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 280 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 280 "zlasr.f"
			    i__4 = j + i__ * a_dim1;
#line 280 "zlasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 280 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 280 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 281 "zlasr.f"
/* L10: */
#line 281 "zlasr.f"
			}
#line 282 "zlasr.f"
		    }
#line 283 "zlasr.f"
/* L20: */
#line 283 "zlasr.f"
		}
#line 284 "zlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 285 "zlasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 286 "zlasr.f"
		    ctemp = c__[j];
#line 287 "zlasr.f"
		    stemp = s[j];
#line 288 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 289 "zlasr.f"
			i__1 = *n;
#line 289 "zlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "zlasr.f"
			    i__2 = j + 1 + i__ * a_dim1;
#line 290 "zlasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 291 "zlasr.f"
			    i__2 = j + 1 + i__ * a_dim1;
#line 291 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 291 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 291 "zlasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 291 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 291 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 292 "zlasr.f"
			    i__2 = j + i__ * a_dim1;
#line 292 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 292 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 292 "zlasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 292 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 292 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 293 "zlasr.f"
/* L30: */
#line 293 "zlasr.f"
			}
#line 294 "zlasr.f"
		    }
#line 295 "zlasr.f"
/* L40: */
#line 295 "zlasr.f"
		}
#line 296 "zlasr.f"
	    }
#line 297 "zlasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 298 "zlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 299 "zlasr.f"
		i__1 = *m;
#line 299 "zlasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 300 "zlasr.f"
		    ctemp = c__[j - 1];
#line 301 "zlasr.f"
		    stemp = s[j - 1];
#line 302 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 303 "zlasr.f"
			i__2 = *n;
#line 303 "zlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 304 "zlasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 305 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 305 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 305 "zlasr.f"
			    i__4 = i__ * a_dim1 + 1;
#line 305 "zlasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 305 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 305 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 306 "zlasr.f"
			    i__3 = i__ * a_dim1 + 1;
#line 306 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 306 "zlasr.f"
			    i__4 = i__ * a_dim1 + 1;
#line 306 "zlasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 306 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 306 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 307 "zlasr.f"
/* L50: */
#line 307 "zlasr.f"
			}
#line 308 "zlasr.f"
		    }
#line 309 "zlasr.f"
/* L60: */
#line 309 "zlasr.f"
		}
#line 310 "zlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 311 "zlasr.f"
		for (j = *m; j >= 2; --j) {
#line 312 "zlasr.f"
		    ctemp = c__[j - 1];
#line 313 "zlasr.f"
		    stemp = s[j - 1];
#line 314 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 315 "zlasr.f"
			i__1 = *n;
#line 315 "zlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "zlasr.f"
			    i__2 = j + i__ * a_dim1;
#line 316 "zlasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 317 "zlasr.f"
			    i__2 = j + i__ * a_dim1;
#line 317 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 317 "zlasr.f"
			    i__3 = i__ * a_dim1 + 1;
#line 317 "zlasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 317 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 317 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 318 "zlasr.f"
			    i__2 = i__ * a_dim1 + 1;
#line 318 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 318 "zlasr.f"
			    i__3 = i__ * a_dim1 + 1;
#line 318 "zlasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 318 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 318 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 319 "zlasr.f"
/* L70: */
#line 319 "zlasr.f"
			}
#line 320 "zlasr.f"
		    }
#line 321 "zlasr.f"
/* L80: */
#line 321 "zlasr.f"
		}
#line 322 "zlasr.f"
	    }
#line 323 "zlasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 324 "zlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 325 "zlasr.f"
		i__1 = *m - 1;
#line 325 "zlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "zlasr.f"
		    ctemp = c__[j];
#line 327 "zlasr.f"
		    stemp = s[j];
#line 328 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 329 "zlasr.f"
			i__2 = *n;
#line 329 "zlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 330 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 330 "zlasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 331 "zlasr.f"
			    i__3 = j + i__ * a_dim1;
#line 331 "zlasr.f"
			    i__4 = *m + i__ * a_dim1;
#line 331 "zlasr.f"
			    z__2.r = stemp * a[i__4].r, z__2.i = stemp * a[
				    i__4].i;
#line 331 "zlasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 331 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 331 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 332 "zlasr.f"
			    i__3 = *m + i__ * a_dim1;
#line 332 "zlasr.f"
			    i__4 = *m + i__ * a_dim1;
#line 332 "zlasr.f"
			    z__2.r = ctemp * a[i__4].r, z__2.i = ctemp * a[
				    i__4].i;
#line 332 "zlasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 332 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 332 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 333 "zlasr.f"
/* L90: */
#line 333 "zlasr.f"
			}
#line 334 "zlasr.f"
		    }
#line 335 "zlasr.f"
/* L100: */
#line 335 "zlasr.f"
		}
#line 336 "zlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 337 "zlasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 338 "zlasr.f"
		    ctemp = c__[j];
#line 339 "zlasr.f"
		    stemp = s[j];
#line 340 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 341 "zlasr.f"
			i__1 = *n;
#line 341 "zlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 342 "zlasr.f"
			    i__2 = j + i__ * a_dim1;
#line 342 "zlasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 343 "zlasr.f"
			    i__2 = j + i__ * a_dim1;
#line 343 "zlasr.f"
			    i__3 = *m + i__ * a_dim1;
#line 343 "zlasr.f"
			    z__2.r = stemp * a[i__3].r, z__2.i = stemp * a[
				    i__3].i;
#line 343 "zlasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 343 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 343 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 344 "zlasr.f"
			    i__2 = *m + i__ * a_dim1;
#line 344 "zlasr.f"
			    i__3 = *m + i__ * a_dim1;
#line 344 "zlasr.f"
			    z__2.r = ctemp * a[i__3].r, z__2.i = ctemp * a[
				    i__3].i;
#line 344 "zlasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 344 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 344 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 345 "zlasr.f"
/* L110: */
#line 345 "zlasr.f"
			}
#line 346 "zlasr.f"
		    }
#line 347 "zlasr.f"
/* L120: */
#line 347 "zlasr.f"
		}
#line 348 "zlasr.f"
	    }
#line 349 "zlasr.f"
	}
#line 350 "zlasr.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form A * P**T */

#line 354 "zlasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 355 "zlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 356 "zlasr.f"
		i__1 = *n - 1;
#line 356 "zlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 357 "zlasr.f"
		    ctemp = c__[j];
#line 358 "zlasr.f"
		    stemp = s[j];
#line 359 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 360 "zlasr.f"
			i__2 = *m;
#line 360 "zlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 361 "zlasr.f"
			    i__3 = i__ + (j + 1) * a_dim1;
#line 361 "zlasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 362 "zlasr.f"
			    i__3 = i__ + (j + 1) * a_dim1;
#line 362 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 362 "zlasr.f"
			    i__4 = i__ + j * a_dim1;
#line 362 "zlasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 362 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 362 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 363 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 363 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 363 "zlasr.f"
			    i__4 = i__ + j * a_dim1;
#line 363 "zlasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 363 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 363 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 364 "zlasr.f"
/* L130: */
#line 364 "zlasr.f"
			}
#line 365 "zlasr.f"
		    }
#line 366 "zlasr.f"
/* L140: */
#line 366 "zlasr.f"
		}
#line 367 "zlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 368 "zlasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 369 "zlasr.f"
		    ctemp = c__[j];
#line 370 "zlasr.f"
		    stemp = s[j];
#line 371 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 372 "zlasr.f"
			i__1 = *m;
#line 372 "zlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 373 "zlasr.f"
			    i__2 = i__ + (j + 1) * a_dim1;
#line 373 "zlasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 374 "zlasr.f"
			    i__2 = i__ + (j + 1) * a_dim1;
#line 374 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 374 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 374 "zlasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 374 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 374 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 375 "zlasr.f"
			    i__2 = i__ + j * a_dim1;
#line 375 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 375 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 375 "zlasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 375 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 375 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 376 "zlasr.f"
/* L150: */
#line 376 "zlasr.f"
			}
#line 377 "zlasr.f"
		    }
#line 378 "zlasr.f"
/* L160: */
#line 378 "zlasr.f"
		}
#line 379 "zlasr.f"
	    }
#line 380 "zlasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 381 "zlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 382 "zlasr.f"
		i__1 = *n;
#line 382 "zlasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 383 "zlasr.f"
		    ctemp = c__[j - 1];
#line 384 "zlasr.f"
		    stemp = s[j - 1];
#line 385 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 386 "zlasr.f"
			i__2 = *m;
#line 386 "zlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 387 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 387 "zlasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 388 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 388 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 388 "zlasr.f"
			    i__4 = i__ + a_dim1;
#line 388 "zlasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 388 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 388 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 389 "zlasr.f"
			    i__3 = i__ + a_dim1;
#line 389 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 389 "zlasr.f"
			    i__4 = i__ + a_dim1;
#line 389 "zlasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 389 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 389 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 390 "zlasr.f"
/* L170: */
#line 390 "zlasr.f"
			}
#line 391 "zlasr.f"
		    }
#line 392 "zlasr.f"
/* L180: */
#line 392 "zlasr.f"
		}
#line 393 "zlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 394 "zlasr.f"
		for (j = *n; j >= 2; --j) {
#line 395 "zlasr.f"
		    ctemp = c__[j - 1];
#line 396 "zlasr.f"
		    stemp = s[j - 1];
#line 397 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 398 "zlasr.f"
			i__1 = *m;
#line 398 "zlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 399 "zlasr.f"
			    i__2 = i__ + j * a_dim1;
#line 399 "zlasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 400 "zlasr.f"
			    i__2 = i__ + j * a_dim1;
#line 400 "zlasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 400 "zlasr.f"
			    i__3 = i__ + a_dim1;
#line 400 "zlasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 400 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 400 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 401 "zlasr.f"
			    i__2 = i__ + a_dim1;
#line 401 "zlasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 401 "zlasr.f"
			    i__3 = i__ + a_dim1;
#line 401 "zlasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 401 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 401 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 402 "zlasr.f"
/* L190: */
#line 402 "zlasr.f"
			}
#line 403 "zlasr.f"
		    }
#line 404 "zlasr.f"
/* L200: */
#line 404 "zlasr.f"
		}
#line 405 "zlasr.f"
	    }
#line 406 "zlasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 407 "zlasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 408 "zlasr.f"
		i__1 = *n - 1;
#line 408 "zlasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 409 "zlasr.f"
		    ctemp = c__[j];
#line 410 "zlasr.f"
		    stemp = s[j];
#line 411 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 412 "zlasr.f"
			i__2 = *m;
#line 412 "zlasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 413 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 413 "zlasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 414 "zlasr.f"
			    i__3 = i__ + j * a_dim1;
#line 414 "zlasr.f"
			    i__4 = i__ + *n * a_dim1;
#line 414 "zlasr.f"
			    z__2.r = stemp * a[i__4].r, z__2.i = stemp * a[
				    i__4].i;
#line 414 "zlasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 414 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 414 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 415 "zlasr.f"
			    i__3 = i__ + *n * a_dim1;
#line 415 "zlasr.f"
			    i__4 = i__ + *n * a_dim1;
#line 415 "zlasr.f"
			    z__2.r = ctemp * a[i__4].r, z__2.i = ctemp * a[
				    i__4].i;
#line 415 "zlasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 415 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 415 "zlasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 416 "zlasr.f"
/* L210: */
#line 416 "zlasr.f"
			}
#line 417 "zlasr.f"
		    }
#line 418 "zlasr.f"
/* L220: */
#line 418 "zlasr.f"
		}
#line 419 "zlasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 420 "zlasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 421 "zlasr.f"
		    ctemp = c__[j];
#line 422 "zlasr.f"
		    stemp = s[j];
#line 423 "zlasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 424 "zlasr.f"
			i__1 = *m;
#line 424 "zlasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 425 "zlasr.f"
			    i__2 = i__ + j * a_dim1;
#line 425 "zlasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 426 "zlasr.f"
			    i__2 = i__ + j * a_dim1;
#line 426 "zlasr.f"
			    i__3 = i__ + *n * a_dim1;
#line 426 "zlasr.f"
			    z__2.r = stemp * a[i__3].r, z__2.i = stemp * a[
				    i__3].i;
#line 426 "zlasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 426 "zlasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 426 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 427 "zlasr.f"
			    i__2 = i__ + *n * a_dim1;
#line 427 "zlasr.f"
			    i__3 = i__ + *n * a_dim1;
#line 427 "zlasr.f"
			    z__2.r = ctemp * a[i__3].r, z__2.i = ctemp * a[
				    i__3].i;
#line 427 "zlasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 427 "zlasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 427 "zlasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 428 "zlasr.f"
/* L230: */
#line 428 "zlasr.f"
			}
#line 429 "zlasr.f"
		    }
#line 430 "zlasr.f"
/* L240: */
#line 430 "zlasr.f"
		}
#line 431 "zlasr.f"
	    }
#line 432 "zlasr.f"
	}
#line 433 "zlasr.f"
    }

#line 435 "zlasr.f"
    return 0;

/*     End of ZLASR */

} /* zlasr_ */

