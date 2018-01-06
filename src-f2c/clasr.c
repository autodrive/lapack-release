#line 1 "clasr.f"
/* clasr.f -- translated by f2c (version 20100827).
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

#line 1 "clasr.f"
/* > \brief \b CLASR applies a sequence of plane rotations to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasr.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasr.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasr.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, PIVOT, SIDE */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ), S( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASR applies a sequence of real plane rotations to a complex matrix */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clasr_(char *side, char *pivot, char *direct, integer *m,
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 242 "clasr.f"
    /* Parameter adjustments */
#line 242 "clasr.f"
    --c__;
#line 242 "clasr.f"
    --s;
#line 242 "clasr.f"
    a_dim1 = *lda;
#line 242 "clasr.f"
    a_offset = 1 + a_dim1;
#line 242 "clasr.f"
    a -= a_offset;
#line 242 "clasr.f"

#line 242 "clasr.f"
    /* Function Body */
#line 242 "clasr.f"
    info = 0;
#line 243 "clasr.f"
    if (! (lsame_(side, "L", (ftnlen)1, (ftnlen)1) || lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1))) {
#line 244 "clasr.f"
	info = 1;
#line 245 "clasr.f"
    } else if (! (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1) || lsame_(pivot, 
	    "T", (ftnlen)1, (ftnlen)1) || lsame_(pivot, "B", (ftnlen)1, (
	    ftnlen)1))) {
#line 247 "clasr.f"
	info = 2;
#line 248 "clasr.f"
    } else if (! (lsame_(direct, "F", (ftnlen)1, (ftnlen)1) || lsame_(direct, 
	    "B", (ftnlen)1, (ftnlen)1))) {
#line 250 "clasr.f"
	info = 3;
#line 251 "clasr.f"
    } else if (*m < 0) {
#line 252 "clasr.f"
	info = 4;
#line 253 "clasr.f"
    } else if (*n < 0) {
#line 254 "clasr.f"
	info = 5;
#line 255 "clasr.f"
    } else if (*lda < max(1,*m)) {
#line 256 "clasr.f"
	info = 9;
#line 257 "clasr.f"
    }
#line 258 "clasr.f"
    if (info != 0) {
#line 259 "clasr.f"
	xerbla_("CLASR ", &info, (ftnlen)6);
#line 260 "clasr.f"
	return 0;
#line 261 "clasr.f"
    }

/*     Quick return if possible */

#line 265 "clasr.f"
    if (*m == 0 || *n == 0) {
#line 265 "clasr.f"
	return 0;
#line 265 "clasr.f"
    }
#line 267 "clasr.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  P * A */

#line 271 "clasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 272 "clasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 273 "clasr.f"
		i__1 = *m - 1;
#line 273 "clasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "clasr.f"
		    ctemp = c__[j];
#line 275 "clasr.f"
		    stemp = s[j];
#line 276 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 277 "clasr.f"
			i__2 = *n;
#line 277 "clasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "clasr.f"
			    i__3 = j + 1 + i__ * a_dim1;
#line 278 "clasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 279 "clasr.f"
			    i__3 = j + 1 + i__ * a_dim1;
#line 279 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 279 "clasr.f"
			    i__4 = j + i__ * a_dim1;
#line 279 "clasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 279 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 279 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 280 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 280 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 280 "clasr.f"
			    i__4 = j + i__ * a_dim1;
#line 280 "clasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 280 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 280 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 281 "clasr.f"
/* L10: */
#line 281 "clasr.f"
			}
#line 282 "clasr.f"
		    }
#line 283 "clasr.f"
/* L20: */
#line 283 "clasr.f"
		}
#line 284 "clasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 285 "clasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 286 "clasr.f"
		    ctemp = c__[j];
#line 287 "clasr.f"
		    stemp = s[j];
#line 288 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 289 "clasr.f"
			i__1 = *n;
#line 289 "clasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "clasr.f"
			    i__2 = j + 1 + i__ * a_dim1;
#line 290 "clasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 291 "clasr.f"
			    i__2 = j + 1 + i__ * a_dim1;
#line 291 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 291 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 291 "clasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 291 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 291 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 292 "clasr.f"
			    i__2 = j + i__ * a_dim1;
#line 292 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 292 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 292 "clasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 292 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 292 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 293 "clasr.f"
/* L30: */
#line 293 "clasr.f"
			}
#line 294 "clasr.f"
		    }
#line 295 "clasr.f"
/* L40: */
#line 295 "clasr.f"
		}
#line 296 "clasr.f"
	    }
#line 297 "clasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 298 "clasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 299 "clasr.f"
		i__1 = *m;
#line 299 "clasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 300 "clasr.f"
		    ctemp = c__[j - 1];
#line 301 "clasr.f"
		    stemp = s[j - 1];
#line 302 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 303 "clasr.f"
			i__2 = *n;
#line 303 "clasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 304 "clasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 305 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 305 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 305 "clasr.f"
			    i__4 = i__ * a_dim1 + 1;
#line 305 "clasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 305 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 305 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 306 "clasr.f"
			    i__3 = i__ * a_dim1 + 1;
#line 306 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 306 "clasr.f"
			    i__4 = i__ * a_dim1 + 1;
#line 306 "clasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 306 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 306 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 307 "clasr.f"
/* L50: */
#line 307 "clasr.f"
			}
#line 308 "clasr.f"
		    }
#line 309 "clasr.f"
/* L60: */
#line 309 "clasr.f"
		}
#line 310 "clasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 311 "clasr.f"
		for (j = *m; j >= 2; --j) {
#line 312 "clasr.f"
		    ctemp = c__[j - 1];
#line 313 "clasr.f"
		    stemp = s[j - 1];
#line 314 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 315 "clasr.f"
			i__1 = *n;
#line 315 "clasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "clasr.f"
			    i__2 = j + i__ * a_dim1;
#line 316 "clasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 317 "clasr.f"
			    i__2 = j + i__ * a_dim1;
#line 317 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 317 "clasr.f"
			    i__3 = i__ * a_dim1 + 1;
#line 317 "clasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 317 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 317 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 318 "clasr.f"
			    i__2 = i__ * a_dim1 + 1;
#line 318 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 318 "clasr.f"
			    i__3 = i__ * a_dim1 + 1;
#line 318 "clasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 318 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 318 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 319 "clasr.f"
/* L70: */
#line 319 "clasr.f"
			}
#line 320 "clasr.f"
		    }
#line 321 "clasr.f"
/* L80: */
#line 321 "clasr.f"
		}
#line 322 "clasr.f"
	    }
#line 323 "clasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 324 "clasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 325 "clasr.f"
		i__1 = *m - 1;
#line 325 "clasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "clasr.f"
		    ctemp = c__[j];
#line 327 "clasr.f"
		    stemp = s[j];
#line 328 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 329 "clasr.f"
			i__2 = *n;
#line 329 "clasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 330 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 330 "clasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 331 "clasr.f"
			    i__3 = j + i__ * a_dim1;
#line 331 "clasr.f"
			    i__4 = *m + i__ * a_dim1;
#line 331 "clasr.f"
			    z__2.r = stemp * a[i__4].r, z__2.i = stemp * a[
				    i__4].i;
#line 331 "clasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 331 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 331 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 332 "clasr.f"
			    i__3 = *m + i__ * a_dim1;
#line 332 "clasr.f"
			    i__4 = *m + i__ * a_dim1;
#line 332 "clasr.f"
			    z__2.r = ctemp * a[i__4].r, z__2.i = ctemp * a[
				    i__4].i;
#line 332 "clasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 332 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 332 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 333 "clasr.f"
/* L90: */
#line 333 "clasr.f"
			}
#line 334 "clasr.f"
		    }
#line 335 "clasr.f"
/* L100: */
#line 335 "clasr.f"
		}
#line 336 "clasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 337 "clasr.f"
		for (j = *m - 1; j >= 1; --j) {
#line 338 "clasr.f"
		    ctemp = c__[j];
#line 339 "clasr.f"
		    stemp = s[j];
#line 340 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 341 "clasr.f"
			i__1 = *n;
#line 341 "clasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 342 "clasr.f"
			    i__2 = j + i__ * a_dim1;
#line 342 "clasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 343 "clasr.f"
			    i__2 = j + i__ * a_dim1;
#line 343 "clasr.f"
			    i__3 = *m + i__ * a_dim1;
#line 343 "clasr.f"
			    z__2.r = stemp * a[i__3].r, z__2.i = stemp * a[
				    i__3].i;
#line 343 "clasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 343 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 343 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 344 "clasr.f"
			    i__2 = *m + i__ * a_dim1;
#line 344 "clasr.f"
			    i__3 = *m + i__ * a_dim1;
#line 344 "clasr.f"
			    z__2.r = ctemp * a[i__3].r, z__2.i = ctemp * a[
				    i__3].i;
#line 344 "clasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 344 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 344 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 345 "clasr.f"
/* L110: */
#line 345 "clasr.f"
			}
#line 346 "clasr.f"
		    }
#line 347 "clasr.f"
/* L120: */
#line 347 "clasr.f"
		}
#line 348 "clasr.f"
	    }
#line 349 "clasr.f"
	}
#line 350 "clasr.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*        Form A * P**T */

#line 354 "clasr.f"
	if (lsame_(pivot, "V", (ftnlen)1, (ftnlen)1)) {
#line 355 "clasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 356 "clasr.f"
		i__1 = *n - 1;
#line 356 "clasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 357 "clasr.f"
		    ctemp = c__[j];
#line 358 "clasr.f"
		    stemp = s[j];
#line 359 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 360 "clasr.f"
			i__2 = *m;
#line 360 "clasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 361 "clasr.f"
			    i__3 = i__ + (j + 1) * a_dim1;
#line 361 "clasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 362 "clasr.f"
			    i__3 = i__ + (j + 1) * a_dim1;
#line 362 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 362 "clasr.f"
			    i__4 = i__ + j * a_dim1;
#line 362 "clasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 362 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 362 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 363 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 363 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 363 "clasr.f"
			    i__4 = i__ + j * a_dim1;
#line 363 "clasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 363 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 363 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 364 "clasr.f"
/* L130: */
#line 364 "clasr.f"
			}
#line 365 "clasr.f"
		    }
#line 366 "clasr.f"
/* L140: */
#line 366 "clasr.f"
		}
#line 367 "clasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 368 "clasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 369 "clasr.f"
		    ctemp = c__[j];
#line 370 "clasr.f"
		    stemp = s[j];
#line 371 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 372 "clasr.f"
			i__1 = *m;
#line 372 "clasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 373 "clasr.f"
			    i__2 = i__ + (j + 1) * a_dim1;
#line 373 "clasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 374 "clasr.f"
			    i__2 = i__ + (j + 1) * a_dim1;
#line 374 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 374 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 374 "clasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 374 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 374 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 375 "clasr.f"
			    i__2 = i__ + j * a_dim1;
#line 375 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 375 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 375 "clasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 375 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 375 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 376 "clasr.f"
/* L150: */
#line 376 "clasr.f"
			}
#line 377 "clasr.f"
		    }
#line 378 "clasr.f"
/* L160: */
#line 378 "clasr.f"
		}
#line 379 "clasr.f"
	    }
#line 380 "clasr.f"
	} else if (lsame_(pivot, "T", (ftnlen)1, (ftnlen)1)) {
#line 381 "clasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 382 "clasr.f"
		i__1 = *n;
#line 382 "clasr.f"
		for (j = 2; j <= i__1; ++j) {
#line 383 "clasr.f"
		    ctemp = c__[j - 1];
#line 384 "clasr.f"
		    stemp = s[j - 1];
#line 385 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 386 "clasr.f"
			i__2 = *m;
#line 386 "clasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 387 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 387 "clasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 388 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 388 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 388 "clasr.f"
			    i__4 = i__ + a_dim1;
#line 388 "clasr.f"
			    z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[
				    i__4].i;
#line 388 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 388 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 389 "clasr.f"
			    i__3 = i__ + a_dim1;
#line 389 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 389 "clasr.f"
			    i__4 = i__ + a_dim1;
#line 389 "clasr.f"
			    z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[
				    i__4].i;
#line 389 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 389 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 390 "clasr.f"
/* L170: */
#line 390 "clasr.f"
			}
#line 391 "clasr.f"
		    }
#line 392 "clasr.f"
/* L180: */
#line 392 "clasr.f"
		}
#line 393 "clasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 394 "clasr.f"
		for (j = *n; j >= 2; --j) {
#line 395 "clasr.f"
		    ctemp = c__[j - 1];
#line 396 "clasr.f"
		    stemp = s[j - 1];
#line 397 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 398 "clasr.f"
			i__1 = *m;
#line 398 "clasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 399 "clasr.f"
			    i__2 = i__ + j * a_dim1;
#line 399 "clasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 400 "clasr.f"
			    i__2 = i__ + j * a_dim1;
#line 400 "clasr.f"
			    z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
#line 400 "clasr.f"
			    i__3 = i__ + a_dim1;
#line 400 "clasr.f"
			    z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[
				    i__3].i;
#line 400 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 400 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 401 "clasr.f"
			    i__2 = i__ + a_dim1;
#line 401 "clasr.f"
			    z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
#line 401 "clasr.f"
			    i__3 = i__ + a_dim1;
#line 401 "clasr.f"
			    z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[
				    i__3].i;
#line 401 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 401 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 402 "clasr.f"
/* L190: */
#line 402 "clasr.f"
			}
#line 403 "clasr.f"
		    }
#line 404 "clasr.f"
/* L200: */
#line 404 "clasr.f"
		}
#line 405 "clasr.f"
	    }
#line 406 "clasr.f"
	} else if (lsame_(pivot, "B", (ftnlen)1, (ftnlen)1)) {
#line 407 "clasr.f"
	    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 408 "clasr.f"
		i__1 = *n - 1;
#line 408 "clasr.f"
		for (j = 1; j <= i__1; ++j) {
#line 409 "clasr.f"
		    ctemp = c__[j];
#line 410 "clasr.f"
		    stemp = s[j];
#line 411 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 412 "clasr.f"
			i__2 = *m;
#line 412 "clasr.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 413 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 413 "clasr.f"
			    temp.r = a[i__3].r, temp.i = a[i__3].i;
#line 414 "clasr.f"
			    i__3 = i__ + j * a_dim1;
#line 414 "clasr.f"
			    i__4 = i__ + *n * a_dim1;
#line 414 "clasr.f"
			    z__2.r = stemp * a[i__4].r, z__2.i = stemp * a[
				    i__4].i;
#line 414 "clasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 414 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 414 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 415 "clasr.f"
			    i__3 = i__ + *n * a_dim1;
#line 415 "clasr.f"
			    i__4 = i__ + *n * a_dim1;
#line 415 "clasr.f"
			    z__2.r = ctemp * a[i__4].r, z__2.i = ctemp * a[
				    i__4].i;
#line 415 "clasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 415 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 415 "clasr.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 416 "clasr.f"
/* L210: */
#line 416 "clasr.f"
			}
#line 417 "clasr.f"
		    }
#line 418 "clasr.f"
/* L220: */
#line 418 "clasr.f"
		}
#line 419 "clasr.f"
	    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 420 "clasr.f"
		for (j = *n - 1; j >= 1; --j) {
#line 421 "clasr.f"
		    ctemp = c__[j];
#line 422 "clasr.f"
		    stemp = s[j];
#line 423 "clasr.f"
		    if (ctemp != 1. || stemp != 0.) {
#line 424 "clasr.f"
			i__1 = *m;
#line 424 "clasr.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 425 "clasr.f"
			    i__2 = i__ + j * a_dim1;
#line 425 "clasr.f"
			    temp.r = a[i__2].r, temp.i = a[i__2].i;
#line 426 "clasr.f"
			    i__2 = i__ + j * a_dim1;
#line 426 "clasr.f"
			    i__3 = i__ + *n * a_dim1;
#line 426 "clasr.f"
			    z__2.r = stemp * a[i__3].r, z__2.i = stemp * a[
				    i__3].i;
#line 426 "clasr.f"
			    z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
#line 426 "clasr.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 426 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 427 "clasr.f"
			    i__2 = i__ + *n * a_dim1;
#line 427 "clasr.f"
			    i__3 = i__ + *n * a_dim1;
#line 427 "clasr.f"
			    z__2.r = ctemp * a[i__3].r, z__2.i = ctemp * a[
				    i__3].i;
#line 427 "clasr.f"
			    z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
#line 427 "clasr.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 427 "clasr.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 428 "clasr.f"
/* L230: */
#line 428 "clasr.f"
			}
#line 429 "clasr.f"
		    }
#line 430 "clasr.f"
/* L240: */
#line 430 "clasr.f"
		}
#line 431 "clasr.f"
	    }
#line 432 "clasr.f"
	}
#line 433 "clasr.f"
    }

#line 435 "clasr.f"
    return 0;

/*     End of CLASR */

} /* clasr_ */

