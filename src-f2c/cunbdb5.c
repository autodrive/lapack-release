#line 1 "cunbdb5.f"
/* cunbdb5.f -- translated by f2c (version 20100827).
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

#line 1 "cunbdb5.f"
/* > \brief \b CUNBDB5 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNBDB5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb5
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb5
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb5
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
/*                           LDQ2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, */
/*      $                   N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > CUNBDB5 orthogonalizes the column vector */
/* >      X = [ X1 ] */
/* >          [ X2 ] */
/* > with respect to the columns of */
/* >      Q = [ Q1 ] . */
/* >          [ Q2 ] */
/* > The columns of Q must be orthonormal. */
/* > */
/* > If the projection is zero according to Kahan's "twice is enough" */
/* > criterion, then some other vector from the orthogonal complement */
/* > is returned. This vector is chosen in an arbitrary but deterministic */
/* > way. */
/* > */
/* >\endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M1 */
/* > \verbatim */
/* >          M1 is INTEGER */
/* >           The dimension of X1 and the number of rows in Q1. 0 <= M1. */
/* > \endverbatim */
/* > */
/* > \param[in] M2 */
/* > \verbatim */
/* >          M2 is INTEGER */
/* >           The dimension of X2 and the number of rows in Q2. 0 <= M2. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           The number of columns in Q1 and Q2. 0 <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X1 */
/* > \verbatim */
/* >          X1 is COMPLEX array, dimension (M1) */
/* >           On entry, the top part of the vector to be orthogonalized. */
/* >           On exit, the top part of the projected vector. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX1 */
/* > \verbatim */
/* >          INCX1 is INTEGER */
/* >           Increment for entries of X1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X2 */
/* > \verbatim */
/* >          X2 is COMPLEX array, dimension (M2) */
/* >           On entry, the bottom part of the vector to be */
/* >           orthogonalized. On exit, the bottom part of the projected */
/* >           vector. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX2 */
/* > \verbatim */
/* >          INCX2 is INTEGER */
/* >           Increment for entries of X2. */
/* > \endverbatim */
/* > */
/* > \param[in] Q1 */
/* > \verbatim */
/* >          Q1 is COMPLEX array, dimension (LDQ1, N) */
/* >           The top part of the orthonormal basis matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ1 */
/* > \verbatim */
/* >          LDQ1 is INTEGER */
/* >           The leading dimension of Q1. LDQ1 >= M1. */
/* > \endverbatim */
/* > */
/* > \param[in] Q2 */
/* > \verbatim */
/* >          Q2 is COMPLEX array, dimension (LDQ2, N) */
/* >           The bottom part of the orthonormal basis matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ2 */
/* > \verbatim */
/* >          LDQ2 is INTEGER */
/* >           The leading dimension of Q2. LDQ2 >= M2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >           The dimension of the array WORK. LWORK >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:  successful exit. */
/* >           < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date July 2012 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunbdb5_(integer *m1, integer *m2, integer *n, 
	doublecomplex *x1, integer *incx1, doublecomplex *x2, integer *incx2, 
	doublecomplex *q1, integer *ldq1, doublecomplex *q2, integer *ldq2, 
	doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, childinfo;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cunbdb6_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *)
	    ;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     July 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Function .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input arguments */

#line 195 "cunbdb5.f"
    /* Parameter adjustments */
#line 195 "cunbdb5.f"
    --x1;
#line 195 "cunbdb5.f"
    --x2;
#line 195 "cunbdb5.f"
    q1_dim1 = *ldq1;
#line 195 "cunbdb5.f"
    q1_offset = 1 + q1_dim1;
#line 195 "cunbdb5.f"
    q1 -= q1_offset;
#line 195 "cunbdb5.f"
    q2_dim1 = *ldq2;
#line 195 "cunbdb5.f"
    q2_offset = 1 + q2_dim1;
#line 195 "cunbdb5.f"
    q2 -= q2_offset;
#line 195 "cunbdb5.f"
    --work;
#line 195 "cunbdb5.f"

#line 195 "cunbdb5.f"
    /* Function Body */
#line 195 "cunbdb5.f"
    *info = 0;
#line 196 "cunbdb5.f"
    if (*m1 < 0) {
#line 197 "cunbdb5.f"
	*info = -1;
#line 198 "cunbdb5.f"
    } else if (*m2 < 0) {
#line 199 "cunbdb5.f"
	*info = -2;
#line 200 "cunbdb5.f"
    } else if (*n < 0) {
#line 201 "cunbdb5.f"
	*info = -3;
#line 202 "cunbdb5.f"
    } else if (*incx1 < 1) {
#line 203 "cunbdb5.f"
	*info = -5;
#line 204 "cunbdb5.f"
    } else if (*incx2 < 1) {
#line 205 "cunbdb5.f"
	*info = -7;
#line 206 "cunbdb5.f"
    } else if (*ldq1 < max(1,*m1)) {
#line 207 "cunbdb5.f"
	*info = -9;
#line 208 "cunbdb5.f"
    } else if (*ldq2 < max(1,*m2)) {
#line 209 "cunbdb5.f"
	*info = -11;
#line 210 "cunbdb5.f"
    } else if (*lwork < *n) {
#line 211 "cunbdb5.f"
	*info = -13;
#line 212 "cunbdb5.f"
    }

#line 214 "cunbdb5.f"
    if (*info != 0) {
#line 215 "cunbdb5.f"
	i__1 = -(*info);
#line 215 "cunbdb5.f"
	xerbla_("CUNBDB5", &i__1, (ftnlen)7);
#line 216 "cunbdb5.f"
	return 0;
#line 217 "cunbdb5.f"
    }

/*     Project X onto the orthogonal complement of Q */

#line 221 "cunbdb5.f"
    cunbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], ldq1, &
	    q2[q2_offset], ldq2, &work[1], lwork, &childinfo);

/*     If the projection is nonzero, then return */

#line 226 "cunbdb5.f"
    d__1 = scnrm2_(m1, &x1[1], incx1);
#line 226 "cunbdb5.f"
    d__2 = scnrm2_(m2, &x2[1], incx2);
#line 226 "cunbdb5.f"
    if (d__1 != 0. || d__2 != 0.) {
#line 228 "cunbdb5.f"
	return 0;
#line 229 "cunbdb5.f"
    }

/*     Project each standard basis vector e_1,...,e_M1 in turn, stopping */
/*     when a nonzero projection is found */

#line 234 "cunbdb5.f"
    i__1 = *m1;
#line 234 "cunbdb5.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "cunbdb5.f"
	i__2 = *m1;
#line 235 "cunbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 236 "cunbdb5.f"
	    i__3 = j;
#line 236 "cunbdb5.f"
	    x1[i__3].r = 0., x1[i__3].i = 0.;
#line 237 "cunbdb5.f"
	}
#line 238 "cunbdb5.f"
	i__2 = i__;
#line 238 "cunbdb5.f"
	x1[i__2].r = 1., x1[i__2].i = 0.;
#line 239 "cunbdb5.f"
	i__2 = *m2;
#line 239 "cunbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 240 "cunbdb5.f"
	    i__3 = j;
#line 240 "cunbdb5.f"
	    x2[i__3].r = 0., x2[i__3].i = 0.;
#line 241 "cunbdb5.f"
	}
#line 242 "cunbdb5.f"
	cunbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], 
		ldq1, &q2[q2_offset], ldq2, &work[1], lwork, &childinfo);
#line 244 "cunbdb5.f"
	d__1 = scnrm2_(m1, &x1[1], incx1);
#line 244 "cunbdb5.f"
	d__2 = scnrm2_(m2, &x2[1], incx2);
#line 244 "cunbdb5.f"
	if (d__1 != 0. || d__2 != 0.) {
#line 246 "cunbdb5.f"
	    return 0;
#line 247 "cunbdb5.f"
	}
#line 248 "cunbdb5.f"
    }

/*     Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn, */
/*     stopping when a nonzero projection is found */

#line 253 "cunbdb5.f"
    i__1 = *m2;
#line 253 "cunbdb5.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 254 "cunbdb5.f"
	i__2 = *m1;
#line 254 "cunbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 255 "cunbdb5.f"
	    i__3 = j;
#line 255 "cunbdb5.f"
	    x1[i__3].r = 0., x1[i__3].i = 0.;
#line 256 "cunbdb5.f"
	}
#line 257 "cunbdb5.f"
	i__2 = *m2;
#line 257 "cunbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 258 "cunbdb5.f"
	    i__3 = j;
#line 258 "cunbdb5.f"
	    x2[i__3].r = 0., x2[i__3].i = 0.;
#line 259 "cunbdb5.f"
	}
#line 260 "cunbdb5.f"
	i__2 = i__;
#line 260 "cunbdb5.f"
	x2[i__2].r = 1., x2[i__2].i = 0.;
#line 261 "cunbdb5.f"
	cunbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], 
		ldq1, &q2[q2_offset], ldq2, &work[1], lwork, &childinfo);
#line 263 "cunbdb5.f"
	d__1 = scnrm2_(m1, &x1[1], incx1);
#line 263 "cunbdb5.f"
	d__2 = scnrm2_(m2, &x2[1], incx2);
#line 263 "cunbdb5.f"
	if (d__1 != 0. || d__2 != 0.) {
#line 265 "cunbdb5.f"
	    return 0;
#line 266 "cunbdb5.f"
	}
#line 267 "cunbdb5.f"
    }

#line 269 "cunbdb5.f"
    return 0;

/*     End of CUNBDB5 */

} /* cunbdb5_ */

