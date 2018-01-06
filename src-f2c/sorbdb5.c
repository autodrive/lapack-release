#line 1 "sorbdb5.f"
/* sorbdb5.f -- translated by f2c (version 20100827).
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

#line 1 "sorbdb5.f"
/* > \brief \b SORBDB5 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORBDB5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb5
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb5
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb5
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
/*                           LDQ2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, */
/*      $                   N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > SORBDB5 orthogonalizes the column vector */
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
/* >          X1 is REAL array, dimension (M1) */
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
/* >          X2 is REAL array, dimension (M2) */
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
/* >          Q1 is REAL array, dimension (LDQ1, N) */
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
/* >          Q2 is REAL array, dimension (LDQ2, N) */
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
/* >          WORK is REAL array, dimension (LWORK) */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorbdb5_(integer *m1, integer *m2, integer *n, 
	doublereal *x1, integer *incx1, doublereal *x2, integer *incx2, 
	doublereal *q1, integer *ldq1, doublereal *q2, integer *ldq2, 
	doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, childinfo;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), sorbdb6_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);


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

#line 195 "sorbdb5.f"
    /* Parameter adjustments */
#line 195 "sorbdb5.f"
    --x1;
#line 195 "sorbdb5.f"
    --x2;
#line 195 "sorbdb5.f"
    q1_dim1 = *ldq1;
#line 195 "sorbdb5.f"
    q1_offset = 1 + q1_dim1;
#line 195 "sorbdb5.f"
    q1 -= q1_offset;
#line 195 "sorbdb5.f"
    q2_dim1 = *ldq2;
#line 195 "sorbdb5.f"
    q2_offset = 1 + q2_dim1;
#line 195 "sorbdb5.f"
    q2 -= q2_offset;
#line 195 "sorbdb5.f"
    --work;
#line 195 "sorbdb5.f"

#line 195 "sorbdb5.f"
    /* Function Body */
#line 195 "sorbdb5.f"
    *info = 0;
#line 196 "sorbdb5.f"
    if (*m1 < 0) {
#line 197 "sorbdb5.f"
	*info = -1;
#line 198 "sorbdb5.f"
    } else if (*m2 < 0) {
#line 199 "sorbdb5.f"
	*info = -2;
#line 200 "sorbdb5.f"
    } else if (*n < 0) {
#line 201 "sorbdb5.f"
	*info = -3;
#line 202 "sorbdb5.f"
    } else if (*incx1 < 1) {
#line 203 "sorbdb5.f"
	*info = -5;
#line 204 "sorbdb5.f"
    } else if (*incx2 < 1) {
#line 205 "sorbdb5.f"
	*info = -7;
#line 206 "sorbdb5.f"
    } else if (*ldq1 < max(1,*m1)) {
#line 207 "sorbdb5.f"
	*info = -9;
#line 208 "sorbdb5.f"
    } else if (*ldq2 < max(1,*m2)) {
#line 209 "sorbdb5.f"
	*info = -11;
#line 210 "sorbdb5.f"
    } else if (*lwork < *n) {
#line 211 "sorbdb5.f"
	*info = -13;
#line 212 "sorbdb5.f"
    }

#line 214 "sorbdb5.f"
    if (*info != 0) {
#line 215 "sorbdb5.f"
	i__1 = -(*info);
#line 215 "sorbdb5.f"
	xerbla_("SORBDB5", &i__1, (ftnlen)7);
#line 216 "sorbdb5.f"
	return 0;
#line 217 "sorbdb5.f"
    }

/*     Project X onto the orthogonal complement of Q */

#line 221 "sorbdb5.f"
    sorbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], ldq1, &
	    q2[q2_offset], ldq2, &work[1], lwork, &childinfo);

/*     If the projection is nonzero, then return */

#line 226 "sorbdb5.f"
    if (snrm2_(m1, &x1[1], incx1) != 0. || snrm2_(m2, &x2[1], incx2) != 0.) {
#line 228 "sorbdb5.f"
	return 0;
#line 229 "sorbdb5.f"
    }

/*     Project each standard basis vector e_1,...,e_M1 in turn, stopping */
/*     when a nonzero projection is found */

#line 234 "sorbdb5.f"
    i__1 = *m1;
#line 234 "sorbdb5.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "sorbdb5.f"
	i__2 = *m1;
#line 235 "sorbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 236 "sorbdb5.f"
	    x1[j] = 0.;
#line 237 "sorbdb5.f"
	}
#line 238 "sorbdb5.f"
	x1[i__] = 1.;
#line 239 "sorbdb5.f"
	i__2 = *m2;
#line 239 "sorbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 240 "sorbdb5.f"
	    x2[j] = 0.;
#line 241 "sorbdb5.f"
	}
#line 242 "sorbdb5.f"
	sorbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], 
		ldq1, &q2[q2_offset], ldq2, &work[1], lwork, &childinfo);
#line 244 "sorbdb5.f"
	if (snrm2_(m1, &x1[1], incx1) != 0. || snrm2_(m2, &x2[1], incx2) != 
		0.) {
#line 246 "sorbdb5.f"
	    return 0;
#line 247 "sorbdb5.f"
	}
#line 248 "sorbdb5.f"
    }

/*     Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn, */
/*     stopping when a nonzero projection is found */

#line 253 "sorbdb5.f"
    i__1 = *m2;
#line 253 "sorbdb5.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 254 "sorbdb5.f"
	i__2 = *m1;
#line 254 "sorbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 255 "sorbdb5.f"
	    x1[j] = 0.;
#line 256 "sorbdb5.f"
	}
#line 257 "sorbdb5.f"
	i__2 = *m2;
#line 257 "sorbdb5.f"
	for (j = 1; j <= i__2; ++j) {
#line 258 "sorbdb5.f"
	    x2[j] = 0.;
#line 259 "sorbdb5.f"
	}
#line 260 "sorbdb5.f"
	x2[i__] = 1.;
#line 261 "sorbdb5.f"
	sorbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], 
		ldq1, &q2[q2_offset], ldq2, &work[1], lwork, &childinfo);
#line 263 "sorbdb5.f"
	if (snrm2_(m1, &x1[1], incx1) != 0. || snrm2_(m2, &x2[1], incx2) != 
		0.) {
#line 265 "sorbdb5.f"
	    return 0;
#line 266 "sorbdb5.f"
	}
#line 267 "sorbdb5.f"
    }

#line 269 "sorbdb5.f"
    return 0;

/*     End of SORBDB5 */

} /* sorbdb5_ */

