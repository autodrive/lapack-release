#line 1 "sorbdb6.f"
/* sorbdb6.f -- translated by f2c (version 20100827).
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

#line 1 "sorbdb6.f"
/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static doublereal c_b12 = -1.;

/* > \brief \b SORBDB6 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORBDB6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb6
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb6
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb6
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
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
/* > SORBDB6 orthogonalizes the column vector */
/* >      X = [ X1 ] */
/* >          [ X2 ] */
/* > with respect to the columns of */
/* >      Q = [ Q1 ] . */
/* >          [ Q2 ] */
/* > The columns of Q must be orthonormal. */
/* > */
/* > If the projection is zero according to Kahan's "twice is enough" */
/* > criterion, then the zero vector is returned. */
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
/* Subroutine */ int sorbdb6_(integer *m1, integer *m2, integer *n, 
	doublereal *x1, integer *incx1, doublereal *x2, integer *incx2, 
	doublereal *q1, integer *ldq1, doublereal *q2, integer *ldq2, 
	doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal scl1, scl2, ssq1, ssq2;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), slassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal normsq1, normsq2;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. Intrinsic Function .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input arguments */

#line 193 "sorbdb6.f"
    /* Parameter adjustments */
#line 193 "sorbdb6.f"
    --x1;
#line 193 "sorbdb6.f"
    --x2;
#line 193 "sorbdb6.f"
    q1_dim1 = *ldq1;
#line 193 "sorbdb6.f"
    q1_offset = 1 + q1_dim1;
#line 193 "sorbdb6.f"
    q1 -= q1_offset;
#line 193 "sorbdb6.f"
    q2_dim1 = *ldq2;
#line 193 "sorbdb6.f"
    q2_offset = 1 + q2_dim1;
#line 193 "sorbdb6.f"
    q2 -= q2_offset;
#line 193 "sorbdb6.f"
    --work;
#line 193 "sorbdb6.f"

#line 193 "sorbdb6.f"
    /* Function Body */
#line 193 "sorbdb6.f"
    *info = 0;
#line 194 "sorbdb6.f"
    if (*m1 < 0) {
#line 195 "sorbdb6.f"
	*info = -1;
#line 196 "sorbdb6.f"
    } else if (*m2 < 0) {
#line 197 "sorbdb6.f"
	*info = -2;
#line 198 "sorbdb6.f"
    } else if (*n < 0) {
#line 199 "sorbdb6.f"
	*info = -3;
#line 200 "sorbdb6.f"
    } else if (*incx1 < 1) {
#line 201 "sorbdb6.f"
	*info = -5;
#line 202 "sorbdb6.f"
    } else if (*incx2 < 1) {
#line 203 "sorbdb6.f"
	*info = -7;
#line 204 "sorbdb6.f"
    } else if (*ldq1 < max(1,*m1)) {
#line 205 "sorbdb6.f"
	*info = -9;
#line 206 "sorbdb6.f"
    } else if (*ldq2 < max(1,*m2)) {
#line 207 "sorbdb6.f"
	*info = -11;
#line 208 "sorbdb6.f"
    } else if (*lwork < *n) {
#line 209 "sorbdb6.f"
	*info = -13;
#line 210 "sorbdb6.f"
    }

#line 212 "sorbdb6.f"
    if (*info != 0) {
#line 213 "sorbdb6.f"
	i__1 = -(*info);
#line 213 "sorbdb6.f"
	xerbla_("SORBDB6", &i__1, (ftnlen)7);
#line 214 "sorbdb6.f"
	return 0;
#line 215 "sorbdb6.f"
    }

/*     First, project X onto the orthogonal complement of Q's column */
/*     space */

#line 220 "sorbdb6.f"
    scl1 = 0.;
#line 221 "sorbdb6.f"
    ssq1 = 1.;
#line 222 "sorbdb6.f"
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
#line 223 "sorbdb6.f"
    scl2 = 0.;
#line 224 "sorbdb6.f"
    ssq2 = 1.;
#line 225 "sorbdb6.f"
    slassq_(m2, &x2[1], incx2, &scl2, &ssq2);
/* Computing 2nd power */
#line 226 "sorbdb6.f"
    d__1 = scl1;
/* Computing 2nd power */
#line 226 "sorbdb6.f"
    d__2 = scl2;
#line 226 "sorbdb6.f"
    normsq1 = d__1 * d__1 * ssq1 + d__2 * d__2 * ssq2;

#line 228 "sorbdb6.f"
    if (*m1 == 0) {
#line 229 "sorbdb6.f"
	i__1 = *n;
#line 229 "sorbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 230 "sorbdb6.f"
	    work[i__] = 0.;
#line 231 "sorbdb6.f"
	}
#line 232 "sorbdb6.f"
    } else {
#line 233 "sorbdb6.f"
	sgemv_("C", m1, n, &c_b4, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b5, 
		&work[1], &c__1, (ftnlen)1);
#line 235 "sorbdb6.f"
    }

#line 237 "sorbdb6.f"
    sgemv_("C", m2, n, &c_b4, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b4, &
	    work[1], &c__1, (ftnlen)1);

#line 239 "sorbdb6.f"
    sgemv_("N", m1, n, &c_b12, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b4, &
	    x1[1], incx1, (ftnlen)1);
#line 241 "sorbdb6.f"
    sgemv_("N", m2, n, &c_b12, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b4, &
	    x2[1], incx2, (ftnlen)1);

#line 244 "sorbdb6.f"
    scl1 = 0.;
#line 245 "sorbdb6.f"
    ssq1 = 1.;
#line 246 "sorbdb6.f"
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
#line 247 "sorbdb6.f"
    scl2 = 0.;
#line 248 "sorbdb6.f"
    ssq2 = 1.;
#line 249 "sorbdb6.f"
    slassq_(m2, &x2[1], incx2, &scl2, &ssq2);
/* Computing 2nd power */
#line 250 "sorbdb6.f"
    d__1 = scl1;
/* Computing 2nd power */
#line 250 "sorbdb6.f"
    d__2 = scl2;
#line 250 "sorbdb6.f"
    normsq2 = d__1 * d__1 * ssq1 + d__2 * d__2 * ssq2;

/*     If projection is sufficiently large in norm, then stop. */
/*     If projection is zero, then stop. */
/*     Otherwise, project again. */

#line 256 "sorbdb6.f"
    if (normsq2 >= normsq1 * .01) {
#line 257 "sorbdb6.f"
	return 0;
#line 258 "sorbdb6.f"
    }

#line 260 "sorbdb6.f"
    if (normsq2 == 0.) {
#line 261 "sorbdb6.f"
	return 0;
#line 262 "sorbdb6.f"
    }

#line 264 "sorbdb6.f"
    normsq1 = normsq2;

#line 266 "sorbdb6.f"
    i__1 = *n;
#line 266 "sorbdb6.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "sorbdb6.f"
	work[i__] = 0.;
#line 268 "sorbdb6.f"
    }

#line 270 "sorbdb6.f"
    if (*m1 == 0) {
#line 271 "sorbdb6.f"
	i__1 = *n;
#line 271 "sorbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 272 "sorbdb6.f"
	    work[i__] = 0.;
#line 273 "sorbdb6.f"
	}
#line 274 "sorbdb6.f"
    } else {
#line 275 "sorbdb6.f"
	sgemv_("C", m1, n, &c_b4, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b5, 
		&work[1], &c__1, (ftnlen)1);
#line 277 "sorbdb6.f"
    }

#line 279 "sorbdb6.f"
    sgemv_("C", m2, n, &c_b4, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b4, &
	    work[1], &c__1, (ftnlen)1);

#line 281 "sorbdb6.f"
    sgemv_("N", m1, n, &c_b12, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b4, &
	    x1[1], incx1, (ftnlen)1);
#line 283 "sorbdb6.f"
    sgemv_("N", m2, n, &c_b12, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b4, &
	    x2[1], incx2, (ftnlen)1);

#line 286 "sorbdb6.f"
    scl1 = 0.;
#line 287 "sorbdb6.f"
    ssq1 = 1.;
#line 288 "sorbdb6.f"
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
#line 289 "sorbdb6.f"
    scl2 = 0.;
#line 290 "sorbdb6.f"
    ssq2 = 1.;
#line 291 "sorbdb6.f"
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
/* Computing 2nd power */
#line 292 "sorbdb6.f"
    d__1 = scl1;
/* Computing 2nd power */
#line 292 "sorbdb6.f"
    d__2 = scl2;
#line 292 "sorbdb6.f"
    normsq2 = d__1 * d__1 * ssq1 + d__2 * d__2 * ssq2;

/*     If second projection is sufficiently large in norm, then do */
/*     nothing more. Alternatively, if it shrunk significantly, then */
/*     truncate it to zero. */

#line 298 "sorbdb6.f"
    if (normsq2 < normsq1 * .01) {
#line 299 "sorbdb6.f"
	i__1 = *m1;
#line 299 "sorbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 300 "sorbdb6.f"
	    x1[i__] = 0.;
#line 301 "sorbdb6.f"
	}
#line 302 "sorbdb6.f"
	i__1 = *m2;
#line 302 "sorbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 303 "sorbdb6.f"
	    x2[i__] = 0.;
#line 304 "sorbdb6.f"
	}
#line 305 "sorbdb6.f"
    }

#line 307 "sorbdb6.f"
    return 0;

/*     End of SORBDB6 */

} /* sorbdb6_ */

