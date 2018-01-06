#line 1 "cunbdb6.f"
/* cunbdb6.f -- translated by f2c (version 20100827).
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

#line 1 "cunbdb6.f"
/* Table of constant values */

static doublecomplex c_b1 = {-1.,0.};
static doublecomplex c_b2 = {1.,0.};
static doublecomplex c_b3 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CUNBDB6 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNBDB6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb6
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb6
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb6
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
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
/* > CUNBDB6 orthogonalizes the column vector */
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
/* Subroutine */ int cunbdb6_(integer *m1, integer *m2, integer *n, 
	doublecomplex *x1, integer *incx1, doublecomplex *x2, integer *incx2, 
	doublecomplex *q1, integer *ldq1, doublecomplex *q2, integer *ldq2, 
	doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal scl1, scl2, ssq1, ssq2;
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), classq_(integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *);
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

#line 194 "cunbdb6.f"
    /* Parameter adjustments */
#line 194 "cunbdb6.f"
    --x1;
#line 194 "cunbdb6.f"
    --x2;
#line 194 "cunbdb6.f"
    q1_dim1 = *ldq1;
#line 194 "cunbdb6.f"
    q1_offset = 1 + q1_dim1;
#line 194 "cunbdb6.f"
    q1 -= q1_offset;
#line 194 "cunbdb6.f"
    q2_dim1 = *ldq2;
#line 194 "cunbdb6.f"
    q2_offset = 1 + q2_dim1;
#line 194 "cunbdb6.f"
    q2 -= q2_offset;
#line 194 "cunbdb6.f"
    --work;
#line 194 "cunbdb6.f"

#line 194 "cunbdb6.f"
    /* Function Body */
#line 194 "cunbdb6.f"
    *info = 0;
#line 195 "cunbdb6.f"
    if (*m1 < 0) {
#line 196 "cunbdb6.f"
	*info = -1;
#line 197 "cunbdb6.f"
    } else if (*m2 < 0) {
#line 198 "cunbdb6.f"
	*info = -2;
#line 199 "cunbdb6.f"
    } else if (*n < 0) {
#line 200 "cunbdb6.f"
	*info = -3;
#line 201 "cunbdb6.f"
    } else if (*incx1 < 1) {
#line 202 "cunbdb6.f"
	*info = -5;
#line 203 "cunbdb6.f"
    } else if (*incx2 < 1) {
#line 204 "cunbdb6.f"
	*info = -7;
#line 205 "cunbdb6.f"
    } else if (*ldq1 < max(1,*m1)) {
#line 206 "cunbdb6.f"
	*info = -9;
#line 207 "cunbdb6.f"
    } else if (*ldq2 < max(1,*m2)) {
#line 208 "cunbdb6.f"
	*info = -11;
#line 209 "cunbdb6.f"
    } else if (*lwork < *n) {
#line 210 "cunbdb6.f"
	*info = -13;
#line 211 "cunbdb6.f"
    }

#line 213 "cunbdb6.f"
    if (*info != 0) {
#line 214 "cunbdb6.f"
	i__1 = -(*info);
#line 214 "cunbdb6.f"
	xerbla_("CUNBDB6", &i__1, (ftnlen)7);
#line 215 "cunbdb6.f"
	return 0;
#line 216 "cunbdb6.f"
    }

/*     First, project X onto the orthogonal complement of Q's column */
/*     space */

#line 221 "cunbdb6.f"
    scl1 = 0.;
#line 222 "cunbdb6.f"
    ssq1 = 1.;
#line 223 "cunbdb6.f"
    classq_(m1, &x1[1], incx1, &scl1, &ssq1);
#line 224 "cunbdb6.f"
    scl2 = 0.;
#line 225 "cunbdb6.f"
    ssq2 = 1.;
#line 226 "cunbdb6.f"
    classq_(m2, &x2[1], incx2, &scl2, &ssq2);
/* Computing 2nd power */
#line 227 "cunbdb6.f"
    d__1 = scl1;
/* Computing 2nd power */
#line 227 "cunbdb6.f"
    d__2 = scl2;
#line 227 "cunbdb6.f"
    normsq1 = d__1 * d__1 * ssq1 + d__2 * d__2 * ssq2;

#line 229 "cunbdb6.f"
    if (*m1 == 0) {
#line 230 "cunbdb6.f"
	i__1 = *n;
#line 230 "cunbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 231 "cunbdb6.f"
	    i__2 = i__;
#line 231 "cunbdb6.f"
	    work[i__2].r = 0., work[i__2].i = 0.;
#line 232 "cunbdb6.f"
	}
#line 233 "cunbdb6.f"
    } else {
#line 234 "cunbdb6.f"
	cgemv_("C", m1, n, &c_b2, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b3, 
		&work[1], &c__1, (ftnlen)1);
#line 236 "cunbdb6.f"
    }

#line 238 "cunbdb6.f"
    cgemv_("C", m2, n, &c_b2, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b2, &
	    work[1], &c__1, (ftnlen)1);

#line 240 "cunbdb6.f"
    cgemv_("N", m1, n, &c_b1, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b2, &
	    x1[1], incx1, (ftnlen)1);
#line 242 "cunbdb6.f"
    cgemv_("N", m2, n, &c_b1, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b2, &
	    x2[1], incx2, (ftnlen)1);

#line 245 "cunbdb6.f"
    scl1 = 0.;
#line 246 "cunbdb6.f"
    ssq1 = 1.;
#line 247 "cunbdb6.f"
    classq_(m1, &x1[1], incx1, &scl1, &ssq1);
#line 248 "cunbdb6.f"
    scl2 = 0.;
#line 249 "cunbdb6.f"
    ssq2 = 1.;
#line 250 "cunbdb6.f"
    classq_(m2, &x2[1], incx2, &scl2, &ssq2);
/* Computing 2nd power */
#line 251 "cunbdb6.f"
    d__1 = scl1;
/* Computing 2nd power */
#line 251 "cunbdb6.f"
    d__2 = scl2;
#line 251 "cunbdb6.f"
    normsq2 = d__1 * d__1 * ssq1 + d__2 * d__2 * ssq2;

/*     If projection is sufficiently large in norm, then stop. */
/*     If projection is zero, then stop. */
/*     Otherwise, project again. */

#line 257 "cunbdb6.f"
    if (normsq2 >= normsq1 * .01) {
#line 258 "cunbdb6.f"
	return 0;
#line 259 "cunbdb6.f"
    }

#line 261 "cunbdb6.f"
    if (normsq2 == 0.) {
#line 262 "cunbdb6.f"
	return 0;
#line 263 "cunbdb6.f"
    }

#line 265 "cunbdb6.f"
    normsq1 = normsq2;

#line 267 "cunbdb6.f"
    i__1 = *n;
#line 267 "cunbdb6.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 268 "cunbdb6.f"
	i__2 = i__;
#line 268 "cunbdb6.f"
	work[i__2].r = 0., work[i__2].i = 0.;
#line 269 "cunbdb6.f"
    }

#line 271 "cunbdb6.f"
    if (*m1 == 0) {
#line 272 "cunbdb6.f"
	i__1 = *n;
#line 272 "cunbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "cunbdb6.f"
	    i__2 = i__;
#line 273 "cunbdb6.f"
	    work[i__2].r = 0., work[i__2].i = 0.;
#line 274 "cunbdb6.f"
	}
#line 275 "cunbdb6.f"
    } else {
#line 276 "cunbdb6.f"
	cgemv_("C", m1, n, &c_b2, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b3, 
		&work[1], &c__1, (ftnlen)1);
#line 278 "cunbdb6.f"
    }

#line 280 "cunbdb6.f"
    cgemv_("C", m2, n, &c_b2, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b2, &
	    work[1], &c__1, (ftnlen)1);

#line 282 "cunbdb6.f"
    cgemv_("N", m1, n, &c_b1, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b2, &
	    x1[1], incx1, (ftnlen)1);
#line 284 "cunbdb6.f"
    cgemv_("N", m2, n, &c_b1, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b2, &
	    x2[1], incx2, (ftnlen)1);

#line 287 "cunbdb6.f"
    scl1 = 0.;
#line 288 "cunbdb6.f"
    ssq1 = 1.;
#line 289 "cunbdb6.f"
    classq_(m1, &x1[1], incx1, &scl1, &ssq1);
#line 290 "cunbdb6.f"
    scl2 = 0.;
#line 291 "cunbdb6.f"
    ssq2 = 1.;
#line 292 "cunbdb6.f"
    classq_(m1, &x1[1], incx1, &scl1, &ssq1);
/* Computing 2nd power */
#line 293 "cunbdb6.f"
    d__1 = scl1;
/* Computing 2nd power */
#line 293 "cunbdb6.f"
    d__2 = scl2;
#line 293 "cunbdb6.f"
    normsq2 = d__1 * d__1 * ssq1 + d__2 * d__2 * ssq2;

/*     If second projection is sufficiently large in norm, then do */
/*     nothing more. Alternatively, if it shrunk significantly, then */
/*     truncate it to zero. */

#line 299 "cunbdb6.f"
    if (normsq2 < normsq1 * .01) {
#line 300 "cunbdb6.f"
	i__1 = *m1;
#line 300 "cunbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 301 "cunbdb6.f"
	    i__2 = i__;
#line 301 "cunbdb6.f"
	    x1[i__2].r = 0., x1[i__2].i = 0.;
#line 302 "cunbdb6.f"
	}
#line 303 "cunbdb6.f"
	i__1 = *m2;
#line 303 "cunbdb6.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 304 "cunbdb6.f"
	    i__2 = i__;
#line 304 "cunbdb6.f"
	    x2[i__2].r = 0., x2[i__2].i = 0.;
#line 305 "cunbdb6.f"
	}
#line 306 "cunbdb6.f"
    }

#line 308 "cunbdb6.f"
    return 0;

/*     End of CUNBDB6 */

} /* cunbdb6_ */

