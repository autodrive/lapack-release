#line 1 "dorcsd2by1.f"
/* dorcsd2by1.f -- translated by f2c (version 20100827).
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

#line 1 "dorcsd2by1.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__1 = 1;
static logical c_false = FALSE_;

/* > \brief \b DORCSD2BY1 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORCSD2BY1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorcsd2
by1.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorcsd2
by1.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorcsd2
by1.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, */
/*                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T, */
/*                              LDV1T, WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU1, JOBU2, JOBV1T */
/*       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, */
/*      $                   M, P, Q */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   THETA(*) */
/*       DOUBLE PRECISION   U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), */
/*      $                   X11(LDX11,*), X21(LDX21,*) */
/*       INTEGER            IWORK(*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > Purpose: */
/* > ======== */
/* > */
/* > DORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with */
/* > orthonormal columns that has been partitioned into a 2-by-1 block */
/* > structure: */
/* > */
/* >                                [  I  0  0 ] */
/* >                                [  0  C  0 ] */
/* >          [ X11 ]   [ U1 |    ] [  0  0  0 ] */
/* >      X = [-----] = [---------] [----------] V1**T . */
/* >          [ X21 ]   [    | U2 ] [  0  0  0 ] */
/* >                                [  0  S  0 ] */
/* >                                [  0  0  I ] */
/* > */
/* > X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P, */
/* > (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R */
/* > nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which */
/* > R = MIN(P,M-P,Q,M-Q). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU1 */
/* > \verbatim */
/* >          JOBU1 is CHARACTER */
/* >          = 'Y':      U1 is computed; */
/* >          otherwise:  U1 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU2 */
/* > \verbatim */
/* >          JOBU2 is CHARACTER */
/* >          = 'Y':      U2 is computed; */
/* >          otherwise:  U2 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV1T */
/* > \verbatim */
/* >          JOBV1T is CHARACTER */
/* >          = 'Y':      V1T is computed; */
/* >          otherwise:  V1T is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows in X. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >          The number of rows in X11. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is INTEGER */
/* >          The number of columns in X11 and X21. 0 <= Q <= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* >          X11 is DOUBLE PRECISION array, dimension (LDX11,Q) */
/* >          On entry, part of the orthogonal matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* >          LDX11 is INTEGER */
/* >          The leading dimension of X11. LDX11 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* >          X21 is DOUBLE PRECISION array, dimension (LDX21,Q) */
/* >          On entry, part of the orthogonal matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* >          LDX21 is INTEGER */
/* >          The leading dimension of X21. LDX21 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* >          THETA is DOUBLE PRECISION array, dimension (R), in which R = */
/* >          MIN(P,M-P,Q,M-Q). */
/* >          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* >          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* >          U1 is DOUBLE PRECISION array, dimension (P) */
/* >          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU1 */
/* > \verbatim */
/* >          LDU1 is INTEGER */
/* >          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >= */
/* >          MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* >          U2 is DOUBLE PRECISION array, dimension (M-P) */
/* >          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal */
/* >          matrix U2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* >          LDU2 is INTEGER */
/* >          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >= */
/* >          MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] V1T */
/* > \verbatim */
/* >          V1T is DOUBLE PRECISION array, dimension (Q) */
/* >          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal */
/* >          matrix V1**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV1T */
/* > \verbatim */
/* >          LDV1T is INTEGER */
/* >          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >= */
/* >          MAX(1,Q). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* >          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1), */
/* >          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* >          define the matrix in intermediate bidiagonal-block form */
/* >          remaining after nonconvergence. INFO specifies the number */
/* >          of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the work array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  DBBCSD did not converge. See the description of WORK */
/* >                above for details. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date July 2012 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorcsd2by1_(char *jobu1, char *jobu2, char *jobv1t, 
	integer *m, integer *p, integer *q, doublereal *x11, integer *ldx11, 
	doublereal *x21, integer *ldx21, doublereal *theta, doublereal *u1, 
	integer *ldu1, doublereal *u2, integer *ldu2, doublereal *v1t, 
	integer *ldv1t, doublereal *work, integer *lwork, integer *iwork, 
	integer *info, ftnlen jobu1_len, ftnlen jobu2_len, ftnlen jobv1t_len)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, 
	    x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6, i__7, i__8, i__9;

    /* Local variables */
    static integer lworkmin, lworkopt, i__, j, r__, childinfo, lorglqmin, 
	    lorgqrmin, lorglqopt, lorgqropt, ib11d, ib11e, ib12d, ib12e, 
	    ib21d, ib21e, ib22d, ib22e, iphi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer itaup1, itaup2, itauq1;
    static logical wantu1, wantu2;
    extern /* Subroutine */ int dbbcsd_();
    static integer ibbcsd, lbbcsd, iorbdb, lorbdb;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlapmr_(logical *, integer *, 
	    integer *, doublereal *, integer *, integer *), dlapmt_(logical *,
	     integer *, integer *, doublereal *, integer *, integer *), 
	    dorglq_();
    static integer iorglq;
    extern /* Subroutine */ int dorgqr_();
    static integer lorglq, iorgqr, lorgqr;
    extern /* Subroutine */ int dorbdb1_(), dorbdb2_(), dorbdb3_(), dorbdb4_()
	    ;
    static logical lquery, wantv1t;


/*  -- LAPACK computational routine (3.5.0) -- */
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

#line 285 "dorcsd2by1.f"
    /* Parameter adjustments */
#line 285 "dorcsd2by1.f"
    x11_dim1 = *ldx11;
#line 285 "dorcsd2by1.f"
    x11_offset = 1 + x11_dim1;
#line 285 "dorcsd2by1.f"
    x11 -= x11_offset;
#line 285 "dorcsd2by1.f"
    x21_dim1 = *ldx21;
#line 285 "dorcsd2by1.f"
    x21_offset = 1 + x21_dim1;
#line 285 "dorcsd2by1.f"
    x21 -= x21_offset;
#line 285 "dorcsd2by1.f"
    --theta;
#line 285 "dorcsd2by1.f"
    u1_dim1 = *ldu1;
#line 285 "dorcsd2by1.f"
    u1_offset = 1 + u1_dim1;
#line 285 "dorcsd2by1.f"
    u1 -= u1_offset;
#line 285 "dorcsd2by1.f"
    u2_dim1 = *ldu2;
#line 285 "dorcsd2by1.f"
    u2_offset = 1 + u2_dim1;
#line 285 "dorcsd2by1.f"
    u2 -= u2_offset;
#line 285 "dorcsd2by1.f"
    v1t_dim1 = *ldv1t;
#line 285 "dorcsd2by1.f"
    v1t_offset = 1 + v1t_dim1;
#line 285 "dorcsd2by1.f"
    v1t -= v1t_offset;
#line 285 "dorcsd2by1.f"
    --work;
#line 285 "dorcsd2by1.f"
    --iwork;
#line 285 "dorcsd2by1.f"

#line 285 "dorcsd2by1.f"
    /* Function Body */
#line 285 "dorcsd2by1.f"
    *info = 0;
#line 286 "dorcsd2by1.f"
    wantu1 = lsame_(jobu1, "Y", (ftnlen)1, (ftnlen)1);
#line 287 "dorcsd2by1.f"
    wantu2 = lsame_(jobu2, "Y", (ftnlen)1, (ftnlen)1);
#line 288 "dorcsd2by1.f"
    wantv1t = lsame_(jobv1t, "Y", (ftnlen)1, (ftnlen)1);
#line 289 "dorcsd2by1.f"
    lquery = *lwork == -1;

#line 291 "dorcsd2by1.f"
    if (*m < 0) {
#line 292 "dorcsd2by1.f"
	*info = -4;
#line 293 "dorcsd2by1.f"
    } else if (*p < 0 || *p > *m) {
#line 294 "dorcsd2by1.f"
	*info = -5;
#line 295 "dorcsd2by1.f"
    } else if (*q < 0 || *q > *m) {
#line 296 "dorcsd2by1.f"
	*info = -6;
#line 297 "dorcsd2by1.f"
    } else if (*ldx11 < max(1,*p)) {
#line 298 "dorcsd2by1.f"
	*info = -8;
#line 299 "dorcsd2by1.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 299 "dorcsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 299 "dorcsd2by1.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 300 "dorcsd2by1.f"
	    *info = -10;
#line 301 "dorcsd2by1.f"
	} else if (wantu1 && *ldu1 < *p) {
#line 302 "dorcsd2by1.f"
	    *info = -13;
#line 303 "dorcsd2by1.f"
	} else if (wantu2 && *ldu2 < *m - *p) {
#line 304 "dorcsd2by1.f"
	    *info = -15;
#line 305 "dorcsd2by1.f"
	} else if (wantv1t && *ldv1t < *q) {
#line 306 "dorcsd2by1.f"
	    *info = -17;
#line 307 "dorcsd2by1.f"
	}
#line 307 "dorcsd2by1.f"
    }

/* Computing MIN */
#line 309 "dorcsd2by1.f"
    i__1 = *p, i__2 = *m - *p, i__1 = min(i__1,i__2), i__1 = min(i__1,*q), 
	    i__2 = *m - *q;
#line 309 "dorcsd2by1.f"
    r__ = min(i__1,i__2);

/*     Compute workspace */

/*       WORK layout: */
/*     |-------------------------------------------------------| */
/*     | LWORKOPT (1)                                          | */
/*     |-------------------------------------------------------| */
/*     | PHI (MAX(1,R-1))                                      | */
/*     |-------------------------------------------------------| */
/*     | TAUP1 (MAX(1,P))                        | B11D (R)    | */
/*     | TAUP2 (MAX(1,M-P))                      | B11E (R-1)  | */
/*     | TAUQ1 (MAX(1,Q))                        | B12D (R)    | */
/*     |-----------------------------------------| B12E (R-1)  | */
/*     | DORBDB WORK | DORGQR WORK | DORGLQ WORK | B21D (R)    | */
/*     |             |             |             | B21E (R-1)  | */
/*     |             |             |             | B22D (R)    | */
/*     |             |             |             | B22E (R-1)  | */
/*     |             |             |             | DBBCSD WORK | */
/*     |-------------------------------------------------------| */

#line 330 "dorcsd2by1.f"
    if (*info == 0) {
#line 331 "dorcsd2by1.f"
	iphi = 2;
/* Computing MAX */
#line 332 "dorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 332 "dorcsd2by1.f"
	ib11d = iphi + max(i__1,i__2);
#line 333 "dorcsd2by1.f"
	ib11e = ib11d + max(1,r__);
/* Computing MAX */
#line 334 "dorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 334 "dorcsd2by1.f"
	ib12d = ib11e + max(i__1,i__2);
#line 335 "dorcsd2by1.f"
	ib12e = ib12d + max(1,r__);
/* Computing MAX */
#line 336 "dorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 336 "dorcsd2by1.f"
	ib21d = ib12e + max(i__1,i__2);
#line 337 "dorcsd2by1.f"
	ib21e = ib21d + max(1,r__);
/* Computing MAX */
#line 338 "dorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 338 "dorcsd2by1.f"
	ib22d = ib21e + max(i__1,i__2);
#line 339 "dorcsd2by1.f"
	ib22e = ib22d + max(1,r__);
/* Computing MAX */
#line 340 "dorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 340 "dorcsd2by1.f"
	ibbcsd = ib22e + max(i__1,i__2);
/* Computing MAX */
#line 341 "dorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 341 "dorcsd2by1.f"
	itaup1 = iphi + max(i__1,i__2);
#line 342 "dorcsd2by1.f"
	itaup2 = itaup1 + max(1,*p);
/* Computing MAX */
#line 343 "dorcsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 343 "dorcsd2by1.f"
	itauq1 = itaup2 + max(i__1,i__2);
#line 344 "dorcsd2by1.f"
	iorbdb = itauq1 + max(1,*q);
#line 345 "dorcsd2by1.f"
	iorgqr = itauq1 + max(1,*q);
#line 346 "dorcsd2by1.f"
	iorglq = itauq1 + max(1,*q);
#line 347 "dorcsd2by1.f"
	if (r__ == *q) {
#line 348 "dorcsd2by1.f"
	    dorbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 350 "dorcsd2by1.f"
	    lorbdb = (integer) work[1];
#line 351 "dorcsd2by1.f"
	    if (*p >= *m - *p) {
#line 352 "dorcsd2by1.f"
		dorgqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 354 "dorcsd2by1.f"
		lorgqrmin = max(1,*p);
#line 355 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 356 "dorcsd2by1.f"
	    } else {
#line 357 "dorcsd2by1.f"
		i__1 = *m - *p;
#line 357 "dorcsd2by1.f"
		i__2 = *m - *p;
#line 357 "dorcsd2by1.f"
		dorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 359 "dorcsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 359 "dorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 360 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 361 "dorcsd2by1.f"
	    }
/* Computing MAX */
#line 362 "dorcsd2by1.f"
	    i__2 = 0, i__3 = *q - 1;
#line 362 "dorcsd2by1.f"
	    i__1 = max(i__2,i__3);
/* Computing MAX */
#line 362 "dorcsd2by1.f"
	    i__5 = 0, i__6 = *q - 1;
#line 362 "dorcsd2by1.f"
	    i__4 = max(i__5,i__6);
/* Computing MAX */
#line 362 "dorcsd2by1.f"
	    i__8 = 0, i__9 = *q - 1;
#line 362 "dorcsd2by1.f"
	    i__7 = max(i__8,i__9);
#line 362 "dorcsd2by1.f"
	    dorglq_(&i__1, &i__4, &i__7, &v1t[v1t_offset], ldv1t, &c__0, &
		    work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 364 "dorcsd2by1.f"
	    i__1 = 1, i__2 = *q - 1;
#line 364 "dorcsd2by1.f"
	    lorglqmin = max(i__1,i__2);
#line 365 "dorcsd2by1.f"
	    lorglqopt = (integer) work[1];
#line 366 "dorcsd2by1.f"
	    dbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &c__0,
		     &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		    v1t_offset], ldv1t, &c__0, &c__1, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &work[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 369 "dorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 370 "dorcsd2by1.f"
	} else if (r__ == *p) {
#line 371 "dorcsd2by1.f"
	    dorbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 373 "dorcsd2by1.f"
	    lorbdb = (integer) work[1];
#line 374 "dorcsd2by1.f"
	    if (*p - 1 >= *m - *p) {
#line 375 "dorcsd2by1.f"
		i__1 = *p - 1;
#line 375 "dorcsd2by1.f"
		i__2 = *p - 1;
#line 375 "dorcsd2by1.f"
		i__3 = *p - 1;
#line 375 "dorcsd2by1.f"
		dorgqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 377 "dorcsd2by1.f"
		i__1 = 1, i__2 = *p - 1;
#line 377 "dorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 378 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 379 "dorcsd2by1.f"
	    } else {
#line 380 "dorcsd2by1.f"
		i__1 = *m - *p;
#line 380 "dorcsd2by1.f"
		i__2 = *m - *p;
#line 380 "dorcsd2by1.f"
		dorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 382 "dorcsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 382 "dorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 383 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 384 "dorcsd2by1.f"
	    }
#line 385 "dorcsd2by1.f"
	    dorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 387 "dorcsd2by1.f"
	    lorglqmin = max(1,*q);
#line 388 "dorcsd2by1.f"
	    lorglqopt = (integer) work[1];
#line 389 "dorcsd2by1.f"
	    dbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &c__0,
		     &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		    ldu1, &u2[u2_offset], ldu2, &c__0, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &work[1], &c_n1, &childinfo, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 392 "dorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 393 "dorcsd2by1.f"
	} else if (r__ == *m - *p) {
#line 394 "dorcsd2by1.f"
	    dorbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 396 "dorcsd2by1.f"
	    lorbdb = (integer) work[1];
#line 397 "dorcsd2by1.f"
	    if (*p >= *m - *p - 1) {
#line 398 "dorcsd2by1.f"
		dorgqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 400 "dorcsd2by1.f"
		lorgqrmin = max(1,*p);
#line 401 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 402 "dorcsd2by1.f"
	    } else {
#line 403 "dorcsd2by1.f"
		i__1 = *m - *p - 1;
#line 403 "dorcsd2by1.f"
		i__2 = *m - *p - 1;
#line 403 "dorcsd2by1.f"
		i__3 = *m - *p - 1;
#line 403 "dorcsd2by1.f"
		dorgqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 405 "dorcsd2by1.f"
		i__1 = 1, i__2 = *m - *p - 1;
#line 405 "dorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 406 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 407 "dorcsd2by1.f"
	    }
#line 408 "dorcsd2by1.f"
	    dorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 410 "dorcsd2by1.f"
	    lorglqmin = max(1,*q);
#line 411 "dorcsd2by1.f"
	    lorglqopt = (integer) work[1];
#line 412 "dorcsd2by1.f"
	    i__1 = *m - *q;
#line 412 "dorcsd2by1.f"
	    i__2 = *m - *p;
#line 412 "dorcsd2by1.f"
	    dbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1]
		    , &c__0, &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		    u2_offset], ldu2, &u1[u1_offset], ldu1, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &work[1], &c_n1, 
		    &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 416 "dorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 417 "dorcsd2by1.f"
	} else {
#line 418 "dorcsd2by1.f"
	    dorbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &c__0, &
		    work[1], &c_n1, &childinfo);
#line 420 "dorcsd2by1.f"
	    lorbdb = *m + (integer) work[1];
#line 421 "dorcsd2by1.f"
	    if (*p >= *m - *p) {
#line 422 "dorcsd2by1.f"
		i__1 = *m - *q;
#line 422 "dorcsd2by1.f"
		dorgqr_(p, p, &i__1, &u1[u1_offset], ldu1, &c__0, &work[1], &
			c_n1, &childinfo);
#line 424 "dorcsd2by1.f"
		lorgqrmin = max(1,*p);
#line 425 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 426 "dorcsd2by1.f"
	    } else {
#line 427 "dorcsd2by1.f"
		i__1 = *m - *p;
#line 427 "dorcsd2by1.f"
		i__2 = *m - *p;
#line 427 "dorcsd2by1.f"
		i__3 = *m - *q;
#line 427 "dorcsd2by1.f"
		dorgqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &c__0, &
			work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 429 "dorcsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 429 "dorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 430 "dorcsd2by1.f"
		lorgqropt = (integer) work[1];
#line 431 "dorcsd2by1.f"
	    }
#line 432 "dorcsd2by1.f"
	    dorglq_(q, q, q, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &c_n1, 
		    &childinfo);
#line 434 "dorcsd2by1.f"
	    lorglqmin = max(1,*q);
#line 435 "dorcsd2by1.f"
	    lorglqopt = (integer) work[1];
#line 436 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 436 "dorcsd2by1.f"
	    i__2 = *m - *q;
#line 436 "dorcsd2by1.f"
	    dbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1]
		    , &c__0, &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		    c__0, &c__1, &v1t[v1t_offset], ldv1t, &c__0, &c__0, &c__0,
		     &c__0, &c__0, &c__0, &c__0, &c__0, &work[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 440 "dorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 441 "dorcsd2by1.f"
	}
/* Computing MAX */
#line 442 "dorcsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqrmin - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqmin - 1, i__1 = max(i__1,
		i__2), i__2 = ibbcsd + lbbcsd - 1;
#line 442 "dorcsd2by1.f"
	lworkmin = max(i__1,i__2);
/* Computing MAX */
#line 446 "dorcsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqropt - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqopt - 1, i__1 = max(i__1,
		i__2), i__2 = ibbcsd + lbbcsd - 1;
#line 446 "dorcsd2by1.f"
	lworkopt = max(i__1,i__2);
#line 450 "dorcsd2by1.f"
	work[1] = (doublereal) lworkopt;
#line 451 "dorcsd2by1.f"
	if (*lwork < lworkmin && ! lquery) {
#line 452 "dorcsd2by1.f"
	    *info = -19;
#line 453 "dorcsd2by1.f"
	}
#line 454 "dorcsd2by1.f"
    }
#line 455 "dorcsd2by1.f"
    if (*info != 0) {
#line 456 "dorcsd2by1.f"
	i__1 = -(*info);
#line 456 "dorcsd2by1.f"
	xerbla_("DORCSD2BY1", &i__1, (ftnlen)10);
#line 457 "dorcsd2by1.f"
	return 0;
#line 458 "dorcsd2by1.f"
    } else if (lquery) {
#line 459 "dorcsd2by1.f"
	return 0;
#line 460 "dorcsd2by1.f"
    }
#line 461 "dorcsd2by1.f"
    lorgqr = *lwork - iorgqr + 1;
#line 462 "dorcsd2by1.f"
    lorglq = *lwork - iorglq + 1;

/*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q, */
/*     in which R = MIN(P,M-P,Q,M-Q) */

#line 467 "dorcsd2by1.f"
    if (r__ == *q) {

/*        Case 1: R = Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 473 "dorcsd2by1.f"
	dorbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 479 "dorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 480 "dorcsd2by1.f"
	    dlacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 481 "dorcsd2by1.f"
	    dorgqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 483 "dorcsd2by1.f"
	}
#line 484 "dorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 485 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 485 "dorcsd2by1.f"
	    dlacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 486 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 486 "dorcsd2by1.f"
	    i__2 = *m - *p;
#line 486 "dorcsd2by1.f"
	    dorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 488 "dorcsd2by1.f"
	}
#line 489 "dorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 490 "dorcsd2by1.f"
	    v1t[v1t_dim1 + 1] = 1.;
#line 491 "dorcsd2by1.f"
	    i__1 = *q;
#line 491 "dorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 492 "dorcsd2by1.f"
		v1t[j * v1t_dim1 + 1] = 0.;
#line 493 "dorcsd2by1.f"
		v1t[j + v1t_dim1] = 0.;
#line 494 "dorcsd2by1.f"
	    }
#line 495 "dorcsd2by1.f"
	    i__1 = *q - 1;
#line 495 "dorcsd2by1.f"
	    i__2 = *q - 1;
#line 495 "dorcsd2by1.f"
	    dlacpy_("U", &i__1, &i__2, &x21[(x21_dim1 << 1) + 1], ldx21, &v1t[
		    (v1t_dim1 << 1) + 2], ldv1t, (ftnlen)1);
#line 497 "dorcsd2by1.f"
	    i__1 = *q - 1;
#line 497 "dorcsd2by1.f"
	    i__2 = *q - 1;
#line 497 "dorcsd2by1.f"
	    i__3 = *q - 1;
#line 497 "dorcsd2by1.f"
	    dorglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &
		    work[itauq1], &work[iorglq], &lorglq, &childinfo);
#line 499 "dorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 503 "dorcsd2by1.f"
	dbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &work[
		iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		v1t_offset], ldv1t, &c__0, &c__1, &work[ib11d], &work[ib11e], 
		&work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e], &work[
		ib22d], &work[ib22e], &work[ibbcsd], &lbbcsd, &childinfo, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Permute rows and columns to place zero submatrices in */
/*        preferred positions */

#line 513 "dorcsd2by1.f"
	if (*q > 0 && wantu2) {
#line 514 "dorcsd2by1.f"
	    i__1 = *q;
#line 514 "dorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 515 "dorcsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 516 "dorcsd2by1.f"
	    }
#line 517 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 517 "dorcsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 518 "dorcsd2by1.f"
		iwork[i__] = i__ - *q;
#line 519 "dorcsd2by1.f"
	    }
#line 520 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 520 "dorcsd2by1.f"
	    i__2 = *m - *p;
#line 520 "dorcsd2by1.f"
	    dlapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 521 "dorcsd2by1.f"
	}
#line 522 "dorcsd2by1.f"
    } else if (r__ == *p) {

/*        Case 2: R = P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 528 "dorcsd2by1.f"
	dorbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 534 "dorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 535 "dorcsd2by1.f"
	    u1[u1_dim1 + 1] = 1.;
#line 536 "dorcsd2by1.f"
	    i__1 = *p;
#line 536 "dorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 537 "dorcsd2by1.f"
		u1[j * u1_dim1 + 1] = 0.;
#line 538 "dorcsd2by1.f"
		u1[j + u1_dim1] = 0.;
#line 539 "dorcsd2by1.f"
	    }
#line 540 "dorcsd2by1.f"
	    i__1 = *p - 1;
#line 540 "dorcsd2by1.f"
	    i__2 = *p - 1;
#line 540 "dorcsd2by1.f"
	    dlacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 541 "dorcsd2by1.f"
	    i__1 = *p - 1;
#line 541 "dorcsd2by1.f"
	    i__2 = *p - 1;
#line 541 "dorcsd2by1.f"
	    i__3 = *p - 1;
#line 541 "dorcsd2by1.f"
	    dorgqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &work[
		    itaup1], &work[iorgqr], &lorgqr, &childinfo);
#line 543 "dorcsd2by1.f"
	}
#line 544 "dorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 545 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 545 "dorcsd2by1.f"
	    dlacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 546 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 546 "dorcsd2by1.f"
	    i__2 = *m - *p;
#line 546 "dorcsd2by1.f"
	    dorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 548 "dorcsd2by1.f"
	}
#line 549 "dorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 550 "dorcsd2by1.f"
	    dlacpy_("U", p, q, &x11[x11_offset], ldx11, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 551 "dorcsd2by1.f"
	    dorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 553 "dorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 557 "dorcsd2by1.f"
	dbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &work[
		iphi], &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		ldu1, &u2[u2_offset], ldu2, &work[ib11d], &work[ib11e], &work[
		ib12d], &work[ib12e], &work[ib21d], &work[ib21e], &work[ib22d]
		, &work[ib22e], &work[ibbcsd], &lbbcsd, &childinfo, (ftnlen)1,
		 (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 567 "dorcsd2by1.f"
	if (*q > 0 && wantu2) {
#line 568 "dorcsd2by1.f"
	    i__1 = *q;
#line 568 "dorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 569 "dorcsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 570 "dorcsd2by1.f"
	    }
#line 571 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 571 "dorcsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 572 "dorcsd2by1.f"
		iwork[i__] = i__ - *q;
#line 573 "dorcsd2by1.f"
	    }
#line 574 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 574 "dorcsd2by1.f"
	    i__2 = *m - *p;
#line 574 "dorcsd2by1.f"
	    dlapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 575 "dorcsd2by1.f"
	}
#line 576 "dorcsd2by1.f"
    } else if (r__ == *m - *p) {

/*        Case 3: R = M-P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 582 "dorcsd2by1.f"
	dorbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 588 "dorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 589 "dorcsd2by1.f"
	    dlacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 590 "dorcsd2by1.f"
	    dorgqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 592 "dorcsd2by1.f"
	}
#line 593 "dorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 594 "dorcsd2by1.f"
	    u2[u2_dim1 + 1] = 1.;
#line 595 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 595 "dorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 596 "dorcsd2by1.f"
		u2[j * u2_dim1 + 1] = 0.;
#line 597 "dorcsd2by1.f"
		u2[j + u2_dim1] = 0.;
#line 598 "dorcsd2by1.f"
	    }
#line 599 "dorcsd2by1.f"
	    i__1 = *m - *p - 1;
#line 599 "dorcsd2by1.f"
	    i__2 = *m - *p - 1;
#line 599 "dorcsd2by1.f"
	    dlacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 601 "dorcsd2by1.f"
	    i__1 = *m - *p - 1;
#line 601 "dorcsd2by1.f"
	    i__2 = *m - *p - 1;
#line 601 "dorcsd2by1.f"
	    i__3 = *m - *p - 1;
#line 601 "dorcsd2by1.f"
	    dorgqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &work[
		    itaup2], &work[iorgqr], &lorgqr, &childinfo);
#line 603 "dorcsd2by1.f"
	}
#line 604 "dorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 605 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 605 "dorcsd2by1.f"
	    dlacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 606 "dorcsd2by1.f"
	    dorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 608 "dorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 612 "dorcsd2by1.f"
	i__1 = *m - *q;
#line 612 "dorcsd2by1.f"
	i__2 = *m - *p;
#line 612 "dorcsd2by1.f"
	dbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1], &
		work[iphi], &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		u2_offset], ldu2, &u1[u1_offset], ldu1, &work[ib11d], &work[
		ib11e], &work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e]
		, &work[ib22d], &work[ib22e], &work[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 622 "dorcsd2by1.f"
	if (*q > r__) {
#line 623 "dorcsd2by1.f"
	    i__1 = r__;
#line 623 "dorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 624 "dorcsd2by1.f"
		iwork[i__] = *q - r__ + i__;
#line 625 "dorcsd2by1.f"
	    }
#line 626 "dorcsd2by1.f"
	    i__1 = *q;
#line 626 "dorcsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 627 "dorcsd2by1.f"
		iwork[i__] = i__ - r__;
#line 628 "dorcsd2by1.f"
	    }
#line 629 "dorcsd2by1.f"
	    if (wantu1) {
#line 630 "dorcsd2by1.f"
		dlapmt_(&c_false, p, q, &u1[u1_offset], ldu1, &iwork[1]);
#line 631 "dorcsd2by1.f"
	    }
#line 632 "dorcsd2by1.f"
	    if (wantv1t) {
#line 633 "dorcsd2by1.f"
		dlapmr_(&c_false, q, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 634 "dorcsd2by1.f"
	    }
#line 635 "dorcsd2by1.f"
	}
#line 636 "dorcsd2by1.f"
    } else {

/*        Case 4: R = M-Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 642 "dorcsd2by1.f"
	i__1 = lorbdb - *m;
#line 642 "dorcsd2by1.f"
	dorbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &work[iorbdb + *m], &i__1, &childinfo)
		;

/*        Accumulate Householder reflectors */

#line 649 "dorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 650 "dorcsd2by1.f"
	    dcopy_(p, &work[iorbdb], &c__1, &u1[u1_offset], &c__1);
#line 651 "dorcsd2by1.f"
	    i__1 = *p;
#line 651 "dorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 652 "dorcsd2by1.f"
		u1[j * u1_dim1 + 1] = 0.;
#line 653 "dorcsd2by1.f"
	    }
#line 654 "dorcsd2by1.f"
	    i__1 = *p - 1;
#line 654 "dorcsd2by1.f"
	    i__2 = *m - *q - 1;
#line 654 "dorcsd2by1.f"
	    dlacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 656 "dorcsd2by1.f"
	    i__1 = *m - *q;
#line 656 "dorcsd2by1.f"
	    dorgqr_(p, p, &i__1, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 658 "dorcsd2by1.f"
	}
#line 659 "dorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 660 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 660 "dorcsd2by1.f"
	    dcopy_(&i__1, &work[iorbdb + *p], &c__1, &u2[u2_offset], &c__1);
#line 661 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 661 "dorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 662 "dorcsd2by1.f"
		u2[j * u2_dim1 + 1] = 0.;
#line 663 "dorcsd2by1.f"
	    }
#line 664 "dorcsd2by1.f"
	    i__1 = *m - *p - 1;
#line 664 "dorcsd2by1.f"
	    i__2 = *m - *q - 1;
#line 664 "dorcsd2by1.f"
	    dlacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 666 "dorcsd2by1.f"
	    i__1 = *m - *p;
#line 666 "dorcsd2by1.f"
	    i__2 = *m - *p;
#line 666 "dorcsd2by1.f"
	    i__3 = *m - *q;
#line 666 "dorcsd2by1.f"
	    dorgqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &work[itaup2], 
		    &work[iorgqr], &lorgqr, &childinfo);
#line 668 "dorcsd2by1.f"
	}
#line 669 "dorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 670 "dorcsd2by1.f"
	    i__1 = *m - *q;
#line 670 "dorcsd2by1.f"
	    dlacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 671 "dorcsd2by1.f"
	    i__1 = *p - (*m - *q);
#line 671 "dorcsd2by1.f"
	    i__2 = *q - (*m - *q);
#line 671 "dorcsd2by1.f"
	    dlacpy_("U", &i__1, &i__2, &x11[*m - *q + 1 + (*m - *q + 1) * 
		    x11_dim1], ldx11, &v1t[*m - *q + 1 + (*m - *q + 1) * 
		    v1t_dim1], ldv1t, (ftnlen)1);
#line 673 "dorcsd2by1.f"
	    i__1 = -(*p) + *q;
#line 673 "dorcsd2by1.f"
	    i__2 = *q - *p;
#line 673 "dorcsd2by1.f"
	    dlacpy_("U", &i__1, &i__2, &x21[*m - *q + 1 + (*p + 1) * x21_dim1]
		    , ldx21, &v1t[*p + 1 + (*p + 1) * v1t_dim1], ldv1t, (
		    ftnlen)1);
#line 675 "dorcsd2by1.f"
	    dorglq_(q, q, q, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 677 "dorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 681 "dorcsd2by1.f"
	i__1 = *m - *p;
#line 681 "dorcsd2by1.f"
	i__2 = *m - *q;
#line 681 "dorcsd2by1.f"
	dbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1], &
		work[iphi], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &c__0,
		 &c__1, &v1t[v1t_offset], ldv1t, &work[ib11d], &work[ib11e], &
		work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e], &work[
		ib22d], &work[ib22e], &work[ibbcsd], &lbbcsd, &childinfo, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 691 "dorcsd2by1.f"
	if (*p > r__) {
#line 692 "dorcsd2by1.f"
	    i__1 = r__;
#line 692 "dorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 693 "dorcsd2by1.f"
		iwork[i__] = *p - r__ + i__;
#line 694 "dorcsd2by1.f"
	    }
#line 695 "dorcsd2by1.f"
	    i__1 = *p;
#line 695 "dorcsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 696 "dorcsd2by1.f"
		iwork[i__] = i__ - r__;
#line 697 "dorcsd2by1.f"
	    }
#line 698 "dorcsd2by1.f"
	    if (wantu1) {
#line 699 "dorcsd2by1.f"
		dlapmt_(&c_false, p, p, &u1[u1_offset], ldu1, &iwork[1]);
#line 700 "dorcsd2by1.f"
	    }
#line 701 "dorcsd2by1.f"
	    if (wantv1t) {
#line 702 "dorcsd2by1.f"
		dlapmr_(&c_false, p, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 703 "dorcsd2by1.f"
	    }
#line 704 "dorcsd2by1.f"
	}
#line 705 "dorcsd2by1.f"
    }

#line 707 "dorcsd2by1.f"
    return 0;

/*     End of DORCSD2BY1 */

} /* dorcsd2by1_ */

