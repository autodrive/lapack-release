#line 1 "zuncsd2by1.f"
/* zuncsd2by1.f -- translated by f2c (version 20100827).
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

#line 1 "zuncsd2by1.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__1 = 1;
static logical c_false = FALSE_;

/* > \brief \b ZUNCSD2BY1 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNCSD2BY1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zuncsd2
by1.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zuncsd2
by1.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zuncsd2
by1.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, */
/*                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T, */
/*                              LDV1T, WORK, LWORK, RWORK, LRWORK, IWORK, */
/*                              INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU1, JOBU2, JOBV1T */
/*       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, */
/*      $                   M, P, Q */
/*       INTEGER            LRWORK, LRWORKMIN, LRWORKOPT */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK(*) */
/*       DOUBLE PRECISION   THETA(*) */
/*       COMPLEX*16         U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), */
/*      $                   X11(LDX11,*), X21(LDX21,*) */
/*       INTEGER            IWORK(*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > ZUNCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with */
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
/* > X11 is P-by-Q. The unitary matrices U1, U2, and V1 are P-by-P, */
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
/* >          X11 is COMPLEX*16 array, dimension (LDX11,Q) */
/* >          On entry, part of the unitary matrix whose CSD is desired. */
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
/* >          X21 is COMPLEX*16 array, dimension (LDX21,Q) */
/* >          On entry, part of the unitary matrix whose CSD is desired. */
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
/* >          U1 is COMPLEX*16 array, dimension (P) */
/* >          If JOBU1 = 'Y', U1 contains the P-by-P unitary matrix U1. */
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
/* >          U2 is COMPLEX*16 array, dimension (M-P) */
/* >          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) unitary */
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
/* >          V1T is COMPLEX*16 array, dimension (Q) */
/* >          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* >          If INFO > 0 on exit, RWORK(2:R) contains the values PHI(1), */
/* >          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* >          define the matrix in intermediate bidiagonal-block form */
/* >          remaining after nonconvergence. INFO specifies the number */
/* >          of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of the array RWORK. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the RWORK array, returns */
/* >          this value as the first entry of the work array, and no error */
/* >          message related to LRWORK is issued by XERBLA. */
/* > \endverbatim */

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
/* >          > 0:  ZBBCSD did not converge. See the description of WORK */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zuncsd2by1_(char *jobu1, char *jobu2, char *jobv1t, 
	integer *m, integer *p, integer *q, doublecomplex *x11, integer *
	ldx11, doublecomplex *x21, integer *ldx21, doublereal *theta, 
	doublecomplex *u1, integer *ldu1, doublecomplex *u2, integer *ldu2, 
	doublecomplex *v1t, integer *ldv1t, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *
	info, ftnlen jobu1_len, ftnlen jobu2_len, ftnlen jobv1t_len)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, 
	    x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6, i__7, i__8, i__9;

    /* Local variables */
    static integer lworkmin, lworkopt, i__, j, r__, childinfo, lorglqmin, 
	    lorgqrmin, lorglqopt, lrworkmin, lorgqropt, lrworkopt, ib11d, 
	    ib11e, ib12d, ib12e, ib21d, ib21e, ib22d, ib22e, iphi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer itaup1, itaup2, itauq1;
    static logical wantu1, wantu2;
    static integer ibbcsd, lbbcsd, iorbdb, lorbdb;
    extern /* Subroutine */ int zbbcsd_(), xerbla_(char *, integer *, ftnlen);
    static integer iorglq, lorglq;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer iorgqr;
    extern /* Subroutine */ int zlapmr_(logical *, integer *, integer *, 
	    doublecomplex *, integer *, integer *);
    static integer lorgqr;
    extern /* Subroutine */ int zlapmt_(logical *, integer *, integer *, 
	    doublecomplex *, integer *, integer *);
    static logical lquery;
    extern /* Subroutine */ int zunglq_(), zungqr_(), zunbdb1_(), zunbdb2_(), 
	    zunbdb3_(), zunbdb4_();
    static logical wantv1t;


/*  -- LAPACK computational routine (version 3.6.0) -- */
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

#line 306 "zuncsd2by1.f"
    /* Parameter adjustments */
#line 306 "zuncsd2by1.f"
    x11_dim1 = *ldx11;
#line 306 "zuncsd2by1.f"
    x11_offset = 1 + x11_dim1;
#line 306 "zuncsd2by1.f"
    x11 -= x11_offset;
#line 306 "zuncsd2by1.f"
    x21_dim1 = *ldx21;
#line 306 "zuncsd2by1.f"
    x21_offset = 1 + x21_dim1;
#line 306 "zuncsd2by1.f"
    x21 -= x21_offset;
#line 306 "zuncsd2by1.f"
    --theta;
#line 306 "zuncsd2by1.f"
    u1_dim1 = *ldu1;
#line 306 "zuncsd2by1.f"
    u1_offset = 1 + u1_dim1;
#line 306 "zuncsd2by1.f"
    u1 -= u1_offset;
#line 306 "zuncsd2by1.f"
    u2_dim1 = *ldu2;
#line 306 "zuncsd2by1.f"
    u2_offset = 1 + u2_dim1;
#line 306 "zuncsd2by1.f"
    u2 -= u2_offset;
#line 306 "zuncsd2by1.f"
    v1t_dim1 = *ldv1t;
#line 306 "zuncsd2by1.f"
    v1t_offset = 1 + v1t_dim1;
#line 306 "zuncsd2by1.f"
    v1t -= v1t_offset;
#line 306 "zuncsd2by1.f"
    --work;
#line 306 "zuncsd2by1.f"
    --rwork;
#line 306 "zuncsd2by1.f"
    --iwork;
#line 306 "zuncsd2by1.f"

#line 306 "zuncsd2by1.f"
    /* Function Body */
#line 306 "zuncsd2by1.f"
    *info = 0;
#line 307 "zuncsd2by1.f"
    wantu1 = lsame_(jobu1, "Y", (ftnlen)1, (ftnlen)1);
#line 308 "zuncsd2by1.f"
    wantu2 = lsame_(jobu2, "Y", (ftnlen)1, (ftnlen)1);
#line 309 "zuncsd2by1.f"
    wantv1t = lsame_(jobv1t, "Y", (ftnlen)1, (ftnlen)1);
#line 310 "zuncsd2by1.f"
    lquery = *lwork == -1;

#line 312 "zuncsd2by1.f"
    if (*m < 0) {
#line 313 "zuncsd2by1.f"
	*info = -4;
#line 314 "zuncsd2by1.f"
    } else if (*p < 0 || *p > *m) {
#line 315 "zuncsd2by1.f"
	*info = -5;
#line 316 "zuncsd2by1.f"
    } else if (*q < 0 || *q > *m) {
#line 317 "zuncsd2by1.f"
	*info = -6;
#line 318 "zuncsd2by1.f"
    } else if (*ldx11 < max(1,*p)) {
#line 319 "zuncsd2by1.f"
	*info = -8;
#line 320 "zuncsd2by1.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 320 "zuncsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 320 "zuncsd2by1.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 321 "zuncsd2by1.f"
	    *info = -10;
#line 322 "zuncsd2by1.f"
	} else if (wantu1 && *ldu1 < *p) {
#line 323 "zuncsd2by1.f"
	    *info = -13;
#line 324 "zuncsd2by1.f"
	} else if (wantu2 && *ldu2 < *m - *p) {
#line 325 "zuncsd2by1.f"
	    *info = -15;
#line 326 "zuncsd2by1.f"
	} else if (wantv1t && *ldv1t < *q) {
#line 327 "zuncsd2by1.f"
	    *info = -17;
#line 328 "zuncsd2by1.f"
	}
#line 328 "zuncsd2by1.f"
    }

/* Computing MIN */
#line 330 "zuncsd2by1.f"
    i__1 = *p, i__2 = *m - *p, i__1 = min(i__1,i__2), i__1 = min(i__1,*q), 
	    i__2 = *m - *q;
#line 330 "zuncsd2by1.f"
    r__ = min(i__1,i__2);

/*     Compute workspace */

/*       WORK layout: */
/*     |-----------------------------------------| */
/*     | LWORKOPT (1)                            | */
/*     |-----------------------------------------| */
/*     | TAUP1 (MAX(1,P))                        | */
/*     | TAUP2 (MAX(1,M-P))                      | */
/*     | TAUQ1 (MAX(1,Q))                        | */
/*     |-----------------------------------------| */
/*     | ZUNBDB WORK | ZUNGQR WORK | ZUNGLQ WORK | */
/*     |             |             |             | */
/*     |             |             |             | */
/*     |             |             |             | */
/*     |             |             |             | */
/*     |-----------------------------------------| */
/*       RWORK layout: */
/*     |------------------| */
/*     | LRWORKOPT (1)    | */
/*     |------------------| */
/*     | PHI (MAX(1,R-1)) | */
/*     |------------------| */
/*     | B11D (R)         | */
/*     | B11E (R-1)       | */
/*     | B12D (R)         | */
/*     | B12E (R-1)       | */
/*     | B21D (R)         | */
/*     | B21E (R-1)       | */
/*     | B22D (R)         | */
/*     | B22E (R-1)       | */
/*     | ZBBCSD RWORK     | */
/*     |------------------| */

#line 365 "zuncsd2by1.f"
    if (*info == 0) {
#line 366 "zuncsd2by1.f"
	iphi = 2;
/* Computing MAX */
#line 367 "zuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 367 "zuncsd2by1.f"
	ib11d = iphi + max(i__1,i__2);
#line 368 "zuncsd2by1.f"
	ib11e = ib11d + max(1,r__);
/* Computing MAX */
#line 369 "zuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 369 "zuncsd2by1.f"
	ib12d = ib11e + max(i__1,i__2);
#line 370 "zuncsd2by1.f"
	ib12e = ib12d + max(1,r__);
/* Computing MAX */
#line 371 "zuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 371 "zuncsd2by1.f"
	ib21d = ib12e + max(i__1,i__2);
#line 372 "zuncsd2by1.f"
	ib21e = ib21d + max(1,r__);
/* Computing MAX */
#line 373 "zuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 373 "zuncsd2by1.f"
	ib22d = ib21e + max(i__1,i__2);
#line 374 "zuncsd2by1.f"
	ib22e = ib22d + max(1,r__);
/* Computing MAX */
#line 375 "zuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 375 "zuncsd2by1.f"
	ibbcsd = ib22e + max(i__1,i__2);
#line 376 "zuncsd2by1.f"
	itaup1 = 2;
#line 377 "zuncsd2by1.f"
	itaup2 = itaup1 + max(1,*p);
/* Computing MAX */
#line 378 "zuncsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 378 "zuncsd2by1.f"
	itauq1 = itaup2 + max(i__1,i__2);
#line 379 "zuncsd2by1.f"
	iorbdb = itauq1 + max(1,*q);
#line 380 "zuncsd2by1.f"
	iorgqr = itauq1 + max(1,*q);
#line 381 "zuncsd2by1.f"
	iorglq = itauq1 + max(1,*q);
#line 382 "zuncsd2by1.f"
	if (r__ == *q) {
#line 383 "zuncsd2by1.f"
	    zunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 385 "zuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 386 "zuncsd2by1.f"
	    if (*p >= *m - *p) {
#line 387 "zuncsd2by1.f"
		zungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 389 "zuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 390 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 391 "zuncsd2by1.f"
	    } else {
#line 392 "zuncsd2by1.f"
		i__1 = *m - *p;
#line 392 "zuncsd2by1.f"
		i__2 = *m - *p;
#line 392 "zuncsd2by1.f"
		zungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 394 "zuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 394 "zuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 395 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 396 "zuncsd2by1.f"
	    }
/* Computing MAX */
#line 397 "zuncsd2by1.f"
	    i__2 = 0, i__3 = *q - 1;
#line 397 "zuncsd2by1.f"
	    i__1 = max(i__2,i__3);
/* Computing MAX */
#line 397 "zuncsd2by1.f"
	    i__5 = 0, i__6 = *q - 1;
#line 397 "zuncsd2by1.f"
	    i__4 = max(i__5,i__6);
/* Computing MAX */
#line 397 "zuncsd2by1.f"
	    i__8 = 0, i__9 = *q - 1;
#line 397 "zuncsd2by1.f"
	    i__7 = max(i__8,i__9);
#line 397 "zuncsd2by1.f"
	    zunglq_(&i__1, &i__4, &i__7, &v1t[v1t_offset], ldv1t, &c__0, &
		    work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 399 "zuncsd2by1.f"
	    i__1 = 1, i__2 = *q - 1;
#line 399 "zuncsd2by1.f"
	    lorglqmin = max(i__1,i__2);
#line 400 "zuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 401 "zuncsd2by1.f"
	    zbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &c__0,
		     &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		    v1t_offset], ldv1t, &c__0, &c__1, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 404 "zuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 405 "zuncsd2by1.f"
	} else if (r__ == *p) {
#line 406 "zuncsd2by1.f"
	    zunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 408 "zuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 409 "zuncsd2by1.f"
	    if (*p - 1 >= *m - *p) {
#line 410 "zuncsd2by1.f"
		i__1 = *p - 1;
#line 410 "zuncsd2by1.f"
		i__2 = *p - 1;
#line 410 "zuncsd2by1.f"
		i__3 = *p - 1;
#line 410 "zuncsd2by1.f"
		zungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 412 "zuncsd2by1.f"
		i__1 = 1, i__2 = *p - 1;
#line 412 "zuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 413 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 414 "zuncsd2by1.f"
	    } else {
#line 415 "zuncsd2by1.f"
		i__1 = *m - *p;
#line 415 "zuncsd2by1.f"
		i__2 = *m - *p;
#line 415 "zuncsd2by1.f"
		zungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 417 "zuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 417 "zuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 418 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 419 "zuncsd2by1.f"
	    }
#line 420 "zuncsd2by1.f"
	    zunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 422 "zuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 423 "zuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 424 "zuncsd2by1.f"
	    zbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &c__0,
		     &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		    ldu1, &u2[u2_offset], ldu2, &c__0, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &childinfo, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 427 "zuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 428 "zuncsd2by1.f"
	} else if (r__ == *m - *p) {
#line 429 "zuncsd2by1.f"
	    zunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 431 "zuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 432 "zuncsd2by1.f"
	    if (*p >= *m - *p - 1) {
#line 433 "zuncsd2by1.f"
		zungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 435 "zuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 436 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 437 "zuncsd2by1.f"
	    } else {
#line 438 "zuncsd2by1.f"
		i__1 = *m - *p - 1;
#line 438 "zuncsd2by1.f"
		i__2 = *m - *p - 1;
#line 438 "zuncsd2by1.f"
		i__3 = *m - *p - 1;
#line 438 "zuncsd2by1.f"
		zungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 440 "zuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p - 1;
#line 440 "zuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 441 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 442 "zuncsd2by1.f"
	    }
#line 443 "zuncsd2by1.f"
	    zunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 445 "zuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 446 "zuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 447 "zuncsd2by1.f"
	    i__1 = *m - *q;
#line 447 "zuncsd2by1.f"
	    i__2 = *m - *p;
#line 447 "zuncsd2by1.f"
	    zbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1]
		    , &c__0, &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		    u2_offset], ldu2, &u1[u1_offset], ldu1, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1,
		     &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 451 "zuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 452 "zuncsd2by1.f"
	} else {
#line 453 "zuncsd2by1.f"
	    zunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &c__0, &
		    work[1], &c_n1, &childinfo);
#line 455 "zuncsd2by1.f"
	    lorbdb = *m + (integer) work[1].r;
#line 456 "zuncsd2by1.f"
	    if (*p >= *m - *p) {
#line 457 "zuncsd2by1.f"
		i__1 = *m - *q;
#line 457 "zuncsd2by1.f"
		zungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &c__0, &work[1], &
			c_n1, &childinfo);
#line 459 "zuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 460 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 461 "zuncsd2by1.f"
	    } else {
#line 462 "zuncsd2by1.f"
		i__1 = *m - *p;
#line 462 "zuncsd2by1.f"
		i__2 = *m - *p;
#line 462 "zuncsd2by1.f"
		i__3 = *m - *q;
#line 462 "zuncsd2by1.f"
		zungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &c__0, &
			work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 464 "zuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 464 "zuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 465 "zuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 466 "zuncsd2by1.f"
	    }
#line 467 "zuncsd2by1.f"
	    zunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &c_n1, 
		    &childinfo);
#line 469 "zuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 470 "zuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 471 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 471 "zuncsd2by1.f"
	    i__2 = *m - *q;
#line 471 "zuncsd2by1.f"
	    zbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1]
		    , &c__0, &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		    c__0, &c__1, &v1t[v1t_offset], ldv1t, &c__0, &c__0, &c__0,
		     &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 475 "zuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 476 "zuncsd2by1.f"
	}
#line 477 "zuncsd2by1.f"
	lrworkmin = ibbcsd + lbbcsd - 1;
#line 478 "zuncsd2by1.f"
	lrworkopt = lrworkmin;
#line 479 "zuncsd2by1.f"
	rwork[1] = (doublereal) lrworkopt;
/* Computing MAX */
#line 480 "zuncsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqrmin - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqmin - 1;
#line 480 "zuncsd2by1.f"
	lworkmin = max(i__1,i__2);
/* Computing MAX */
#line 483 "zuncsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqropt - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqopt - 1;
#line 483 "zuncsd2by1.f"
	lworkopt = max(i__1,i__2);
#line 486 "zuncsd2by1.f"
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
#line 487 "zuncsd2by1.f"
	if (*lwork < lworkmin && ! lquery) {
#line 488 "zuncsd2by1.f"
	    *info = -19;
#line 489 "zuncsd2by1.f"
	}
#line 490 "zuncsd2by1.f"
    }
#line 491 "zuncsd2by1.f"
    if (*info != 0) {
#line 492 "zuncsd2by1.f"
	i__1 = -(*info);
#line 492 "zuncsd2by1.f"
	xerbla_("ZUNCSD2BY1", &i__1, (ftnlen)10);
#line 493 "zuncsd2by1.f"
	return 0;
#line 494 "zuncsd2by1.f"
    } else if (lquery) {
#line 495 "zuncsd2by1.f"
	return 0;
#line 496 "zuncsd2by1.f"
    }
#line 497 "zuncsd2by1.f"
    lorgqr = *lwork - iorgqr + 1;
#line 498 "zuncsd2by1.f"
    lorglq = *lwork - iorglq + 1;

/*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q, */
/*     in which R = MIN(P,M-P,Q,M-Q) */

#line 503 "zuncsd2by1.f"
    if (r__ == *q) {

/*        Case 1: R = Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 509 "zuncsd2by1.f"
	zunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 515 "zuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 516 "zuncsd2by1.f"
	    zlacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 517 "zuncsd2by1.f"
	    zungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 519 "zuncsd2by1.f"
	}
#line 520 "zuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 521 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 521 "zuncsd2by1.f"
	    zlacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 522 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 522 "zuncsd2by1.f"
	    i__2 = *m - *p;
#line 522 "zuncsd2by1.f"
	    zungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 524 "zuncsd2by1.f"
	}
#line 525 "zuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 526 "zuncsd2by1.f"
	    i__1 = v1t_dim1 + 1;
#line 526 "zuncsd2by1.f"
	    v1t[i__1].r = 1., v1t[i__1].i = 0.;
#line 527 "zuncsd2by1.f"
	    i__1 = *q;
#line 527 "zuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 528 "zuncsd2by1.f"
		i__2 = j * v1t_dim1 + 1;
#line 528 "zuncsd2by1.f"
		v1t[i__2].r = 0., v1t[i__2].i = 0.;
#line 529 "zuncsd2by1.f"
		i__2 = j + v1t_dim1;
#line 529 "zuncsd2by1.f"
		v1t[i__2].r = 0., v1t[i__2].i = 0.;
#line 530 "zuncsd2by1.f"
	    }
#line 531 "zuncsd2by1.f"
	    i__1 = *q - 1;
#line 531 "zuncsd2by1.f"
	    i__2 = *q - 1;
#line 531 "zuncsd2by1.f"
	    zlacpy_("U", &i__1, &i__2, &x21[(x21_dim1 << 1) + 1], ldx21, &v1t[
		    (v1t_dim1 << 1) + 2], ldv1t, (ftnlen)1);
#line 533 "zuncsd2by1.f"
	    i__1 = *q - 1;
#line 533 "zuncsd2by1.f"
	    i__2 = *q - 1;
#line 533 "zuncsd2by1.f"
	    i__3 = *q - 1;
#line 533 "zuncsd2by1.f"
	    zunglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &
		    work[itauq1], &work[iorglq], &lorglq, &childinfo);
#line 535 "zuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 539 "zuncsd2by1.f"
	zbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &rwork[
		iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		v1t_offset], ldv1t, &c__0, &c__1, &rwork[ib11d], &rwork[ib11e]
		, &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], 
		&rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place zero submatrices in */
/*        preferred positions */

#line 549 "zuncsd2by1.f"
	if (*q > 0 && wantu2) {
#line 550 "zuncsd2by1.f"
	    i__1 = *q;
#line 550 "zuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 551 "zuncsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 552 "zuncsd2by1.f"
	    }
#line 553 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 553 "zuncsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 554 "zuncsd2by1.f"
		iwork[i__] = i__ - *q;
#line 555 "zuncsd2by1.f"
	    }
#line 556 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 556 "zuncsd2by1.f"
	    i__2 = *m - *p;
#line 556 "zuncsd2by1.f"
	    zlapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 557 "zuncsd2by1.f"
	}
#line 558 "zuncsd2by1.f"
    } else if (r__ == *p) {

/*        Case 2: R = P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 564 "zuncsd2by1.f"
	zunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 570 "zuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 571 "zuncsd2by1.f"
	    i__1 = u1_dim1 + 1;
#line 571 "zuncsd2by1.f"
	    u1[i__1].r = 1., u1[i__1].i = 0.;
#line 572 "zuncsd2by1.f"
	    i__1 = *p;
#line 572 "zuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 573 "zuncsd2by1.f"
		i__2 = j * u1_dim1 + 1;
#line 573 "zuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 574 "zuncsd2by1.f"
		i__2 = j + u1_dim1;
#line 574 "zuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 575 "zuncsd2by1.f"
	    }
#line 576 "zuncsd2by1.f"
	    i__1 = *p - 1;
#line 576 "zuncsd2by1.f"
	    i__2 = *p - 1;
#line 576 "zuncsd2by1.f"
	    zlacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 577 "zuncsd2by1.f"
	    i__1 = *p - 1;
#line 577 "zuncsd2by1.f"
	    i__2 = *p - 1;
#line 577 "zuncsd2by1.f"
	    i__3 = *p - 1;
#line 577 "zuncsd2by1.f"
	    zungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &work[
		    itaup1], &work[iorgqr], &lorgqr, &childinfo);
#line 579 "zuncsd2by1.f"
	}
#line 580 "zuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 581 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 581 "zuncsd2by1.f"
	    zlacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 582 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 582 "zuncsd2by1.f"
	    i__2 = *m - *p;
#line 582 "zuncsd2by1.f"
	    zungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 584 "zuncsd2by1.f"
	}
#line 585 "zuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 586 "zuncsd2by1.f"
	    zlacpy_("U", p, q, &x11[x11_offset], ldx11, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 587 "zuncsd2by1.f"
	    zunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 589 "zuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 593 "zuncsd2by1.f"
	zbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &rwork[
		iphi], &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		ldu1, &u2[u2_offset], ldu2, &rwork[ib11d], &rwork[ib11e], &
		rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], &
		rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 603 "zuncsd2by1.f"
	if (*q > 0 && wantu2) {
#line 604 "zuncsd2by1.f"
	    i__1 = *q;
#line 604 "zuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 605 "zuncsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 606 "zuncsd2by1.f"
	    }
#line 607 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 607 "zuncsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 608 "zuncsd2by1.f"
		iwork[i__] = i__ - *q;
#line 609 "zuncsd2by1.f"
	    }
#line 610 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 610 "zuncsd2by1.f"
	    i__2 = *m - *p;
#line 610 "zuncsd2by1.f"
	    zlapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 611 "zuncsd2by1.f"
	}
#line 612 "zuncsd2by1.f"
    } else if (r__ == *m - *p) {

/*        Case 3: R = M-P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 618 "zuncsd2by1.f"
	zunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 624 "zuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 625 "zuncsd2by1.f"
	    zlacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 626 "zuncsd2by1.f"
	    zungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 628 "zuncsd2by1.f"
	}
#line 629 "zuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 630 "zuncsd2by1.f"
	    i__1 = u2_dim1 + 1;
#line 630 "zuncsd2by1.f"
	    u2[i__1].r = 1., u2[i__1].i = 0.;
#line 631 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 631 "zuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 632 "zuncsd2by1.f"
		i__2 = j * u2_dim1 + 1;
#line 632 "zuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 633 "zuncsd2by1.f"
		i__2 = j + u2_dim1;
#line 633 "zuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 634 "zuncsd2by1.f"
	    }
#line 635 "zuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 635 "zuncsd2by1.f"
	    i__2 = *m - *p - 1;
#line 635 "zuncsd2by1.f"
	    zlacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 637 "zuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 637 "zuncsd2by1.f"
	    i__2 = *m - *p - 1;
#line 637 "zuncsd2by1.f"
	    i__3 = *m - *p - 1;
#line 637 "zuncsd2by1.f"
	    zungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &work[
		    itaup2], &work[iorgqr], &lorgqr, &childinfo);
#line 639 "zuncsd2by1.f"
	}
#line 640 "zuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 641 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 641 "zuncsd2by1.f"
	    zlacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 642 "zuncsd2by1.f"
	    zunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 644 "zuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 648 "zuncsd2by1.f"
	i__1 = *m - *q;
#line 648 "zuncsd2by1.f"
	i__2 = *m - *p;
#line 648 "zuncsd2by1.f"
	zbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1], &
		rwork[iphi], &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		u2_offset], ldu2, &u1[u1_offset], ldu1, &rwork[ib11d], &rwork[
		ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[
		ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd,
		 &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 658 "zuncsd2by1.f"
	if (*q > r__) {
#line 659 "zuncsd2by1.f"
	    i__1 = r__;
#line 659 "zuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 660 "zuncsd2by1.f"
		iwork[i__] = *q - r__ + i__;
#line 661 "zuncsd2by1.f"
	    }
#line 662 "zuncsd2by1.f"
	    i__1 = *q;
#line 662 "zuncsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 663 "zuncsd2by1.f"
		iwork[i__] = i__ - r__;
#line 664 "zuncsd2by1.f"
	    }
#line 665 "zuncsd2by1.f"
	    if (wantu1) {
#line 666 "zuncsd2by1.f"
		zlapmt_(&c_false, p, q, &u1[u1_offset], ldu1, &iwork[1]);
#line 667 "zuncsd2by1.f"
	    }
#line 668 "zuncsd2by1.f"
	    if (wantv1t) {
#line 669 "zuncsd2by1.f"
		zlapmr_(&c_false, q, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 670 "zuncsd2by1.f"
	    }
#line 671 "zuncsd2by1.f"
	}
#line 672 "zuncsd2by1.f"
    } else {

/*        Case 4: R = M-Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 678 "zuncsd2by1.f"
	i__1 = lorbdb - *m;
#line 678 "zuncsd2by1.f"
	zunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &work[iorbdb + *m], &i__1, &childinfo)
		;

/*        Accumulate Householder reflectors */

#line 685 "zuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 686 "zuncsd2by1.f"
	    zcopy_(p, &work[iorbdb], &c__1, &u1[u1_offset], &c__1);
#line 687 "zuncsd2by1.f"
	    i__1 = *p;
#line 687 "zuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 688 "zuncsd2by1.f"
		i__2 = j * u1_dim1 + 1;
#line 688 "zuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 689 "zuncsd2by1.f"
	    }
#line 690 "zuncsd2by1.f"
	    i__1 = *p - 1;
#line 690 "zuncsd2by1.f"
	    i__2 = *m - *q - 1;
#line 690 "zuncsd2by1.f"
	    zlacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 692 "zuncsd2by1.f"
	    i__1 = *m - *q;
#line 692 "zuncsd2by1.f"
	    zungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 694 "zuncsd2by1.f"
	}
#line 695 "zuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 696 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 696 "zuncsd2by1.f"
	    zcopy_(&i__1, &work[iorbdb + *p], &c__1, &u2[u2_offset], &c__1);
#line 697 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 697 "zuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 698 "zuncsd2by1.f"
		i__2 = j * u2_dim1 + 1;
#line 698 "zuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 699 "zuncsd2by1.f"
	    }
#line 700 "zuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 700 "zuncsd2by1.f"
	    i__2 = *m - *q - 1;
#line 700 "zuncsd2by1.f"
	    zlacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 702 "zuncsd2by1.f"
	    i__1 = *m - *p;
#line 702 "zuncsd2by1.f"
	    i__2 = *m - *p;
#line 702 "zuncsd2by1.f"
	    i__3 = *m - *q;
#line 702 "zuncsd2by1.f"
	    zungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &work[itaup2], 
		    &work[iorgqr], &lorgqr, &childinfo);
#line 704 "zuncsd2by1.f"
	}
#line 705 "zuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 706 "zuncsd2by1.f"
	    i__1 = *m - *q;
#line 706 "zuncsd2by1.f"
	    zlacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 707 "zuncsd2by1.f"
	    i__1 = *p - (*m - *q);
#line 707 "zuncsd2by1.f"
	    i__2 = *q - (*m - *q);
#line 707 "zuncsd2by1.f"
	    zlacpy_("U", &i__1, &i__2, &x11[*m - *q + 1 + (*m - *q + 1) * 
		    x11_dim1], ldx11, &v1t[*m - *q + 1 + (*m - *q + 1) * 
		    v1t_dim1], ldv1t, (ftnlen)1);
#line 709 "zuncsd2by1.f"
	    i__1 = -(*p) + *q;
#line 709 "zuncsd2by1.f"
	    i__2 = *q - *p;
#line 709 "zuncsd2by1.f"
	    zlacpy_("U", &i__1, &i__2, &x21[*m - *q + 1 + (*p + 1) * x21_dim1]
		    , ldx21, &v1t[*p + 1 + (*p + 1) * v1t_dim1], ldv1t, (
		    ftnlen)1);
#line 711 "zuncsd2by1.f"
	    zunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 713 "zuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 717 "zuncsd2by1.f"
	i__1 = *m - *p;
#line 717 "zuncsd2by1.f"
	i__2 = *m - *q;
#line 717 "zuncsd2by1.f"
	zbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1], &
		rwork[iphi], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		c__0, &c__1, &v1t[v1t_offset], ldv1t, &rwork[ib11d], &rwork[
		ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[
		ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd,
		 &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 727 "zuncsd2by1.f"
	if (*p > r__) {
#line 728 "zuncsd2by1.f"
	    i__1 = r__;
#line 728 "zuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 729 "zuncsd2by1.f"
		iwork[i__] = *p - r__ + i__;
#line 730 "zuncsd2by1.f"
	    }
#line 731 "zuncsd2by1.f"
	    i__1 = *p;
#line 731 "zuncsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 732 "zuncsd2by1.f"
		iwork[i__] = i__ - r__;
#line 733 "zuncsd2by1.f"
	    }
#line 734 "zuncsd2by1.f"
	    if (wantu1) {
#line 735 "zuncsd2by1.f"
		zlapmt_(&c_false, p, p, &u1[u1_offset], ldu1, &iwork[1]);
#line 736 "zuncsd2by1.f"
	    }
#line 737 "zuncsd2by1.f"
	    if (wantv1t) {
#line 738 "zuncsd2by1.f"
		zlapmr_(&c_false, p, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 739 "zuncsd2by1.f"
	    }
#line 740 "zuncsd2by1.f"
	}
#line 741 "zuncsd2by1.f"
    }

#line 743 "zuncsd2by1.f"
    return 0;

/*     End of ZUNCSD2BY1 */

} /* zuncsd2by1_ */

