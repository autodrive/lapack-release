#line 1 "cuncsd2by1.f"
/* cuncsd2by1.f -- translated by f2c (version 20100827).
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

#line 1 "cuncsd2by1.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__1 = 1;
static logical c_false = FALSE_;

/* > \brief \b CUNCSD2BY1 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNCSD2BY1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cuncsd2
by1.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cuncsd2
by1.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cuncsd2
by1.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, */
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
/*       REAL               RWORK(*) */
/*       REAL               THETA(*) */
/*       COMPLEX            U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), */
/*      $                   X11(LDX11,*), X21(LDX21,*) */
/*       INTEGER            IWORK(*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > CUNCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with */
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
/* > */
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
/* >          X11 is COMPLEX array, dimension (LDX11,Q) */
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
/* >          X21 is COMPLEX array, dimension (LDX21,Q) */
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
/* >          THETA is REAL array, dimension (R), in which R = */
/* >          MIN(P,M-P,Q,M-Q). */
/* >          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* >          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* >          U1 is COMPLEX array, dimension (P) */
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
/* >          U2 is COMPLEX array, dimension (M-P) */
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
/* >          V1T is COMPLEX array, dimension (Q) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
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
/* >          > 0:  CBBCSD did not converge. See the description of WORK */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cuncsd2by1_(char *jobu1, char *jobu2, char *jobv1t, 
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
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer itaup1, itaup2, itauq1;
    static logical wantu1, wantu2;
    extern /* Subroutine */ int cbbcsd_();
    static integer ibbcsd, lbbcsd, iorbdb, lorbdb;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), clapmr_(logical *, integer *, 
	    integer *, doublecomplex *, integer *, integer *), clapmt_(
	    logical *, integer *, integer *, doublecomplex *, integer *, 
	    integer *), cunglq_();
    static integer iorglq;
    extern /* Subroutine */ int cungqr_();
    static integer lorglq, iorgqr, lorgqr;
    extern /* Subroutine */ int cunbdb1_(), cunbdb2_();
    static logical lquery;
    extern /* Subroutine */ int cunbdb3_(), cunbdb4_();
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

#line 307 "cuncsd2by1.f"
    /* Parameter adjustments */
#line 307 "cuncsd2by1.f"
    x11_dim1 = *ldx11;
#line 307 "cuncsd2by1.f"
    x11_offset = 1 + x11_dim1;
#line 307 "cuncsd2by1.f"
    x11 -= x11_offset;
#line 307 "cuncsd2by1.f"
    x21_dim1 = *ldx21;
#line 307 "cuncsd2by1.f"
    x21_offset = 1 + x21_dim1;
#line 307 "cuncsd2by1.f"
    x21 -= x21_offset;
#line 307 "cuncsd2by1.f"
    --theta;
#line 307 "cuncsd2by1.f"
    u1_dim1 = *ldu1;
#line 307 "cuncsd2by1.f"
    u1_offset = 1 + u1_dim1;
#line 307 "cuncsd2by1.f"
    u1 -= u1_offset;
#line 307 "cuncsd2by1.f"
    u2_dim1 = *ldu2;
#line 307 "cuncsd2by1.f"
    u2_offset = 1 + u2_dim1;
#line 307 "cuncsd2by1.f"
    u2 -= u2_offset;
#line 307 "cuncsd2by1.f"
    v1t_dim1 = *ldv1t;
#line 307 "cuncsd2by1.f"
    v1t_offset = 1 + v1t_dim1;
#line 307 "cuncsd2by1.f"
    v1t -= v1t_offset;
#line 307 "cuncsd2by1.f"
    --work;
#line 307 "cuncsd2by1.f"
    --rwork;
#line 307 "cuncsd2by1.f"
    --iwork;
#line 307 "cuncsd2by1.f"

#line 307 "cuncsd2by1.f"
    /* Function Body */
#line 307 "cuncsd2by1.f"
    *info = 0;
#line 308 "cuncsd2by1.f"
    wantu1 = lsame_(jobu1, "Y", (ftnlen)1, (ftnlen)1);
#line 309 "cuncsd2by1.f"
    wantu2 = lsame_(jobu2, "Y", (ftnlen)1, (ftnlen)1);
#line 310 "cuncsd2by1.f"
    wantv1t = lsame_(jobv1t, "Y", (ftnlen)1, (ftnlen)1);
#line 311 "cuncsd2by1.f"
    lquery = *lwork == -1;

#line 313 "cuncsd2by1.f"
    if (*m < 0) {
#line 314 "cuncsd2by1.f"
	*info = -4;
#line 315 "cuncsd2by1.f"
    } else if (*p < 0 || *p > *m) {
#line 316 "cuncsd2by1.f"
	*info = -5;
#line 317 "cuncsd2by1.f"
    } else if (*q < 0 || *q > *m) {
#line 318 "cuncsd2by1.f"
	*info = -6;
#line 319 "cuncsd2by1.f"
    } else if (*ldx11 < max(1,*p)) {
#line 320 "cuncsd2by1.f"
	*info = -8;
#line 321 "cuncsd2by1.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 321 "cuncsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 321 "cuncsd2by1.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 322 "cuncsd2by1.f"
	    *info = -10;
#line 323 "cuncsd2by1.f"
	} else if (wantu1 && *ldu1 < *p) {
#line 324 "cuncsd2by1.f"
	    *info = -13;
#line 325 "cuncsd2by1.f"
	} else if (wantu2 && *ldu2 < *m - *p) {
#line 326 "cuncsd2by1.f"
	    *info = -15;
#line 327 "cuncsd2by1.f"
	} else if (wantv1t && *ldv1t < *q) {
#line 328 "cuncsd2by1.f"
	    *info = -17;
#line 329 "cuncsd2by1.f"
	}
#line 329 "cuncsd2by1.f"
    }

/* Computing MIN */
#line 331 "cuncsd2by1.f"
    i__1 = *p, i__2 = *m - *p, i__1 = min(i__1,i__2), i__1 = min(i__1,*q), 
	    i__2 = *m - *q;
#line 331 "cuncsd2by1.f"
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
/*     | CUNBDB WORK | CUNGQR WORK | CUNGLQ WORK | */
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
/*     | CBBCSD RWORK     | */
/*     |------------------| */

#line 366 "cuncsd2by1.f"
    if (*info == 0) {
#line 367 "cuncsd2by1.f"
	iphi = 2;
/* Computing MAX */
#line 368 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 368 "cuncsd2by1.f"
	ib11d = iphi + max(i__1,i__2);
#line 369 "cuncsd2by1.f"
	ib11e = ib11d + max(1,r__);
/* Computing MAX */
#line 370 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 370 "cuncsd2by1.f"
	ib12d = ib11e + max(i__1,i__2);
#line 371 "cuncsd2by1.f"
	ib12e = ib12d + max(1,r__);
/* Computing MAX */
#line 372 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 372 "cuncsd2by1.f"
	ib21d = ib12e + max(i__1,i__2);
#line 373 "cuncsd2by1.f"
	ib21e = ib21d + max(1,r__);
/* Computing MAX */
#line 374 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 374 "cuncsd2by1.f"
	ib22d = ib21e + max(i__1,i__2);
#line 375 "cuncsd2by1.f"
	ib22e = ib22d + max(1,r__);
/* Computing MAX */
#line 376 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 376 "cuncsd2by1.f"
	ibbcsd = ib22e + max(i__1,i__2);
#line 377 "cuncsd2by1.f"
	itaup1 = 2;
#line 378 "cuncsd2by1.f"
	itaup2 = itaup1 + max(1,*p);
/* Computing MAX */
#line 379 "cuncsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 379 "cuncsd2by1.f"
	itauq1 = itaup2 + max(i__1,i__2);
#line 380 "cuncsd2by1.f"
	iorbdb = itauq1 + max(1,*q);
#line 381 "cuncsd2by1.f"
	iorgqr = itauq1 + max(1,*q);
#line 382 "cuncsd2by1.f"
	iorglq = itauq1 + max(1,*q);
#line 383 "cuncsd2by1.f"
	if (r__ == *q) {
#line 384 "cuncsd2by1.f"
	    cunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 386 "cuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 387 "cuncsd2by1.f"
	    if (*p >= *m - *p) {
#line 388 "cuncsd2by1.f"
		cungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 390 "cuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 391 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 392 "cuncsd2by1.f"
	    } else {
#line 393 "cuncsd2by1.f"
		i__1 = *m - *p;
#line 393 "cuncsd2by1.f"
		i__2 = *m - *p;
#line 393 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 395 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 395 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 396 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 397 "cuncsd2by1.f"
	    }
/* Computing MAX */
#line 398 "cuncsd2by1.f"
	    i__2 = 0, i__3 = *q - 1;
#line 398 "cuncsd2by1.f"
	    i__1 = max(i__2,i__3);
/* Computing MAX */
#line 398 "cuncsd2by1.f"
	    i__5 = 0, i__6 = *q - 1;
#line 398 "cuncsd2by1.f"
	    i__4 = max(i__5,i__6);
/* Computing MAX */
#line 398 "cuncsd2by1.f"
	    i__8 = 0, i__9 = *q - 1;
#line 398 "cuncsd2by1.f"
	    i__7 = max(i__8,i__9);
#line 398 "cuncsd2by1.f"
	    cunglq_(&i__1, &i__4, &i__7, &v1t[v1t_offset], ldv1t, &c__0, &
		    work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 400 "cuncsd2by1.f"
	    i__1 = 1, i__2 = *q - 1;
#line 400 "cuncsd2by1.f"
	    lorglqmin = max(i__1,i__2);
#line 401 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 402 "cuncsd2by1.f"
	    cbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &c__0,
		     &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		    v1t_offset], ldv1t, &c__0, &c__1, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 405 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 406 "cuncsd2by1.f"
	} else if (r__ == *p) {
#line 407 "cuncsd2by1.f"
	    cunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 409 "cuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 410 "cuncsd2by1.f"
	    if (*p - 1 >= *m - *p) {
#line 411 "cuncsd2by1.f"
		i__1 = *p - 1;
#line 411 "cuncsd2by1.f"
		i__2 = *p - 1;
#line 411 "cuncsd2by1.f"
		i__3 = *p - 1;
#line 411 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 413 "cuncsd2by1.f"
		i__1 = 1, i__2 = *p - 1;
#line 413 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 414 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 415 "cuncsd2by1.f"
	    } else {
#line 416 "cuncsd2by1.f"
		i__1 = *m - *p;
#line 416 "cuncsd2by1.f"
		i__2 = *m - *p;
#line 416 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 418 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 418 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 419 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 420 "cuncsd2by1.f"
	    }
#line 421 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 423 "cuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 424 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 425 "cuncsd2by1.f"
	    cbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &c__0,
		     &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		    ldu1, &u2[u2_offset], ldu2, &c__0, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &childinfo, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 428 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 429 "cuncsd2by1.f"
	} else if (r__ == *m - *p) {
#line 430 "cuncsd2by1.f"
	    cunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 432 "cuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 433 "cuncsd2by1.f"
	    if (*p >= *m - *p - 1) {
#line 434 "cuncsd2by1.f"
		cungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 436 "cuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 437 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 438 "cuncsd2by1.f"
	    } else {
#line 439 "cuncsd2by1.f"
		i__1 = *m - *p - 1;
#line 439 "cuncsd2by1.f"
		i__2 = *m - *p - 1;
#line 439 "cuncsd2by1.f"
		i__3 = *m - *p - 1;
#line 439 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 441 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p - 1;
#line 441 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 442 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 443 "cuncsd2by1.f"
	    }
#line 444 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 446 "cuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 447 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 448 "cuncsd2by1.f"
	    i__1 = *m - *q;
#line 448 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 448 "cuncsd2by1.f"
	    cbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1]
		    , &c__0, &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		    u2_offset], ldu2, &u1[u1_offset], ldu1, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1,
		     &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 452 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 453 "cuncsd2by1.f"
	} else {
#line 454 "cuncsd2by1.f"
	    cunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &c__0, &
		    work[1], &c_n1, &childinfo);
#line 456 "cuncsd2by1.f"
	    lorbdb = *m + (integer) work[1].r;
#line 457 "cuncsd2by1.f"
	    if (*p >= *m - *p) {
#line 458 "cuncsd2by1.f"
		i__1 = *m - *q;
#line 458 "cuncsd2by1.f"
		cungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &c__0, &work[1], &
			c_n1, &childinfo);
#line 460 "cuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 461 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 462 "cuncsd2by1.f"
	    } else {
#line 463 "cuncsd2by1.f"
		i__1 = *m - *p;
#line 463 "cuncsd2by1.f"
		i__2 = *m - *p;
#line 463 "cuncsd2by1.f"
		i__3 = *m - *q;
#line 463 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &c__0, &
			work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 465 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 465 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 466 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 467 "cuncsd2by1.f"
	    }
#line 468 "cuncsd2by1.f"
	    cunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &c_n1, 
		    &childinfo);
#line 470 "cuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 471 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 472 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 472 "cuncsd2by1.f"
	    i__2 = *m - *q;
#line 472 "cuncsd2by1.f"
	    cbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1]
		    , &c__0, &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		    c__0, &c__1, &v1t[v1t_offset], ldv1t, &c__0, &c__0, &c__0,
		     &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 476 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 477 "cuncsd2by1.f"
	}
#line 478 "cuncsd2by1.f"
	lrworkmin = ibbcsd + lbbcsd - 1;
#line 479 "cuncsd2by1.f"
	lrworkopt = lrworkmin;
#line 480 "cuncsd2by1.f"
	rwork[1] = (doublereal) lrworkopt;
/* Computing MAX */
#line 481 "cuncsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqrmin - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqmin - 1;
#line 481 "cuncsd2by1.f"
	lworkmin = max(i__1,i__2);
/* Computing MAX */
#line 484 "cuncsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqropt - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqopt - 1;
#line 484 "cuncsd2by1.f"
	lworkopt = max(i__1,i__2);
#line 487 "cuncsd2by1.f"
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
#line 488 "cuncsd2by1.f"
	if (*lwork < lworkmin && ! lquery) {
#line 489 "cuncsd2by1.f"
	    *info = -19;
#line 490 "cuncsd2by1.f"
	}
#line 491 "cuncsd2by1.f"
    }
#line 492 "cuncsd2by1.f"
    if (*info != 0) {
#line 493 "cuncsd2by1.f"
	i__1 = -(*info);
#line 493 "cuncsd2by1.f"
	xerbla_("CUNCSD2BY1", &i__1, (ftnlen)10);
#line 494 "cuncsd2by1.f"
	return 0;
#line 495 "cuncsd2by1.f"
    } else if (lquery) {
#line 496 "cuncsd2by1.f"
	return 0;
#line 497 "cuncsd2by1.f"
    }
#line 498 "cuncsd2by1.f"
    lorgqr = *lwork - iorgqr + 1;
#line 499 "cuncsd2by1.f"
    lorglq = *lwork - iorglq + 1;

/*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q, */
/*     in which R = MIN(P,M-P,Q,M-Q) */

#line 504 "cuncsd2by1.f"
    if (r__ == *q) {

/*        Case 1: R = Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 510 "cuncsd2by1.f"
	cunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 516 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 517 "cuncsd2by1.f"
	    clacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 518 "cuncsd2by1.f"
	    cungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 520 "cuncsd2by1.f"
	}
#line 521 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 522 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 522 "cuncsd2by1.f"
	    clacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 523 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 523 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 523 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 525 "cuncsd2by1.f"
	}
#line 526 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 527 "cuncsd2by1.f"
	    i__1 = v1t_dim1 + 1;
#line 527 "cuncsd2by1.f"
	    v1t[i__1].r = 1., v1t[i__1].i = 0.;
#line 528 "cuncsd2by1.f"
	    i__1 = *q;
#line 528 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 529 "cuncsd2by1.f"
		i__2 = j * v1t_dim1 + 1;
#line 529 "cuncsd2by1.f"
		v1t[i__2].r = 0., v1t[i__2].i = 0.;
#line 530 "cuncsd2by1.f"
		i__2 = j + v1t_dim1;
#line 530 "cuncsd2by1.f"
		v1t[i__2].r = 0., v1t[i__2].i = 0.;
#line 531 "cuncsd2by1.f"
	    }
#line 532 "cuncsd2by1.f"
	    i__1 = *q - 1;
#line 532 "cuncsd2by1.f"
	    i__2 = *q - 1;
#line 532 "cuncsd2by1.f"
	    clacpy_("U", &i__1, &i__2, &x21[(x21_dim1 << 1) + 1], ldx21, &v1t[
		    (v1t_dim1 << 1) + 2], ldv1t, (ftnlen)1);
#line 534 "cuncsd2by1.f"
	    i__1 = *q - 1;
#line 534 "cuncsd2by1.f"
	    i__2 = *q - 1;
#line 534 "cuncsd2by1.f"
	    i__3 = *q - 1;
#line 534 "cuncsd2by1.f"
	    cunglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &
		    work[itauq1], &work[iorglq], &lorglq, &childinfo);
#line 536 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 540 "cuncsd2by1.f"
	cbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &rwork[
		iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		v1t_offset], ldv1t, &c__0, &c__1, &rwork[ib11d], &rwork[ib11e]
		, &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], 
		&rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place zero submatrices in */
/*        preferred positions */

#line 550 "cuncsd2by1.f"
	if (*q > 0 && wantu2) {
#line 551 "cuncsd2by1.f"
	    i__1 = *q;
#line 551 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 552 "cuncsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 553 "cuncsd2by1.f"
	    }
#line 554 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 554 "cuncsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 555 "cuncsd2by1.f"
		iwork[i__] = i__ - *q;
#line 556 "cuncsd2by1.f"
	    }
#line 557 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 557 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 557 "cuncsd2by1.f"
	    clapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 558 "cuncsd2by1.f"
	}
#line 559 "cuncsd2by1.f"
    } else if (r__ == *p) {

/*        Case 2: R = P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 565 "cuncsd2by1.f"
	cunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 571 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 572 "cuncsd2by1.f"
	    i__1 = u1_dim1 + 1;
#line 572 "cuncsd2by1.f"
	    u1[i__1].r = 1., u1[i__1].i = 0.;
#line 573 "cuncsd2by1.f"
	    i__1 = *p;
#line 573 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 574 "cuncsd2by1.f"
		i__2 = j * u1_dim1 + 1;
#line 574 "cuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 575 "cuncsd2by1.f"
		i__2 = j + u1_dim1;
#line 575 "cuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 576 "cuncsd2by1.f"
	    }
#line 577 "cuncsd2by1.f"
	    i__1 = *p - 1;
#line 577 "cuncsd2by1.f"
	    i__2 = *p - 1;
#line 577 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 578 "cuncsd2by1.f"
	    i__1 = *p - 1;
#line 578 "cuncsd2by1.f"
	    i__2 = *p - 1;
#line 578 "cuncsd2by1.f"
	    i__3 = *p - 1;
#line 578 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &work[
		    itaup1], &work[iorgqr], &lorgqr, &childinfo);
#line 580 "cuncsd2by1.f"
	}
#line 581 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 582 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 582 "cuncsd2by1.f"
	    clacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 583 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 583 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 583 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 585 "cuncsd2by1.f"
	}
#line 586 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 587 "cuncsd2by1.f"
	    clacpy_("U", p, q, &x11[x11_offset], ldx11, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 588 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 590 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 594 "cuncsd2by1.f"
	cbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &rwork[
		iphi], &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		ldu1, &u2[u2_offset], ldu2, &rwork[ib11d], &rwork[ib11e], &
		rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], &
		rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 604 "cuncsd2by1.f"
	if (*q > 0 && wantu2) {
#line 605 "cuncsd2by1.f"
	    i__1 = *q;
#line 605 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 606 "cuncsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 607 "cuncsd2by1.f"
	    }
#line 608 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 608 "cuncsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 609 "cuncsd2by1.f"
		iwork[i__] = i__ - *q;
#line 610 "cuncsd2by1.f"
	    }
#line 611 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 611 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 611 "cuncsd2by1.f"
	    clapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 612 "cuncsd2by1.f"
	}
#line 613 "cuncsd2by1.f"
    } else if (r__ == *m - *p) {

/*        Case 3: R = M-P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 619 "cuncsd2by1.f"
	cunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 625 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 626 "cuncsd2by1.f"
	    clacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 627 "cuncsd2by1.f"
	    cungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 629 "cuncsd2by1.f"
	}
#line 630 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 631 "cuncsd2by1.f"
	    i__1 = u2_dim1 + 1;
#line 631 "cuncsd2by1.f"
	    u2[i__1].r = 1., u2[i__1].i = 0.;
#line 632 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 632 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 633 "cuncsd2by1.f"
		i__2 = j * u2_dim1 + 1;
#line 633 "cuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 634 "cuncsd2by1.f"
		i__2 = j + u2_dim1;
#line 634 "cuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 635 "cuncsd2by1.f"
	    }
#line 636 "cuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 636 "cuncsd2by1.f"
	    i__2 = *m - *p - 1;
#line 636 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 638 "cuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 638 "cuncsd2by1.f"
	    i__2 = *m - *p - 1;
#line 638 "cuncsd2by1.f"
	    i__3 = *m - *p - 1;
#line 638 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &work[
		    itaup2], &work[iorgqr], &lorgqr, &childinfo);
#line 640 "cuncsd2by1.f"
	}
#line 641 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 642 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 642 "cuncsd2by1.f"
	    clacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 643 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 645 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 649 "cuncsd2by1.f"
	i__1 = *m - *q;
#line 649 "cuncsd2by1.f"
	i__2 = *m - *p;
#line 649 "cuncsd2by1.f"
	cbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1], &
		rwork[iphi], &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		u2_offset], ldu2, &u1[u1_offset], ldu1, &rwork[ib11d], &rwork[
		ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[
		ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd,
		 &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 659 "cuncsd2by1.f"
	if (*q > r__) {
#line 660 "cuncsd2by1.f"
	    i__1 = r__;
#line 660 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 661 "cuncsd2by1.f"
		iwork[i__] = *q - r__ + i__;
#line 662 "cuncsd2by1.f"
	    }
#line 663 "cuncsd2by1.f"
	    i__1 = *q;
#line 663 "cuncsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 664 "cuncsd2by1.f"
		iwork[i__] = i__ - r__;
#line 665 "cuncsd2by1.f"
	    }
#line 666 "cuncsd2by1.f"
	    if (wantu1) {
#line 667 "cuncsd2by1.f"
		clapmt_(&c_false, p, q, &u1[u1_offset], ldu1, &iwork[1]);
#line 668 "cuncsd2by1.f"
	    }
#line 669 "cuncsd2by1.f"
	    if (wantv1t) {
#line 670 "cuncsd2by1.f"
		clapmr_(&c_false, q, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 671 "cuncsd2by1.f"
	    }
#line 672 "cuncsd2by1.f"
	}
#line 673 "cuncsd2by1.f"
    } else {

/*        Case 4: R = M-Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 679 "cuncsd2by1.f"
	i__1 = lorbdb - *m;
#line 679 "cuncsd2by1.f"
	cunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &work[iorbdb + *m], &i__1, &childinfo)
		;

/*        Accumulate Householder reflectors */

#line 686 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 687 "cuncsd2by1.f"
	    ccopy_(p, &work[iorbdb], &c__1, &u1[u1_offset], &c__1);
#line 688 "cuncsd2by1.f"
	    i__1 = *p;
#line 688 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 689 "cuncsd2by1.f"
		i__2 = j * u1_dim1 + 1;
#line 689 "cuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 690 "cuncsd2by1.f"
	    }
#line 691 "cuncsd2by1.f"
	    i__1 = *p - 1;
#line 691 "cuncsd2by1.f"
	    i__2 = *m - *q - 1;
#line 691 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 693 "cuncsd2by1.f"
	    i__1 = *m - *q;
#line 693 "cuncsd2by1.f"
	    cungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 695 "cuncsd2by1.f"
	}
#line 696 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 697 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 697 "cuncsd2by1.f"
	    ccopy_(&i__1, &work[iorbdb + *p], &c__1, &u2[u2_offset], &c__1);
#line 698 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 698 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 699 "cuncsd2by1.f"
		i__2 = j * u2_dim1 + 1;
#line 699 "cuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 700 "cuncsd2by1.f"
	    }
#line 701 "cuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 701 "cuncsd2by1.f"
	    i__2 = *m - *q - 1;
#line 701 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 703 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 703 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 703 "cuncsd2by1.f"
	    i__3 = *m - *q;
#line 703 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &work[itaup2], 
		    &work[iorgqr], &lorgqr, &childinfo);
#line 705 "cuncsd2by1.f"
	}
#line 706 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 707 "cuncsd2by1.f"
	    i__1 = *m - *q;
#line 707 "cuncsd2by1.f"
	    clacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 708 "cuncsd2by1.f"
	    i__1 = *p - (*m - *q);
#line 708 "cuncsd2by1.f"
	    i__2 = *q - (*m - *q);
#line 708 "cuncsd2by1.f"
	    clacpy_("U", &i__1, &i__2, &x11[*m - *q + 1 + (*m - *q + 1) * 
		    x11_dim1], ldx11, &v1t[*m - *q + 1 + (*m - *q + 1) * 
		    v1t_dim1], ldv1t, (ftnlen)1);
#line 710 "cuncsd2by1.f"
	    i__1 = -(*p) + *q;
#line 710 "cuncsd2by1.f"
	    i__2 = *q - *p;
#line 710 "cuncsd2by1.f"
	    clacpy_("U", &i__1, &i__2, &x21[*m - *q + 1 + (*p + 1) * x21_dim1]
		    , ldx21, &v1t[*p + 1 + (*p + 1) * v1t_dim1], ldv1t, (
		    ftnlen)1);
#line 712 "cuncsd2by1.f"
	    cunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 714 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 718 "cuncsd2by1.f"
	i__1 = *m - *p;
#line 718 "cuncsd2by1.f"
	i__2 = *m - *q;
#line 718 "cuncsd2by1.f"
	cbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1], &
		rwork[iphi], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		c__0, &c__1, &v1t[v1t_offset], ldv1t, &rwork[ib11d], &rwork[
		ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[
		ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd,
		 &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 728 "cuncsd2by1.f"
	if (*p > r__) {
#line 729 "cuncsd2by1.f"
	    i__1 = r__;
#line 729 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 730 "cuncsd2by1.f"
		iwork[i__] = *p - r__ + i__;
#line 731 "cuncsd2by1.f"
	    }
#line 732 "cuncsd2by1.f"
	    i__1 = *p;
#line 732 "cuncsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 733 "cuncsd2by1.f"
		iwork[i__] = i__ - r__;
#line 734 "cuncsd2by1.f"
	    }
#line 735 "cuncsd2by1.f"
	    if (wantu1) {
#line 736 "cuncsd2by1.f"
		clapmt_(&c_false, p, p, &u1[u1_offset], ldu1, &iwork[1]);
#line 737 "cuncsd2by1.f"
	    }
#line 738 "cuncsd2by1.f"
	    if (wantv1t) {
#line 739 "cuncsd2by1.f"
		clapmr_(&c_false, p, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 740 "cuncsd2by1.f"
	    }
#line 741 "cuncsd2by1.f"
	}
#line 742 "cuncsd2by1.f"
    }

#line 744 "cuncsd2by1.f"
    return 0;

/*     End of CUNCSD2BY1 */

} /* cuncsd2by1_ */

