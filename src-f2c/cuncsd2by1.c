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
/* > X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P, */
/* > (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are */
/* > R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in */
/* > which R = MIN(P,M-P,Q,M-Q). */
/* > */
/* >\endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU1 */
/* > \verbatim */
/* >          JOBU1 is CHARACTER */
/* >           = 'Y':      U1 is computed; */
/* >           otherwise:  U1 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU2 */
/* > \verbatim */
/* >          JOBU2 is CHARACTER */
/* >           = 'Y':      U2 is computed; */
/* >           otherwise:  U2 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV1T */
/* > \verbatim */
/* >          JOBV1T is CHARACTER */
/* >           = 'Y':      V1T is computed; */
/* >           otherwise:  V1T is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           The number of rows and columns in X. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >           The number of rows in X11 and X12. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is INTEGER */
/* >           The number of columns in X11 and X21. 0 <= Q <= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* >          X11 is COMPLEX array, dimension (LDX11,Q) */
/* >           On entry, part of the unitary matrix whose CSD is */
/* >           desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* >          LDX11 is INTEGER */
/* >           The leading dimension of X11. LDX11 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* >          X21 is COMPLEX array, dimension (LDX21,Q) */
/* >           On entry, part of the unitary matrix whose CSD is */
/* >           desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* >          LDX21 is INTEGER */
/* >           The leading dimension of X21. LDX21 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* >          THETA is COMPLEX array, dimension (R), in which R = */
/* >           MIN(P,M-P,Q,M-Q). */
/* >           C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* >           S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* >          U1 is COMPLEX array, dimension (P) */
/* >           If JOBU1 = 'Y', U1 contains the P-by-P unitary matrix U1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU1 */
/* > \verbatim */
/* >          LDU1 is INTEGER */
/* >           The leading dimension of U1. If JOBU1 = 'Y', LDU1 >= */
/* >           MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* >          U2 is COMPLEX array, dimension (M-P) */
/* >           If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) unitary */
/* >           matrix U2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* >          LDU2 is INTEGER */
/* >           The leading dimension of U2. If JOBU2 = 'Y', LDU2 >= */
/* >           MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] V1T */
/* > \verbatim */
/* >          V1T is COMPLEX array, dimension (Q) */
/* >           If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary */
/* >           matrix V1**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV1T */
/* > \verbatim */
/* >          LDV1T is INTEGER */
/* >           The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >= */
/* >           MAX(1,Q). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >           On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* >           If INFO > 0 on exit, WORK(2:R) contains the values PHI(1), */
/* >           ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* >           define the matrix in intermediate bidiagonal-block form */
/* >           remaining after nonconvergence. INFO specifies the number */
/* >           of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >           The dimension of the array WORK. */
/* > \endverbatim */
/* > \verbatim */
/* >           If LWORK = -1, then a workspace query is assumed; the routine */
/* >           only calculates the optimal size of the WORK array, returns */
/* >           this value as the first entry of the work array, and no error */
/* >           message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* >           On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* >           If INFO > 0 on exit, RWORK(2:R) contains the values PHI(1), */
/* >           ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* >           define the matrix in intermediate bidiagonal-block form */
/* >           remaining after nonconvergence. INFO specifies the number */
/* >           of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >           The dimension of the array RWORK. */
/* > */
/* >           If LRWORK = -1, then a workspace query is assumed; the routine */
/* >           only calculates the optimal size of the RWORK array, returns */
/* >           this value as the first entry of the work array, and no error */
/* >           message related to LRWORK is issued by XERBLA. */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q)) */
/* > \endverbatim */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:  successful exit. */
/* >           < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >           > 0:  CBBCSD did not converge. See the description of WORK */
/* >                above for details. */
/* > \endverbatim */

/* >  \par References: */
/* >  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */

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

#line 315 "cuncsd2by1.f"
    /* Parameter adjustments */
#line 315 "cuncsd2by1.f"
    x11_dim1 = *ldx11;
#line 315 "cuncsd2by1.f"
    x11_offset = 1 + x11_dim1;
#line 315 "cuncsd2by1.f"
    x11 -= x11_offset;
#line 315 "cuncsd2by1.f"
    x21_dim1 = *ldx21;
#line 315 "cuncsd2by1.f"
    x21_offset = 1 + x21_dim1;
#line 315 "cuncsd2by1.f"
    x21 -= x21_offset;
#line 315 "cuncsd2by1.f"
    --theta;
#line 315 "cuncsd2by1.f"
    u1_dim1 = *ldu1;
#line 315 "cuncsd2by1.f"
    u1_offset = 1 + u1_dim1;
#line 315 "cuncsd2by1.f"
    u1 -= u1_offset;
#line 315 "cuncsd2by1.f"
    u2_dim1 = *ldu2;
#line 315 "cuncsd2by1.f"
    u2_offset = 1 + u2_dim1;
#line 315 "cuncsd2by1.f"
    u2 -= u2_offset;
#line 315 "cuncsd2by1.f"
    v1t_dim1 = *ldv1t;
#line 315 "cuncsd2by1.f"
    v1t_offset = 1 + v1t_dim1;
#line 315 "cuncsd2by1.f"
    v1t -= v1t_offset;
#line 315 "cuncsd2by1.f"
    --work;
#line 315 "cuncsd2by1.f"
    --rwork;
#line 315 "cuncsd2by1.f"
    --iwork;
#line 315 "cuncsd2by1.f"

#line 315 "cuncsd2by1.f"
    /* Function Body */
#line 315 "cuncsd2by1.f"
    *info = 0;
#line 316 "cuncsd2by1.f"
    wantu1 = lsame_(jobu1, "Y", (ftnlen)1, (ftnlen)1);
#line 317 "cuncsd2by1.f"
    wantu2 = lsame_(jobu2, "Y", (ftnlen)1, (ftnlen)1);
#line 318 "cuncsd2by1.f"
    wantv1t = lsame_(jobv1t, "Y", (ftnlen)1, (ftnlen)1);
#line 319 "cuncsd2by1.f"
    lquery = *lwork == -1;

#line 321 "cuncsd2by1.f"
    if (*m < 0) {
#line 322 "cuncsd2by1.f"
	*info = -4;
#line 323 "cuncsd2by1.f"
    } else if (*p < 0 || *p > *m) {
#line 324 "cuncsd2by1.f"
	*info = -5;
#line 325 "cuncsd2by1.f"
    } else if (*q < 0 || *q > *m) {
#line 326 "cuncsd2by1.f"
	*info = -6;
#line 327 "cuncsd2by1.f"
    } else if (*ldx11 < max(1,*p)) {
#line 328 "cuncsd2by1.f"
	*info = -8;
#line 329 "cuncsd2by1.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 329 "cuncsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 329 "cuncsd2by1.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 330 "cuncsd2by1.f"
	    *info = -10;
#line 331 "cuncsd2by1.f"
	} else if (wantu1 && *ldu1 < *p) {
#line 332 "cuncsd2by1.f"
	    *info = -13;
#line 333 "cuncsd2by1.f"
	} else if (wantu2 && *ldu2 < *m - *p) {
#line 334 "cuncsd2by1.f"
	    *info = -15;
#line 335 "cuncsd2by1.f"
	} else if (wantv1t && *ldv1t < *q) {
#line 336 "cuncsd2by1.f"
	    *info = -17;
#line 337 "cuncsd2by1.f"
	}
#line 337 "cuncsd2by1.f"
    }

/* Computing MIN */
#line 339 "cuncsd2by1.f"
    i__1 = *p, i__2 = *m - *p, i__1 = min(i__1,i__2), i__1 = min(i__1,*q), 
	    i__2 = *m - *q;
#line 339 "cuncsd2by1.f"
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

#line 374 "cuncsd2by1.f"
    if (*info == 0) {
#line 375 "cuncsd2by1.f"
	iphi = 2;
/* Computing MAX */
#line 376 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 376 "cuncsd2by1.f"
	ib11d = iphi + max(i__1,i__2);
#line 377 "cuncsd2by1.f"
	ib11e = ib11d + max(1,r__);
/* Computing MAX */
#line 378 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 378 "cuncsd2by1.f"
	ib12d = ib11e + max(i__1,i__2);
#line 379 "cuncsd2by1.f"
	ib12e = ib12d + max(1,r__);
/* Computing MAX */
#line 380 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 380 "cuncsd2by1.f"
	ib21d = ib12e + max(i__1,i__2);
#line 381 "cuncsd2by1.f"
	ib21e = ib21d + max(1,r__);
/* Computing MAX */
#line 382 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 382 "cuncsd2by1.f"
	ib22d = ib21e + max(i__1,i__2);
#line 383 "cuncsd2by1.f"
	ib22e = ib22d + max(1,r__);
/* Computing MAX */
#line 384 "cuncsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 384 "cuncsd2by1.f"
	ibbcsd = ib22e + max(i__1,i__2);
#line 385 "cuncsd2by1.f"
	itaup1 = 2;
#line 386 "cuncsd2by1.f"
	itaup2 = itaup1 + max(1,*p);
/* Computing MAX */
#line 387 "cuncsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 387 "cuncsd2by1.f"
	itauq1 = itaup2 + max(i__1,i__2);
#line 388 "cuncsd2by1.f"
	iorbdb = itauq1 + max(1,*q);
#line 389 "cuncsd2by1.f"
	iorgqr = itauq1 + max(1,*q);
#line 390 "cuncsd2by1.f"
	iorglq = itauq1 + max(1,*q);
#line 391 "cuncsd2by1.f"
	if (r__ == *q) {
#line 392 "cuncsd2by1.f"
	    cunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 394 "cuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 395 "cuncsd2by1.f"
	    if (*p >= *m - *p) {
#line 396 "cuncsd2by1.f"
		cungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 398 "cuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 399 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 400 "cuncsd2by1.f"
	    } else {
#line 401 "cuncsd2by1.f"
		i__1 = *m - *p;
#line 401 "cuncsd2by1.f"
		i__2 = *m - *p;
#line 401 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 403 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 403 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 404 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 405 "cuncsd2by1.f"
	    }
/* Computing MAX */
#line 406 "cuncsd2by1.f"
	    i__2 = 0, i__3 = *q - 1;
#line 406 "cuncsd2by1.f"
	    i__1 = max(i__2,i__3);
/* Computing MAX */
#line 406 "cuncsd2by1.f"
	    i__5 = 0, i__6 = *q - 1;
#line 406 "cuncsd2by1.f"
	    i__4 = max(i__5,i__6);
/* Computing MAX */
#line 406 "cuncsd2by1.f"
	    i__8 = 0, i__9 = *q - 1;
#line 406 "cuncsd2by1.f"
	    i__7 = max(i__8,i__9);
#line 406 "cuncsd2by1.f"
	    cunglq_(&i__1, &i__4, &i__7, &v1t[v1t_offset], ldv1t, &c__0, &
		    work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 408 "cuncsd2by1.f"
	    i__1 = 1, i__2 = *q - 1;
#line 408 "cuncsd2by1.f"
	    lorglqmin = max(i__1,i__2);
#line 409 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 410 "cuncsd2by1.f"
	    cbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &c__0,
		     &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		    v1t_offset], ldv1t, &c__0, &c__1, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 413 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 414 "cuncsd2by1.f"
	} else if (r__ == *p) {
#line 415 "cuncsd2by1.f"
	    cunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 417 "cuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 418 "cuncsd2by1.f"
	    if (*p - 1 >= *m - *p) {
#line 419 "cuncsd2by1.f"
		i__1 = *p - 1;
#line 419 "cuncsd2by1.f"
		i__2 = *p - 1;
#line 419 "cuncsd2by1.f"
		i__3 = *p - 1;
#line 419 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 421 "cuncsd2by1.f"
		i__1 = 1, i__2 = *p - 1;
#line 421 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 422 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 423 "cuncsd2by1.f"
	    } else {
#line 424 "cuncsd2by1.f"
		i__1 = *m - *p;
#line 424 "cuncsd2by1.f"
		i__2 = *m - *p;
#line 424 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1]
			, &c_n1, &childinfo);
/* Computing MAX */
#line 426 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 426 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 427 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 428 "cuncsd2by1.f"
	    }
#line 429 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 431 "cuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 432 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 433 "cuncsd2by1.f"
	    cbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &c__0,
		     &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		    ldu1, &u2[u2_offset], ldu2, &c__0, &c__0, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &childinfo, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 436 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 437 "cuncsd2by1.f"
	} else if (r__ == *m - *p) {
#line 438 "cuncsd2by1.f"
	    cunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 440 "cuncsd2by1.f"
	    lorbdb = (integer) work[1].r;
#line 441 "cuncsd2by1.f"
	    if (*p >= *m - *p - 1) {
#line 442 "cuncsd2by1.f"
		cungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1,
			 &childinfo);
#line 444 "cuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 445 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 446 "cuncsd2by1.f"
	    } else {
#line 447 "cuncsd2by1.f"
		i__1 = *m - *p - 1;
#line 447 "cuncsd2by1.f"
		i__2 = *m - *p - 1;
#line 447 "cuncsd2by1.f"
		i__3 = *m - *p - 1;
#line 447 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &
			c__0, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 449 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p - 1;
#line 449 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 450 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 451 "cuncsd2by1.f"
	    }
#line 452 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &
		    c_n1, &childinfo);
#line 454 "cuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 455 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 456 "cuncsd2by1.f"
	    i__1 = *m - *q;
#line 456 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 456 "cuncsd2by1.f"
	    cbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1]
		    , &c__0, &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		    u2_offset], ldu2, &u1[u1_offset], ldu1, &c__0, &c__0, &
		    c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1,
		     &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 460 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 461 "cuncsd2by1.f"
	} else {
#line 462 "cuncsd2by1.f"
	    cunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &c__0, &
		    work[1], &c_n1, &childinfo);
#line 464 "cuncsd2by1.f"
	    lorbdb = *m + (integer) work[1].r;
#line 465 "cuncsd2by1.f"
	    if (*p >= *m - *p) {
#line 466 "cuncsd2by1.f"
		i__1 = *m - *q;
#line 466 "cuncsd2by1.f"
		cungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &c__0, &work[1], &
			c_n1, &childinfo);
#line 468 "cuncsd2by1.f"
		lorgqrmin = max(1,*p);
#line 469 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 470 "cuncsd2by1.f"
	    } else {
#line 471 "cuncsd2by1.f"
		i__1 = *m - *p;
#line 471 "cuncsd2by1.f"
		i__2 = *m - *p;
#line 471 "cuncsd2by1.f"
		i__3 = *m - *q;
#line 471 "cuncsd2by1.f"
		cungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &c__0, &
			work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 473 "cuncsd2by1.f"
		i__1 = 1, i__2 = *m - *p;
#line 473 "cuncsd2by1.f"
		lorgqrmin = max(i__1,i__2);
#line 474 "cuncsd2by1.f"
		lorgqropt = (integer) work[1].r;
#line 475 "cuncsd2by1.f"
	    }
#line 476 "cuncsd2by1.f"
	    cunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &c_n1, 
		    &childinfo);
#line 478 "cuncsd2by1.f"
	    lorglqmin = max(1,*q);
#line 479 "cuncsd2by1.f"
	    lorglqopt = (integer) work[1].r;
#line 480 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 480 "cuncsd2by1.f"
	    i__2 = *m - *q;
#line 480 "cuncsd2by1.f"
	    cbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1]
		    , &c__0, &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		    c__0, &c__1, &v1t[v1t_offset], ldv1t, &c__0, &c__0, &c__0,
		     &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &
		    childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 484 "cuncsd2by1.f"
	    lbbcsd = (integer) rwork[1];
#line 485 "cuncsd2by1.f"
	}
#line 486 "cuncsd2by1.f"
	lrworkmin = ibbcsd + lbbcsd - 1;
#line 487 "cuncsd2by1.f"
	lrworkopt = lrworkmin;
#line 488 "cuncsd2by1.f"
	rwork[1] = (doublereal) lrworkopt;
/* Computing MAX */
#line 489 "cuncsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqrmin - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqmin - 1;
#line 489 "cuncsd2by1.f"
	lworkmin = max(i__1,i__2);
/* Computing MAX */
#line 492 "cuncsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqropt - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqopt - 1;
#line 492 "cuncsd2by1.f"
	lworkopt = max(i__1,i__2);
#line 495 "cuncsd2by1.f"
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
#line 496 "cuncsd2by1.f"
	if (*lwork < lworkmin && ! lquery) {
#line 497 "cuncsd2by1.f"
	    *info = -19;
#line 498 "cuncsd2by1.f"
	}
#line 499 "cuncsd2by1.f"
    }
#line 500 "cuncsd2by1.f"
    if (*info != 0) {
#line 501 "cuncsd2by1.f"
	i__1 = -(*info);
#line 501 "cuncsd2by1.f"
	xerbla_("CUNCSD2BY1", &i__1, (ftnlen)10);
#line 502 "cuncsd2by1.f"
	return 0;
#line 503 "cuncsd2by1.f"
    } else if (lquery) {
#line 504 "cuncsd2by1.f"
	return 0;
#line 505 "cuncsd2by1.f"
    }
#line 506 "cuncsd2by1.f"
    lorgqr = *lwork - iorgqr + 1;
#line 507 "cuncsd2by1.f"
    lorglq = *lwork - iorglq + 1;

/*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q, */
/*     in which R = MIN(P,M-P,Q,M-Q) */

#line 512 "cuncsd2by1.f"
    if (r__ == *q) {

/*        Case 1: R = Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 518 "cuncsd2by1.f"
	cunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 524 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 525 "cuncsd2by1.f"
	    clacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 526 "cuncsd2by1.f"
	    cungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 528 "cuncsd2by1.f"
	}
#line 529 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 530 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 530 "cuncsd2by1.f"
	    clacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 531 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 531 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 531 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 533 "cuncsd2by1.f"
	}
#line 534 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 535 "cuncsd2by1.f"
	    i__1 = v1t_dim1 + 1;
#line 535 "cuncsd2by1.f"
	    v1t[i__1].r = 1., v1t[i__1].i = 0.;
#line 536 "cuncsd2by1.f"
	    i__1 = *q;
#line 536 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 537 "cuncsd2by1.f"
		i__2 = j * v1t_dim1 + 1;
#line 537 "cuncsd2by1.f"
		v1t[i__2].r = 0., v1t[i__2].i = 0.;
#line 538 "cuncsd2by1.f"
		i__2 = j + v1t_dim1;
#line 538 "cuncsd2by1.f"
		v1t[i__2].r = 0., v1t[i__2].i = 0.;
#line 539 "cuncsd2by1.f"
	    }
#line 540 "cuncsd2by1.f"
	    i__1 = *q - 1;
#line 540 "cuncsd2by1.f"
	    i__2 = *q - 1;
#line 540 "cuncsd2by1.f"
	    clacpy_("U", &i__1, &i__2, &x21[(x21_dim1 << 1) + 1], ldx21, &v1t[
		    (v1t_dim1 << 1) + 2], ldv1t, (ftnlen)1);
#line 542 "cuncsd2by1.f"
	    i__1 = *q - 1;
#line 542 "cuncsd2by1.f"
	    i__2 = *q - 1;
#line 542 "cuncsd2by1.f"
	    i__3 = *q - 1;
#line 542 "cuncsd2by1.f"
	    cunglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &
		    work[itauq1], &work[iorglq], &lorglq, &childinfo);
#line 544 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 548 "cuncsd2by1.f"
	cbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &rwork[
		iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		v1t_offset], ldv1t, &c__0, &c__1, &rwork[ib11d], &rwork[ib11e]
		, &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], 
		&rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place zero submatrices in */
/*        preferred positions */

#line 558 "cuncsd2by1.f"
	if (*q > 0 && wantu2) {
#line 559 "cuncsd2by1.f"
	    i__1 = *q;
#line 559 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 560 "cuncsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 561 "cuncsd2by1.f"
	    }
#line 562 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 562 "cuncsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 563 "cuncsd2by1.f"
		iwork[i__] = i__ - *q;
#line 564 "cuncsd2by1.f"
	    }
#line 565 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 565 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 565 "cuncsd2by1.f"
	    clapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 566 "cuncsd2by1.f"
	}
#line 567 "cuncsd2by1.f"
    } else if (r__ == *p) {

/*        Case 2: R = P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 573 "cuncsd2by1.f"
	cunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 579 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 580 "cuncsd2by1.f"
	    i__1 = u1_dim1 + 1;
#line 580 "cuncsd2by1.f"
	    u1[i__1].r = 1., u1[i__1].i = 0.;
#line 581 "cuncsd2by1.f"
	    i__1 = *p;
#line 581 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 582 "cuncsd2by1.f"
		i__2 = j * u1_dim1 + 1;
#line 582 "cuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 583 "cuncsd2by1.f"
		i__2 = j + u1_dim1;
#line 583 "cuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 584 "cuncsd2by1.f"
	    }
#line 585 "cuncsd2by1.f"
	    i__1 = *p - 1;
#line 585 "cuncsd2by1.f"
	    i__2 = *p - 1;
#line 585 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 586 "cuncsd2by1.f"
	    i__1 = *p - 1;
#line 586 "cuncsd2by1.f"
	    i__2 = *p - 1;
#line 586 "cuncsd2by1.f"
	    i__3 = *p - 1;
#line 586 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &work[
		    itaup1], &work[iorgqr], &lorgqr, &childinfo);
#line 588 "cuncsd2by1.f"
	}
#line 589 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 590 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 590 "cuncsd2by1.f"
	    clacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 591 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 591 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 591 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 593 "cuncsd2by1.f"
	}
#line 594 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 595 "cuncsd2by1.f"
	    clacpy_("U", p, q, &x11[x11_offset], ldx11, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 596 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 598 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 602 "cuncsd2by1.f"
	cbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &rwork[
		iphi], &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], 
		ldu1, &u2[u2_offset], ldu2, &rwork[ib11d], &rwork[ib11e], &
		rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], &
		rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 612 "cuncsd2by1.f"
	if (*q > 0 && wantu2) {
#line 613 "cuncsd2by1.f"
	    i__1 = *q;
#line 613 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 614 "cuncsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 615 "cuncsd2by1.f"
	    }
#line 616 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 616 "cuncsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 617 "cuncsd2by1.f"
		iwork[i__] = i__ - *q;
#line 618 "cuncsd2by1.f"
	    }
#line 619 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 619 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 619 "cuncsd2by1.f"
	    clapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 620 "cuncsd2by1.f"
	}
#line 621 "cuncsd2by1.f"
    } else if (r__ == *m - *p) {

/*        Case 3: R = M-P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 627 "cuncsd2by1.f"
	cunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 633 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 634 "cuncsd2by1.f"
	    clacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 635 "cuncsd2by1.f"
	    cungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 637 "cuncsd2by1.f"
	}
#line 638 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 639 "cuncsd2by1.f"
	    i__1 = u2_dim1 + 1;
#line 639 "cuncsd2by1.f"
	    u2[i__1].r = 1., u2[i__1].i = 0.;
#line 640 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 640 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 641 "cuncsd2by1.f"
		i__2 = j * u2_dim1 + 1;
#line 641 "cuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 642 "cuncsd2by1.f"
		i__2 = j + u2_dim1;
#line 642 "cuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 643 "cuncsd2by1.f"
	    }
#line 644 "cuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 644 "cuncsd2by1.f"
	    i__2 = *m - *p - 1;
#line 644 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 646 "cuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 646 "cuncsd2by1.f"
	    i__2 = *m - *p - 1;
#line 646 "cuncsd2by1.f"
	    i__3 = *m - *p - 1;
#line 646 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &work[
		    itaup2], &work[iorgqr], &lorgqr, &childinfo);
#line 648 "cuncsd2by1.f"
	}
#line 649 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 650 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 650 "cuncsd2by1.f"
	    clacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 651 "cuncsd2by1.f"
	    cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 653 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 657 "cuncsd2by1.f"
	i__1 = *m - *q;
#line 657 "cuncsd2by1.f"
	i__2 = *m - *p;
#line 657 "cuncsd2by1.f"
	cbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1], &
		rwork[iphi], &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		u2_offset], ldu2, &u1[u1_offset], ldu1, &rwork[ib11d], &rwork[
		ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[
		ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd,
		 &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 667 "cuncsd2by1.f"
	if (*q > r__) {
#line 668 "cuncsd2by1.f"
	    i__1 = r__;
#line 668 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 669 "cuncsd2by1.f"
		iwork[i__] = *q - r__ + i__;
#line 670 "cuncsd2by1.f"
	    }
#line 671 "cuncsd2by1.f"
	    i__1 = *q;
#line 671 "cuncsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 672 "cuncsd2by1.f"
		iwork[i__] = i__ - r__;
#line 673 "cuncsd2by1.f"
	    }
#line 674 "cuncsd2by1.f"
	    if (wantu1) {
#line 675 "cuncsd2by1.f"
		clapmt_(&c_false, p, q, &u1[u1_offset], ldu1, &iwork[1]);
#line 676 "cuncsd2by1.f"
	    }
#line 677 "cuncsd2by1.f"
	    if (wantv1t) {
#line 678 "cuncsd2by1.f"
		clapmr_(&c_false, q, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 679 "cuncsd2by1.f"
	    }
#line 680 "cuncsd2by1.f"
	}
#line 681 "cuncsd2by1.f"
    } else {

/*        Case 4: R = M-Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 687 "cuncsd2by1.f"
	i__1 = lorbdb - *m;
#line 687 "cuncsd2by1.f"
	cunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &work[iorbdb + *m], &i__1, &childinfo)
		;

/*        Accumulate Householder reflectors */

#line 694 "cuncsd2by1.f"
	if (wantu1 && *p > 0) {
#line 695 "cuncsd2by1.f"
	    ccopy_(p, &work[iorbdb], &c__1, &u1[u1_offset], &c__1);
#line 696 "cuncsd2by1.f"
	    i__1 = *p;
#line 696 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 697 "cuncsd2by1.f"
		i__2 = j * u1_dim1 + 1;
#line 697 "cuncsd2by1.f"
		u1[i__2].r = 0., u1[i__2].i = 0.;
#line 698 "cuncsd2by1.f"
	    }
#line 699 "cuncsd2by1.f"
	    i__1 = *p - 1;
#line 699 "cuncsd2by1.f"
	    i__2 = *m - *q - 1;
#line 699 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 701 "cuncsd2by1.f"
	    i__1 = *m - *q;
#line 701 "cuncsd2by1.f"
	    cungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 703 "cuncsd2by1.f"
	}
#line 704 "cuncsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 705 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 705 "cuncsd2by1.f"
	    ccopy_(&i__1, &work[iorbdb + *p], &c__1, &u2[u2_offset], &c__1);
#line 706 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 706 "cuncsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 707 "cuncsd2by1.f"
		i__2 = j * u2_dim1 + 1;
#line 707 "cuncsd2by1.f"
		u2[i__2].r = 0., u2[i__2].i = 0.;
#line 708 "cuncsd2by1.f"
	    }
#line 709 "cuncsd2by1.f"
	    i__1 = *m - *p - 1;
#line 709 "cuncsd2by1.f"
	    i__2 = *m - *q - 1;
#line 709 "cuncsd2by1.f"
	    clacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 711 "cuncsd2by1.f"
	    i__1 = *m - *p;
#line 711 "cuncsd2by1.f"
	    i__2 = *m - *p;
#line 711 "cuncsd2by1.f"
	    i__3 = *m - *q;
#line 711 "cuncsd2by1.f"
	    cungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &work[itaup2], 
		    &work[iorgqr], &lorgqr, &childinfo);
#line 713 "cuncsd2by1.f"
	}
#line 714 "cuncsd2by1.f"
	if (wantv1t && *q > 0) {
#line 715 "cuncsd2by1.f"
	    i__1 = *m - *q;
#line 715 "cuncsd2by1.f"
	    clacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 716 "cuncsd2by1.f"
	    i__1 = *p - (*m - *q);
#line 716 "cuncsd2by1.f"
	    i__2 = *q - (*m - *q);
#line 716 "cuncsd2by1.f"
	    clacpy_("U", &i__1, &i__2, &x11[*m - *q + 1 + (*m - *q + 1) * 
		    x11_dim1], ldx11, &v1t[*m - *q + 1 + (*m - *q + 1) * 
		    v1t_dim1], ldv1t, (ftnlen)1);
#line 718 "cuncsd2by1.f"
	    i__1 = -(*p) + *q;
#line 718 "cuncsd2by1.f"
	    i__2 = *q - *p;
#line 718 "cuncsd2by1.f"
	    clacpy_("U", &i__1, &i__2, &x21[*m - *q + 1 + (*p + 1) * x21_dim1]
		    , ldx21, &v1t[*p + 1 + (*p + 1) * v1t_dim1], ldv1t, (
		    ftnlen)1);
#line 720 "cuncsd2by1.f"
	    cunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 722 "cuncsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 726 "cuncsd2by1.f"
	i__1 = *m - *p;
#line 726 "cuncsd2by1.f"
	i__2 = *m - *q;
#line 726 "cuncsd2by1.f"
	cbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1], &
		rwork[iphi], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, &
		c__0, &c__1, &v1t[v1t_offset], ldv1t, &rwork[ib11d], &rwork[
		ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[
		ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd,
		 &childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 736 "cuncsd2by1.f"
	if (*p > r__) {
#line 737 "cuncsd2by1.f"
	    i__1 = r__;
#line 737 "cuncsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 738 "cuncsd2by1.f"
		iwork[i__] = *p - r__ + i__;
#line 739 "cuncsd2by1.f"
	    }
#line 740 "cuncsd2by1.f"
	    i__1 = *p;
#line 740 "cuncsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 741 "cuncsd2by1.f"
		iwork[i__] = i__ - r__;
#line 742 "cuncsd2by1.f"
	    }
#line 743 "cuncsd2by1.f"
	    if (wantu1) {
#line 744 "cuncsd2by1.f"
		clapmt_(&c_false, p, p, &u1[u1_offset], ldu1, &iwork[1]);
#line 745 "cuncsd2by1.f"
	    }
#line 746 "cuncsd2by1.f"
	    if (wantv1t) {
#line 747 "cuncsd2by1.f"
		clapmr_(&c_false, p, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 748 "cuncsd2by1.f"
	    }
#line 749 "cuncsd2by1.f"
	}
#line 750 "cuncsd2by1.f"
    }

#line 752 "cuncsd2by1.f"
    return 0;

/*     End of CUNCSD2BY1 */

} /* cuncsd2by1_ */

