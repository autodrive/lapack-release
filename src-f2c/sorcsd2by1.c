#line 1 "sorcsd2by1.f"
/* sorcsd2by1.f -- translated by f2c (version 20100827).
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

#line 1 "sorcsd2by1.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;
static logical c_false = FALSE_;

/* > \brief \b SORCSD2BY1 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORCSD2BY1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorcsd2
by1.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorcsd2
by1.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorcsd2
by1.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, */
/*                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T, */
/*                              LDV1T, WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU1, JOBU2, JOBV1T */
/*       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, */
/*      $                   M, P, Q */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               THETA(*) */
/*       REAL               U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), */
/*      $                   X11(LDX11,*), X21(LDX21,*) */
/*       INTEGER            IWORK(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* >\verbatim */
/* > */
/* > SORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with */
/* > orthonormal columns that has been partitioned into a 2-by-1 block */
/* > structure: */
/* > */
/* >                                [  I1 0  0 ] */
/* >                                [  0  C  0 ] */
/* >          [ X11 ]   [ U1 |    ] [  0  0  0 ] */
/* >      X = [-----] = [---------] [----------] V1**T . */
/* >          [ X21 ]   [    | U2 ] [  0  0  0 ] */
/* >                                [  0  S  0 ] */
/* >                                [  0  0  I2] */
/* > */
/* > X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P, */
/* > (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R */
/* > nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which */
/* > R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a */
/* > K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0). */
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
/* >          X11 is REAL array, dimension (LDX11,Q) */
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
/* >          X21 is REAL array, dimension (LDX21,Q) */
/* >          On entry, part of the orthogonal matrix whose CSD is desired. */
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
/* >          THETA is REAL array, dimension (R), in which R = */
/* >          MIN(P,M-P,Q,M-Q). */
/* >          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* >          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* >          U1 is REAL array, dimension (P) */
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
/* >          U2 is REAL array, dimension (M-P) */
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
/* >          V1T is REAL array, dimension (Q) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* >          > 0:  SBBCSD did not converge. See the description of WORK */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorcsd2by1_(char *jobu1, char *jobu2, char *jobv1t, 
	integer *m, integer *p, integer *q, doublereal *x11, integer *ldx11, 
	doublereal *x21, integer *ldx21, doublereal *theta, doublereal *u1, 
	integer *ldu1, doublereal *u2, integer *ldu2, doublereal *v1t, 
	integer *ldv1t, doublereal *work, integer *lwork, integer *iwork, 
	integer *info, ftnlen jobu1_len, ftnlen jobu2_len, ftnlen jobv1t_len)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, 
	    x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer lworkmin, lworkopt, i__, j, r__, childinfo, lorglqmin, 
	    lorgqrmin, lorglqopt, lorgqropt;
    static doublereal dum1[1], dum2[1]	/* was [1][1] */;
    static integer ib11d, ib11e, ib12d, ib12e, ib21d, ib21e, ib22d, ib22e, 
	    iphi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer itaup1, itaup2, itauq1;
    static logical wantu1, wantu2;
    static integer ibbcsd, lbbcsd;
    extern /* Subroutine */ int sbbcsd_(char *, char *, char *, char *, char *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer iorbdb, lorbdb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
    static integer iorglq;
    extern /* Subroutine */ int slapmr_(logical *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    static integer lorglq;
    extern /* Subroutine */ int slapmt_(logical *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    static integer iorgqr, lorgqr;
    extern /* Subroutine */ int sorglq_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), sorgqr_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *);
    static logical lquery;
    extern /* Subroutine */ int sorbdb1_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), sorbdb2_(integer *, integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), sorbdb3_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    sorbdb4_(integer *, integer *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, integer *);
    static logical wantv1t;


/*  -- LAPACK computational routine (version 3.7.1) -- */
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Function .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input arguments */

#line 287 "sorcsd2by1.f"
    /* Parameter adjustments */
#line 287 "sorcsd2by1.f"
    x11_dim1 = *ldx11;
#line 287 "sorcsd2by1.f"
    x11_offset = 1 + x11_dim1;
#line 287 "sorcsd2by1.f"
    x11 -= x11_offset;
#line 287 "sorcsd2by1.f"
    x21_dim1 = *ldx21;
#line 287 "sorcsd2by1.f"
    x21_offset = 1 + x21_dim1;
#line 287 "sorcsd2by1.f"
    x21 -= x21_offset;
#line 287 "sorcsd2by1.f"
    --theta;
#line 287 "sorcsd2by1.f"
    u1_dim1 = *ldu1;
#line 287 "sorcsd2by1.f"
    u1_offset = 1 + u1_dim1;
#line 287 "sorcsd2by1.f"
    u1 -= u1_offset;
#line 287 "sorcsd2by1.f"
    u2_dim1 = *ldu2;
#line 287 "sorcsd2by1.f"
    u2_offset = 1 + u2_dim1;
#line 287 "sorcsd2by1.f"
    u2 -= u2_offset;
#line 287 "sorcsd2by1.f"
    v1t_dim1 = *ldv1t;
#line 287 "sorcsd2by1.f"
    v1t_offset = 1 + v1t_dim1;
#line 287 "sorcsd2by1.f"
    v1t -= v1t_offset;
#line 287 "sorcsd2by1.f"
    --work;
#line 287 "sorcsd2by1.f"
    --iwork;
#line 287 "sorcsd2by1.f"

#line 287 "sorcsd2by1.f"
    /* Function Body */
#line 287 "sorcsd2by1.f"
    *info = 0;
#line 288 "sorcsd2by1.f"
    wantu1 = lsame_(jobu1, "Y", (ftnlen)1, (ftnlen)1);
#line 289 "sorcsd2by1.f"
    wantu2 = lsame_(jobu2, "Y", (ftnlen)1, (ftnlen)1);
#line 290 "sorcsd2by1.f"
    wantv1t = lsame_(jobv1t, "Y", (ftnlen)1, (ftnlen)1);
#line 291 "sorcsd2by1.f"
    lquery = *lwork == -1;

#line 293 "sorcsd2by1.f"
    if (*m < 0) {
#line 294 "sorcsd2by1.f"
	*info = -4;
#line 295 "sorcsd2by1.f"
    } else if (*p < 0 || *p > *m) {
#line 296 "sorcsd2by1.f"
	*info = -5;
#line 297 "sorcsd2by1.f"
    } else if (*q < 0 || *q > *m) {
#line 298 "sorcsd2by1.f"
	*info = -6;
#line 299 "sorcsd2by1.f"
    } else if (*ldx11 < max(1,*p)) {
#line 300 "sorcsd2by1.f"
	*info = -8;
#line 301 "sorcsd2by1.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 301 "sorcsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 301 "sorcsd2by1.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 302 "sorcsd2by1.f"
	    *info = -10;
#line 303 "sorcsd2by1.f"
	} else if (wantu1 && *ldu1 < max(1,*p)) {
#line 304 "sorcsd2by1.f"
	    *info = -13;
#line 305 "sorcsd2by1.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 305 "sorcsd2by1.f"
	    i__1 = 1, i__2 = *m - *p;
#line 305 "sorcsd2by1.f"
	    if (wantu2 && *ldu2 < max(i__1,i__2)) {
#line 306 "sorcsd2by1.f"
		*info = -15;
#line 307 "sorcsd2by1.f"
	    } else if (wantv1t && *ldv1t < max(1,*q)) {
#line 308 "sorcsd2by1.f"
		*info = -17;
#line 309 "sorcsd2by1.f"
	    }
#line 309 "sorcsd2by1.f"
	}
#line 309 "sorcsd2by1.f"
    }

/* Computing MIN */
#line 311 "sorcsd2by1.f"
    i__1 = *p, i__2 = *m - *p, i__1 = min(i__1,i__2), i__1 = min(i__1,*q), 
	    i__2 = *m - *q;
#line 311 "sorcsd2by1.f"
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
/*     | SORBDB WORK | SORGQR WORK | SORGLQ WORK | B21D (R)    | */
/*     |             |             |             | B21E (R-1)  | */
/*     |             |             |             | B22D (R)    | */
/*     |             |             |             | B22E (R-1)  | */
/*     |             |             |             | SBBCSD WORK | */
/*     |-------------------------------------------------------| */

#line 332 "sorcsd2by1.f"
    if (*info == 0) {
#line 333 "sorcsd2by1.f"
	iphi = 2;
/* Computing MAX */
#line 334 "sorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 334 "sorcsd2by1.f"
	ib11d = iphi + max(i__1,i__2);
#line 335 "sorcsd2by1.f"
	ib11e = ib11d + max(1,r__);
/* Computing MAX */
#line 336 "sorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 336 "sorcsd2by1.f"
	ib12d = ib11e + max(i__1,i__2);
#line 337 "sorcsd2by1.f"
	ib12e = ib12d + max(1,r__);
/* Computing MAX */
#line 338 "sorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 338 "sorcsd2by1.f"
	ib21d = ib12e + max(i__1,i__2);
#line 339 "sorcsd2by1.f"
	ib21e = ib21d + max(1,r__);
/* Computing MAX */
#line 340 "sorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 340 "sorcsd2by1.f"
	ib22d = ib21e + max(i__1,i__2);
#line 341 "sorcsd2by1.f"
	ib22e = ib22d + max(1,r__);
/* Computing MAX */
#line 342 "sorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 342 "sorcsd2by1.f"
	ibbcsd = ib22e + max(i__1,i__2);
/* Computing MAX */
#line 343 "sorcsd2by1.f"
	i__1 = 1, i__2 = r__ - 1;
#line 343 "sorcsd2by1.f"
	itaup1 = iphi + max(i__1,i__2);
#line 344 "sorcsd2by1.f"
	itaup2 = itaup1 + max(1,*p);
/* Computing MAX */
#line 345 "sorcsd2by1.f"
	i__1 = 1, i__2 = *m - *p;
#line 345 "sorcsd2by1.f"
	itauq1 = itaup2 + max(i__1,i__2);
#line 346 "sorcsd2by1.f"
	iorbdb = itauq1 + max(1,*q);
#line 347 "sorcsd2by1.f"
	iorgqr = itauq1 + max(1,*q);
#line 348 "sorcsd2by1.f"
	iorglq = itauq1 + max(1,*q);
#line 349 "sorcsd2by1.f"
	lorgqrmin = 1;
#line 350 "sorcsd2by1.f"
	lorgqropt = 1;
#line 351 "sorcsd2by1.f"
	lorglqmin = 1;
#line 352 "sorcsd2by1.f"
	lorglqopt = 1;
#line 353 "sorcsd2by1.f"
	if (r__ == *q) {
#line 354 "sorcsd2by1.f"
	    sorbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], dum1, dum1, dum1, dum1, &work[1], &c_n1,
		     &childinfo);
#line 357 "sorcsd2by1.f"
	    lorbdb = (integer) work[1];
#line 358 "sorcsd2by1.f"
	    if (wantu1 && *p > 0) {
#line 359 "sorcsd2by1.f"
		sorgqr_(p, p, q, &u1[u1_offset], ldu1, dum1, &work[1], &c_n1, 
			&childinfo);
#line 361 "sorcsd2by1.f"
		lorgqrmin = max(lorgqrmin,*p);
/* Computing MAX */
#line 362 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 362 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 363 "sorcsd2by1.f"
	    }
#line 364 "sorcsd2by1.f"
	    if (wantu2 && *m - *p > 0) {
#line 365 "sorcsd2by1.f"
		i__1 = *m - *p;
#line 365 "sorcsd2by1.f"
		i__2 = *m - *p;
#line 365 "sorcsd2by1.f"
		sorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, dum1, &work[1],
			 &c_n1, &childinfo);
/* Computing MAX */
#line 367 "sorcsd2by1.f"
		i__1 = lorgqrmin, i__2 = *m - *p;
#line 367 "sorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
/* Computing MAX */
#line 368 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 368 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 369 "sorcsd2by1.f"
	    }
#line 370 "sorcsd2by1.f"
	    if (wantv1t && *q > 0) {
#line 371 "sorcsd2by1.f"
		i__1 = *q - 1;
#line 371 "sorcsd2by1.f"
		i__2 = *q - 1;
#line 371 "sorcsd2by1.f"
		i__3 = *q - 1;
#line 371 "sorcsd2by1.f"
		sorglq_(&i__1, &i__2, &i__3, &v1t[v1t_offset], ldv1t, dum1, &
			work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 373 "sorcsd2by1.f"
		i__1 = lorglqmin, i__2 = *q - 1;
#line 373 "sorcsd2by1.f"
		lorglqmin = max(i__1,i__2);
/* Computing MAX */
#line 374 "sorcsd2by1.f"
		i__1 = lorglqopt, i__2 = (integer) work[1];
#line 374 "sorcsd2by1.f"
		lorglqopt = max(i__1,i__2);
#line 375 "sorcsd2by1.f"
	    }
#line 376 "sorcsd2by1.f"
	    sbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], dum1, 
		    &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		    v1t_offset], ldv1t, dum2, &c__1, dum1, dum1, dum1, dum1, 
		    dum1, dum1, dum1, dum1, &work[1], &c_n1, &childinfo, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 381 "sorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 382 "sorcsd2by1.f"
	} else if (r__ == *p) {
#line 383 "sorcsd2by1.f"
	    sorbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], dum1, dum1, dum1, dum1, &work[1], &c_n1,
		     &childinfo);
#line 386 "sorcsd2by1.f"
	    lorbdb = (integer) work[1];
#line 387 "sorcsd2by1.f"
	    if (wantu1 && *p > 0) {
#line 388 "sorcsd2by1.f"
		i__1 = *p - 1;
#line 388 "sorcsd2by1.f"
		i__2 = *p - 1;
#line 388 "sorcsd2by1.f"
		i__3 = *p - 1;
#line 388 "sorcsd2by1.f"
		sorgqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, 
			dum1, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 390 "sorcsd2by1.f"
		i__1 = lorgqrmin, i__2 = *p - 1;
#line 390 "sorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
/* Computing MAX */
#line 391 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 391 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 392 "sorcsd2by1.f"
	    }
#line 393 "sorcsd2by1.f"
	    if (wantu2 && *m - *p > 0) {
#line 394 "sorcsd2by1.f"
		i__1 = *m - *p;
#line 394 "sorcsd2by1.f"
		i__2 = *m - *p;
#line 394 "sorcsd2by1.f"
		sorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, dum1, &work[1],
			 &c_n1, &childinfo);
/* Computing MAX */
#line 396 "sorcsd2by1.f"
		i__1 = lorgqrmin, i__2 = *m - *p;
#line 396 "sorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
/* Computing MAX */
#line 397 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 397 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 398 "sorcsd2by1.f"
	    }
#line 399 "sorcsd2by1.f"
	    if (wantv1t && *q > 0) {
#line 400 "sorcsd2by1.f"
		sorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, dum1, &work[1], &
			c_n1, &childinfo);
#line 402 "sorcsd2by1.f"
		lorglqmin = max(lorglqmin,*q);
/* Computing MAX */
#line 403 "sorcsd2by1.f"
		i__1 = lorglqopt, i__2 = (integer) work[1];
#line 403 "sorcsd2by1.f"
		lorglqopt = max(i__1,i__2);
#line 404 "sorcsd2by1.f"
	    }
#line 405 "sorcsd2by1.f"
	    sbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], dum1, 
		    &v1t[v1t_offset], ldv1t, dum2, &c__1, &u1[u1_offset], 
		    ldu1, &u2[u2_offset], ldu2, dum1, dum1, dum1, dum1, dum1, 
		    dum1, dum1, dum1, &work[1], &c_n1, &childinfo, (ftnlen)1, 
		    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 410 "sorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 411 "sorcsd2by1.f"
	} else if (r__ == *m - *p) {
#line 412 "sorcsd2by1.f"
	    sorbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], dum1, dum1, dum1, dum1, &work[1], &c_n1,
		     &childinfo);
#line 415 "sorcsd2by1.f"
	    lorbdb = (integer) work[1];
#line 416 "sorcsd2by1.f"
	    if (wantu1 && *p > 0) {
#line 417 "sorcsd2by1.f"
		sorgqr_(p, p, q, &u1[u1_offset], ldu1, dum1, &work[1], &c_n1, 
			&childinfo);
#line 419 "sorcsd2by1.f"
		lorgqrmin = max(lorgqrmin,*p);
/* Computing MAX */
#line 420 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 420 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 421 "sorcsd2by1.f"
	    }
#line 422 "sorcsd2by1.f"
	    if (wantu2 && *m - *p > 0) {
#line 423 "sorcsd2by1.f"
		i__1 = *m - *p - 1;
#line 423 "sorcsd2by1.f"
		i__2 = *m - *p - 1;
#line 423 "sorcsd2by1.f"
		i__3 = *m - *p - 1;
#line 423 "sorcsd2by1.f"
		sorgqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, 
			dum1, &work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 425 "sorcsd2by1.f"
		i__1 = lorgqrmin, i__2 = *m - *p - 1;
#line 425 "sorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
/* Computing MAX */
#line 426 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 426 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 427 "sorcsd2by1.f"
	    }
#line 428 "sorcsd2by1.f"
	    if (wantv1t && *q > 0) {
#line 429 "sorcsd2by1.f"
		sorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, dum1, &work[1], &
			c_n1, &childinfo);
#line 431 "sorcsd2by1.f"
		lorglqmin = max(lorglqmin,*q);
/* Computing MAX */
#line 432 "sorcsd2by1.f"
		i__1 = lorglqopt, i__2 = (integer) work[1];
#line 432 "sorcsd2by1.f"
		lorglqopt = max(i__1,i__2);
#line 433 "sorcsd2by1.f"
	    }
#line 434 "sorcsd2by1.f"
	    i__1 = *m - *q;
#line 434 "sorcsd2by1.f"
	    i__2 = *m - *p;
#line 434 "sorcsd2by1.f"
	    sbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1]
		    , dum1, dum2, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		    u2_offset], ldu2, &u1[u1_offset], ldu1, dum1, dum1, dum1, 
		    dum1, dum1, dum1, dum1, dum1, &work[1], &c_n1, &childinfo,
		     (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 439 "sorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 440 "sorcsd2by1.f"
	} else {
#line 441 "sorcsd2by1.f"
	    sorbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], 
		    ldx21, &theta[1], dum1, dum1, dum1, dum1, dum1, &work[1], 
		    &c_n1, &childinfo);
#line 444 "sorcsd2by1.f"
	    lorbdb = *m + (integer) work[1];
#line 445 "sorcsd2by1.f"
	    if (wantu1 && *p > 0) {
#line 446 "sorcsd2by1.f"
		i__1 = *m - *q;
#line 446 "sorcsd2by1.f"
		sorgqr_(p, p, &i__1, &u1[u1_offset], ldu1, dum1, &work[1], &
			c_n1, &childinfo);
#line 448 "sorcsd2by1.f"
		lorgqrmin = max(lorgqrmin,*p);
/* Computing MAX */
#line 449 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 449 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 450 "sorcsd2by1.f"
	    }
#line 451 "sorcsd2by1.f"
	    if (wantu2 && *m - *p > 0) {
#line 452 "sorcsd2by1.f"
		i__1 = *m - *p;
#line 452 "sorcsd2by1.f"
		i__2 = *m - *p;
#line 452 "sorcsd2by1.f"
		i__3 = *m - *q;
#line 452 "sorcsd2by1.f"
		sorgqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, dum1, &
			work[1], &c_n1, &childinfo);
/* Computing MAX */
#line 454 "sorcsd2by1.f"
		i__1 = lorgqrmin, i__2 = *m - *p;
#line 454 "sorcsd2by1.f"
		lorgqrmin = max(i__1,i__2);
/* Computing MAX */
#line 455 "sorcsd2by1.f"
		i__1 = lorgqropt, i__2 = (integer) work[1];
#line 455 "sorcsd2by1.f"
		lorgqropt = max(i__1,i__2);
#line 456 "sorcsd2by1.f"
	    }
#line 457 "sorcsd2by1.f"
	    if (wantv1t && *q > 0) {
#line 458 "sorcsd2by1.f"
		sorglq_(q, q, q, &v1t[v1t_offset], ldv1t, dum1, &work[1], &
			c_n1, &childinfo);
#line 460 "sorcsd2by1.f"
		lorglqmin = max(lorglqmin,*q);
/* Computing MAX */
#line 461 "sorcsd2by1.f"
		i__1 = lorglqopt, i__2 = (integer) work[1];
#line 461 "sorcsd2by1.f"
		lorglqopt = max(i__1,i__2);
#line 462 "sorcsd2by1.f"
	    }
#line 463 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 463 "sorcsd2by1.f"
	    i__2 = *m - *q;
#line 463 "sorcsd2by1.f"
	    sbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1]
		    , dum1, &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, dum2, 
		    &c__1, &v1t[v1t_offset], ldv1t, dum1, dum1, dum1, dum1, 
		    dum1, dum1, dum1, dum1, &work[1], &c_n1, &childinfo, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 468 "sorcsd2by1.f"
	    lbbcsd = (integer) work[1];
#line 469 "sorcsd2by1.f"
	}
/* Computing MAX */
#line 470 "sorcsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqrmin - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqmin - 1, i__1 = max(i__1,
		i__2), i__2 = ibbcsd + lbbcsd - 1;
#line 470 "sorcsd2by1.f"
	lworkmin = max(i__1,i__2);
/* Computing MAX */
#line 474 "sorcsd2by1.f"
	i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqropt - 1, i__1 = max(
		i__1,i__2), i__2 = iorglq + lorglqopt - 1, i__1 = max(i__1,
		i__2), i__2 = ibbcsd + lbbcsd - 1;
#line 474 "sorcsd2by1.f"
	lworkopt = max(i__1,i__2);
#line 478 "sorcsd2by1.f"
	work[1] = (doublereal) lworkopt;
#line 479 "sorcsd2by1.f"
	if (*lwork < lworkmin && ! lquery) {
#line 480 "sorcsd2by1.f"
	    *info = -19;
#line 481 "sorcsd2by1.f"
	}
#line 482 "sorcsd2by1.f"
    }
#line 483 "sorcsd2by1.f"
    if (*info != 0) {
#line 484 "sorcsd2by1.f"
	i__1 = -(*info);
#line 484 "sorcsd2by1.f"
	xerbla_("SORCSD2BY1", &i__1, (ftnlen)10);
#line 485 "sorcsd2by1.f"
	return 0;
#line 486 "sorcsd2by1.f"
    } else if (lquery) {
#line 487 "sorcsd2by1.f"
	return 0;
#line 488 "sorcsd2by1.f"
    }
#line 489 "sorcsd2by1.f"
    lorgqr = *lwork - iorgqr + 1;
#line 490 "sorcsd2by1.f"
    lorglq = *lwork - iorglq + 1;

/*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q, */
/*     in which R = MIN(P,M-P,Q,M-Q) */

#line 495 "sorcsd2by1.f"
    if (r__ == *q) {

/*        Case 1: R = Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 501 "sorcsd2by1.f"
	sorbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 507 "sorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 508 "sorcsd2by1.f"
	    slacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 509 "sorcsd2by1.f"
	    sorgqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 511 "sorcsd2by1.f"
	}
#line 512 "sorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 513 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 513 "sorcsd2by1.f"
	    slacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 514 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 514 "sorcsd2by1.f"
	    i__2 = *m - *p;
#line 514 "sorcsd2by1.f"
	    sorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 516 "sorcsd2by1.f"
	}
#line 517 "sorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 518 "sorcsd2by1.f"
	    v1t[v1t_dim1 + 1] = 1.;
#line 519 "sorcsd2by1.f"
	    i__1 = *q;
#line 519 "sorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 520 "sorcsd2by1.f"
		v1t[j * v1t_dim1 + 1] = 0.;
#line 521 "sorcsd2by1.f"
		v1t[j + v1t_dim1] = 0.;
#line 522 "sorcsd2by1.f"
	    }
#line 523 "sorcsd2by1.f"
	    i__1 = *q - 1;
#line 523 "sorcsd2by1.f"
	    i__2 = *q - 1;
#line 523 "sorcsd2by1.f"
	    slacpy_("U", &i__1, &i__2, &x21[(x21_dim1 << 1) + 1], ldx21, &v1t[
		    (v1t_dim1 << 1) + 2], ldv1t, (ftnlen)1);
#line 525 "sorcsd2by1.f"
	    i__1 = *q - 1;
#line 525 "sorcsd2by1.f"
	    i__2 = *q - 1;
#line 525 "sorcsd2by1.f"
	    i__3 = *q - 1;
#line 525 "sorcsd2by1.f"
	    sorglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, &
		    work[itauq1], &work[iorglq], &lorglq, &childinfo);
#line 527 "sorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 531 "sorcsd2by1.f"
	sbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &work[
		iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[
		v1t_offset], ldv1t, dum2, &c__1, &work[ib11d], &work[ib11e], &
		work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e], &work[
		ib22d], &work[ib22e], &work[ibbcsd], &lbbcsd, &childinfo, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Permute rows and columns to place zero submatrices in */
/*        preferred positions */

#line 541 "sorcsd2by1.f"
	if (*q > 0 && wantu2) {
#line 542 "sorcsd2by1.f"
	    i__1 = *q;
#line 542 "sorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 543 "sorcsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 544 "sorcsd2by1.f"
	    }
#line 545 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 545 "sorcsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 546 "sorcsd2by1.f"
		iwork[i__] = i__ - *q;
#line 547 "sorcsd2by1.f"
	    }
#line 548 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 548 "sorcsd2by1.f"
	    i__2 = *m - *p;
#line 548 "sorcsd2by1.f"
	    slapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 549 "sorcsd2by1.f"
	}
#line 550 "sorcsd2by1.f"
    } else if (r__ == *p) {

/*        Case 2: R = P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 556 "sorcsd2by1.f"
	sorbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 562 "sorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 563 "sorcsd2by1.f"
	    u1[u1_dim1 + 1] = 1.;
#line 564 "sorcsd2by1.f"
	    i__1 = *p;
#line 564 "sorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 565 "sorcsd2by1.f"
		u1[j * u1_dim1 + 1] = 0.;
#line 566 "sorcsd2by1.f"
		u1[j + u1_dim1] = 0.;
#line 567 "sorcsd2by1.f"
	    }
#line 568 "sorcsd2by1.f"
	    i__1 = *p - 1;
#line 568 "sorcsd2by1.f"
	    i__2 = *p - 1;
#line 568 "sorcsd2by1.f"
	    slacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 569 "sorcsd2by1.f"
	    i__1 = *p - 1;
#line 569 "sorcsd2by1.f"
	    i__2 = *p - 1;
#line 569 "sorcsd2by1.f"
	    i__3 = *p - 1;
#line 569 "sorcsd2by1.f"
	    sorgqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &work[
		    itaup1], &work[iorgqr], &lorgqr, &childinfo);
#line 571 "sorcsd2by1.f"
	}
#line 572 "sorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 573 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 573 "sorcsd2by1.f"
	    slacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], 
		    ldu2, (ftnlen)1);
#line 574 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 574 "sorcsd2by1.f"
	    i__2 = *m - *p;
#line 574 "sorcsd2by1.f"
	    sorgqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], &
		    work[iorgqr], &lorgqr, &childinfo);
#line 576 "sorcsd2by1.f"
	}
#line 577 "sorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 578 "sorcsd2by1.f"
	    slacpy_("U", p, q, &x11[x11_offset], ldx11, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 579 "sorcsd2by1.f"
	    sorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 581 "sorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 585 "sorcsd2by1.f"
	sbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &work[
		iphi], &v1t[v1t_offset], ldv1t, dum1, &c__1, &u1[u1_offset], 
		ldu1, &u2[u2_offset], ldu2, &work[ib11d], &work[ib11e], &work[
		ib12d], &work[ib12e], &work[ib21d], &work[ib21e], &work[ib22d]
		, &work[ib22e], &work[ibbcsd], &lbbcsd, &childinfo, (ftnlen)1,
		 (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 595 "sorcsd2by1.f"
	if (*q > 0 && wantu2) {
#line 596 "sorcsd2by1.f"
	    i__1 = *q;
#line 596 "sorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 597 "sorcsd2by1.f"
		iwork[i__] = *m - *p - *q + i__;
#line 598 "sorcsd2by1.f"
	    }
#line 599 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 599 "sorcsd2by1.f"
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
#line 600 "sorcsd2by1.f"
		iwork[i__] = i__ - *q;
#line 601 "sorcsd2by1.f"
	    }
#line 602 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 602 "sorcsd2by1.f"
	    i__2 = *m - *p;
#line 602 "sorcsd2by1.f"
	    slapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
#line 603 "sorcsd2by1.f"
	}
#line 604 "sorcsd2by1.f"
    } else if (r__ == *m - *p) {

/*        Case 3: R = M-P */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 610 "sorcsd2by1.f"
	sorbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &lorbdb, &childinfo);

/*        Accumulate Householder reflectors */

#line 616 "sorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 617 "sorcsd2by1.f"
	    slacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1, 
		    (ftnlen)1);
#line 618 "sorcsd2by1.f"
	    sorgqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 620 "sorcsd2by1.f"
	}
#line 621 "sorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 622 "sorcsd2by1.f"
	    u2[u2_dim1 + 1] = 1.;
#line 623 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 623 "sorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 624 "sorcsd2by1.f"
		u2[j * u2_dim1 + 1] = 0.;
#line 625 "sorcsd2by1.f"
		u2[j + u2_dim1] = 0.;
#line 626 "sorcsd2by1.f"
	    }
#line 627 "sorcsd2by1.f"
	    i__1 = *m - *p - 1;
#line 627 "sorcsd2by1.f"
	    i__2 = *m - *p - 1;
#line 627 "sorcsd2by1.f"
	    slacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 629 "sorcsd2by1.f"
	    i__1 = *m - *p - 1;
#line 629 "sorcsd2by1.f"
	    i__2 = *m - *p - 1;
#line 629 "sorcsd2by1.f"
	    i__3 = *m - *p - 1;
#line 629 "sorcsd2by1.f"
	    sorgqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &work[
		    itaup2], &work[iorgqr], &lorgqr, &childinfo);
#line 631 "sorcsd2by1.f"
	}
#line 632 "sorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 633 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 633 "sorcsd2by1.f"
	    slacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 634 "sorcsd2by1.f"
	    sorglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 636 "sorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 640 "sorcsd2by1.f"
	i__1 = *m - *q;
#line 640 "sorcsd2by1.f"
	i__2 = *m - *p;
#line 640 "sorcsd2by1.f"
	sbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1], &
		work[iphi], dum1, &c__1, &v1t[v1t_offset], ldv1t, &u2[
		u2_offset], ldu2, &u1[u1_offset], ldu1, &work[ib11d], &work[
		ib11e], &work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e]
		, &work[ib22d], &work[ib22e], &work[ibbcsd], &lbbcsd, &
		childinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 650 "sorcsd2by1.f"
	if (*q > r__) {
#line 651 "sorcsd2by1.f"
	    i__1 = r__;
#line 651 "sorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 652 "sorcsd2by1.f"
		iwork[i__] = *q - r__ + i__;
#line 653 "sorcsd2by1.f"
	    }
#line 654 "sorcsd2by1.f"
	    i__1 = *q;
#line 654 "sorcsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 655 "sorcsd2by1.f"
		iwork[i__] = i__ - r__;
#line 656 "sorcsd2by1.f"
	    }
#line 657 "sorcsd2by1.f"
	    if (wantu1) {
#line 658 "sorcsd2by1.f"
		slapmt_(&c_false, p, q, &u1[u1_offset], ldu1, &iwork[1]);
#line 659 "sorcsd2by1.f"
	    }
#line 660 "sorcsd2by1.f"
	    if (wantv1t) {
#line 661 "sorcsd2by1.f"
		slapmr_(&c_false, q, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 662 "sorcsd2by1.f"
	    }
#line 663 "sorcsd2by1.f"
	}
#line 664 "sorcsd2by1.f"
    } else {

/*        Case 4: R = M-Q */

/*        Simultaneously bidiagonalize X11 and X21 */

#line 670 "sorcsd2by1.f"
	i__1 = lorbdb - *m;
#line 670 "sorcsd2by1.f"
	sorbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &
		theta[1], &work[iphi], &work[itaup1], &work[itaup2], &work[
		itauq1], &work[iorbdb], &work[iorbdb + *m], &i__1, &childinfo)
		;

/*        Accumulate Householder reflectors */

#line 677 "sorcsd2by1.f"
	if (wantu1 && *p > 0) {
#line 678 "sorcsd2by1.f"
	    scopy_(p, &work[iorbdb], &c__1, &u1[u1_offset], &c__1);
#line 679 "sorcsd2by1.f"
	    i__1 = *p;
#line 679 "sorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 680 "sorcsd2by1.f"
		u1[j * u1_dim1 + 1] = 0.;
#line 681 "sorcsd2by1.f"
	    }
#line 682 "sorcsd2by1.f"
	    i__1 = *p - 1;
#line 682 "sorcsd2by1.f"
	    i__2 = *m - *q - 1;
#line 682 "sorcsd2by1.f"
	    slacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[(
		    u1_dim1 << 1) + 2], ldu1, (ftnlen)1);
#line 684 "sorcsd2by1.f"
	    i__1 = *m - *q;
#line 684 "sorcsd2by1.f"
	    sorgqr_(p, p, &i__1, &u1[u1_offset], ldu1, &work[itaup1], &work[
		    iorgqr], &lorgqr, &childinfo);
#line 686 "sorcsd2by1.f"
	}
#line 687 "sorcsd2by1.f"
	if (wantu2 && *m - *p > 0) {
#line 688 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 688 "sorcsd2by1.f"
	    scopy_(&i__1, &work[iorbdb + *p], &c__1, &u2[u2_offset], &c__1);
#line 689 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 689 "sorcsd2by1.f"
	    for (j = 2; j <= i__1; ++j) {
#line 690 "sorcsd2by1.f"
		u2[j * u2_dim1 + 1] = 0.;
#line 691 "sorcsd2by1.f"
	    }
#line 692 "sorcsd2by1.f"
	    i__1 = *m - *p - 1;
#line 692 "sorcsd2by1.f"
	    i__2 = *m - *q - 1;
#line 692 "sorcsd2by1.f"
	    slacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[(
		    u2_dim1 << 1) + 2], ldu2, (ftnlen)1);
#line 694 "sorcsd2by1.f"
	    i__1 = *m - *p;
#line 694 "sorcsd2by1.f"
	    i__2 = *m - *p;
#line 694 "sorcsd2by1.f"
	    i__3 = *m - *q;
#line 694 "sorcsd2by1.f"
	    sorgqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &work[itaup2], 
		    &work[iorgqr], &lorgqr, &childinfo);
#line 696 "sorcsd2by1.f"
	}
#line 697 "sorcsd2by1.f"
	if (wantv1t && *q > 0) {
#line 698 "sorcsd2by1.f"
	    i__1 = *m - *q;
#line 698 "sorcsd2by1.f"
	    slacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], 
		    ldv1t, (ftnlen)1);
#line 699 "sorcsd2by1.f"
	    i__1 = *p - (*m - *q);
#line 699 "sorcsd2by1.f"
	    i__2 = *q - (*m - *q);
#line 699 "sorcsd2by1.f"
	    slacpy_("U", &i__1, &i__2, &x11[*m - *q + 1 + (*m - *q + 1) * 
		    x11_dim1], ldx11, &v1t[*m - *q + 1 + (*m - *q + 1) * 
		    v1t_dim1], ldv1t, (ftnlen)1);
#line 701 "sorcsd2by1.f"
	    i__1 = -(*p) + *q;
#line 701 "sorcsd2by1.f"
	    i__2 = *q - *p;
#line 701 "sorcsd2by1.f"
	    slacpy_("U", &i__1, &i__2, &x21[*m - *q + 1 + (*p + 1) * x21_dim1]
		    , ldx21, &v1t[*p + 1 + (*p + 1) * v1t_dim1], ldv1t, (
		    ftnlen)1);
#line 703 "sorcsd2by1.f"
	    sorglq_(q, q, q, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[
		    iorglq], &lorglq, &childinfo);
#line 705 "sorcsd2by1.f"
	}

/*        Simultaneously diagonalize X11 and X21. */

#line 709 "sorcsd2by1.f"
	i__1 = *m - *p;
#line 709 "sorcsd2by1.f"
	i__2 = *m - *q;
#line 709 "sorcsd2by1.f"
	sbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1], &
		work[iphi], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, dum1, 
		&c__1, &v1t[v1t_offset], ldv1t, &work[ib11d], &work[ib11e], &
		work[ib12d], &work[ib12e], &work[ib21d], &work[ib21e], &work[
		ib22d], &work[ib22e], &work[ibbcsd], &lbbcsd, &childinfo, (
		ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Permute rows and columns to place identity submatrices in */
/*        preferred positions */

#line 719 "sorcsd2by1.f"
	if (*p > r__) {
#line 720 "sorcsd2by1.f"
	    i__1 = r__;
#line 720 "sorcsd2by1.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 721 "sorcsd2by1.f"
		iwork[i__] = *p - r__ + i__;
#line 722 "sorcsd2by1.f"
	    }
#line 723 "sorcsd2by1.f"
	    i__1 = *p;
#line 723 "sorcsd2by1.f"
	    for (i__ = r__ + 1; i__ <= i__1; ++i__) {
#line 724 "sorcsd2by1.f"
		iwork[i__] = i__ - r__;
#line 725 "sorcsd2by1.f"
	    }
#line 726 "sorcsd2by1.f"
	    if (wantu1) {
#line 727 "sorcsd2by1.f"
		slapmt_(&c_false, p, p, &u1[u1_offset], ldu1, &iwork[1]);
#line 728 "sorcsd2by1.f"
	    }
#line 729 "sorcsd2by1.f"
	    if (wantv1t) {
#line 730 "sorcsd2by1.f"
		slapmr_(&c_false, p, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
#line 731 "sorcsd2by1.f"
	    }
#line 732 "sorcsd2by1.f"
	}
#line 733 "sorcsd2by1.f"
    }

#line 735 "sorcsd2by1.f"
    return 0;

/*     End of SORCSD2BY1 */

} /* sorcsd2by1_ */

