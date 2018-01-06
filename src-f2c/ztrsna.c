#line 1 "ztrsna.f"
/* ztrsna.f -- translated by f2c (version 20100827).
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

#line 1 "ztrsna.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTRSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTRSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                          LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       DOUBLE PRECISION   RWORK( * ), S( * ), SEP( * ) */
/*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or right eigenvectors of a complex upper triangular */
/* > matrix T (or of any matrix Q*T*Q**H with Q unitary). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies whether condition numbers are required for */
/* >          eigenvalues (S) or eigenvectors (SEP): */
/* >          = 'E': for eigenvalues only (S); */
/* >          = 'V': for eigenvectors only (SEP); */
/* >          = 'B': for both eigenvalues and eigenvectors (S and SEP). */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* >          HOWMNY is CHARACTER*1 */
/* >          = 'A': compute condition numbers for all eigenpairs; */
/* >          = 'S': compute condition numbers for selected eigenpairs */
/* >                 specified by the array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY = 'S', SELECT specifies the eigenpairs for which */
/* >          condition numbers are required. To select condition numbers */
/* >          for the j-th eigenpair, SELECT(j) must be set to .TRUE.. */
/* >          If HOWMNY = 'A', SELECT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,N) */
/* >          The upper triangular matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,M) */
/* >          If JOB = 'E' or 'B', VL must contain left eigenvectors of T */
/* >          (or of any Q*T*Q**H with Q unitary), corresponding to the */
/* >          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* >          must be stored in consecutive columns of VL, as returned by */
/* >          ZHSEIN or ZTREVC. */
/* >          If JOB = 'V', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL. */
/* >          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] VR */
/* > \verbatim */
/* >          VR is COMPLEX*16 array, dimension (LDVR,M) */
/* >          If JOB = 'E' or 'B', VR must contain right eigenvectors of T */
/* >          (or of any Q*T*Q**H with Q unitary), corresponding to the */
/* >          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* >          must be stored in consecutive columns of VR, as returned by */
/* >          ZHSEIN or ZTREVC. */
/* >          If JOB = 'V', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR. */
/* >          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (MM) */
/* >          If JOB = 'E' or 'B', the reciprocal condition numbers of the */
/* >          selected eigenvalues, stored in consecutive elements of the */
/* >          array. Thus S(j), SEP(j), and the j-th columns of VL and VR */
/* >          all correspond to the same eigenpair (but not in general the */
/* >          j-th eigenpair, unless all eigenpairs are selected). */
/* >          If JOB = 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SEP */
/* > \verbatim */
/* >          SEP is DOUBLE PRECISION array, dimension (MM) */
/* >          If JOB = 'V' or 'B', the estimated reciprocal condition */
/* >          numbers of the selected eigenvectors, stored in consecutive */
/* >          elements of the array. */
/* >          If JOB = 'E', SEP is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* >          MM is INTEGER */
/* >          The number of elements in the arrays S (if JOB = 'E' or 'B') */
/* >           and/or SEP (if JOB = 'V' or 'B'). MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of elements of the arrays S and/or SEP actually */
/* >          used to store the estimated condition numbers. */
/* >          If HOWMNY = 'A', M is set to N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (LDWORK,N+6) */
/* >          If JOB = 'E', WORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* >          LDWORK is INTEGER */
/* >          The leading dimension of the array WORK. */
/* >          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
/* >          If JOB = 'E', RWORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The reciprocal of the condition number of an eigenvalue lambda is */
/* >  defined as */
/* > */
/* >          S(lambda) = |v**H*u| / (norm(u)*norm(v)) */
/* > */
/* >  where u and v are the right and left eigenvectors of T corresponding */
/* >  to lambda; v**H denotes the conjugate transpose of v, and norm(u) */
/* >  denotes the Euclidean norm. These reciprocal condition numbers always */
/* >  lie between zero (very badly conditioned) and one (very well */
/* >  conditioned). If n = 1, S(lambda) is defined to be 1. */
/* > */
/* >  An approximate error bound for a computed eigenvalue W(i) is given by */
/* > */
/* >                      EPS * norm(T) / S(i) */
/* > */
/* >  where EPS is the machine precision. */
/* > */
/* >  The reciprocal of the condition number of the right eigenvector u */
/* >  corresponding to lambda is defined as follows. Suppose */
/* > */
/* >              T = ( lambda  c  ) */
/* >                  (   0    T22 ) */
/* > */
/* >  Then the reciprocal condition number is */
/* > */
/* >          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I ) */
/* > */
/* >  where sigma-min denotes the smallest singular value. We approximate */
/* >  the smallest singular value by the reciprocal of an estimate of the */
/* >  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is */
/* >  defined to be abs(T(1,1)). */
/* > */
/* >  An approximate error bound for a computed right eigenvector VR(i) */
/* >  is given by */
/* > */
/* >                      EPS * norm(T) / SEP(i) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztrsna_(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, 
	doublereal *sep, integer *mm, integer *m, doublecomplex *work, 
	integer *ldwork, doublereal *rwork, integer *info, ftnlen job_len, 
	ftnlen howmny_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, 
	    work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, k, ks, ix;
    static doublereal eps, est;
    static integer kase, ierr;
    static doublecomplex prod;
    static doublereal lnrm, rnrm, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublecomplex dummy[1];
    static logical wants;
    static doublereal xnorm;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical wantbh;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical somcon;
    extern /* Subroutine */ int zdrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static char normin[1];
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal smlnum;
    static logical wantsp;
    extern /* Subroutine */ int zlatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), ztrexc_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

#line 310 "ztrsna.f"
    /* Parameter adjustments */
#line 310 "ztrsna.f"
    --select;
#line 310 "ztrsna.f"
    t_dim1 = *ldt;
#line 310 "ztrsna.f"
    t_offset = 1 + t_dim1;
#line 310 "ztrsna.f"
    t -= t_offset;
#line 310 "ztrsna.f"
    vl_dim1 = *ldvl;
#line 310 "ztrsna.f"
    vl_offset = 1 + vl_dim1;
#line 310 "ztrsna.f"
    vl -= vl_offset;
#line 310 "ztrsna.f"
    vr_dim1 = *ldvr;
#line 310 "ztrsna.f"
    vr_offset = 1 + vr_dim1;
#line 310 "ztrsna.f"
    vr -= vr_offset;
#line 310 "ztrsna.f"
    --s;
#line 310 "ztrsna.f"
    --sep;
#line 310 "ztrsna.f"
    work_dim1 = *ldwork;
#line 310 "ztrsna.f"
    work_offset = 1 + work_dim1;
#line 310 "ztrsna.f"
    work -= work_offset;
#line 310 "ztrsna.f"
    --rwork;
#line 310 "ztrsna.f"

#line 310 "ztrsna.f"
    /* Function Body */
#line 310 "ztrsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 311 "ztrsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 312 "ztrsna.f"
    wantsp = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 314 "ztrsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of eigenpairs for which condition numbers are */
/*     to be computed. */

#line 319 "ztrsna.f"
    if (somcon) {
#line 320 "ztrsna.f"
	*m = 0;
#line 321 "ztrsna.f"
	i__1 = *n;
#line 321 "ztrsna.f"
	for (j = 1; j <= i__1; ++j) {
#line 322 "ztrsna.f"
	    if (select[j]) {
#line 322 "ztrsna.f"
		++(*m);
#line 322 "ztrsna.f"
	    }
#line 324 "ztrsna.f"
/* L10: */
#line 324 "ztrsna.f"
	}
#line 325 "ztrsna.f"
    } else {
#line 326 "ztrsna.f"
	*m = *n;
#line 327 "ztrsna.f"
    }

#line 329 "ztrsna.f"
    *info = 0;
#line 330 "ztrsna.f"
    if (! wants && ! wantsp) {
#line 331 "ztrsna.f"
	*info = -1;
#line 332 "ztrsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 333 "ztrsna.f"
	*info = -2;
#line 334 "ztrsna.f"
    } else if (*n < 0) {
#line 335 "ztrsna.f"
	*info = -4;
#line 336 "ztrsna.f"
    } else if (*ldt < max(1,*n)) {
#line 337 "ztrsna.f"
	*info = -6;
#line 338 "ztrsna.f"
    } else if (*ldvl < 1 || wants && *ldvl < *n) {
#line 339 "ztrsna.f"
	*info = -8;
#line 340 "ztrsna.f"
    } else if (*ldvr < 1 || wants && *ldvr < *n) {
#line 341 "ztrsna.f"
	*info = -10;
#line 342 "ztrsna.f"
    } else if (*mm < *m) {
#line 343 "ztrsna.f"
	*info = -13;
#line 344 "ztrsna.f"
    } else if (*ldwork < 1 || wantsp && *ldwork < *n) {
#line 345 "ztrsna.f"
	*info = -16;
#line 346 "ztrsna.f"
    }
#line 347 "ztrsna.f"
    if (*info != 0) {
#line 348 "ztrsna.f"
	i__1 = -(*info);
#line 348 "ztrsna.f"
	xerbla_("ZTRSNA", &i__1, (ftnlen)6);
#line 349 "ztrsna.f"
	return 0;
#line 350 "ztrsna.f"
    }

/*     Quick return if possible */

#line 354 "ztrsna.f"
    if (*n == 0) {
#line 354 "ztrsna.f"
	return 0;
#line 354 "ztrsna.f"
    }

#line 357 "ztrsna.f"
    if (*n == 1) {
#line 358 "ztrsna.f"
	if (somcon) {
#line 359 "ztrsna.f"
	    if (! select[1]) {
#line 359 "ztrsna.f"
		return 0;
#line 359 "ztrsna.f"
	    }
#line 361 "ztrsna.f"
	}
#line 362 "ztrsna.f"
	if (wants) {
#line 362 "ztrsna.f"
	    s[1] = 1.;
#line 362 "ztrsna.f"
	}
#line 364 "ztrsna.f"
	if (wantsp) {
#line 364 "ztrsna.f"
	    sep[1] = z_abs(&t[t_dim1 + 1]);
#line 364 "ztrsna.f"
	}
#line 366 "ztrsna.f"
	return 0;
#line 367 "ztrsna.f"
    }

/*     Get machine constants */

#line 371 "ztrsna.f"
    eps = dlamch_("P", (ftnlen)1);
#line 372 "ztrsna.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 373 "ztrsna.f"
    bignum = 1. / smlnum;
#line 374 "ztrsna.f"
    dlabad_(&smlnum, &bignum);

#line 376 "ztrsna.f"
    ks = 1;
#line 377 "ztrsna.f"
    i__1 = *n;
#line 377 "ztrsna.f"
    for (k = 1; k <= i__1; ++k) {

#line 379 "ztrsna.f"
	if (somcon) {
#line 380 "ztrsna.f"
	    if (! select[k]) {
#line 380 "ztrsna.f"
		goto L50;
#line 380 "ztrsna.f"
	    }
#line 382 "ztrsna.f"
	}

#line 384 "ztrsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 389 "ztrsna.f"
	    zdotc_(&z__1, n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * vl_dim1 + 
		    1], &c__1);
#line 389 "ztrsna.f"
	    prod.r = z__1.r, prod.i = z__1.i;
#line 390 "ztrsna.f"
	    rnrm = dznrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 391 "ztrsna.f"
	    lnrm = dznrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 392 "ztrsna.f"
	    s[ks] = z_abs(&prod) / (rnrm * lnrm);

#line 394 "ztrsna.f"
	}

#line 396 "ztrsna.f"
	if (wantsp) {

/*           Estimate the reciprocal condition number of the k-th */
/*           eigenvector. */

/*           Copy the matrix T to the array WORK and swap the k-th */
/*           diagonal element to the (1,1) position. */

#line 404 "ztrsna.f"
	    zlacpy_("Full", n, n, &t[t_offset], ldt, &work[work_offset], 
		    ldwork, (ftnlen)4);
#line 405 "ztrsna.f"
	    ztrexc_("No Q", n, &work[work_offset], ldwork, dummy, &c__1, &k, &
		    c__1, &ierr, (ftnlen)4);

/*           Form  C = T22 - lambda*I in WORK(2:N,2:N). */

#line 409 "ztrsna.f"
	    i__2 = *n;
#line 409 "ztrsna.f"
	    for (i__ = 2; i__ <= i__2; ++i__) {
#line 410 "ztrsna.f"
		i__3 = i__ + i__ * work_dim1;
#line 410 "ztrsna.f"
		i__4 = i__ + i__ * work_dim1;
#line 410 "ztrsna.f"
		i__5 = work_dim1 + 1;
#line 410 "ztrsna.f"
		z__1.r = work[i__4].r - work[i__5].r, z__1.i = work[i__4].i - 
			work[i__5].i;
#line 410 "ztrsna.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 411 "ztrsna.f"
/* L20: */
#line 411 "ztrsna.f"
	    }

/*           Estimate a lower bound for the 1-norm of inv(C**H). The 1st */
/*           and (N+1)th columns of WORK are used to store work vectors. */

#line 416 "ztrsna.f"
	    sep[ks] = 0.;
#line 417 "ztrsna.f"
	    est = 0.;
#line 418 "ztrsna.f"
	    kase = 0;
#line 419 "ztrsna.f"
	    *(unsigned char *)normin = 'N';
#line 420 "ztrsna.f"
L30:
#line 421 "ztrsna.f"
	    i__2 = *n - 1;
#line 421 "ztrsna.f"
	    zlacn2_(&i__2, &work[(*n + 1) * work_dim1 + 1], &work[work_offset]
		    , &est, &kase, isave);

#line 423 "ztrsna.f"
	    if (kase != 0) {
#line 424 "ztrsna.f"
		if (kase == 1) {

/*                 Solve C**H*x = scale*b */

#line 428 "ztrsna.f"
		    i__2 = *n - 1;
#line 428 "ztrsna.f"
		    zlatrs_("Upper", "Conjugate transpose", "Nonunit", normin,
			     &i__2, &work[(work_dim1 << 1) + 2], ldwork, &
			    work[work_offset], &scale, &rwork[1], &ierr, (
			    ftnlen)5, (ftnlen)19, (ftnlen)7, (ftnlen)1);
#line 431 "ztrsna.f"
		} else {

/*                 Solve C*x = scale*b */

#line 435 "ztrsna.f"
		    i__2 = *n - 1;
#line 435 "ztrsna.f"
		    zlatrs_("Upper", "No transpose", "Nonunit", normin, &i__2,
			     &work[(work_dim1 << 1) + 2], ldwork, &work[
			    work_offset], &scale, &rwork[1], &ierr, (ftnlen)5,
			     (ftnlen)12, (ftnlen)7, (ftnlen)1);
#line 438 "ztrsna.f"
		}
#line 439 "ztrsna.f"
		*(unsigned char *)normin = 'Y';
#line 440 "ztrsna.f"
		if (scale != 1.) {

/*                 Multiply by 1/SCALE if doing so will not cause */
/*                 overflow. */

#line 445 "ztrsna.f"
		    i__2 = *n - 1;
#line 445 "ztrsna.f"
		    ix = izamax_(&i__2, &work[work_offset], &c__1);
#line 446 "ztrsna.f"
		    i__2 = ix + work_dim1;
#line 446 "ztrsna.f"
		    xnorm = (d__1 = work[i__2].r, abs(d__1)) + (d__2 = d_imag(
			    &work[ix + work_dim1]), abs(d__2));
#line 447 "ztrsna.f"
		    if (scale < xnorm * smlnum || scale == 0.) {
#line 447 "ztrsna.f"
			goto L40;
#line 447 "ztrsna.f"
		    }
#line 449 "ztrsna.f"
		    zdrscl_(n, &scale, &work[work_offset], &c__1);
#line 450 "ztrsna.f"
		}
#line 451 "ztrsna.f"
		goto L30;
#line 452 "ztrsna.f"
	    }

#line 454 "ztrsna.f"
	    sep[ks] = 1. / max(est,smlnum);
#line 455 "ztrsna.f"
	}

#line 457 "ztrsna.f"
L40:
#line 458 "ztrsna.f"
	++ks;
#line 459 "ztrsna.f"
L50:
#line 459 "ztrsna.f"
	;
#line 459 "ztrsna.f"
    }
#line 460 "ztrsna.f"
    return 0;

/*     End of ZTRSNA */

} /* ztrsna_ */

