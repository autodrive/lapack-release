#line 1 "strsna.f"
/* strsna.f -- translated by f2c (version 20100827).
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

#line 1 "strsna.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_true = TRUE_;
static logical c_false = FALSE_;

/* > \brief \b STRSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                          LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       REAL               S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or right eigenvectors of a real upper */
/* > quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q */
/* > orthogonal). */
/* > */
/* > T must be in Schur canonical form (as returned by SHSEQR), that is, */
/* > block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
/* > 2-by-2 diagonal block has its diagonal elements equal and its */
/* > off-diagonal elements of opposite sign. */
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
/* >          for the eigenpair corresponding to a real eigenvalue w(j), */
/* >          SELECT(j) must be set to .TRUE.. To select condition numbers */
/* >          corresponding to a complex conjugate pair of eigenvalues w(j) */
/* >          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be */
/* >          set to .TRUE.. */
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
/* >          T is REAL array, dimension (LDT,N) */
/* >          The upper quasi-triangular matrix T, in Schur canonical form. */
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
/* >          VL is REAL array, dimension (LDVL,M) */
/* >          If JOB = 'E' or 'B', VL must contain left eigenvectors of T */
/* >          (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
/* >          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* >          must be stored in consecutive columns of VL, as returned by */
/* >          SHSEIN or STREVC. */
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
/* >          VR is REAL array, dimension (LDVR,M) */
/* >          If JOB = 'E' or 'B', VR must contain right eigenvectors of T */
/* >          (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
/* >          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* >          must be stored in consecutive columns of VR, as returned by */
/* >          SHSEIN or STREVC. */
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
/* >          S is REAL array, dimension (MM) */
/* >          If JOB = 'E' or 'B', the reciprocal condition numbers of the */
/* >          selected eigenvalues, stored in consecutive elements of the */
/* >          array. For a complex conjugate pair of eigenvalues two */
/* >          consecutive elements of S are set to the same value. Thus */
/* >          S(j), SEP(j), and the j-th columns of VL and VR all */
/* >          correspond to the same eigenpair (but not in general the */
/* >          j-th eigenpair, unless all eigenpairs are selected). */
/* >          If JOB = 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SEP */
/* > \verbatim */
/* >          SEP is REAL array, dimension (MM) */
/* >          If JOB = 'V' or 'B', the estimated reciprocal condition */
/* >          numbers of the selected eigenvectors, stored in consecutive */
/* >          elements of the array. For a complex eigenvector two */
/* >          consecutive elements of SEP are set to the same value. If */
/* >          the eigenvalues cannot be reordered to compute SEP(j), SEP(j) */
/* >          is set to 0; this can only occur when the true value would be */
/* >          very small anyway. */
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
/* >          WORK is REAL array, dimension (LDWORK,N+6) */
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
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (2*(N-1)) */
/* >          If JOB = 'E', IWORK is not referenced. */
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

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The reciprocal of the condition number of an eigenvalue lambda is */
/* >  defined as */
/* > */
/* >          S(lambda) = |v**T*u| / (norm(u)*norm(v)) */
/* > */
/* >  where u and v are the right and left eigenvectors of T corresponding */
/* >  to lambda; v**T denotes the transpose of v, and norm(u) */
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
/* Subroutine */ int strsna_(char *job, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep, 
	integer *mm, integer *m, doublereal *work, integer *ldwork, integer *
	iwork, integer *info, ftnlen job_len, ftnlen howmny_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, 
	    work_dim1, work_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, n2;
    static doublereal cs;
    static integer nn, ks;
    static doublereal sn, mu, eps, est;
    static integer kase;
    static doublereal cond;
    static logical pair;
    static integer ierr;
    static doublereal dumm, prod;
    static integer ifst;
    static doublereal lnrm;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ilst;
    static doublereal rnrm, prod1, prod2;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal scale, delta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical wants;
    static doublereal dummy[1];
    extern /* Subroutine */ int slacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal slapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical wantbh;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static logical somcon;
    extern /* Subroutine */ int slaqtr_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *), strexc_(char *, integer *
	    , doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static doublereal smlnum;
    static logical wantsp;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

#line 315 "strsna.f"
    /* Parameter adjustments */
#line 315 "strsna.f"
    --select;
#line 315 "strsna.f"
    t_dim1 = *ldt;
#line 315 "strsna.f"
    t_offset = 1 + t_dim1;
#line 315 "strsna.f"
    t -= t_offset;
#line 315 "strsna.f"
    vl_dim1 = *ldvl;
#line 315 "strsna.f"
    vl_offset = 1 + vl_dim1;
#line 315 "strsna.f"
    vl -= vl_offset;
#line 315 "strsna.f"
    vr_dim1 = *ldvr;
#line 315 "strsna.f"
    vr_offset = 1 + vr_dim1;
#line 315 "strsna.f"
    vr -= vr_offset;
#line 315 "strsna.f"
    --s;
#line 315 "strsna.f"
    --sep;
#line 315 "strsna.f"
    work_dim1 = *ldwork;
#line 315 "strsna.f"
    work_offset = 1 + work_dim1;
#line 315 "strsna.f"
    work -= work_offset;
#line 315 "strsna.f"
    --iwork;
#line 315 "strsna.f"

#line 315 "strsna.f"
    /* Function Body */
#line 315 "strsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 316 "strsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 317 "strsna.f"
    wantsp = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 319 "strsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 321 "strsna.f"
    *info = 0;
#line 322 "strsna.f"
    if (! wants && ! wantsp) {
#line 323 "strsna.f"
	*info = -1;
#line 324 "strsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 325 "strsna.f"
	*info = -2;
#line 326 "strsna.f"
    } else if (*n < 0) {
#line 327 "strsna.f"
	*info = -4;
#line 328 "strsna.f"
    } else if (*ldt < max(1,*n)) {
#line 329 "strsna.f"
	*info = -6;
#line 330 "strsna.f"
    } else if (*ldvl < 1 || wants && *ldvl < *n) {
#line 331 "strsna.f"
	*info = -8;
#line 332 "strsna.f"
    } else if (*ldvr < 1 || wants && *ldvr < *n) {
#line 333 "strsna.f"
	*info = -10;
#line 334 "strsna.f"
    } else {

/*        Set M to the number of eigenpairs for which condition numbers */
/*        are required, and test MM. */

#line 339 "strsna.f"
	if (somcon) {
#line 340 "strsna.f"
	    *m = 0;
#line 341 "strsna.f"
	    pair = FALSE_;
#line 342 "strsna.f"
	    i__1 = *n;
#line 342 "strsna.f"
	    for (k = 1; k <= i__1; ++k) {
#line 343 "strsna.f"
		if (pair) {
#line 344 "strsna.f"
		    pair = FALSE_;
#line 345 "strsna.f"
		} else {
#line 346 "strsna.f"
		    if (k < *n) {
#line 347 "strsna.f"
			if (t[k + 1 + k * t_dim1] == 0.) {
#line 348 "strsna.f"
			    if (select[k]) {
#line 348 "strsna.f"
				++(*m);
#line 348 "strsna.f"
			    }
#line 350 "strsna.f"
			} else {
#line 351 "strsna.f"
			    pair = TRUE_;
#line 352 "strsna.f"
			    if (select[k] || select[k + 1]) {
#line 352 "strsna.f"
				*m += 2;
#line 352 "strsna.f"
			    }
#line 354 "strsna.f"
			}
#line 355 "strsna.f"
		    } else {
#line 356 "strsna.f"
			if (select[*n]) {
#line 356 "strsna.f"
			    ++(*m);
#line 356 "strsna.f"
			}
#line 358 "strsna.f"
		    }
#line 359 "strsna.f"
		}
#line 360 "strsna.f"
/* L10: */
#line 360 "strsna.f"
	    }
#line 361 "strsna.f"
	} else {
#line 362 "strsna.f"
	    *m = *n;
#line 363 "strsna.f"
	}

#line 365 "strsna.f"
	if (*mm < *m) {
#line 366 "strsna.f"
	    *info = -13;
#line 367 "strsna.f"
	} else if (*ldwork < 1 || wantsp && *ldwork < *n) {
#line 368 "strsna.f"
	    *info = -16;
#line 369 "strsna.f"
	}
#line 370 "strsna.f"
    }
#line 371 "strsna.f"
    if (*info != 0) {
#line 372 "strsna.f"
	i__1 = -(*info);
#line 372 "strsna.f"
	xerbla_("STRSNA", &i__1, (ftnlen)6);
#line 373 "strsna.f"
	return 0;
#line 374 "strsna.f"
    }

/*     Quick return if possible */

#line 378 "strsna.f"
    if (*n == 0) {
#line 378 "strsna.f"
	return 0;
#line 378 "strsna.f"
    }

#line 381 "strsna.f"
    if (*n == 1) {
#line 382 "strsna.f"
	if (somcon) {
#line 383 "strsna.f"
	    if (! select[1]) {
#line 383 "strsna.f"
		return 0;
#line 383 "strsna.f"
	    }
#line 385 "strsna.f"
	}
#line 386 "strsna.f"
	if (wants) {
#line 386 "strsna.f"
	    s[1] = 1.;
#line 386 "strsna.f"
	}
#line 388 "strsna.f"
	if (wantsp) {
#line 388 "strsna.f"
	    sep[1] = (d__1 = t[t_dim1 + 1], abs(d__1));
#line 388 "strsna.f"
	}
#line 390 "strsna.f"
	return 0;
#line 391 "strsna.f"
    }

/*     Get machine constants */

#line 395 "strsna.f"
    eps = slamch_("P", (ftnlen)1);
#line 396 "strsna.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 397 "strsna.f"
    bignum = 1. / smlnum;
#line 398 "strsna.f"
    slabad_(&smlnum, &bignum);

#line 400 "strsna.f"
    ks = 0;
#line 401 "strsna.f"
    pair = FALSE_;
#line 402 "strsna.f"
    i__1 = *n;
#line 402 "strsna.f"
    for (k = 1; k <= i__1; ++k) {

/*        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block. */

#line 406 "strsna.f"
	if (pair) {
#line 407 "strsna.f"
	    pair = FALSE_;
#line 408 "strsna.f"
	    goto L60;
#line 409 "strsna.f"
	} else {
#line 410 "strsna.f"
	    if (k < *n) {
#line 410 "strsna.f"
		pair = t[k + 1 + k * t_dim1] != 0.;
#line 410 "strsna.f"
	    }
#line 412 "strsna.f"
	}

/*        Determine whether condition numbers are required for the k-th */
/*        eigenpair. */

#line 417 "strsna.f"
	if (somcon) {
#line 418 "strsna.f"
	    if (pair) {
#line 419 "strsna.f"
		if (! select[k] && ! select[k + 1]) {
#line 419 "strsna.f"
		    goto L60;
#line 419 "strsna.f"
		}
#line 421 "strsna.f"
	    } else {
#line 422 "strsna.f"
		if (! select[k]) {
#line 422 "strsna.f"
		    goto L60;
#line 422 "strsna.f"
		}
#line 424 "strsna.f"
	    }
#line 425 "strsna.f"
	}

#line 427 "strsna.f"
	++ks;

#line 429 "strsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 434 "strsna.f"
	    if (! pair) {

/*              Real eigenvalue. */

#line 438 "strsna.f"
		prod = sdot_(n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * 
			vl_dim1 + 1], &c__1);
#line 439 "strsna.f"
		rnrm = snrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 440 "strsna.f"
		lnrm = snrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 441 "strsna.f"
		s[ks] = abs(prod) / (rnrm * lnrm);
#line 442 "strsna.f"
	    } else {

/*              Complex eigenvalue. */

#line 446 "strsna.f"
		prod1 = sdot_(n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * 
			vl_dim1 + 1], &c__1);
#line 447 "strsna.f"
		prod1 += sdot_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1, &vl[(ks 
			+ 1) * vl_dim1 + 1], &c__1);
#line 449 "strsna.f"
		prod2 = sdot_(n, &vl[ks * vl_dim1 + 1], &c__1, &vr[(ks + 1) * 
			vr_dim1 + 1], &c__1);
#line 450 "strsna.f"
		prod2 -= sdot_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1, &vr[ks *
			 vr_dim1 + 1], &c__1);
#line 452 "strsna.f"
		d__1 = snrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 452 "strsna.f"
		d__2 = snrm2_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
#line 452 "strsna.f"
		rnrm = slapy2_(&d__1, &d__2);
#line 454 "strsna.f"
		d__1 = snrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 454 "strsna.f"
		d__2 = snrm2_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
#line 454 "strsna.f"
		lnrm = slapy2_(&d__1, &d__2);
#line 456 "strsna.f"
		cond = slapy2_(&prod1, &prod2) / (rnrm * lnrm);
#line 457 "strsna.f"
		s[ks] = cond;
#line 458 "strsna.f"
		s[ks + 1] = cond;
#line 459 "strsna.f"
	    }
#line 460 "strsna.f"
	}

#line 462 "strsna.f"
	if (wantsp) {

/*           Estimate the reciprocal condition number of the k-th */
/*           eigenvector. */

/*           Copy the matrix T to the array WORK and swap the diagonal */
/*           block beginning at T(k,k) to the (1,1) position. */

#line 470 "strsna.f"
	    slacpy_("Full", n, n, &t[t_offset], ldt, &work[work_offset], 
		    ldwork, (ftnlen)4);
#line 471 "strsna.f"
	    ifst = k;
#line 472 "strsna.f"
	    ilst = 1;
#line 473 "strsna.f"
	    strexc_("No Q", n, &work[work_offset], ldwork, dummy, &c__1, &
		    ifst, &ilst, &work[(*n + 1) * work_dim1 + 1], &ierr, (
		    ftnlen)4);

#line 476 "strsna.f"
	    if (ierr == 1 || ierr == 2) {

/*              Could not swap because blocks not well separated */

#line 480 "strsna.f"
		scale = 1.;
#line 481 "strsna.f"
		est = bignum;
#line 482 "strsna.f"
	    } else {

/*              Reordering successful */

#line 486 "strsna.f"
		if (work[work_dim1 + 2] == 0.) {

/*                 Form C = T22 - lambda*I in WORK(2:N,2:N). */

#line 490 "strsna.f"
		    i__2 = *n;
#line 490 "strsna.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 491 "strsna.f"
			work[i__ + i__ * work_dim1] -= work[work_dim1 + 1];
#line 492 "strsna.f"
/* L20: */
#line 492 "strsna.f"
		    }
#line 493 "strsna.f"
		    n2 = 1;
#line 494 "strsna.f"
		    nn = *n - 1;
#line 495 "strsna.f"
		} else {

/*                 Triangularize the 2 by 2 block by unitary */
/*                 transformation U = [  cs   i*ss ] */
/*                                    [ i*ss   cs  ]. */
/*                 such that the (1,1) position of WORK is complex */
/*                 eigenvalue lambda with positive imaginary part. (2,2) */
/*                 position of WORK is the complex eigenvalue lambda */
/*                 with negative imaginary  part. */

#line 505 "strsna.f"
		    mu = sqrt((d__1 = work[(work_dim1 << 1) + 1], abs(d__1))) 
			    * sqrt((d__2 = work[work_dim1 + 2], abs(d__2)));
#line 507 "strsna.f"
		    delta = slapy2_(&mu, &work[work_dim1 + 2]);
#line 508 "strsna.f"
		    cs = mu / delta;
#line 509 "strsna.f"
		    sn = -work[work_dim1 + 2] / delta;

/*                 Form */

/*                 C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ] */
/*                                          [   mu                     ] */
/*                                          [         ..               ] */
/*                                          [             ..           ] */
/*                                          [                  mu      ] */
/*                 where C**T is transpose of matrix C, */
/*                 and RWORK is stored starting in the N+1-st column of */
/*                 WORK. */

#line 522 "strsna.f"
		    i__2 = *n;
#line 522 "strsna.f"
		    for (j = 3; j <= i__2; ++j) {
#line 523 "strsna.f"
			work[j * work_dim1 + 2] = cs * work[j * work_dim1 + 2]
				;
#line 524 "strsna.f"
			work[j + j * work_dim1] -= work[work_dim1 + 1];
#line 525 "strsna.f"
/* L30: */
#line 525 "strsna.f"
		    }
#line 526 "strsna.f"
		    work[(work_dim1 << 1) + 2] = 0.;

#line 528 "strsna.f"
		    work[(*n + 1) * work_dim1 + 1] = mu * 2.;
#line 529 "strsna.f"
		    i__2 = *n - 1;
#line 529 "strsna.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 530 "strsna.f"
			work[i__ + (*n + 1) * work_dim1] = sn * work[(i__ + 1)
				 * work_dim1 + 1];
#line 531 "strsna.f"
/* L40: */
#line 531 "strsna.f"
		    }
#line 532 "strsna.f"
		    n2 = 2;
#line 533 "strsna.f"
		    nn = *n - 1 << 1;
#line 534 "strsna.f"
		}

/*              Estimate norm(inv(C**T)) */

#line 538 "strsna.f"
		est = 0.;
#line 539 "strsna.f"
		kase = 0;
#line 540 "strsna.f"
L50:
#line 541 "strsna.f"
		slacn2_(&nn, &work[(*n + 2) * work_dim1 + 1], &work[(*n + 4) *
			 work_dim1 + 1], &iwork[1], &est, &kase, isave);
#line 543 "strsna.f"
		if (kase != 0) {
#line 544 "strsna.f"
		    if (kase == 1) {
#line 545 "strsna.f"
			if (n2 == 1) {

/*                       Real eigenvalue: solve C**T*x = scale*c. */

#line 549 "strsna.f"
			    i__2 = *n - 1;
#line 549 "strsna.f"
			    slaqtr_(&c_true, &c_true, &i__2, &work[(work_dim1 
				    << 1) + 2], ldwork, dummy, &dumm, &scale, 
				    &work[(*n + 4) * work_dim1 + 1], &work[(*
				    n + 6) * work_dim1 + 1], &ierr);
#line 553 "strsna.f"
			} else {

/*                       Complex eigenvalue: solve */
/*                       C**T*(p+iq) = scale*(c+id) in real arithmetic. */

#line 558 "strsna.f"
			    i__2 = *n - 1;
#line 558 "strsna.f"
			    slaqtr_(&c_true, &c_false, &i__2, &work[(
				    work_dim1 << 1) + 2], ldwork, &work[(*n + 
				    1) * work_dim1 + 1], &mu, &scale, &work[(*
				    n + 4) * work_dim1 + 1], &work[(*n + 6) * 
				    work_dim1 + 1], &ierr);
#line 562 "strsna.f"
			}
#line 563 "strsna.f"
		    } else {
#line 564 "strsna.f"
			if (n2 == 1) {

/*                       Real eigenvalue: solve C*x = scale*c. */

#line 568 "strsna.f"
			    i__2 = *n - 1;
#line 568 "strsna.f"
			    slaqtr_(&c_false, &c_true, &i__2, &work[(
				    work_dim1 << 1) + 2], ldwork, dummy, &
				    dumm, &scale, &work[(*n + 4) * work_dim1 
				    + 1], &work[(*n + 6) * work_dim1 + 1], &
				    ierr);
#line 572 "strsna.f"
			} else {

/*                       Complex eigenvalue: solve */
/*                       C*(p+iq) = scale*(c+id) in real arithmetic. */

#line 577 "strsna.f"
			    i__2 = *n - 1;
#line 577 "strsna.f"
			    slaqtr_(&c_false, &c_false, &i__2, &work[(
				    work_dim1 << 1) + 2], ldwork, &work[(*n + 
				    1) * work_dim1 + 1], &mu, &scale, &work[(*
				    n + 4) * work_dim1 + 1], &work[(*n + 6) * 
				    work_dim1 + 1], &ierr);

#line 583 "strsna.f"
			}
#line 584 "strsna.f"
		    }

#line 586 "strsna.f"
		    goto L50;
#line 587 "strsna.f"
		}
#line 588 "strsna.f"
	    }

#line 590 "strsna.f"
	    sep[ks] = scale / max(est,smlnum);
#line 591 "strsna.f"
	    if (pair) {
#line 591 "strsna.f"
		sep[ks + 1] = sep[ks];
#line 591 "strsna.f"
	    }
#line 593 "strsna.f"
	}

#line 595 "strsna.f"
	if (pair) {
#line 595 "strsna.f"
	    ++ks;
#line 595 "strsna.f"
	}

#line 598 "strsna.f"
L60:
#line 598 "strsna.f"
	;
#line 598 "strsna.f"
    }
#line 599 "strsna.f"
    return 0;

/*     End of STRSNA */

} /* strsna_ */

