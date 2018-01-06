#line 1 "dtrsna.f"
/* dtrsna.f -- translated by f2c (version 20100827).
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

#line 1 "dtrsna.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_true = TRUE_;
static logical c_false = FALSE_;

/* > \brief \b DTRSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                          LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or right eigenvectors of a real upper */
/* > quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q */
/* > orthogonal). */
/* > */
/* > T must be in Schur canonical form (as returned by DHSEQR), that is, */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
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
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,M) */
/* >          If JOB = 'E' or 'B', VL must contain left eigenvectors of T */
/* >          (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
/* >          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* >          must be stored in consecutive columns of VL, as returned by */
/* >          DHSEIN or DTREVC. */
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
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,M) */
/* >          If JOB = 'E' or 'B', VR must contain right eigenvectors of T */
/* >          (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
/* >          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
/* >          must be stored in consecutive columns of VR, as returned by */
/* >          DHSEIN or DTREVC. */
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
/* >          SEP is DOUBLE PRECISION array, dimension (MM) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (LDWORK,N+6) */
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

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dtrsna_(char *job, char *howmny, logical *select, 
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
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical pair;
    static integer ierr;
    static doublereal dumm, prod;
    static integer ifst;
    static doublereal lnrm;
    static integer ilst;
    static doublereal rnrm;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal prod1, prod2, scale, delta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical wants;
    static doublereal dummy[1];
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical wantbh;
    extern /* Subroutine */ int dlaqtr_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *), dtrexc_(char *, integer *
	    , doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static logical somcon;
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

#line 315 "dtrsna.f"
    /* Parameter adjustments */
#line 315 "dtrsna.f"
    --select;
#line 315 "dtrsna.f"
    t_dim1 = *ldt;
#line 315 "dtrsna.f"
    t_offset = 1 + t_dim1;
#line 315 "dtrsna.f"
    t -= t_offset;
#line 315 "dtrsna.f"
    vl_dim1 = *ldvl;
#line 315 "dtrsna.f"
    vl_offset = 1 + vl_dim1;
#line 315 "dtrsna.f"
    vl -= vl_offset;
#line 315 "dtrsna.f"
    vr_dim1 = *ldvr;
#line 315 "dtrsna.f"
    vr_offset = 1 + vr_dim1;
#line 315 "dtrsna.f"
    vr -= vr_offset;
#line 315 "dtrsna.f"
    --s;
#line 315 "dtrsna.f"
    --sep;
#line 315 "dtrsna.f"
    work_dim1 = *ldwork;
#line 315 "dtrsna.f"
    work_offset = 1 + work_dim1;
#line 315 "dtrsna.f"
    work -= work_offset;
#line 315 "dtrsna.f"
    --iwork;
#line 315 "dtrsna.f"

#line 315 "dtrsna.f"
    /* Function Body */
#line 315 "dtrsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 316 "dtrsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 317 "dtrsna.f"
    wantsp = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 319 "dtrsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 321 "dtrsna.f"
    *info = 0;
#line 322 "dtrsna.f"
    if (! wants && ! wantsp) {
#line 323 "dtrsna.f"
	*info = -1;
#line 324 "dtrsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 325 "dtrsna.f"
	*info = -2;
#line 326 "dtrsna.f"
    } else if (*n < 0) {
#line 327 "dtrsna.f"
	*info = -4;
#line 328 "dtrsna.f"
    } else if (*ldt < max(1,*n)) {
#line 329 "dtrsna.f"
	*info = -6;
#line 330 "dtrsna.f"
    } else if (*ldvl < 1 || wants && *ldvl < *n) {
#line 331 "dtrsna.f"
	*info = -8;
#line 332 "dtrsna.f"
    } else if (*ldvr < 1 || wants && *ldvr < *n) {
#line 333 "dtrsna.f"
	*info = -10;
#line 334 "dtrsna.f"
    } else {

/*        Set M to the number of eigenpairs for which condition numbers */
/*        are required, and test MM. */

#line 339 "dtrsna.f"
	if (somcon) {
#line 340 "dtrsna.f"
	    *m = 0;
#line 341 "dtrsna.f"
	    pair = FALSE_;
#line 342 "dtrsna.f"
	    i__1 = *n;
#line 342 "dtrsna.f"
	    for (k = 1; k <= i__1; ++k) {
#line 343 "dtrsna.f"
		if (pair) {
#line 344 "dtrsna.f"
		    pair = FALSE_;
#line 345 "dtrsna.f"
		} else {
#line 346 "dtrsna.f"
		    if (k < *n) {
#line 347 "dtrsna.f"
			if (t[k + 1 + k * t_dim1] == 0.) {
#line 348 "dtrsna.f"
			    if (select[k]) {
#line 348 "dtrsna.f"
				++(*m);
#line 348 "dtrsna.f"
			    }
#line 350 "dtrsna.f"
			} else {
#line 351 "dtrsna.f"
			    pair = TRUE_;
#line 352 "dtrsna.f"
			    if (select[k] || select[k + 1]) {
#line 352 "dtrsna.f"
				*m += 2;
#line 352 "dtrsna.f"
			    }
#line 354 "dtrsna.f"
			}
#line 355 "dtrsna.f"
		    } else {
#line 356 "dtrsna.f"
			if (select[*n]) {
#line 356 "dtrsna.f"
			    ++(*m);
#line 356 "dtrsna.f"
			}
#line 358 "dtrsna.f"
		    }
#line 359 "dtrsna.f"
		}
#line 360 "dtrsna.f"
/* L10: */
#line 360 "dtrsna.f"
	    }
#line 361 "dtrsna.f"
	} else {
#line 362 "dtrsna.f"
	    *m = *n;
#line 363 "dtrsna.f"
	}

#line 365 "dtrsna.f"
	if (*mm < *m) {
#line 366 "dtrsna.f"
	    *info = -13;
#line 367 "dtrsna.f"
	} else if (*ldwork < 1 || wantsp && *ldwork < *n) {
#line 368 "dtrsna.f"
	    *info = -16;
#line 369 "dtrsna.f"
	}
#line 370 "dtrsna.f"
    }
#line 371 "dtrsna.f"
    if (*info != 0) {
#line 372 "dtrsna.f"
	i__1 = -(*info);
#line 372 "dtrsna.f"
	xerbla_("DTRSNA", &i__1, (ftnlen)6);
#line 373 "dtrsna.f"
	return 0;
#line 374 "dtrsna.f"
    }

/*     Quick return if possible */

#line 378 "dtrsna.f"
    if (*n == 0) {
#line 378 "dtrsna.f"
	return 0;
#line 378 "dtrsna.f"
    }

#line 381 "dtrsna.f"
    if (*n == 1) {
#line 382 "dtrsna.f"
	if (somcon) {
#line 383 "dtrsna.f"
	    if (! select[1]) {
#line 383 "dtrsna.f"
		return 0;
#line 383 "dtrsna.f"
	    }
#line 385 "dtrsna.f"
	}
#line 386 "dtrsna.f"
	if (wants) {
#line 386 "dtrsna.f"
	    s[1] = 1.;
#line 386 "dtrsna.f"
	}
#line 388 "dtrsna.f"
	if (wantsp) {
#line 388 "dtrsna.f"
	    sep[1] = (d__1 = t[t_dim1 + 1], abs(d__1));
#line 388 "dtrsna.f"
	}
#line 390 "dtrsna.f"
	return 0;
#line 391 "dtrsna.f"
    }

/*     Get machine constants */

#line 395 "dtrsna.f"
    eps = dlamch_("P", (ftnlen)1);
#line 396 "dtrsna.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 397 "dtrsna.f"
    bignum = 1. / smlnum;
#line 398 "dtrsna.f"
    dlabad_(&smlnum, &bignum);

#line 400 "dtrsna.f"
    ks = 0;
#line 401 "dtrsna.f"
    pair = FALSE_;
#line 402 "dtrsna.f"
    i__1 = *n;
#line 402 "dtrsna.f"
    for (k = 1; k <= i__1; ++k) {

/*        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block. */

#line 406 "dtrsna.f"
	if (pair) {
#line 407 "dtrsna.f"
	    pair = FALSE_;
#line 408 "dtrsna.f"
	    goto L60;
#line 409 "dtrsna.f"
	} else {
#line 410 "dtrsna.f"
	    if (k < *n) {
#line 410 "dtrsna.f"
		pair = t[k + 1 + k * t_dim1] != 0.;
#line 410 "dtrsna.f"
	    }
#line 412 "dtrsna.f"
	}

/*        Determine whether condition numbers are required for the k-th */
/*        eigenpair. */

#line 417 "dtrsna.f"
	if (somcon) {
#line 418 "dtrsna.f"
	    if (pair) {
#line 419 "dtrsna.f"
		if (! select[k] && ! select[k + 1]) {
#line 419 "dtrsna.f"
		    goto L60;
#line 419 "dtrsna.f"
		}
#line 421 "dtrsna.f"
	    } else {
#line 422 "dtrsna.f"
		if (! select[k]) {
#line 422 "dtrsna.f"
		    goto L60;
#line 422 "dtrsna.f"
		}
#line 424 "dtrsna.f"
	    }
#line 425 "dtrsna.f"
	}

#line 427 "dtrsna.f"
	++ks;

#line 429 "dtrsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 434 "dtrsna.f"
	    if (! pair) {

/*              Real eigenvalue. */

#line 438 "dtrsna.f"
		prod = ddot_(n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * 
			vl_dim1 + 1], &c__1);
#line 439 "dtrsna.f"
		rnrm = dnrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 440 "dtrsna.f"
		lnrm = dnrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 441 "dtrsna.f"
		s[ks] = abs(prod) / (rnrm * lnrm);
#line 442 "dtrsna.f"
	    } else {

/*              Complex eigenvalue. */

#line 446 "dtrsna.f"
		prod1 = ddot_(n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * 
			vl_dim1 + 1], &c__1);
#line 447 "dtrsna.f"
		prod1 += ddot_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1, &vl[(ks 
			+ 1) * vl_dim1 + 1], &c__1);
#line 449 "dtrsna.f"
		prod2 = ddot_(n, &vl[ks * vl_dim1 + 1], &c__1, &vr[(ks + 1) * 
			vr_dim1 + 1], &c__1);
#line 450 "dtrsna.f"
		prod2 -= ddot_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1, &vr[ks *
			 vr_dim1 + 1], &c__1);
#line 452 "dtrsna.f"
		d__1 = dnrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 452 "dtrsna.f"
		d__2 = dnrm2_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
#line 452 "dtrsna.f"
		rnrm = dlapy2_(&d__1, &d__2);
#line 454 "dtrsna.f"
		d__1 = dnrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 454 "dtrsna.f"
		d__2 = dnrm2_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
#line 454 "dtrsna.f"
		lnrm = dlapy2_(&d__1, &d__2);
#line 456 "dtrsna.f"
		cond = dlapy2_(&prod1, &prod2) / (rnrm * lnrm);
#line 457 "dtrsna.f"
		s[ks] = cond;
#line 458 "dtrsna.f"
		s[ks + 1] = cond;
#line 459 "dtrsna.f"
	    }
#line 460 "dtrsna.f"
	}

#line 462 "dtrsna.f"
	if (wantsp) {

/*           Estimate the reciprocal condition number of the k-th */
/*           eigenvector. */

/*           Copy the matrix T to the array WORK and swap the diagonal */
/*           block beginning at T(k,k) to the (1,1) position. */

#line 470 "dtrsna.f"
	    dlacpy_("Full", n, n, &t[t_offset], ldt, &work[work_offset], 
		    ldwork, (ftnlen)4);
#line 471 "dtrsna.f"
	    ifst = k;
#line 472 "dtrsna.f"
	    ilst = 1;
#line 473 "dtrsna.f"
	    dtrexc_("No Q", n, &work[work_offset], ldwork, dummy, &c__1, &
		    ifst, &ilst, &work[(*n + 1) * work_dim1 + 1], &ierr, (
		    ftnlen)4);

#line 476 "dtrsna.f"
	    if (ierr == 1 || ierr == 2) {

/*              Could not swap because blocks not well separated */

#line 480 "dtrsna.f"
		scale = 1.;
#line 481 "dtrsna.f"
		est = bignum;
#line 482 "dtrsna.f"
	    } else {

/*              Reordering successful */

#line 486 "dtrsna.f"
		if (work[work_dim1 + 2] == 0.) {

/*                 Form C = T22 - lambda*I in WORK(2:N,2:N). */

#line 490 "dtrsna.f"
		    i__2 = *n;
#line 490 "dtrsna.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 491 "dtrsna.f"
			work[i__ + i__ * work_dim1] -= work[work_dim1 + 1];
#line 492 "dtrsna.f"
/* L20: */
#line 492 "dtrsna.f"
		    }
#line 493 "dtrsna.f"
		    n2 = 1;
#line 494 "dtrsna.f"
		    nn = *n - 1;
#line 495 "dtrsna.f"
		} else {

/*                 Triangularize the 2 by 2 block by unitary */
/*                 transformation U = [  cs   i*ss ] */
/*                                    [ i*ss   cs  ]. */
/*                 such that the (1,1) position of WORK is complex */
/*                 eigenvalue lambda with positive imaginary part. (2,2) */
/*                 position of WORK is the complex eigenvalue lambda */
/*                 with negative imaginary  part. */

#line 505 "dtrsna.f"
		    mu = sqrt((d__1 = work[(work_dim1 << 1) + 1], abs(d__1))) 
			    * sqrt((d__2 = work[work_dim1 + 2], abs(d__2)));
#line 507 "dtrsna.f"
		    delta = dlapy2_(&mu, &work[work_dim1 + 2]);
#line 508 "dtrsna.f"
		    cs = mu / delta;
#line 509 "dtrsna.f"
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

#line 522 "dtrsna.f"
		    i__2 = *n;
#line 522 "dtrsna.f"
		    for (j = 3; j <= i__2; ++j) {
#line 523 "dtrsna.f"
			work[j * work_dim1 + 2] = cs * work[j * work_dim1 + 2]
				;
#line 524 "dtrsna.f"
			work[j + j * work_dim1] -= work[work_dim1 + 1];
#line 525 "dtrsna.f"
/* L30: */
#line 525 "dtrsna.f"
		    }
#line 526 "dtrsna.f"
		    work[(work_dim1 << 1) + 2] = 0.;

#line 528 "dtrsna.f"
		    work[(*n + 1) * work_dim1 + 1] = mu * 2.;
#line 529 "dtrsna.f"
		    i__2 = *n - 1;
#line 529 "dtrsna.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 530 "dtrsna.f"
			work[i__ + (*n + 1) * work_dim1] = sn * work[(i__ + 1)
				 * work_dim1 + 1];
#line 531 "dtrsna.f"
/* L40: */
#line 531 "dtrsna.f"
		    }
#line 532 "dtrsna.f"
		    n2 = 2;
#line 533 "dtrsna.f"
		    nn = *n - 1 << 1;
#line 534 "dtrsna.f"
		}

/*              Estimate norm(inv(C**T)) */

#line 538 "dtrsna.f"
		est = 0.;
#line 539 "dtrsna.f"
		kase = 0;
#line 540 "dtrsna.f"
L50:
#line 541 "dtrsna.f"
		dlacn2_(&nn, &work[(*n + 2) * work_dim1 + 1], &work[(*n + 4) *
			 work_dim1 + 1], &iwork[1], &est, &kase, isave);
#line 543 "dtrsna.f"
		if (kase != 0) {
#line 544 "dtrsna.f"
		    if (kase == 1) {
#line 545 "dtrsna.f"
			if (n2 == 1) {

/*                       Real eigenvalue: solve C**T*x = scale*c. */

#line 549 "dtrsna.f"
			    i__2 = *n - 1;
#line 549 "dtrsna.f"
			    dlaqtr_(&c_true, &c_true, &i__2, &work[(work_dim1 
				    << 1) + 2], ldwork, dummy, &dumm, &scale, 
				    &work[(*n + 4) * work_dim1 + 1], &work[(*
				    n + 6) * work_dim1 + 1], &ierr);
#line 553 "dtrsna.f"
			} else {

/*                       Complex eigenvalue: solve */
/*                       C**T*(p+iq) = scale*(c+id) in real arithmetic. */

#line 558 "dtrsna.f"
			    i__2 = *n - 1;
#line 558 "dtrsna.f"
			    dlaqtr_(&c_true, &c_false, &i__2, &work[(
				    work_dim1 << 1) + 2], ldwork, &work[(*n + 
				    1) * work_dim1 + 1], &mu, &scale, &work[(*
				    n + 4) * work_dim1 + 1], &work[(*n + 6) * 
				    work_dim1 + 1], &ierr);
#line 562 "dtrsna.f"
			}
#line 563 "dtrsna.f"
		    } else {
#line 564 "dtrsna.f"
			if (n2 == 1) {

/*                       Real eigenvalue: solve C*x = scale*c. */

#line 568 "dtrsna.f"
			    i__2 = *n - 1;
#line 568 "dtrsna.f"
			    dlaqtr_(&c_false, &c_true, &i__2, &work[(
				    work_dim1 << 1) + 2], ldwork, dummy, &
				    dumm, &scale, &work[(*n + 4) * work_dim1 
				    + 1], &work[(*n + 6) * work_dim1 + 1], &
				    ierr);
#line 572 "dtrsna.f"
			} else {

/*                       Complex eigenvalue: solve */
/*                       C*(p+iq) = scale*(c+id) in real arithmetic. */

#line 577 "dtrsna.f"
			    i__2 = *n - 1;
#line 577 "dtrsna.f"
			    dlaqtr_(&c_false, &c_false, &i__2, &work[(
				    work_dim1 << 1) + 2], ldwork, &work[(*n + 
				    1) * work_dim1 + 1], &mu, &scale, &work[(*
				    n + 4) * work_dim1 + 1], &work[(*n + 6) * 
				    work_dim1 + 1], &ierr);

#line 583 "dtrsna.f"
			}
#line 584 "dtrsna.f"
		    }

#line 586 "dtrsna.f"
		    goto L50;
#line 587 "dtrsna.f"
		}
#line 588 "dtrsna.f"
	    }

#line 590 "dtrsna.f"
	    sep[ks] = scale / max(est,smlnum);
#line 591 "dtrsna.f"
	    if (pair) {
#line 591 "dtrsna.f"
		sep[ks + 1] = sep[ks];
#line 591 "dtrsna.f"
	    }
#line 593 "dtrsna.f"
	}

#line 595 "dtrsna.f"
	if (pair) {
#line 595 "dtrsna.f"
	    ++ks;
#line 595 "dtrsna.f"
	}

#line 598 "dtrsna.f"
L60:
#line 598 "dtrsna.f"
	;
#line 598 "dtrsna.f"
    }
#line 599 "dtrsna.f"
    return 0;

/*     End of DTRSNA */

} /* dtrsna_ */

