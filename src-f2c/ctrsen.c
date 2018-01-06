#line 1 "ctrsen.f"
/* ctrsen.f -- translated by f2c (version 20100827).
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

#line 1 "ctrsen.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief \b CTRSEN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTRSEN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrsen.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrsen.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrsen.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S, */
/*                          SEP, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, JOB */
/*       INTEGER            INFO, LDQ, LDT, LWORK, M, N */
/*       REAL               S, SEP */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       COMPLEX            Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRSEN reorders the Schur factorization of a complex matrix */
/* > A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in */
/* > the leading positions on the diagonal of the upper triangular matrix */
/* > T, and the leading columns of Q form an orthonormal basis of the */
/* > corresponding right invariant subspace. */
/* > */
/* > Optionally the routine computes the reciprocal condition numbers of */
/* > the cluster of eigenvalues and/or the invariant subspace. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies whether condition numbers are required for the */
/* >          cluster of eigenvalues (S) or the invariant subspace (SEP): */
/* >          = 'N': none; */
/* >          = 'E': for eigenvalues only (S); */
/* >          = 'V': for invariant subspace only (SEP); */
/* >          = 'B': for both eigenvalues and invariant subspace (S and */
/* >                 SEP). */
/* > \endverbatim */
/* > */
/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'V': update the matrix Q of Schur vectors; */
/* >          = 'N': do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          SELECT specifies the eigenvalues in the selected cluster. To */
/* >          select the j-th eigenvalue, SELECT(j) must be set to .TRUE.. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,N) */
/* >          On entry, the upper triangular matrix T. */
/* >          On exit, T is overwritten by the reordered matrix T, with the */
/* >          selected eigenvalues as the leading diagonal elements. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX array, dimension (LDQ,N) */
/* >          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
/* >          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
/* >          unitary transformation matrix which reorders T; the leading M */
/* >          columns of Q form an orthonormal basis for the specified */
/* >          invariant subspace. */
/* >          If COMPQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. */
/* >          LDQ >= 1; and if COMPQ = 'V', LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (N) */
/* >          The reordered eigenvalues of T, in the same order as they */
/* >          appear on the diagonal of T. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The dimension of the specified invariant subspace. */
/* >          0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL */
/* >          If JOB = 'E' or 'B', S is a lower bound on the reciprocal */
/* >          condition number for the selected cluster of eigenvalues. */
/* >          S cannot underestimate the true reciprocal condition number */
/* >          by more than a factor of sqrt(N). If M = 0 or N, S = 1. */
/* >          If JOB = 'N' or 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SEP */
/* > \verbatim */
/* >          SEP is REAL */
/* >          If JOB = 'V' or 'B', SEP is the estimated reciprocal */
/* >          condition number of the specified invariant subspace. If */
/* >          M = 0 or N, SEP = norm(T). */
/* >          If JOB = 'N' or 'E', SEP is not referenced. */
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
/* >          If JOB = 'N', LWORK >= 1; */
/* >          if JOB = 'E', LWORK = max(1,M*(N-M)); */
/* >          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  CTRSEN first collects the selected eigenvalues by computing a unitary */
/* >  transformation Z to move them to the top left corner of T. In other */
/* >  words, the selected eigenvalues are the eigenvalues of T11 in: */
/* > */
/* >          Z**H * T * Z = ( T11 T12 ) n1 */
/* >                         (  0  T22 ) n2 */
/* >                            n1  n2 */
/* > */
/* >  where N = n1+n2. The first */
/* >  n1 columns of Z span the specified invariant subspace of T. */
/* > */
/* >  If T has been obtained from the Schur factorization of a matrix */
/* >  A = Q*T*Q**H, then the reordered Schur factorization of A is given by */
/* >  A = (Q*Z)*(Z**H*T*Z)*(Q*Z)**H, and the first n1 columns of Q*Z span the */
/* >  corresponding invariant subspace of A. */
/* > */
/* >  The reciprocal condition number of the average of the eigenvalues of */
/* >  T11 may be returned in S. S lies between 0 (very badly conditioned) */
/* >  and 1 (very well conditioned). It is computed as follows. First we */
/* >  compute R so that */
/* > */
/* >                         P = ( I  R ) n1 */
/* >                             ( 0  0 ) n2 */
/* >                               n1 n2 */
/* > */
/* >  is the projector on the invariant subspace associated with T11. */
/* >  R is the solution of the Sylvester equation: */
/* > */
/* >                        T11*R - R*T22 = T12. */
/* > */
/* >  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote */
/* >  the two-norm of M. Then S is computed as the lower bound */
/* > */
/* >                      (1 + F-norm(R)**2)**(-1/2) */
/* > */
/* >  on the reciprocal of 2-norm(P), the true reciprocal condition number. */
/* >  S cannot underestimate 1 / 2-norm(P) by more than a factor of */
/* >  sqrt(N). */
/* > */
/* >  An approximate error bound for the computed average of the */
/* >  eigenvalues of T11 is */
/* > */
/* >                         EPS * norm(T) / S */
/* > */
/* >  where EPS is the machine precision. */
/* > */
/* >  The reciprocal condition number of the right invariant subspace */
/* >  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP. */
/* >  SEP is defined as the separation of T11 and T22: */
/* > */
/* >                     sep( T11, T22 ) = sigma-min( C ) */
/* > */
/* >  where sigma-min(C) is the smallest singular value of the */
/* >  n1*n2-by-n1*n2 matrix */
/* > */
/* >     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) ) */
/* > */
/* >  I(m) is an m by m identity matrix, and kprod denotes the Kronecker */
/* >  product. We estimate sigma-min(C) by the reciprocal of an estimate of */
/* >  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C) */
/* >  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2). */
/* > */
/* >  When SEP is small, small changes in T can cause large changes in */
/* >  the invariant subspace. An approximate bound on the maximum angular */
/* >  error in the computed right invariant subspace is */
/* > */
/* >                      EPS * norm(T) / SEP */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctrsen_(char *job, char *compq, logical *select, integer 
	*n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, 
	doublecomplex *w, integer *m, doublereal *s, doublereal *sep, 
	doublecomplex *work, integer *lwork, integer *info, ftnlen job_len, 
	ftnlen compq_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, n1, n2, nn, ks;
    static doublereal est;
    static integer kase, ierr;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3], lwmin;
    static logical wantq, wants;
    static doublereal rnorm;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    static doublereal rwork[1];
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical wantbh;
    extern /* Subroutine */ int ctrexc_(char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, integer *, 
	    integer *, ftnlen);
    static logical wantsp;
    extern /* Subroutine */ int ctrsyl_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static logical lquery;


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

/*     Decode and test the input parameters. */

#line 312 "ctrsen.f"
    /* Parameter adjustments */
#line 312 "ctrsen.f"
    --select;
#line 312 "ctrsen.f"
    t_dim1 = *ldt;
#line 312 "ctrsen.f"
    t_offset = 1 + t_dim1;
#line 312 "ctrsen.f"
    t -= t_offset;
#line 312 "ctrsen.f"
    q_dim1 = *ldq;
#line 312 "ctrsen.f"
    q_offset = 1 + q_dim1;
#line 312 "ctrsen.f"
    q -= q_offset;
#line 312 "ctrsen.f"
    --w;
#line 312 "ctrsen.f"
    --work;
#line 312 "ctrsen.f"

#line 312 "ctrsen.f"
    /* Function Body */
#line 312 "ctrsen.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 313 "ctrsen.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 314 "ctrsen.f"
    wantsp = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;
#line 315 "ctrsen.f"
    wantq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of selected eigenvalues. */

#line 319 "ctrsen.f"
    *m = 0;
#line 320 "ctrsen.f"
    i__1 = *n;
#line 320 "ctrsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 321 "ctrsen.f"
	if (select[k]) {
#line 321 "ctrsen.f"
	    ++(*m);
#line 321 "ctrsen.f"
	}
#line 323 "ctrsen.f"
/* L10: */
#line 323 "ctrsen.f"
    }

#line 325 "ctrsen.f"
    n1 = *m;
#line 326 "ctrsen.f"
    n2 = *n - *m;
#line 327 "ctrsen.f"
    nn = n1 * n2;

#line 329 "ctrsen.f"
    *info = 0;
#line 330 "ctrsen.f"
    lquery = *lwork == -1;

#line 332 "ctrsen.f"
    if (wantsp) {
/* Computing MAX */
#line 333 "ctrsen.f"
	i__1 = 1, i__2 = nn << 1;
#line 333 "ctrsen.f"
	lwmin = max(i__1,i__2);
#line 334 "ctrsen.f"
    } else if (lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 335 "ctrsen.f"
	lwmin = 1;
#line 336 "ctrsen.f"
    } else if (lsame_(job, "E", (ftnlen)1, (ftnlen)1)) {
#line 337 "ctrsen.f"
	lwmin = max(1,nn);
#line 338 "ctrsen.f"
    }

#line 340 "ctrsen.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! wants && ! wantsp) {
#line 342 "ctrsen.f"
	*info = -1;
#line 343 "ctrsen.f"
    } else if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 344 "ctrsen.f"
	*info = -2;
#line 345 "ctrsen.f"
    } else if (*n < 0) {
#line 346 "ctrsen.f"
	*info = -4;
#line 347 "ctrsen.f"
    } else if (*ldt < max(1,*n)) {
#line 348 "ctrsen.f"
	*info = -6;
#line 349 "ctrsen.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 350 "ctrsen.f"
	*info = -8;
#line 351 "ctrsen.f"
    } else if (*lwork < lwmin && ! lquery) {
#line 352 "ctrsen.f"
	*info = -14;
#line 353 "ctrsen.f"
    }

#line 355 "ctrsen.f"
    if (*info == 0) {
#line 356 "ctrsen.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 357 "ctrsen.f"
    }

#line 359 "ctrsen.f"
    if (*info != 0) {
#line 360 "ctrsen.f"
	i__1 = -(*info);
#line 360 "ctrsen.f"
	xerbla_("CTRSEN", &i__1, (ftnlen)6);
#line 361 "ctrsen.f"
	return 0;
#line 362 "ctrsen.f"
    } else if (lquery) {
#line 363 "ctrsen.f"
	return 0;
#line 364 "ctrsen.f"
    }

/*     Quick return if possible */

#line 368 "ctrsen.f"
    if (*m == *n || *m == 0) {
#line 369 "ctrsen.f"
	if (wants) {
#line 369 "ctrsen.f"
	    *s = 1.;
#line 369 "ctrsen.f"
	}
#line 371 "ctrsen.f"
	if (wantsp) {
#line 371 "ctrsen.f"
	    *sep = clange_("1", n, n, &t[t_offset], ldt, rwork, (ftnlen)1);
#line 371 "ctrsen.f"
	}
#line 373 "ctrsen.f"
	goto L40;
#line 374 "ctrsen.f"
    }

/*     Collect the selected eigenvalues at the top left corner of T. */

#line 378 "ctrsen.f"
    ks = 0;
#line 379 "ctrsen.f"
    i__1 = *n;
#line 379 "ctrsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 380 "ctrsen.f"
	if (select[k]) {
#line 381 "ctrsen.f"
	    ++ks;

/*           Swap the K-th eigenvalue to position KS. */

#line 385 "ctrsen.f"
	    if (k != ks) {
#line 385 "ctrsen.f"
		ctrexc_(compq, n, &t[t_offset], ldt, &q[q_offset], ldq, &k, &
			ks, &ierr, (ftnlen)1);
#line 385 "ctrsen.f"
	    }
#line 387 "ctrsen.f"
	}
#line 388 "ctrsen.f"
/* L20: */
#line 388 "ctrsen.f"
    }

#line 390 "ctrsen.f"
    if (wants) {

/*        Solve the Sylvester equation for R: */

/*           T11*R - R*T22 = scale*T12 */

#line 396 "ctrsen.f"
	clacpy_("F", &n1, &n2, &t[(n1 + 1) * t_dim1 + 1], ldt, &work[1], &n1, 
		(ftnlen)1);
#line 397 "ctrsen.f"
	ctrsyl_("N", "N", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 1 + (n1 
		+ 1) * t_dim1], ldt, &work[1], &n1, &scale, &ierr, (ftnlen)1, 
		(ftnlen)1);

/*        Estimate the reciprocal of the condition number of the cluster */
/*        of eigenvalues. */

#line 403 "ctrsen.f"
	rnorm = clange_("F", &n1, &n2, &work[1], &n1, rwork, (ftnlen)1);
#line 404 "ctrsen.f"
	if (rnorm == 0.) {
#line 405 "ctrsen.f"
	    *s = 1.;
#line 406 "ctrsen.f"
	} else {
#line 407 "ctrsen.f"
	    *s = scale / (sqrt(scale * scale / rnorm + rnorm) * sqrt(rnorm));
#line 409 "ctrsen.f"
	}
#line 410 "ctrsen.f"
    }

#line 412 "ctrsen.f"
    if (wantsp) {

/*        Estimate sep(T11,T22). */

#line 416 "ctrsen.f"
	est = 0.;
#line 417 "ctrsen.f"
	kase = 0;
#line 418 "ctrsen.f"
L30:
#line 419 "ctrsen.f"
	clacn2_(&nn, &work[nn + 1], &work[1], &est, &kase, isave);
#line 420 "ctrsen.f"
	if (kase != 0) {
#line 421 "ctrsen.f"
	    if (kase == 1) {

/*              Solve T11*R - R*T22 = scale*X. */

#line 425 "ctrsen.f"
		ctrsyl_("N", "N", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 
			1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &
			ierr, (ftnlen)1, (ftnlen)1);
#line 428 "ctrsen.f"
	    } else {

/*              Solve T11**H*R - R*T22**H = scale*X. */

#line 432 "ctrsen.f"
		ctrsyl_("C", "C", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 
			1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &
			ierr, (ftnlen)1, (ftnlen)1);
#line 435 "ctrsen.f"
	    }
#line 436 "ctrsen.f"
	    goto L30;
#line 437 "ctrsen.f"
	}

#line 439 "ctrsen.f"
	*sep = scale / est;
#line 440 "ctrsen.f"
    }

#line 442 "ctrsen.f"
L40:

/*     Copy reordered eigenvalues to W. */

#line 446 "ctrsen.f"
    i__1 = *n;
#line 446 "ctrsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 447 "ctrsen.f"
	i__2 = k;
#line 447 "ctrsen.f"
	i__3 = k + k * t_dim1;
#line 447 "ctrsen.f"
	w[i__2].r = t[i__3].r, w[i__2].i = t[i__3].i;
#line 448 "ctrsen.f"
/* L50: */
#line 448 "ctrsen.f"
    }

#line 450 "ctrsen.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 452 "ctrsen.f"
    return 0;

/*     End of CTRSEN */

} /* ctrsen_ */

