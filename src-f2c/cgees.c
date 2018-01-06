#line 1 "cgees.f"
/* cgees.f -- translated by f2c (version 20100827).
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

#line 1 "cgees.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> CGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgees.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgees.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgees.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, */
/*                         LDVS, WORK, LWORK, RWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVS, SORT */
/*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * ) */
/*       .. */
/*       .. Function Arguments .. */
/*       LOGICAL            SELECT */
/*       EXTERNAL           SELECT */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEES computes for an N-by-N complex nonsymmetric matrix A, the */
/* > eigenvalues, the Schur form T, and, optionally, the matrix of Schur */
/* > vectors Z.  This gives the Schur factorization A = Z*T*(Z**H). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > Schur form so that selected eigenvalues are at the top left. */
/* > The leading columns of Z then form an orthonormal basis for the */
/* > invariant subspace corresponding to the selected eigenvalues. */
/* > */
/* > A complex matrix is in Schur form if it is upper triangular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBVS */
/* > \verbatim */
/* >          JOBVS is CHARACTER*1 */
/* >          = 'N': Schur vectors are not computed; */
/* >          = 'V': Schur vectors are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] SORT */
/* > \verbatim */
/* >          SORT is CHARACTER*1 */
/* >          Specifies whether or not to order the eigenvalues on the */
/* >          diagonal of the Schur form. */
/* >          = 'N': Eigenvalues are not ordered: */
/* >          = 'S': Eigenvalues are ordered (see SELECT). */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is a LOGICAL FUNCTION of one COMPLEX argument */
/* >          SELECT must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'S', SELECT is used to select eigenvalues to order */
/* >          to the top left of the Schur form. */
/* >          IF SORT = 'N', SELECT is not referenced. */
/* >          The eigenvalue W(j) is selected if SELECT(W(j)) is true. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A has been overwritten by its Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SDIM */
/* > \verbatim */
/* >          SDIM is INTEGER */
/* >          If SORT = 'N', SDIM = 0. */
/* >          If SORT = 'S', SDIM = number of eigenvalues for which */
/* >                         SELECT is true. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (N) */
/* >          W contains the computed eigenvalues, in the same order that */
/* >          they appear on the diagonal of the output Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* >          VS is COMPLEX array, dimension (LDVS,N) */
/* >          If JOBVS = 'V', VS contains the unitary matrix Z of Schur */
/* >          vectors. */
/* >          If JOBVS = 'N', VS is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVS */
/* > \verbatim */
/* >          LDVS is INTEGER */
/* >          The leading dimension of the array VS.  LDVS >= 1; if */
/* >          JOBVS = 'V', LDVS >= N. */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BWORK */
/* > \verbatim */
/* >          BWORK is LOGICAL array, dimension (N) */
/* >          Not referenced if SORT = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0: if INFO = i, and i is */
/* >               <= N:  the QR algorithm failed to compute all the */
/* >                      eigenvalues; elements 1:ILO-1 and i+1:N of W */
/* >                      contain those eigenvalues which have converged; */
/* >                      if JOBVS = 'V', VS contains the matrix which */
/* >                      reduces A to its partially converged Schur form. */
/* >               = N+1: the eigenvalues could not be reordered because */
/* >                      some eigenvalues were too close to separate (the */
/* >                      problem is very ill-conditioned); */
/* >               = N+2: after reordering, roundoff changed values of */
/* >                      some complex eigenvalues so that leading */
/* >                      eigenvalues in the Schur form no longer satisfy */
/* >                      SELECT = .TRUE..  This could also be caused by */
/* >                      underflow due to scaling. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGEeigen */

/*  ===================================================================== */
/* Subroutine */ int cgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, 
	doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork,
	 doublereal *rwork, logical *bwork, integer *info, ftnlen jobvs_len, 
	ftnlen sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s;
    static integer ihi, ilo;
    static doublereal dum[1], eps, sep;
    static integer ibal;
    static doublereal anrm;
    static integer ierr, itau, iwrk, icond, ieval;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cgebak_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen), cgebal_(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen), slabad_(doublereal *, doublereal *);
    static logical scalea;
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int cgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int chseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), cunghr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), ctrsen_(char *, char *, logical *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static integer hswork;
    static logical wantst, lquery, wantvs;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Function Arguments .. */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 251 "cgees.f"
    /* Parameter adjustments */
#line 251 "cgees.f"
    a_dim1 = *lda;
#line 251 "cgees.f"
    a_offset = 1 + a_dim1;
#line 251 "cgees.f"
    a -= a_offset;
#line 251 "cgees.f"
    --w;
#line 251 "cgees.f"
    vs_dim1 = *ldvs;
#line 251 "cgees.f"
    vs_offset = 1 + vs_dim1;
#line 251 "cgees.f"
    vs -= vs_offset;
#line 251 "cgees.f"
    --work;
#line 251 "cgees.f"
    --rwork;
#line 251 "cgees.f"
    --bwork;
#line 251 "cgees.f"

#line 251 "cgees.f"
    /* Function Body */
#line 251 "cgees.f"
    *info = 0;
#line 252 "cgees.f"
    lquery = *lwork == -1;
#line 253 "cgees.f"
    wantvs = lsame_(jobvs, "V", (ftnlen)1, (ftnlen)1);
#line 254 "cgees.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 255 "cgees.f"
    if (! wantvs && ! lsame_(jobvs, "N", (ftnlen)1, (ftnlen)1)) {
#line 256 "cgees.f"
	*info = -1;
#line 257 "cgees.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 258 "cgees.f"
	*info = -2;
#line 259 "cgees.f"
    } else if (*n < 0) {
#line 260 "cgees.f"
	*info = -4;
#line 261 "cgees.f"
    } else if (*lda < max(1,*n)) {
#line 262 "cgees.f"
	*info = -6;
#line 263 "cgees.f"
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
#line 264 "cgees.f"
	*info = -10;
#line 265 "cgees.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by CHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 278 "cgees.f"
    if (*info == 0) {
#line 279 "cgees.f"
	if (*n == 0) {
#line 280 "cgees.f"
	    minwrk = 1;
#line 281 "cgees.f"
	    maxwrk = 1;
#line 282 "cgees.f"
	} else {
#line 283 "cgees.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "CGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);
#line 284 "cgees.f"
	    minwrk = *n << 1;

#line 286 "cgees.f"
	    chseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &w[1], &vs[
		    vs_offset], ldvs, &work[1], &c_n1, &ieval, (ftnlen)1, (
		    ftnlen)1);
#line 288 "cgees.f"
	    hswork = (integer) work[1].r;

#line 290 "cgees.f"
	    if (! wantvs) {
#line 291 "cgees.f"
		maxwrk = max(maxwrk,hswork);
#line 292 "cgees.f"
	    } else {
/* Computing MAX */
#line 293 "cgees.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 293 "cgees.f"
		maxwrk = max(i__1,i__2);
#line 295 "cgees.f"
		maxwrk = max(maxwrk,hswork);
#line 296 "cgees.f"
	    }
#line 297 "cgees.f"
	}
#line 298 "cgees.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 300 "cgees.f"
	if (*lwork < minwrk && ! lquery) {
#line 301 "cgees.f"
	    *info = -12;
#line 302 "cgees.f"
	}
#line 303 "cgees.f"
    }

#line 305 "cgees.f"
    if (*info != 0) {
#line 306 "cgees.f"
	i__1 = -(*info);
#line 306 "cgees.f"
	xerbla_("CGEES ", &i__1, (ftnlen)6);
#line 307 "cgees.f"
	return 0;
#line 308 "cgees.f"
    } else if (lquery) {
#line 309 "cgees.f"
	return 0;
#line 310 "cgees.f"
    }

/*     Quick return if possible */

#line 314 "cgees.f"
    if (*n == 0) {
#line 315 "cgees.f"
	*sdim = 0;
#line 316 "cgees.f"
	return 0;
#line 317 "cgees.f"
    }

/*     Get machine constants */

#line 321 "cgees.f"
    eps = slamch_("P", (ftnlen)1);
#line 322 "cgees.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 323 "cgees.f"
    bignum = 1. / smlnum;
#line 324 "cgees.f"
    slabad_(&smlnum, &bignum);
#line 325 "cgees.f"
    smlnum = sqrt(smlnum) / eps;
#line 326 "cgees.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 330 "cgees.f"
    anrm = clange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 331 "cgees.f"
    scalea = FALSE_;
#line 332 "cgees.f"
    if (anrm > 0. && anrm < smlnum) {
#line 333 "cgees.f"
	scalea = TRUE_;
#line 334 "cgees.f"
	cscale = smlnum;
#line 335 "cgees.f"
    } else if (anrm > bignum) {
#line 336 "cgees.f"
	scalea = TRUE_;
#line 337 "cgees.f"
	cscale = bignum;
#line 338 "cgees.f"
    }
#line 339 "cgees.f"
    if (scalea) {
#line 339 "cgees.f"
	clascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 339 "cgees.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (CWorkspace: none) */
/*     (RWorkspace: need N) */

#line 346 "cgees.f"
    ibal = 1;
#line 347 "cgees.f"
    cgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 353 "cgees.f"
    itau = 1;
#line 354 "cgees.f"
    iwrk = *n + itau;
#line 355 "cgees.f"
    i__1 = *lwork - iwrk + 1;
#line 355 "cgees.f"
    cgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 358 "cgees.f"
    if (wantvs) {

/*        Copy Householder vectors to VS */

#line 362 "cgees.f"
	clacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs, (ftnlen)1)
		;

/*        Generate unitary matrix in VS */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 368 "cgees.f"
	i__1 = *lwork - iwrk + 1;
#line 368 "cgees.f"
	cunghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk],
		 &i__1, &ierr);
#line 370 "cgees.f"
    }

#line 372 "cgees.f"
    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired */
/*     (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*     (RWorkspace: none) */

#line 378 "cgees.f"
    iwrk = itau;
#line 379 "cgees.f"
    i__1 = *lwork - iwrk + 1;
#line 379 "cgees.f"
    chseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vs[
	    vs_offset], ldvs, &work[iwrk], &i__1, &ieval, (ftnlen)1, (ftnlen)
	    1);
#line 381 "cgees.f"
    if (ieval > 0) {
#line 381 "cgees.f"
	*info = ieval;
#line 381 "cgees.f"
    }

/*     Sort eigenvalues if desired */

#line 386 "cgees.f"
    if (wantst && *info == 0) {
#line 387 "cgees.f"
	if (scalea) {
#line 387 "cgees.f"
	    clascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &w[1], n, &
		    ierr, (ftnlen)1);
#line 387 "cgees.f"
	}
#line 389 "cgees.f"
	i__1 = *n;
#line 389 "cgees.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 390 "cgees.f"
	    bwork[i__] = (*select)(&w[i__]);
#line 391 "cgees.f"
/* L10: */
#line 391 "cgees.f"
	}

/*        Reorder eigenvalues and transform Schur vectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: none) */

#line 397 "cgees.f"
	i__1 = *lwork - iwrk + 1;
#line 397 "cgees.f"
	ctrsen_("N", jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset], 
		ldvs, &w[1], sdim, &s, &sep, &work[iwrk], &i__1, &icond, (
		ftnlen)1, (ftnlen)1);
#line 399 "cgees.f"
    }

#line 401 "cgees.f"
    if (wantvs) {

/*        Undo balancing */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 407 "cgees.f"
	cgebak_("P", "R", n, &ilo, &ihi, &rwork[ibal], n, &vs[vs_offset], 
		ldvs, &ierr, (ftnlen)1, (ftnlen)1);
#line 409 "cgees.f"
    }

#line 411 "cgees.f"
    if (scalea) {

/*        Undo scaling for the Schur form of A */

#line 415 "cgees.f"
	clascl_("U", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 416 "cgees.f"
	i__1 = *lda + 1;
#line 416 "cgees.f"
	ccopy_(n, &a[a_offset], &i__1, &w[1], &c__1);
#line 417 "cgees.f"
    }

#line 419 "cgees.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 420 "cgees.f"
    return 0;

/*     End of CGEES */

} /* cgees_ */

