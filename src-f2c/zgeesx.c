#line 1 "zgeesx.f"
/* zgeesx.f -- translated by f2c (version 20100827).
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

#line 1 "zgeesx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W, */
/*                          VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK, */
/*                          BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVS, SENSE, SORT */
/*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM */
/*       DOUBLE PRECISION   RCONDE, RCONDV */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * ) */
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
/* > ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the */
/* > eigenvalues, the Schur form T, and, optionally, the matrix of Schur */
/* > vectors Z.  This gives the Schur factorization A = Z*T*(Z**H). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > Schur form so that selected eigenvalues are at the top left; */
/* > computes a reciprocal condition number for the average of the */
/* > selected eigenvalues (RCONDE); and computes a reciprocal condition */
/* > number for the right invariant subspace corresponding to the */
/* > selected eigenvalues (RCONDV).  The leading columns of Z form an */
/* > orthonormal basis for this invariant subspace. */
/* > */
/* > For further explanation of the reciprocal condition numbers RCONDE */
/* > and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where */
/* > these quantities are called s and sep respectively). */
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
/* >          = 'N': Eigenvalues are not ordered; */
/* >          = 'S': Eigenvalues are ordered (see SELECT). */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is a LOGICAL FUNCTION of one COMPLEX*16 argument */
/* >          SELECT must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'S', SELECT is used to select eigenvalues to order */
/* >          to the top left of the Schur form. */
/* >          If SORT = 'N', SELECT is not referenced. */
/* >          An eigenvalue W(j) is selected if SELECT(W(j)) is true. */
/* > \endverbatim */
/* > */
/* > \param[in] SENSE */
/* > \verbatim */
/* >          SENSE is CHARACTER*1 */
/* >          Determines which reciprocal condition numbers are computed. */
/* >          = 'N': None are computed; */
/* >          = 'E': Computed for average of selected eigenvalues only; */
/* >          = 'V': Computed for selected right invariant subspace only; */
/* >          = 'B': Computed for both. */
/* >          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'. */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A is overwritten by its Schur form T. */
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
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          W contains the computed eigenvalues, in the same order */
/* >          that they appear on the diagonal of the output Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* >          VS is COMPLEX*16 array, dimension (LDVS,N) */
/* >          If JOBVS = 'V', VS contains the unitary matrix Z of Schur */
/* >          vectors. */
/* >          If JOBVS = 'N', VS is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVS */
/* > \verbatim */
/* >          LDVS is INTEGER */
/* >          The leading dimension of the array VS.  LDVS >= 1, and if */
/* >          JOBVS = 'V', LDVS >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDE */
/* > \verbatim */
/* >          RCONDE is DOUBLE PRECISION */
/* >          If SENSE = 'E' or 'B', RCONDE contains the reciprocal */
/* >          condition number for the average of the selected eigenvalues. */
/* >          Not referenced if SENSE = 'N' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is DOUBLE PRECISION */
/* >          If SENSE = 'V' or 'B', RCONDV contains the reciprocal */
/* >          condition number for the selected right invariant subspace. */
/* >          Not referenced if SENSE = 'N' or 'E'. */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
/* >          Also, if SENSE = 'E' or 'V' or 'B', LWORK >= 2*SDIM*(N-SDIM), */
/* >          where SDIM is the number of selected eigenvalues computed by */
/* >          this routine.  Note that 2*SDIM*(N-SDIM) <= N*N/2. Note also */
/* >          that an error is only returned if LWORK < max(1,2*N), but if */
/* >          SENSE = 'E' or 'V' or 'B' this may not be large enough. */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates upper bound on the optimal size of the */
/* >          array WORK, returns this value as the first entry of the WORK */
/* >          array, and no error message related to LWORK is issued by */
/* >          XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
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
/* >             <= N: the QR algorithm failed to compute all the */
/* >                   eigenvalues; elements 1:ILO-1 and i+1:N of W */
/* >                   contain those eigenvalues which have converged; if */
/* >                   JOBVS = 'V', VS contains the transformation which */
/* >                   reduces A to its partially converged Schur form. */
/* >             = N+1: the eigenvalues could not be reordered because some */
/* >                   eigenvalues were too close to separate (the problem */
/* >                   is very ill-conditioned); */
/* >             = N+2: after reordering, roundoff changed values of some */
/* >                   complex eigenvalues so that leading eigenvalues in */
/* >                   the Schur form no longer satisfy SELECT=.TRUE.  This */
/* >                   could also be caused by underflow due to scaling. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, 
	doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *
	rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, 
	doublereal *rwork, logical *bwork, integer *info, ftnlen jobvs_len, 
	ftnlen sort_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ihi, ilo;
    static doublereal dum[1], eps;
    static integer ibal;
    static doublereal anrm;
    static integer ierr, itau, iwrk, lwrk, icond, ieval;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    static logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), zgebak_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen), zgebal_(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static logical wantsb, wantse;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    extern /* Subroutine */ int zhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer hswork;
    extern /* Subroutine */ int zunghr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static logical wantst, lquery, wantsv, wantvs;
    extern /* Subroutine */ int ztrsen_(char *, char *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 295 "zgeesx.f"
    /* Parameter adjustments */
#line 295 "zgeesx.f"
    a_dim1 = *lda;
#line 295 "zgeesx.f"
    a_offset = 1 + a_dim1;
#line 295 "zgeesx.f"
    a -= a_offset;
#line 295 "zgeesx.f"
    --w;
#line 295 "zgeesx.f"
    vs_dim1 = *ldvs;
#line 295 "zgeesx.f"
    vs_offset = 1 + vs_dim1;
#line 295 "zgeesx.f"
    vs -= vs_offset;
#line 295 "zgeesx.f"
    --work;
#line 295 "zgeesx.f"
    --rwork;
#line 295 "zgeesx.f"
    --bwork;
#line 295 "zgeesx.f"

#line 295 "zgeesx.f"
    /* Function Body */
#line 295 "zgeesx.f"
    *info = 0;
#line 296 "zgeesx.f"
    wantvs = lsame_(jobvs, "V", (ftnlen)1, (ftnlen)1);
#line 297 "zgeesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 298 "zgeesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 299 "zgeesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 300 "zgeesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 301 "zgeesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 302 "zgeesx.f"
    lquery = *lwork == -1;

#line 304 "zgeesx.f"
    if (! wantvs && ! lsame_(jobvs, "N", (ftnlen)1, (ftnlen)1)) {
#line 305 "zgeesx.f"
	*info = -1;
#line 306 "zgeesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 307 "zgeesx.f"
	*info = -2;
#line 308 "zgeesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 310 "zgeesx.f"
	*info = -4;
#line 311 "zgeesx.f"
    } else if (*n < 0) {
#line 312 "zgeesx.f"
	*info = -5;
#line 313 "zgeesx.f"
    } else if (*lda < max(1,*n)) {
#line 314 "zgeesx.f"
	*info = -7;
#line 315 "zgeesx.f"
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
#line 316 "zgeesx.f"
	*info = -11;
#line 317 "zgeesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of real workspace needed at that point in the */
/*       code, as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by ZHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case. */
/*       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed */
/*       depends on SDIM, which is computed by the routine ZTRSEN later */
/*       in the code.) */

#line 333 "zgeesx.f"
    if (*info == 0) {
#line 334 "zgeesx.f"
	if (*n == 0) {
#line 335 "zgeesx.f"
	    minwrk = 1;
#line 336 "zgeesx.f"
	    lwrk = 1;
#line 337 "zgeesx.f"
	} else {
#line 338 "zgeesx.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);
#line 339 "zgeesx.f"
	    minwrk = *n << 1;

#line 341 "zgeesx.f"
	    zhseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &w[1], &vs[
		    vs_offset], ldvs, &work[1], &c_n1, &ieval, (ftnlen)1, (
		    ftnlen)1);
#line 343 "zgeesx.f"
	    hswork = (integer) work[1].r;

#line 345 "zgeesx.f"
	    if (! wantvs) {
#line 346 "zgeesx.f"
		maxwrk = max(maxwrk,hswork);
#line 347 "zgeesx.f"
	    } else {
/* Computing MAX */
#line 348 "zgeesx.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 348 "zgeesx.f"
		maxwrk = max(i__1,i__2);
#line 350 "zgeesx.f"
		maxwrk = max(maxwrk,hswork);
#line 351 "zgeesx.f"
	    }
#line 352 "zgeesx.f"
	    lwrk = maxwrk;
#line 353 "zgeesx.f"
	    if (! wantsn) {
/* Computing MAX */
#line 353 "zgeesx.f"
		i__1 = lwrk, i__2 = *n * *n / 2;
#line 353 "zgeesx.f"
		lwrk = max(i__1,i__2);
#line 353 "zgeesx.f"
	    }
#line 355 "zgeesx.f"
	}
#line 356 "zgeesx.f"
	work[1].r = (doublereal) lwrk, work[1].i = 0.;

#line 358 "zgeesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 359 "zgeesx.f"
	    *info = -15;
#line 360 "zgeesx.f"
	}
#line 361 "zgeesx.f"
    }

#line 363 "zgeesx.f"
    if (*info != 0) {
#line 364 "zgeesx.f"
	i__1 = -(*info);
#line 364 "zgeesx.f"
	xerbla_("ZGEESX", &i__1, (ftnlen)6);
#line 365 "zgeesx.f"
	return 0;
#line 366 "zgeesx.f"
    } else if (lquery) {
#line 367 "zgeesx.f"
	return 0;
#line 368 "zgeesx.f"
    }

/*     Quick return if possible */

#line 372 "zgeesx.f"
    if (*n == 0) {
#line 373 "zgeesx.f"
	*sdim = 0;
#line 374 "zgeesx.f"
	return 0;
#line 375 "zgeesx.f"
    }

/*     Get machine constants */

#line 379 "zgeesx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 380 "zgeesx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 381 "zgeesx.f"
    bignum = 1. / smlnum;
#line 382 "zgeesx.f"
    dlabad_(&smlnum, &bignum);
#line 383 "zgeesx.f"
    smlnum = sqrt(smlnum) / eps;
#line 384 "zgeesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 388 "zgeesx.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 389 "zgeesx.f"
    scalea = FALSE_;
#line 390 "zgeesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 391 "zgeesx.f"
	scalea = TRUE_;
#line 392 "zgeesx.f"
	cscale = smlnum;
#line 393 "zgeesx.f"
    } else if (anrm > bignum) {
#line 394 "zgeesx.f"
	scalea = TRUE_;
#line 395 "zgeesx.f"
	cscale = bignum;
#line 396 "zgeesx.f"
    }
#line 397 "zgeesx.f"
    if (scalea) {
#line 397 "zgeesx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 397 "zgeesx.f"
    }


/*     Permute the matrix to make it more nearly triangular */
/*     (CWorkspace: none) */
/*     (RWorkspace: need N) */

#line 405 "zgeesx.f"
    ibal = 1;
#line 406 "zgeesx.f"
    zgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 412 "zgeesx.f"
    itau = 1;
#line 413 "zgeesx.f"
    iwrk = *n + itau;
#line 414 "zgeesx.f"
    i__1 = *lwork - iwrk + 1;
#line 414 "zgeesx.f"
    zgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 417 "zgeesx.f"
    if (wantvs) {

/*        Copy Householder vectors to VS */

#line 421 "zgeesx.f"
	zlacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs, (ftnlen)1)
		;

/*        Generate unitary matrix in VS */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 427 "zgeesx.f"
	i__1 = *lwork - iwrk + 1;
#line 427 "zgeesx.f"
	zunghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk],
		 &i__1, &ierr);
#line 429 "zgeesx.f"
    }

#line 431 "zgeesx.f"
    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired */
/*     (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*     (RWorkspace: none) */

#line 437 "zgeesx.f"
    iwrk = itau;
#line 438 "zgeesx.f"
    i__1 = *lwork - iwrk + 1;
#line 438 "zgeesx.f"
    zhseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vs[
	    vs_offset], ldvs, &work[iwrk], &i__1, &ieval, (ftnlen)1, (ftnlen)
	    1);
#line 440 "zgeesx.f"
    if (ieval > 0) {
#line 440 "zgeesx.f"
	*info = ieval;
#line 440 "zgeesx.f"
    }

/*     Sort eigenvalues if desired */

#line 445 "zgeesx.f"
    if (wantst && *info == 0) {
#line 446 "zgeesx.f"
	if (scalea) {
#line 446 "zgeesx.f"
	    zlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &w[1], n, &
		    ierr, (ftnlen)1);
#line 446 "zgeesx.f"
	}
#line 448 "zgeesx.f"
	i__1 = *n;
#line 448 "zgeesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 449 "zgeesx.f"
	    bwork[i__] = (*select)(&w[i__]);
#line 450 "zgeesx.f"
/* L10: */
#line 450 "zgeesx.f"
	}

/*        Reorder eigenvalues, transform Schur vectors, and compute */
/*        reciprocal condition numbers */
/*        (CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM) */
/*                     otherwise, need none ) */
/*        (RWorkspace: none) */

#line 458 "zgeesx.f"
	i__1 = *lwork - iwrk + 1;
#line 458 "zgeesx.f"
	ztrsen_(sense, jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset],
		 ldvs, &w[1], sdim, rconde, rcondv, &work[iwrk], &i__1, &
		icond, (ftnlen)1, (ftnlen)1);
#line 461 "zgeesx.f"
	if (! wantsn) {
/* Computing MAX */
#line 461 "zgeesx.f"
	    i__1 = maxwrk, i__2 = (*sdim << 1) * (*n - *sdim);
#line 461 "zgeesx.f"
	    maxwrk = max(i__1,i__2);
#line 461 "zgeesx.f"
	}
#line 463 "zgeesx.f"
	if (icond == -14) {

/*           Not enough complex workspace */

#line 467 "zgeesx.f"
	    *info = -15;
#line 468 "zgeesx.f"
	}
#line 469 "zgeesx.f"
    }

#line 471 "zgeesx.f"
    if (wantvs) {

/*        Undo balancing */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 477 "zgeesx.f"
	zgebak_("P", "R", n, &ilo, &ihi, &rwork[ibal], n, &vs[vs_offset], 
		ldvs, &ierr, (ftnlen)1, (ftnlen)1);
#line 479 "zgeesx.f"
    }

#line 481 "zgeesx.f"
    if (scalea) {

/*        Undo scaling for the Schur form of A */

#line 485 "zgeesx.f"
	zlascl_("U", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 486 "zgeesx.f"
	i__1 = *lda + 1;
#line 486 "zgeesx.f"
	zcopy_(n, &a[a_offset], &i__1, &w[1], &c__1);
#line 487 "zgeesx.f"
	if ((wantsv || wantsb) && *info == 0) {
#line 488 "zgeesx.f"
	    dum[0] = *rcondv;
#line 489 "zgeesx.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &
		    c__1, &ierr, (ftnlen)1);
#line 490 "zgeesx.f"
	    *rcondv = dum[0];
#line 491 "zgeesx.f"
	}
#line 492 "zgeesx.f"
    }

#line 494 "zgeesx.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 495 "zgeesx.f"
    return 0;

/*     End of ZGEESX */

} /* zgeesx_ */

