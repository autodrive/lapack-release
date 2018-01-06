#line 1 "zgees.f"
/* zgees.f -- translated by f2c (version 20100827).
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

#line 1 "zgees.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> ZGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgees.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgees.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgees.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, */
/*                         LDVS, WORK, LWORK, RWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVS, SORT */
/*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM */
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
/* > ZGEES computes for an N-by-N complex nonsymmetric matrix A, the */
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
/* >          SELECT is a LOGICAL FUNCTION of one COMPLEX*16 argument */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          W contains the computed eigenvalues, in the same order that */
/* >          they appear on the diagonal of the output Schur form T. */
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
/* >          The leading dimension of the array VS.  LDVS >= 1; if */
/* >          JOBVS = 'V', LDVS >= N. */
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

/* > \date November 2011 */

/* > \ingroup complex16GEeigen */

/*  ===================================================================== */
/* Subroutine */ int zgees_(char *jobvs, char *sort, L_fp select, integer *n, 
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
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    static logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern /* Subroutine */ int zgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen, ftnlen), zgebal_(char *, integer *, 
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
	     integer *, integer *, ftnlen), zlacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int zhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer hswork;
    extern /* Subroutine */ int zunghr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static logical wantst, lquery, wantvs;
    extern /* Subroutine */ int ztrsen_(char *, char *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 251 "zgees.f"
    /* Parameter adjustments */
#line 251 "zgees.f"
    a_dim1 = *lda;
#line 251 "zgees.f"
    a_offset = 1 + a_dim1;
#line 251 "zgees.f"
    a -= a_offset;
#line 251 "zgees.f"
    --w;
#line 251 "zgees.f"
    vs_dim1 = *ldvs;
#line 251 "zgees.f"
    vs_offset = 1 + vs_dim1;
#line 251 "zgees.f"
    vs -= vs_offset;
#line 251 "zgees.f"
    --work;
#line 251 "zgees.f"
    --rwork;
#line 251 "zgees.f"
    --bwork;
#line 251 "zgees.f"

#line 251 "zgees.f"
    /* Function Body */
#line 251 "zgees.f"
    *info = 0;
#line 252 "zgees.f"
    lquery = *lwork == -1;
#line 253 "zgees.f"
    wantvs = lsame_(jobvs, "V", (ftnlen)1, (ftnlen)1);
#line 254 "zgees.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 255 "zgees.f"
    if (! wantvs && ! lsame_(jobvs, "N", (ftnlen)1, (ftnlen)1)) {
#line 256 "zgees.f"
	*info = -1;
#line 257 "zgees.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 258 "zgees.f"
	*info = -2;
#line 259 "zgees.f"
    } else if (*n < 0) {
#line 260 "zgees.f"
	*info = -4;
#line 261 "zgees.f"
    } else if (*lda < max(1,*n)) {
#line 262 "zgees.f"
	*info = -6;
#line 263 "zgees.f"
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
#line 264 "zgees.f"
	*info = -10;
#line 265 "zgees.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace to real */
/*       workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by ZHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 278 "zgees.f"
    if (*info == 0) {
#line 279 "zgees.f"
	if (*n == 0) {
#line 280 "zgees.f"
	    minwrk = 1;
#line 281 "zgees.f"
	    maxwrk = 1;
#line 282 "zgees.f"
	} else {
#line 283 "zgees.f"
	    maxwrk = *n + *n * ilaenv_(&c__1, "ZGEHRD", " ", n, &c__1, n, &
		    c__0, (ftnlen)6, (ftnlen)1);
#line 284 "zgees.f"
	    minwrk = *n << 1;

#line 286 "zgees.f"
	    zhseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &w[1], &vs[
		    vs_offset], ldvs, &work[1], &c_n1, &ieval, (ftnlen)1, (
		    ftnlen)1);
#line 288 "zgees.f"
	    hswork = (integer) work[1].r;

#line 290 "zgees.f"
	    if (! wantvs) {
#line 291 "zgees.f"
		maxwrk = max(maxwrk,hswork);
#line 292 "zgees.f"
	    } else {
/* Computing MAX */
#line 293 "zgees.f"
		i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_(&c__1, "ZUNGHR",
			 " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
#line 293 "zgees.f"
		maxwrk = max(i__1,i__2);
#line 295 "zgees.f"
		maxwrk = max(maxwrk,hswork);
#line 296 "zgees.f"
	    }
#line 297 "zgees.f"
	}
#line 298 "zgees.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 300 "zgees.f"
	if (*lwork < minwrk && ! lquery) {
#line 301 "zgees.f"
	    *info = -12;
#line 302 "zgees.f"
	}
#line 303 "zgees.f"
    }

#line 305 "zgees.f"
    if (*info != 0) {
#line 306 "zgees.f"
	i__1 = -(*info);
#line 306 "zgees.f"
	xerbla_("ZGEES ", &i__1, (ftnlen)6);
#line 307 "zgees.f"
	return 0;
#line 308 "zgees.f"
    } else if (lquery) {
#line 309 "zgees.f"
	return 0;
#line 310 "zgees.f"
    }

/*     Quick return if possible */

#line 314 "zgees.f"
    if (*n == 0) {
#line 315 "zgees.f"
	*sdim = 0;
#line 316 "zgees.f"
	return 0;
#line 317 "zgees.f"
    }

/*     Get machine constants */

#line 321 "zgees.f"
    eps = dlamch_("P", (ftnlen)1);
#line 322 "zgees.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 323 "zgees.f"
    bignum = 1. / smlnum;
#line 324 "zgees.f"
    dlabad_(&smlnum, &bignum);
#line 325 "zgees.f"
    smlnum = sqrt(smlnum) / eps;
#line 326 "zgees.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 330 "zgees.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 331 "zgees.f"
    scalea = FALSE_;
#line 332 "zgees.f"
    if (anrm > 0. && anrm < smlnum) {
#line 333 "zgees.f"
	scalea = TRUE_;
#line 334 "zgees.f"
	cscale = smlnum;
#line 335 "zgees.f"
    } else if (anrm > bignum) {
#line 336 "zgees.f"
	scalea = TRUE_;
#line 337 "zgees.f"
	cscale = bignum;
#line 338 "zgees.f"
    }
#line 339 "zgees.f"
    if (scalea) {
#line 339 "zgees.f"
	zlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 339 "zgees.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (CWorkspace: none) */
/*     (RWorkspace: need N) */

#line 346 "zgees.f"
    ibal = 1;
#line 347 "zgees.f"
    zgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (CWorkspace: need 2*N, prefer N+N*NB) */
/*     (RWorkspace: none) */

#line 353 "zgees.f"
    itau = 1;
#line 354 "zgees.f"
    iwrk = *n + itau;
#line 355 "zgees.f"
    i__1 = *lwork - iwrk + 1;
#line 355 "zgees.f"
    zgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 358 "zgees.f"
    if (wantvs) {

/*        Copy Householder vectors to VS */

#line 362 "zgees.f"
	zlacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs, (ftnlen)1)
		;

/*        Generate unitary matrix in VS */
/*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 368 "zgees.f"
	i__1 = *lwork - iwrk + 1;
#line 368 "zgees.f"
	zunghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk],
		 &i__1, &ierr);
#line 370 "zgees.f"
    }

#line 372 "zgees.f"
    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired */
/*     (CWorkspace: need 1, prefer HSWORK (see comments) ) */
/*     (RWorkspace: none) */

#line 378 "zgees.f"
    iwrk = itau;
#line 379 "zgees.f"
    i__1 = *lwork - iwrk + 1;
#line 379 "zgees.f"
    zhseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vs[
	    vs_offset], ldvs, &work[iwrk], &i__1, &ieval, (ftnlen)1, (ftnlen)
	    1);
#line 381 "zgees.f"
    if (ieval > 0) {
#line 381 "zgees.f"
	*info = ieval;
#line 381 "zgees.f"
    }

/*     Sort eigenvalues if desired */

#line 386 "zgees.f"
    if (wantst && *info == 0) {
#line 387 "zgees.f"
	if (scalea) {
#line 387 "zgees.f"
	    zlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &w[1], n, &
		    ierr, (ftnlen)1);
#line 387 "zgees.f"
	}
#line 389 "zgees.f"
	i__1 = *n;
#line 389 "zgees.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 390 "zgees.f"
	    bwork[i__] = (*select)(&w[i__]);
#line 391 "zgees.f"
/* L10: */
#line 391 "zgees.f"
	}

/*        Reorder eigenvalues and transform Schur vectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: none) */

#line 397 "zgees.f"
	i__1 = *lwork - iwrk + 1;
#line 397 "zgees.f"
	ztrsen_("N", jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset], 
		ldvs, &w[1], sdim, &s, &sep, &work[iwrk], &i__1, &icond, (
		ftnlen)1, (ftnlen)1);
#line 399 "zgees.f"
    }

#line 401 "zgees.f"
    if (wantvs) {

/*        Undo balancing */
/*        (CWorkspace: none) */
/*        (RWorkspace: need N) */

#line 407 "zgees.f"
	zgebak_("P", "R", n, &ilo, &ihi, &rwork[ibal], n, &vs[vs_offset], 
		ldvs, &ierr, (ftnlen)1, (ftnlen)1);
#line 409 "zgees.f"
    }

#line 411 "zgees.f"
    if (scalea) {

/*        Undo scaling for the Schur form of A */

#line 415 "zgees.f"
	zlascl_("U", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 416 "zgees.f"
	i__1 = *lda + 1;
#line 416 "zgees.f"
	zcopy_(n, &a[a_offset], &i__1, &w[1], &c__1);
#line 417 "zgees.f"
    }

#line 419 "zgees.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 420 "zgees.f"
    return 0;

/*     End of ZGEES */

} /* zgees_ */

