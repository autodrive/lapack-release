#line 1 "dgees.f"
/* dgees.f -- translated by f2c (version 20100827).
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

#line 1 "dgees.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> DGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f
or GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgees.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgees.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgees.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI, */
/*                         VS, LDVS, WORK, LWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVS, SORT */
/*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ), */
/*      $                   WR( * ) */
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
/* > DGEES computes for an N-by-N real nonsymmetric matrix A, the */
/* > eigenvalues, the real Schur form T, and, optionally, the matrix of */
/* > Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > real Schur form so that selected eigenvalues are at the top left. */
/* > The leading columns of Z then form an orthonormal basis for the */
/* > invariant subspace corresponding to the selected eigenvalues. */
/* > */
/* > A matrix is in real Schur form if it is upper quasi-triangular with */
/* > 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the */
/* > form */
/* >         [  a  b  ] */
/* >         [  c  a  ] */
/* > */
/* > where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc). */
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
/* >          SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments */
/* >          SELECT must be declared EXTERNAL in the calling subroutine. */
/* >          If SORT = 'S', SELECT is used to select eigenvalues to sort */
/* >          to the top left of the Schur form. */
/* >          If SORT = 'N', SELECT is not referenced. */
/* >          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if */
/* >          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex */
/* >          conjugate pair of eigenvalues is selected, then both complex */
/* >          eigenvalues are selected. */
/* >          Note that a selected complex eigenvalue may no longer */
/* >          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since */
/* >          ordering may change the value of complex eigenvalues */
/* >          (especially if the eigenvalue is ill-conditioned); in this */
/* >          case INFO is set to N+2 (see INFO below). */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A has been overwritten by its real Schur form T. */
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
/* >          If SORT = 'S', SDIM = number of eigenvalues (after sorting) */
/* >                         for which SELECT is true. (Complex conjugate */
/* >                         pairs for which SELECT is true for either */
/* >                         eigenvalue count as 2.) */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* >          WR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION array, dimension (N) */
/* >          WR and WI contain the real and imaginary parts, */
/* >          respectively, of the computed eigenvalues in the same order */
/* >          that they appear on the diagonal of the output Schur form T. */
/* >          Complex conjugate pairs of eigenvalues will appear */
/* >          consecutively with the eigenvalue having the positive */
/* >          imaginary part first. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* >          VS is DOUBLE PRECISION array, dimension (LDVS,N) */
/* >          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) contains the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,3*N). */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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
/* >                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI */
/* >                   contain those eigenvalues which have converged; if */
/* >                   JOBVS = 'V', VS contains the matrix which reduces A */
/* >                   to its partially converged Schur form. */
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

/* > \date November 2011 */

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	doublereal *a, integer *lda, integer *sdim, doublereal *wr, 
	doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, 
	integer *lwork, logical *bwork, integer *info, ftnlen jobvs_len, 
	ftnlen sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s;
    static integer i1, i2, ip, ihi, ilo;
    static doublereal dum[1], eps, sep;
    static integer ibal;
    static doublereal anrm;
    static integer idum[1], ierr, itau, iwrk, inxt, icond, ieval;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical cursl;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dgebal_(char *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static logical lst2sl, scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dhseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtrsen_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen, ftnlen);
    static logical lastsl;
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static integer hswork;
    static logical wantst, lquery, wantvs;


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

#line 272 "dgees.f"
    /* Parameter adjustments */
#line 272 "dgees.f"
    a_dim1 = *lda;
#line 272 "dgees.f"
    a_offset = 1 + a_dim1;
#line 272 "dgees.f"
    a -= a_offset;
#line 272 "dgees.f"
    --wr;
#line 272 "dgees.f"
    --wi;
#line 272 "dgees.f"
    vs_dim1 = *ldvs;
#line 272 "dgees.f"
    vs_offset = 1 + vs_dim1;
#line 272 "dgees.f"
    vs -= vs_offset;
#line 272 "dgees.f"
    --work;
#line 272 "dgees.f"
    --bwork;
#line 272 "dgees.f"

#line 272 "dgees.f"
    /* Function Body */
#line 272 "dgees.f"
    *info = 0;
#line 273 "dgees.f"
    lquery = *lwork == -1;
#line 274 "dgees.f"
    wantvs = lsame_(jobvs, "V", (ftnlen)1, (ftnlen)1);
#line 275 "dgees.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 276 "dgees.f"
    if (! wantvs && ! lsame_(jobvs, "N", (ftnlen)1, (ftnlen)1)) {
#line 277 "dgees.f"
	*info = -1;
#line 278 "dgees.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 279 "dgees.f"
	*info = -2;
#line 280 "dgees.f"
    } else if (*n < 0) {
#line 281 "dgees.f"
	*info = -4;
#line 282 "dgees.f"
    } else if (*lda < max(1,*n)) {
#line 283 "dgees.f"
	*info = -6;
#line 284 "dgees.f"
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
#line 285 "dgees.f"
	*info = -11;
#line 286 "dgees.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by DHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

#line 298 "dgees.f"
    if (*info == 0) {
#line 299 "dgees.f"
	if (*n == 0) {
#line 300 "dgees.f"
	    minwrk = 1;
#line 301 "dgees.f"
	    maxwrk = 1;
#line 302 "dgees.f"
	} else {
#line 303 "dgees.f"
	    maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "DGEHRD", " ", n, &c__1, 
		    n, &c__0, (ftnlen)6, (ftnlen)1);
#line 304 "dgees.f"
	    minwrk = *n * 3;

#line 306 "dgees.f"
	    dhseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1]
		    , &vs[vs_offset], ldvs, &work[1], &c_n1, &ieval, (ftnlen)
		    1, (ftnlen)1);
#line 308 "dgees.f"
	    hswork = (integer) work[1];

#line 310 "dgees.f"
	    if (! wantvs) {
/* Computing MAX */
#line 311 "dgees.f"
		i__1 = maxwrk, i__2 = *n + hswork;
#line 311 "dgees.f"
		maxwrk = max(i__1,i__2);
#line 312 "dgees.f"
	    } else {
/* Computing MAX */
#line 313 "dgees.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"DORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
#line 313 "dgees.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 315 "dgees.f"
		i__1 = maxwrk, i__2 = *n + hswork;
#line 315 "dgees.f"
		maxwrk = max(i__1,i__2);
#line 316 "dgees.f"
	    }
#line 317 "dgees.f"
	}
#line 318 "dgees.f"
	work[1] = (doublereal) maxwrk;

#line 320 "dgees.f"
	if (*lwork < minwrk && ! lquery) {
#line 321 "dgees.f"
	    *info = -13;
#line 322 "dgees.f"
	}
#line 323 "dgees.f"
    }

#line 325 "dgees.f"
    if (*info != 0) {
#line 326 "dgees.f"
	i__1 = -(*info);
#line 326 "dgees.f"
	xerbla_("DGEES ", &i__1, (ftnlen)6);
#line 327 "dgees.f"
	return 0;
#line 328 "dgees.f"
    } else if (lquery) {
#line 329 "dgees.f"
	return 0;
#line 330 "dgees.f"
    }

/*     Quick return if possible */

#line 334 "dgees.f"
    if (*n == 0) {
#line 335 "dgees.f"
	*sdim = 0;
#line 336 "dgees.f"
	return 0;
#line 337 "dgees.f"
    }

/*     Get machine constants */

#line 341 "dgees.f"
    eps = dlamch_("P", (ftnlen)1);
#line 342 "dgees.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 343 "dgees.f"
    bignum = 1. / smlnum;
#line 344 "dgees.f"
    dlabad_(&smlnum, &bignum);
#line 345 "dgees.f"
    smlnum = sqrt(smlnum) / eps;
#line 346 "dgees.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 350 "dgees.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 351 "dgees.f"
    scalea = FALSE_;
#line 352 "dgees.f"
    if (anrm > 0. && anrm < smlnum) {
#line 353 "dgees.f"
	scalea = TRUE_;
#line 354 "dgees.f"
	cscale = smlnum;
#line 355 "dgees.f"
    } else if (anrm > bignum) {
#line 356 "dgees.f"
	scalea = TRUE_;
#line 357 "dgees.f"
	cscale = bignum;
#line 358 "dgees.f"
    }
#line 359 "dgees.f"
    if (scalea) {
#line 359 "dgees.f"
	dlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 359 "dgees.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (Workspace: need N) */

#line 365 "dgees.f"
    ibal = 1;
#line 366 "dgees.f"
    dgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (Workspace: need 3*N, prefer 2*N+N*NB) */

#line 371 "dgees.f"
    itau = *n + ibal;
#line 372 "dgees.f"
    iwrk = *n + itau;
#line 373 "dgees.f"
    i__1 = *lwork - iwrk + 1;
#line 373 "dgees.f"
    dgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 376 "dgees.f"
    if (wantvs) {

/*        Copy Householder vectors to VS */

#line 380 "dgees.f"
	dlacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VS */
/*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

#line 385 "dgees.f"
	i__1 = *lwork - iwrk + 1;
#line 385 "dgees.f"
	dorghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk],
		 &i__1, &ierr);
#line 387 "dgees.f"
    }

#line 389 "dgees.f"
    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired */
/*     (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

#line 394 "dgees.f"
    iwrk = itau;
#line 395 "dgees.f"
    i__1 = *lwork - iwrk + 1;
#line 395 "dgees.f"
    dhseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vs[
	    vs_offset], ldvs, &work[iwrk], &i__1, &ieval, (ftnlen)1, (ftnlen)
	    1);
#line 397 "dgees.f"
    if (ieval > 0) {
#line 397 "dgees.f"
	*info = ieval;
#line 397 "dgees.f"
    }

/*     Sort eigenvalues if desired */

#line 402 "dgees.f"
    if (wantst && *info == 0) {
#line 403 "dgees.f"
	if (scalea) {
#line 404 "dgees.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wr[1], n, &
		    ierr, (ftnlen)1);
#line 405 "dgees.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wi[1], n, &
		    ierr, (ftnlen)1);
#line 406 "dgees.f"
	}
#line 407 "dgees.f"
	i__1 = *n;
#line 407 "dgees.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 408 "dgees.f"
	    bwork[i__] = (*select)(&wr[i__], &wi[i__]);
#line 409 "dgees.f"
/* L10: */
#line 409 "dgees.f"
	}

/*        Reorder eigenvalues and transform Schur vectors */
/*        (Workspace: none needed) */

#line 414 "dgees.f"
	i__1 = *lwork - iwrk + 1;
#line 414 "dgees.f"
	dtrsen_("N", jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset], 
		ldvs, &wr[1], &wi[1], sdim, &s, &sep, &work[iwrk], &i__1, 
		idum, &c__1, &icond, (ftnlen)1, (ftnlen)1);
#line 417 "dgees.f"
	if (icond > 0) {
#line 417 "dgees.f"
	    *info = *n + icond;
#line 417 "dgees.f"
	}
#line 419 "dgees.f"
    }

#line 421 "dgees.f"
    if (wantvs) {

/*        Undo balancing */
/*        (Workspace: need N) */

#line 426 "dgees.f"
	dgebak_("P", "R", n, &ilo, &ihi, &work[ibal], n, &vs[vs_offset], ldvs,
		 &ierr, (ftnlen)1, (ftnlen)1);
#line 428 "dgees.f"
    }

#line 430 "dgees.f"
    if (scalea) {

/*        Undo scaling for the Schur form of A */

#line 434 "dgees.f"
	dlascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 435 "dgees.f"
	i__1 = *lda + 1;
#line 435 "dgees.f"
	dcopy_(n, &a[a_offset], &i__1, &wr[1], &c__1);
#line 436 "dgees.f"
	if (cscale == smlnum) {

/*           If scaling back towards underflow, adjust WI if an */
/*           offdiagonal element of a 2-by-2 block in the Schur form */
/*           underflows. */

#line 442 "dgees.f"
	    if (ieval > 0) {
#line 443 "dgees.f"
		i1 = ieval + 1;
#line 444 "dgees.f"
		i2 = ihi - 1;
#line 445 "dgees.f"
		i__1 = ilo - 1;
/* Computing MAX */
#line 445 "dgees.f"
		i__3 = ilo - 1;
#line 445 "dgees.f"
		i__2 = max(i__3,1);
#line 445 "dgees.f"
		dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[
			1], &i__2, &ierr, (ftnlen)1);
#line 447 "dgees.f"
	    } else if (wantst) {
#line 448 "dgees.f"
		i1 = 1;
#line 449 "dgees.f"
		i2 = *n - 1;
#line 450 "dgees.f"
	    } else {
#line 451 "dgees.f"
		i1 = ilo;
#line 452 "dgees.f"
		i2 = ihi - 1;
#line 453 "dgees.f"
	    }
#line 454 "dgees.f"
	    inxt = i1 - 1;
#line 455 "dgees.f"
	    i__1 = i2;
#line 455 "dgees.f"
	    for (i__ = i1; i__ <= i__1; ++i__) {
#line 456 "dgees.f"
		if (i__ < inxt) {
#line 456 "dgees.f"
		    goto L20;
#line 456 "dgees.f"
		}
#line 458 "dgees.f"
		if (wi[i__] == 0.) {
#line 459 "dgees.f"
		    inxt = i__ + 1;
#line 460 "dgees.f"
		} else {
#line 461 "dgees.f"
		    if (a[i__ + 1 + i__ * a_dim1] == 0.) {
#line 462 "dgees.f"
			wi[i__] = 0.;
#line 463 "dgees.f"
			wi[i__ + 1] = 0.;
#line 464 "dgees.f"
		    } else if (a[i__ + 1 + i__ * a_dim1] != 0. && a[i__ + (
			    i__ + 1) * a_dim1] == 0.) {
#line 466 "dgees.f"
			wi[i__] = 0.;
#line 467 "dgees.f"
			wi[i__ + 1] = 0.;
#line 468 "dgees.f"
			if (i__ > 1) {
#line 468 "dgees.f"
			    i__2 = i__ - 1;
#line 468 "dgees.f"
			    dswap_(&i__2, &a[i__ * a_dim1 + 1], &c__1, &a[(
				    i__ + 1) * a_dim1 + 1], &c__1);
#line 468 "dgees.f"
			}
#line 470 "dgees.f"
			if (*n > i__ + 1) {
#line 470 "dgees.f"
			    i__2 = *n - i__ - 1;
#line 470 "dgees.f"
			    dswap_(&i__2, &a[i__ + (i__ + 2) * a_dim1], lda, &
				    a[i__ + 1 + (i__ + 2) * a_dim1], lda);
#line 470 "dgees.f"
			}
#line 473 "dgees.f"
			if (wantvs) {
#line 474 "dgees.f"
			    dswap_(n, &vs[i__ * vs_dim1 + 1], &c__1, &vs[(i__ 
				    + 1) * vs_dim1 + 1], &c__1);
#line 475 "dgees.f"
			}
#line 476 "dgees.f"
			a[i__ + (i__ + 1) * a_dim1] = a[i__ + 1 + i__ * 
				a_dim1];
#line 477 "dgees.f"
			a[i__ + 1 + i__ * a_dim1] = 0.;
#line 478 "dgees.f"
		    }
#line 479 "dgees.f"
		    inxt = i__ + 2;
#line 480 "dgees.f"
		}
#line 481 "dgees.f"
L20:
#line 481 "dgees.f"
		;
#line 481 "dgees.f"
	    }
#line 482 "dgees.f"
	}

/*        Undo scaling for the imaginary part of the eigenvalues */

#line 486 "dgees.f"
	i__1 = *n - ieval;
/* Computing MAX */
#line 486 "dgees.f"
	i__3 = *n - ieval;
#line 486 "dgees.f"
	i__2 = max(i__3,1);
#line 486 "dgees.f"
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[ieval + 
		1], &i__2, &ierr, (ftnlen)1);
#line 488 "dgees.f"
    }

#line 490 "dgees.f"
    if (wantst && *info == 0) {

/*        Check if reordering successful */

#line 494 "dgees.f"
	lastsl = TRUE_;
#line 495 "dgees.f"
	lst2sl = TRUE_;
#line 496 "dgees.f"
	*sdim = 0;
#line 497 "dgees.f"
	ip = 0;
#line 498 "dgees.f"
	i__1 = *n;
#line 498 "dgees.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 499 "dgees.f"
	    cursl = (*select)(&wr[i__], &wi[i__]);
#line 500 "dgees.f"
	    if (wi[i__] == 0.) {
#line 501 "dgees.f"
		if (cursl) {
#line 501 "dgees.f"
		    ++(*sdim);
#line 501 "dgees.f"
		}
#line 503 "dgees.f"
		ip = 0;
#line 504 "dgees.f"
		if (cursl && ! lastsl) {
#line 504 "dgees.f"
		    *info = *n + 2;
#line 504 "dgees.f"
		}
#line 506 "dgees.f"
	    } else {
#line 507 "dgees.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 511 "dgees.f"
		    cursl = cursl || lastsl;
#line 512 "dgees.f"
		    lastsl = cursl;
#line 513 "dgees.f"
		    if (cursl) {
#line 513 "dgees.f"
			*sdim += 2;
#line 513 "dgees.f"
		    }
#line 515 "dgees.f"
		    ip = -1;
#line 516 "dgees.f"
		    if (cursl && ! lst2sl) {
#line 516 "dgees.f"
			*info = *n + 2;
#line 516 "dgees.f"
		    }
#line 518 "dgees.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 522 "dgees.f"
		    ip = 1;
#line 523 "dgees.f"
		}
#line 524 "dgees.f"
	    }
#line 525 "dgees.f"
	    lst2sl = lastsl;
#line 526 "dgees.f"
	    lastsl = cursl;
#line 527 "dgees.f"
/* L30: */
#line 527 "dgees.f"
	}
#line 528 "dgees.f"
    }

#line 530 "dgees.f"
    work[1] = (doublereal) maxwrk;
#line 531 "dgees.f"
    return 0;

/*     End of DGEES */

} /* dgees_ */

