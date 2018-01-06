#line 1 "dgeesx.f"
/* dgeesx.f -- translated by f2c (version 20100827).
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

#line 1 "dgeesx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> DGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, */
/*                          WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, */
/*                          IWORK, LIWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVS, SENSE, SORT */
/*       INTEGER            INFO, LDA, LDVS, LIWORK, LWORK, N, SDIM */
/*       DOUBLE PRECISION   RCONDE, RCONDV */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
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
/* > DGEESX computes for an N-by-N real nonsymmetric matrix A, the */
/* > eigenvalues, the real Schur form T, and, optionally, the matrix of */
/* > Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > real Schur form so that selected eigenvalues are at the top left; */
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
/* > A real matrix is in real Schur form if it is upper quasi-triangular */
/* > with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in */
/* > the form */
/* >           [  a  b  ] */
/* >           [  c  a  ] */
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
/* >          SELECT(WR(j),WI(j)) is true; i.e., if either one of a */
/* >          complex conjugate pair of eigenvalues is selected, then both */
/* >          are.  Note that a selected complex eigenvalue may no longer */
/* >          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since */
/* >          ordering may change the value of complex eigenvalues */
/* >          (especially if the eigenvalue is ill-conditioned); in this */
/* >          case INFO may be set to N+3 (see INFO below). */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the N-by-N matrix A. */
/* >          On exit, A is overwritten by its real Schur form T. */
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
/* >          WR and WI contain the real and imaginary parts, respectively, */
/* >          of the computed eigenvalues, in the same order that they */
/* >          appear on the diagonal of the output Schur form T.  Complex */
/* >          conjugate pairs of eigenvalues appear consecutively with the */
/* >          eigenvalue having the positive imaginary part first. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,3*N). */
/* >          Also, if SENSE = 'E' or 'V' or 'B', */
/* >          LWORK >= N+2*SDIM*(N-SDIM), where SDIM is the number of */
/* >          selected eigenvalues computed by this routine.  Note that */
/* >          N+2*SDIM*(N-SDIM) <= N+N*N/2. Note also that an error is only */
/* >          returned if LWORK < max(1,3*N), but if SENSE = 'E' or 'V' or */
/* >          'B' this may not be large enough. */
/* >          For good performance, LWORK must generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates upper bounds on the optimal sizes of the */
/* >          arrays WORK and IWORK, returns these values as the first */
/* >          entries of the WORK and IWORK arrays, and no error messages */
/* >          related to LWORK or LIWORK are issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          LIWORK >= 1; if SENSE = 'V' or 'B', LIWORK >= SDIM*(N-SDIM). */
/* >          Note that SDIM*(N-SDIM) <= N*N/4. Note also that an error is */
/* >          only returned if LIWORK < 1, but if SENSE = 'V' or 'B' this */
/* >          may not be large enough. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates upper bounds on the optimal sizes of */
/* >          the arrays WORK and IWORK, returns these values as the first */
/* >          entries of the WORK and IWORK arrays, and no error messages */
/* >          related to LWORK or LIWORK are issued by XERBLA. */
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

/* > \ingroup doubleGEeigen */

/*  ===================================================================== */
/* Subroutine */ int dgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublereal *a, integer *lda, integer *sdim, 
	doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, 
	doublereal *rconde, doublereal *rcondv, doublereal *work, integer *
	lwork, integer *iwork, integer *liwork, logical *bwork, integer *info,
	 ftnlen jobvs_len, ftnlen sort_len, ftnlen sense_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, i1, i2, ip, ihi, ilo;
    static doublereal dum[1], eps;
    static integer ibal;
    static doublereal anrm;
    static integer ierr, itau, iwrk, lwrk, inxt, icond, ieval;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical cursl;
    static integer liwrk;
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
	    ftnlen, ftnlen);
    static logical wantsb;
    extern /* Subroutine */ int dtrsen_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen, ftnlen);
    static logical wantse, lastsl;
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    static integer hswork;
    static logical wantst, lquery, wantsv, wantvs;


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

#line 339 "dgeesx.f"
    /* Parameter adjustments */
#line 339 "dgeesx.f"
    a_dim1 = *lda;
#line 339 "dgeesx.f"
    a_offset = 1 + a_dim1;
#line 339 "dgeesx.f"
    a -= a_offset;
#line 339 "dgeesx.f"
    --wr;
#line 339 "dgeesx.f"
    --wi;
#line 339 "dgeesx.f"
    vs_dim1 = *ldvs;
#line 339 "dgeesx.f"
    vs_offset = 1 + vs_dim1;
#line 339 "dgeesx.f"
    vs -= vs_offset;
#line 339 "dgeesx.f"
    --work;
#line 339 "dgeesx.f"
    --iwork;
#line 339 "dgeesx.f"
    --bwork;
#line 339 "dgeesx.f"

#line 339 "dgeesx.f"
    /* Function Body */
#line 339 "dgeesx.f"
    *info = 0;
#line 340 "dgeesx.f"
    wantvs = lsame_(jobvs, "V", (ftnlen)1, (ftnlen)1);
#line 341 "dgeesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 342 "dgeesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 343 "dgeesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 344 "dgeesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 345 "dgeesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 346 "dgeesx.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 348 "dgeesx.f"
    if (! wantvs && ! lsame_(jobvs, "N", (ftnlen)1, (ftnlen)1)) {
#line 349 "dgeesx.f"
	*info = -1;
#line 350 "dgeesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 351 "dgeesx.f"
	*info = -2;
#line 352 "dgeesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 354 "dgeesx.f"
	*info = -4;
#line 355 "dgeesx.f"
    } else if (*n < 0) {
#line 356 "dgeesx.f"
	*info = -5;
#line 357 "dgeesx.f"
    } else if (*lda < max(1,*n)) {
#line 358 "dgeesx.f"
	*info = -7;
#line 359 "dgeesx.f"
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
#line 360 "dgeesx.f"
	*info = -12;
#line 361 "dgeesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "RWorkspace:" describe the */
/*       minimal amount of real workspace needed at that point in the */
/*       code, as well as the preferred amount for good performance. */
/*       IWorkspace refers to integer workspace. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by DHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case. */
/*       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed */
/*       depends on SDIM, which is computed by the routine DTRSEN later */
/*       in the code.) */

#line 377 "dgeesx.f"
    if (*info == 0) {
#line 378 "dgeesx.f"
	liwrk = 1;
#line 379 "dgeesx.f"
	if (*n == 0) {
#line 380 "dgeesx.f"
	    minwrk = 1;
#line 381 "dgeesx.f"
	    lwrk = 1;
#line 382 "dgeesx.f"
	} else {
#line 383 "dgeesx.f"
	    maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "DGEHRD", " ", n, &c__1, 
		    n, &c__0, (ftnlen)6, (ftnlen)1);
#line 384 "dgeesx.f"
	    minwrk = *n * 3;

#line 386 "dgeesx.f"
	    dhseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1]
		    , &vs[vs_offset], ldvs, &work[1], &c_n1, &ieval, (ftnlen)
		    1, (ftnlen)1);
#line 388 "dgeesx.f"
	    hswork = (integer) work[1];

#line 390 "dgeesx.f"
	    if (! wantvs) {
/* Computing MAX */
#line 391 "dgeesx.f"
		i__1 = maxwrk, i__2 = *n + hswork;
#line 391 "dgeesx.f"
		maxwrk = max(i__1,i__2);
#line 392 "dgeesx.f"
	    } else {
/* Computing MAX */
#line 393 "dgeesx.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"DORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
#line 393 "dgeesx.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 395 "dgeesx.f"
		i__1 = maxwrk, i__2 = *n + hswork;
#line 395 "dgeesx.f"
		maxwrk = max(i__1,i__2);
#line 396 "dgeesx.f"
	    }
#line 397 "dgeesx.f"
	    lwrk = maxwrk;
#line 398 "dgeesx.f"
	    if (! wantsn) {
/* Computing MAX */
#line 398 "dgeesx.f"
		i__1 = lwrk, i__2 = *n + *n * *n / 2;
#line 398 "dgeesx.f"
		lwrk = max(i__1,i__2);
#line 398 "dgeesx.f"
	    }
#line 400 "dgeesx.f"
	    if (wantsv || wantsb) {
#line 400 "dgeesx.f"
		liwrk = *n * *n / 4;
#line 400 "dgeesx.f"
	    }
#line 402 "dgeesx.f"
	}
#line 403 "dgeesx.f"
	iwork[1] = liwrk;
#line 404 "dgeesx.f"
	work[1] = (doublereal) lwrk;

#line 406 "dgeesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 407 "dgeesx.f"
	    *info = -16;
#line 408 "dgeesx.f"
	} else if (*liwork < 1 && ! lquery) {
#line 409 "dgeesx.f"
	    *info = -18;
#line 410 "dgeesx.f"
	}
#line 411 "dgeesx.f"
    }

#line 413 "dgeesx.f"
    if (*info != 0) {
#line 414 "dgeesx.f"
	i__1 = -(*info);
#line 414 "dgeesx.f"
	xerbla_("DGEESX", &i__1, (ftnlen)6);
#line 415 "dgeesx.f"
	return 0;
#line 416 "dgeesx.f"
    } else if (lquery) {
#line 417 "dgeesx.f"
	return 0;
#line 418 "dgeesx.f"
    }

/*     Quick return if possible */

#line 422 "dgeesx.f"
    if (*n == 0) {
#line 423 "dgeesx.f"
	*sdim = 0;
#line 424 "dgeesx.f"
	return 0;
#line 425 "dgeesx.f"
    }

/*     Get machine constants */

#line 429 "dgeesx.f"
    eps = dlamch_("P", (ftnlen)1);
#line 430 "dgeesx.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 431 "dgeesx.f"
    bignum = 1. / smlnum;
#line 432 "dgeesx.f"
    dlabad_(&smlnum, &bignum);
#line 433 "dgeesx.f"
    smlnum = sqrt(smlnum) / eps;
#line 434 "dgeesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 438 "dgeesx.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 439 "dgeesx.f"
    scalea = FALSE_;
#line 440 "dgeesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 441 "dgeesx.f"
	scalea = TRUE_;
#line 442 "dgeesx.f"
	cscale = smlnum;
#line 443 "dgeesx.f"
    } else if (anrm > bignum) {
#line 444 "dgeesx.f"
	scalea = TRUE_;
#line 445 "dgeesx.f"
	cscale = bignum;
#line 446 "dgeesx.f"
    }
#line 447 "dgeesx.f"
    if (scalea) {
#line 447 "dgeesx.f"
	dlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 447 "dgeesx.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (RWorkspace: need N) */

#line 453 "dgeesx.f"
    ibal = 1;
#line 454 "dgeesx.f"
    dgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (RWorkspace: need 3*N, prefer 2*N+N*NB) */

#line 459 "dgeesx.f"
    itau = *n + ibal;
#line 460 "dgeesx.f"
    iwrk = *n + itau;
#line 461 "dgeesx.f"
    i__1 = *lwork - iwrk + 1;
#line 461 "dgeesx.f"
    dgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 464 "dgeesx.f"
    if (wantvs) {

/*        Copy Householder vectors to VS */

#line 468 "dgeesx.f"
	dlacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VS */
/*        (RWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

#line 473 "dgeesx.f"
	i__1 = *lwork - iwrk + 1;
#line 473 "dgeesx.f"
	dorghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk],
		 &i__1, &ierr);
#line 475 "dgeesx.f"
    }

#line 477 "dgeesx.f"
    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired */
/*     (RWorkspace: need N+1, prefer N+HSWORK (see comments) ) */

#line 482 "dgeesx.f"
    iwrk = itau;
#line 483 "dgeesx.f"
    i__1 = *lwork - iwrk + 1;
#line 483 "dgeesx.f"
    dhseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vs[
	    vs_offset], ldvs, &work[iwrk], &i__1, &ieval, (ftnlen)1, (ftnlen)
	    1);
#line 485 "dgeesx.f"
    if (ieval > 0) {
#line 485 "dgeesx.f"
	*info = ieval;
#line 485 "dgeesx.f"
    }

/*     Sort eigenvalues if desired */

#line 490 "dgeesx.f"
    if (wantst && *info == 0) {
#line 491 "dgeesx.f"
	if (scalea) {
#line 492 "dgeesx.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wr[1], n, &
		    ierr, (ftnlen)1);
#line 493 "dgeesx.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wi[1], n, &
		    ierr, (ftnlen)1);
#line 494 "dgeesx.f"
	}
#line 495 "dgeesx.f"
	i__1 = *n;
#line 495 "dgeesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 496 "dgeesx.f"
	    bwork[i__] = (*select)(&wr[i__], &wi[i__]);
#line 497 "dgeesx.f"
/* L10: */
#line 497 "dgeesx.f"
	}

/*        Reorder eigenvalues, transform Schur vectors, and compute */
/*        reciprocal condition numbers */
/*        (RWorkspace: if SENSE is not 'N', need N+2*SDIM*(N-SDIM) */
/*                     otherwise, need N ) */
/*        (IWorkspace: if SENSE is 'V' or 'B', need SDIM*(N-SDIM) */
/*                     otherwise, need 0 ) */

#line 506 "dgeesx.f"
	i__1 = *lwork - iwrk + 1;
#line 506 "dgeesx.f"
	dtrsen_(sense, jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset],
		 ldvs, &wr[1], &wi[1], sdim, rconde, rcondv, &work[iwrk], &
		i__1, &iwork[1], liwork, &icond, (ftnlen)1, (ftnlen)1);
#line 509 "dgeesx.f"
	if (! wantsn) {
/* Computing MAX */
#line 509 "dgeesx.f"
	    i__1 = maxwrk, i__2 = *n + (*sdim << 1) * (*n - *sdim);
#line 509 "dgeesx.f"
	    maxwrk = max(i__1,i__2);
#line 509 "dgeesx.f"
	}
#line 511 "dgeesx.f"
	if (icond == -15) {

/*           Not enough real workspace */

#line 515 "dgeesx.f"
	    *info = -16;
#line 516 "dgeesx.f"
	} else if (icond == -17) {

/*           Not enough integer workspace */

#line 520 "dgeesx.f"
	    *info = -18;
#line 521 "dgeesx.f"
	} else if (icond > 0) {

/*           DTRSEN failed to reorder or to restore standard Schur form */

#line 525 "dgeesx.f"
	    *info = icond + *n;
#line 526 "dgeesx.f"
	}
#line 527 "dgeesx.f"
    }

#line 529 "dgeesx.f"
    if (wantvs) {

/*        Undo balancing */
/*        (RWorkspace: need N) */

#line 534 "dgeesx.f"
	dgebak_("P", "R", n, &ilo, &ihi, &work[ibal], n, &vs[vs_offset], ldvs,
		 &ierr, (ftnlen)1, (ftnlen)1);
#line 536 "dgeesx.f"
    }

#line 538 "dgeesx.f"
    if (scalea) {

/*        Undo scaling for the Schur form of A */

#line 542 "dgeesx.f"
	dlascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 543 "dgeesx.f"
	i__1 = *lda + 1;
#line 543 "dgeesx.f"
	dcopy_(n, &a[a_offset], &i__1, &wr[1], &c__1);
#line 544 "dgeesx.f"
	if ((wantsv || wantsb) && *info == 0) {
#line 545 "dgeesx.f"
	    dum[0] = *rcondv;
#line 546 "dgeesx.f"
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &
		    c__1, &ierr, (ftnlen)1);
#line 547 "dgeesx.f"
	    *rcondv = dum[0];
#line 548 "dgeesx.f"
	}
#line 549 "dgeesx.f"
	if (cscale == smlnum) {

/*           If scaling back towards underflow, adjust WI if an */
/*           offdiagonal element of a 2-by-2 block in the Schur form */
/*           underflows. */

#line 555 "dgeesx.f"
	    if (ieval > 0) {
#line 556 "dgeesx.f"
		i1 = ieval + 1;
#line 557 "dgeesx.f"
		i2 = ihi - 1;
#line 558 "dgeesx.f"
		i__1 = ilo - 1;
#line 558 "dgeesx.f"
		dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[
			1], n, &ierr, (ftnlen)1);
#line 560 "dgeesx.f"
	    } else if (wantst) {
#line 561 "dgeesx.f"
		i1 = 1;
#line 562 "dgeesx.f"
		i2 = *n - 1;
#line 563 "dgeesx.f"
	    } else {
#line 564 "dgeesx.f"
		i1 = ilo;
#line 565 "dgeesx.f"
		i2 = ihi - 1;
#line 566 "dgeesx.f"
	    }
#line 567 "dgeesx.f"
	    inxt = i1 - 1;
#line 568 "dgeesx.f"
	    i__1 = i2;
#line 568 "dgeesx.f"
	    for (i__ = i1; i__ <= i__1; ++i__) {
#line 569 "dgeesx.f"
		if (i__ < inxt) {
#line 569 "dgeesx.f"
		    goto L20;
#line 569 "dgeesx.f"
		}
#line 571 "dgeesx.f"
		if (wi[i__] == 0.) {
#line 572 "dgeesx.f"
		    inxt = i__ + 1;
#line 573 "dgeesx.f"
		} else {
#line 574 "dgeesx.f"
		    if (a[i__ + 1 + i__ * a_dim1] == 0.) {
#line 575 "dgeesx.f"
			wi[i__] = 0.;
#line 576 "dgeesx.f"
			wi[i__ + 1] = 0.;
#line 577 "dgeesx.f"
		    } else if (a[i__ + 1 + i__ * a_dim1] != 0. && a[i__ + (
			    i__ + 1) * a_dim1] == 0.) {
#line 579 "dgeesx.f"
			wi[i__] = 0.;
#line 580 "dgeesx.f"
			wi[i__ + 1] = 0.;
#line 581 "dgeesx.f"
			if (i__ > 1) {
#line 581 "dgeesx.f"
			    i__2 = i__ - 1;
#line 581 "dgeesx.f"
			    dswap_(&i__2, &a[i__ * a_dim1 + 1], &c__1, &a[(
				    i__ + 1) * a_dim1 + 1], &c__1);
#line 581 "dgeesx.f"
			}
#line 583 "dgeesx.f"
			if (*n > i__ + 1) {
#line 583 "dgeesx.f"
			    i__2 = *n - i__ - 1;
#line 583 "dgeesx.f"
			    dswap_(&i__2, &a[i__ + (i__ + 2) * a_dim1], lda, &
				    a[i__ + 1 + (i__ + 2) * a_dim1], lda);
#line 583 "dgeesx.f"
			}
#line 586 "dgeesx.f"
			dswap_(n, &vs[i__ * vs_dim1 + 1], &c__1, &vs[(i__ + 1)
				 * vs_dim1 + 1], &c__1);
#line 587 "dgeesx.f"
			a[i__ + (i__ + 1) * a_dim1] = a[i__ + 1 + i__ * 
				a_dim1];
#line 588 "dgeesx.f"
			a[i__ + 1 + i__ * a_dim1] = 0.;
#line 589 "dgeesx.f"
		    }
#line 590 "dgeesx.f"
		    inxt = i__ + 2;
#line 591 "dgeesx.f"
		}
#line 592 "dgeesx.f"
L20:
#line 592 "dgeesx.f"
		;
#line 592 "dgeesx.f"
	    }
#line 593 "dgeesx.f"
	}
#line 594 "dgeesx.f"
	i__1 = *n - ieval;
/* Computing MAX */
#line 594 "dgeesx.f"
	i__3 = *n - ieval;
#line 594 "dgeesx.f"
	i__2 = max(i__3,1);
#line 594 "dgeesx.f"
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[ieval + 
		1], &i__2, &ierr, (ftnlen)1);
#line 596 "dgeesx.f"
    }

#line 598 "dgeesx.f"
    if (wantst && *info == 0) {

/*        Check if reordering successful */

#line 602 "dgeesx.f"
	lastsl = TRUE_;
#line 603 "dgeesx.f"
	lst2sl = TRUE_;
#line 604 "dgeesx.f"
	*sdim = 0;
#line 605 "dgeesx.f"
	ip = 0;
#line 606 "dgeesx.f"
	i__1 = *n;
#line 606 "dgeesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 607 "dgeesx.f"
	    cursl = (*select)(&wr[i__], &wi[i__]);
#line 608 "dgeesx.f"
	    if (wi[i__] == 0.) {
#line 609 "dgeesx.f"
		if (cursl) {
#line 609 "dgeesx.f"
		    ++(*sdim);
#line 609 "dgeesx.f"
		}
#line 611 "dgeesx.f"
		ip = 0;
#line 612 "dgeesx.f"
		if (cursl && ! lastsl) {
#line 612 "dgeesx.f"
		    *info = *n + 2;
#line 612 "dgeesx.f"
		}
#line 614 "dgeesx.f"
	    } else {
#line 615 "dgeesx.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 619 "dgeesx.f"
		    cursl = cursl || lastsl;
#line 620 "dgeesx.f"
		    lastsl = cursl;
#line 621 "dgeesx.f"
		    if (cursl) {
#line 621 "dgeesx.f"
			*sdim += 2;
#line 621 "dgeesx.f"
		    }
#line 623 "dgeesx.f"
		    ip = -1;
#line 624 "dgeesx.f"
		    if (cursl && ! lst2sl) {
#line 624 "dgeesx.f"
			*info = *n + 2;
#line 624 "dgeesx.f"
		    }
#line 626 "dgeesx.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 630 "dgeesx.f"
		    ip = 1;
#line 631 "dgeesx.f"
		}
#line 632 "dgeesx.f"
	    }
#line 633 "dgeesx.f"
	    lst2sl = lastsl;
#line 634 "dgeesx.f"
	    lastsl = cursl;
#line 635 "dgeesx.f"
/* L30: */
#line 635 "dgeesx.f"
	}
#line 636 "dgeesx.f"
    }

#line 638 "dgeesx.f"
    work[1] = (doublereal) maxwrk;
#line 639 "dgeesx.f"
    if (wantsv || wantsb) {
/* Computing MAX */
#line 640 "dgeesx.f"
	i__1 = 1, i__2 = *sdim * (*n - *sdim);
#line 640 "dgeesx.f"
	iwork[1] = max(i__1,i__2);
#line 641 "dgeesx.f"
    } else {
#line 642 "dgeesx.f"
	iwork[1] = 1;
#line 643 "dgeesx.f"
    }

#line 645 "dgeesx.f"
    return 0;

/*     End of DGEESX */

} /* dgeesx_ */

