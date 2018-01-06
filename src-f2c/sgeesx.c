#line 1 "sgeesx.f"
/* sgeesx.f -- translated by f2c (version 20100827).
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

#line 1 "sgeesx.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief <b> SGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors 
for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEESX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeesx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeesx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeesx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, */
/*                          WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK, */
/*                          IWORK, LIWORK, BWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVS, SENSE, SORT */
/*       INTEGER            INFO, LDA, LDVS, LIWORK, LWORK, N, SDIM */
/*       REAL               RCONDE, RCONDV */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            BWORK( * ) */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ), */
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
/* > SGEESX computes for an N-by-N real nonsymmetric matrix A, the */
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
/* >          SELECT is a LOGICAL FUNCTION of two REAL arguments */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (N) */
/* >          WR and WI contain the real and imaginary parts, respectively, */
/* >          of the computed eigenvalues, in the same order that they */
/* >          appear on the diagonal of the output Schur form T.  Complex */
/* >          conjugate pairs of eigenvalues appear consecutively with the */
/* >          eigenvalue having the positive imaginary part first. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* >          VS is REAL array, dimension (LDVS,N) */
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
/* >          RCONDE is REAL */
/* >          If SENSE = 'E' or 'B', RCONDE contains the reciprocal */
/* >          condition number for the average of the selected eigenvalues. */
/* >          Not referenced if SENSE = 'N' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* >          RCONDV is REAL */
/* >          If SENSE = 'V' or 'B', RCONDV contains the reciprocal */
/* >          condition number for the selected right invariant subspace. */
/* >          Not referenced if SENSE = 'N' or 'E'. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup realGEeigen */

/*  ===================================================================== */
/* Subroutine */ int sgeesx_(char *jobvs, char *sort, L_fp select, char *
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
    static logical cursl;
    static integer liwrk;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lst2sl;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    static logical scalea;
    static doublereal cscale;
    extern /* Subroutine */ int sgebak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sgebal_(char *, integer *, doublereal 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int sgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), slacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen);
    static logical wantsb, wantse, lastsl;
    extern /* Subroutine */ int sorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), shseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static logical wantsn;
    static doublereal smlnum;
    static integer hswork;
    extern /* Subroutine */ int strsen_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen, ftnlen);
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

#line 339 "sgeesx.f"
    /* Parameter adjustments */
#line 339 "sgeesx.f"
    a_dim1 = *lda;
#line 339 "sgeesx.f"
    a_offset = 1 + a_dim1;
#line 339 "sgeesx.f"
    a -= a_offset;
#line 339 "sgeesx.f"
    --wr;
#line 339 "sgeesx.f"
    --wi;
#line 339 "sgeesx.f"
    vs_dim1 = *ldvs;
#line 339 "sgeesx.f"
    vs_offset = 1 + vs_dim1;
#line 339 "sgeesx.f"
    vs -= vs_offset;
#line 339 "sgeesx.f"
    --work;
#line 339 "sgeesx.f"
    --iwork;
#line 339 "sgeesx.f"
    --bwork;
#line 339 "sgeesx.f"

#line 339 "sgeesx.f"
    /* Function Body */
#line 339 "sgeesx.f"
    *info = 0;
#line 340 "sgeesx.f"
    wantvs = lsame_(jobvs, "V", (ftnlen)1, (ftnlen)1);
#line 341 "sgeesx.f"
    wantst = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
#line 342 "sgeesx.f"
    wantsn = lsame_(sense, "N", (ftnlen)1, (ftnlen)1);
#line 343 "sgeesx.f"
    wantse = lsame_(sense, "E", (ftnlen)1, (ftnlen)1);
#line 344 "sgeesx.f"
    wantsv = lsame_(sense, "V", (ftnlen)1, (ftnlen)1);
#line 345 "sgeesx.f"
    wantsb = lsame_(sense, "B", (ftnlen)1, (ftnlen)1);
#line 346 "sgeesx.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 348 "sgeesx.f"
    if (! wantvs && ! lsame_(jobvs, "N", (ftnlen)1, (ftnlen)1)) {
#line 349 "sgeesx.f"
	*info = -1;
#line 350 "sgeesx.f"
    } else if (! wantst && ! lsame_(sort, "N", (ftnlen)1, (ftnlen)1)) {
#line 351 "sgeesx.f"
	*info = -2;
#line 352 "sgeesx.f"
    } else if (! (wantsn || wantse || wantsv || wantsb) || ! wantst && ! 
	    wantsn) {
#line 354 "sgeesx.f"
	*info = -4;
#line 355 "sgeesx.f"
    } else if (*n < 0) {
#line 356 "sgeesx.f"
	*info = -5;
#line 357 "sgeesx.f"
    } else if (*lda < max(1,*n)) {
#line 358 "sgeesx.f"
	*info = -7;
#line 359 "sgeesx.f"
    } else if (*ldvs < 1 || wantvs && *ldvs < *n) {
#line 360 "sgeesx.f"
	*info = -12;
#line 361 "sgeesx.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "RWorkspace:" describe the */
/*       minimal amount of real workspace needed at that point in the */
/*       code, as well as the preferred amount for good performance. */
/*       IWorkspace refers to integer workspace. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSWORK refers to the workspace preferred by SHSEQR, as */
/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
/*       the worst case. */
/*       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed */
/*       depends on SDIM, which is computed by the routine STRSEN later */
/*       in the code.) */

#line 377 "sgeesx.f"
    if (*info == 0) {
#line 378 "sgeesx.f"
	liwrk = 1;
#line 379 "sgeesx.f"
	if (*n == 0) {
#line 380 "sgeesx.f"
	    minwrk = 1;
#line 381 "sgeesx.f"
	    lwrk = 1;
#line 382 "sgeesx.f"
	} else {
#line 383 "sgeesx.f"
	    maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "SGEHRD", " ", n, &c__1, 
		    n, &c__0, (ftnlen)6, (ftnlen)1);
#line 384 "sgeesx.f"
	    minwrk = *n * 3;

#line 386 "sgeesx.f"
	    shseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1]
		    , &vs[vs_offset], ldvs, &work[1], &c_n1, &ieval, (ftnlen)
		    1, (ftnlen)1);
#line 388 "sgeesx.f"
	    hswork = (integer) work[1];

#line 390 "sgeesx.f"
	    if (! wantvs) {
/* Computing MAX */
#line 391 "sgeesx.f"
		i__1 = maxwrk, i__2 = *n + hswork;
#line 391 "sgeesx.f"
		maxwrk = max(i__1,i__2);
#line 392 "sgeesx.f"
	    } else {
/* Computing MAX */
#line 393 "sgeesx.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"SORGHR", " ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)
			1);
#line 393 "sgeesx.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 395 "sgeesx.f"
		i__1 = maxwrk, i__2 = *n + hswork;
#line 395 "sgeesx.f"
		maxwrk = max(i__1,i__2);
#line 396 "sgeesx.f"
	    }
#line 397 "sgeesx.f"
	    lwrk = maxwrk;
#line 398 "sgeesx.f"
	    if (! wantsn) {
/* Computing MAX */
#line 398 "sgeesx.f"
		i__1 = lwrk, i__2 = *n + *n * *n / 2;
#line 398 "sgeesx.f"
		lwrk = max(i__1,i__2);
#line 398 "sgeesx.f"
	    }
#line 400 "sgeesx.f"
	    if (wantsv || wantsb) {
#line 400 "sgeesx.f"
		liwrk = *n * *n / 4;
#line 400 "sgeesx.f"
	    }
#line 402 "sgeesx.f"
	}
#line 403 "sgeesx.f"
	iwork[1] = liwrk;
#line 404 "sgeesx.f"
	work[1] = (doublereal) lwrk;

#line 406 "sgeesx.f"
	if (*lwork < minwrk && ! lquery) {
#line 407 "sgeesx.f"
	    *info = -16;
#line 408 "sgeesx.f"
	} else if (*liwork < 1 && ! lquery) {
#line 409 "sgeesx.f"
	    *info = -18;
#line 410 "sgeesx.f"
	}
#line 411 "sgeesx.f"
    }

#line 413 "sgeesx.f"
    if (*info != 0) {
#line 414 "sgeesx.f"
	i__1 = -(*info);
#line 414 "sgeesx.f"
	xerbla_("SGEESX", &i__1, (ftnlen)6);
#line 415 "sgeesx.f"
	return 0;
#line 416 "sgeesx.f"
    } else if (lquery) {
#line 417 "sgeesx.f"
	return 0;
#line 418 "sgeesx.f"
    }

/*     Quick return if possible */

#line 422 "sgeesx.f"
    if (*n == 0) {
#line 423 "sgeesx.f"
	*sdim = 0;
#line 424 "sgeesx.f"
	return 0;
#line 425 "sgeesx.f"
    }

/*     Get machine constants */

#line 429 "sgeesx.f"
    eps = slamch_("P", (ftnlen)1);
#line 430 "sgeesx.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 431 "sgeesx.f"
    bignum = 1. / smlnum;
#line 432 "sgeesx.f"
    slabad_(&smlnum, &bignum);
#line 433 "sgeesx.f"
    smlnum = sqrt(smlnum) / eps;
#line 434 "sgeesx.f"
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 438 "sgeesx.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
#line 439 "sgeesx.f"
    scalea = FALSE_;
#line 440 "sgeesx.f"
    if (anrm > 0. && anrm < smlnum) {
#line 441 "sgeesx.f"
	scalea = TRUE_;
#line 442 "sgeesx.f"
	cscale = smlnum;
#line 443 "sgeesx.f"
    } else if (anrm > bignum) {
#line 444 "sgeesx.f"
	scalea = TRUE_;
#line 445 "sgeesx.f"
	cscale = bignum;
#line 446 "sgeesx.f"
    }
#line 447 "sgeesx.f"
    if (scalea) {
#line 447 "sgeesx.f"
	slascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 447 "sgeesx.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     (RWorkspace: need N) */

#line 453 "sgeesx.f"
    ibal = 1;
#line 454 "sgeesx.f"
    sgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form */
/*     (RWorkspace: need 3*N, prefer 2*N+N*NB) */

#line 459 "sgeesx.f"
    itau = *n + ibal;
#line 460 "sgeesx.f"
    iwrk = *n + itau;
#line 461 "sgeesx.f"
    i__1 = *lwork - iwrk + 1;
#line 461 "sgeesx.f"
    sgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1,
	     &ierr);

#line 464 "sgeesx.f"
    if (wantvs) {

/*        Copy Householder vectors to VS */

#line 468 "sgeesx.f"
	slacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs, (ftnlen)1)
		;

/*        Generate orthogonal matrix in VS */
/*        (RWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

#line 473 "sgeesx.f"
	i__1 = *lwork - iwrk + 1;
#line 473 "sgeesx.f"
	sorghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk],
		 &i__1, &ierr);
#line 475 "sgeesx.f"
    }

#line 477 "sgeesx.f"
    *sdim = 0;

/*     Perform QR iteration, accumulating Schur vectors in VS if desired */
/*     (RWorkspace: need N+1, prefer N+HSWORK (see comments) ) */

#line 482 "sgeesx.f"
    iwrk = itau;
#line 483 "sgeesx.f"
    i__1 = *lwork - iwrk + 1;
#line 483 "sgeesx.f"
    shseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vs[
	    vs_offset], ldvs, &work[iwrk], &i__1, &ieval, (ftnlen)1, (ftnlen)
	    1);
#line 485 "sgeesx.f"
    if (ieval > 0) {
#line 485 "sgeesx.f"
	*info = ieval;
#line 485 "sgeesx.f"
    }

/*     Sort eigenvalues if desired */

#line 490 "sgeesx.f"
    if (wantst && *info == 0) {
#line 491 "sgeesx.f"
	if (scalea) {
#line 492 "sgeesx.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wr[1], n, &
		    ierr, (ftnlen)1);
#line 493 "sgeesx.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wi[1], n, &
		    ierr, (ftnlen)1);
#line 494 "sgeesx.f"
	}
#line 495 "sgeesx.f"
	i__1 = *n;
#line 495 "sgeesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 496 "sgeesx.f"
	    bwork[i__] = (*select)(&wr[i__], &wi[i__]);
#line 497 "sgeesx.f"
/* L10: */
#line 497 "sgeesx.f"
	}

/*        Reorder eigenvalues, transform Schur vectors, and compute */
/*        reciprocal condition numbers */
/*        (RWorkspace: if SENSE is not 'N', need N+2*SDIM*(N-SDIM) */
/*                     otherwise, need N ) */
/*        (IWorkspace: if SENSE is 'V' or 'B', need SDIM*(N-SDIM) */
/*                     otherwise, need 0 ) */

#line 506 "sgeesx.f"
	i__1 = *lwork - iwrk + 1;
#line 506 "sgeesx.f"
	strsen_(sense, jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset],
		 ldvs, &wr[1], &wi[1], sdim, rconde, rcondv, &work[iwrk], &
		i__1, &iwork[1], liwork, &icond, (ftnlen)1, (ftnlen)1);
#line 509 "sgeesx.f"
	if (! wantsn) {
/* Computing MAX */
#line 509 "sgeesx.f"
	    i__1 = maxwrk, i__2 = *n + (*sdim << 1) * (*n - *sdim);
#line 509 "sgeesx.f"
	    maxwrk = max(i__1,i__2);
#line 509 "sgeesx.f"
	}
#line 511 "sgeesx.f"
	if (icond == -15) {

/*           Not enough real workspace */

#line 515 "sgeesx.f"
	    *info = -16;
#line 516 "sgeesx.f"
	} else if (icond == -17) {

/*           Not enough integer workspace */

#line 520 "sgeesx.f"
	    *info = -18;
#line 521 "sgeesx.f"
	} else if (icond > 0) {

/*           STRSEN failed to reorder or to restore standard Schur form */

#line 525 "sgeesx.f"
	    *info = icond + *n;
#line 526 "sgeesx.f"
	}
#line 527 "sgeesx.f"
    }

#line 529 "sgeesx.f"
    if (wantvs) {

/*        Undo balancing */
/*        (RWorkspace: need N) */

#line 534 "sgeesx.f"
	sgebak_("P", "R", n, &ilo, &ihi, &work[ibal], n, &vs[vs_offset], ldvs,
		 &ierr, (ftnlen)1, (ftnlen)1);
#line 536 "sgeesx.f"
    }

#line 538 "sgeesx.f"
    if (scalea) {

/*        Undo scaling for the Schur form of A */

#line 542 "sgeesx.f"
	slascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
#line 543 "sgeesx.f"
	i__1 = *lda + 1;
#line 543 "sgeesx.f"
	scopy_(n, &a[a_offset], &i__1, &wr[1], &c__1);
#line 544 "sgeesx.f"
	if ((wantsv || wantsb) && *info == 0) {
#line 545 "sgeesx.f"
	    dum[0] = *rcondv;
#line 546 "sgeesx.f"
	    slascl_("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &
		    c__1, &ierr, (ftnlen)1);
#line 547 "sgeesx.f"
	    *rcondv = dum[0];
#line 548 "sgeesx.f"
	}
#line 549 "sgeesx.f"
	if (cscale == smlnum) {

/*           If scaling back towards underflow, adjust WI if an */
/*           offdiagonal element of a 2-by-2 block in the Schur form */
/*           underflows. */

#line 555 "sgeesx.f"
	    if (ieval > 0) {
#line 556 "sgeesx.f"
		i1 = ieval + 1;
#line 557 "sgeesx.f"
		i2 = ihi - 1;
#line 558 "sgeesx.f"
		i__1 = ilo - 1;
#line 558 "sgeesx.f"
		slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[
			1], n, &ierr, (ftnlen)1);
#line 560 "sgeesx.f"
	    } else if (wantst) {
#line 561 "sgeesx.f"
		i1 = 1;
#line 562 "sgeesx.f"
		i2 = *n - 1;
#line 563 "sgeesx.f"
	    } else {
#line 564 "sgeesx.f"
		i1 = ilo;
#line 565 "sgeesx.f"
		i2 = ihi - 1;
#line 566 "sgeesx.f"
	    }
#line 567 "sgeesx.f"
	    inxt = i1 - 1;
#line 568 "sgeesx.f"
	    i__1 = i2;
#line 568 "sgeesx.f"
	    for (i__ = i1; i__ <= i__1; ++i__) {
#line 569 "sgeesx.f"
		if (i__ < inxt) {
#line 569 "sgeesx.f"
		    goto L20;
#line 569 "sgeesx.f"
		}
#line 571 "sgeesx.f"
		if (wi[i__] == 0.) {
#line 572 "sgeesx.f"
		    inxt = i__ + 1;
#line 573 "sgeesx.f"
		} else {
#line 574 "sgeesx.f"
		    if (a[i__ + 1 + i__ * a_dim1] == 0.) {
#line 575 "sgeesx.f"
			wi[i__] = 0.;
#line 576 "sgeesx.f"
			wi[i__ + 1] = 0.;
#line 577 "sgeesx.f"
		    } else if (a[i__ + 1 + i__ * a_dim1] != 0. && a[i__ + (
			    i__ + 1) * a_dim1] == 0.) {
#line 579 "sgeesx.f"
			wi[i__] = 0.;
#line 580 "sgeesx.f"
			wi[i__ + 1] = 0.;
#line 581 "sgeesx.f"
			if (i__ > 1) {
#line 581 "sgeesx.f"
			    i__2 = i__ - 1;
#line 581 "sgeesx.f"
			    sswap_(&i__2, &a[i__ * a_dim1 + 1], &c__1, &a[(
				    i__ + 1) * a_dim1 + 1], &c__1);
#line 581 "sgeesx.f"
			}
#line 583 "sgeesx.f"
			if (*n > i__ + 1) {
#line 583 "sgeesx.f"
			    i__2 = *n - i__ - 1;
#line 583 "sgeesx.f"
			    sswap_(&i__2, &a[i__ + (i__ + 2) * a_dim1], lda, &
				    a[i__ + 1 + (i__ + 2) * a_dim1], lda);
#line 583 "sgeesx.f"
			}
#line 586 "sgeesx.f"
			sswap_(n, &vs[i__ * vs_dim1 + 1], &c__1, &vs[(i__ + 1)
				 * vs_dim1 + 1], &c__1);
#line 587 "sgeesx.f"
			a[i__ + (i__ + 1) * a_dim1] = a[i__ + 1 + i__ * 
				a_dim1];
#line 588 "sgeesx.f"
			a[i__ + 1 + i__ * a_dim1] = 0.;
#line 589 "sgeesx.f"
		    }
#line 590 "sgeesx.f"
		    inxt = i__ + 2;
#line 591 "sgeesx.f"
		}
#line 592 "sgeesx.f"
L20:
#line 592 "sgeesx.f"
		;
#line 592 "sgeesx.f"
	    }
#line 593 "sgeesx.f"
	}
#line 594 "sgeesx.f"
	i__1 = *n - ieval;
/* Computing MAX */
#line 594 "sgeesx.f"
	i__3 = *n - ieval;
#line 594 "sgeesx.f"
	i__2 = max(i__3,1);
#line 594 "sgeesx.f"
	slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[ieval + 
		1], &i__2, &ierr, (ftnlen)1);
#line 596 "sgeesx.f"
    }

#line 598 "sgeesx.f"
    if (wantst && *info == 0) {

/*        Check if reordering successful */

#line 602 "sgeesx.f"
	lastsl = TRUE_;
#line 603 "sgeesx.f"
	lst2sl = TRUE_;
#line 604 "sgeesx.f"
	*sdim = 0;
#line 605 "sgeesx.f"
	ip = 0;
#line 606 "sgeesx.f"
	i__1 = *n;
#line 606 "sgeesx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 607 "sgeesx.f"
	    cursl = (*select)(&wr[i__], &wi[i__]);
#line 608 "sgeesx.f"
	    if (wi[i__] == 0.) {
#line 609 "sgeesx.f"
		if (cursl) {
#line 609 "sgeesx.f"
		    ++(*sdim);
#line 609 "sgeesx.f"
		}
#line 611 "sgeesx.f"
		ip = 0;
#line 612 "sgeesx.f"
		if (cursl && ! lastsl) {
#line 612 "sgeesx.f"
		    *info = *n + 2;
#line 612 "sgeesx.f"
		}
#line 614 "sgeesx.f"
	    } else {
#line 615 "sgeesx.f"
		if (ip == 1) {

/*                 Last eigenvalue of conjugate pair */

#line 619 "sgeesx.f"
		    cursl = cursl || lastsl;
#line 620 "sgeesx.f"
		    lastsl = cursl;
#line 621 "sgeesx.f"
		    if (cursl) {
#line 621 "sgeesx.f"
			*sdim += 2;
#line 621 "sgeesx.f"
		    }
#line 623 "sgeesx.f"
		    ip = -1;
#line 624 "sgeesx.f"
		    if (cursl && ! lst2sl) {
#line 624 "sgeesx.f"
			*info = *n + 2;
#line 624 "sgeesx.f"
		    }
#line 626 "sgeesx.f"
		} else {

/*                 First eigenvalue of conjugate pair */

#line 630 "sgeesx.f"
		    ip = 1;
#line 631 "sgeesx.f"
		}
#line 632 "sgeesx.f"
	    }
#line 633 "sgeesx.f"
	    lst2sl = lastsl;
#line 634 "sgeesx.f"
	    lastsl = cursl;
#line 635 "sgeesx.f"
/* L30: */
#line 635 "sgeesx.f"
	}
#line 636 "sgeesx.f"
    }

#line 638 "sgeesx.f"
    work[1] = (doublereal) maxwrk;
#line 639 "sgeesx.f"
    if (wantsv || wantsb) {
#line 640 "sgeesx.f"
	iwork[1] = *sdim * (*n - *sdim);
#line 641 "sgeesx.f"
    } else {
#line 642 "sgeesx.f"
	iwork[1] = 1;
#line 643 "sgeesx.f"
    }

#line 645 "sgeesx.f"
    return 0;

/*     End of SGEESX */

} /* sgeesx_ */

