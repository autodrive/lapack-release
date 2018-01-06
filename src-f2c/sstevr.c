#line 1 "sstevr.f"
/* sstevr.f -- translated by f2c (version 20100827).
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

#line 1 "sstevr.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;

/* > \brief <b> SSTEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEVR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, */
/*                          M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N */
/*       REAL               ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEVR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix T.  Eigenvalues and */
/* > eigenvectors can be selected by specifying either a range of values */
/* > or a range of indices for the desired eigenvalues. */
/* > */
/* > Whenever possible, SSTEVR calls SSTEMR to compute the */
/* > eigenspectrum using Relatively Robust Representations.  SSTEMR */
/* > computes eigenvalues by the dqds algorithm, while orthogonal */
/* > eigenvectors are computed from various "good" L D L^T representations */
/* > (also known as Relatively Robust Representations). Gram-Schmidt */
/* > orthogonalization is avoided as far as possible. More specifically, */
/* > the various steps of the algorithm are as follows. For the i-th */
/* > unreduced block of T, */
/* >    (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T */
/* >         is a relatively robust representation, */
/* >    (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high */
/* >        relative accuracy by the dqds algorithm, */
/* >    (c) If there is a cluster of close eigenvalues, "choose" sigma_i */
/* >        close to the cluster, and go to step (a), */
/* >    (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T, */
/* >        compute the corresponding eigenvector by forming a */
/* >        rank-revealing twisted factorization. */
/* > The desired accuracy of the output can be specified by the input */
/* > parameter ABSTOL. */
/* > */
/* > For more details, see "A new O(n^2) algorithm for the symmetric */
/* > tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon, */
/* > Computer Science Division Technical Report No. UCB//CSD-97-971, */
/* > UC Berkeley, May 1997. */
/* > */
/* > */
/* > Note 1 : SSTEVR calls SSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > SSTEVR calls SSTEBZ and SSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of SSTEMR may create NaNs and infinities and */
/* > hence may abort due to a floating point exception in environments */
/* > which do not handle NaNs and infinities in the ieee standard default */
/* > manner. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only; */
/* >          = 'V':  Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RANGE */
/* > \verbatim */
/* >          RANGE is CHARACTER*1 */
/* >          = 'A': all eigenvalues will be found. */
/* >          = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* >                 will be found. */
/* >          = 'I': the IL-th through IU-th eigenvalues will be found. */
/* >          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and */
/* >          SSTEIN are called */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, D may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (max(1,N-1)) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A in elements 1 to N-1 of E. */
/* >          On exit, E may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is REAL */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          smallest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the */
/* >          largest eigenvalue to be returned. */
/* >          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* >          ABSTOL is REAL */
/* >          The absolute error tolerance for the eigenvalues. */
/* >          An approximate eigenvalue is accepted as converged */
/* >          when it is determined to lie in an interval [a,b] */
/* >          of width less than or equal to */
/* > */
/* >                  ABSTOL + EPS *   max( |a|,|b| ) , */
/* > */
/* >          where EPS is the machine precision.  If ABSTOL is less than */
/* >          or equal to zero, then  EPS*|T|  will be used in its place, */
/* >          where |T| is the 1-norm of the tridiagonal matrix obtained */
/* >          by reducing A to tridiagonal form. */
/* > */
/* >          See "Computing Small Singular Values of Bidiagonal Matrices */
/* >          with Guaranteed High Relative Accuracy," by Demmel and */
/* >          Kahan, LAPACK Working Note #3. */
/* > */
/* >          If high relative accuracy is important, set ABSTOL to */
/* >          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that */
/* >          eigenvalues are computed to high relative accuracy when */
/* >          possible in future releases.  The current code does not */
/* >          make any guarantees about high relative accuracy, but */
/* >          future releases will. See J. Barlow and J. Demmel, */
/* >          "Computing Accurate Eigensystems of Scaled Diagonally */
/* >          Dominant Matrices", LAPACK Working Note #7, for a discussion */
/* >          of which matrices define their eigenvalues to high relative */
/* >          accuracy. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The total number of eigenvalues found.  0 <= M <= N. */
/* >          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, max(1,M) ) */
/* >          If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* >          contain the orthonormal eigenvectors of the matrix A */
/* >          corresponding to the selected eigenvalues, with the i-th */
/* >          column of Z holding the eigenvector associated with W(i). */
/* >          Note: the user must ensure that at least max(1,M) columns are */
/* >          supplied in the array Z; if RANGE = 'V', the exact value of M */
/* >          is not known in advance and an upper bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* >          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) ) */
/* >          The support of the eigenvectors in Z, i.e., the indices */
/* >          indicating the nonzero elements in Z. The i-th eigenvector */
/* >          is nonzero only in elements ISUPPZ( 2*i-1 ) through */
/* >          ISUPPZ( 2*i ). */
/* >          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal (and */
/* >          minimal) LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= 20*N. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal sizes of the WORK and IWORK */
/* >          arrays, returns these values as the first entries of the WORK */
/* >          and IWORK arrays, and no error message related to LWORK or */
/* >          LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal (and */
/* >          minimal) LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK.  LIWORK >= 10*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK and IWORK arrays, and no error message related to */
/* >          LWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  Internal error */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup realOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Inderjit Dhillon, IBM Almaden, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Ken Stanley, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Jason Riedy, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sstevr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen range_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, jj;
    static doublereal eps, vll, vuu, tmp1;
    static integer imax;
    static doublereal rmin, rmax;
    static logical test;
    static doublereal tnrm, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static char order[1];
    static integer lwmin;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical wantz, alleig, indeig;
    static integer iscale, ieeeok, indibl, indifl;
    static logical valeig;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indisp, indiwo, liwmin;
    static logical tryrac;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int sstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *);
    static integer nsplit;
    extern /* Subroutine */ int sstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int sstemr_(char *, char *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, logical *, doublereal *, integer *, integer *, integer 
	    *, integer *, ftnlen, ftnlen);
    static logical lquery;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */


/*     Test the input parameters. */

#line 357 "sstevr.f"
    /* Parameter adjustments */
#line 357 "sstevr.f"
    --d__;
#line 357 "sstevr.f"
    --e;
#line 357 "sstevr.f"
    --w;
#line 357 "sstevr.f"
    z_dim1 = *ldz;
#line 357 "sstevr.f"
    z_offset = 1 + z_dim1;
#line 357 "sstevr.f"
    z__ -= z_offset;
#line 357 "sstevr.f"
    --isuppz;
#line 357 "sstevr.f"
    --work;
#line 357 "sstevr.f"
    --iwork;
#line 357 "sstevr.f"

#line 357 "sstevr.f"
    /* Function Body */
#line 357 "sstevr.f"
    ieeeok = ilaenv_(&c__10, "SSTEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 359 "sstevr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 360 "sstevr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 361 "sstevr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 362 "sstevr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 364 "sstevr.f"
    lquery = *lwork == -1 || *liwork == -1;
/* Computing MAX */
#line 365 "sstevr.f"
    i__1 = 1, i__2 = *n * 20;
#line 365 "sstevr.f"
    lwmin = max(i__1,i__2);
/* Computing MAX */
#line 366 "sstevr.f"
    i__1 = 1, i__2 = *n * 10;
#line 366 "sstevr.f"
    liwmin = max(i__1,i__2);


#line 369 "sstevr.f"
    *info = 0;
#line 370 "sstevr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 371 "sstevr.f"
	*info = -1;
#line 372 "sstevr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 373 "sstevr.f"
	*info = -2;
#line 374 "sstevr.f"
    } else if (*n < 0) {
#line 375 "sstevr.f"
	*info = -3;
#line 376 "sstevr.f"
    } else {
#line 377 "sstevr.f"
	if (valeig) {
#line 378 "sstevr.f"
	    if (*n > 0 && *vu <= *vl) {
#line 378 "sstevr.f"
		*info = -7;
#line 378 "sstevr.f"
	    }
#line 380 "sstevr.f"
	} else if (indeig) {
#line 381 "sstevr.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 382 "sstevr.f"
		*info = -8;
#line 383 "sstevr.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 384 "sstevr.f"
		*info = -9;
#line 385 "sstevr.f"
	    }
#line 386 "sstevr.f"
	}
#line 387 "sstevr.f"
    }
#line 388 "sstevr.f"
    if (*info == 0) {
#line 389 "sstevr.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 390 "sstevr.f"
	    *info = -14;
#line 391 "sstevr.f"
	}
#line 392 "sstevr.f"
    }

#line 394 "sstevr.f"
    if (*info == 0) {
#line 395 "sstevr.f"
	work[1] = (doublereal) lwmin;
#line 396 "sstevr.f"
	iwork[1] = liwmin;

#line 398 "sstevr.f"
	if (*lwork < lwmin && ! lquery) {
#line 399 "sstevr.f"
	    *info = -17;
#line 400 "sstevr.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 401 "sstevr.f"
	    *info = -19;
#line 402 "sstevr.f"
	}
#line 403 "sstevr.f"
    }

#line 405 "sstevr.f"
    if (*info != 0) {
#line 406 "sstevr.f"
	i__1 = -(*info);
#line 406 "sstevr.f"
	xerbla_("SSTEVR", &i__1, (ftnlen)6);
#line 407 "sstevr.f"
	return 0;
#line 408 "sstevr.f"
    } else if (lquery) {
#line 409 "sstevr.f"
	return 0;
#line 410 "sstevr.f"
    }

/*     Quick return if possible */

#line 414 "sstevr.f"
    *m = 0;
#line 415 "sstevr.f"
    if (*n == 0) {
#line 415 "sstevr.f"
	return 0;
#line 415 "sstevr.f"
    }

#line 418 "sstevr.f"
    if (*n == 1) {
#line 419 "sstevr.f"
	if (alleig || indeig) {
#line 420 "sstevr.f"
	    *m = 1;
#line 421 "sstevr.f"
	    w[1] = d__[1];
#line 422 "sstevr.f"
	} else {
#line 423 "sstevr.f"
	    if (*vl < d__[1] && *vu >= d__[1]) {
#line 424 "sstevr.f"
		*m = 1;
#line 425 "sstevr.f"
		w[1] = d__[1];
#line 426 "sstevr.f"
	    }
#line 427 "sstevr.f"
	}
#line 428 "sstevr.f"
	if (wantz) {
#line 428 "sstevr.f"
	    z__[z_dim1 + 1] = 1.;
#line 428 "sstevr.f"
	}
#line 430 "sstevr.f"
	return 0;
#line 431 "sstevr.f"
    }

/*     Get machine constants. */

#line 435 "sstevr.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 436 "sstevr.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 437 "sstevr.f"
    smlnum = safmin / eps;
#line 438 "sstevr.f"
    bignum = 1. / smlnum;
#line 439 "sstevr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 440 "sstevr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 440 "sstevr.f"
    rmax = min(d__1,d__2);


/*     Scale matrix to allowable range, if necessary. */

#line 445 "sstevr.f"
    iscale = 0;
#line 446 "sstevr.f"
    if (valeig) {
#line 447 "sstevr.f"
	vll = *vl;
#line 448 "sstevr.f"
	vuu = *vu;
#line 449 "sstevr.f"
    }

#line 451 "sstevr.f"
    tnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 452 "sstevr.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 453 "sstevr.f"
	iscale = 1;
#line 454 "sstevr.f"
	sigma = rmin / tnrm;
#line 455 "sstevr.f"
    } else if (tnrm > rmax) {
#line 456 "sstevr.f"
	iscale = 1;
#line 457 "sstevr.f"
	sigma = rmax / tnrm;
#line 458 "sstevr.f"
    }
#line 459 "sstevr.f"
    if (iscale == 1) {
#line 460 "sstevr.f"
	sscal_(n, &sigma, &d__[1], &c__1);
#line 461 "sstevr.f"
	i__1 = *n - 1;
#line 461 "sstevr.f"
	sscal_(&i__1, &sigma, &e[1], &c__1);
#line 462 "sstevr.f"
	if (valeig) {
#line 463 "sstevr.f"
	    vll = *vl * sigma;
#line 464 "sstevr.f"
	    vuu = *vu * sigma;
#line 465 "sstevr.f"
	}
#line 466 "sstevr.f"
    }
/*     Initialize indices into workspaces.  Note: These indices are used only */
/*     if SSTERF or SSTEMR fail. */
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 473 "sstevr.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 476 "sstevr.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     SSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 481 "sstevr.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 483 "sstevr.f"
    indiwo = indisp + *n;

/*     If all eigenvalues are desired, then */
/*     call SSTERF or SSTEMR.  If this fails for some eigenvalue, then */
/*     try SSTEBZ. */


#line 490 "sstevr.f"
    test = FALSE_;
#line 491 "sstevr.f"
    if (indeig) {
#line 492 "sstevr.f"
	if (*il == 1 && *iu == *n) {
#line 493 "sstevr.f"
	    test = TRUE_;
#line 494 "sstevr.f"
	}
#line 495 "sstevr.f"
    }
#line 496 "sstevr.f"
    if ((alleig || test) && ieeeok == 1) {
#line 497 "sstevr.f"
	i__1 = *n - 1;
#line 497 "sstevr.f"
	scopy_(&i__1, &e[1], &c__1, &work[1], &c__1);
#line 498 "sstevr.f"
	if (! wantz) {
#line 499 "sstevr.f"
	    scopy_(n, &d__[1], &c__1, &w[1], &c__1);
#line 500 "sstevr.f"
	    ssterf_(n, &w[1], &work[1], info);
#line 501 "sstevr.f"
	} else {
#line 502 "sstevr.f"
	    scopy_(n, &d__[1], &c__1, &work[*n + 1], &c__1);
#line 503 "sstevr.f"
	    if (*abstol <= *n * 2. * eps) {
#line 504 "sstevr.f"
		tryrac = TRUE_;
#line 505 "sstevr.f"
	    } else {
#line 506 "sstevr.f"
		tryrac = FALSE_;
#line 507 "sstevr.f"
	    }
#line 508 "sstevr.f"
	    i__1 = *lwork - (*n << 1);
#line 508 "sstevr.f"
	    sstemr_(jobz, "A", n, &work[*n + 1], &work[1], vl, vu, il, iu, m, 
		    &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &work[
		    (*n << 1) + 1], &i__1, &iwork[1], liwork, info, (ftnlen)1,
		     (ftnlen)1);

#line 512 "sstevr.f"
	}
#line 513 "sstevr.f"
	if (*info == 0) {
#line 514 "sstevr.f"
	    *m = *n;
#line 515 "sstevr.f"
	    goto L10;
#line 516 "sstevr.f"
	}
#line 517 "sstevr.f"
	*info = 0;
#line 518 "sstevr.f"
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */

#line 522 "sstevr.f"
    if (wantz) {
#line 523 "sstevr.f"
	*(unsigned char *)order = 'B';
#line 524 "sstevr.f"
    } else {
#line 525 "sstevr.f"
	*(unsigned char *)order = 'E';
#line 526 "sstevr.f"
    }
#line 528 "sstevr.f"
    sstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, &
	    nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[1], &iwork[
	    indiwo], info, (ftnlen)1, (ftnlen)1);

#line 532 "sstevr.f"
    if (wantz) {
#line 533 "sstevr.f"
	sstein_(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], &
		z__[z_offset], ldz, &work[1], &iwork[indiwo], &iwork[indifl], 
		info);
#line 536 "sstevr.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 540 "sstevr.f"
L10:
#line 541 "sstevr.f"
    if (iscale == 1) {
#line 542 "sstevr.f"
	if (*info == 0) {
#line 543 "sstevr.f"
	    imax = *m;
#line 544 "sstevr.f"
	} else {
#line 545 "sstevr.f"
	    imax = *info - 1;
#line 546 "sstevr.f"
	}
#line 547 "sstevr.f"
	d__1 = 1. / sigma;
#line 547 "sstevr.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 548 "sstevr.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 553 "sstevr.f"
    if (wantz) {
#line 554 "sstevr.f"
	i__1 = *m - 1;
#line 554 "sstevr.f"
	for (j = 1; j <= i__1; ++j) {
#line 555 "sstevr.f"
	    i__ = 0;
#line 556 "sstevr.f"
	    tmp1 = w[j];
#line 557 "sstevr.f"
	    i__2 = *m;
#line 557 "sstevr.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 558 "sstevr.f"
		if (w[jj] < tmp1) {
#line 559 "sstevr.f"
		    i__ = jj;
#line 560 "sstevr.f"
		    tmp1 = w[jj];
#line 561 "sstevr.f"
		}
#line 562 "sstevr.f"
/* L20: */
#line 562 "sstevr.f"
	    }

#line 564 "sstevr.f"
	    if (i__ != 0) {
#line 565 "sstevr.f"
		w[i__] = w[j];
#line 566 "sstevr.f"
		w[j] = tmp1;
#line 567 "sstevr.f"
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 568 "sstevr.f"
	    }
#line 569 "sstevr.f"
/* L30: */
#line 569 "sstevr.f"
	}
#line 570 "sstevr.f"
    }

/*      Causes problems with tests 19 & 20: */
/*      IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002 */


#line 576 "sstevr.f"
    work[1] = (doublereal) lwmin;
#line 577 "sstevr.f"
    iwork[1] = liwmin;
#line 578 "sstevr.f"
    return 0;

/*     End of SSTEVR */

} /* sstevr_ */

