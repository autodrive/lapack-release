#line 1 "dstevr.f"
/* dstevr.f -- translated by f2c (version 20100827).
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

#line 1 "dstevr.f"
/* Table of constant values */

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;

/* > \brief <b> DSTEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEVR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstevr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstevr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstevr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, */
/*                          M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, RANGE */
/*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N */
/*       DOUBLE PRECISION   ABSTOL, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISUPPZ( * ), IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEVR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix T.  Eigenvalues and */
/* > eigenvectors can be selected by specifying either a range of values */
/* > or a range of indices for the desired eigenvalues. */
/* > */
/* > Whenever possible, DSTEVR calls DSTEMR to compute the */
/* > eigenspectrum using Relatively Robust Representations.  DSTEMR */
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
/* > Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of DSTEMR may create NaNs and infinities and */
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
/* >          For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and */
/* >          DSTEIN are called */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, D may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (max(1,N-1)) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A in elements 1 to N-1 of E. */
/* >          On exit, E may be multiplied by a constant factor chosen */
/* >          to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for eigenvalues. VL < VU. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
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
/* >          ABSTOL is DOUBLE PRECISION */
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
/* >          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          The first M elements contain the selected eigenvalues in */
/* >          ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal (and */
/* >          minimal) LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,20*N). */
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
/* >          The dimension of the array IWORK.  LIWORK >= max(1,10*N). */
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

/* > \ingroup doubleOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Inderjit Dhillon, IBM Almaden, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */
/* >     Ken Stanley, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dstevr_(char *jobz, char *range, integer *n, doublereal *
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
    static doublereal tnrm;
    static integer itmp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char order[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer lwmin;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static logical alleig, indeig;
    static integer iscale, ieeeok, indibl, indifl;
    static logical valeig;
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    static integer indisp;
    extern /* Subroutine */ int dstein_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    dsterf_(integer *, doublereal *, doublereal *, integer *);
    static integer indiwo;
    extern /* Subroutine */ int dstebz_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dstemr_(char *, char *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen);
    static integer liwmin;
    static logical tryrac;
    static integer nsplit;
    static doublereal smlnum;
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

#line 356 "dstevr.f"
    /* Parameter adjustments */
#line 356 "dstevr.f"
    --d__;
#line 356 "dstevr.f"
    --e;
#line 356 "dstevr.f"
    --w;
#line 356 "dstevr.f"
    z_dim1 = *ldz;
#line 356 "dstevr.f"
    z_offset = 1 + z_dim1;
#line 356 "dstevr.f"
    z__ -= z_offset;
#line 356 "dstevr.f"
    --isuppz;
#line 356 "dstevr.f"
    --work;
#line 356 "dstevr.f"
    --iwork;
#line 356 "dstevr.f"

#line 356 "dstevr.f"
    /* Function Body */
#line 356 "dstevr.f"
    ieeeok = ilaenv_(&c__10, "DSTEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 358 "dstevr.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 359 "dstevr.f"
    alleig = lsame_(range, "A", (ftnlen)1, (ftnlen)1);
#line 360 "dstevr.f"
    valeig = lsame_(range, "V", (ftnlen)1, (ftnlen)1);
#line 361 "dstevr.f"
    indeig = lsame_(range, "I", (ftnlen)1, (ftnlen)1);

#line 363 "dstevr.f"
    lquery = *lwork == -1 || *liwork == -1;
/* Computing MAX */
#line 364 "dstevr.f"
    i__1 = 1, i__2 = *n * 20;
#line 364 "dstevr.f"
    lwmin = max(i__1,i__2);
/* Computing MAX */
#line 365 "dstevr.f"
    i__1 = 1, i__2 = *n * 10;
#line 365 "dstevr.f"
    liwmin = max(i__1,i__2);


#line 368 "dstevr.f"
    *info = 0;
#line 369 "dstevr.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 370 "dstevr.f"
	*info = -1;
#line 371 "dstevr.f"
    } else if (! (alleig || valeig || indeig)) {
#line 372 "dstevr.f"
	*info = -2;
#line 373 "dstevr.f"
    } else if (*n < 0) {
#line 374 "dstevr.f"
	*info = -3;
#line 375 "dstevr.f"
    } else {
#line 376 "dstevr.f"
	if (valeig) {
#line 377 "dstevr.f"
	    if (*n > 0 && *vu <= *vl) {
#line 377 "dstevr.f"
		*info = -7;
#line 377 "dstevr.f"
	    }
#line 379 "dstevr.f"
	} else if (indeig) {
#line 380 "dstevr.f"
	    if (*il < 1 || *il > max(1,*n)) {
#line 381 "dstevr.f"
		*info = -8;
#line 382 "dstevr.f"
	    } else if (*iu < min(*n,*il) || *iu > *n) {
#line 383 "dstevr.f"
		*info = -9;
#line 384 "dstevr.f"
	    }
#line 385 "dstevr.f"
	}
#line 386 "dstevr.f"
    }
#line 387 "dstevr.f"
    if (*info == 0) {
#line 388 "dstevr.f"
	if (*ldz < 1 || wantz && *ldz < *n) {
#line 389 "dstevr.f"
	    *info = -14;
#line 390 "dstevr.f"
	}
#line 391 "dstevr.f"
    }

#line 393 "dstevr.f"
    if (*info == 0) {
#line 394 "dstevr.f"
	work[1] = (doublereal) lwmin;
#line 395 "dstevr.f"
	iwork[1] = liwmin;

#line 397 "dstevr.f"
	if (*lwork < lwmin && ! lquery) {
#line 398 "dstevr.f"
	    *info = -17;
#line 399 "dstevr.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 400 "dstevr.f"
	    *info = -19;
#line 401 "dstevr.f"
	}
#line 402 "dstevr.f"
    }

#line 404 "dstevr.f"
    if (*info != 0) {
#line 405 "dstevr.f"
	i__1 = -(*info);
#line 405 "dstevr.f"
	xerbla_("DSTEVR", &i__1, (ftnlen)6);
#line 406 "dstevr.f"
	return 0;
#line 407 "dstevr.f"
    } else if (lquery) {
#line 408 "dstevr.f"
	return 0;
#line 409 "dstevr.f"
    }

/*     Quick return if possible */

#line 413 "dstevr.f"
    *m = 0;
#line 414 "dstevr.f"
    if (*n == 0) {
#line 414 "dstevr.f"
	return 0;
#line 414 "dstevr.f"
    }

#line 417 "dstevr.f"
    if (*n == 1) {
#line 418 "dstevr.f"
	if (alleig || indeig) {
#line 419 "dstevr.f"
	    *m = 1;
#line 420 "dstevr.f"
	    w[1] = d__[1];
#line 421 "dstevr.f"
	} else {
#line 422 "dstevr.f"
	    if (*vl < d__[1] && *vu >= d__[1]) {
#line 423 "dstevr.f"
		*m = 1;
#line 424 "dstevr.f"
		w[1] = d__[1];
#line 425 "dstevr.f"
	    }
#line 426 "dstevr.f"
	}
#line 427 "dstevr.f"
	if (wantz) {
#line 427 "dstevr.f"
	    z__[z_dim1 + 1] = 1.;
#line 427 "dstevr.f"
	}
#line 429 "dstevr.f"
	return 0;
#line 430 "dstevr.f"
    }

/*     Get machine constants. */

#line 434 "dstevr.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 435 "dstevr.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 436 "dstevr.f"
    smlnum = safmin / eps;
#line 437 "dstevr.f"
    bignum = 1. / smlnum;
#line 438 "dstevr.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 439 "dstevr.f"
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
#line 439 "dstevr.f"
    rmax = min(d__1,d__2);


/*     Scale matrix to allowable range, if necessary. */

#line 444 "dstevr.f"
    iscale = 0;
#line 445 "dstevr.f"
    if (valeig) {
#line 446 "dstevr.f"
	vll = *vl;
#line 447 "dstevr.f"
	vuu = *vu;
#line 448 "dstevr.f"
    }

#line 450 "dstevr.f"
    tnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 451 "dstevr.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 452 "dstevr.f"
	iscale = 1;
#line 453 "dstevr.f"
	sigma = rmin / tnrm;
#line 454 "dstevr.f"
    } else if (tnrm > rmax) {
#line 455 "dstevr.f"
	iscale = 1;
#line 456 "dstevr.f"
	sigma = rmax / tnrm;
#line 457 "dstevr.f"
    }
#line 458 "dstevr.f"
    if (iscale == 1) {
#line 459 "dstevr.f"
	dscal_(n, &sigma, &d__[1], &c__1);
#line 460 "dstevr.f"
	i__1 = *n - 1;
#line 460 "dstevr.f"
	dscal_(&i__1, &sigma, &e[1], &c__1);
#line 461 "dstevr.f"
	if (valeig) {
#line 462 "dstevr.f"
	    vll = *vl * sigma;
#line 463 "dstevr.f"
	    vuu = *vu * sigma;
#line 464 "dstevr.f"
	}
#line 465 "dstevr.f"
    }
/*     Initialize indices into workspaces.  Note: These indices are used only */
/*     if DSTERF or DSTEMR fail. */
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
#line 472 "dstevr.f"
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and */
/*     stores the starting and finishing indices of each block. */
#line 475 "dstevr.f"
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     DSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
#line 480 "dstevr.f"
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
#line 482 "dstevr.f"
    indiwo = indisp + *n;

/*     If all eigenvalues are desired, then */
/*     call DSTERF or DSTEMR.  If this fails for some eigenvalue, then */
/*     try DSTEBZ. */


#line 489 "dstevr.f"
    test = FALSE_;
#line 490 "dstevr.f"
    if (indeig) {
#line 491 "dstevr.f"
	if (*il == 1 && *iu == *n) {
#line 492 "dstevr.f"
	    test = TRUE_;
#line 493 "dstevr.f"
	}
#line 494 "dstevr.f"
    }
#line 495 "dstevr.f"
    if ((alleig || test) && ieeeok == 1) {
#line 496 "dstevr.f"
	i__1 = *n - 1;
#line 496 "dstevr.f"
	dcopy_(&i__1, &e[1], &c__1, &work[1], &c__1);
#line 497 "dstevr.f"
	if (! wantz) {
#line 498 "dstevr.f"
	    dcopy_(n, &d__[1], &c__1, &w[1], &c__1);
#line 499 "dstevr.f"
	    dsterf_(n, &w[1], &work[1], info);
#line 500 "dstevr.f"
	} else {
#line 501 "dstevr.f"
	    dcopy_(n, &d__[1], &c__1, &work[*n + 1], &c__1);
#line 502 "dstevr.f"
	    if (*abstol <= *n * 2. * eps) {
#line 503 "dstevr.f"
		tryrac = TRUE_;
#line 504 "dstevr.f"
	    } else {
#line 505 "dstevr.f"
		tryrac = FALSE_;
#line 506 "dstevr.f"
	    }
#line 507 "dstevr.f"
	    i__1 = *lwork - (*n << 1);
#line 507 "dstevr.f"
	    dstemr_(jobz, "A", n, &work[*n + 1], &work[1], vl, vu, il, iu, m, 
		    &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &work[
		    (*n << 1) + 1], &i__1, &iwork[1], liwork, info, (ftnlen)1,
		     (ftnlen)1);

#line 511 "dstevr.f"
	}
#line 512 "dstevr.f"
	if (*info == 0) {
#line 513 "dstevr.f"
	    *m = *n;
#line 514 "dstevr.f"
	    goto L10;
#line 515 "dstevr.f"
	}
#line 516 "dstevr.f"
	*info = 0;
#line 517 "dstevr.f"
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN. */

#line 521 "dstevr.f"
    if (wantz) {
#line 522 "dstevr.f"
	*(unsigned char *)order = 'B';
#line 523 "dstevr.f"
    } else {
#line 524 "dstevr.f"
	*(unsigned char *)order = 'E';
#line 525 "dstevr.f"
    }
#line 527 "dstevr.f"
    dstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, &
	    nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[1], &iwork[
	    indiwo], info, (ftnlen)1, (ftnlen)1);

#line 531 "dstevr.f"
    if (wantz) {
#line 532 "dstevr.f"
	dstein_(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], &
		z__[z_offset], ldz, &work[1], &iwork[indiwo], &iwork[indifl], 
		info);
#line 535 "dstevr.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 539 "dstevr.f"
L10:
#line 540 "dstevr.f"
    if (iscale == 1) {
#line 541 "dstevr.f"
	if (*info == 0) {
#line 542 "dstevr.f"
	    imax = *m;
#line 543 "dstevr.f"
	} else {
#line 544 "dstevr.f"
	    imax = *info - 1;
#line 545 "dstevr.f"
	}
#line 546 "dstevr.f"
	d__1 = 1. / sigma;
#line 546 "dstevr.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 547 "dstevr.f"
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

#line 552 "dstevr.f"
    if (wantz) {
#line 553 "dstevr.f"
	i__1 = *m - 1;
#line 553 "dstevr.f"
	for (j = 1; j <= i__1; ++j) {
#line 554 "dstevr.f"
	    i__ = 0;
#line 555 "dstevr.f"
	    tmp1 = w[j];
#line 556 "dstevr.f"
	    i__2 = *m;
#line 556 "dstevr.f"
	    for (jj = j + 1; jj <= i__2; ++jj) {
#line 557 "dstevr.f"
		if (w[jj] < tmp1) {
#line 558 "dstevr.f"
		    i__ = jj;
#line 559 "dstevr.f"
		    tmp1 = w[jj];
#line 560 "dstevr.f"
		}
#line 561 "dstevr.f"
/* L20: */
#line 561 "dstevr.f"
	    }

#line 563 "dstevr.f"
	    if (i__ != 0) {
#line 564 "dstevr.f"
		itmp1 = iwork[i__];
#line 565 "dstevr.f"
		w[i__] = w[j];
#line 566 "dstevr.f"
		iwork[i__] = iwork[j];
#line 567 "dstevr.f"
		w[j] = tmp1;
#line 568 "dstevr.f"
		iwork[j] = itmp1;
#line 569 "dstevr.f"
		dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
#line 570 "dstevr.f"
	    }
#line 571 "dstevr.f"
/* L30: */
#line 571 "dstevr.f"
	}
#line 572 "dstevr.f"
    }

/*      Causes problems with tests 19 & 20: */
/*      IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002 */


#line 578 "dstevr.f"
    work[1] = (doublereal) lwmin;
#line 579 "dstevr.f"
    iwork[1] = liwmin;
#line 580 "dstevr.f"
    return 0;

/*     End of DSTEVR */

} /* dstevr_ */

