#line 1 "zhbgvd.f"
/* zhbgvd.f -- translated by f2c (version 20100827).
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

#line 1 "zhbgvd.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b ZHBGVD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHBGVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, */
/*                          Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LRWORK, */
/*      $                   LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBGVD computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite banded eigenproblem, of */
/* > the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian */
/* > and banded, and B is also positive definite.  If eigenvectors are */
/* > desired, it uses a divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
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
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangles of A and B are stored; */
/* >          = 'L':  Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KA */
/* > \verbatim */
/* >          KA is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'. KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of superdiagonals of the matrix B if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'. KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first ka+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). */
/* > */
/* >          On exit, the contents of AB are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BB */
/* > \verbatim */
/* >          BB is COMPLEX*16 array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**H*S, as returned by ZPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* >          LDBB is INTEGER */
/* >          The leading dimension of the array BB.  LDBB >= KB+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors, with the i-th column of Z holding the */
/* >          eigenvector associated with W(i). The eigenvectors are */
/* >          normalized so that Z**H*B*Z = I. */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO=0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK >= N. */
/* >          If JOBZ = 'V' and N > 1, LWORK >= 2*N**2. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
/* >          On exit, if INFO=0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of array RWORK. */
/* >          If N <= 1,               LRWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LRWORK >= N. */
/* >          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO=0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of array IWORK. */
/* >          If JOBZ = 'N' or N <= 1, LIWORK >= 1. */
/* >          If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, and i is: */
/* >             <= N:  the algorithm failed to converge: */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF */
/* >                    returned INFO = i: B is not positive definite. */
/* >                    The factorization of B could not be completed and */
/* >                    no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complex16OTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int zhbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info, ftnlen 
	jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;

    /* Local variables */
    static integer inde;
    static char vect[1];
    static integer llwk2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical upper;
    static integer llrwk;
    static logical wantz;
    static integer indwk2;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dsterf_(
	    integer *, doublereal *, doublereal *, integer *), zstedc_(char *,
	     integer *, doublereal *, doublereal *, doublecomplex *, integer *
	    , doublecomplex *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen), zhbtrd_(char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int zhbgst_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublereal *, 
	    integer *, ftnlen, ftnlen), zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer lrwmin;
    extern /* Subroutine */ int zpbstf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static logical lquery;


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 297 "zhbgvd.f"
    /* Parameter adjustments */
#line 297 "zhbgvd.f"
    ab_dim1 = *ldab;
#line 297 "zhbgvd.f"
    ab_offset = 1 + ab_dim1;
#line 297 "zhbgvd.f"
    ab -= ab_offset;
#line 297 "zhbgvd.f"
    bb_dim1 = *ldbb;
#line 297 "zhbgvd.f"
    bb_offset = 1 + bb_dim1;
#line 297 "zhbgvd.f"
    bb -= bb_offset;
#line 297 "zhbgvd.f"
    --w;
#line 297 "zhbgvd.f"
    z_dim1 = *ldz;
#line 297 "zhbgvd.f"
    z_offset = 1 + z_dim1;
#line 297 "zhbgvd.f"
    z__ -= z_offset;
#line 297 "zhbgvd.f"
    --work;
#line 297 "zhbgvd.f"
    --rwork;
#line 297 "zhbgvd.f"
    --iwork;
#line 297 "zhbgvd.f"

#line 297 "zhbgvd.f"
    /* Function Body */
#line 297 "zhbgvd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 298 "zhbgvd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 299 "zhbgvd.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 301 "zhbgvd.f"
    *info = 0;
#line 302 "zhbgvd.f"
    if (*n <= 1) {
#line 303 "zhbgvd.f"
	lwmin = *n + 1;
#line 304 "zhbgvd.f"
	lrwmin = *n + 1;
#line 305 "zhbgvd.f"
	liwmin = 1;
#line 306 "zhbgvd.f"
    } else if (wantz) {
/* Computing 2nd power */
#line 307 "zhbgvd.f"
	i__1 = *n;
#line 307 "zhbgvd.f"
	lwmin = i__1 * i__1 << 1;
/* Computing 2nd power */
#line 308 "zhbgvd.f"
	i__1 = *n;
#line 308 "zhbgvd.f"
	lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 309 "zhbgvd.f"
	liwmin = *n * 5 + 3;
#line 310 "zhbgvd.f"
    } else {
#line 311 "zhbgvd.f"
	lwmin = *n;
#line 312 "zhbgvd.f"
	lrwmin = *n;
#line 313 "zhbgvd.f"
	liwmin = 1;
#line 314 "zhbgvd.f"
    }
#line 315 "zhbgvd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 316 "zhbgvd.f"
	*info = -1;
#line 317 "zhbgvd.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 318 "zhbgvd.f"
	*info = -2;
#line 319 "zhbgvd.f"
    } else if (*n < 0) {
#line 320 "zhbgvd.f"
	*info = -3;
#line 321 "zhbgvd.f"
    } else if (*ka < 0) {
#line 322 "zhbgvd.f"
	*info = -4;
#line 323 "zhbgvd.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 324 "zhbgvd.f"
	*info = -5;
#line 325 "zhbgvd.f"
    } else if (*ldab < *ka + 1) {
#line 326 "zhbgvd.f"
	*info = -7;
#line 327 "zhbgvd.f"
    } else if (*ldbb < *kb + 1) {
#line 328 "zhbgvd.f"
	*info = -9;
#line 329 "zhbgvd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 330 "zhbgvd.f"
	*info = -12;
#line 331 "zhbgvd.f"
    }

#line 333 "zhbgvd.f"
    if (*info == 0) {
#line 334 "zhbgvd.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 335 "zhbgvd.f"
	rwork[1] = (doublereal) lrwmin;
#line 336 "zhbgvd.f"
	iwork[1] = liwmin;

#line 338 "zhbgvd.f"
	if (*lwork < lwmin && ! lquery) {
#line 339 "zhbgvd.f"
	    *info = -14;
#line 340 "zhbgvd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 341 "zhbgvd.f"
	    *info = -16;
#line 342 "zhbgvd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 343 "zhbgvd.f"
	    *info = -18;
#line 344 "zhbgvd.f"
	}
#line 345 "zhbgvd.f"
    }

#line 347 "zhbgvd.f"
    if (*info != 0) {
#line 348 "zhbgvd.f"
	i__1 = -(*info);
#line 348 "zhbgvd.f"
	xerbla_("ZHBGVD", &i__1, (ftnlen)6);
#line 349 "zhbgvd.f"
	return 0;
#line 350 "zhbgvd.f"
    } else if (lquery) {
#line 351 "zhbgvd.f"
	return 0;
#line 352 "zhbgvd.f"
    }

/*     Quick return if possible */

#line 356 "zhbgvd.f"
    if (*n == 0) {
#line 356 "zhbgvd.f"
	return 0;
#line 356 "zhbgvd.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 361 "zhbgvd.f"
    zpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 362 "zhbgvd.f"
    if (*info != 0) {
#line 363 "zhbgvd.f"
	*info = *n + *info;
#line 364 "zhbgvd.f"
	return 0;
#line 365 "zhbgvd.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 369 "zhbgvd.f"
    inde = 1;
#line 370 "zhbgvd.f"
    indwrk = inde + *n;
#line 371 "zhbgvd.f"
    indwk2 = *n * *n + 1;
#line 372 "zhbgvd.f"
    llwk2 = *lwork - indwk2 + 2;
#line 373 "zhbgvd.f"
    llrwk = *lrwork - indwrk + 2;
#line 374 "zhbgvd.f"
    zhbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &z__[z_offset], ldz, &work[1], &rwork[indwrk], &iinfo, (ftnlen)1,
	     (ftnlen)1);

/*     Reduce Hermitian band matrix to tridiagonal form. */

#line 379 "zhbgvd.f"
    if (wantz) {
#line 380 "zhbgvd.f"
	*(unsigned char *)vect = 'U';
#line 381 "zhbgvd.f"
    } else {
#line 382 "zhbgvd.f"
	*(unsigned char *)vect = 'N';
#line 383 "zhbgvd.f"
    }
#line 384 "zhbgvd.f"
    zhbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &w[1], &rwork[inde], &
	    z__[z_offset], ldz, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC. */

#line 389 "zhbgvd.f"
    if (! wantz) {
#line 390 "zhbgvd.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 391 "zhbgvd.f"
    } else {
#line 392 "zhbgvd.f"
	zstedc_("I", n, &w[1], &rwork[inde], &work[1], n, &work[indwk2], &
		llwk2, &rwork[indwrk], &llrwk, &iwork[1], liwork, info, (
		ftnlen)1);
#line 395 "zhbgvd.f"
	zgemm_("N", "N", n, n, n, &c_b1, &z__[z_offset], ldz, &work[1], n, &
		c_b2, &work[indwk2], n, (ftnlen)1, (ftnlen)1);
#line 397 "zhbgvd.f"
	zlacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz, (ftnlen)1);
#line 398 "zhbgvd.f"
    }

#line 400 "zhbgvd.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 401 "zhbgvd.f"
    rwork[1] = (doublereal) lrwmin;
#line 402 "zhbgvd.f"
    iwork[1] = liwmin;
#line 403 "zhbgvd.f"
    return 0;

/*     End of ZHBGVD */

} /* zhbgvd_ */

