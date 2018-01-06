#line 1 "cheevd.f"
/* cheevd.f -- translated by f2c (version 20100827).
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

#line 1 "cheevd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b18 = 1.;

/* > \brief <b> CHEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, */
/*                          LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEVD computes all eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix A.  If eigenvectors are desired, it uses a */
/* > divide and conquer algorithm. */
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
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* >          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* >          orthonormal eigenvectors of the matrix A. */
/* >          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') */
/* >          or the upper triangle (if UPLO='U') of A, including the */
/* >          diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
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
/* >          The length of the array WORK. */
/* >          If N <= 1,                LWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1. */
/* >          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2. */
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
/* >          RWORK is REAL array, */
/* >                                         dimension (LRWORK) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of the array RWORK. */
/* >          If N <= 1,                LRWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LRWORK must be at least N. */
/* >          If JOBZ  = 'V' and N > 1, LRWORK must be at least */
/* >                         1 + 5*N + 2*N**2. */
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
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          If N <= 1,                LIWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
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
/* >          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed */
/* >                to converge; i off-diagonal elements of an intermediate */
/* >                tridiagonal form did not converge to zero; */
/* >                if INFO = i and JOBZ = 'V', then the algorithm failed */
/* >                to compute an eigenvalue while working on the submatrix */
/* >                lying in rows and columns INFO/(N+1) through */
/* >                mod(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexHEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  Modified description of INFO. Sven, 16 Feb 05. */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cheevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax;
    static integer lopt;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lwmin, liopt;
    static logical lower;
    static integer llrwk, lropt;
    static logical wantz;
    static integer indwk2, llwrk2;
    extern doublereal clanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), cstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int chetrd_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen), clacpy_(char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau, indrwk, indwrk, liwmin;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer lrwmin;
    extern /* Subroutine */ int cunmtr_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen, ftnlen);
    static integer llwork;
    static doublereal smlnum;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 256 "cheevd.f"
    /* Parameter adjustments */
#line 256 "cheevd.f"
    a_dim1 = *lda;
#line 256 "cheevd.f"
    a_offset = 1 + a_dim1;
#line 256 "cheevd.f"
    a -= a_offset;
#line 256 "cheevd.f"
    --w;
#line 256 "cheevd.f"
    --work;
#line 256 "cheevd.f"
    --rwork;
#line 256 "cheevd.f"
    --iwork;
#line 256 "cheevd.f"

#line 256 "cheevd.f"
    /* Function Body */
#line 256 "cheevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 257 "cheevd.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 258 "cheevd.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 260 "cheevd.f"
    *info = 0;
#line 261 "cheevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 262 "cheevd.f"
	*info = -1;
#line 263 "cheevd.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 264 "cheevd.f"
	*info = -2;
#line 265 "cheevd.f"
    } else if (*n < 0) {
#line 266 "cheevd.f"
	*info = -3;
#line 267 "cheevd.f"
    } else if (*lda < max(1,*n)) {
#line 268 "cheevd.f"
	*info = -5;
#line 269 "cheevd.f"
    }

#line 271 "cheevd.f"
    if (*info == 0) {
#line 272 "cheevd.f"
	if (*n <= 1) {
#line 273 "cheevd.f"
	    lwmin = 1;
#line 274 "cheevd.f"
	    lrwmin = 1;
#line 275 "cheevd.f"
	    liwmin = 1;
#line 276 "cheevd.f"
	    lopt = lwmin;
#line 277 "cheevd.f"
	    lropt = lrwmin;
#line 278 "cheevd.f"
	    liopt = liwmin;
#line 279 "cheevd.f"
	} else {
#line 280 "cheevd.f"
	    if (wantz) {
#line 281 "cheevd.f"
		lwmin = (*n << 1) + *n * *n;
/* Computing 2nd power */
#line 282 "cheevd.f"
		i__1 = *n;
#line 282 "cheevd.f"
		lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 283 "cheevd.f"
		liwmin = *n * 5 + 3;
#line 284 "cheevd.f"
	    } else {
#line 285 "cheevd.f"
		lwmin = *n + 1;
#line 286 "cheevd.f"
		lrwmin = *n;
#line 287 "cheevd.f"
		liwmin = 1;
#line 288 "cheevd.f"
	    }
/* Computing MAX */
#line 289 "cheevd.f"
	    i__1 = lwmin, i__2 = *n + ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1,
		     &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 289 "cheevd.f"
	    lopt = max(i__1,i__2);
#line 291 "cheevd.f"
	    lropt = lrwmin;
#line 292 "cheevd.f"
	    liopt = liwmin;
#line 293 "cheevd.f"
	}
#line 294 "cheevd.f"
	work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 295 "cheevd.f"
	rwork[1] = (doublereal) lropt;
#line 296 "cheevd.f"
	iwork[1] = liopt;

#line 298 "cheevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 299 "cheevd.f"
	    *info = -8;
#line 300 "cheevd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 301 "cheevd.f"
	    *info = -10;
#line 302 "cheevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 303 "cheevd.f"
	    *info = -12;
#line 304 "cheevd.f"
	}
#line 305 "cheevd.f"
    }

#line 307 "cheevd.f"
    if (*info != 0) {
#line 308 "cheevd.f"
	i__1 = -(*info);
#line 308 "cheevd.f"
	xerbla_("CHEEVD", &i__1, (ftnlen)6);
#line 309 "cheevd.f"
	return 0;
#line 310 "cheevd.f"
    } else if (lquery) {
#line 311 "cheevd.f"
	return 0;
#line 312 "cheevd.f"
    }

/*     Quick return if possible */

#line 316 "cheevd.f"
    if (*n == 0) {
#line 316 "cheevd.f"
	return 0;
#line 316 "cheevd.f"
    }

#line 319 "cheevd.f"
    if (*n == 1) {
#line 320 "cheevd.f"
	i__1 = a_dim1 + 1;
#line 320 "cheevd.f"
	w[1] = a[i__1].r;
#line 321 "cheevd.f"
	if (wantz) {
#line 321 "cheevd.f"
	    i__1 = a_dim1 + 1;
#line 321 "cheevd.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 321 "cheevd.f"
	}
#line 323 "cheevd.f"
	return 0;
#line 324 "cheevd.f"
    }

/*     Get machine constants. */

#line 328 "cheevd.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 329 "cheevd.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 330 "cheevd.f"
    smlnum = safmin / eps;
#line 331 "cheevd.f"
    bignum = 1. / smlnum;
#line 332 "cheevd.f"
    rmin = sqrt(smlnum);
#line 333 "cheevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 337 "cheevd.f"
    anrm = clanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 338 "cheevd.f"
    iscale = 0;
#line 339 "cheevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 340 "cheevd.f"
	iscale = 1;
#line 341 "cheevd.f"
	sigma = rmin / anrm;
#line 342 "cheevd.f"
    } else if (anrm > rmax) {
#line 343 "cheevd.f"
	iscale = 1;
#line 344 "cheevd.f"
	sigma = rmax / anrm;
#line 345 "cheevd.f"
    }
#line 346 "cheevd.f"
    if (iscale == 1) {
#line 346 "cheevd.f"
	clascl_(uplo, &c__0, &c__0, &c_b18, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 346 "cheevd.f"
    }

/*     Call CHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 351 "cheevd.f"
    inde = 1;
#line 352 "cheevd.f"
    indtau = 1;
#line 353 "cheevd.f"
    indwrk = indtau + *n;
#line 354 "cheevd.f"
    indrwk = inde + *n;
#line 355 "cheevd.f"
    indwk2 = indwrk + *n * *n;
#line 356 "cheevd.f"
    llwork = *lwork - indwrk + 1;
#line 357 "cheevd.f"
    llwrk2 = *lwork - indwk2 + 1;
#line 358 "cheevd.f"
    llrwk = *lrwork - indrwk + 1;
#line 359 "cheevd.f"
    chetrd_(uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     CSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call CUNMTR to multiply it to the */
/*     Householder transformations represented as Householder vectors in */
/*     A. */

#line 368 "cheevd.f"
    if (! wantz) {
#line 369 "cheevd.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 370 "cheevd.f"
    } else {
#line 371 "cheevd.f"
	cstedc_("I", n, &w[1], &rwork[inde], &work[indwrk], n, &work[indwk2], 
		&llwrk2, &rwork[indrwk], &llrwk, &iwork[1], liwork, info, (
		ftnlen)1);
#line 374 "cheevd.f"
	cunmtr_("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[
		indwrk], n, &work[indwk2], &llwrk2, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 376 "cheevd.f"
	clacpy_("A", n, n, &work[indwrk], n, &a[a_offset], lda, (ftnlen)1);
#line 377 "cheevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 381 "cheevd.f"
    if (iscale == 1) {
#line 382 "cheevd.f"
	if (*info == 0) {
#line 383 "cheevd.f"
	    imax = *n;
#line 384 "cheevd.f"
	} else {
#line 385 "cheevd.f"
	    imax = *info - 1;
#line 386 "cheevd.f"
	}
#line 387 "cheevd.f"
	d__1 = 1. / sigma;
#line 387 "cheevd.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 388 "cheevd.f"
    }

#line 390 "cheevd.f"
    work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 391 "cheevd.f"
    rwork[1] = (doublereal) lropt;
#line 392 "cheevd.f"
    iwork[1] = liopt;

#line 394 "cheevd.f"
    return 0;

/*     End of CHEEVD */

} /* cheevd_ */

