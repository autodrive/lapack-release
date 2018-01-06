#line 1 "chpevd.f"
/* chpevd.f -- translated by f2c (version 20100827).
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

#line 1 "chpevd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CHPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, */
/*                          RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian matrix A in packed storage.  If eigenvectors are */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, AP is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the diagonal */
/* >          and first superdiagonal of the tridiagonal matrix T overwrite */
/* >          the corresponding elements of A, and if UPLO = 'L', the */
/* >          diagonal and first subdiagonal of T overwrite the */
/* >          corresponding elements of A. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* >          eigenvectors of the matrix A, with the i-th column of Z */
/* >          holding the eigenvector associated with W(i). */
/* >          If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of array WORK. */
/* >          If N <= 1,               LWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be at least N. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the required sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* >          On exit, if INFO = 0, RWORK(1) returns the required LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of array RWORK. */
/* >          If N <= 1,               LRWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LRWORK must be at least N. */
/* >          If JOBZ = 'V' and N > 1, LRWORK must be at least */
/* >                    1 + 5*N + 2*N**2. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the required LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of array IWORK. */
/* >          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the algorithm failed to converge; i */
/* >                off-diagonal elements of an intermediate tridiagonal */
/* >                form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int chpevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info, ftnlen 
	jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lwmin, llrwk, llwrk;
    static logical wantz;
    static integer iscale;
    extern doublereal clanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int cstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int chptrd_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, ftnlen);
    static integer indrwk, indwrk, liwmin;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer lrwmin;
    extern /* Subroutine */ int cupmtr_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);
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

#line 249 "chpevd.f"
    /* Parameter adjustments */
#line 249 "chpevd.f"
    --ap;
#line 249 "chpevd.f"
    --w;
#line 249 "chpevd.f"
    z_dim1 = *ldz;
#line 249 "chpevd.f"
    z_offset = 1 + z_dim1;
#line 249 "chpevd.f"
    z__ -= z_offset;
#line 249 "chpevd.f"
    --work;
#line 249 "chpevd.f"
    --rwork;
#line 249 "chpevd.f"
    --iwork;
#line 249 "chpevd.f"

#line 249 "chpevd.f"
    /* Function Body */
#line 249 "chpevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 250 "chpevd.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 252 "chpevd.f"
    *info = 0;
#line 253 "chpevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 254 "chpevd.f"
	*info = -1;
#line 255 "chpevd.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 257 "chpevd.f"
	*info = -2;
#line 258 "chpevd.f"
    } else if (*n < 0) {
#line 259 "chpevd.f"
	*info = -3;
#line 260 "chpevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 261 "chpevd.f"
	*info = -7;
#line 262 "chpevd.f"
    }

#line 264 "chpevd.f"
    if (*info == 0) {
#line 265 "chpevd.f"
	if (*n <= 1) {
#line 266 "chpevd.f"
	    lwmin = 1;
#line 267 "chpevd.f"
	    liwmin = 1;
#line 268 "chpevd.f"
	    lrwmin = 1;
#line 269 "chpevd.f"
	} else {
#line 270 "chpevd.f"
	    if (wantz) {
#line 271 "chpevd.f"
		lwmin = *n << 1;
/* Computing 2nd power */
#line 272 "chpevd.f"
		i__1 = *n;
#line 272 "chpevd.f"
		lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 273 "chpevd.f"
		liwmin = *n * 5 + 3;
#line 274 "chpevd.f"
	    } else {
#line 275 "chpevd.f"
		lwmin = *n;
#line 276 "chpevd.f"
		lrwmin = *n;
#line 277 "chpevd.f"
		liwmin = 1;
#line 278 "chpevd.f"
	    }
#line 279 "chpevd.f"
	}
#line 280 "chpevd.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 281 "chpevd.f"
	rwork[1] = (doublereal) lrwmin;
#line 282 "chpevd.f"
	iwork[1] = liwmin;

#line 284 "chpevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 285 "chpevd.f"
	    *info = -9;
#line 286 "chpevd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 287 "chpevd.f"
	    *info = -11;
#line 288 "chpevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 289 "chpevd.f"
	    *info = -13;
#line 290 "chpevd.f"
	}
#line 291 "chpevd.f"
    }

#line 293 "chpevd.f"
    if (*info != 0) {
#line 294 "chpevd.f"
	i__1 = -(*info);
#line 294 "chpevd.f"
	xerbla_("CHPEVD", &i__1, (ftnlen)6);
#line 295 "chpevd.f"
	return 0;
#line 296 "chpevd.f"
    } else if (lquery) {
#line 297 "chpevd.f"
	return 0;
#line 298 "chpevd.f"
    }

/*     Quick return if possible */

#line 302 "chpevd.f"
    if (*n == 0) {
#line 302 "chpevd.f"
	return 0;
#line 302 "chpevd.f"
    }

#line 305 "chpevd.f"
    if (*n == 1) {
#line 306 "chpevd.f"
	w[1] = ap[1].r;
#line 307 "chpevd.f"
	if (wantz) {
#line 307 "chpevd.f"
	    i__1 = z_dim1 + 1;
#line 307 "chpevd.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 307 "chpevd.f"
	}
#line 309 "chpevd.f"
	return 0;
#line 310 "chpevd.f"
    }

/*     Get machine constants. */

#line 314 "chpevd.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 315 "chpevd.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 316 "chpevd.f"
    smlnum = safmin / eps;
#line 317 "chpevd.f"
    bignum = 1. / smlnum;
#line 318 "chpevd.f"
    rmin = sqrt(smlnum);
#line 319 "chpevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 323 "chpevd.f"
    anrm = clanhp_("M", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);
#line 324 "chpevd.f"
    iscale = 0;
#line 325 "chpevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 326 "chpevd.f"
	iscale = 1;
#line 327 "chpevd.f"
	sigma = rmin / anrm;
#line 328 "chpevd.f"
    } else if (anrm > rmax) {
#line 329 "chpevd.f"
	iscale = 1;
#line 330 "chpevd.f"
	sigma = rmax / anrm;
#line 331 "chpevd.f"
    }
#line 332 "chpevd.f"
    if (iscale == 1) {
#line 333 "chpevd.f"
	i__1 = *n * (*n + 1) / 2;
#line 333 "chpevd.f"
	csscal_(&i__1, &sigma, &ap[1], &c__1);
#line 334 "chpevd.f"
    }

/*     Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form. */

#line 338 "chpevd.f"
    inde = 1;
#line 339 "chpevd.f"
    indtau = 1;
#line 340 "chpevd.f"
    indrwk = inde + *n;
#line 341 "chpevd.f"
    indwrk = indtau + *n;
#line 342 "chpevd.f"
    llwrk = *lwork - indwrk + 1;
#line 343 "chpevd.f"
    llrwk = *lrwork - indrwk + 1;
#line 344 "chpevd.f"
    chptrd_(uplo, n, &ap[1], &w[1], &rwork[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     CUPGTR to generate the orthogonal matrix, then call CSTEDC. */

#line 350 "chpevd.f"
    if (! wantz) {
#line 351 "chpevd.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 352 "chpevd.f"
    } else {
#line 353 "chpevd.f"
	cstedc_("I", n, &w[1], &rwork[inde], &z__[z_offset], ldz, &work[
		indwrk], &llwrk, &rwork[indrwk], &llrwk, &iwork[1], liwork, 
		info, (ftnlen)1);
#line 356 "chpevd.f"
	cupmtr_("L", uplo, "N", n, n, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 358 "chpevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 362 "chpevd.f"
    if (iscale == 1) {
#line 363 "chpevd.f"
	if (*info == 0) {
#line 364 "chpevd.f"
	    imax = *n;
#line 365 "chpevd.f"
	} else {
#line 366 "chpevd.f"
	    imax = *info - 1;
#line 367 "chpevd.f"
	}
#line 368 "chpevd.f"
	d__1 = 1. / sigma;
#line 368 "chpevd.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 369 "chpevd.f"
    }

#line 371 "chpevd.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 372 "chpevd.f"
    rwork[1] = (doublereal) lrwmin;
#line 373 "chpevd.f"
    iwork[1] = liwmin;
#line 374 "chpevd.f"
    return 0;

/*     End of CHPEVD */

} /* chpevd_ */

