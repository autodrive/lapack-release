#line 1 "dsyevd.f"
/* dsyevd.f -- translated by f2c (version 20100827).
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

#line 1 "dsyevd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b17 = 1.;

/* > \brief <b> DSYEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEVD computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric matrix A. If eigenvectors are desired, it uses a */
/* > divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > */
/* > Because of large use of BLAS of level 3, DSYEVD needs N**2 more */
/* > workspace than DSYEVX. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, */
/* >                                         dimension (LWORK) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least */
/* >                                                1 + 6*N + 2*N**2. */
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

/* > \date December 2016 */

/* > \ingroup doubleSYeigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* >  Modified by Francoise Tisseur, University of Tennessee \n */
/* >  Modified description of INFO. Sven, 16 Feb 05. \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
	a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *iwork, integer *liwork, integer *info, ftnlen jobz_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm, rmin, rmax;
    static integer lopt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, lwmin, liopt;
    static logical lower, wantz;
    static integer indwk2, llwrk2;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int dormtr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), dsytrd_(char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen);
    static integer llwork;
    static doublereal smlnum;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 233 "dsyevd.f"
    /* Parameter adjustments */
#line 233 "dsyevd.f"
    a_dim1 = *lda;
#line 233 "dsyevd.f"
    a_offset = 1 + a_dim1;
#line 233 "dsyevd.f"
    a -= a_offset;
#line 233 "dsyevd.f"
    --w;
#line 233 "dsyevd.f"
    --work;
#line 233 "dsyevd.f"
    --iwork;
#line 233 "dsyevd.f"

#line 233 "dsyevd.f"
    /* Function Body */
#line 233 "dsyevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 234 "dsyevd.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 235 "dsyevd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 237 "dsyevd.f"
    *info = 0;
#line 238 "dsyevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 239 "dsyevd.f"
	*info = -1;
#line 240 "dsyevd.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 241 "dsyevd.f"
	*info = -2;
#line 242 "dsyevd.f"
    } else if (*n < 0) {
#line 243 "dsyevd.f"
	*info = -3;
#line 244 "dsyevd.f"
    } else if (*lda < max(1,*n)) {
#line 245 "dsyevd.f"
	*info = -5;
#line 246 "dsyevd.f"
    }

#line 248 "dsyevd.f"
    if (*info == 0) {
#line 249 "dsyevd.f"
	if (*n <= 1) {
#line 250 "dsyevd.f"
	    liwmin = 1;
#line 251 "dsyevd.f"
	    lwmin = 1;
#line 252 "dsyevd.f"
	    lopt = lwmin;
#line 253 "dsyevd.f"
	    liopt = liwmin;
#line 254 "dsyevd.f"
	} else {
#line 255 "dsyevd.f"
	    if (wantz) {
#line 256 "dsyevd.f"
		liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 257 "dsyevd.f"
		i__1 = *n;
#line 257 "dsyevd.f"
		lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
#line 258 "dsyevd.f"
	    } else {
#line 259 "dsyevd.f"
		liwmin = 1;
#line 260 "dsyevd.f"
		lwmin = (*n << 1) + 1;
#line 261 "dsyevd.f"
	    }
/* Computing MAX */
#line 262 "dsyevd.f"
	    i__1 = lwmin, i__2 = (*n << 1) + ilaenv_(&c__1, "DSYTRD", uplo, n,
		     &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 262 "dsyevd.f"
	    lopt = max(i__1,i__2);
#line 264 "dsyevd.f"
	    liopt = liwmin;
#line 265 "dsyevd.f"
	}
#line 266 "dsyevd.f"
	work[1] = (doublereal) lopt;
#line 267 "dsyevd.f"
	iwork[1] = liopt;

#line 269 "dsyevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 270 "dsyevd.f"
	    *info = -8;
#line 271 "dsyevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 272 "dsyevd.f"
	    *info = -10;
#line 273 "dsyevd.f"
	}
#line 274 "dsyevd.f"
    }

#line 276 "dsyevd.f"
    if (*info != 0) {
#line 277 "dsyevd.f"
	i__1 = -(*info);
#line 277 "dsyevd.f"
	xerbla_("DSYEVD", &i__1, (ftnlen)6);
#line 278 "dsyevd.f"
	return 0;
#line 279 "dsyevd.f"
    } else if (lquery) {
#line 280 "dsyevd.f"
	return 0;
#line 281 "dsyevd.f"
    }

/*     Quick return if possible */

#line 285 "dsyevd.f"
    if (*n == 0) {
#line 285 "dsyevd.f"
	return 0;
#line 285 "dsyevd.f"
    }

#line 288 "dsyevd.f"
    if (*n == 1) {
#line 289 "dsyevd.f"
	w[1] = a[a_dim1 + 1];
#line 290 "dsyevd.f"
	if (wantz) {
#line 290 "dsyevd.f"
	    a[a_dim1 + 1] = 1.;
#line 290 "dsyevd.f"
	}
#line 292 "dsyevd.f"
	return 0;
#line 293 "dsyevd.f"
    }

/*     Get machine constants. */

#line 297 "dsyevd.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 298 "dsyevd.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 299 "dsyevd.f"
    smlnum = safmin / eps;
#line 300 "dsyevd.f"
    bignum = 1. / smlnum;
#line 301 "dsyevd.f"
    rmin = sqrt(smlnum);
#line 302 "dsyevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 306 "dsyevd.f"
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 307 "dsyevd.f"
    iscale = 0;
#line 308 "dsyevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 309 "dsyevd.f"
	iscale = 1;
#line 310 "dsyevd.f"
	sigma = rmin / anrm;
#line 311 "dsyevd.f"
    } else if (anrm > rmax) {
#line 312 "dsyevd.f"
	iscale = 1;
#line 313 "dsyevd.f"
	sigma = rmax / anrm;
#line 314 "dsyevd.f"
    }
#line 315 "dsyevd.f"
    if (iscale == 1) {
#line 315 "dsyevd.f"
	dlascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 315 "dsyevd.f"
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

#line 320 "dsyevd.f"
    inde = 1;
#line 321 "dsyevd.f"
    indtau = inde + *n;
#line 322 "dsyevd.f"
    indwrk = indtau + *n;
#line 323 "dsyevd.f"
    llwork = *lwork - indwrk + 1;
#line 324 "dsyevd.f"
    indwk2 = indwrk + *n * *n;
#line 325 "dsyevd.f"
    llwrk2 = *lwork - indwk2 + 1;

#line 327 "dsyevd.f"
    dsytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call DORMTR to multiply it by the */
/*     Householder transformations stored in A. */

#line 335 "dsyevd.f"
    if (! wantz) {
#line 336 "dsyevd.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 337 "dsyevd.f"
    } else {
#line 338 "dsyevd.f"
	dstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &
		llwrk2, &iwork[1], liwork, info, (ftnlen)1);
#line 340 "dsyevd.f"
	dormtr_("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[
		indwrk], n, &work[indwk2], &llwrk2, &iinfo, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 342 "dsyevd.f"
	dlacpy_("A", n, n, &work[indwrk], n, &a[a_offset], lda, (ftnlen)1);
#line 343 "dsyevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 347 "dsyevd.f"
    if (iscale == 1) {
#line 347 "dsyevd.f"
	d__1 = 1. / sigma;
#line 347 "dsyevd.f"
	dscal_(n, &d__1, &w[1], &c__1);
#line 347 "dsyevd.f"
    }

#line 350 "dsyevd.f"
    work[1] = (doublereal) lopt;
#line 351 "dsyevd.f"
    iwork[1] = liopt;

#line 353 "dsyevd.f"
    return 0;

/*     End of DSYEVD */

} /* dsyevd_ */

