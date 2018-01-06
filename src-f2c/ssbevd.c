#line 1 "ssbevd.f"
/* ssbevd.f -- translated by f2c (version 20100827).
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

#line 1 "ssbevd.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static doublereal c_b18 = 0.;
static integer c__1 = 1;

/* > \brief <b> SSBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSBEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/*                          LWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A. If eigenvectors are desired, it uses */
/* > a divide and conquer algorithm. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, AB is overwritten by values generated during the */
/* >          reduction to tridiagonal form.  If UPLO = 'U', the first */
/* >          superdiagonal and the diagonal of the tridiagonal matrix T */
/* >          are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* >          the diagonal and first subdiagonal of T are returned in the */
/* >          first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD + 1. */
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
/* >          Z is REAL array, dimension (LDZ, N) */
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
/* >          WORK is REAL array, */
/* >                                         dimension (LWORK) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          IF N <= 1,                LWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N. */
/* >          If JOBZ  = 'V' and N > 2, LWORK must be at least */
/* >                         ( 1 + 5*N + 2*N**2 ). */
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
/* >          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N. */
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

/* > \date December 2016 */

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int ssbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps;
    static integer inde;
    static doublereal anrm, rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical lower, wantz;
    static integer indwk2, llwrk2, iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal slansb_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), sstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen), slacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int ssbtrd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen, ftnlen), ssterf_(
	    integer *, doublereal *, doublereal *, integer *);
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

#line 239 "ssbevd.f"
    /* Parameter adjustments */
#line 239 "ssbevd.f"
    ab_dim1 = *ldab;
#line 239 "ssbevd.f"
    ab_offset = 1 + ab_dim1;
#line 239 "ssbevd.f"
    ab -= ab_offset;
#line 239 "ssbevd.f"
    --w;
#line 239 "ssbevd.f"
    z_dim1 = *ldz;
#line 239 "ssbevd.f"
    z_offset = 1 + z_dim1;
#line 239 "ssbevd.f"
    z__ -= z_offset;
#line 239 "ssbevd.f"
    --work;
#line 239 "ssbevd.f"
    --iwork;
#line 239 "ssbevd.f"

#line 239 "ssbevd.f"
    /* Function Body */
#line 239 "ssbevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 240 "ssbevd.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 241 "ssbevd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 243 "ssbevd.f"
    *info = 0;
#line 244 "ssbevd.f"
    if (*n <= 1) {
#line 245 "ssbevd.f"
	liwmin = 1;
#line 246 "ssbevd.f"
	lwmin = 1;
#line 247 "ssbevd.f"
    } else {
#line 248 "ssbevd.f"
	if (wantz) {
#line 249 "ssbevd.f"
	    liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 250 "ssbevd.f"
	    i__1 = *n;
#line 250 "ssbevd.f"
	    lwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 251 "ssbevd.f"
	} else {
#line 252 "ssbevd.f"
	    liwmin = 1;
#line 253 "ssbevd.f"
	    lwmin = *n << 1;
#line 254 "ssbevd.f"
	}
#line 255 "ssbevd.f"
    }
#line 256 "ssbevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 257 "ssbevd.f"
	*info = -1;
#line 258 "ssbevd.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 259 "ssbevd.f"
	*info = -2;
#line 260 "ssbevd.f"
    } else if (*n < 0) {
#line 261 "ssbevd.f"
	*info = -3;
#line 262 "ssbevd.f"
    } else if (*kd < 0) {
#line 263 "ssbevd.f"
	*info = -4;
#line 264 "ssbevd.f"
    } else if (*ldab < *kd + 1) {
#line 265 "ssbevd.f"
	*info = -6;
#line 266 "ssbevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 267 "ssbevd.f"
	*info = -9;
#line 268 "ssbevd.f"
    }

#line 270 "ssbevd.f"
    if (*info == 0) {
#line 271 "ssbevd.f"
	work[1] = (doublereal) lwmin;
#line 272 "ssbevd.f"
	iwork[1] = liwmin;

#line 274 "ssbevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 275 "ssbevd.f"
	    *info = -11;
#line 276 "ssbevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 277 "ssbevd.f"
	    *info = -13;
#line 278 "ssbevd.f"
	}
#line 279 "ssbevd.f"
    }

#line 281 "ssbevd.f"
    if (*info != 0) {
#line 282 "ssbevd.f"
	i__1 = -(*info);
#line 282 "ssbevd.f"
	xerbla_("SSBEVD", &i__1, (ftnlen)6);
#line 283 "ssbevd.f"
	return 0;
#line 284 "ssbevd.f"
    } else if (lquery) {
#line 285 "ssbevd.f"
	return 0;
#line 286 "ssbevd.f"
    }

/*     Quick return if possible */

#line 290 "ssbevd.f"
    if (*n == 0) {
#line 290 "ssbevd.f"
	return 0;
#line 290 "ssbevd.f"
    }

#line 293 "ssbevd.f"
    if (*n == 1) {
#line 294 "ssbevd.f"
	w[1] = ab[ab_dim1 + 1];
#line 295 "ssbevd.f"
	if (wantz) {
#line 295 "ssbevd.f"
	    z__[z_dim1 + 1] = 1.;
#line 295 "ssbevd.f"
	}
#line 297 "ssbevd.f"
	return 0;
#line 298 "ssbevd.f"
    }

/*     Get machine constants. */

#line 302 "ssbevd.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 303 "ssbevd.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 304 "ssbevd.f"
    smlnum = safmin / eps;
#line 305 "ssbevd.f"
    bignum = 1. / smlnum;
#line 306 "ssbevd.f"
    rmin = sqrt(smlnum);
#line 307 "ssbevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 311 "ssbevd.f"
    anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);
#line 312 "ssbevd.f"
    iscale = 0;
#line 313 "ssbevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 314 "ssbevd.f"
	iscale = 1;
#line 315 "ssbevd.f"
	sigma = rmin / anrm;
#line 316 "ssbevd.f"
    } else if (anrm > rmax) {
#line 317 "ssbevd.f"
	iscale = 1;
#line 318 "ssbevd.f"
	sigma = rmax / anrm;
#line 319 "ssbevd.f"
    }
#line 320 "ssbevd.f"
    if (iscale == 1) {
#line 321 "ssbevd.f"
	if (lower) {
#line 322 "ssbevd.f"
	    slascl_("B", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 323 "ssbevd.f"
	} else {
#line 324 "ssbevd.f"
	    slascl_("Q", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, 
		    info, (ftnlen)1);
#line 325 "ssbevd.f"
	}
#line 326 "ssbevd.f"
    }

/*     Call SSBTRD to reduce symmetric band matrix to tridiagonal form. */

#line 330 "ssbevd.f"
    inde = 1;
#line 331 "ssbevd.f"
    indwrk = inde + *n;
#line 332 "ssbevd.f"
    indwk2 = indwrk + *n * *n;
#line 333 "ssbevd.f"
    llwrk2 = *lwork - indwk2 + 1;
#line 334 "ssbevd.f"
    ssbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[
	    z_offset], ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, call SSTEDC. */

#line 339 "ssbevd.f"
    if (! wantz) {
#line 340 "ssbevd.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 341 "ssbevd.f"
    } else {
#line 342 "ssbevd.f"
	sstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &
		llwrk2, &iwork[1], liwork, info, (ftnlen)1);
#line 344 "ssbevd.f"
	sgemm_("N", "N", n, n, n, &c_b11, &z__[z_offset], ldz, &work[indwrk], 
		n, &c_b18, &work[indwk2], n, (ftnlen)1, (ftnlen)1);
#line 346 "ssbevd.f"
	slacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz, (ftnlen)1);
#line 347 "ssbevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 351 "ssbevd.f"
    if (iscale == 1) {
#line 351 "ssbevd.f"
	d__1 = 1. / sigma;
#line 351 "ssbevd.f"
	sscal_(n, &d__1, &w[1], &c__1);
#line 351 "ssbevd.f"
    }

#line 354 "ssbevd.f"
    work[1] = (doublereal) lwmin;
#line 355 "ssbevd.f"
    iwork[1] = liwmin;
#line 356 "ssbevd.f"
    return 0;

/*     End of SSBEVD */

} /* ssbevd_ */

