#line 1 "ssyev.f"
/* ssyev.f -- translated by f2c (version 20100827).
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

#line 1 "ssyev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b17 = 1.;

/* > \brief <b> SSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEV computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric matrix A. */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
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
/* >          The length of the array WORK.  LWORK >= max(1,3*N-1). */
/* >          For optimal efficiency, LWORK >= (NB+2)*N, */
/* >          where NB is the blocksize for SSYTRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \ingroup realSYeigen */

/*  ===================================================================== */
/* Subroutine */ int ssyev_(char *jobz, char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer nb;
    static doublereal eps;
    static integer inde;
    static doublereal anrm;
    static integer imax;
    static doublereal rmin, rmax, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lower, wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer indtau, indwrk;
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    extern doublereal slansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer llwork;
    static doublereal smlnum;
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sorgtr_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), ssteqr_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ssytrd_(char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);


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

#line 178 "ssyev.f"
    /* Parameter adjustments */
#line 178 "ssyev.f"
    a_dim1 = *lda;
#line 178 "ssyev.f"
    a_offset = 1 + a_dim1;
#line 178 "ssyev.f"
    a -= a_offset;
#line 178 "ssyev.f"
    --w;
#line 178 "ssyev.f"
    --work;
#line 178 "ssyev.f"

#line 178 "ssyev.f"
    /* Function Body */
#line 178 "ssyev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 179 "ssyev.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 180 "ssyev.f"
    lquery = *lwork == -1;

#line 182 "ssyev.f"
    *info = 0;
#line 183 "ssyev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 184 "ssyev.f"
	*info = -1;
#line 185 "ssyev.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 186 "ssyev.f"
	*info = -2;
#line 187 "ssyev.f"
    } else if (*n < 0) {
#line 188 "ssyev.f"
	*info = -3;
#line 189 "ssyev.f"
    } else if (*lda < max(1,*n)) {
#line 190 "ssyev.f"
	*info = -5;
#line 191 "ssyev.f"
    }

#line 193 "ssyev.f"
    if (*info == 0) {
#line 194 "ssyev.f"
	nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 195 "ssyev.f"
	i__1 = 1, i__2 = (nb + 2) * *n;
#line 195 "ssyev.f"
	lwkopt = max(i__1,i__2);
#line 196 "ssyev.f"
	work[1] = (doublereal) lwkopt;

/* Computing MAX */
#line 198 "ssyev.f"
	i__1 = 1, i__2 = *n * 3 - 1;
#line 198 "ssyev.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 198 "ssyev.f"
	    *info = -8;
#line 198 "ssyev.f"
	}
#line 200 "ssyev.f"
    }

#line 202 "ssyev.f"
    if (*info != 0) {
#line 203 "ssyev.f"
	i__1 = -(*info);
#line 203 "ssyev.f"
	xerbla_("SSYEV ", &i__1, (ftnlen)6);
#line 204 "ssyev.f"
	return 0;
#line 205 "ssyev.f"
    } else if (lquery) {
#line 206 "ssyev.f"
	return 0;
#line 207 "ssyev.f"
    }

/*     Quick return if possible */

#line 211 "ssyev.f"
    if (*n == 0) {
#line 212 "ssyev.f"
	return 0;
#line 213 "ssyev.f"
    }

#line 215 "ssyev.f"
    if (*n == 1) {
#line 216 "ssyev.f"
	w[1] = a[a_dim1 + 1];
#line 217 "ssyev.f"
	work[1] = 2.;
#line 218 "ssyev.f"
	if (wantz) {
#line 218 "ssyev.f"
	    a[a_dim1 + 1] = 1.;
#line 218 "ssyev.f"
	}
#line 220 "ssyev.f"
	return 0;
#line 221 "ssyev.f"
    }

/*     Get machine constants. */

#line 225 "ssyev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 226 "ssyev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 227 "ssyev.f"
    smlnum = safmin / eps;
#line 228 "ssyev.f"
    bignum = 1. / smlnum;
#line 229 "ssyev.f"
    rmin = sqrt(smlnum);
#line 230 "ssyev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 234 "ssyev.f"
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);
#line 235 "ssyev.f"
    iscale = 0;
#line 236 "ssyev.f"
    if (anrm > 0. && anrm < rmin) {
#line 237 "ssyev.f"
	iscale = 1;
#line 238 "ssyev.f"
	sigma = rmin / anrm;
#line 239 "ssyev.f"
    } else if (anrm > rmax) {
#line 240 "ssyev.f"
	iscale = 1;
#line 241 "ssyev.f"
	sigma = rmax / anrm;
#line 242 "ssyev.f"
    }
#line 243 "ssyev.f"
    if (iscale == 1) {
#line 243 "ssyev.f"
	slascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 243 "ssyev.f"
    }

/*     Call SSYTRD to reduce symmetric matrix to tridiagonal form. */

#line 248 "ssyev.f"
    inde = 1;
#line 249 "ssyev.f"
    indtau = inde + *n;
#line 250 "ssyev.f"
    indwrk = indtau + *n;
#line 251 "ssyev.f"
    llwork = *lwork - indwrk + 1;
#line 252 "ssyev.f"
    ssytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     SORGTR to generate the orthogonal matrix, then call SSTEQR. */

#line 258 "ssyev.f"
    if (! wantz) {
#line 259 "ssyev.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 260 "ssyev.f"
    } else {
#line 261 "ssyev.f"
	sorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo, (ftnlen)1);
#line 263 "ssyev.f"
	ssteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
		 info, (ftnlen)1);
#line 265 "ssyev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 269 "ssyev.f"
    if (iscale == 1) {
#line 270 "ssyev.f"
	if (*info == 0) {
#line 271 "ssyev.f"
	    imax = *n;
#line 272 "ssyev.f"
	} else {
#line 273 "ssyev.f"
	    imax = *info - 1;
#line 274 "ssyev.f"
	}
#line 275 "ssyev.f"
	d__1 = 1. / sigma;
#line 275 "ssyev.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 276 "ssyev.f"
    }

/*     Set WORK(1) to optimal workspace size. */

#line 280 "ssyev.f"
    work[1] = (doublereal) lwkopt;

#line 282 "ssyev.f"
    return 0;

/*     End of SSYEV */

} /* ssyev_ */

