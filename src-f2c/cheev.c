#line 1 "cheev.f"
/* cheev.f -- translated by f2c (version 20100827).
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

#line 1 "cheev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b18 = 1.;

/* > \brief <b> CHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEV computes all eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix A. */
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
/* >          The length of the array WORK.  LWORK >= max(1,2*N-1). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the blocksize for CHETRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (max(1, 3*N-2)) */
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

/* > \ingroup complexHEeigen */

/*  ===================================================================== */
/* Subroutine */ int cheev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
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
    extern doublereal clanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int chetrd_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau, indwrk;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), cungtr_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen), 
	    ssterf_(integer *, doublereal *, doublereal *, integer *);
    static integer llwork;
    static doublereal smlnum;
    static integer lwkopt;
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

#line 189 "cheev.f"
    /* Parameter adjustments */
#line 189 "cheev.f"
    a_dim1 = *lda;
#line 189 "cheev.f"
    a_offset = 1 + a_dim1;
#line 189 "cheev.f"
    a -= a_offset;
#line 189 "cheev.f"
    --w;
#line 189 "cheev.f"
    --work;
#line 189 "cheev.f"
    --rwork;
#line 189 "cheev.f"

#line 189 "cheev.f"
    /* Function Body */
#line 189 "cheev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 190 "cheev.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 191 "cheev.f"
    lquery = *lwork == -1;

#line 193 "cheev.f"
    *info = 0;
#line 194 "cheev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 195 "cheev.f"
	*info = -1;
#line 196 "cheev.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 197 "cheev.f"
	*info = -2;
#line 198 "cheev.f"
    } else if (*n < 0) {
#line 199 "cheev.f"
	*info = -3;
#line 200 "cheev.f"
    } else if (*lda < max(1,*n)) {
#line 201 "cheev.f"
	*info = -5;
#line 202 "cheev.f"
    }

#line 204 "cheev.f"
    if (*info == 0) {
#line 205 "cheev.f"
	nb = ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 206 "cheev.f"
	i__1 = 1, i__2 = (nb + 1) * *n;
#line 206 "cheev.f"
	lwkopt = max(i__1,i__2);
#line 207 "cheev.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/* Computing MAX */
#line 209 "cheev.f"
	i__1 = 1, i__2 = (*n << 1) - 1;
#line 209 "cheev.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 209 "cheev.f"
	    *info = -8;
#line 209 "cheev.f"
	}
#line 211 "cheev.f"
    }

#line 213 "cheev.f"
    if (*info != 0) {
#line 214 "cheev.f"
	i__1 = -(*info);
#line 214 "cheev.f"
	xerbla_("CHEEV ", &i__1, (ftnlen)6);
#line 215 "cheev.f"
	return 0;
#line 216 "cheev.f"
    } else if (lquery) {
#line 217 "cheev.f"
	return 0;
#line 218 "cheev.f"
    }

/*     Quick return if possible */

#line 222 "cheev.f"
    if (*n == 0) {
#line 223 "cheev.f"
	return 0;
#line 224 "cheev.f"
    }

#line 226 "cheev.f"
    if (*n == 1) {
#line 227 "cheev.f"
	i__1 = a_dim1 + 1;
#line 227 "cheev.f"
	w[1] = a[i__1].r;
#line 228 "cheev.f"
	work[1].r = 1., work[1].i = 0.;
#line 229 "cheev.f"
	if (wantz) {
#line 229 "cheev.f"
	    i__1 = a_dim1 + 1;
#line 229 "cheev.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 229 "cheev.f"
	}
#line 231 "cheev.f"
	return 0;
#line 232 "cheev.f"
    }

/*     Get machine constants. */

#line 236 "cheev.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 237 "cheev.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 238 "cheev.f"
    smlnum = safmin / eps;
#line 239 "cheev.f"
    bignum = 1. / smlnum;
#line 240 "cheev.f"
    rmin = sqrt(smlnum);
#line 241 "cheev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 245 "cheev.f"
    anrm = clanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 246 "cheev.f"
    iscale = 0;
#line 247 "cheev.f"
    if (anrm > 0. && anrm < rmin) {
#line 248 "cheev.f"
	iscale = 1;
#line 249 "cheev.f"
	sigma = rmin / anrm;
#line 250 "cheev.f"
    } else if (anrm > rmax) {
#line 251 "cheev.f"
	iscale = 1;
#line 252 "cheev.f"
	sigma = rmax / anrm;
#line 253 "cheev.f"
    }
#line 254 "cheev.f"
    if (iscale == 1) {
#line 254 "cheev.f"
	clascl_(uplo, &c__0, &c__0, &c_b18, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 254 "cheev.f"
    }

/*     Call CHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 259 "cheev.f"
    inde = 1;
#line 260 "cheev.f"
    indtau = 1;
#line 261 "cheev.f"
    indwrk = indtau + *n;
#line 262 "cheev.f"
    llwork = *lwork - indwrk + 1;
#line 263 "cheev.f"
    chetrd_(uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     CUNGTR to generate the unitary matrix, then call CSTEQR. */

#line 269 "cheev.f"
    if (! wantz) {
#line 270 "cheev.f"
	ssterf_(n, &w[1], &rwork[inde], info);
#line 271 "cheev.f"
    } else {
#line 272 "cheev.f"
	cungtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo, (ftnlen)1);
#line 274 "cheev.f"
	indwrk = inde + *n;
#line 275 "cheev.f"
	csteqr_(jobz, n, &w[1], &rwork[inde], &a[a_offset], lda, &rwork[
		indwrk], info, (ftnlen)1);
#line 277 "cheev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 281 "cheev.f"
    if (iscale == 1) {
#line 282 "cheev.f"
	if (*info == 0) {
#line 283 "cheev.f"
	    imax = *n;
#line 284 "cheev.f"
	} else {
#line 285 "cheev.f"
	    imax = *info - 1;
#line 286 "cheev.f"
	}
#line 287 "cheev.f"
	d__1 = 1. / sigma;
#line 287 "cheev.f"
	sscal_(&imax, &d__1, &w[1], &c__1);
#line 288 "cheev.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 292 "cheev.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 294 "cheev.f"
    return 0;

/*     End of CHEEV */

} /* cheev_ */

