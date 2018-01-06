#line 1 "zheev.f"
/* zheev.f -- translated by f2c (version 20100827).
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

#line 1 "zheev.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b18 = 1.;

/* > \brief <b> ZHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEV computes all eigenvalues and, optionally, eigenvectors of a */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,2*N-1). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the blocksize for ZHETRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2)) */
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

/* > \ingroup complex16HEeigen */

/*  ===================================================================== */
/* Subroutine */ int zheev_(char *jobz, char *uplo, integer *n, doublecomplex 
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
    static doublereal rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), zlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen);
    static integer indwrk;
    extern /* Subroutine */ int zhetrd_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer llwork;
    static doublereal smlnum;
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen), zungtr_(char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen);


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

#line 189 "zheev.f"
    /* Parameter adjustments */
#line 189 "zheev.f"
    a_dim1 = *lda;
#line 189 "zheev.f"
    a_offset = 1 + a_dim1;
#line 189 "zheev.f"
    a -= a_offset;
#line 189 "zheev.f"
    --w;
#line 189 "zheev.f"
    --work;
#line 189 "zheev.f"
    --rwork;
#line 189 "zheev.f"

#line 189 "zheev.f"
    /* Function Body */
#line 189 "zheev.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 190 "zheev.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 191 "zheev.f"
    lquery = *lwork == -1;

#line 193 "zheev.f"
    *info = 0;
#line 194 "zheev.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 195 "zheev.f"
	*info = -1;
#line 196 "zheev.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 197 "zheev.f"
	*info = -2;
#line 198 "zheev.f"
    } else if (*n < 0) {
#line 199 "zheev.f"
	*info = -3;
#line 200 "zheev.f"
    } else if (*lda < max(1,*n)) {
#line 201 "zheev.f"
	*info = -5;
#line 202 "zheev.f"
    }

#line 204 "zheev.f"
    if (*info == 0) {
#line 205 "zheev.f"
	nb = ilaenv_(&c__1, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 206 "zheev.f"
	i__1 = 1, i__2 = (nb + 1) * *n;
#line 206 "zheev.f"
	lwkopt = max(i__1,i__2);
#line 207 "zheev.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/* Computing MAX */
#line 209 "zheev.f"
	i__1 = 1, i__2 = (*n << 1) - 1;
#line 209 "zheev.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 209 "zheev.f"
	    *info = -8;
#line 209 "zheev.f"
	}
#line 211 "zheev.f"
    }

#line 213 "zheev.f"
    if (*info != 0) {
#line 214 "zheev.f"
	i__1 = -(*info);
#line 214 "zheev.f"
	xerbla_("ZHEEV ", &i__1, (ftnlen)6);
#line 215 "zheev.f"
	return 0;
#line 216 "zheev.f"
    } else if (lquery) {
#line 217 "zheev.f"
	return 0;
#line 218 "zheev.f"
    }

/*     Quick return if possible */

#line 222 "zheev.f"
    if (*n == 0) {
#line 223 "zheev.f"
	return 0;
#line 224 "zheev.f"
    }

#line 226 "zheev.f"
    if (*n == 1) {
#line 227 "zheev.f"
	i__1 = a_dim1 + 1;
#line 227 "zheev.f"
	w[1] = a[i__1].r;
#line 228 "zheev.f"
	work[1].r = 1., work[1].i = 0.;
#line 229 "zheev.f"
	if (wantz) {
#line 229 "zheev.f"
	    i__1 = a_dim1 + 1;
#line 229 "zheev.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 229 "zheev.f"
	}
#line 231 "zheev.f"
	return 0;
#line 232 "zheev.f"
    }

/*     Get machine constants. */

#line 236 "zheev.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 237 "zheev.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 238 "zheev.f"
    smlnum = safmin / eps;
#line 239 "zheev.f"
    bignum = 1. / smlnum;
#line 240 "zheev.f"
    rmin = sqrt(smlnum);
#line 241 "zheev.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 245 "zheev.f"
    anrm = zlanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);
#line 246 "zheev.f"
    iscale = 0;
#line 247 "zheev.f"
    if (anrm > 0. && anrm < rmin) {
#line 248 "zheev.f"
	iscale = 1;
#line 249 "zheev.f"
	sigma = rmin / anrm;
#line 250 "zheev.f"
    } else if (anrm > rmax) {
#line 251 "zheev.f"
	iscale = 1;
#line 252 "zheev.f"
	sigma = rmax / anrm;
#line 253 "zheev.f"
    }
#line 254 "zheev.f"
    if (iscale == 1) {
#line 254 "zheev.f"
	zlascl_(uplo, &c__0, &c__0, &c_b18, &sigma, n, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 254 "zheev.f"
    }

/*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form. */

#line 259 "zheev.f"
    inde = 1;
#line 260 "zheev.f"
    indtau = 1;
#line 261 "zheev.f"
    indwrk = indtau + *n;
#line 262 "zheev.f"
    llwork = *lwork - indwrk + 1;
#line 263 "zheev.f"
    zhetrd_(uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     ZUNGTR to generate the unitary matrix, then call ZSTEQR. */

#line 269 "zheev.f"
    if (! wantz) {
#line 270 "zheev.f"
	dsterf_(n, &w[1], &rwork[inde], info);
#line 271 "zheev.f"
    } else {
#line 272 "zheev.f"
	zungtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo, (ftnlen)1);
#line 274 "zheev.f"
	indwrk = inde + *n;
#line 275 "zheev.f"
	zsteqr_(jobz, n, &w[1], &rwork[inde], &a[a_offset], lda, &rwork[
		indwrk], info, (ftnlen)1);
#line 277 "zheev.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 281 "zheev.f"
    if (iscale == 1) {
#line 282 "zheev.f"
	if (*info == 0) {
#line 283 "zheev.f"
	    imax = *n;
#line 284 "zheev.f"
	} else {
#line 285 "zheev.f"
	    imax = *info - 1;
#line 286 "zheev.f"
	}
#line 287 "zheev.f"
	d__1 = 1. / sigma;
#line 287 "zheev.f"
	dscal_(&imax, &d__1, &w[1], &c__1);
#line 288 "zheev.f"
    }

/*     Set WORK(1) to optimal complex workspace size. */

#line 292 "zheev.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 294 "zheev.f"
    return 0;

/*     End of ZHEEV */

} /* zheev_ */

