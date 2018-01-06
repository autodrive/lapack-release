#line 1 "sspevd.f"
/* sspevd.f -- translated by f2c (version 20100827).
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

#line 1 "sspevd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SSPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, */
/*                          IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPEVD computes all the eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A in packed storage. If eigenvectors are */
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
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK must be at least 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be at least */
/* >                                                 1 + 6*N + N**2. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the required sizes of the WORK and IWORK */
/* >          arrays, returns these values as the first entries of the WORK */
/* >          and IWORK arrays, and no error message related to LWORK or */
/* >          LIWORK is issued by XERBLA. */
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
/* >          The dimension of the array IWORK. */
/* >          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK and IWORK arrays, and no error message related to */
/* >          LWORK or LIWORK is issued by XERBLA. */
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

/* > \date December 2016 */

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int sspevd_(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
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
	    integer *);
    static integer lwmin;
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer indtau;
    extern /* Subroutine */ int sstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, ftnlen);
    static integer indwrk, liwmin;
    extern doublereal slansp_(char *, char *, integer *, doublereal *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer llwork;
    static doublereal smlnum;
    extern /* Subroutine */ int ssptrd_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int sopmtr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);


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

#line 223 "sspevd.f"
    /* Parameter adjustments */
#line 223 "sspevd.f"
    --ap;
#line 223 "sspevd.f"
    --w;
#line 223 "sspevd.f"
    z_dim1 = *ldz;
#line 223 "sspevd.f"
    z_offset = 1 + z_dim1;
#line 223 "sspevd.f"
    z__ -= z_offset;
#line 223 "sspevd.f"
    --work;
#line 223 "sspevd.f"
    --iwork;
#line 223 "sspevd.f"

#line 223 "sspevd.f"
    /* Function Body */
#line 223 "sspevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 224 "sspevd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 226 "sspevd.f"
    *info = 0;
#line 227 "sspevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 228 "sspevd.f"
	*info = -1;
#line 229 "sspevd.f"
    } else if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1))) {
#line 231 "sspevd.f"
	*info = -2;
#line 232 "sspevd.f"
    } else if (*n < 0) {
#line 233 "sspevd.f"
	*info = -3;
#line 234 "sspevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 235 "sspevd.f"
	*info = -7;
#line 236 "sspevd.f"
    }

#line 238 "sspevd.f"
    if (*info == 0) {
#line 239 "sspevd.f"
	if (*n <= 1) {
#line 240 "sspevd.f"
	    liwmin = 1;
#line 241 "sspevd.f"
	    lwmin = 1;
#line 242 "sspevd.f"
	} else {
#line 243 "sspevd.f"
	    if (wantz) {
#line 244 "sspevd.f"
		liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 245 "sspevd.f"
		i__1 = *n;
#line 245 "sspevd.f"
		lwmin = *n * 6 + 1 + i__1 * i__1;
#line 246 "sspevd.f"
	    } else {
#line 247 "sspevd.f"
		liwmin = 1;
#line 248 "sspevd.f"
		lwmin = *n << 1;
#line 249 "sspevd.f"
	    }
#line 250 "sspevd.f"
	}
#line 251 "sspevd.f"
	iwork[1] = liwmin;
#line 252 "sspevd.f"
	work[1] = (doublereal) lwmin;

#line 254 "sspevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 255 "sspevd.f"
	    *info = -9;
#line 256 "sspevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 257 "sspevd.f"
	    *info = -11;
#line 258 "sspevd.f"
	}
#line 259 "sspevd.f"
    }

#line 261 "sspevd.f"
    if (*info != 0) {
#line 262 "sspevd.f"
	i__1 = -(*info);
#line 262 "sspevd.f"
	xerbla_("SSPEVD", &i__1, (ftnlen)6);
#line 263 "sspevd.f"
	return 0;
#line 264 "sspevd.f"
    } else if (lquery) {
#line 265 "sspevd.f"
	return 0;
#line 266 "sspevd.f"
    }

/*     Quick return if possible */

#line 270 "sspevd.f"
    if (*n == 0) {
#line 270 "sspevd.f"
	return 0;
#line 270 "sspevd.f"
    }

#line 273 "sspevd.f"
    if (*n == 1) {
#line 274 "sspevd.f"
	w[1] = ap[1];
#line 275 "sspevd.f"
	if (wantz) {
#line 275 "sspevd.f"
	    z__[z_dim1 + 1] = 1.;
#line 275 "sspevd.f"
	}
#line 277 "sspevd.f"
	return 0;
#line 278 "sspevd.f"
    }

/*     Get machine constants. */

#line 282 "sspevd.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 283 "sspevd.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 284 "sspevd.f"
    smlnum = safmin / eps;
#line 285 "sspevd.f"
    bignum = 1. / smlnum;
#line 286 "sspevd.f"
    rmin = sqrt(smlnum);
#line 287 "sspevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 291 "sspevd.f"
    anrm = slansp_("M", uplo, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)1);
#line 292 "sspevd.f"
    iscale = 0;
#line 293 "sspevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 294 "sspevd.f"
	iscale = 1;
#line 295 "sspevd.f"
	sigma = rmin / anrm;
#line 296 "sspevd.f"
    } else if (anrm > rmax) {
#line 297 "sspevd.f"
	iscale = 1;
#line 298 "sspevd.f"
	sigma = rmax / anrm;
#line 299 "sspevd.f"
    }
#line 300 "sspevd.f"
    if (iscale == 1) {
#line 301 "sspevd.f"
	i__1 = *n * (*n + 1) / 2;
#line 301 "sspevd.f"
	sscal_(&i__1, &sigma, &ap[1], &c__1);
#line 302 "sspevd.f"
    }

/*     Call SSPTRD to reduce symmetric packed matrix to tridiagonal form. */

#line 306 "sspevd.f"
    inde = 1;
#line 307 "sspevd.f"
    indtau = inde + *n;
#line 308 "sspevd.f"
    ssptrd_(uplo, n, &ap[1], &w[1], &work[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call */
/*     SSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call SOPMTR to multiply it by the */
/*     Householder transformations represented in AP. */

#line 315 "sspevd.f"
    if (! wantz) {
#line 316 "sspevd.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 317 "sspevd.f"
    } else {
#line 318 "sspevd.f"
	indwrk = indtau + *n;
#line 319 "sspevd.f"
	llwork = *lwork - indwrk + 1;
#line 320 "sspevd.f"
	sstedc_("I", n, &w[1], &work[inde], &z__[z_offset], ldz, &work[indwrk]
		, &llwork, &iwork[1], liwork, info, (ftnlen)1);
#line 322 "sspevd.f"
	sopmtr_("L", uplo, "N", n, n, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 324 "sspevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 328 "sspevd.f"
    if (iscale == 1) {
#line 328 "sspevd.f"
	d__1 = 1. / sigma;
#line 328 "sspevd.f"
	sscal_(n, &d__1, &w[1], &c__1);
#line 328 "sspevd.f"
    }

#line 331 "sspevd.f"
    work[1] = (doublereal) lwmin;
#line 332 "sspevd.f"
    iwork[1] = liwmin;
#line 333 "sspevd.f"
    return 0;

/*     End of SSPEVD */

} /* sspevd_ */

