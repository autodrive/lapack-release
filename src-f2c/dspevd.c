#line 1 "dspevd.f"
/* dspevd.f -- translated by f2c (version 20100827).
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

#line 1 "dspevd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DSPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, */
/*                          IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPEVD computes all the eigenvalues and, optionally, eigenvectors */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \date June 2017 */

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dspevd_(char *jobz, char *uplo, integer *n, doublereal *
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
    static doublereal anrm, rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo, lwmin;
    static logical wantz;
    extern doublereal dlamch_(char *, ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlansp_(char *, char *, integer *, doublereal *, 
	    doublereal *, ftnlen, ftnlen);
    static integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int dsptrd_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dopmtr_(char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer llwork;
    static doublereal smlnum;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 223 "dspevd.f"
    /* Parameter adjustments */
#line 223 "dspevd.f"
    --ap;
#line 223 "dspevd.f"
    --w;
#line 223 "dspevd.f"
    z_dim1 = *ldz;
#line 223 "dspevd.f"
    z_offset = 1 + z_dim1;
#line 223 "dspevd.f"
    z__ -= z_offset;
#line 223 "dspevd.f"
    --work;
#line 223 "dspevd.f"
    --iwork;
#line 223 "dspevd.f"

#line 223 "dspevd.f"
    /* Function Body */
#line 223 "dspevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 224 "dspevd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 226 "dspevd.f"
    *info = 0;
#line 227 "dspevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 228 "dspevd.f"
	*info = -1;
#line 229 "dspevd.f"
    } else if (! (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1))) {
#line 231 "dspevd.f"
	*info = -2;
#line 232 "dspevd.f"
    } else if (*n < 0) {
#line 233 "dspevd.f"
	*info = -3;
#line 234 "dspevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 235 "dspevd.f"
	*info = -7;
#line 236 "dspevd.f"
    }

#line 238 "dspevd.f"
    if (*info == 0) {
#line 239 "dspevd.f"
	if (*n <= 1) {
#line 240 "dspevd.f"
	    liwmin = 1;
#line 241 "dspevd.f"
	    lwmin = 1;
#line 242 "dspevd.f"
	} else {
#line 243 "dspevd.f"
	    if (wantz) {
#line 244 "dspevd.f"
		liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 245 "dspevd.f"
		i__1 = *n;
#line 245 "dspevd.f"
		lwmin = *n * 6 + 1 + i__1 * i__1;
#line 246 "dspevd.f"
	    } else {
#line 247 "dspevd.f"
		liwmin = 1;
#line 248 "dspevd.f"
		lwmin = *n << 1;
#line 249 "dspevd.f"
	    }
#line 250 "dspevd.f"
	}
#line 251 "dspevd.f"
	iwork[1] = liwmin;
#line 252 "dspevd.f"
	work[1] = (doublereal) lwmin;

#line 254 "dspevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 255 "dspevd.f"
	    *info = -9;
#line 256 "dspevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 257 "dspevd.f"
	    *info = -11;
#line 258 "dspevd.f"
	}
#line 259 "dspevd.f"
    }

#line 261 "dspevd.f"
    if (*info != 0) {
#line 262 "dspevd.f"
	i__1 = -(*info);
#line 262 "dspevd.f"
	xerbla_("DSPEVD", &i__1, (ftnlen)6);
#line 263 "dspevd.f"
	return 0;
#line 264 "dspevd.f"
    } else if (lquery) {
#line 265 "dspevd.f"
	return 0;
#line 266 "dspevd.f"
    }

/*     Quick return if possible */

#line 270 "dspevd.f"
    if (*n == 0) {
#line 270 "dspevd.f"
	return 0;
#line 270 "dspevd.f"
    }

#line 273 "dspevd.f"
    if (*n == 1) {
#line 274 "dspevd.f"
	w[1] = ap[1];
#line 275 "dspevd.f"
	if (wantz) {
#line 275 "dspevd.f"
	    z__[z_dim1 + 1] = 1.;
#line 275 "dspevd.f"
	}
#line 277 "dspevd.f"
	return 0;
#line 278 "dspevd.f"
    }

/*     Get machine constants. */

#line 282 "dspevd.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 283 "dspevd.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 284 "dspevd.f"
    smlnum = safmin / eps;
#line 285 "dspevd.f"
    bignum = 1. / smlnum;
#line 286 "dspevd.f"
    rmin = sqrt(smlnum);
#line 287 "dspevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 291 "dspevd.f"
    anrm = dlansp_("M", uplo, n, &ap[1], &work[1], (ftnlen)1, (ftnlen)1);
#line 292 "dspevd.f"
    iscale = 0;
#line 293 "dspevd.f"
    if (anrm > 0. && anrm < rmin) {
#line 294 "dspevd.f"
	iscale = 1;
#line 295 "dspevd.f"
	sigma = rmin / anrm;
#line 296 "dspevd.f"
    } else if (anrm > rmax) {
#line 297 "dspevd.f"
	iscale = 1;
#line 298 "dspevd.f"
	sigma = rmax / anrm;
#line 299 "dspevd.f"
    }
#line 300 "dspevd.f"
    if (iscale == 1) {
#line 301 "dspevd.f"
	i__1 = *n * (*n + 1) / 2;
#line 301 "dspevd.f"
	dscal_(&i__1, &sigma, &ap[1], &c__1);
#line 302 "dspevd.f"
    }

/*     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form. */

#line 306 "dspevd.f"
    inde = 1;
#line 307 "dspevd.f"
    indtau = inde + *n;
#line 308 "dspevd.f"
    dsptrd_(uplo, n, &ap[1], &w[1], &work[inde], &work[indtau], &iinfo, (
	    ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call DOPMTR to multiply it by the */
/*     Householder transformations represented in AP. */

#line 315 "dspevd.f"
    if (! wantz) {
#line 316 "dspevd.f"
	dsterf_(n, &w[1], &work[inde], info);
#line 317 "dspevd.f"
    } else {
#line 318 "dspevd.f"
	indwrk = indtau + *n;
#line 319 "dspevd.f"
	llwork = *lwork - indwrk + 1;
#line 320 "dspevd.f"
	dstedc_("I", n, &w[1], &work[inde], &z__[z_offset], ldz, &work[indwrk]
		, &llwork, &iwork[1], liwork, info, (ftnlen)1);
#line 322 "dspevd.f"
	dopmtr_("L", uplo, "N", n, n, &ap[1], &work[indtau], &z__[z_offset], 
		ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 324 "dspevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 328 "dspevd.f"
    if (iscale == 1) {
#line 328 "dspevd.f"
	d__1 = 1. / sigma;
#line 328 "dspevd.f"
	dscal_(n, &d__1, &w[1], &c__1);
#line 328 "dspevd.f"
    }

#line 331 "dspevd.f"
    work[1] = (doublereal) lwmin;
#line 332 "dspevd.f"
    iwork[1] = liwmin;
#line 333 "dspevd.f"
    return 0;

/*     End of DSPEVD */

} /* dspevd_ */

