#line 1 "sstevd.f"
/* sstevd.f -- translated by f2c (version 20100827).
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

#line 1 "sstevd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SSTEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER 
matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ */
/*       INTEGER            INFO, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEVD computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric tridiagonal matrix. If eigenvectors are desired, it */
/* > uses a divide and conquer algorithm. */
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
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A, stored in elements 1 to N-1 of E. */
/* >          On exit, the contents of E are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* >          eigenvectors of the matrix A, with the i-th column of Z */
/* >          holding the eigenvector associated with D(i). */
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
/* >          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1 then LWORK must be at least */
/* >                         ( 1 + 4*N + N**2 ). */
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
/* >          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N. */
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
/* >                off-diagonal elements of E did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int sstevd_(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal eps, rmin, rmax, tnrm, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer lwmin;
    static logical wantz;
    static integer iscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int sstedc_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, ftnlen);
    static integer liwmin;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
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

#line 207 "sstevd.f"
    /* Parameter adjustments */
#line 207 "sstevd.f"
    --d__;
#line 207 "sstevd.f"
    --e;
#line 207 "sstevd.f"
    z_dim1 = *ldz;
#line 207 "sstevd.f"
    z_offset = 1 + z_dim1;
#line 207 "sstevd.f"
    z__ -= z_offset;
#line 207 "sstevd.f"
    --work;
#line 207 "sstevd.f"
    --iwork;
#line 207 "sstevd.f"

#line 207 "sstevd.f"
    /* Function Body */
#line 207 "sstevd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 208 "sstevd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 210 "sstevd.f"
    *info = 0;
#line 211 "sstevd.f"
    liwmin = 1;
#line 212 "sstevd.f"
    lwmin = 1;
#line 213 "sstevd.f"
    if (*n > 1 && wantz) {
/* Computing 2nd power */
#line 214 "sstevd.f"
	i__1 = *n;
#line 214 "sstevd.f"
	lwmin = (*n << 2) + 1 + i__1 * i__1;
#line 215 "sstevd.f"
	liwmin = *n * 5 + 3;
#line 216 "sstevd.f"
    }

#line 218 "sstevd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 219 "sstevd.f"
	*info = -1;
#line 220 "sstevd.f"
    } else if (*n < 0) {
#line 221 "sstevd.f"
	*info = -2;
#line 222 "sstevd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 223 "sstevd.f"
	*info = -6;
#line 224 "sstevd.f"
    }

#line 226 "sstevd.f"
    if (*info == 0) {
#line 227 "sstevd.f"
	work[1] = (doublereal) lwmin;
#line 228 "sstevd.f"
	iwork[1] = liwmin;

#line 230 "sstevd.f"
	if (*lwork < lwmin && ! lquery) {
#line 231 "sstevd.f"
	    *info = -8;
#line 232 "sstevd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 233 "sstevd.f"
	    *info = -10;
#line 234 "sstevd.f"
	}
#line 235 "sstevd.f"
    }

#line 237 "sstevd.f"
    if (*info != 0) {
#line 238 "sstevd.f"
	i__1 = -(*info);
#line 238 "sstevd.f"
	xerbla_("SSTEVD", &i__1, (ftnlen)6);
#line 239 "sstevd.f"
	return 0;
#line 240 "sstevd.f"
    } else if (lquery) {
#line 241 "sstevd.f"
	return 0;
#line 242 "sstevd.f"
    }

/*     Quick return if possible */

#line 246 "sstevd.f"
    if (*n == 0) {
#line 246 "sstevd.f"
	return 0;
#line 246 "sstevd.f"
    }

#line 249 "sstevd.f"
    if (*n == 1) {
#line 250 "sstevd.f"
	if (wantz) {
#line 250 "sstevd.f"
	    z__[z_dim1 + 1] = 1.;
#line 250 "sstevd.f"
	}
#line 252 "sstevd.f"
	return 0;
#line 253 "sstevd.f"
    }

/*     Get machine constants. */

#line 257 "sstevd.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 258 "sstevd.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 259 "sstevd.f"
    smlnum = safmin / eps;
#line 260 "sstevd.f"
    bignum = 1. / smlnum;
#line 261 "sstevd.f"
    rmin = sqrt(smlnum);
#line 262 "sstevd.f"
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

#line 266 "sstevd.f"
    iscale = 0;
#line 267 "sstevd.f"
    tnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 268 "sstevd.f"
    if (tnrm > 0. && tnrm < rmin) {
#line 269 "sstevd.f"
	iscale = 1;
#line 270 "sstevd.f"
	sigma = rmin / tnrm;
#line 271 "sstevd.f"
    } else if (tnrm > rmax) {
#line 272 "sstevd.f"
	iscale = 1;
#line 273 "sstevd.f"
	sigma = rmax / tnrm;
#line 274 "sstevd.f"
    }
#line 275 "sstevd.f"
    if (iscale == 1) {
#line 276 "sstevd.f"
	sscal_(n, &sigma, &d__[1], &c__1);
#line 277 "sstevd.f"
	i__1 = *n - 1;
#line 277 "sstevd.f"
	sscal_(&i__1, &sigma, &e[1], &c__1);
#line 278 "sstevd.f"
    }

/*     For eigenvalues only, call SSTERF.  For eigenvalues and */
/*     eigenvectors, call SSTEDC. */

#line 283 "sstevd.f"
    if (! wantz) {
#line 284 "sstevd.f"
	ssterf_(n, &d__[1], &e[1], info);
#line 285 "sstevd.f"
    } else {
#line 286 "sstevd.f"
	sstedc_("I", n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], lwork, 
		&iwork[1], liwork, info, (ftnlen)1);
#line 288 "sstevd.f"
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

#line 292 "sstevd.f"
    if (iscale == 1) {
#line 292 "sstevd.f"
	d__1 = 1. / sigma;
#line 292 "sstevd.f"
	sscal_(n, &d__1, &d__[1], &c__1);
#line 292 "sstevd.f"
    }

#line 295 "sstevd.f"
    work[1] = (doublereal) lwmin;
#line 296 "sstevd.f"
    iwork[1] = liwmin;

#line 298 "sstevd.f"
    return 0;

/*     End of SSTEVD */

} /* sstevd_ */

