#line 1 "zstedc.f"
/* zstedc.f -- translated by f2c (version 20100827).
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

#line 1 "zstedc.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;

/* > \brief \b ZSTEDC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSTEDC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstedc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstedc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstedc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, */
/*                          LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * ) */
/*       COMPLEX*16         WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSTEDC computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > The eigenvectors of a full or band complex Hermitian matrix can also */
/* > be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this */
/* > matrix to tridiagonal form. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none.  See DLAED3 for details. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only. */
/* >          = 'I':  Compute eigenvectors of tridiagonal matrix also. */
/* >          = 'V':  Compute eigenvectors of original Hermitian matrix */
/* >                  also.  On entry, Z contains the unitary matrix used */
/* >                  to reduce the original matrix to tridiagonal form. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The dimension of the symmetric tridiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the subdiagonal elements of the tridiagonal matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ,N) */
/* >          On entry, if COMPZ = 'V', then Z contains the unitary */
/* >          matrix used in the reduction to tridiagonal form. */
/* >          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the */
/* >          orthonormal eigenvectors of the original Hermitian matrix, */
/* >          and if COMPZ = 'I', Z contains the orthonormal eigenvectors */
/* >          of the symmetric tridiagonal matrix. */
/* >          If  COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1. */
/* >          If eigenvectors are desired, then LDZ >= max(1,N). */
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
/* >          The dimension of the array WORK. */
/* >          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1. */
/* >          If COMPZ = 'V' and N > 1, LWORK must be at least N*N. */
/* >          Note that for COMPZ = 'V', then if N is less than or */
/* >          equal to the minimum divide size, usually 25, then LWORK need */
/* >          only be 1. */
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
/* >          RWORK is DOUBLE PRECISION array, */
/* >                                         dimension (LRWORK) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of the array RWORK. */
/* >          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1. */
/* >          If COMPZ = 'V' and N > 1, LRWORK must be at least */
/* >                         1 + 3*N + 2*N*lg N + 4*N**2 , */
/* >                         where lg( N ) = smallest integer k such */
/* >                         that 2**k >= N. */
/* >          If COMPZ = 'I' and N > 1, LRWORK must be at least */
/* >                         1 + 4*N + 2*N**2 . */
/* >          Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* >          equal to the minimum divide size, usually 25, then LRWORK */
/* >          need only be max(1,2*(N-1)). */
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
/* >          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* >          If COMPZ = 'V' or N > 1,  LIWORK must be at least */
/* >                                    6 + 6*N + 5*N*lg N. */
/* >          If COMPZ = 'I' or N > 1,  LIWORK must be at least */
/* >                                    3 + 5*N . */
/* >          Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* >          equal to the minimum divide size, usually 25, then LIWORK */
/* >          need only be 1. */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  The algorithm failed to compute an eigenvalue while */
/* >                working on the submatrix lying in rows and columns */
/* >                INFO/(N+1) through mod(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int zstedc_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen compz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal p;
    static integer ii, ll, lgn;
    static doublereal eps, tiny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lwmin, start;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlaed0_(integer *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen), dlaset_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer finish;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), zlacrm_(integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    static integer liwmin, icompz;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), zlacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen);
    static doublereal orgnrm;
    static integer lrwmin;
    static logical lquery;
    static integer smlsiz;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 260 "zstedc.f"
    /* Parameter adjustments */
#line 260 "zstedc.f"
    --d__;
#line 260 "zstedc.f"
    --e;
#line 260 "zstedc.f"
    z_dim1 = *ldz;
#line 260 "zstedc.f"
    z_offset = 1 + z_dim1;
#line 260 "zstedc.f"
    z__ -= z_offset;
#line 260 "zstedc.f"
    --work;
#line 260 "zstedc.f"
    --rwork;
#line 260 "zstedc.f"
    --iwork;
#line 260 "zstedc.f"

#line 260 "zstedc.f"
    /* Function Body */
#line 260 "zstedc.f"
    *info = 0;
#line 261 "zstedc.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 263 "zstedc.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 264 "zstedc.f"
	icompz = 0;
#line 265 "zstedc.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 266 "zstedc.f"
	icompz = 1;
#line 267 "zstedc.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 268 "zstedc.f"
	icompz = 2;
#line 269 "zstedc.f"
    } else {
#line 270 "zstedc.f"
	icompz = -1;
#line 271 "zstedc.f"
    }
#line 272 "zstedc.f"
    if (icompz < 0) {
#line 273 "zstedc.f"
	*info = -1;
#line 274 "zstedc.f"
    } else if (*n < 0) {
#line 275 "zstedc.f"
	*info = -2;
#line 276 "zstedc.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 278 "zstedc.f"
	*info = -6;
#line 279 "zstedc.f"
    }

#line 281 "zstedc.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 285 "zstedc.f"
	smlsiz = ilaenv_(&c__9, "ZSTEDC", " ", &c__0, &c__0, &c__0, &c__0, (
		ftnlen)6, (ftnlen)1);
#line 286 "zstedc.f"
	if (*n <= 1 || icompz == 0) {
#line 287 "zstedc.f"
	    lwmin = 1;
#line 288 "zstedc.f"
	    liwmin = 1;
#line 289 "zstedc.f"
	    lrwmin = 1;
#line 290 "zstedc.f"
	} else if (*n <= smlsiz) {
#line 291 "zstedc.f"
	    lwmin = 1;
#line 292 "zstedc.f"
	    liwmin = 1;
#line 293 "zstedc.f"
	    lrwmin = *n - 1 << 1;
#line 294 "zstedc.f"
	} else if (icompz == 1) {
#line 295 "zstedc.f"
	    lgn = (integer) (log((doublereal) (*n)) / log(2.));
#line 296 "zstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 296 "zstedc.f"
		++lgn;
#line 296 "zstedc.f"
	    }
#line 298 "zstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 298 "zstedc.f"
		++lgn;
#line 298 "zstedc.f"
	    }
#line 300 "zstedc.f"
	    lwmin = *n * *n;
/* Computing 2nd power */
#line 301 "zstedc.f"
	    i__1 = *n;
#line 301 "zstedc.f"
	    lrwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
#line 302 "zstedc.f"
	    liwmin = *n * 6 + 6 + *n * 5 * lgn;
#line 303 "zstedc.f"
	} else if (icompz == 2) {
#line 304 "zstedc.f"
	    lwmin = 1;
/* Computing 2nd power */
#line 305 "zstedc.f"
	    i__1 = *n;
#line 305 "zstedc.f"
	    lrwmin = (*n << 2) + 1 + (i__1 * i__1 << 1);
#line 306 "zstedc.f"
	    liwmin = *n * 5 + 3;
#line 307 "zstedc.f"
	}
#line 308 "zstedc.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 309 "zstedc.f"
	rwork[1] = (doublereal) lrwmin;
#line 310 "zstedc.f"
	iwork[1] = liwmin;

#line 312 "zstedc.f"
	if (*lwork < lwmin && ! lquery) {
#line 313 "zstedc.f"
	    *info = -8;
#line 314 "zstedc.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 315 "zstedc.f"
	    *info = -10;
#line 316 "zstedc.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 317 "zstedc.f"
	    *info = -12;
#line 318 "zstedc.f"
	}
#line 319 "zstedc.f"
    }

#line 321 "zstedc.f"
    if (*info != 0) {
#line 322 "zstedc.f"
	i__1 = -(*info);
#line 322 "zstedc.f"
	xerbla_("ZSTEDC", &i__1, (ftnlen)6);
#line 323 "zstedc.f"
	return 0;
#line 324 "zstedc.f"
    } else if (lquery) {
#line 325 "zstedc.f"
	return 0;
#line 326 "zstedc.f"
    }

/*     Quick return if possible */

#line 330 "zstedc.f"
    if (*n == 0) {
#line 330 "zstedc.f"
	return 0;
#line 330 "zstedc.f"
    }
#line 332 "zstedc.f"
    if (*n == 1) {
#line 333 "zstedc.f"
	if (icompz != 0) {
#line 333 "zstedc.f"
	    i__1 = z_dim1 + 1;
#line 333 "zstedc.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 333 "zstedc.f"
	}
#line 335 "zstedc.f"
	return 0;
#line 336 "zstedc.f"
    }

/*     If the following conditional clause is removed, then the routine */
/*     will use the Divide and Conquer routine to compute only the */
/*     eigenvalues, which requires (3N + 3N**2) real workspace and */
/*     (2 + 5N + 2N lg(N)) integer workspace. */
/*     Since on many architectures DSTERF is much faster than any other */
/*     algorithm for finding eigenvalues only, it is used here */
/*     as the default. If the conditional clause is removed, then */
/*     information on the size of workspace needs to be changed. */

/*     If COMPZ = 'N', use DSTERF to compute the eigenvalues. */

#line 349 "zstedc.f"
    if (icompz == 0) {
#line 350 "zstedc.f"
	dsterf_(n, &d__[1], &e[1], info);
#line 351 "zstedc.f"
	goto L70;
#line 352 "zstedc.f"
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then */
/*     solve the problem with another solver. */

#line 357 "zstedc.f"
    if (*n <= smlsiz) {

#line 359 "zstedc.f"
	zsteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &rwork[1], 
		info, (ftnlen)1);

#line 361 "zstedc.f"
    } else {

/*        If COMPZ = 'I', we simply call DSTEDC instead. */

#line 365 "zstedc.f"
	if (icompz == 2) {
#line 366 "zstedc.f"
	    dlaset_("Full", n, n, &c_b17, &c_b18, &rwork[1], n, (ftnlen)4);
#line 367 "zstedc.f"
	    ll = *n * *n + 1;
#line 368 "zstedc.f"
	    i__1 = *lrwork - ll + 1;
#line 368 "zstedc.f"
	    dstedc_("I", n, &d__[1], &e[1], &rwork[1], n, &rwork[ll], &i__1, &
		    iwork[1], liwork, info, (ftnlen)1);
#line 370 "zstedc.f"
	    i__1 = *n;
#line 370 "zstedc.f"
	    for (j = 1; j <= i__1; ++j) {
#line 371 "zstedc.f"
		i__2 = *n;
#line 371 "zstedc.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 372 "zstedc.f"
		    i__3 = i__ + j * z_dim1;
#line 372 "zstedc.f"
		    i__4 = (j - 1) * *n + i__;
#line 372 "zstedc.f"
		    z__[i__3].r = rwork[i__4], z__[i__3].i = 0.;
#line 373 "zstedc.f"
/* L10: */
#line 373 "zstedc.f"
		}
#line 374 "zstedc.f"
/* L20: */
#line 374 "zstedc.f"
	    }
#line 375 "zstedc.f"
	    goto L70;
#line 376 "zstedc.f"
	}

/*        From now on, only option left to be handled is COMPZ = 'V', */
/*        i.e. ICOMPZ = 1. */

/*        Scale. */

#line 383 "zstedc.f"
	orgnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 384 "zstedc.f"
	if (orgnrm == 0.) {
#line 384 "zstedc.f"
	    goto L70;
#line 384 "zstedc.f"
	}

#line 387 "zstedc.f"
	eps = dlamch_("Epsilon", (ftnlen)7);

#line 389 "zstedc.f"
	start = 1;

/*        while ( START <= N ) */

#line 393 "zstedc.f"
L30:
#line 394 "zstedc.f"
	if (start <= *n) {

/*           Let FINISH be the position of the next subdiagonal entry */
/*           such that E( FINISH ) <= TINY or FINISH = N if no such */
/*           subdiagonal exists.  The matrix identified by the elements */
/*           between START and FINISH constitutes an independent */
/*           sub-problem. */

#line 402 "zstedc.f"
	    finish = start;
#line 403 "zstedc.f"
L40:
#line 404 "zstedc.f"
	    if (finish < *n) {
#line 405 "zstedc.f"
		tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) * sqrt((
			d__2 = d__[finish + 1], abs(d__2)));
#line 407 "zstedc.f"
		if ((d__1 = e[finish], abs(d__1)) > tiny) {
#line 408 "zstedc.f"
		    ++finish;
#line 409 "zstedc.f"
		    goto L40;
#line 410 "zstedc.f"
		}
#line 411 "zstedc.f"
	    }

/*           (Sub) Problem determined.  Compute its size and solve it. */

#line 415 "zstedc.f"
	    m = finish - start + 1;
#line 416 "zstedc.f"
	    if (m > smlsiz) {

/*              Scale. */

#line 420 "zstedc.f"
		orgnrm = dlanst_("M", &m, &d__[start], &e[start], (ftnlen)1);
#line 421 "zstedc.f"
		dlascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);
#line 423 "zstedc.f"
		i__1 = m - 1;
#line 423 "zstedc.f"
		i__2 = m - 1;
#line 423 "zstedc.f"
		dlascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[
			start], &i__2, info, (ftnlen)1);

#line 426 "zstedc.f"
		zlaed0_(n, &m, &d__[start], &e[start], &z__[start * z_dim1 + 
			1], ldz, &work[1], n, &rwork[1], &iwork[1], info);
#line 428 "zstedc.f"
		if (*info > 0) {
#line 429 "zstedc.f"
		    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info %
			     (m + 1) + start - 1;
#line 431 "zstedc.f"
		    goto L70;
#line 432 "zstedc.f"
		}

/*              Scale back. */

#line 436 "zstedc.f"
		dlascl_("G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);

#line 439 "zstedc.f"
	    } else {
#line 440 "zstedc.f"
		dsteqr_("I", &m, &d__[start], &e[start], &rwork[1], &m, &
			rwork[m * m + 1], info, (ftnlen)1);
#line 442 "zstedc.f"
		zlacrm_(n, &m, &z__[start * z_dim1 + 1], ldz, &rwork[1], &m, &
			work[1], n, &rwork[m * m + 1]);
#line 444 "zstedc.f"
		zlacpy_("A", n, &m, &work[1], n, &z__[start * z_dim1 + 1], 
			ldz, (ftnlen)1);
#line 445 "zstedc.f"
		if (*info > 0) {
#line 446 "zstedc.f"
		    *info = start * (*n + 1) + finish;
#line 447 "zstedc.f"
		    goto L70;
#line 448 "zstedc.f"
		}
#line 449 "zstedc.f"
	    }

#line 451 "zstedc.f"
	    start = finish + 1;
#line 452 "zstedc.f"
	    goto L30;
#line 453 "zstedc.f"
	}

/*        endwhile */


/*        Use Selection Sort to minimize swaps of eigenvectors */

#line 460 "zstedc.f"
	i__1 = *n;
#line 460 "zstedc.f"
	for (ii = 2; ii <= i__1; ++ii) {
#line 461 "zstedc.f"
	    i__ = ii - 1;
#line 462 "zstedc.f"
	    k = i__;
#line 463 "zstedc.f"
	    p = d__[i__];
#line 464 "zstedc.f"
	    i__2 = *n;
#line 464 "zstedc.f"
	    for (j = ii; j <= i__2; ++j) {
#line 465 "zstedc.f"
		if (d__[j] < p) {
#line 466 "zstedc.f"
		    k = j;
#line 467 "zstedc.f"
		    p = d__[j];
#line 468 "zstedc.f"
		}
#line 469 "zstedc.f"
/* L50: */
#line 469 "zstedc.f"
	    }
#line 470 "zstedc.f"
	    if (k != i__) {
#line 471 "zstedc.f"
		d__[k] = d__[i__];
#line 472 "zstedc.f"
		d__[i__] = p;
#line 473 "zstedc.f"
		zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
#line 474 "zstedc.f"
	    }
#line 475 "zstedc.f"
/* L60: */
#line 475 "zstedc.f"
	}
#line 476 "zstedc.f"
    }

#line 478 "zstedc.f"
L70:
#line 479 "zstedc.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 480 "zstedc.f"
    rwork[1] = (doublereal) lrwmin;
#line 481 "zstedc.f"
    iwork[1] = liwmin;

#line 483 "zstedc.f"
    return 0;

/*     End of ZSTEDC */

} /* zstedc_ */

