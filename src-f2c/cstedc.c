#line 1 "cstedc.f"
/* cstedc.f -- translated by f2c (version 20100827).
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

#line 1 "cstedc.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;

/* > \brief \b CSTEDC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSTEDC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cstedc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cstedc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cstedc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, */
/*                          LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), E( * ), RWORK( * ) */
/*       COMPLEX            WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSTEDC computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > The eigenvectors of a full or band complex Hermitian matrix can also */
/* > be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this */
/* > matrix to tridiagonal form. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none.  See SLAED3 for details. */
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
/* >          D is REAL array, dimension (N) */
/* >          On entry, the diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the subdiagonal elements of the tridiagonal matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int cstedc_(char *compz, integer *n, doublereal *d__, 
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
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer lwmin;
    extern /* Subroutine */ int claed0_(integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *);
    static integer start;
    extern /* Subroutine */ int clacrm_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer finish;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), sstedc_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen), slaset_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static integer liwmin, icompz;
    extern /* Subroutine */ int csteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     ftnlen);
    static doublereal orgnrm;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *);
    static integer lrwmin;
    static logical lquery;
    static integer smlsiz;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

#line 259 "cstedc.f"
    /* Parameter adjustments */
#line 259 "cstedc.f"
    --d__;
#line 259 "cstedc.f"
    --e;
#line 259 "cstedc.f"
    z_dim1 = *ldz;
#line 259 "cstedc.f"
    z_offset = 1 + z_dim1;
#line 259 "cstedc.f"
    z__ -= z_offset;
#line 259 "cstedc.f"
    --work;
#line 259 "cstedc.f"
    --rwork;
#line 259 "cstedc.f"
    --iwork;
#line 259 "cstedc.f"

#line 259 "cstedc.f"
    /* Function Body */
#line 259 "cstedc.f"
    *info = 0;
#line 260 "cstedc.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 262 "cstedc.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 263 "cstedc.f"
	icompz = 0;
#line 264 "cstedc.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 265 "cstedc.f"
	icompz = 1;
#line 266 "cstedc.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 267 "cstedc.f"
	icompz = 2;
#line 268 "cstedc.f"
    } else {
#line 269 "cstedc.f"
	icompz = -1;
#line 270 "cstedc.f"
    }
#line 271 "cstedc.f"
    if (icompz < 0) {
#line 272 "cstedc.f"
	*info = -1;
#line 273 "cstedc.f"
    } else if (*n < 0) {
#line 274 "cstedc.f"
	*info = -2;
#line 275 "cstedc.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 277 "cstedc.f"
	*info = -6;
#line 278 "cstedc.f"
    }

#line 280 "cstedc.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 284 "cstedc.f"
	smlsiz = ilaenv_(&c__9, "CSTEDC", " ", &c__0, &c__0, &c__0, &c__0, (
		ftnlen)6, (ftnlen)1);
#line 285 "cstedc.f"
	if (*n <= 1 || icompz == 0) {
#line 286 "cstedc.f"
	    lwmin = 1;
#line 287 "cstedc.f"
	    liwmin = 1;
#line 288 "cstedc.f"
	    lrwmin = 1;
#line 289 "cstedc.f"
	} else if (*n <= smlsiz) {
#line 290 "cstedc.f"
	    lwmin = 1;
#line 291 "cstedc.f"
	    liwmin = 1;
#line 292 "cstedc.f"
	    lrwmin = *n - 1 << 1;
#line 293 "cstedc.f"
	} else if (icompz == 1) {
#line 294 "cstedc.f"
	    lgn = (integer) (log((doublereal) (*n)) / log(2.));
#line 295 "cstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 295 "cstedc.f"
		++lgn;
#line 295 "cstedc.f"
	    }
#line 297 "cstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 297 "cstedc.f"
		++lgn;
#line 297 "cstedc.f"
	    }
#line 299 "cstedc.f"
	    lwmin = *n * *n;
/* Computing 2nd power */
#line 300 "cstedc.f"
	    i__1 = *n;
#line 300 "cstedc.f"
	    lrwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
#line 301 "cstedc.f"
	    liwmin = *n * 6 + 6 + *n * 5 * lgn;
#line 302 "cstedc.f"
	} else if (icompz == 2) {
#line 303 "cstedc.f"
	    lwmin = 1;
/* Computing 2nd power */
#line 304 "cstedc.f"
	    i__1 = *n;
#line 304 "cstedc.f"
	    lrwmin = (*n << 2) + 1 + (i__1 * i__1 << 1);
#line 305 "cstedc.f"
	    liwmin = *n * 5 + 3;
#line 306 "cstedc.f"
	}
#line 307 "cstedc.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 308 "cstedc.f"
	rwork[1] = (doublereal) lrwmin;
#line 309 "cstedc.f"
	iwork[1] = liwmin;

#line 311 "cstedc.f"
	if (*lwork < lwmin && ! lquery) {
#line 312 "cstedc.f"
	    *info = -8;
#line 313 "cstedc.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 314 "cstedc.f"
	    *info = -10;
#line 315 "cstedc.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 316 "cstedc.f"
	    *info = -12;
#line 317 "cstedc.f"
	}
#line 318 "cstedc.f"
    }

#line 320 "cstedc.f"
    if (*info != 0) {
#line 321 "cstedc.f"
	i__1 = -(*info);
#line 321 "cstedc.f"
	xerbla_("CSTEDC", &i__1, (ftnlen)6);
#line 322 "cstedc.f"
	return 0;
#line 323 "cstedc.f"
    } else if (lquery) {
#line 324 "cstedc.f"
	return 0;
#line 325 "cstedc.f"
    }

/*     Quick return if possible */

#line 329 "cstedc.f"
    if (*n == 0) {
#line 329 "cstedc.f"
	return 0;
#line 329 "cstedc.f"
    }
#line 331 "cstedc.f"
    if (*n == 1) {
#line 332 "cstedc.f"
	if (icompz != 0) {
#line 332 "cstedc.f"
	    i__1 = z_dim1 + 1;
#line 332 "cstedc.f"
	    z__[i__1].r = 1., z__[i__1].i = 0.;
#line 332 "cstedc.f"
	}
#line 334 "cstedc.f"
	return 0;
#line 335 "cstedc.f"
    }

/*     If the following conditional clause is removed, then the routine */
/*     will use the Divide and Conquer routine to compute only the */
/*     eigenvalues, which requires (3N + 3N**2) real workspace and */
/*     (2 + 5N + 2N lg(N)) integer workspace. */
/*     Since on many architectures SSTERF is much faster than any other */
/*     algorithm for finding eigenvalues only, it is used here */
/*     as the default. If the conditional clause is removed, then */
/*     information on the size of workspace needs to be changed. */

/*     If COMPZ = 'N', use SSTERF to compute the eigenvalues. */

#line 348 "cstedc.f"
    if (icompz == 0) {
#line 349 "cstedc.f"
	ssterf_(n, &d__[1], &e[1], info);
#line 350 "cstedc.f"
	goto L70;
#line 351 "cstedc.f"
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then */
/*     solve the problem with another solver. */

#line 356 "cstedc.f"
    if (*n <= smlsiz) {

#line 358 "cstedc.f"
	csteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &rwork[1], 
		info, (ftnlen)1);

#line 360 "cstedc.f"
    } else {

/*        If COMPZ = 'I', we simply call SSTEDC instead. */

#line 364 "cstedc.f"
	if (icompz == 2) {
#line 365 "cstedc.f"
	    slaset_("Full", n, n, &c_b17, &c_b18, &rwork[1], n, (ftnlen)4);
#line 366 "cstedc.f"
	    ll = *n * *n + 1;
#line 367 "cstedc.f"
	    i__1 = *lrwork - ll + 1;
#line 367 "cstedc.f"
	    sstedc_("I", n, &d__[1], &e[1], &rwork[1], n, &rwork[ll], &i__1, &
		    iwork[1], liwork, info, (ftnlen)1);
#line 369 "cstedc.f"
	    i__1 = *n;
#line 369 "cstedc.f"
	    for (j = 1; j <= i__1; ++j) {
#line 370 "cstedc.f"
		i__2 = *n;
#line 370 "cstedc.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 371 "cstedc.f"
		    i__3 = i__ + j * z_dim1;
#line 371 "cstedc.f"
		    i__4 = (j - 1) * *n + i__;
#line 371 "cstedc.f"
		    z__[i__3].r = rwork[i__4], z__[i__3].i = 0.;
#line 372 "cstedc.f"
/* L10: */
#line 372 "cstedc.f"
		}
#line 373 "cstedc.f"
/* L20: */
#line 373 "cstedc.f"
	    }
#line 374 "cstedc.f"
	    goto L70;
#line 375 "cstedc.f"
	}

/*        From now on, only option left to be handled is COMPZ = 'V', */
/*        i.e. ICOMPZ = 1. */

/*        Scale. */

#line 382 "cstedc.f"
	orgnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 383 "cstedc.f"
	if (orgnrm == 0.) {
#line 383 "cstedc.f"
	    goto L70;
#line 383 "cstedc.f"
	}

#line 386 "cstedc.f"
	eps = slamch_("Epsilon", (ftnlen)7);

#line 388 "cstedc.f"
	start = 1;

/*        while ( START <= N ) */

#line 392 "cstedc.f"
L30:
#line 393 "cstedc.f"
	if (start <= *n) {

/*           Let FINISH be the position of the next subdiagonal entry */
/*           such that E( FINISH ) <= TINY or FINISH = N if no such */
/*           subdiagonal exists.  The matrix identified by the elements */
/*           between START and FINISH constitutes an independent */
/*           sub-problem. */

#line 401 "cstedc.f"
	    finish = start;
#line 402 "cstedc.f"
L40:
#line 403 "cstedc.f"
	    if (finish < *n) {
#line 404 "cstedc.f"
		tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) * sqrt((
			d__2 = d__[finish + 1], abs(d__2)));
#line 406 "cstedc.f"
		if ((d__1 = e[finish], abs(d__1)) > tiny) {
#line 407 "cstedc.f"
		    ++finish;
#line 408 "cstedc.f"
		    goto L40;
#line 409 "cstedc.f"
		}
#line 410 "cstedc.f"
	    }

/*           (Sub) Problem determined.  Compute its size and solve it. */

#line 414 "cstedc.f"
	    m = finish - start + 1;
#line 415 "cstedc.f"
	    if (m > smlsiz) {

/*              Scale. */

#line 419 "cstedc.f"
		orgnrm = slanst_("M", &m, &d__[start], &e[start], (ftnlen)1);
#line 420 "cstedc.f"
		slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);
#line 422 "cstedc.f"
		i__1 = m - 1;
#line 422 "cstedc.f"
		i__2 = m - 1;
#line 422 "cstedc.f"
		slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[
			start], &i__2, info, (ftnlen)1);

#line 425 "cstedc.f"
		claed0_(n, &m, &d__[start], &e[start], &z__[start * z_dim1 + 
			1], ldz, &work[1], n, &rwork[1], &iwork[1], info);
#line 427 "cstedc.f"
		if (*info > 0) {
#line 428 "cstedc.f"
		    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info %
			     (m + 1) + start - 1;
#line 430 "cstedc.f"
		    goto L70;
#line 431 "cstedc.f"
		}

/*              Scale back. */

#line 435 "cstedc.f"
		slascl_("G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);

#line 438 "cstedc.f"
	    } else {
#line 439 "cstedc.f"
		ssteqr_("I", &m, &d__[start], &e[start], &rwork[1], &m, &
			rwork[m * m + 1], info, (ftnlen)1);
#line 441 "cstedc.f"
		clacrm_(n, &m, &z__[start * z_dim1 + 1], ldz, &rwork[1], &m, &
			work[1], n, &rwork[m * m + 1]);
#line 443 "cstedc.f"
		clacpy_("A", n, &m, &work[1], n, &z__[start * z_dim1 + 1], 
			ldz, (ftnlen)1);
#line 444 "cstedc.f"
		if (*info > 0) {
#line 445 "cstedc.f"
		    *info = start * (*n + 1) + finish;
#line 446 "cstedc.f"
		    goto L70;
#line 447 "cstedc.f"
		}
#line 448 "cstedc.f"
	    }

#line 450 "cstedc.f"
	    start = finish + 1;
#line 451 "cstedc.f"
	    goto L30;
#line 452 "cstedc.f"
	}

/*        endwhile */


/*        Use Selection Sort to minimize swaps of eigenvectors */

#line 459 "cstedc.f"
	i__1 = *n;
#line 459 "cstedc.f"
	for (ii = 2; ii <= i__1; ++ii) {
#line 460 "cstedc.f"
	    i__ = ii - 1;
#line 461 "cstedc.f"
	    k = i__;
#line 462 "cstedc.f"
	    p = d__[i__];
#line 463 "cstedc.f"
	    i__2 = *n;
#line 463 "cstedc.f"
	    for (j = ii; j <= i__2; ++j) {
#line 464 "cstedc.f"
		if (d__[j] < p) {
#line 465 "cstedc.f"
		    k = j;
#line 466 "cstedc.f"
		    p = d__[j];
#line 467 "cstedc.f"
		}
#line 468 "cstedc.f"
/* L50: */
#line 468 "cstedc.f"
	    }
#line 469 "cstedc.f"
	    if (k != i__) {
#line 470 "cstedc.f"
		d__[k] = d__[i__];
#line 471 "cstedc.f"
		d__[i__] = p;
#line 472 "cstedc.f"
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
#line 473 "cstedc.f"
	    }
#line 474 "cstedc.f"
/* L60: */
#line 474 "cstedc.f"
	}
#line 475 "cstedc.f"
    }

#line 477 "cstedc.f"
L70:
#line 478 "cstedc.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 479 "cstedc.f"
    rwork[1] = (doublereal) lrwmin;
#line 480 "cstedc.f"
    iwork[1] = liwmin;

#line 482 "cstedc.f"
    return 0;

/*     End of CSTEDC */

} /* cstedc_ */

