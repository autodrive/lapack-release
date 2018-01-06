#line 1 "sstedc.f"
/* sstedc.f -- translated by f2c (version 20100827).
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

#line 1 "sstedc.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;

/* > \brief \b SSTEBZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTEDC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstedc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstedc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstedc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
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
/* > SSTEDC computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > The eigenvectors of a full or band real symmetric matrix can also be */
/* > found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this */
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
/* >          = 'V':  Compute eigenvectors of original dense symmetric */
/* >                  matrix also.  On entry, Z contains the orthogonal */
/* >                  matrix used to reduce the original matrix to */
/* >                  tridiagonal form. */
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
/* >          Z is REAL array, dimension (LDZ,N) */
/* >          On entry, if COMPZ = 'V', then Z contains the orthogonal */
/* >          matrix used in the reduction to tridiagonal form. */
/* >          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the */
/* >          orthonormal eigenvectors of the original symmetric matrix, */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1. */
/* >          If COMPZ = 'V' and N > 1 then LWORK must be at least */
/* >                         ( 1 + 3*N + 2*N*lg N + 4*N**2 ), */
/* >                         where lg( N ) = smallest integer k such */
/* >                         that 2**k >= N. */
/* >          If COMPZ = 'I' and N > 1 then LWORK must be at least */
/* >                         ( 1 + 4*N + N**2 ). */
/* >          Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* >          equal to the minimum divide size, usually 25, then LWORK need */
/* >          only be max(1,2*(N-1)). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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
/* >          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1. */
/* >          If COMPZ = 'V' and N > 1 then LIWORK must be at least */
/* >                         ( 6 + 6*N + 5*N*lg N ). */
/* >          If COMPZ = 'I' and N > 1 then LIWORK must be at least */
/* >                         ( 3 + 5*N ). */
/* >          Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* >          equal to the minimum divide size, usually 25, then LIWORK */
/* >          need only be 1. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of the IWORK array, */
/* >          returns this value as the first entry of the IWORK array, and */
/* >          no error message related to LIWORK is issued by XERBLA. */
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

/* > \date November 2011 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* >  Modified by Francoise Tisseur, University of Tennessee */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sstedc_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen compz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal p;
    static integer ii, lgn;
    static doublereal eps, tiny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer lwmin, start;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slaed0_(integer *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, doublereal *
	    , integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer finish;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), slacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static integer liwmin, icompz;
    static doublereal orgnrm;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int ssterf_(integer *, doublereal *, doublereal *,
	     integer *), slasrt_(char *, integer *, doublereal *, integer *, 
	    ftnlen);
    static logical lquery;
    static integer smlsiz;
    extern /* Subroutine */ int ssteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer storez, strtrw;


/*  -- LAPACK computational routine (version 3.4.0) -- */
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

#line 234 "sstedc.f"
    /* Parameter adjustments */
#line 234 "sstedc.f"
    --d__;
#line 234 "sstedc.f"
    --e;
#line 234 "sstedc.f"
    z_dim1 = *ldz;
#line 234 "sstedc.f"
    z_offset = 1 + z_dim1;
#line 234 "sstedc.f"
    z__ -= z_offset;
#line 234 "sstedc.f"
    --work;
#line 234 "sstedc.f"
    --iwork;
#line 234 "sstedc.f"

#line 234 "sstedc.f"
    /* Function Body */
#line 234 "sstedc.f"
    *info = 0;
#line 235 "sstedc.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 237 "sstedc.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "sstedc.f"
	icompz = 0;
#line 239 "sstedc.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 240 "sstedc.f"
	icompz = 1;
#line 241 "sstedc.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 242 "sstedc.f"
	icompz = 2;
#line 243 "sstedc.f"
    } else {
#line 244 "sstedc.f"
	icompz = -1;
#line 245 "sstedc.f"
    }
#line 246 "sstedc.f"
    if (icompz < 0) {
#line 247 "sstedc.f"
	*info = -1;
#line 248 "sstedc.f"
    } else if (*n < 0) {
#line 249 "sstedc.f"
	*info = -2;
#line 250 "sstedc.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 252 "sstedc.f"
	*info = -6;
#line 253 "sstedc.f"
    }

#line 255 "sstedc.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 259 "sstedc.f"
	smlsiz = ilaenv_(&c__9, "SSTEDC", " ", &c__0, &c__0, &c__0, &c__0, (
		ftnlen)6, (ftnlen)1);
#line 260 "sstedc.f"
	if (*n <= 1 || icompz == 0) {
#line 261 "sstedc.f"
	    liwmin = 1;
#line 262 "sstedc.f"
	    lwmin = 1;
#line 263 "sstedc.f"
	} else if (*n <= smlsiz) {
#line 264 "sstedc.f"
	    liwmin = 1;
#line 265 "sstedc.f"
	    lwmin = *n - 1 << 1;
#line 266 "sstedc.f"
	} else {
#line 267 "sstedc.f"
	    lgn = (integer) (log((doublereal) (*n)) / log(2.));
#line 268 "sstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 268 "sstedc.f"
		++lgn;
#line 268 "sstedc.f"
	    }
#line 270 "sstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 270 "sstedc.f"
		++lgn;
#line 270 "sstedc.f"
	    }
#line 272 "sstedc.f"
	    if (icompz == 1) {
/* Computing 2nd power */
#line 273 "sstedc.f"
		i__1 = *n;
#line 273 "sstedc.f"
		lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
#line 274 "sstedc.f"
		liwmin = *n * 6 + 6 + *n * 5 * lgn;
#line 275 "sstedc.f"
	    } else if (icompz == 2) {
/* Computing 2nd power */
#line 276 "sstedc.f"
		i__1 = *n;
#line 276 "sstedc.f"
		lwmin = (*n << 2) + 1 + i__1 * i__1;
#line 277 "sstedc.f"
		liwmin = *n * 5 + 3;
#line 278 "sstedc.f"
	    }
#line 279 "sstedc.f"
	}
#line 280 "sstedc.f"
	work[1] = (doublereal) lwmin;
#line 281 "sstedc.f"
	iwork[1] = liwmin;

#line 283 "sstedc.f"
	if (*lwork < lwmin && ! lquery) {
#line 284 "sstedc.f"
	    *info = -8;
#line 285 "sstedc.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 286 "sstedc.f"
	    *info = -10;
#line 287 "sstedc.f"
	}
#line 288 "sstedc.f"
    }

#line 290 "sstedc.f"
    if (*info != 0) {
#line 291 "sstedc.f"
	i__1 = -(*info);
#line 291 "sstedc.f"
	xerbla_("SSTEDC", &i__1, (ftnlen)6);
#line 292 "sstedc.f"
	return 0;
#line 293 "sstedc.f"
    } else if (lquery) {
#line 294 "sstedc.f"
	return 0;
#line 295 "sstedc.f"
    }

/*     Quick return if possible */

#line 299 "sstedc.f"
    if (*n == 0) {
#line 299 "sstedc.f"
	return 0;
#line 299 "sstedc.f"
    }
#line 301 "sstedc.f"
    if (*n == 1) {
#line 302 "sstedc.f"
	if (icompz != 0) {
#line 302 "sstedc.f"
	    z__[z_dim1 + 1] = 1.;
#line 302 "sstedc.f"
	}
#line 304 "sstedc.f"
	return 0;
#line 305 "sstedc.f"
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

#line 318 "sstedc.f"
    if (icompz == 0) {
#line 319 "sstedc.f"
	ssterf_(n, &d__[1], &e[1], info);
#line 320 "sstedc.f"
	goto L50;
#line 321 "sstedc.f"
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then */
/*     solve the problem with another solver. */

#line 326 "sstedc.f"
    if (*n <= smlsiz) {

#line 328 "sstedc.f"
	ssteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info,
		 (ftnlen)1);

#line 330 "sstedc.f"
    } else {

/*        If COMPZ = 'V', the Z matrix must be stored elsewhere for later */
/*        use. */

#line 335 "sstedc.f"
	if (icompz == 1) {
#line 336 "sstedc.f"
	    storez = *n * *n + 1;
#line 337 "sstedc.f"
	} else {
#line 338 "sstedc.f"
	    storez = 1;
#line 339 "sstedc.f"
	}

#line 341 "sstedc.f"
	if (icompz == 2) {
#line 342 "sstedc.f"
	    slaset_("Full", n, n, &c_b17, &c_b18, &z__[z_offset], ldz, (
		    ftnlen)4);
#line 343 "sstedc.f"
	}

/*        Scale. */

#line 347 "sstedc.f"
	orgnrm = slanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 348 "sstedc.f"
	if (orgnrm == 0.) {
#line 348 "sstedc.f"
	    goto L50;
#line 348 "sstedc.f"
	}

#line 351 "sstedc.f"
	eps = slamch_("Epsilon", (ftnlen)7);

#line 353 "sstedc.f"
	start = 1;

/*        while ( START <= N ) */

#line 357 "sstedc.f"
L10:
#line 358 "sstedc.f"
	if (start <= *n) {

/*           Let FINISH be the position of the next subdiagonal entry */
/*           such that E( FINISH ) <= TINY or FINISH = N if no such */
/*           subdiagonal exists.  The matrix identified by the elements */
/*           between START and FINISH constitutes an independent */
/*           sub-problem. */

#line 366 "sstedc.f"
	    finish = start;
#line 367 "sstedc.f"
L20:
#line 368 "sstedc.f"
	    if (finish < *n) {
#line 369 "sstedc.f"
		tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) * sqrt((
			d__2 = d__[finish + 1], abs(d__2)));
#line 371 "sstedc.f"
		if ((d__1 = e[finish], abs(d__1)) > tiny) {
#line 372 "sstedc.f"
		    ++finish;
#line 373 "sstedc.f"
		    goto L20;
#line 374 "sstedc.f"
		}
#line 375 "sstedc.f"
	    }

/*           (Sub) Problem determined.  Compute its size and solve it. */

#line 379 "sstedc.f"
	    m = finish - start + 1;
#line 380 "sstedc.f"
	    if (m == 1) {
#line 381 "sstedc.f"
		start = finish + 1;
#line 382 "sstedc.f"
		goto L10;
#line 383 "sstedc.f"
	    }
#line 384 "sstedc.f"
	    if (m > smlsiz) {

/*              Scale. */

#line 388 "sstedc.f"
		orgnrm = slanst_("M", &m, &d__[start], &e[start], (ftnlen)1);
#line 389 "sstedc.f"
		slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);
#line 391 "sstedc.f"
		i__1 = m - 1;
#line 391 "sstedc.f"
		i__2 = m - 1;
#line 391 "sstedc.f"
		slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[
			start], &i__2, info, (ftnlen)1);

#line 394 "sstedc.f"
		if (icompz == 1) {
#line 395 "sstedc.f"
		    strtrw = 1;
#line 396 "sstedc.f"
		} else {
#line 397 "sstedc.f"
		    strtrw = start;
#line 398 "sstedc.f"
		}
#line 399 "sstedc.f"
		slaed0_(&icompz, n, &m, &d__[start], &e[start], &z__[strtrw + 
			start * z_dim1], ldz, &work[1], n, &work[storez], &
			iwork[1], info);
#line 402 "sstedc.f"
		if (*info != 0) {
#line 403 "sstedc.f"
		    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info %
			     (m + 1) + start - 1;
#line 405 "sstedc.f"
		    goto L50;
#line 406 "sstedc.f"
		}

/*              Scale back. */

#line 410 "sstedc.f"
		slascl_("G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);

#line 413 "sstedc.f"
	    } else {
#line 414 "sstedc.f"
		if (icompz == 1) {

/*                 Since QR won't update a Z matrix which is larger than */
/*                 the length of D, we must solve the sub-problem in a */
/*                 workspace and then multiply back into Z. */

#line 420 "sstedc.f"
		    ssteqr_("I", &m, &d__[start], &e[start], &work[1], &m, &
			    work[m * m + 1], info, (ftnlen)1);
#line 422 "sstedc.f"
		    slacpy_("A", n, &m, &z__[start * z_dim1 + 1], ldz, &work[
			    storez], n, (ftnlen)1);
#line 424 "sstedc.f"
		    sgemm_("N", "N", n, &m, &m, &c_b18, &work[storez], n, &
			    work[1], &m, &c_b17, &z__[start * z_dim1 + 1], 
			    ldz, (ftnlen)1, (ftnlen)1);
#line 427 "sstedc.f"
		} else if (icompz == 2) {
#line 428 "sstedc.f"
		    ssteqr_("I", &m, &d__[start], &e[start], &z__[start + 
			    start * z_dim1], ldz, &work[1], info, (ftnlen)1);
#line 430 "sstedc.f"
		} else {
#line 431 "sstedc.f"
		    ssterf_(&m, &d__[start], &e[start], info);
#line 432 "sstedc.f"
		}
#line 433 "sstedc.f"
		if (*info != 0) {
#line 434 "sstedc.f"
		    *info = start * (*n + 1) + finish;
#line 435 "sstedc.f"
		    goto L50;
#line 436 "sstedc.f"
		}
#line 437 "sstedc.f"
	    }

#line 439 "sstedc.f"
	    start = finish + 1;
#line 440 "sstedc.f"
	    goto L10;
#line 441 "sstedc.f"
	}

/*        endwhile */

/*        If the problem split any number of times, then the eigenvalues */
/*        will not be properly ordered.  Here we permute the eigenvalues */
/*        (and the associated eigenvectors) into ascending order. */

#line 449 "sstedc.f"
	if (m != *n) {
#line 450 "sstedc.f"
	    if (icompz == 0) {

/*              Use Quick Sort */

#line 454 "sstedc.f"
		slasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 456 "sstedc.f"
	    } else {

/*              Use Selection Sort to minimize swaps of eigenvectors */

#line 460 "sstedc.f"
		i__1 = *n;
#line 460 "sstedc.f"
		for (ii = 2; ii <= i__1; ++ii) {
#line 461 "sstedc.f"
		    i__ = ii - 1;
#line 462 "sstedc.f"
		    k = i__;
#line 463 "sstedc.f"
		    p = d__[i__];
#line 464 "sstedc.f"
		    i__2 = *n;
#line 464 "sstedc.f"
		    for (j = ii; j <= i__2; ++j) {
#line 465 "sstedc.f"
			if (d__[j] < p) {
#line 466 "sstedc.f"
			    k = j;
#line 467 "sstedc.f"
			    p = d__[j];
#line 468 "sstedc.f"
			}
#line 469 "sstedc.f"
/* L30: */
#line 469 "sstedc.f"
		    }
#line 470 "sstedc.f"
		    if (k != i__) {
#line 471 "sstedc.f"
			d__[k] = d__[i__];
#line 472 "sstedc.f"
			d__[i__] = p;
#line 473 "sstedc.f"
			sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * 
				z_dim1 + 1], &c__1);
#line 474 "sstedc.f"
		    }
#line 475 "sstedc.f"
/* L40: */
#line 475 "sstedc.f"
		}
#line 476 "sstedc.f"
	    }
#line 477 "sstedc.f"
	}
#line 478 "sstedc.f"
    }

#line 480 "sstedc.f"
L50:
#line 481 "sstedc.f"
    work[1] = (doublereal) lwmin;
#line 482 "sstedc.f"
    iwork[1] = liwmin;

#line 484 "sstedc.f"
    return 0;

/*     End of SSTEDC */

} /* sstedc_ */

