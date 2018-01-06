#line 1 "dstedc.f"
/* dstedc.f -- translated by f2c (version 20100827).
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

#line 1 "dstedc.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;

/* > \brief \b DSTEDC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEDC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstedc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstedc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstedc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, */
/*                          LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEDC computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > The eigenvectors of a full or band real symmetric matrix can also be */
/* > found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \date June 2017 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* >  Modified by Francoise Tisseur, University of Tennessee */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dstedc_(char *compz, integer *n, doublereal *d__, 
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
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lwmin;
    extern /* Subroutine */ int dlaed0_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static integer start;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer finish;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dlasrt_(char *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer liwmin, icompz;
    extern /* Subroutine */ int dsteqr_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static doublereal orgnrm;
    static logical lquery;
    static integer smlsiz, storez, strtrw;


/*  -- LAPACK computational routine (version 3.7.1) -- */
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

#line 234 "dstedc.f"
    /* Parameter adjustments */
#line 234 "dstedc.f"
    --d__;
#line 234 "dstedc.f"
    --e;
#line 234 "dstedc.f"
    z_dim1 = *ldz;
#line 234 "dstedc.f"
    z_offset = 1 + z_dim1;
#line 234 "dstedc.f"
    z__ -= z_offset;
#line 234 "dstedc.f"
    --work;
#line 234 "dstedc.f"
    --iwork;
#line 234 "dstedc.f"

#line 234 "dstedc.f"
    /* Function Body */
#line 234 "dstedc.f"
    *info = 0;
#line 235 "dstedc.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 237 "dstedc.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 238 "dstedc.f"
	icompz = 0;
#line 239 "dstedc.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 240 "dstedc.f"
	icompz = 1;
#line 241 "dstedc.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 242 "dstedc.f"
	icompz = 2;
#line 243 "dstedc.f"
    } else {
#line 244 "dstedc.f"
	icompz = -1;
#line 245 "dstedc.f"
    }
#line 246 "dstedc.f"
    if (icompz < 0) {
#line 247 "dstedc.f"
	*info = -1;
#line 248 "dstedc.f"
    } else if (*n < 0) {
#line 249 "dstedc.f"
	*info = -2;
#line 250 "dstedc.f"
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
#line 252 "dstedc.f"
	*info = -6;
#line 253 "dstedc.f"
    }

#line 255 "dstedc.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 259 "dstedc.f"
	smlsiz = ilaenv_(&c__9, "DSTEDC", " ", &c__0, &c__0, &c__0, &c__0, (
		ftnlen)6, (ftnlen)1);
#line 260 "dstedc.f"
	if (*n <= 1 || icompz == 0) {
#line 261 "dstedc.f"
	    liwmin = 1;
#line 262 "dstedc.f"
	    lwmin = 1;
#line 263 "dstedc.f"
	} else if (*n <= smlsiz) {
#line 264 "dstedc.f"
	    liwmin = 1;
#line 265 "dstedc.f"
	    lwmin = *n - 1 << 1;
#line 266 "dstedc.f"
	} else {
#line 267 "dstedc.f"
	    lgn = (integer) (log((doublereal) (*n)) / log(2.));
#line 268 "dstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 268 "dstedc.f"
		++lgn;
#line 268 "dstedc.f"
	    }
#line 270 "dstedc.f"
	    if (pow_ii(&c__2, &lgn) < *n) {
#line 270 "dstedc.f"
		++lgn;
#line 270 "dstedc.f"
	    }
#line 272 "dstedc.f"
	    if (icompz == 1) {
/* Computing 2nd power */
#line 273 "dstedc.f"
		i__1 = *n;
#line 273 "dstedc.f"
		lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
#line 274 "dstedc.f"
		liwmin = *n * 6 + 6 + *n * 5 * lgn;
#line 275 "dstedc.f"
	    } else if (icompz == 2) {
/* Computing 2nd power */
#line 276 "dstedc.f"
		i__1 = *n;
#line 276 "dstedc.f"
		lwmin = (*n << 2) + 1 + i__1 * i__1;
#line 277 "dstedc.f"
		liwmin = *n * 5 + 3;
#line 278 "dstedc.f"
	    }
#line 279 "dstedc.f"
	}
#line 280 "dstedc.f"
	work[1] = (doublereal) lwmin;
#line 281 "dstedc.f"
	iwork[1] = liwmin;

#line 283 "dstedc.f"
	if (*lwork < lwmin && ! lquery) {
#line 284 "dstedc.f"
	    *info = -8;
#line 285 "dstedc.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 286 "dstedc.f"
	    *info = -10;
#line 287 "dstedc.f"
	}
#line 288 "dstedc.f"
    }

#line 290 "dstedc.f"
    if (*info != 0) {
#line 291 "dstedc.f"
	i__1 = -(*info);
#line 291 "dstedc.f"
	xerbla_("DSTEDC", &i__1, (ftnlen)6);
#line 292 "dstedc.f"
	return 0;
#line 293 "dstedc.f"
    } else if (lquery) {
#line 294 "dstedc.f"
	return 0;
#line 295 "dstedc.f"
    }

/*     Quick return if possible */

#line 299 "dstedc.f"
    if (*n == 0) {
#line 299 "dstedc.f"
	return 0;
#line 299 "dstedc.f"
    }
#line 301 "dstedc.f"
    if (*n == 1) {
#line 302 "dstedc.f"
	if (icompz != 0) {
#line 302 "dstedc.f"
	    z__[z_dim1 + 1] = 1.;
#line 302 "dstedc.f"
	}
#line 304 "dstedc.f"
	return 0;
#line 305 "dstedc.f"
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

#line 318 "dstedc.f"
    if (icompz == 0) {
#line 319 "dstedc.f"
	dsterf_(n, &d__[1], &e[1], info);
#line 320 "dstedc.f"
	goto L50;
#line 321 "dstedc.f"
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then */
/*     solve the problem with another solver. */

#line 326 "dstedc.f"
    if (*n <= smlsiz) {

#line 328 "dstedc.f"
	dsteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info,
		 (ftnlen)1);

#line 330 "dstedc.f"
    } else {

/*        If COMPZ = 'V', the Z matrix must be stored elsewhere for later */
/*        use. */

#line 335 "dstedc.f"
	if (icompz == 1) {
#line 336 "dstedc.f"
	    storez = *n * *n + 1;
#line 337 "dstedc.f"
	} else {
#line 338 "dstedc.f"
	    storez = 1;
#line 339 "dstedc.f"
	}

#line 341 "dstedc.f"
	if (icompz == 2) {
#line 342 "dstedc.f"
	    dlaset_("Full", n, n, &c_b17, &c_b18, &z__[z_offset], ldz, (
		    ftnlen)4);
#line 343 "dstedc.f"
	}

/*        Scale. */

#line 347 "dstedc.f"
	orgnrm = dlanst_("M", n, &d__[1], &e[1], (ftnlen)1);
#line 348 "dstedc.f"
	if (orgnrm == 0.) {
#line 348 "dstedc.f"
	    goto L50;
#line 348 "dstedc.f"
	}

#line 351 "dstedc.f"
	eps = dlamch_("Epsilon", (ftnlen)7);

#line 353 "dstedc.f"
	start = 1;

/*        while ( START <= N ) */

#line 357 "dstedc.f"
L10:
#line 358 "dstedc.f"
	if (start <= *n) {

/*           Let FINISH be the position of the next subdiagonal entry */
/*           such that E( FINISH ) <= TINY or FINISH = N if no such */
/*           subdiagonal exists.  The matrix identified by the elements */
/*           between START and FINISH constitutes an independent */
/*           sub-problem. */

#line 366 "dstedc.f"
	    finish = start;
#line 367 "dstedc.f"
L20:
#line 368 "dstedc.f"
	    if (finish < *n) {
#line 369 "dstedc.f"
		tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) * sqrt((
			d__2 = d__[finish + 1], abs(d__2)));
#line 371 "dstedc.f"
		if ((d__1 = e[finish], abs(d__1)) > tiny) {
#line 372 "dstedc.f"
		    ++finish;
#line 373 "dstedc.f"
		    goto L20;
#line 374 "dstedc.f"
		}
#line 375 "dstedc.f"
	    }

/*           (Sub) Problem determined.  Compute its size and solve it. */

#line 379 "dstedc.f"
	    m = finish - start + 1;
#line 380 "dstedc.f"
	    if (m == 1) {
#line 381 "dstedc.f"
		start = finish + 1;
#line 382 "dstedc.f"
		goto L10;
#line 383 "dstedc.f"
	    }
#line 384 "dstedc.f"
	    if (m > smlsiz) {

/*              Scale. */

#line 388 "dstedc.f"
		orgnrm = dlanst_("M", &m, &d__[start], &e[start], (ftnlen)1);
#line 389 "dstedc.f"
		dlascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);
#line 391 "dstedc.f"
		i__1 = m - 1;
#line 391 "dstedc.f"
		i__2 = m - 1;
#line 391 "dstedc.f"
		dlascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[
			start], &i__2, info, (ftnlen)1);

#line 394 "dstedc.f"
		if (icompz == 1) {
#line 395 "dstedc.f"
		    strtrw = 1;
#line 396 "dstedc.f"
		} else {
#line 397 "dstedc.f"
		    strtrw = start;
#line 398 "dstedc.f"
		}
#line 399 "dstedc.f"
		dlaed0_(&icompz, n, &m, &d__[start], &e[start], &z__[strtrw + 
			start * z_dim1], ldz, &work[1], n, &work[storez], &
			iwork[1], info);
#line 402 "dstedc.f"
		if (*info != 0) {
#line 403 "dstedc.f"
		    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info %
			     (m + 1) + start - 1;
#line 405 "dstedc.f"
		    goto L50;
#line 406 "dstedc.f"
		}

/*              Scale back. */

#line 410 "dstedc.f"
		dlascl_("G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[
			start], &m, info, (ftnlen)1);

#line 413 "dstedc.f"
	    } else {
#line 414 "dstedc.f"
		if (icompz == 1) {

/*                 Since QR won't update a Z matrix which is larger than */
/*                 the length of D, we must solve the sub-problem in a */
/*                 workspace and then multiply back into Z. */

#line 420 "dstedc.f"
		    dsteqr_("I", &m, &d__[start], &e[start], &work[1], &m, &
			    work[m * m + 1], info, (ftnlen)1);
#line 422 "dstedc.f"
		    dlacpy_("A", n, &m, &z__[start * z_dim1 + 1], ldz, &work[
			    storez], n, (ftnlen)1);
#line 424 "dstedc.f"
		    dgemm_("N", "N", n, &m, &m, &c_b18, &work[storez], n, &
			    work[1], &m, &c_b17, &z__[start * z_dim1 + 1], 
			    ldz, (ftnlen)1, (ftnlen)1);
#line 427 "dstedc.f"
		} else if (icompz == 2) {
#line 428 "dstedc.f"
		    dsteqr_("I", &m, &d__[start], &e[start], &z__[start + 
			    start * z_dim1], ldz, &work[1], info, (ftnlen)1);
#line 430 "dstedc.f"
		} else {
#line 431 "dstedc.f"
		    dsterf_(&m, &d__[start], &e[start], info);
#line 432 "dstedc.f"
		}
#line 433 "dstedc.f"
		if (*info != 0) {
#line 434 "dstedc.f"
		    *info = start * (*n + 1) + finish;
#line 435 "dstedc.f"
		    goto L50;
#line 436 "dstedc.f"
		}
#line 437 "dstedc.f"
	    }

#line 439 "dstedc.f"
	    start = finish + 1;
#line 440 "dstedc.f"
	    goto L10;
#line 441 "dstedc.f"
	}

/*        endwhile */

#line 445 "dstedc.f"
	if (icompz == 0) {

/*          Use Quick Sort */

#line 449 "dstedc.f"
	    dlasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 451 "dstedc.f"
	} else {

/*          Use Selection Sort to minimize swaps of eigenvectors */

#line 455 "dstedc.f"
	    i__1 = *n;
#line 455 "dstedc.f"
	    for (ii = 2; ii <= i__1; ++ii) {
#line 456 "dstedc.f"
		i__ = ii - 1;
#line 457 "dstedc.f"
		k = i__;
#line 458 "dstedc.f"
		p = d__[i__];
#line 459 "dstedc.f"
		i__2 = *n;
#line 459 "dstedc.f"
		for (j = ii; j <= i__2; ++j) {
#line 460 "dstedc.f"
		    if (d__[j] < p) {
#line 461 "dstedc.f"
			k = j;
#line 462 "dstedc.f"
			p = d__[j];
#line 463 "dstedc.f"
		    }
#line 464 "dstedc.f"
/* L30: */
#line 464 "dstedc.f"
		}
#line 465 "dstedc.f"
		if (k != i__) {
#line 466 "dstedc.f"
		    d__[k] = d__[i__];
#line 467 "dstedc.f"
		    d__[i__] = p;
#line 468 "dstedc.f"
		    dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 
			    + 1], &c__1);
#line 469 "dstedc.f"
		}
#line 470 "dstedc.f"
/* L40: */
#line 470 "dstedc.f"
	    }
#line 471 "dstedc.f"
	}
#line 472 "dstedc.f"
    }

#line 474 "dstedc.f"
L50:
#line 475 "dstedc.f"
    work[1] = (doublereal) lwmin;
#line 476 "dstedc.f"
    iwork[1] = liwmin;

#line 478 "dstedc.f"
    return 0;

/*     End of DSTEDC */

} /* dstedc_ */

