#line 1 "ssbgvd.f"
/* ssbgvd.f -- translated by f2c (version 20100827).
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

#line 1 "ssbgvd.f"
/* Table of constant values */

static doublereal c_b12 = 1.;
static doublereal c_b13 = 0.;

/* > \brief \b SSBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSBGVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbgvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbgvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbgvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, */
/*                          Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AB( LDAB, * ), BB( LDBB, * ), W( * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBGVD computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite banded eigenproblem, of the */
/* > form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and */
/* > banded, and B is also positive definite.  If eigenvectors are */
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
/* >          = 'U':  Upper triangles of A and B are stored; */
/* >          = 'L':  Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KA */
/* > \verbatim */
/* >          KA is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of superdiagonals of the matrix B if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first ka+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). */
/* > */
/* >          On exit, the contents of AB are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BB */
/* > \verbatim */
/* >          BB is REAL array, dimension (LDBB, N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix B, stored in the first kb+1 rows of the array.  The */
/* >          j-th column of B is stored in the j-th column of the array BB */
/* >          as follows: */
/* >          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j; */
/* >          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb). */
/* > */
/* >          On exit, the factor S from the split Cholesky factorization */
/* >          B = S**T*S, as returned by SPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* >          LDBB is INTEGER */
/* >          The leading dimension of the array BB.  LDBB >= KB+1. */
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
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors, with the i-th column of Z holding the */
/* >          eigenvector associated with W(i).  The eigenvectors are */
/* >          normalized so Z**T*B*Z = I. */
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
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK >= 3*N. */
/* >          If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2. */
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
/* >          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          If JOBZ  = 'N' or N <= 1, LIWORK >= 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. */
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
/* >          > 0:  if INFO = i, and i is: */
/* >             <= N:  the algorithm failed to converge: */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then SPBSTF */
/* >                    returned INFO = i: B is not positive definite. */
/* >                    The factorization of B could not be completed and */
/* >                    no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int ssbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;

    /* Local variables */
    static integer inde;
    static char vect[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer lwmin;
    static logical upper, wantz;
    static integer indwk2, llwrk2;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), sstedc_(
	    char *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen), slacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    static integer indwrk, liwmin;
    extern /* Subroutine */ int spbstf_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), ssbtrd_(char *, char 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), ssbgst_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen),
	     ssterf_(integer *, doublereal *, doublereal *, integer *);
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 269 "ssbgvd.f"
    /* Parameter adjustments */
#line 269 "ssbgvd.f"
    ab_dim1 = *ldab;
#line 269 "ssbgvd.f"
    ab_offset = 1 + ab_dim1;
#line 269 "ssbgvd.f"
    ab -= ab_offset;
#line 269 "ssbgvd.f"
    bb_dim1 = *ldbb;
#line 269 "ssbgvd.f"
    bb_offset = 1 + bb_dim1;
#line 269 "ssbgvd.f"
    bb -= bb_offset;
#line 269 "ssbgvd.f"
    --w;
#line 269 "ssbgvd.f"
    z_dim1 = *ldz;
#line 269 "ssbgvd.f"
    z_offset = 1 + z_dim1;
#line 269 "ssbgvd.f"
    z__ -= z_offset;
#line 269 "ssbgvd.f"
    --work;
#line 269 "ssbgvd.f"
    --iwork;
#line 269 "ssbgvd.f"

#line 269 "ssbgvd.f"
    /* Function Body */
#line 269 "ssbgvd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 270 "ssbgvd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 271 "ssbgvd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 273 "ssbgvd.f"
    *info = 0;
#line 274 "ssbgvd.f"
    if (*n <= 1) {
#line 275 "ssbgvd.f"
	liwmin = 1;
#line 276 "ssbgvd.f"
	lwmin = 1;
#line 277 "ssbgvd.f"
    } else if (wantz) {
#line 278 "ssbgvd.f"
	liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 279 "ssbgvd.f"
	i__1 = *n;
#line 279 "ssbgvd.f"
	lwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 280 "ssbgvd.f"
    } else {
#line 281 "ssbgvd.f"
	liwmin = 1;
#line 282 "ssbgvd.f"
	lwmin = *n << 1;
#line 283 "ssbgvd.f"
    }

#line 285 "ssbgvd.f"
    if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 286 "ssbgvd.f"
	*info = -1;
#line 287 "ssbgvd.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 288 "ssbgvd.f"
	*info = -2;
#line 289 "ssbgvd.f"
    } else if (*n < 0) {
#line 290 "ssbgvd.f"
	*info = -3;
#line 291 "ssbgvd.f"
    } else if (*ka < 0) {
#line 292 "ssbgvd.f"
	*info = -4;
#line 293 "ssbgvd.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 294 "ssbgvd.f"
	*info = -5;
#line 295 "ssbgvd.f"
    } else if (*ldab < *ka + 1) {
#line 296 "ssbgvd.f"
	*info = -7;
#line 297 "ssbgvd.f"
    } else if (*ldbb < *kb + 1) {
#line 298 "ssbgvd.f"
	*info = -9;
#line 299 "ssbgvd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 300 "ssbgvd.f"
	*info = -12;
#line 301 "ssbgvd.f"
    }

#line 303 "ssbgvd.f"
    if (*info == 0) {
#line 304 "ssbgvd.f"
	work[1] = (doublereal) lwmin;
#line 305 "ssbgvd.f"
	iwork[1] = liwmin;

#line 307 "ssbgvd.f"
	if (*lwork < lwmin && ! lquery) {
#line 308 "ssbgvd.f"
	    *info = -14;
#line 309 "ssbgvd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 310 "ssbgvd.f"
	    *info = -16;
#line 311 "ssbgvd.f"
	}
#line 312 "ssbgvd.f"
    }

#line 314 "ssbgvd.f"
    if (*info != 0) {
#line 315 "ssbgvd.f"
	i__1 = -(*info);
#line 315 "ssbgvd.f"
	xerbla_("SSBGVD", &i__1, (ftnlen)6);
#line 316 "ssbgvd.f"
	return 0;
#line 317 "ssbgvd.f"
    } else if (lquery) {
#line 318 "ssbgvd.f"
	return 0;
#line 319 "ssbgvd.f"
    }

/*     Quick return if possible */

#line 323 "ssbgvd.f"
    if (*n == 0) {
#line 323 "ssbgvd.f"
	return 0;
#line 323 "ssbgvd.f"
    }

/*     Form a split Cholesky factorization of B. */

#line 328 "ssbgvd.f"
    spbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info, (ftnlen)1);
#line 329 "ssbgvd.f"
    if (*info != 0) {
#line 330 "ssbgvd.f"
	*info = *n + *info;
#line 331 "ssbgvd.f"
	return 0;
#line 332 "ssbgvd.f"
    }

/*     Transform problem to standard eigenvalue problem. */

#line 336 "ssbgvd.f"
    inde = 1;
#line 337 "ssbgvd.f"
    indwrk = inde + *n;
#line 338 "ssbgvd.f"
    indwk2 = indwrk + *n * *n;
#line 339 "ssbgvd.f"
    llwrk2 = *lwork - indwk2 + 1;
#line 340 "ssbgvd.f"
    ssbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb,
	     &z__[z_offset], ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1)
	    ;

/*     Reduce to tridiagonal form. */

#line 345 "ssbgvd.f"
    if (wantz) {
#line 346 "ssbgvd.f"
	*(unsigned char *)vect = 'U';
#line 347 "ssbgvd.f"
    } else {
#line 348 "ssbgvd.f"
	*(unsigned char *)vect = 'N';
#line 349 "ssbgvd.f"
    }
#line 350 "ssbgvd.f"
    ssbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[
	    z_offset], ldz, &work[indwrk], &iinfo, (ftnlen)1, (ftnlen)1);

/*     For eigenvalues only, call SSTERF. For eigenvectors, call SSTEDC. */

#line 355 "ssbgvd.f"
    if (! wantz) {
#line 356 "ssbgvd.f"
	ssterf_(n, &w[1], &work[inde], info);
#line 357 "ssbgvd.f"
    } else {
#line 358 "ssbgvd.f"
	sstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &
		llwrk2, &iwork[1], liwork, info, (ftnlen)1);
#line 360 "ssbgvd.f"
	sgemm_("N", "N", n, n, n, &c_b12, &z__[z_offset], ldz, &work[indwrk], 
		n, &c_b13, &work[indwk2], n, (ftnlen)1, (ftnlen)1);
#line 362 "ssbgvd.f"
	slacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz, (ftnlen)1);
#line 363 "ssbgvd.f"
    }

#line 365 "ssbgvd.f"
    work[1] = (doublereal) lwmin;
#line 366 "ssbgvd.f"
    iwork[1] = liwmin;

#line 368 "ssbgvd.f"
    return 0;

/*     End of SSBGVD */

} /* ssbgvd_ */

