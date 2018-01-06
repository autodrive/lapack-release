#line 1 "chpgvd.f"
/* chpgvd.f -- translated by f2c (version 20100827).
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

#line 1 "chpgvd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHPGVD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPGVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpgvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpgvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpgvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, */
/*                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), W( * ) */
/*       COMPLEX            AP( * ), BP( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPGVD computes all the eigenvalues and, optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and */
/* > B are assumed to be Hermitian, stored in packed format, and B is also */
/* > positive definite. */
/* > If eigenvectors are desired, it uses a divide and conquer algorithm. */
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

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          Specifies the problem type to be solved: */
/* >          = 1:  A*x = (lambda)*B*x */
/* >          = 2:  A*B*x = (lambda)*x */
/* >          = 3:  B*A*x = (lambda)*x */
/* > \endverbatim */
/* > */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the contents of AP are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BP */
/* > \verbatim */
/* >          BP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          B, packed columnwise in a linear array.  The j-th column of B */
/* >          is stored in the array BP as follows: */
/* >          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the triangular factor U or L from the Cholesky */
/* >          factorization B = U**H*U or B = L*L**H, in the same storage */
/* >          format as B. */
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
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* >          eigenvectors.  The eigenvectors are normalized as follows: */
/* >          if ITYPE = 1 or 2, Z**H*B*Z = I; */
/* >          if ITYPE = 3, Z**H*inv(B)*Z = I. */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of array WORK. */
/* >          If N <= 1,               LWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK >= N. */
/* >          If JOBZ = 'V' and N > 1, LWORK >= 2*N. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the required sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* >          On exit, if INFO = 0, RWORK(1) returns the required LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of array RWORK. */
/* >          If N <= 1,               LRWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LRWORK >= N. */
/* >          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
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
/* >          The dimension of array IWORK. */
/* >          If JOBZ  = 'N' or N <= 1, LIWORK >= 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the required sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  CPPTRF or CHPEVD returned an error code: */
/* >             <= N:  if INFO = i, CHPEVD failed to converge; */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not convergeto zero; */
/* >             > N:   if INFO = N + i, for 1 <= i <= n, then the leading */
/* >                    minor of order i of B is not positive definite. */
/* >                    The factorization of B could not be completed and */
/* >                    no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complexOTHEReigen */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int chpgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lwmin;
    static char trans[1];
    extern /* Subroutine */ int ctpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical wantz;
    extern /* Subroutine */ int chpevd_(char *, char *, integer *, 
	    doublecomplex *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), chpgst_(integer *, char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), cpptrf_(char *, integer *, 
	    doublecomplex *, integer *, ftnlen);
    static integer liwmin, lrwmin;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 270 "chpgvd.f"
    /* Parameter adjustments */
#line 270 "chpgvd.f"
    --ap;
#line 270 "chpgvd.f"
    --bp;
#line 270 "chpgvd.f"
    --w;
#line 270 "chpgvd.f"
    z_dim1 = *ldz;
#line 270 "chpgvd.f"
    z_offset = 1 + z_dim1;
#line 270 "chpgvd.f"
    z__ -= z_offset;
#line 270 "chpgvd.f"
    --work;
#line 270 "chpgvd.f"
    --rwork;
#line 270 "chpgvd.f"
    --iwork;
#line 270 "chpgvd.f"

#line 270 "chpgvd.f"
    /* Function Body */
#line 270 "chpgvd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 271 "chpgvd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 272 "chpgvd.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 274 "chpgvd.f"
    *info = 0;
#line 275 "chpgvd.f"
    if (*itype < 1 || *itype > 3) {
#line 276 "chpgvd.f"
	*info = -1;
#line 277 "chpgvd.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 278 "chpgvd.f"
	*info = -2;
#line 279 "chpgvd.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 280 "chpgvd.f"
	*info = -3;
#line 281 "chpgvd.f"
    } else if (*n < 0) {
#line 282 "chpgvd.f"
	*info = -4;
#line 283 "chpgvd.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 284 "chpgvd.f"
	*info = -9;
#line 285 "chpgvd.f"
    }

#line 287 "chpgvd.f"
    if (*info == 0) {
#line 288 "chpgvd.f"
	if (*n <= 1) {
#line 289 "chpgvd.f"
	    lwmin = 1;
#line 290 "chpgvd.f"
	    liwmin = 1;
#line 291 "chpgvd.f"
	    lrwmin = 1;
#line 292 "chpgvd.f"
	} else {
#line 293 "chpgvd.f"
	    if (wantz) {
#line 294 "chpgvd.f"
		lwmin = *n << 1;
/* Computing 2nd power */
#line 295 "chpgvd.f"
		i__1 = *n;
#line 295 "chpgvd.f"
		lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
#line 296 "chpgvd.f"
		liwmin = *n * 5 + 3;
#line 297 "chpgvd.f"
	    } else {
#line 298 "chpgvd.f"
		lwmin = *n;
#line 299 "chpgvd.f"
		lrwmin = *n;
#line 300 "chpgvd.f"
		liwmin = 1;
#line 301 "chpgvd.f"
	    }
#line 302 "chpgvd.f"
	}

#line 304 "chpgvd.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 305 "chpgvd.f"
	rwork[1] = (doublereal) lrwmin;
#line 306 "chpgvd.f"
	iwork[1] = liwmin;
#line 307 "chpgvd.f"
	if (*lwork < lwmin && ! lquery) {
#line 308 "chpgvd.f"
	    *info = -11;
#line 309 "chpgvd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 310 "chpgvd.f"
	    *info = -13;
#line 311 "chpgvd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 312 "chpgvd.f"
	    *info = -15;
#line 313 "chpgvd.f"
	}
#line 314 "chpgvd.f"
    }

#line 316 "chpgvd.f"
    if (*info != 0) {
#line 317 "chpgvd.f"
	i__1 = -(*info);
#line 317 "chpgvd.f"
	xerbla_("CHPGVD", &i__1, (ftnlen)6);
#line 318 "chpgvd.f"
	return 0;
#line 319 "chpgvd.f"
    } else if (lquery) {
#line 320 "chpgvd.f"
	return 0;
#line 321 "chpgvd.f"
    }

/*     Quick return if possible */

#line 325 "chpgvd.f"
    if (*n == 0) {
#line 325 "chpgvd.f"
	return 0;
#line 325 "chpgvd.f"
    }

/*     Form a Cholesky factorization of B. */

#line 330 "chpgvd.f"
    cpptrf_(uplo, n, &bp[1], info, (ftnlen)1);
#line 331 "chpgvd.f"
    if (*info != 0) {
#line 332 "chpgvd.f"
	*info = *n + *info;
#line 333 "chpgvd.f"
	return 0;
#line 334 "chpgvd.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 338 "chpgvd.f"
    chpgst_(itype, uplo, n, &ap[1], &bp[1], info, (ftnlen)1);
#line 339 "chpgvd.f"
    chpevd_(jobz, uplo, n, &ap[1], &w[1], &z__[z_offset], ldz, &work[1], 
	    lwork, &rwork[1], lrwork, &iwork[1], liwork, info, (ftnlen)1, (
	    ftnlen)1);
/* Computing MAX */
#line 341 "chpgvd.f"
    d__1 = (doublereal) lwmin, d__2 = work[1].r;
#line 341 "chpgvd.f"
    lwmin = (integer) max(d__1,d__2);
/* Computing MAX */
#line 342 "chpgvd.f"
    d__1 = (doublereal) lrwmin;
#line 342 "chpgvd.f"
    lrwmin = (integer) max(d__1,rwork[1]);
/* Computing MAX */
#line 343 "chpgvd.f"
    d__1 = (doublereal) liwmin, d__2 = (doublereal) iwork[1];
#line 343 "chpgvd.f"
    liwmin = (integer) max(d__1,d__2);

#line 345 "chpgvd.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 349 "chpgvd.f"
	neig = *n;
#line 350 "chpgvd.f"
	if (*info > 0) {
#line 350 "chpgvd.f"
	    neig = *info - 1;
#line 350 "chpgvd.f"
	}
#line 352 "chpgvd.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y */

#line 357 "chpgvd.f"
	    if (upper) {
#line 358 "chpgvd.f"
		*(unsigned char *)trans = 'N';
#line 359 "chpgvd.f"
	    } else {
#line 360 "chpgvd.f"
		*(unsigned char *)trans = 'C';
#line 361 "chpgvd.f"
	    }

#line 363 "chpgvd.f"
	    i__1 = neig;
#line 363 "chpgvd.f"
	    for (j = 1; j <= i__1; ++j) {
#line 364 "chpgvd.f"
		ctpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 366 "chpgvd.f"
/* L10: */
#line 366 "chpgvd.f"
	    }

#line 368 "chpgvd.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H *y */

#line 373 "chpgvd.f"
	    if (upper) {
#line 374 "chpgvd.f"
		*(unsigned char *)trans = 'C';
#line 375 "chpgvd.f"
	    } else {
#line 376 "chpgvd.f"
		*(unsigned char *)trans = 'N';
#line 377 "chpgvd.f"
	    }

#line 379 "chpgvd.f"
	    i__1 = neig;
#line 379 "chpgvd.f"
	    for (j = 1; j <= i__1; ++j) {
#line 380 "chpgvd.f"
		ctpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 382 "chpgvd.f"
/* L20: */
#line 382 "chpgvd.f"
	    }
#line 383 "chpgvd.f"
	}
#line 384 "chpgvd.f"
    }

#line 386 "chpgvd.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 387 "chpgvd.f"
    rwork[1] = (doublereal) lrwmin;
#line 388 "chpgvd.f"
    iwork[1] = liwmin;
#line 389 "chpgvd.f"
    return 0;

/*     End of CHPGVD */

} /* chpgvd_ */

