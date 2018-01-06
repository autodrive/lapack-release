#line 1 "zhegvd.f"
/* zhegvd.f -- translated by f2c (version 20100827).
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

#line 1 "zhegvd.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZHEGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEGVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, */
/*                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEGVD computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and */
/* > B are assumed to be Hermitian and B is also positive definite. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* > */
/* >          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* >          matrix Z of eigenvectors.  The eigenvectors are normalized */
/* >          as follows: */
/* >          if ITYPE = 1 or 2, Z**H*B*Z = I; */
/* >          if ITYPE = 3, Z**H*inv(B)*Z = I. */
/* >          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U') */
/* >          or the lower triangle (if UPLO='L') of A, including the */
/* >          diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
/* >          On entry, the Hermitian matrix B.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of B contains the */
/* >          upper triangular part of the matrix B.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of B contains */
/* >          the lower triangular part of the matrix B. */
/* > */
/* >          On exit, if INFO <= N, the part of B containing the matrix is */
/* >          overwritten by the triangular factor U or L from the Cholesky */
/* >          factorization B = U**H*U or B = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
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
/* >          The length of the array WORK. */
/* >          If N <= 1,                LWORK >= 1. */
/* >          If JOBZ  = 'N' and N > 1, LWORK >= N + 1. */
/* >          If JOBZ  = 'V' and N > 1, LWORK >= 2*N + N**2. */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of the array RWORK. */
/* >          If N <= 1,                LRWORK >= 1. */
/* >          If JOBZ  = 'N' and N > 1, LRWORK >= N. */
/* >          If JOBZ  = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2. */
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
/* >          If N <= 1,                LIWORK >= 1. */
/* >          If JOBZ  = 'N' and N > 1, LIWORK >= 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  ZPOTRF or ZHEEVD returned an error code: */
/* >             <= N:  if INFO = i and JOBZ = 'N', then the algorithm */
/* >                    failed to converge; i off-diagonal elements of an */
/* >                    intermediate tridiagonal form did not converge to */
/* >                    zero; */
/* >                    if INFO = i and JOBZ = 'V', then the algorithm */
/* >                    failed to compute an eigenvalue while working on */
/* >                    the submatrix lying in rows and columns INFO/(N+1) */
/* >                    through mod(INFO,N+1); */
/* >             > N:   if INFO = N + i, for 1 <= i <= N, then the leading */
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

/* > \date November 2011 */

/* > \ingroup complex16HEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified so that no backsubstitution is performed if ZHEEVD fails to */
/* >  converge (NEIG in old code could be greater than N causing out of */
/* >  bounds reference to A - reported by Ralf Meyer).  Also corrected the */
/* >  description of INFO and the test on ITYPE. Sven, 16 Feb 05. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhegvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
	 integer *lrwork, integer *iwork, integer *liwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer lopt;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lwmin;
    static char trans[1];
    static integer liopt;
    static logical upper;
    static integer lropt;
    static logical wantz;
    extern /* Subroutine */ int ztrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    ztrsm_(char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), zheevd_(char *, char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen);
    static integer liwmin;
    extern /* Subroutine */ int zhegst_(integer *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);
    static integer lrwmin;
    static logical lquery;
    extern /* Subroutine */ int zpotrf_(char *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen);


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

#line 292 "zhegvd.f"
    /* Parameter adjustments */
#line 292 "zhegvd.f"
    a_dim1 = *lda;
#line 292 "zhegvd.f"
    a_offset = 1 + a_dim1;
#line 292 "zhegvd.f"
    a -= a_offset;
#line 292 "zhegvd.f"
    b_dim1 = *ldb;
#line 292 "zhegvd.f"
    b_offset = 1 + b_dim1;
#line 292 "zhegvd.f"
    b -= b_offset;
#line 292 "zhegvd.f"
    --w;
#line 292 "zhegvd.f"
    --work;
#line 292 "zhegvd.f"
    --rwork;
#line 292 "zhegvd.f"
    --iwork;
#line 292 "zhegvd.f"

#line 292 "zhegvd.f"
    /* Function Body */
#line 292 "zhegvd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 293 "zhegvd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 294 "zhegvd.f"
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

#line 296 "zhegvd.f"
    *info = 0;
#line 297 "zhegvd.f"
    if (*n <= 1) {
#line 298 "zhegvd.f"
	lwmin = 1;
#line 299 "zhegvd.f"
	lrwmin = 1;
#line 300 "zhegvd.f"
	liwmin = 1;
#line 301 "zhegvd.f"
    } else if (wantz) {
#line 302 "zhegvd.f"
	lwmin = (*n << 1) + *n * *n;
#line 303 "zhegvd.f"
	lrwmin = *n * 5 + 1 + (*n << 1) * *n;
#line 304 "zhegvd.f"
	liwmin = *n * 5 + 3;
#line 305 "zhegvd.f"
    } else {
#line 306 "zhegvd.f"
	lwmin = *n + 1;
#line 307 "zhegvd.f"
	lrwmin = *n;
#line 308 "zhegvd.f"
	liwmin = 1;
#line 309 "zhegvd.f"
    }
#line 310 "zhegvd.f"
    lopt = lwmin;
#line 311 "zhegvd.f"
    lropt = lrwmin;
#line 312 "zhegvd.f"
    liopt = liwmin;
#line 313 "zhegvd.f"
    if (*itype < 1 || *itype > 3) {
#line 314 "zhegvd.f"
	*info = -1;
#line 315 "zhegvd.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 316 "zhegvd.f"
	*info = -2;
#line 317 "zhegvd.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 318 "zhegvd.f"
	*info = -3;
#line 319 "zhegvd.f"
    } else if (*n < 0) {
#line 320 "zhegvd.f"
	*info = -4;
#line 321 "zhegvd.f"
    } else if (*lda < max(1,*n)) {
#line 322 "zhegvd.f"
	*info = -6;
#line 323 "zhegvd.f"
    } else if (*ldb < max(1,*n)) {
#line 324 "zhegvd.f"
	*info = -8;
#line 325 "zhegvd.f"
    }

#line 327 "zhegvd.f"
    if (*info == 0) {
#line 328 "zhegvd.f"
	work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 329 "zhegvd.f"
	rwork[1] = (doublereal) lropt;
#line 330 "zhegvd.f"
	iwork[1] = liopt;

#line 332 "zhegvd.f"
	if (*lwork < lwmin && ! lquery) {
#line 333 "zhegvd.f"
	    *info = -11;
#line 334 "zhegvd.f"
	} else if (*lrwork < lrwmin && ! lquery) {
#line 335 "zhegvd.f"
	    *info = -13;
#line 336 "zhegvd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 337 "zhegvd.f"
	    *info = -15;
#line 338 "zhegvd.f"
	}
#line 339 "zhegvd.f"
    }

#line 341 "zhegvd.f"
    if (*info != 0) {
#line 342 "zhegvd.f"
	i__1 = -(*info);
#line 342 "zhegvd.f"
	xerbla_("ZHEGVD", &i__1, (ftnlen)6);
#line 343 "zhegvd.f"
	return 0;
#line 344 "zhegvd.f"
    } else if (lquery) {
#line 345 "zhegvd.f"
	return 0;
#line 346 "zhegvd.f"
    }

/*     Quick return if possible */

#line 350 "zhegvd.f"
    if (*n == 0) {
#line 350 "zhegvd.f"
	return 0;
#line 350 "zhegvd.f"
    }

/*     Form a Cholesky factorization of B. */

#line 355 "zhegvd.f"
    zpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 356 "zhegvd.f"
    if (*info != 0) {
#line 357 "zhegvd.f"
	*info = *n + *info;
#line 358 "zhegvd.f"
	return 0;
#line 359 "zhegvd.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 363 "zhegvd.f"
    zhegst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 364 "zhegvd.f"
    zheevd_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, &rwork[
	    1], lrwork, &iwork[1], liwork, info, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 366 "zhegvd.f"
    d__1 = (doublereal) lopt, d__2 = work[1].r;
#line 366 "zhegvd.f"
    lopt = (integer) max(d__1,d__2);
/* Computing MAX */
#line 367 "zhegvd.f"
    d__1 = (doublereal) lropt;
#line 367 "zhegvd.f"
    lropt = (integer) max(d__1,rwork[1]);
/* Computing MAX */
#line 368 "zhegvd.f"
    d__1 = (doublereal) liopt, d__2 = (doublereal) iwork[1];
#line 368 "zhegvd.f"
    liopt = (integer) max(d__1,d__2);

#line 370 "zhegvd.f"
    if (wantz && *info == 0) {

/*        Backtransform eigenvectors to the original problem. */

#line 374 "zhegvd.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y */

#line 379 "zhegvd.f"
	    if (upper) {
#line 380 "zhegvd.f"
		*(unsigned char *)trans = 'N';
#line 381 "zhegvd.f"
	    } else {
#line 382 "zhegvd.f"
		*(unsigned char *)trans = 'C';
#line 383 "zhegvd.f"
	    }

#line 385 "zhegvd.f"
	    ztrsm_("Left", uplo, trans, "Non-unit", n, n, &c_b1, &b[b_offset],
		     ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (ftnlen)1, 
		    (ftnlen)8);

#line 388 "zhegvd.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**H *y */

#line 393 "zhegvd.f"
	    if (upper) {
#line 394 "zhegvd.f"
		*(unsigned char *)trans = 'C';
#line 395 "zhegvd.f"
	    } else {
#line 396 "zhegvd.f"
		*(unsigned char *)trans = 'N';
#line 397 "zhegvd.f"
	    }

#line 399 "zhegvd.f"
	    ztrmm_("Left", uplo, trans, "Non-unit", n, n, &c_b1, &b[b_offset],
		     ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (ftnlen)1, 
		    (ftnlen)8);
#line 401 "zhegvd.f"
	}
#line 402 "zhegvd.f"
    }

#line 404 "zhegvd.f"
    work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 405 "zhegvd.f"
    rwork[1] = (doublereal) lropt;
#line 406 "zhegvd.f"
    iwork[1] = liopt;

#line 408 "zhegvd.f"
    return 0;

/*     End of ZHEGVD */

} /* zhegvd_ */

