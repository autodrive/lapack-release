#line 1 "dsygvd.f"
/* dsygvd.f -- translated by f2c (version 20100827).
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

#line 1 "dsygvd.f"
/* Table of constant values */

static doublereal c_b11 = 1.;

/* > \brief \b DSYGVD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYGVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, */
/*                          LWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYGVD computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and */
/* > B are assumed to be symmetric and B is also positive definite. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* > */
/* >          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* >          matrix Z of eigenvectors.  The eigenvectors are normalized */
/* >          as follows: */
/* >          if ITYPE = 1 or 2, Z**T*B*Z = I; */
/* >          if ITYPE = 3, Z**T*inv(B)*Z = I. */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          On entry, the symmetric matrix B.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of B contains the */
/* >          upper triangular part of the matrix B.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of B contains */
/* >          the lower triangular part of the matrix B. */
/* > */
/* >          On exit, if INFO <= N, the part of B containing the matrix is */
/* >          overwritten by the triangular factor U or L from the Cholesky */
/* >          factorization B = U**T*U or B = L*L**T. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= 1,               LWORK >= 1. */
/* >          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1. */
/* >          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2. */
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
/* >          If N <= 1,                LIWORK >= 1. */
/* >          If JOBZ  = 'N' and N > 1, LIWORK >= 1. */
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
/* >          > 0:  DPOTRF or DSYEVD returned an error code: */
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

/* > \date December 2016 */

/* > \ingroup doubleSYeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified so that no backsubstitution is performed if DSYEVD fails to */
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
/* Subroutine */ int dsygvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer lopt;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer lwmin;
    static char trans[1];
    static integer liopt;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper, wantz;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpotrf_(
	    char *, integer *, doublereal *, integer *, integer *, ftnlen);
    static integer liwmin;
    extern /* Subroutine */ int dsyevd_(char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen), dsygst_(integer *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical lquery;


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

#line 269 "dsygvd.f"
    /* Parameter adjustments */
#line 269 "dsygvd.f"
    a_dim1 = *lda;
#line 269 "dsygvd.f"
    a_offset = 1 + a_dim1;
#line 269 "dsygvd.f"
    a -= a_offset;
#line 269 "dsygvd.f"
    b_dim1 = *ldb;
#line 269 "dsygvd.f"
    b_offset = 1 + b_dim1;
#line 269 "dsygvd.f"
    b -= b_offset;
#line 269 "dsygvd.f"
    --w;
#line 269 "dsygvd.f"
    --work;
#line 269 "dsygvd.f"
    --iwork;
#line 269 "dsygvd.f"

#line 269 "dsygvd.f"
    /* Function Body */
#line 269 "dsygvd.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 270 "dsygvd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 271 "dsygvd.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 273 "dsygvd.f"
    *info = 0;
#line 274 "dsygvd.f"
    if (*n <= 1) {
#line 275 "dsygvd.f"
	liwmin = 1;
#line 276 "dsygvd.f"
	lwmin = 1;
#line 277 "dsygvd.f"
    } else if (wantz) {
#line 278 "dsygvd.f"
	liwmin = *n * 5 + 3;
/* Computing 2nd power */
#line 279 "dsygvd.f"
	i__1 = *n;
#line 279 "dsygvd.f"
	lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
#line 280 "dsygvd.f"
    } else {
#line 281 "dsygvd.f"
	liwmin = 1;
#line 282 "dsygvd.f"
	lwmin = (*n << 1) + 1;
#line 283 "dsygvd.f"
    }
#line 284 "dsygvd.f"
    lopt = lwmin;
#line 285 "dsygvd.f"
    liopt = liwmin;
#line 286 "dsygvd.f"
    if (*itype < 1 || *itype > 3) {
#line 287 "dsygvd.f"
	*info = -1;
#line 288 "dsygvd.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 289 "dsygvd.f"
	*info = -2;
#line 290 "dsygvd.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 291 "dsygvd.f"
	*info = -3;
#line 292 "dsygvd.f"
    } else if (*n < 0) {
#line 293 "dsygvd.f"
	*info = -4;
#line 294 "dsygvd.f"
    } else if (*lda < max(1,*n)) {
#line 295 "dsygvd.f"
	*info = -6;
#line 296 "dsygvd.f"
    } else if (*ldb < max(1,*n)) {
#line 297 "dsygvd.f"
	*info = -8;
#line 298 "dsygvd.f"
    }

#line 300 "dsygvd.f"
    if (*info == 0) {
#line 301 "dsygvd.f"
	work[1] = (doublereal) lopt;
#line 302 "dsygvd.f"
	iwork[1] = liopt;

#line 304 "dsygvd.f"
	if (*lwork < lwmin && ! lquery) {
#line 305 "dsygvd.f"
	    *info = -11;
#line 306 "dsygvd.f"
	} else if (*liwork < liwmin && ! lquery) {
#line 307 "dsygvd.f"
	    *info = -13;
#line 308 "dsygvd.f"
	}
#line 309 "dsygvd.f"
    }

#line 311 "dsygvd.f"
    if (*info != 0) {
#line 312 "dsygvd.f"
	i__1 = -(*info);
#line 312 "dsygvd.f"
	xerbla_("DSYGVD", &i__1, (ftnlen)6);
#line 313 "dsygvd.f"
	return 0;
#line 314 "dsygvd.f"
    } else if (lquery) {
#line 315 "dsygvd.f"
	return 0;
#line 316 "dsygvd.f"
    }

/*     Quick return if possible */

#line 320 "dsygvd.f"
    if (*n == 0) {
#line 320 "dsygvd.f"
	return 0;
#line 320 "dsygvd.f"
    }

/*     Form a Cholesky factorization of B. */

#line 325 "dsygvd.f"
    dpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 326 "dsygvd.f"
    if (*info != 0) {
#line 327 "dsygvd.f"
	*info = *n + *info;
#line 328 "dsygvd.f"
	return 0;
#line 329 "dsygvd.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 333 "dsygvd.f"
    dsygst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 334 "dsygvd.f"
    dsyevd_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, &iwork[
	    1], liwork, info, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 336 "dsygvd.f"
    d__1 = (doublereal) lopt;
#line 336 "dsygvd.f"
    lopt = (integer) max(d__1,work[1]);
/* Computing MAX */
#line 337 "dsygvd.f"
    d__1 = (doublereal) liopt, d__2 = (doublereal) iwork[1];
#line 337 "dsygvd.f"
    liopt = (integer) max(d__1,d__2);

#line 339 "dsygvd.f"
    if (wantz && *info == 0) {

/*        Backtransform eigenvectors to the original problem. */

#line 343 "dsygvd.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 348 "dsygvd.f"
	    if (upper) {
#line 349 "dsygvd.f"
		*(unsigned char *)trans = 'N';
#line 350 "dsygvd.f"
	    } else {
#line 351 "dsygvd.f"
		*(unsigned char *)trans = 'T';
#line 352 "dsygvd.f"
	    }

#line 354 "dsygvd.f"
	    dtrsm_("Left", uplo, trans, "Non-unit", n, n, &c_b11, &b[b_offset]
		    , ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (ftnlen)1,
		     (ftnlen)8);

#line 357 "dsygvd.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 362 "dsygvd.f"
	    if (upper) {
#line 363 "dsygvd.f"
		*(unsigned char *)trans = 'T';
#line 364 "dsygvd.f"
	    } else {
#line 365 "dsygvd.f"
		*(unsigned char *)trans = 'N';
#line 366 "dsygvd.f"
	    }

#line 368 "dsygvd.f"
	    dtrmm_("Left", uplo, trans, "Non-unit", n, n, &c_b11, &b[b_offset]
		    , ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (ftnlen)1,
		     (ftnlen)8);
#line 370 "dsygvd.f"
	}
#line 371 "dsygvd.f"
    }

#line 373 "dsygvd.f"
    work[1] = (doublereal) lopt;
#line 374 "dsygvd.f"
    iwork[1] = liopt;

#line 376 "dsygvd.f"
    return 0;

/*     End of DSYGVD */

} /* dsygvd_ */

