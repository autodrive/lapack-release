#line 1 "ssygv_2stage.f"
/* ssygv_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "ssygv_2stage.f"
/* Table of constant values */

static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static doublereal c_b26 = 1.;

/* > \brief \b SSYGV_2STAGE */

/*  @generated from dsygv_2stage.f, fortran d -> s, Sun Nov  6 12:54:29 2016 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYGV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygv_2
stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygv_2
stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygv_2
stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYGV_2STAGE( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, */
/*                                WORK, LWORK, INFO ) */

/*       IMPLICIT NONE */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYGV_2STAGE computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be symmetric and B is also */
/* > positive definite. */
/* > This routine use the 2stage technique for the reduction to tridiagonal */
/* > which showed higher performance on recent architecture and for large */
/*  sizes N>2000. */
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
/* >                  Not available in this release. */
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
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          B is REAL array, dimension (LDB, N) */
/* >          On entry, the symmetric positive definite matrix B. */
/* >          If UPLO = 'U', the leading N-by-N upper triangular part of B */
/* >          contains the upper triangular part of the matrix B. */
/* >          If UPLO = 'L', the leading N-by-N lower triangular part of B */
/* >          contains the lower triangular part of the matrix B. */
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
/* >          W is REAL array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
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
/* >          The length of the array WORK. LWORK >= 1, when N <= 1; */
/* >          otherwise */
/* >          If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* >                                   LWORK = MAX(1, dimension) where */
/* >                                   dimension = max(stage1,stage2) + (KD+1)*N + 2*N */
/* >                                             = N*KD + N*max(KD+1,FACTOPTNB) */
/* >                                               + max(2*KD*KD, KD*NTHREADS) */
/* >                                               + (KD+1)*N + 2*N */
/* >                                   where KD is the blocking size of the reduction, */
/* >                                   FACTOPTNB is the blocking used by the QR or LQ */
/* >                                   algorithm, usually FACTOPTNB=128 is a good choice */
/* >                                   NTHREADS is the number of threads used when */
/* >                                   openMP compilation is enabled, otherwise =1. */
/* >          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  SPOTRF or SSYEV returned an error code: */
/* >             <= N:  if INFO = i, SSYEV failed to converge; */
/* >                    i off-diagonal elements of an intermediate */
/* >                    tridiagonal form did not converge to zero; */
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

/* > \ingroup realSYeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  All details about the 2stage techniques are available in: */
/* > */
/* >  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* >  Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* >  using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* >  of 2011 International Conference for High Performance Computing, */
/* >  Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* >  Article 8 , 11 pages. */
/* >  http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* >  A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* >  An improved parallel singular value algorithm and its implementation */
/* >  for multicore hardware, In Proceedings of 2013 International Conference */
/* >  for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* >  Denver, Colorado, USA, 2013. */
/* >  Article 90, 12 pages. */
/* >  http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* >  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* >  A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* >  calculations based on fine-grained memory aware tasks. */
/* >  International Journal of High Performance Computing Applications. */
/* >  Volume 28 Issue 2, Pages 196-209, May 2014. */
/* >  http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int ssygv_2stage__(integer *itype, char *jobz, char *uplo, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static integer ib, kd, neig;
    extern /* Subroutine */ int ssyev_2stage__(char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lhtrd, lwmin;
    static char trans[1];
    static logical upper;
    static integer lwtrd;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical wantz;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int spotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int ssygst_(integer *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);



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

#line 272 "ssygv_2stage.f"
    /* Parameter adjustments */
#line 272 "ssygv_2stage.f"
    a_dim1 = *lda;
#line 272 "ssygv_2stage.f"
    a_offset = 1 + a_dim1;
#line 272 "ssygv_2stage.f"
    a -= a_offset;
#line 272 "ssygv_2stage.f"
    b_dim1 = *ldb;
#line 272 "ssygv_2stage.f"
    b_offset = 1 + b_dim1;
#line 272 "ssygv_2stage.f"
    b -= b_offset;
#line 272 "ssygv_2stage.f"
    --w;
#line 272 "ssygv_2stage.f"
    --work;
#line 272 "ssygv_2stage.f"

#line 272 "ssygv_2stage.f"
    /* Function Body */
#line 272 "ssygv_2stage.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 273 "ssygv_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 274 "ssygv_2stage.f"
    lquery = *lwork == -1;

#line 276 "ssygv_2stage.f"
    *info = 0;
#line 277 "ssygv_2stage.f"
    if (*itype < 1 || *itype > 3) {
#line 278 "ssygv_2stage.f"
	*info = -1;
#line 279 "ssygv_2stage.f"
    } else if (! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 280 "ssygv_2stage.f"
	*info = -2;
#line 281 "ssygv_2stage.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 282 "ssygv_2stage.f"
	*info = -3;
#line 283 "ssygv_2stage.f"
    } else if (*n < 0) {
#line 284 "ssygv_2stage.f"
	*info = -4;
#line 285 "ssygv_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 286 "ssygv_2stage.f"
	*info = -6;
#line 287 "ssygv_2stage.f"
    } else if (*ldb < max(1,*n)) {
#line 288 "ssygv_2stage.f"
	*info = -8;
#line 289 "ssygv_2stage.f"
    }

#line 291 "ssygv_2stage.f"
    if (*info == 0) {
#line 292 "ssygv_2stage.f"
	kd = ilaenv_(&c__17, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 293 "ssygv_2stage.f"
	ib = ilaenv_(&c__18, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 294 "ssygv_2stage.f"
	lhtrd = ilaenv_(&c__19, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 295 "ssygv_2stage.f"
	lwtrd = ilaenv_(&c__20, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1, (
		ftnlen)13, (ftnlen)1);
#line 296 "ssygv_2stage.f"
	lwmin = (*n << 1) + lhtrd + lwtrd;
#line 297 "ssygv_2stage.f"
	work[1] = (doublereal) lwmin;

#line 299 "ssygv_2stage.f"
	if (*lwork < lwmin && ! lquery) {
#line 300 "ssygv_2stage.f"
	    *info = -11;
#line 301 "ssygv_2stage.f"
	}
#line 302 "ssygv_2stage.f"
    }

#line 304 "ssygv_2stage.f"
    if (*info != 0) {
#line 305 "ssygv_2stage.f"
	i__1 = -(*info);
#line 305 "ssygv_2stage.f"
	xerbla_("SSYGV_2STAGE ", &i__1, (ftnlen)13);
#line 306 "ssygv_2stage.f"
	return 0;
#line 307 "ssygv_2stage.f"
    } else if (lquery) {
#line 308 "ssygv_2stage.f"
	return 0;
#line 309 "ssygv_2stage.f"
    }

/*     Quick return if possible */

#line 313 "ssygv_2stage.f"
    if (*n == 0) {
#line 313 "ssygv_2stage.f"
	return 0;
#line 313 "ssygv_2stage.f"
    }

/*     Form a Cholesky factorization of B. */

#line 318 "ssygv_2stage.f"
    spotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 319 "ssygv_2stage.f"
    if (*info != 0) {
#line 320 "ssygv_2stage.f"
	*info = *n + *info;
#line 321 "ssygv_2stage.f"
	return 0;
#line 322 "ssygv_2stage.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 326 "ssygv_2stage.f"
    ssygst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 327 "ssygv_2stage.f"
    ssyev_2stage__(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, 
	    info, (ftnlen)1, (ftnlen)1);

#line 329 "ssygv_2stage.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 333 "ssygv_2stage.f"
	neig = *n;
#line 334 "ssygv_2stage.f"
	if (*info > 0) {
#line 334 "ssygv_2stage.f"
	    neig = *info - 1;
#line 334 "ssygv_2stage.f"
	}
#line 336 "ssygv_2stage.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 341 "ssygv_2stage.f"
	    if (upper) {
#line 342 "ssygv_2stage.f"
		*(unsigned char *)trans = 'N';
#line 343 "ssygv_2stage.f"
	    } else {
#line 344 "ssygv_2stage.f"
		*(unsigned char *)trans = 'T';
#line 345 "ssygv_2stage.f"
	    }

#line 347 "ssygv_2stage.f"
	    strsm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b26, &b[
		    b_offset], ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);

#line 350 "ssygv_2stage.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 355 "ssygv_2stage.f"
	    if (upper) {
#line 356 "ssygv_2stage.f"
		*(unsigned char *)trans = 'T';
#line 357 "ssygv_2stage.f"
	    } else {
#line 358 "ssygv_2stage.f"
		*(unsigned char *)trans = 'N';
#line 359 "ssygv_2stage.f"
	    }

#line 361 "ssygv_2stage.f"
	    strmm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b26, &b[
		    b_offset], ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);
#line 363 "ssygv_2stage.f"
	}
#line 364 "ssygv_2stage.f"
    }

#line 366 "ssygv_2stage.f"
    work[1] = (doublereal) lwmin;
#line 367 "ssygv_2stage.f"
    return 0;

/*     End of SSYGV_2STAGE */

} /* ssygv_2stage__ */

