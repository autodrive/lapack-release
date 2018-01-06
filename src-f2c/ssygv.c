#line 1 "ssygv.f"
/* ssygv.f -- translated by f2c (version 20100827).
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

#line 1 "ssygv.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b16 = 1.;

/* > \brief \b SSYGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, */
/*                         LWORK, INFO ) */

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
/* > SSYGV computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be symmetric and B is also */
/* > positive definite. */
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
/* >          The length of the array WORK.  LWORK >= max(1,3*N-1). */
/* >          For optimal efficiency, LWORK >= (NB+2)*N, */
/* >          where NB is the blocksize for SSYTRD returned by ILAENV. */
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

/* > \date November 2011 */

/* > \ingroup realSYeigen */

/*  ===================================================================== */
/* Subroutine */ int ssygv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *info, 
	ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer nb, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char trans[1];
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical wantz;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), ssyev_(
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), xerbla_(char 
	    *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkmin;
    extern /* Subroutine */ int spotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int ssygst_(integer *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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

#line 217 "ssygv.f"
    /* Parameter adjustments */
#line 217 "ssygv.f"
    a_dim1 = *lda;
#line 217 "ssygv.f"
    a_offset = 1 + a_dim1;
#line 217 "ssygv.f"
    a -= a_offset;
#line 217 "ssygv.f"
    b_dim1 = *ldb;
#line 217 "ssygv.f"
    b_offset = 1 + b_dim1;
#line 217 "ssygv.f"
    b -= b_offset;
#line 217 "ssygv.f"
    --w;
#line 217 "ssygv.f"
    --work;
#line 217 "ssygv.f"

#line 217 "ssygv.f"
    /* Function Body */
#line 217 "ssygv.f"
    wantz = lsame_(jobz, "V", (ftnlen)1, (ftnlen)1);
#line 218 "ssygv.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 219 "ssygv.f"
    lquery = *lwork == -1;

#line 221 "ssygv.f"
    *info = 0;
#line 222 "ssygv.f"
    if (*itype < 1 || *itype > 3) {
#line 223 "ssygv.f"
	*info = -1;
#line 224 "ssygv.f"
    } else if (! (wantz || lsame_(jobz, "N", (ftnlen)1, (ftnlen)1))) {
#line 225 "ssygv.f"
	*info = -2;
#line 226 "ssygv.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 227 "ssygv.f"
	*info = -3;
#line 228 "ssygv.f"
    } else if (*n < 0) {
#line 229 "ssygv.f"
	*info = -4;
#line 230 "ssygv.f"
    } else if (*lda < max(1,*n)) {
#line 231 "ssygv.f"
	*info = -6;
#line 232 "ssygv.f"
    } else if (*ldb < max(1,*n)) {
#line 233 "ssygv.f"
	*info = -8;
#line 234 "ssygv.f"
    }

#line 236 "ssygv.f"
    if (*info == 0) {
/* Computing MAX */
#line 237 "ssygv.f"
	i__1 = 1, i__2 = *n * 3 - 1;
#line 237 "ssygv.f"
	lwkmin = max(i__1,i__2);
#line 238 "ssygv.f"
	nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
#line 239 "ssygv.f"
	i__1 = lwkmin, i__2 = (nb + 2) * *n;
#line 239 "ssygv.f"
	lwkopt = max(i__1,i__2);
#line 240 "ssygv.f"
	work[1] = (doublereal) lwkopt;

#line 242 "ssygv.f"
	if (*lwork < lwkmin && ! lquery) {
#line 243 "ssygv.f"
	    *info = -11;
#line 244 "ssygv.f"
	}
#line 245 "ssygv.f"
    }

#line 247 "ssygv.f"
    if (*info != 0) {
#line 248 "ssygv.f"
	i__1 = -(*info);
#line 248 "ssygv.f"
	xerbla_("SSYGV ", &i__1, (ftnlen)6);
#line 249 "ssygv.f"
	return 0;
#line 250 "ssygv.f"
    } else if (lquery) {
#line 251 "ssygv.f"
	return 0;
#line 252 "ssygv.f"
    }

/*     Quick return if possible */

#line 256 "ssygv.f"
    if (*n == 0) {
#line 256 "ssygv.f"
	return 0;
#line 256 "ssygv.f"
    }

/*     Form a Cholesky factorization of B. */

#line 261 "ssygv.f"
    spotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
#line 262 "ssygv.f"
    if (*info != 0) {
#line 263 "ssygv.f"
	*info = *n + *info;
#line 264 "ssygv.f"
	return 0;
#line 265 "ssygv.f"
    }

/*     Transform problem to standard eigenvalue problem and solve. */

#line 269 "ssygv.f"
    ssygst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
	    ftnlen)1);
#line 270 "ssygv.f"
    ssyev_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, info, (
	    ftnlen)1, (ftnlen)1);

#line 272 "ssygv.f"
    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

#line 276 "ssygv.f"
	neig = *n;
#line 277 "ssygv.f"
	if (*info > 0) {
#line 277 "ssygv.f"
	    neig = *info - 1;
#line 277 "ssygv.f"
	}
#line 279 "ssygv.f"
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x; */
/*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */

#line 284 "ssygv.f"
	    if (upper) {
#line 285 "ssygv.f"
		*(unsigned char *)trans = 'N';
#line 286 "ssygv.f"
	    } else {
#line 287 "ssygv.f"
		*(unsigned char *)trans = 'T';
#line 288 "ssygv.f"
	    }

#line 290 "ssygv.f"
	    strsm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b16, &b[
		    b_offset], ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);

#line 293 "ssygv.f"
	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x; */
/*           backtransform eigenvectors: x = L*y or U**T*y */

#line 298 "ssygv.f"
	    if (upper) {
#line 299 "ssygv.f"
		*(unsigned char *)trans = 'T';
#line 300 "ssygv.f"
	    } else {
#line 301 "ssygv.f"
		*(unsigned char *)trans = 'N';
#line 302 "ssygv.f"
	    }

#line 304 "ssygv.f"
	    strmm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b16, &b[
		    b_offset], ldb, &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);
#line 306 "ssygv.f"
	}
#line 307 "ssygv.f"
    }

#line 309 "ssygv.f"
    work[1] = (doublereal) lwkopt;
#line 310 "ssygv.f"
    return 0;

/*     End of SSYGV */

} /* ssygv_ */

