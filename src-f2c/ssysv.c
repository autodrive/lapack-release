#line 1 "ssysv.f"
/* ssysv.f -- translated by f2c (version 20100827).
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

#line 1 "ssysv.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> SSYSV computes the solution to system of linear equations A * X = B for SY matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, */
/*                         LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYSV computes the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > The diagonal pivoting method is used to factor A as */
/* >    A = U * D * U**T,  if UPLO = 'U', or */
/* >    A = L * D * L**T,  if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then */
/* > used to solve the system of equations A * X = B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of linear equations, i.e., the order of the */
/* >          matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the block diagonal matrix D and the */
/* >          multipliers used to obtain the factor U or L from the */
/* >          factorization A = U*D*U**T or A = L*D*L**T as computed by */
/* >          SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D, as */
/* >          determined by SSYTRF.  If IPIV(k) > 0, then rows and columns */
/* >          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1 */
/* >          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, */
/* >          then rows and columns k-1 and -IPIV(k) were interchanged and */
/* >          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and */
/* >          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and */
/* >          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2 */
/* >          diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the N-by-NRHS right hand side matrix B. */
/* >          On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
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
/* >          The length of WORK.  LWORK >= 1, and for best performance */
/* >          LWORK >= max(1,N*NB), where NB is the optimal blocksize for */
/* >          SSYTRF. */
/* >          for LWORK < N, TRS will be done with Level BLAS 2 */
/* >          for LWORK >= N, TRS will be done with Level BLAS 3 */
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
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization */
/* >               has been completed, but the block diagonal matrix D is */
/* >               exactly singular, so the solution could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realSYsolve */

/*  ===================================================================== */
/* Subroutine */ int ssysv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int ssytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     ssytrs_(char *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), ssytrs2_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 208 "ssysv.f"
    /* Parameter adjustments */
#line 208 "ssysv.f"
    a_dim1 = *lda;
#line 208 "ssysv.f"
    a_offset = 1 + a_dim1;
#line 208 "ssysv.f"
    a -= a_offset;
#line 208 "ssysv.f"
    --ipiv;
#line 208 "ssysv.f"
    b_dim1 = *ldb;
#line 208 "ssysv.f"
    b_offset = 1 + b_dim1;
#line 208 "ssysv.f"
    b -= b_offset;
#line 208 "ssysv.f"
    --work;
#line 208 "ssysv.f"

#line 208 "ssysv.f"
    /* Function Body */
#line 208 "ssysv.f"
    *info = 0;
#line 209 "ssysv.f"
    lquery = *lwork == -1;
#line 210 "ssysv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 211 "ssysv.f"
	*info = -1;
#line 212 "ssysv.f"
    } else if (*n < 0) {
#line 213 "ssysv.f"
	*info = -2;
#line 214 "ssysv.f"
    } else if (*nrhs < 0) {
#line 215 "ssysv.f"
	*info = -3;
#line 216 "ssysv.f"
    } else if (*lda < max(1,*n)) {
#line 217 "ssysv.f"
	*info = -5;
#line 218 "ssysv.f"
    } else if (*ldb < max(1,*n)) {
#line 219 "ssysv.f"
	*info = -8;
#line 220 "ssysv.f"
    } else if (*lwork < 1 && ! lquery) {
#line 221 "ssysv.f"
	*info = -10;
#line 222 "ssysv.f"
    }

#line 224 "ssysv.f"
    if (*info == 0) {
#line 225 "ssysv.f"
	if (*n == 0) {
#line 226 "ssysv.f"
	    lwkopt = 1;
#line 227 "ssysv.f"
	} else {
#line 228 "ssysv.f"
	    ssytrf_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &c_n1, 
		    info, (ftnlen)1);
#line 229 "ssysv.f"
	    lwkopt = (integer) work[1];
#line 230 "ssysv.f"
	}
#line 231 "ssysv.f"
	work[1] = (doublereal) lwkopt;
#line 232 "ssysv.f"
    }

#line 234 "ssysv.f"
    if (*info != 0) {
#line 235 "ssysv.f"
	i__1 = -(*info);
#line 235 "ssysv.f"
	xerbla_("SSYSV ", &i__1, (ftnlen)6);
#line 236 "ssysv.f"
	return 0;
#line 237 "ssysv.f"
    } else if (lquery) {
#line 238 "ssysv.f"
	return 0;
#line 239 "ssysv.f"
    }

/*     Compute the factorization A = U*D*U**T or A = L*D*L**T. */

#line 243 "ssysv.f"
    ssytrf_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info, (
	    ftnlen)1);
#line 244 "ssysv.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 248 "ssysv.f"
	if (*lwork < *n) {

/*        Solve with TRS ( Use Level BLAS 2) */

#line 252 "ssysv.f"
	    ssytrs_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], 
		    ldb, info, (ftnlen)1);

#line 254 "ssysv.f"
	} else {

/*        Solve with TRS2 ( Use Level BLAS 3) */

#line 258 "ssysv.f"
	    ssytrs2_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset],
		     ldb, &work[1], info, (ftnlen)1);

#line 260 "ssysv.f"
	}

#line 262 "ssysv.f"
    }

#line 264 "ssysv.f"
    work[1] = (doublereal) lwkopt;

#line 266 "ssysv.f"
    return 0;

/*     End of SSYSV */

} /* ssysv_ */

