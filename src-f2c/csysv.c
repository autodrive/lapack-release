#line 1 "csysv.f"
/* csysv.f -- translated by f2c (version 20100827).
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

#line 1 "csysv.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> CSYSV computes the solution to system of linear equations A * X = B for SY matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csysv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csysv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csysv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, */
/*                         LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYSV computes the solution to a complex system of linear equations */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          CSYTRF. */
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
/* >          determined by CSYTRF.  If IPIV(k) > 0, then rows and columns */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK.  LWORK >= 1, and for best performance */
/* >          LWORK >= max(1,N*NB), where NB is the optimal blocksize for */
/* >          CSYTRF. */
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

/* > \date November 2011 */

/* > \ingroup complexSYsolve */

/*  ===================================================================== */
/* Subroutine */ int csysv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), csytrf_(
	    char *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int csytrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen), csytrs2_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 208 "csysv.f"
    /* Parameter adjustments */
#line 208 "csysv.f"
    a_dim1 = *lda;
#line 208 "csysv.f"
    a_offset = 1 + a_dim1;
#line 208 "csysv.f"
    a -= a_offset;
#line 208 "csysv.f"
    --ipiv;
#line 208 "csysv.f"
    b_dim1 = *ldb;
#line 208 "csysv.f"
    b_offset = 1 + b_dim1;
#line 208 "csysv.f"
    b -= b_offset;
#line 208 "csysv.f"
    --work;
#line 208 "csysv.f"

#line 208 "csysv.f"
    /* Function Body */
#line 208 "csysv.f"
    *info = 0;
#line 209 "csysv.f"
    lquery = *lwork == -1;
#line 210 "csysv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 211 "csysv.f"
	*info = -1;
#line 212 "csysv.f"
    } else if (*n < 0) {
#line 213 "csysv.f"
	*info = -2;
#line 214 "csysv.f"
    } else if (*nrhs < 0) {
#line 215 "csysv.f"
	*info = -3;
#line 216 "csysv.f"
    } else if (*lda < max(1,*n)) {
#line 217 "csysv.f"
	*info = -5;
#line 218 "csysv.f"
    } else if (*ldb < max(1,*n)) {
#line 219 "csysv.f"
	*info = -8;
#line 220 "csysv.f"
    } else if (*lwork < 1 && ! lquery) {
#line 221 "csysv.f"
	*info = -10;
#line 222 "csysv.f"
    }

#line 224 "csysv.f"
    if (*info == 0) {
#line 225 "csysv.f"
	if (*n == 0) {
#line 226 "csysv.f"
	    lwkopt = 1;
#line 227 "csysv.f"
	} else {
#line 228 "csysv.f"
	    csytrf_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &c_n1, 
		    info, (ftnlen)1);
#line 229 "csysv.f"
	    lwkopt = (integer) work[1].r;
#line 230 "csysv.f"
	}
#line 231 "csysv.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 232 "csysv.f"
    }

#line 234 "csysv.f"
    if (*info != 0) {
#line 235 "csysv.f"
	i__1 = -(*info);
#line 235 "csysv.f"
	xerbla_("CSYSV ", &i__1, (ftnlen)6);
#line 236 "csysv.f"
	return 0;
#line 237 "csysv.f"
    } else if (lquery) {
#line 238 "csysv.f"
	return 0;
#line 239 "csysv.f"
    }

/*     Compute the factorization A = U*D*U**T or A = L*D*L**T. */

#line 243 "csysv.f"
    csytrf_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info, (
	    ftnlen)1);
#line 244 "csysv.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 248 "csysv.f"
	if (*lwork < *n) {

/*        Solve with TRS ( Use Level BLAS 2) */

#line 252 "csysv.f"
	    csytrs_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], 
		    ldb, info, (ftnlen)1);

#line 254 "csysv.f"
	} else {

/*        Solve with TRS2 ( Use Level BLAS 3) */

#line 258 "csysv.f"
	    csytrs2_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset],
		     ldb, &work[1], info, (ftnlen)1);

#line 260 "csysv.f"
	}

#line 262 "csysv.f"
    }

#line 264 "csysv.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 266 "csysv.f"
    return 0;

/*     End of CSYSV */

} /* csysv_ */

