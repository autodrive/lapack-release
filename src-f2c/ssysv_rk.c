#line 1 "ssysv_rk.f"
/* ssysv_rk.f -- translated by f2c (version 20100827).
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

#line 1 "ssysv_rk.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> SSYSV_RK computes the solution to system of linear equations A * X = B for SY matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYSV_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv_r
k.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv_r
k.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv_r
k.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, */
/*                            WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > SSYSV_RK computes the solution to a real system of linear */
/* > equations A * X = B, where A is an N-by-N symmetric matrix */
/* > and X and B are N-by-NRHS matrices. */
/* > */
/* > The bounded Bunch-Kaufman (rook) diagonal pivoting method is used */
/* > to factor A as */
/* >    A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or */
/* >    A = P*L*D*(L**T)*(P**T),  if UPLO = 'L', */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > SSYTRF_RK is called to compute the factorization of a real */
/* > symmetric matrix.  The factored form of A is then used to solve */
/* > the system of equations A * X = B by calling BLAS3 routine SSYTRS_3. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored: */
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
/* >          On entry, the symmetric matrix A. */
/* >            If UPLO = 'U': the leading N-by-N upper triangular part */
/* >            of A contains the upper triangular part of the matrix A, */
/* >            and the strictly lower triangular part of A is not */
/* >            referenced. */
/* > */
/* >            If UPLO = 'L': the leading N-by-N lower triangular part */
/* >            of A contains the lower triangular part of the matrix A, */
/* >            and the strictly upper triangular part of A is not */
/* >            referenced. */
/* > */
/* >          On exit, if INFO = 0, diagonal of the block diagonal */
/* >          matrix D and factors U or L  as computed by SSYTRF_RK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          For more info see the description of DSYTRF_RK routine. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* >          On exit, contains the output computed by the factorization */
/* >          routine DSYTRF_RK, i.e. the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          NOTE: For 1-by-1 diagonal block D(k), where */
/* >          1 <= k <= N, the element E(k) is set to 0 in both */
/* >          UPLO = 'U' or UPLO = 'L' cases. */
/* > */
/* >          For more info see the description of DSYTRF_RK routine. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D, */
/* >          as determined by SSYTRF_RK. */
/* > */
/* >          For more info see the description of DSYTRF_RK routine. */
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
/* >          WORK is REAL array, dimension ( MAX(1,LWORK) ). */
/* >          Work array used in the factorization stage. */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK.  LWORK >= 1. For best performance */
/* >          of factorization stage LWORK >= max(1,N*NB), where NB is */
/* >          the optimal blocksize for DSYTRF_RK. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; */
/* >          the routine only calculates the optimal size of the WORK */
/* >          array for factorization stage, returns this value as */
/* >          the first entry of the WORK array, and no error message */
/* >          related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* > */
/* >          < 0: If INFO = -k, the k-th argument had an illegal value */
/* > */
/* >          > 0: If INFO = k, the matrix A is singular, because: */
/* >                 If UPLO = 'U': column k in the upper */
/* >                 triangular part of A contains all zeros. */
/* >                 If UPLO = 'L': column k in the lower */
/* >                 triangular part of A contains all zeros. */
/* > */
/* >               Therefore D(k,k) is exactly zero, and superdiagonal */
/* >               elements of column k of U (or subdiagonal elements of */
/* >               column k of L ) are all zeros. The factorization has */
/* >               been completed, but the block diagonal matrix D is */
/* >               exactly singular, and division by zero will occur if */
/* >               it is used to solve a system of equations. */
/* > */
/* >               NOTE: INFO only stores the first occurrence of */
/* >               a singularity, any subsequent occurrence of singularity */
/* >               is not stored in INFO even though the factorization */
/* >               always completes. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup singleSYsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int ssysv_rk__(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *e, integer *ipiv, doublereal 
	*b, integer *ldb, doublereal *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int ssytrs_3__(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), ssytrf_rk__(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
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

#line 265 "ssysv_rk.f"
    /* Parameter adjustments */
#line 265 "ssysv_rk.f"
    a_dim1 = *lda;
#line 265 "ssysv_rk.f"
    a_offset = 1 + a_dim1;
#line 265 "ssysv_rk.f"
    a -= a_offset;
#line 265 "ssysv_rk.f"
    --e;
#line 265 "ssysv_rk.f"
    --ipiv;
#line 265 "ssysv_rk.f"
    b_dim1 = *ldb;
#line 265 "ssysv_rk.f"
    b_offset = 1 + b_dim1;
#line 265 "ssysv_rk.f"
    b -= b_offset;
#line 265 "ssysv_rk.f"
    --work;
#line 265 "ssysv_rk.f"

#line 265 "ssysv_rk.f"
    /* Function Body */
#line 265 "ssysv_rk.f"
    *info = 0;
#line 266 "ssysv_rk.f"
    lquery = *lwork == -1;
#line 267 "ssysv_rk.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 268 "ssysv_rk.f"
	*info = -1;
#line 269 "ssysv_rk.f"
    } else if (*n < 0) {
#line 270 "ssysv_rk.f"
	*info = -2;
#line 271 "ssysv_rk.f"
    } else if (*nrhs < 0) {
#line 272 "ssysv_rk.f"
	*info = -3;
#line 273 "ssysv_rk.f"
    } else if (*lda < max(1,*n)) {
#line 274 "ssysv_rk.f"
	*info = -5;
#line 275 "ssysv_rk.f"
    } else if (*ldb < max(1,*n)) {
#line 276 "ssysv_rk.f"
	*info = -9;
#line 277 "ssysv_rk.f"
    } else if (*lwork < 1 && ! lquery) {
#line 278 "ssysv_rk.f"
	*info = -11;
#line 279 "ssysv_rk.f"
    }

#line 281 "ssysv_rk.f"
    if (*info == 0) {
#line 282 "ssysv_rk.f"
	if (*n == 0) {
#line 283 "ssysv_rk.f"
	    lwkopt = 1;
#line 284 "ssysv_rk.f"
	} else {
#line 285 "ssysv_rk.f"
	    ssytrf_rk__(uplo, n, &a[a_offset], lda, &e[1], &ipiv[1], &work[1],
		     &c_n1, info, (ftnlen)1);
#line 286 "ssysv_rk.f"
	    lwkopt = (integer) work[1];
#line 287 "ssysv_rk.f"
	}
#line 288 "ssysv_rk.f"
	work[1] = (doublereal) lwkopt;
#line 289 "ssysv_rk.f"
    }

#line 291 "ssysv_rk.f"
    if (*info != 0) {
#line 292 "ssysv_rk.f"
	i__1 = -(*info);
#line 292 "ssysv_rk.f"
	xerbla_("SSYSV_RK ", &i__1, (ftnlen)9);
#line 293 "ssysv_rk.f"
	return 0;
#line 294 "ssysv_rk.f"
    } else if (lquery) {
#line 295 "ssysv_rk.f"
	return 0;
#line 296 "ssysv_rk.f"
    }

/*     Compute the factorization A = P*U*D*(U**T)*(P**T) or */
/*     A = P*U*D*(U**T)*(P**T). */

#line 301 "ssysv_rk.f"
    ssytrf_rk__(uplo, n, &a[a_offset], lda, &e[1], &ipiv[1], &work[1], lwork, 
	    info, (ftnlen)1);

#line 303 "ssysv_rk.f"
    if (*info == 0) {

/*        Solve the system A*X = B with BLAS3 solver, overwriting B with X. */

#line 307 "ssysv_rk.f"
	ssytrs_3__(uplo, n, nrhs, &a[a_offset], lda, &e[1], &ipiv[1], &b[
		b_offset], ldb, info, (ftnlen)1);

#line 309 "ssysv_rk.f"
    }

#line 311 "ssysv_rk.f"
    work[1] = (doublereal) lwkopt;

#line 313 "ssysv_rk.f"
    return 0;

/*     End of SSYSV_RK */

} /* ssysv_rk__ */

