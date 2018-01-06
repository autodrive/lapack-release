#line 1 "zsysv_aa.f"
/* zsysv_aa.f -- translated by f2c (version 20100827).
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

#line 1 "zsysv_aa.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> ZSYSV_AA computes the solution to system of linear equations A * X = B for SY matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYSV_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysv_a
a.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysv_a
a.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysv_a
a.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYSV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, */
/*                            LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYSV computes the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > Aasen's algorithm is used to factor A as */
/* >    A = U * T * U**T,  if UPLO = 'U', or */
/* >    A = L * T * L**T,  if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is symmetric tridiagonal. The factored */
/* > form of A is then used to solve the system of equations A * X = B. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the tridiagonal matrix T and the */
/* >          multipliers used to obtain the factor U or L from the */
/* >          factorization A = U*T*U**T or A = L*T*L**T as computed by */
/* >          ZSYTRF. */
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
/* >          On exit, it contains the details of the interchanges, i.e., */
/* >          the row and column k of A were interchanged with the */
/* >          row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK.  LWORK >= MAX(1,2*N,3*N-2), and for */
/* >          the best performance, LWORK >= MAX(1,N*NB), where NB is */
/* >          the optimal blocksize for ZSYTRF_AA. */
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

/* > \date November 2017 */

/* > \ingroup complex16SYsolve */

/*  ===================================================================== */
/* Subroutine */ int zsysv_aa__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int zsytrf_aa__(char *, integer *, doublecomplex *
	    , integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), zsytrs_aa__(char *, integer *, integer *, doublecomplex *
	    , integer *, integer *, doublecomplex *, integer *, doublecomplex 
	    *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lwkopt_sytrf__, lwkopt_sytrs__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK driver routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 200 "zsysv_aa.f"
    /* Parameter adjustments */
#line 200 "zsysv_aa.f"
    a_dim1 = *lda;
#line 200 "zsysv_aa.f"
    a_offset = 1 + a_dim1;
#line 200 "zsysv_aa.f"
    a -= a_offset;
#line 200 "zsysv_aa.f"
    --ipiv;
#line 200 "zsysv_aa.f"
    b_dim1 = *ldb;
#line 200 "zsysv_aa.f"
    b_offset = 1 + b_dim1;
#line 200 "zsysv_aa.f"
    b -= b_offset;
#line 200 "zsysv_aa.f"
    --work;
#line 200 "zsysv_aa.f"

#line 200 "zsysv_aa.f"
    /* Function Body */
#line 200 "zsysv_aa.f"
    *info = 0;
#line 201 "zsysv_aa.f"
    lquery = *lwork == -1;
#line 202 "zsysv_aa.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 203 "zsysv_aa.f"
	*info = -1;
#line 204 "zsysv_aa.f"
    } else if (*n < 0) {
#line 205 "zsysv_aa.f"
	*info = -2;
#line 206 "zsysv_aa.f"
    } else if (*nrhs < 0) {
#line 207 "zsysv_aa.f"
	*info = -3;
#line 208 "zsysv_aa.f"
    } else if (*lda < max(1,*n)) {
#line 209 "zsysv_aa.f"
	*info = -5;
#line 210 "zsysv_aa.f"
    } else if (*ldb < max(1,*n)) {
#line 211 "zsysv_aa.f"
	*info = -8;
#line 212 "zsysv_aa.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 212 "zsysv_aa.f"
	i__1 = *n << 1, i__2 = *n * 3 - 2;
#line 212 "zsysv_aa.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 213 "zsysv_aa.f"
	    *info = -10;
#line 214 "zsysv_aa.f"
	}
#line 214 "zsysv_aa.f"
    }

#line 216 "zsysv_aa.f"
    if (*info == 0) {
#line 217 "zsysv_aa.f"
	zsytrf_aa__(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &c_n1, 
		info, (ftnlen)1);
#line 218 "zsysv_aa.f"
	lwkopt_sytrf__ = (integer) work[1].r;
#line 219 "zsysv_aa.f"
	zsytrs_aa__(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], 
		ldb, &work[1], &c_n1, info, (ftnlen)1);
#line 221 "zsysv_aa.f"
	lwkopt_sytrs__ = (integer) work[1].r;
#line 222 "zsysv_aa.f"
	lwkopt = max(lwkopt_sytrf__,lwkopt_sytrs__);
#line 223 "zsysv_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 224 "zsysv_aa.f"
	if (*lwork < lwkopt && ! lquery) {
#line 225 "zsysv_aa.f"
	    *info = -10;
#line 226 "zsysv_aa.f"
	}
#line 227 "zsysv_aa.f"
    }

#line 229 "zsysv_aa.f"
    if (*info != 0) {
#line 230 "zsysv_aa.f"
	i__1 = -(*info);
#line 230 "zsysv_aa.f"
	xerbla_("ZSYSV_AA ", &i__1, (ftnlen)9);
#line 231 "zsysv_aa.f"
	return 0;
#line 232 "zsysv_aa.f"
    } else if (lquery) {
#line 233 "zsysv_aa.f"
	return 0;
#line 234 "zsysv_aa.f"
    }

/*     Compute the factorization A = U*T*U**T or A = L*T*L**T. */

#line 238 "zsysv_aa.f"
    zsytrf_aa__(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info, (
	    ftnlen)1);
#line 239 "zsysv_aa.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 243 "zsysv_aa.f"
	zsytrs_aa__(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], 
		ldb, &work[1], lwork, info, (ftnlen)1);

#line 246 "zsysv_aa.f"
    }

#line 248 "zsysv_aa.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 250 "zsysv_aa.f"
    return 0;

/*     End of ZSYSV_AA */

} /* zsysv_aa__ */

