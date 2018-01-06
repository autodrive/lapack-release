#line 1 "chesv_aa.f"
/* chesv_aa.f -- translated by f2c (version 20100827).
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

#line 1 "chesv_aa.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> CHESV_AA computes the solution to system of linear equations A * X = B for HE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHESV_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chesv_a
a.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chesv_a
a.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chesv_a
a.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHESV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, */
/*                            LWORK, INFO ) */

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
/* > CHESV_AA computes the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > Aasen's algorithm is used to factor A as */
/* >    A = U * T * U**H,  if UPLO = 'U', or */
/* >    A = L * T * L**H,  if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is Hermitian and tridiagonal. The factored form */
/* > of A is then used to solve the system of equations A * X = B. */
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
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the tridiagonal matrix T and the */
/* >          multipliers used to obtain the factor U or L from the */
/* >          factorization A = U*T*U**H or A = L*T*L**H as computed by */
/* >          CHETRF_AA. */
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
/* >          The length of WORK.  LWORK >= MAX(1,2*N,3*N-2), and for best */
/* >          performance LWORK >= MAX(1,N*NB), where NB is the optimal */
/* >          blocksize for CHETRF. */
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

/* > \ingroup complexHEsolve */

/*  ===================================================================== */
/* Subroutine */ int chesv_aa__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int chetrf_aa__(char *, integer *, doublecomplex *
	    , integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), chetrs_aa__(char *, integer *, integer *, doublecomplex *
	    , integer *, integer *, doublecomplex *, integer *, doublecomplex 
	    *, integer *, integer *, ftnlen);
    static integer lwkopt_hetrf__, lwkopt_hetrs__;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
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

#line 200 "chesv_aa.f"
    /* Parameter adjustments */
#line 200 "chesv_aa.f"
    a_dim1 = *lda;
#line 200 "chesv_aa.f"
    a_offset = 1 + a_dim1;
#line 200 "chesv_aa.f"
    a -= a_offset;
#line 200 "chesv_aa.f"
    --ipiv;
#line 200 "chesv_aa.f"
    b_dim1 = *ldb;
#line 200 "chesv_aa.f"
    b_offset = 1 + b_dim1;
#line 200 "chesv_aa.f"
    b -= b_offset;
#line 200 "chesv_aa.f"
    --work;
#line 200 "chesv_aa.f"

#line 200 "chesv_aa.f"
    /* Function Body */
#line 200 "chesv_aa.f"
    *info = 0;
#line 201 "chesv_aa.f"
    lquery = *lwork == -1;
#line 202 "chesv_aa.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 203 "chesv_aa.f"
	*info = -1;
#line 204 "chesv_aa.f"
    } else if (*n < 0) {
#line 205 "chesv_aa.f"
	*info = -2;
#line 206 "chesv_aa.f"
    } else if (*nrhs < 0) {
#line 207 "chesv_aa.f"
	*info = -3;
#line 208 "chesv_aa.f"
    } else if (*lda < max(1,*n)) {
#line 209 "chesv_aa.f"
	*info = -5;
#line 210 "chesv_aa.f"
    } else if (*ldb < max(1,*n)) {
#line 211 "chesv_aa.f"
	*info = -8;
#line 212 "chesv_aa.f"
    }

#line 214 "chesv_aa.f"
    if (*info == 0) {
#line 215 "chesv_aa.f"
	chetrf_aa__(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &c_n1, 
		info, (ftnlen)1);
#line 216 "chesv_aa.f"
	lwkopt_hetrf__ = (integer) work[1].r;
#line 217 "chesv_aa.f"
	chetrs_aa__(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], 
		ldb, &work[1], &c_n1, info, (ftnlen)1);
#line 219 "chesv_aa.f"
	lwkopt_hetrs__ = (integer) work[1].r;
#line 220 "chesv_aa.f"
	lwkopt = max(lwkopt_hetrf__,lwkopt_hetrs__);
#line 221 "chesv_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 222 "chesv_aa.f"
	if (*lwork < lwkopt && ! lquery) {
#line 223 "chesv_aa.f"
	    *info = -10;
#line 224 "chesv_aa.f"
	}
#line 225 "chesv_aa.f"
    }

#line 227 "chesv_aa.f"
    if (*info != 0) {
#line 228 "chesv_aa.f"
	i__1 = -(*info);
#line 228 "chesv_aa.f"
	xerbla_("CHESV_AA ", &i__1, (ftnlen)9);
#line 229 "chesv_aa.f"
	return 0;
#line 230 "chesv_aa.f"
    } else if (lquery) {
#line 231 "chesv_aa.f"
	return 0;
#line 232 "chesv_aa.f"
    }

/*     Compute the factorization A = U*T*U**H or A = L*T*L**H. */

#line 236 "chesv_aa.f"
    chetrf_aa__(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info, (
	    ftnlen)1);
#line 237 "chesv_aa.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 241 "chesv_aa.f"
	chetrs_aa__(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], 
		ldb, &work[1], lwork, info, (ftnlen)1);

#line 244 "chesv_aa.f"
    }

#line 246 "chesv_aa.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 248 "chesv_aa.f"
    return 0;

/*     End of CHESV_AA */

} /* chesv_aa__ */

