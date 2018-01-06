#line 1 "chesv_aa_2stage.f"
/* chesv_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "chesv_aa_2stage.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> CHESV_AA_2STAGE computes the solution to system of linear equations A * X = B for HE matrices
</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHESV_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chesv_a
a_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chesv_a
a_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chesv_a
a_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE CHESV_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, */
/*                                  IPIV, IPIV2, B, LDB, WORK, LWORK, */
/*                                  INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, NRHS, LDA, LTB, LDB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IPIV2( * ) */
/*       COMPLEX            A( LDA, * ), TB( * ), B( LDB, *), WORK( * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHESV_AA_2STAGE computes the solution to a complex system of */
/* > linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > Aasen's 2-stage algorithm is used to factor A as */
/* >    A = U * T * U**H,  if UPLO = 'U', or */
/* >    A = L * T * L**H,  if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is Hermitian and band. The matrix T is */
/* > then LU-factored with partial pivoting. The factored form of A */
/* > is then used to solve the system of equations A * X = B. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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
/* >          The order of the matrix A.  N >= 0. */
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
/* >          On entry, the hermitian matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, L is stored below (or above) the subdiaonal blocks, */
/* >          when UPLO  is 'L' (or 'U'). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TB */
/* > \verbatim */
/* >          TB is COMPLEX array, dimension (LTB) */
/* >          On exit, details of the LU factorization of the band matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LTB */
/* > \verbatim */
/* >          The size of the array TB. LTB >= 4*N, internally */
/* >          used to select NB such that LTB >= (3*NB+1)*N. */
/* > */
/* >          If LTB = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of LTB, */
/* >          returns this value as the first entry of TB, and */
/* >          no error message related to LTB is issued by XERBLA. */
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
/* > \param[out] IPIV2 */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On exit, it contains the details of the interchanges, i.e., */
/* >          the row and column k of T were interchanged with the */
/* >          row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, the solution matrix X. */
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
/* >          WORK is COMPLEX workspace of size LWORK */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          The size of WORK. LWORK >= N, internally used to select NB */
/* >          such that LWORK >= N*NB. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of the WORK array, */
/* >          returns this value as the first entry of the WORK array, and */
/* >          no error message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, band LU factorization failed on i-th column */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int chesv_aa_2stage__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *tb, integer *ltb, 
	integer *ipiv, integer *ipiv2, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int chetrs_aa_2stage__(char *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical tquery, wquery;
    extern /* Subroutine */ int chetrf_aa_2stage__(char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.8.0) -- */
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

#line 224 "chesv_aa_2stage.f"
    /* Parameter adjustments */
#line 224 "chesv_aa_2stage.f"
    a_dim1 = *lda;
#line 224 "chesv_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 224 "chesv_aa_2stage.f"
    a -= a_offset;
#line 224 "chesv_aa_2stage.f"
    --tb;
#line 224 "chesv_aa_2stage.f"
    --ipiv;
#line 224 "chesv_aa_2stage.f"
    --ipiv2;
#line 224 "chesv_aa_2stage.f"
    b_dim1 = *ldb;
#line 224 "chesv_aa_2stage.f"
    b_offset = 1 + b_dim1;
#line 224 "chesv_aa_2stage.f"
    b -= b_offset;
#line 224 "chesv_aa_2stage.f"
    --work;
#line 224 "chesv_aa_2stage.f"

#line 224 "chesv_aa_2stage.f"
    /* Function Body */
#line 224 "chesv_aa_2stage.f"
    *info = 0;
#line 225 "chesv_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 226 "chesv_aa_2stage.f"
    wquery = *lwork == -1;
#line 227 "chesv_aa_2stage.f"
    tquery = *ltb == -1;
#line 228 "chesv_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 229 "chesv_aa_2stage.f"
	*info = -1;
#line 230 "chesv_aa_2stage.f"
    } else if (*n < 0) {
#line 231 "chesv_aa_2stage.f"
	*info = -2;
#line 232 "chesv_aa_2stage.f"
    } else if (*nrhs < 0) {
#line 233 "chesv_aa_2stage.f"
	*info = -3;
#line 234 "chesv_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 235 "chesv_aa_2stage.f"
	*info = -5;
#line 236 "chesv_aa_2stage.f"
    } else if (*ldb < max(1,*n)) {
#line 237 "chesv_aa_2stage.f"
	*info = -11;
#line 238 "chesv_aa_2stage.f"
    }

#line 240 "chesv_aa_2stage.f"
    if (*info == 0) {
#line 241 "chesv_aa_2stage.f"
	chetrf_aa_2stage__(uplo, n, &a[a_offset], lda, &tb[1], &c_n1, &ipiv[1]
		, &ipiv2[1], &work[1], &c_n1, info, (ftnlen)1);
#line 243 "chesv_aa_2stage.f"
	lwkopt = (integer) work[1].r;
#line 244 "chesv_aa_2stage.f"
	if (*ltb < (integer) tb[1].r && ! tquery) {
#line 245 "chesv_aa_2stage.f"
	    *info = -7;
#line 246 "chesv_aa_2stage.f"
	} else if (*lwork < lwkopt && ! wquery) {
#line 247 "chesv_aa_2stage.f"
	    *info = -13;
#line 248 "chesv_aa_2stage.f"
	}
#line 249 "chesv_aa_2stage.f"
    }

#line 251 "chesv_aa_2stage.f"
    if (*info != 0) {
#line 252 "chesv_aa_2stage.f"
	i__1 = -(*info);
#line 252 "chesv_aa_2stage.f"
	xerbla_("CHESV_AA_2STAGE", &i__1, (ftnlen)15);
#line 253 "chesv_aa_2stage.f"
	return 0;
#line 254 "chesv_aa_2stage.f"
    } else if (wquery || tquery) {
#line 255 "chesv_aa_2stage.f"
	return 0;
#line 256 "chesv_aa_2stage.f"
    }


/*     Compute the factorization A = U*T*U**H or A = L*T*L**H. */

#line 261 "chesv_aa_2stage.f"
    chetrf_aa_2stage__(uplo, n, &a[a_offset], lda, &tb[1], ltb, &ipiv[1], &
	    ipiv2[1], &work[1], lwork, info, (ftnlen)1);
#line 263 "chesv_aa_2stage.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 267 "chesv_aa_2stage.f"
	chetrs_aa_2stage__(uplo, n, nrhs, &a[a_offset], lda, &tb[1], ltb, &
		ipiv[1], &ipiv2[1], &b[b_offset], ldb, info, (ftnlen)1);

#line 270 "chesv_aa_2stage.f"
    }

#line 272 "chesv_aa_2stage.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/*     End of CHESV_AA_2STAGE */

#line 276 "chesv_aa_2stage.f"
    return 0;
} /* chesv_aa_2stage__ */

