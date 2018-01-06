#line 1 "zhesv_aa_2stage.f"
/* zhesv_aa_2stage.f -- translated by f2c (version 20100827).
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

#line 1 "zhesv_aa_2stage.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief <b> ZHESV_AA_2STAGE computes the solution to system of linear equations A * X = B for HE matrices
</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHESV_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesv_a
a_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesv_a
a_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesv_a
a_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*      SUBROUTINE ZHESV_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, */
/*                                  IPIV, IPIV2, B, LDB, WORK, LWORK, */
/*                                  INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, NRHS, LDA, LTB, LDB, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IPIV2( * ) */
/*       COMPLEX*16         A( LDA, * ), TB( * ), B( LDB, *), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHESV_AA_2STAGE computes the solution to a complex system of */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          TB is COMPLEX*16 array, dimension (LTB) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* >          WORK is COMPLEX*16 workspace of size LWORK */
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

/* > \ingroup complex16HEsolve */

/*  ===================================================================== */
/* Subroutine */ int zhesv_aa_2stage__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *tb, integer *ltb, 
	integer *ipiv, integer *ipiv2, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int zhetrf_aa_2stage__(char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen), 
	    zhetrs_aa_2stage__(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical tquery, wquery;


/*  -- LAPACK driver routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*     .. Parameters .. */

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

#line 231 "zhesv_aa_2stage.f"
    /* Parameter adjustments */
#line 231 "zhesv_aa_2stage.f"
    a_dim1 = *lda;
#line 231 "zhesv_aa_2stage.f"
    a_offset = 1 + a_dim1;
#line 231 "zhesv_aa_2stage.f"
    a -= a_offset;
#line 231 "zhesv_aa_2stage.f"
    --tb;
#line 231 "zhesv_aa_2stage.f"
    --ipiv;
#line 231 "zhesv_aa_2stage.f"
    --ipiv2;
#line 231 "zhesv_aa_2stage.f"
    b_dim1 = *ldb;
#line 231 "zhesv_aa_2stage.f"
    b_offset = 1 + b_dim1;
#line 231 "zhesv_aa_2stage.f"
    b -= b_offset;
#line 231 "zhesv_aa_2stage.f"
    --work;
#line 231 "zhesv_aa_2stage.f"

#line 231 "zhesv_aa_2stage.f"
    /* Function Body */
#line 231 "zhesv_aa_2stage.f"
    *info = 0;
#line 232 "zhesv_aa_2stage.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 233 "zhesv_aa_2stage.f"
    wquery = *lwork == -1;
#line 234 "zhesv_aa_2stage.f"
    tquery = *ltb == -1;
#line 235 "zhesv_aa_2stage.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 236 "zhesv_aa_2stage.f"
	*info = -1;
#line 237 "zhesv_aa_2stage.f"
    } else if (*n < 0) {
#line 238 "zhesv_aa_2stage.f"
	*info = -2;
#line 239 "zhesv_aa_2stage.f"
    } else if (*nrhs < 0) {
#line 240 "zhesv_aa_2stage.f"
	*info = -3;
#line 241 "zhesv_aa_2stage.f"
    } else if (*lda < max(1,*n)) {
#line 242 "zhesv_aa_2stage.f"
	*info = -5;
#line 243 "zhesv_aa_2stage.f"
    } else if (*ldb < max(1,*n)) {
#line 244 "zhesv_aa_2stage.f"
	*info = -11;
#line 245 "zhesv_aa_2stage.f"
    }

#line 247 "zhesv_aa_2stage.f"
    if (*info == 0) {
#line 248 "zhesv_aa_2stage.f"
	zhetrf_aa_2stage__(uplo, n, &a[a_offset], lda, &tb[1], &c_n1, &ipiv[1]
		, &ipiv2[1], &work[1], &c_n1, info, (ftnlen)1);
#line 250 "zhesv_aa_2stage.f"
	lwkopt = (integer) work[1].r;
#line 251 "zhesv_aa_2stage.f"
	if (*ltb < (integer) tb[1].r && ! tquery) {
#line 252 "zhesv_aa_2stage.f"
	    *info = -7;
#line 253 "zhesv_aa_2stage.f"
	} else if (*lwork < lwkopt && ! wquery) {
#line 254 "zhesv_aa_2stage.f"
	    *info = -13;
#line 255 "zhesv_aa_2stage.f"
	}
#line 256 "zhesv_aa_2stage.f"
    }

#line 258 "zhesv_aa_2stage.f"
    if (*info != 0) {
#line 259 "zhesv_aa_2stage.f"
	i__1 = -(*info);
#line 259 "zhesv_aa_2stage.f"
	xerbla_("ZHESV_AA_2STAGE", &i__1, (ftnlen)15);
#line 260 "zhesv_aa_2stage.f"
	return 0;
#line 261 "zhesv_aa_2stage.f"
    } else if (wquery || tquery) {
#line 262 "zhesv_aa_2stage.f"
	return 0;
#line 263 "zhesv_aa_2stage.f"
    }

/*     Compute the factorization A = U*T*U**H or A = L*T*L**H. */

#line 267 "zhesv_aa_2stage.f"
    zhetrf_aa_2stage__(uplo, n, &a[a_offset], lda, &tb[1], ltb, &ipiv[1], &
	    ipiv2[1], &work[1], lwork, info, (ftnlen)1);
#line 269 "zhesv_aa_2stage.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 273 "zhesv_aa_2stage.f"
	zhetrs_aa_2stage__(uplo, n, nrhs, &a[a_offset], lda, &tb[1], ltb, &
		ipiv[1], &ipiv2[1], &b[b_offset], ldb, info, (ftnlen)1);

#line 276 "zhesv_aa_2stage.f"
    }

#line 278 "zhesv_aa_2stage.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 280 "zhesv_aa_2stage.f"
    return 0;

/*     End of ZHESV_AA_2STAGE */

} /* zhesv_aa_2stage__ */

