#line 1 "zgbsv.f"
/* zgbsv.f -- translated by f2c (version 20100827).
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

#line 1 "zgbsv.f"
/* > \brief <b> ZGBSV computes the solution to system of linear equations A * X = B for GB matrices</b> (simpl
e driver) */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbsv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbsv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbsv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AB( LDAB, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGBSV computes the solution to a complex system of linear equations */
/* > A * X = B, where A is a band matrix of order N with KL subdiagonals */
/* > and KU superdiagonals, and X and B are N-by-NRHS matrices. */
/* > */
/* > The LU decomposition with partial pivoting and row interchanges is */
/* > used to factor A as A = L * U, where L is a product of permutation */
/* > and unit lower triangular matrices with KL subdiagonals, and U is */
/* > upper triangular with KL+KU superdiagonals.  The factored form of A */
/* > is then used to solve the system of equations A * X = B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of linear equations, i.e., the order of the */
/* >          matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On entry, the matrix A in band storage, in rows KL+1 to */
/* >          2*KL+KU+1; rows 1 to KL of the array need not be set. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL) */
/* >          On exit, details of the factorization: U is stored as an */
/* >          upper triangular band matrix with KL+KU superdiagonals in */
/* >          rows 1 to KL+KU+1, and the multipliers used during the */
/* >          factorization are stored in rows KL+KU+2 to 2*KL+KU+1. */
/* >          See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices that define the permutation matrix P; */
/* >          row i of the matrix was interchanged with row IPIV(i). */
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
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization */
/* >                has been completed, but the factor U is exactly */
/* >                singular, and the solution has not been computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GBsolve */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  M = N = 6, KL = 2, KU = 1: */
/* > */
/* >  On entry:                       On exit: */
/* > */
/* >      *    *    *    +    +    +       *    *    *   u14  u25  u36 */
/* >      *    *    +    +    +    +       *    *   u13  u24  u35  u46 */
/* >      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 */
/* >     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 */
/* >     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   * */
/* >     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    * */
/* > */
/* >  Array elements marked * are not used by the routine; elements marked */
/* >  + need not be set on entry, but are required by the routine to store */
/* >  elements of U because of fill-in resulting from the row interchanges. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *
	b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zgbtrf_(
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, integer *, integer *), zgbtrs_(char *, integer *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 190 "zgbsv.f"
    /* Parameter adjustments */
#line 190 "zgbsv.f"
    ab_dim1 = *ldab;
#line 190 "zgbsv.f"
    ab_offset = 1 + ab_dim1;
#line 190 "zgbsv.f"
    ab -= ab_offset;
#line 190 "zgbsv.f"
    --ipiv;
#line 190 "zgbsv.f"
    b_dim1 = *ldb;
#line 190 "zgbsv.f"
    b_offset = 1 + b_dim1;
#line 190 "zgbsv.f"
    b -= b_offset;
#line 190 "zgbsv.f"

#line 190 "zgbsv.f"
    /* Function Body */
#line 190 "zgbsv.f"
    *info = 0;
#line 191 "zgbsv.f"
    if (*n < 0) {
#line 192 "zgbsv.f"
	*info = -1;
#line 193 "zgbsv.f"
    } else if (*kl < 0) {
#line 194 "zgbsv.f"
	*info = -2;
#line 195 "zgbsv.f"
    } else if (*ku < 0) {
#line 196 "zgbsv.f"
	*info = -3;
#line 197 "zgbsv.f"
    } else if (*nrhs < 0) {
#line 198 "zgbsv.f"
	*info = -4;
#line 199 "zgbsv.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 200 "zgbsv.f"
	*info = -6;
#line 201 "zgbsv.f"
    } else if (*ldb < max(*n,1)) {
#line 202 "zgbsv.f"
	*info = -9;
#line 203 "zgbsv.f"
    }
#line 204 "zgbsv.f"
    if (*info != 0) {
#line 205 "zgbsv.f"
	i__1 = -(*info);
#line 205 "zgbsv.f"
	xerbla_("ZGBSV ", &i__1, (ftnlen)6);
#line 206 "zgbsv.f"
	return 0;
#line 207 "zgbsv.f"
    }

/*     Compute the LU factorization of the band matrix A. */

#line 211 "zgbsv.f"
    zgbtrf_(n, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
#line 212 "zgbsv.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 216 "zgbsv.f"
	zgbtrs_("No transpose", n, kl, ku, nrhs, &ab[ab_offset], ldab, &ipiv[
		1], &b[b_offset], ldb, info, (ftnlen)12);
#line 218 "zgbsv.f"
    }
#line 219 "zgbsv.f"
    return 0;

/*     End of ZGBSV */

} /* zgbsv_ */

