#line 1 "cgesv.f"
/* cgesv.f -- translated by f2c (version 20100827).
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

#line 1 "cgesv.f"
/* > \brief <b> CGESV computes the solution to system of linear equations A * X = B for GE matrices (simple dr
iver) </b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESV computes the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* > The LU decomposition with partial pivoting and row interchanges is */
/* > used to factor A as */
/* >    A = P * L * U, */
/* > where P is a permutation matrix, L is unit lower triangular, and U is */
/* > upper triangular.  The factored form of A is then used to solve the */
/* > system of equations A * X = B. */
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
/* >          On entry, the N-by-N coefficient matrix A. */
/* >          On exit, the factors L and U from the factorization */
/* >          A = P*L*U; the unit diagonal elements of L are not stored. */
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
/* >          The pivot indices that define the permutation matrix P; */
/* >          row i of the matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >          On entry, the N-by-NRHS matrix of right hand side matrix B. */
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
/* >                singular, so the solution could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complexGEsolve */

/*  ===================================================================== */
/* Subroutine */ int cgesv_(integer *n, integer *nrhs, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int cgetrf_(integer *, integer *, doublecomplex *,
	     integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen), cgetrs_(char *, integer *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);


/*  -- LAPACK driver routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 150 "cgesv.f"
    /* Parameter adjustments */
#line 150 "cgesv.f"
    a_dim1 = *lda;
#line 150 "cgesv.f"
    a_offset = 1 + a_dim1;
#line 150 "cgesv.f"
    a -= a_offset;
#line 150 "cgesv.f"
    --ipiv;
#line 150 "cgesv.f"
    b_dim1 = *ldb;
#line 150 "cgesv.f"
    b_offset = 1 + b_dim1;
#line 150 "cgesv.f"
    b -= b_offset;
#line 150 "cgesv.f"

#line 150 "cgesv.f"
    /* Function Body */
#line 150 "cgesv.f"
    *info = 0;
#line 151 "cgesv.f"
    if (*n < 0) {
#line 152 "cgesv.f"
	*info = -1;
#line 153 "cgesv.f"
    } else if (*nrhs < 0) {
#line 154 "cgesv.f"
	*info = -2;
#line 155 "cgesv.f"
    } else if (*lda < max(1,*n)) {
#line 156 "cgesv.f"
	*info = -4;
#line 157 "cgesv.f"
    } else if (*ldb < max(1,*n)) {
#line 158 "cgesv.f"
	*info = -7;
#line 159 "cgesv.f"
    }
#line 160 "cgesv.f"
    if (*info != 0) {
#line 161 "cgesv.f"
	i__1 = -(*info);
#line 161 "cgesv.f"
	xerbla_("CGESV ", &i__1, (ftnlen)6);
#line 162 "cgesv.f"
	return 0;
#line 163 "cgesv.f"
    }

/*     Compute the LU factorization of A. */

#line 167 "cgesv.f"
    cgetrf_(n, n, &a[a_offset], lda, &ipiv[1], info);
#line 168 "cgesv.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 172 "cgesv.f"
	cgetrs_("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &b[
		b_offset], ldb, info, (ftnlen)12);
#line 174 "cgesv.f"
    }
#line 175 "cgesv.f"
    return 0;

/*     End of CGESV */

} /* cgesv_ */

