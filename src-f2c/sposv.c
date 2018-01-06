#line 1 "sposv.f"
/* sposv.f -- translated by f2c (version 20100827).
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

#line 1 "sposv.f"
/* > \brief <b> SPOSV computes the solution to system of linear equations A * X = B for PO matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPOSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sposv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sposv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sposv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPOSV computes the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric positive definite matrix and X and B */
/* > are N-by-NRHS matrices. */
/* > */
/* > The Cholesky decomposition is used to factor A as */
/* >    A = U**T* U,  if UPLO = 'U', or */
/* >    A = L * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is a lower triangular */
/* > matrix.  The factored form of A is then used to solve the system of */
/* > equations A * X = B. */
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
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T. */
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
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i of A is not */
/* >                positive definite, so the factorization could not be */
/* >                completed, and the solution has not been computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realPOsolve */

/*  ===================================================================== */
/* Subroutine */ int sposv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), spotrf_(
	    char *, integer *, doublereal *, integer *, integer *, ftnlen), 
	    spotrs_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 162 "sposv.f"
    /* Parameter adjustments */
#line 162 "sposv.f"
    a_dim1 = *lda;
#line 162 "sposv.f"
    a_offset = 1 + a_dim1;
#line 162 "sposv.f"
    a -= a_offset;
#line 162 "sposv.f"
    b_dim1 = *ldb;
#line 162 "sposv.f"
    b_offset = 1 + b_dim1;
#line 162 "sposv.f"
    b -= b_offset;
#line 162 "sposv.f"

#line 162 "sposv.f"
    /* Function Body */
#line 162 "sposv.f"
    *info = 0;
#line 163 "sposv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 164 "sposv.f"
	*info = -1;
#line 165 "sposv.f"
    } else if (*n < 0) {
#line 166 "sposv.f"
	*info = -2;
#line 167 "sposv.f"
    } else if (*nrhs < 0) {
#line 168 "sposv.f"
	*info = -3;
#line 169 "sposv.f"
    } else if (*lda < max(1,*n)) {
#line 170 "sposv.f"
	*info = -5;
#line 171 "sposv.f"
    } else if (*ldb < max(1,*n)) {
#line 172 "sposv.f"
	*info = -7;
#line 173 "sposv.f"
    }
#line 174 "sposv.f"
    if (*info != 0) {
#line 175 "sposv.f"
	i__1 = -(*info);
#line 175 "sposv.f"
	xerbla_("SPOSV ", &i__1, (ftnlen)6);
#line 176 "sposv.f"
	return 0;
#line 177 "sposv.f"
    }

/*     Compute the Cholesky factorization A = U**T*U or A = L*L**T. */

#line 181 "sposv.f"
    spotrf_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
#line 182 "sposv.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 186 "sposv.f"
	spotrs_(uplo, n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info, (
		ftnlen)1);

#line 188 "sposv.f"
    }
#line 189 "sposv.f"
    return 0;

/*     End of SPOSV */

} /* sposv_ */

