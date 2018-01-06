#line 1 "dptsv.f"
/* dptsv.f -- translated by f2c (version 20100827).
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

#line 1 "dptsv.f"
/* > \brief <b> DPTSV computes the solution to system of linear equations A * X = B for PT matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dptsv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dptsv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dptsv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTSV( N, NRHS, D, E, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTSV computes the solution to a real system of linear equations */
/* > A*X = B, where A is an N-by-N symmetric positive definite tridiagonal */
/* > matrix, and X and B are N-by-NRHS matrices. */
/* > */
/* > A is factored as A = L*D*L**T, and the factored form of A is then */
/* > used to solve the system of equations. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A.  On exit, the n diagonal elements of the diagonal matrix */
/* >          D from the factorization A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A.  On exit, the (n-1) subdiagonal elements of the */
/* >          unit bidiagonal factor L from the L*D*L**T factorization of */
/* >          A.  (E can also be regarded as the superdiagonal of the unit */
/* >          bidiagonal factor U from the U**T*D*U factorization of A.) */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* >          > 0:  if INFO = i, the leading minor of order i is not */
/* >                positive definite, and the solution has not been */
/* >                computed.  The factorization has not been completed */
/* >                unless i = N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doublePTsolve */

/*  ===================================================================== */
/* Subroutine */ int dptsv_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpttrf_(
	    integer *, doublereal *, doublereal *, integer *), dpttrs_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);


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

#line 141 "dptsv.f"
    /* Parameter adjustments */
#line 141 "dptsv.f"
    --d__;
#line 141 "dptsv.f"
    --e;
#line 141 "dptsv.f"
    b_dim1 = *ldb;
#line 141 "dptsv.f"
    b_offset = 1 + b_dim1;
#line 141 "dptsv.f"
    b -= b_offset;
#line 141 "dptsv.f"

#line 141 "dptsv.f"
    /* Function Body */
#line 141 "dptsv.f"
    *info = 0;
#line 142 "dptsv.f"
    if (*n < 0) {
#line 143 "dptsv.f"
	*info = -1;
#line 144 "dptsv.f"
    } else if (*nrhs < 0) {
#line 145 "dptsv.f"
	*info = -2;
#line 146 "dptsv.f"
    } else if (*ldb < max(1,*n)) {
#line 147 "dptsv.f"
	*info = -6;
#line 148 "dptsv.f"
    }
#line 149 "dptsv.f"
    if (*info != 0) {
#line 150 "dptsv.f"
	i__1 = -(*info);
#line 150 "dptsv.f"
	xerbla_("DPTSV ", &i__1, (ftnlen)6);
#line 151 "dptsv.f"
	return 0;
#line 152 "dptsv.f"
    }

/*     Compute the L*D*L**T (or U**T*D*U) factorization of A. */

#line 156 "dptsv.f"
    dpttrf_(n, &d__[1], &e[1], info);
#line 157 "dptsv.f"
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

#line 161 "dptsv.f"
	dpttrs_(n, nrhs, &d__[1], &e[1], &b[b_offset], ldb, info);
#line 162 "dptsv.f"
    }
#line 163 "dptsv.f"
    return 0;

/*     End of DPTSV */

} /* dptsv_ */

