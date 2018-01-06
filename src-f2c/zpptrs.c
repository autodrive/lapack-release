#line 1 "zpptrs.f"
/* zpptrs.f -- translated by f2c (version 20100827).
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

#line 1 "zpptrs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZPPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPTRS solves a system of linear equations A*X = B with a Hermitian */
/* > positive definite matrix A in packed storage using the Cholesky */
/* > factorization A = U**H * U or A = L * L**H computed by ZPPTRF. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H * U or A = L * L**H, packed columnwise in a linear */
/* >          array.  The j-th column of U or L is stored in the array AP */
/* >          as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
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
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1;

    /* Local variables */
    static integer i__;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int ztpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

#line 144 "zpptrs.f"
    /* Parameter adjustments */
#line 144 "zpptrs.f"
    --ap;
#line 144 "zpptrs.f"
    b_dim1 = *ldb;
#line 144 "zpptrs.f"
    b_offset = 1 + b_dim1;
#line 144 "zpptrs.f"
    b -= b_offset;
#line 144 "zpptrs.f"

#line 144 "zpptrs.f"
    /* Function Body */
#line 144 "zpptrs.f"
    *info = 0;
#line 145 "zpptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 146 "zpptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 147 "zpptrs.f"
	*info = -1;
#line 148 "zpptrs.f"
    } else if (*n < 0) {
#line 149 "zpptrs.f"
	*info = -2;
#line 150 "zpptrs.f"
    } else if (*nrhs < 0) {
#line 151 "zpptrs.f"
	*info = -3;
#line 152 "zpptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 153 "zpptrs.f"
	*info = -6;
#line 154 "zpptrs.f"
    }
#line 155 "zpptrs.f"
    if (*info != 0) {
#line 156 "zpptrs.f"
	i__1 = -(*info);
#line 156 "zpptrs.f"
	xerbla_("ZPPTRS", &i__1, (ftnlen)6);
#line 157 "zpptrs.f"
	return 0;
#line 158 "zpptrs.f"
    }

/*     Quick return if possible */

#line 162 "zpptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 162 "zpptrs.f"
	return 0;
#line 162 "zpptrs.f"
    }

#line 165 "zpptrs.f"
    if (upper) {

/*        Solve A*X = B where A = U**H * U. */

#line 169 "zpptrs.f"
	i__1 = *nrhs;
#line 169 "zpptrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U**H *X = B, overwriting B with X. */

#line 173 "zpptrs.f"
	    ztpsv_("Upper", "Conjugate transpose", "Non-unit", n, &ap[1], &b[
		    i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)19, (ftnlen)
		    8);

/*           Solve U*X = B, overwriting B with X. */

#line 178 "zpptrs.f"
	    ztpsv_("Upper", "No transpose", "Non-unit", n, &ap[1], &b[i__ * 
		    b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 180 "zpptrs.f"
/* L10: */
#line 180 "zpptrs.f"
	}
#line 181 "zpptrs.f"
    } else {

/*        Solve A*X = B where A = L * L**H. */

#line 185 "zpptrs.f"
	i__1 = *nrhs;
#line 185 "zpptrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve L*Y = B, overwriting B with X. */

#line 189 "zpptrs.f"
	    ztpsv_("Lower", "No transpose", "Non-unit", n, &ap[1], &b[i__ * 
		    b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           Solve L**H *X = Y, overwriting B with X. */

#line 194 "zpptrs.f"
	    ztpsv_("Lower", "Conjugate transpose", "Non-unit", n, &ap[1], &b[
		    i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)19, (ftnlen)
		    8);
#line 196 "zpptrs.f"
/* L20: */
#line 196 "zpptrs.f"
	}
#line 197 "zpptrs.f"
    }

#line 199 "zpptrs.f"
    return 0;

/*     End of ZPPTRS */

} /* zpptrs_ */

