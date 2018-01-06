#line 1 "cpbtrs.f"
/* cpbtrs.f -- translated by f2c (version 20100827).
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

#line 1 "cpbtrs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CPBTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AB( LDAB, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBTRS solves a system of linear equations A*X = B with a Hermitian */
/* > positive definite band matrix A using the Cholesky factorization */
/* > A = U**H*U or A = L*L**H computed by CPBTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangular factor stored in AB; */
/* >          = 'L':  Lower triangular factor stored in AB. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H of the band matrix A, stored in the */
/* >          first KD+1 rows of the array.  The j-th column of U or L is */
/* >          stored in the j-th column of the array AB as follows: */
/* >          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static integer j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctbsv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 157 "cpbtrs.f"
    /* Parameter adjustments */
#line 157 "cpbtrs.f"
    ab_dim1 = *ldab;
#line 157 "cpbtrs.f"
    ab_offset = 1 + ab_dim1;
#line 157 "cpbtrs.f"
    ab -= ab_offset;
#line 157 "cpbtrs.f"
    b_dim1 = *ldb;
#line 157 "cpbtrs.f"
    b_offset = 1 + b_dim1;
#line 157 "cpbtrs.f"
    b -= b_offset;
#line 157 "cpbtrs.f"

#line 157 "cpbtrs.f"
    /* Function Body */
#line 157 "cpbtrs.f"
    *info = 0;
#line 158 "cpbtrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 159 "cpbtrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 160 "cpbtrs.f"
	*info = -1;
#line 161 "cpbtrs.f"
    } else if (*n < 0) {
#line 162 "cpbtrs.f"
	*info = -2;
#line 163 "cpbtrs.f"
    } else if (*kd < 0) {
#line 164 "cpbtrs.f"
	*info = -3;
#line 165 "cpbtrs.f"
    } else if (*nrhs < 0) {
#line 166 "cpbtrs.f"
	*info = -4;
#line 167 "cpbtrs.f"
    } else if (*ldab < *kd + 1) {
#line 168 "cpbtrs.f"
	*info = -6;
#line 169 "cpbtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 170 "cpbtrs.f"
	*info = -8;
#line 171 "cpbtrs.f"
    }
#line 172 "cpbtrs.f"
    if (*info != 0) {
#line 173 "cpbtrs.f"
	i__1 = -(*info);
#line 173 "cpbtrs.f"
	xerbla_("CPBTRS", &i__1, (ftnlen)6);
#line 174 "cpbtrs.f"
	return 0;
#line 175 "cpbtrs.f"
    }

/*     Quick return if possible */

#line 179 "cpbtrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 179 "cpbtrs.f"
	return 0;
#line 179 "cpbtrs.f"
    }

#line 182 "cpbtrs.f"
    if (upper) {

/*        Solve A*X = B where A = U**H *U. */

#line 186 "cpbtrs.f"
	i__1 = *nrhs;
#line 186 "cpbtrs.f"
	for (j = 1; j <= i__1; ++j) {

/*           Solve U**H *X = B, overwriting B with X. */

#line 190 "cpbtrs.f"
	    ctbsv_("Upper", "Conjugate transpose", "Non-unit", n, kd, &ab[
		    ab_offset], ldab, &b[j * b_dim1 + 1], &c__1, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

/*           Solve U*X = B, overwriting B with X. */

#line 195 "cpbtrs.f"
	    ctbsv_("Upper", "No transpose", "Non-unit", n, kd, &ab[ab_offset],
		     ldab, &b[j * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);
#line 197 "cpbtrs.f"
/* L10: */
#line 197 "cpbtrs.f"
	}
#line 198 "cpbtrs.f"
    } else {

/*        Solve A*X = B where A = L*L**H. */

#line 202 "cpbtrs.f"
	i__1 = *nrhs;
#line 202 "cpbtrs.f"
	for (j = 1; j <= i__1; ++j) {

/*           Solve L*X = B, overwriting B with X. */

#line 206 "cpbtrs.f"
	    ctbsv_("Lower", "No transpose", "Non-unit", n, kd, &ab[ab_offset],
		     ldab, &b[j * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);

/*           Solve L**H *X = B, overwriting B with X. */

#line 211 "cpbtrs.f"
	    ctbsv_("Lower", "Conjugate transpose", "Non-unit", n, kd, &ab[
		    ab_offset], ldab, &b[j * b_dim1 + 1], &c__1, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);
#line 213 "cpbtrs.f"
/* L20: */
#line 213 "cpbtrs.f"
	}
#line 214 "cpbtrs.f"
    }

#line 216 "cpbtrs.f"
    return 0;

/*     End of CPBTRS */

} /* cpbtrs_ */

