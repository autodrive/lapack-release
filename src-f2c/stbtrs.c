#line 1 "stbtrs.f"
/* stbtrs.f -- translated by f2c (version 20100827).
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

#line 1 "stbtrs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b STBTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STBTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, */
/*                          LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STBTRS solves a triangular system of the form */
/* > */
/* >    A * X = B  or  A**T * X = B, */
/* > */
/* > where A is a triangular band matrix of order N, and B is an */
/* > N-by NRHS matrix.  A check is made to verify that A is nonsingular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form the system of equations: */
/* >          = 'N':  A * X = B  (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
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
/* >          The number of superdiagonals or subdiagonals of the */
/* >          triangular band matrix A.  KD >= 0. */
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
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          The upper or lower triangular band matrix A, stored in the */
/* >          first kd+1 rows of AB.  The j-th column of A is stored */
/* >          in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* >          If DIAG = 'U', the diagonal elements of A are not referenced */
/* >          and are assumed to be 1. */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, if INFO = 0, the solution matrix X. */
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
/* >          > 0:  if INFO = i, the i-th diagonal element of A is zero, */
/* >                indicating that the matrix is singular and the */
/* >                solutions X have not been computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int stbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	*b, integer *ldb, integer *info, ftnlen uplo_len, ftnlen trans_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static integer j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int stbsv_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical nounit;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
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

#line 186 "stbtrs.f"
    /* Parameter adjustments */
#line 186 "stbtrs.f"
    ab_dim1 = *ldab;
#line 186 "stbtrs.f"
    ab_offset = 1 + ab_dim1;
#line 186 "stbtrs.f"
    ab -= ab_offset;
#line 186 "stbtrs.f"
    b_dim1 = *ldb;
#line 186 "stbtrs.f"
    b_offset = 1 + b_dim1;
#line 186 "stbtrs.f"
    b -= b_offset;
#line 186 "stbtrs.f"

#line 186 "stbtrs.f"
    /* Function Body */
#line 186 "stbtrs.f"
    *info = 0;
#line 187 "stbtrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 188 "stbtrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 189 "stbtrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 190 "stbtrs.f"
	*info = -1;
#line 191 "stbtrs.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 193 "stbtrs.f"
	*info = -2;
#line 194 "stbtrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 195 "stbtrs.f"
	*info = -3;
#line 196 "stbtrs.f"
    } else if (*n < 0) {
#line 197 "stbtrs.f"
	*info = -4;
#line 198 "stbtrs.f"
    } else if (*kd < 0) {
#line 199 "stbtrs.f"
	*info = -5;
#line 200 "stbtrs.f"
    } else if (*nrhs < 0) {
#line 201 "stbtrs.f"
	*info = -6;
#line 202 "stbtrs.f"
    } else if (*ldab < *kd + 1) {
#line 203 "stbtrs.f"
	*info = -8;
#line 204 "stbtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 205 "stbtrs.f"
	*info = -10;
#line 206 "stbtrs.f"
    }
#line 207 "stbtrs.f"
    if (*info != 0) {
#line 208 "stbtrs.f"
	i__1 = -(*info);
#line 208 "stbtrs.f"
	xerbla_("STBTRS", &i__1, (ftnlen)6);
#line 209 "stbtrs.f"
	return 0;
#line 210 "stbtrs.f"
    }

/*     Quick return if possible */

#line 214 "stbtrs.f"
    if (*n == 0) {
#line 214 "stbtrs.f"
	return 0;
#line 214 "stbtrs.f"
    }

/*     Check for singularity. */

#line 219 "stbtrs.f"
    if (nounit) {
#line 220 "stbtrs.f"
	if (upper) {
#line 221 "stbtrs.f"
	    i__1 = *n;
#line 221 "stbtrs.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 222 "stbtrs.f"
		if (ab[*kd + 1 + *info * ab_dim1] == 0.) {
#line 222 "stbtrs.f"
		    return 0;
#line 222 "stbtrs.f"
		}
#line 224 "stbtrs.f"
/* L10: */
#line 224 "stbtrs.f"
	    }
#line 225 "stbtrs.f"
	} else {
#line 226 "stbtrs.f"
	    i__1 = *n;
#line 226 "stbtrs.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 227 "stbtrs.f"
		if (ab[*info * ab_dim1 + 1] == 0.) {
#line 227 "stbtrs.f"
		    return 0;
#line 227 "stbtrs.f"
		}
#line 229 "stbtrs.f"
/* L20: */
#line 229 "stbtrs.f"
	    }
#line 230 "stbtrs.f"
	}
#line 231 "stbtrs.f"
    }
#line 232 "stbtrs.f"
    *info = 0;

/*     Solve A * X = B  or  A**T * X = B. */

#line 236 "stbtrs.f"
    i__1 = *nrhs;
#line 236 "stbtrs.f"
    for (j = 1; j <= i__1; ++j) {
#line 237 "stbtrs.f"
	stbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &b[j * b_dim1 
		+ 1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 238 "stbtrs.f"
/* L30: */
#line 238 "stbtrs.f"
    }

#line 240 "stbtrs.f"
    return 0;

/*     End of STBTRS */

} /* stbtrs_ */

