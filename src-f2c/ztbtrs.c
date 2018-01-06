#line 1 "ztbtrs.f"
/* ztbtrs.f -- translated by f2c (version 20100827).
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

#line 1 "ztbtrs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTBTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTBTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztbtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztbtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztbtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, */
/*                          LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AB( LDAB, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTBTRS solves a triangular system of the form */
/* > */
/* >    A * X = B,  A**T * X = B,  or  A**H * X = B, */
/* > */
/* > where A is a triangular band matrix of order N, and B is an */
/* > N-by-NRHS matrix.  A check is made to verify that A is nonsingular. */
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
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, integer *info, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int ztbsv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
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

#line 186 "ztbtrs.f"
    /* Parameter adjustments */
#line 186 "ztbtrs.f"
    ab_dim1 = *ldab;
#line 186 "ztbtrs.f"
    ab_offset = 1 + ab_dim1;
#line 186 "ztbtrs.f"
    ab -= ab_offset;
#line 186 "ztbtrs.f"
    b_dim1 = *ldb;
#line 186 "ztbtrs.f"
    b_offset = 1 + b_dim1;
#line 186 "ztbtrs.f"
    b -= b_offset;
#line 186 "ztbtrs.f"

#line 186 "ztbtrs.f"
    /* Function Body */
#line 186 "ztbtrs.f"
    *info = 0;
#line 187 "ztbtrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 188 "ztbtrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 189 "ztbtrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 190 "ztbtrs.f"
	*info = -1;
#line 191 "ztbtrs.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 193 "ztbtrs.f"
	*info = -2;
#line 194 "ztbtrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 195 "ztbtrs.f"
	*info = -3;
#line 196 "ztbtrs.f"
    } else if (*n < 0) {
#line 197 "ztbtrs.f"
	*info = -4;
#line 198 "ztbtrs.f"
    } else if (*kd < 0) {
#line 199 "ztbtrs.f"
	*info = -5;
#line 200 "ztbtrs.f"
    } else if (*nrhs < 0) {
#line 201 "ztbtrs.f"
	*info = -6;
#line 202 "ztbtrs.f"
    } else if (*ldab < *kd + 1) {
#line 203 "ztbtrs.f"
	*info = -8;
#line 204 "ztbtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 205 "ztbtrs.f"
	*info = -10;
#line 206 "ztbtrs.f"
    }
#line 207 "ztbtrs.f"
    if (*info != 0) {
#line 208 "ztbtrs.f"
	i__1 = -(*info);
#line 208 "ztbtrs.f"
	xerbla_("ZTBTRS", &i__1, (ftnlen)6);
#line 209 "ztbtrs.f"
	return 0;
#line 210 "ztbtrs.f"
    }

/*     Quick return if possible */

#line 214 "ztbtrs.f"
    if (*n == 0) {
#line 214 "ztbtrs.f"
	return 0;
#line 214 "ztbtrs.f"
    }

/*     Check for singularity. */

#line 219 "ztbtrs.f"
    if (nounit) {
#line 220 "ztbtrs.f"
	if (upper) {
#line 221 "ztbtrs.f"
	    i__1 = *n;
#line 221 "ztbtrs.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 222 "ztbtrs.f"
		i__2 = *kd + 1 + *info * ab_dim1;
#line 222 "ztbtrs.f"
		if (ab[i__2].r == 0. && ab[i__2].i == 0.) {
#line 222 "ztbtrs.f"
		    return 0;
#line 222 "ztbtrs.f"
		}
#line 224 "ztbtrs.f"
/* L10: */
#line 224 "ztbtrs.f"
	    }
#line 225 "ztbtrs.f"
	} else {
#line 226 "ztbtrs.f"
	    i__1 = *n;
#line 226 "ztbtrs.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 227 "ztbtrs.f"
		i__2 = *info * ab_dim1 + 1;
#line 227 "ztbtrs.f"
		if (ab[i__2].r == 0. && ab[i__2].i == 0.) {
#line 227 "ztbtrs.f"
		    return 0;
#line 227 "ztbtrs.f"
		}
#line 229 "ztbtrs.f"
/* L20: */
#line 229 "ztbtrs.f"
	    }
#line 230 "ztbtrs.f"
	}
#line 231 "ztbtrs.f"
    }
#line 232 "ztbtrs.f"
    *info = 0;

/*     Solve A * X = B,  A**T * X = B,  or  A**H * X = B. */

#line 236 "ztbtrs.f"
    i__1 = *nrhs;
#line 236 "ztbtrs.f"
    for (j = 1; j <= i__1; ++j) {
#line 237 "ztbtrs.f"
	ztbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &b[j * b_dim1 
		+ 1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 238 "ztbtrs.f"
/* L30: */
#line 238 "ztbtrs.f"
    }

#line 240 "ztbtrs.f"
    return 0;

/*     End of ZTBTRS */

} /* ztbtrs_ */

