#line 1 "stptrs.f"
/* stptrs.f -- translated by f2c (version 20100827).
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

#line 1 "stptrs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b STPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPTRS solves a triangular system of the form */
/* > */
/* >    A * X = B  or  A**T * X = B, */
/* > */
/* > where A is a triangular matrix of order N stored in packed format, */
/* > and B is an N-by-NRHS matrix.  A check is made to verify that A is */
/* > nonsingular. */
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
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
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
/* Subroutine */ int stptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1;

    /* Local variables */
    static integer j, jc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int stpsv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
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

#line 170 "stptrs.f"
    /* Parameter adjustments */
#line 170 "stptrs.f"
    --ap;
#line 170 "stptrs.f"
    b_dim1 = *ldb;
#line 170 "stptrs.f"
    b_offset = 1 + b_dim1;
#line 170 "stptrs.f"
    b -= b_offset;
#line 170 "stptrs.f"

#line 170 "stptrs.f"
    /* Function Body */
#line 170 "stptrs.f"
    *info = 0;
#line 171 "stptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 172 "stptrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 173 "stptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "stptrs.f"
	*info = -1;
#line 175 "stptrs.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 177 "stptrs.f"
	*info = -2;
#line 178 "stptrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 179 "stptrs.f"
	*info = -3;
#line 180 "stptrs.f"
    } else if (*n < 0) {
#line 181 "stptrs.f"
	*info = -4;
#line 182 "stptrs.f"
    } else if (*nrhs < 0) {
#line 183 "stptrs.f"
	*info = -5;
#line 184 "stptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 185 "stptrs.f"
	*info = -8;
#line 186 "stptrs.f"
    }
#line 187 "stptrs.f"
    if (*info != 0) {
#line 188 "stptrs.f"
	i__1 = -(*info);
#line 188 "stptrs.f"
	xerbla_("STPTRS", &i__1, (ftnlen)6);
#line 189 "stptrs.f"
	return 0;
#line 190 "stptrs.f"
    }

/*     Quick return if possible */

#line 194 "stptrs.f"
    if (*n == 0) {
#line 194 "stptrs.f"
	return 0;
#line 194 "stptrs.f"
    }

/*     Check for singularity. */

#line 199 "stptrs.f"
    if (nounit) {
#line 200 "stptrs.f"
	if (upper) {
#line 201 "stptrs.f"
	    jc = 1;
#line 202 "stptrs.f"
	    i__1 = *n;
#line 202 "stptrs.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 203 "stptrs.f"
		if (ap[jc + *info - 1] == 0.) {
#line 203 "stptrs.f"
		    return 0;
#line 203 "stptrs.f"
		}
#line 205 "stptrs.f"
		jc += *info;
#line 206 "stptrs.f"
/* L10: */
#line 206 "stptrs.f"
	    }
#line 207 "stptrs.f"
	} else {
#line 208 "stptrs.f"
	    jc = 1;
#line 209 "stptrs.f"
	    i__1 = *n;
#line 209 "stptrs.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 210 "stptrs.f"
		if (ap[jc] == 0.) {
#line 210 "stptrs.f"
		    return 0;
#line 210 "stptrs.f"
		}
#line 212 "stptrs.f"
		jc = jc + *n - *info + 1;
#line 213 "stptrs.f"
/* L20: */
#line 213 "stptrs.f"
	    }
#line 214 "stptrs.f"
	}
#line 215 "stptrs.f"
    }
#line 216 "stptrs.f"
    *info = 0;

/*     Solve A * x = b  or  A**T * x = b. */

#line 220 "stptrs.f"
    i__1 = *nrhs;
#line 220 "stptrs.f"
    for (j = 1; j <= i__1; ++j) {
#line 221 "stptrs.f"
	stpsv_(uplo, trans, diag, n, &ap[1], &b[j * b_dim1 + 1], &c__1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 222 "stptrs.f"
/* L30: */
#line 222 "stptrs.f"
    }

#line 224 "stptrs.f"
    return 0;

/*     End of STPTRS */

} /* stptrs_ */

