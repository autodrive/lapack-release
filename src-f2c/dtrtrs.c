#line 1 "dtrtrs.f"
/* dtrtrs.f -- translated by f2c (version 20100827).
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

#line 1 "dtrtrs.f"
/* Table of constant values */

static doublereal c_b12 = 1.;

/* > \brief \b DTRTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRTRS solves a triangular system of the form */
/* > */
/* >    A * X = B  or  A**T * X = B, */
/* > */
/* > where A is a triangular matrix of order N, and B is an N-by-NRHS */
/* > matrix.  A check is made to verify that A is nonsingular. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The triangular matrix A.  If UPLO = 'U', the leading N-by-N */
/* >          upper triangular part of the array A contains the upper */
/* >          triangular matrix, and the strictly lower triangular part of */
/* >          A is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/* >          triangular part of the array A contains the lower triangular */
/* >          matrix, and the strictly upper triangular part of A is not */
/* >          referenced.  If DIAG = 'U', the diagonal elements of A are */
/* >          also not referenced and are assumed to be 1. */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, the i-th diagonal element of A is zero, */
/* >               indicating that the matrix is singular and the solutions */
/* >               X have not been computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, integer *info, ftnlen uplo_len, ftnlen trans_len, ftnlen 
	diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
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

#line 179 "dtrtrs.f"
    /* Parameter adjustments */
#line 179 "dtrtrs.f"
    a_dim1 = *lda;
#line 179 "dtrtrs.f"
    a_offset = 1 + a_dim1;
#line 179 "dtrtrs.f"
    a -= a_offset;
#line 179 "dtrtrs.f"
    b_dim1 = *ldb;
#line 179 "dtrtrs.f"
    b_offset = 1 + b_dim1;
#line 179 "dtrtrs.f"
    b -= b_offset;
#line 179 "dtrtrs.f"

#line 179 "dtrtrs.f"
    /* Function Body */
#line 179 "dtrtrs.f"
    *info = 0;
#line 180 "dtrtrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 181 "dtrtrs.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 182 "dtrtrs.f"
	*info = -1;
#line 183 "dtrtrs.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 185 "dtrtrs.f"
	*info = -2;
#line 186 "dtrtrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 187 "dtrtrs.f"
	*info = -3;
#line 188 "dtrtrs.f"
    } else if (*n < 0) {
#line 189 "dtrtrs.f"
	*info = -4;
#line 190 "dtrtrs.f"
    } else if (*nrhs < 0) {
#line 191 "dtrtrs.f"
	*info = -5;
#line 192 "dtrtrs.f"
    } else if (*lda < max(1,*n)) {
#line 193 "dtrtrs.f"
	*info = -7;
#line 194 "dtrtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 195 "dtrtrs.f"
	*info = -9;
#line 196 "dtrtrs.f"
    }
#line 197 "dtrtrs.f"
    if (*info != 0) {
#line 198 "dtrtrs.f"
	i__1 = -(*info);
#line 198 "dtrtrs.f"
	xerbla_("DTRTRS", &i__1, (ftnlen)6);
#line 199 "dtrtrs.f"
	return 0;
#line 200 "dtrtrs.f"
    }

/*     Quick return if possible */

#line 204 "dtrtrs.f"
    if (*n == 0) {
#line 204 "dtrtrs.f"
	return 0;
#line 204 "dtrtrs.f"
    }

/*     Check for singularity. */

#line 209 "dtrtrs.f"
    if (nounit) {
#line 210 "dtrtrs.f"
	i__1 = *n;
#line 210 "dtrtrs.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 211 "dtrtrs.f"
	    if (a[*info + *info * a_dim1] == 0.) {
#line 211 "dtrtrs.f"
		return 0;
#line 211 "dtrtrs.f"
	    }
#line 213 "dtrtrs.f"
/* L10: */
#line 213 "dtrtrs.f"
	}
#line 214 "dtrtrs.f"
    }
#line 215 "dtrtrs.f"
    *info = 0;

/*     Solve A * x = b  or  A**T * x = b. */

#line 219 "dtrtrs.f"
    dtrsm_("Left", uplo, trans, diag, n, nrhs, &c_b12, &a[a_offset], lda, &b[
	    b_offset], ldb, (ftnlen)4, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 222 "dtrtrs.f"
    return 0;

/*     End of DTRTRS */

} /* dtrtrs_ */

