#line 1 "ztrtrs.f"
/* ztrtrs.f -- translated by f2c (version 20100827).
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

#line 1 "ztrtrs.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};

/* > \brief \b ZTRTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTRTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRTRS solves a triangular system of the form */
/* > */
/* >    A * X = B,  A**T * X = B,  or  A**H * X = B, */
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
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *info, ftnlen uplo_len, ftnlen trans_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical nounit;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 180 "ztrtrs.f"
    /* Parameter adjustments */
#line 180 "ztrtrs.f"
    a_dim1 = *lda;
#line 180 "ztrtrs.f"
    a_offset = 1 + a_dim1;
#line 180 "ztrtrs.f"
    a -= a_offset;
#line 180 "ztrtrs.f"
    b_dim1 = *ldb;
#line 180 "ztrtrs.f"
    b_offset = 1 + b_dim1;
#line 180 "ztrtrs.f"
    b -= b_offset;
#line 180 "ztrtrs.f"

#line 180 "ztrtrs.f"
    /* Function Body */
#line 180 "ztrtrs.f"
    *info = 0;
#line 181 "ztrtrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 182 "ztrtrs.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 183 "ztrtrs.f"
	*info = -1;
#line 184 "ztrtrs.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 186 "ztrtrs.f"
	*info = -2;
#line 187 "ztrtrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 188 "ztrtrs.f"
	*info = -3;
#line 189 "ztrtrs.f"
    } else if (*n < 0) {
#line 190 "ztrtrs.f"
	*info = -4;
#line 191 "ztrtrs.f"
    } else if (*nrhs < 0) {
#line 192 "ztrtrs.f"
	*info = -5;
#line 193 "ztrtrs.f"
    } else if (*lda < max(1,*n)) {
#line 194 "ztrtrs.f"
	*info = -7;
#line 195 "ztrtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 196 "ztrtrs.f"
	*info = -9;
#line 197 "ztrtrs.f"
    }
#line 198 "ztrtrs.f"
    if (*info != 0) {
#line 199 "ztrtrs.f"
	i__1 = -(*info);
#line 199 "ztrtrs.f"
	xerbla_("ZTRTRS", &i__1, (ftnlen)6);
#line 200 "ztrtrs.f"
	return 0;
#line 201 "ztrtrs.f"
    }

/*     Quick return if possible */

#line 205 "ztrtrs.f"
    if (*n == 0) {
#line 205 "ztrtrs.f"
	return 0;
#line 205 "ztrtrs.f"
    }

/*     Check for singularity. */

#line 210 "ztrtrs.f"
    if (nounit) {
#line 211 "ztrtrs.f"
	i__1 = *n;
#line 211 "ztrtrs.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 212 "ztrtrs.f"
	    i__2 = *info + *info * a_dim1;
#line 212 "ztrtrs.f"
	    if (a[i__2].r == 0. && a[i__2].i == 0.) {
#line 212 "ztrtrs.f"
		return 0;
#line 212 "ztrtrs.f"
	    }
#line 214 "ztrtrs.f"
/* L10: */
#line 214 "ztrtrs.f"
	}
#line 215 "ztrtrs.f"
    }
#line 216 "ztrtrs.f"
    *info = 0;

/*     Solve A * x = b,  A**T * x = b,  or  A**H * x = b. */

#line 220 "ztrtrs.f"
    ztrsm_("Left", uplo, trans, diag, n, nrhs, &c_b2, &a[a_offset], lda, &b[
	    b_offset], ldb, (ftnlen)4, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 223 "ztrtrs.f"
    return 0;

/*     End of ZTRTRS */

} /* ztrtrs_ */

