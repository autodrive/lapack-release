#line 1 "zpotrs.f"
/* zpotrs.f -- translated by f2c (version 20100827).
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

#line 1 "zpotrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZPOTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPOTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpotrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpotrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpotrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
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
/* > ZPOTRS solves a system of linear equations A*X = B with a Hermitian */
/* > positive definite matrix A using the Cholesky factorization */
/* > A = U**H * U or A = L * L**H computed by ZPOTRF. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H * U or A = L * L**H, as computed by ZPOTRF. */
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

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpotrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

#line 149 "zpotrs.f"
    /* Parameter adjustments */
#line 149 "zpotrs.f"
    a_dim1 = *lda;
#line 149 "zpotrs.f"
    a_offset = 1 + a_dim1;
#line 149 "zpotrs.f"
    a -= a_offset;
#line 149 "zpotrs.f"
    b_dim1 = *ldb;
#line 149 "zpotrs.f"
    b_offset = 1 + b_dim1;
#line 149 "zpotrs.f"
    b -= b_offset;
#line 149 "zpotrs.f"

#line 149 "zpotrs.f"
    /* Function Body */
#line 149 "zpotrs.f"
    *info = 0;
#line 150 "zpotrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 151 "zpotrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 152 "zpotrs.f"
	*info = -1;
#line 153 "zpotrs.f"
    } else if (*n < 0) {
#line 154 "zpotrs.f"
	*info = -2;
#line 155 "zpotrs.f"
    } else if (*nrhs < 0) {
#line 156 "zpotrs.f"
	*info = -3;
#line 157 "zpotrs.f"
    } else if (*lda < max(1,*n)) {
#line 158 "zpotrs.f"
	*info = -5;
#line 159 "zpotrs.f"
    } else if (*ldb < max(1,*n)) {
#line 160 "zpotrs.f"
	*info = -7;
#line 161 "zpotrs.f"
    }
#line 162 "zpotrs.f"
    if (*info != 0) {
#line 163 "zpotrs.f"
	i__1 = -(*info);
#line 163 "zpotrs.f"
	xerbla_("ZPOTRS", &i__1, (ftnlen)6);
#line 164 "zpotrs.f"
	return 0;
#line 165 "zpotrs.f"
    }

/*     Quick return if possible */

#line 169 "zpotrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 169 "zpotrs.f"
	return 0;
#line 169 "zpotrs.f"
    }

#line 172 "zpotrs.f"
    if (upper) {

/*        Solve A*X = B where A = U**H *U. */

/*        Solve U**H *X = B, overwriting B with X. */

#line 178 "zpotrs.f"
	ztrsm_("Left", "Upper", "Conjugate transpose", "Non-unit", n, nrhs, &
		c_b1, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (
		ftnlen)5, (ftnlen)19, (ftnlen)8);

/*        Solve U*X = B, overwriting B with X. */

#line 183 "zpotrs.f"
	ztrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b1, &
		a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8);
#line 185 "zpotrs.f"
    } else {

/*        Solve A*X = B where A = L*L**H. */

/*        Solve L*X = B, overwriting B with X. */

#line 191 "zpotrs.f"
	ztrsm_("Left", "Lower", "No transpose", "Non-unit", n, nrhs, &c_b1, &
		a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8);

/*        Solve L**H *X = B, overwriting B with X. */

#line 196 "zpotrs.f"
	ztrsm_("Left", "Lower", "Conjugate transpose", "Non-unit", n, nrhs, &
		c_b1, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (
		ftnlen)5, (ftnlen)19, (ftnlen)8);
#line 198 "zpotrs.f"
    }

#line 200 "zpotrs.f"
    return 0;

/*     End of ZPOTRS */

} /* zpotrs_ */

