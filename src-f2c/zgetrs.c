#line 1 "zgetrs.f"
/* zgetrs.f -- translated by f2c (version 20100827).
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

#line 1 "zgetrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZGETRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGETRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGETRS solves a system of linear equations */
/* >    A * X = B,  A**T * X = B,  or  A**H * X = B */
/* > with a general N-by-N matrix A using the LU factorization computed */
/* > by ZGETRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by ZGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices from ZGETRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
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

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgetrs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical notran;
    extern /* Subroutine */ int zlaswp_(integer *, doublecomplex *, integer *,
	     integer *, integer *, integer *, integer *);


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

#line 161 "zgetrs.f"
    /* Parameter adjustments */
#line 161 "zgetrs.f"
    a_dim1 = *lda;
#line 161 "zgetrs.f"
    a_offset = 1 + a_dim1;
#line 161 "zgetrs.f"
    a -= a_offset;
#line 161 "zgetrs.f"
    --ipiv;
#line 161 "zgetrs.f"
    b_dim1 = *ldb;
#line 161 "zgetrs.f"
    b_offset = 1 + b_dim1;
#line 161 "zgetrs.f"
    b -= b_offset;
#line 161 "zgetrs.f"

#line 161 "zgetrs.f"
    /* Function Body */
#line 161 "zgetrs.f"
    *info = 0;
#line 162 "zgetrs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 163 "zgetrs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 165 "zgetrs.f"
	*info = -1;
#line 166 "zgetrs.f"
    } else if (*n < 0) {
#line 167 "zgetrs.f"
	*info = -2;
#line 168 "zgetrs.f"
    } else if (*nrhs < 0) {
#line 169 "zgetrs.f"
	*info = -3;
#line 170 "zgetrs.f"
    } else if (*lda < max(1,*n)) {
#line 171 "zgetrs.f"
	*info = -5;
#line 172 "zgetrs.f"
    } else if (*ldb < max(1,*n)) {
#line 173 "zgetrs.f"
	*info = -8;
#line 174 "zgetrs.f"
    }
#line 175 "zgetrs.f"
    if (*info != 0) {
#line 176 "zgetrs.f"
	i__1 = -(*info);
#line 176 "zgetrs.f"
	xerbla_("ZGETRS", &i__1, (ftnlen)6);
#line 177 "zgetrs.f"
	return 0;
#line 178 "zgetrs.f"
    }

/*     Quick return if possible */

#line 182 "zgetrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 182 "zgetrs.f"
	return 0;
#line 182 "zgetrs.f"
    }

#line 185 "zgetrs.f"
    if (notran) {

/*        Solve A * X = B. */

/*        Apply row interchanges to the right hand sides. */

#line 191 "zgetrs.f"
	zlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);

/*        Solve L*X = B, overwriting B with X. */

#line 195 "zgetrs.f"
	ztrsm_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b1, &a[
		a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);

/*        Solve U*X = B, overwriting B with X. */

#line 200 "zgetrs.f"
	ztrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b1, &
		a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8);
#line 202 "zgetrs.f"
    } else {

/*        Solve A**T * X = B  or A**H * X = B. */

/*        Solve U**T *X = B or U**H *X = B, overwriting B with X. */

#line 208 "zgetrs.f"
	ztrsm_("Left", "Upper", trans, "Non-unit", n, nrhs, &c_b1, &a[
		a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)1, (ftnlen)8);

/*        Solve L**T *X = B, or L**H *X = B overwriting B with X. */

#line 213 "zgetrs.f"
	ztrsm_("Left", "Lower", trans, "Unit", n, nrhs, &c_b1, &a[a_offset], 
		lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)1, (
		ftnlen)4);

/*        Apply row interchanges to the solution vectors. */

#line 218 "zgetrs.f"
	zlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
#line 219 "zgetrs.f"
    }

#line 221 "zgetrs.f"
    return 0;

/*     End of ZGETRS */

} /* zgetrs_ */

