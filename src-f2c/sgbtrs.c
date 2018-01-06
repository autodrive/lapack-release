#line 1 "sgbtrs.f"
/* sgbtrs.f -- translated by f2c (version 20100827).
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

#line 1 "sgbtrs.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b23 = 1.;

/* > \brief \b SGBTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               AB( LDAB, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBTRS solves a system of linear equations */
/* >    A * X = B  or  A**T * X = B */
/* > with a general band matrix A using the LU factorization computed */
/* > by SGBTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations. */
/* >          = 'N':  A * X = B  (No transpose) */
/* >          = 'T':  A**T* X = B  (Transpose) */
/* >          = 'C':  A**T* X = B  (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
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
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by SGBTRF.  U is stored as an upper triangular band */
/* >          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
/* >          the multipliers used during the factorization are stored in */
/* >          rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= N, row i of the matrix was */
/* >          interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv, 
	doublereal *b, integer *ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, kd, lm;
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical lnoti;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), stbsv_(char *, char *, char *, integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;


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

#line 179 "sgbtrs.f"
    /* Parameter adjustments */
#line 179 "sgbtrs.f"
    ab_dim1 = *ldab;
#line 179 "sgbtrs.f"
    ab_offset = 1 + ab_dim1;
#line 179 "sgbtrs.f"
    ab -= ab_offset;
#line 179 "sgbtrs.f"
    --ipiv;
#line 179 "sgbtrs.f"
    b_dim1 = *ldb;
#line 179 "sgbtrs.f"
    b_offset = 1 + b_dim1;
#line 179 "sgbtrs.f"
    b -= b_offset;
#line 179 "sgbtrs.f"

#line 179 "sgbtrs.f"
    /* Function Body */
#line 179 "sgbtrs.f"
    *info = 0;
#line 180 "sgbtrs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 181 "sgbtrs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 183 "sgbtrs.f"
	*info = -1;
#line 184 "sgbtrs.f"
    } else if (*n < 0) {
#line 185 "sgbtrs.f"
	*info = -2;
#line 186 "sgbtrs.f"
    } else if (*kl < 0) {
#line 187 "sgbtrs.f"
	*info = -3;
#line 188 "sgbtrs.f"
    } else if (*ku < 0) {
#line 189 "sgbtrs.f"
	*info = -4;
#line 190 "sgbtrs.f"
    } else if (*nrhs < 0) {
#line 191 "sgbtrs.f"
	*info = -5;
#line 192 "sgbtrs.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 193 "sgbtrs.f"
	*info = -7;
#line 194 "sgbtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 195 "sgbtrs.f"
	*info = -10;
#line 196 "sgbtrs.f"
    }
#line 197 "sgbtrs.f"
    if (*info != 0) {
#line 198 "sgbtrs.f"
	i__1 = -(*info);
#line 198 "sgbtrs.f"
	xerbla_("SGBTRS", &i__1, (ftnlen)6);
#line 199 "sgbtrs.f"
	return 0;
#line 200 "sgbtrs.f"
    }

/*     Quick return if possible */

#line 204 "sgbtrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 204 "sgbtrs.f"
	return 0;
#line 204 "sgbtrs.f"
    }

#line 207 "sgbtrs.f"
    kd = *ku + *kl + 1;
#line 208 "sgbtrs.f"
    lnoti = *kl > 0;

#line 210 "sgbtrs.f"
    if (notran) {

/*        Solve  A*X = B. */

/*        Solve L*X = B, overwriting B with X. */

/*        L is represented as a product of permutations and unit lower */
/*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
/*        where each transformation L(i) is a rank-one modification of */
/*        the identity matrix. */

#line 221 "sgbtrs.f"
	if (lnoti) {
#line 222 "sgbtrs.f"
	    i__1 = *n - 1;
#line 222 "sgbtrs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 223 "sgbtrs.f"
		i__2 = *kl, i__3 = *n - j;
#line 223 "sgbtrs.f"
		lm = min(i__2,i__3);
#line 224 "sgbtrs.f"
		l = ipiv[j];
#line 225 "sgbtrs.f"
		if (l != j) {
#line 225 "sgbtrs.f"
		    sswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 225 "sgbtrs.f"
		}
#line 227 "sgbtrs.f"
		sger_(&lm, nrhs, &c_b7, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
			j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
#line 229 "sgbtrs.f"
/* L10: */
#line 229 "sgbtrs.f"
	    }
#line 230 "sgbtrs.f"
	}

#line 232 "sgbtrs.f"
	i__1 = *nrhs;
#line 232 "sgbtrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U*X = B, overwriting B with X. */

#line 236 "sgbtrs.f"
	    i__2 = *kl + *ku;
#line 236 "sgbtrs.f"
	    stbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[
		    ab_offset], ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, 
		    (ftnlen)12, (ftnlen)8);
#line 238 "sgbtrs.f"
/* L20: */
#line 238 "sgbtrs.f"
	}

#line 240 "sgbtrs.f"
    } else {

/*        Solve A**T*X = B. */

#line 244 "sgbtrs.f"
	i__1 = *nrhs;
#line 244 "sgbtrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U**T*X = B, overwriting B with X. */

#line 248 "sgbtrs.f"
	    i__2 = *kl + *ku;
#line 248 "sgbtrs.f"
	    stbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset],
		     ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9, 
		    (ftnlen)8);
#line 250 "sgbtrs.f"
/* L30: */
#line 250 "sgbtrs.f"
	}

/*        Solve L**T*X = B, overwriting B with X. */

#line 254 "sgbtrs.f"
	if (lnoti) {
#line 255 "sgbtrs.f"
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 256 "sgbtrs.f"
		i__1 = *kl, i__2 = *n - j;
#line 256 "sgbtrs.f"
		lm = min(i__1,i__2);
#line 257 "sgbtrs.f"
		sgemv_("Transpose", &lm, nrhs, &c_b7, &b[j + 1 + b_dim1], ldb,
			 &ab[kd + 1 + j * ab_dim1], &c__1, &c_b23, &b[j + 
			b_dim1], ldb, (ftnlen)9);
#line 259 "sgbtrs.f"
		l = ipiv[j];
#line 260 "sgbtrs.f"
		if (l != j) {
#line 260 "sgbtrs.f"
		    sswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 260 "sgbtrs.f"
		}
#line 262 "sgbtrs.f"
/* L40: */
#line 262 "sgbtrs.f"
	    }
#line 263 "sgbtrs.f"
	}
#line 264 "sgbtrs.f"
    }
#line 265 "sgbtrs.f"
    return 0;

/*     End of SGBTRS */

} /* sgbtrs_ */

