#line 1 "cgbtrs.f"
/* cgbtrs.f -- translated by f2c (version 20100827).
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

#line 1 "cgbtrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CGBTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbtrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbtrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbtrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AB( LDAB, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBTRS solves a system of linear equations */
/* >    A * X = B,  A**T * X = B,  or  A**H * X = B */
/* > with a general band matrix A using the LU factorization computed */
/* > by CGBTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations. */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by CGBTRF.  U is stored as an upper triangular band */
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

/* > \date December 2016 */

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, 
	doublecomplex *b, integer *ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, l, kd, lm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    cgeru_(integer *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    , cswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), ctbsv_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static logical lnoti;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
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

#line 179 "cgbtrs.f"
    /* Parameter adjustments */
#line 179 "cgbtrs.f"
    ab_dim1 = *ldab;
#line 179 "cgbtrs.f"
    ab_offset = 1 + ab_dim1;
#line 179 "cgbtrs.f"
    ab -= ab_offset;
#line 179 "cgbtrs.f"
    --ipiv;
#line 179 "cgbtrs.f"
    b_dim1 = *ldb;
#line 179 "cgbtrs.f"
    b_offset = 1 + b_dim1;
#line 179 "cgbtrs.f"
    b -= b_offset;
#line 179 "cgbtrs.f"

#line 179 "cgbtrs.f"
    /* Function Body */
#line 179 "cgbtrs.f"
    *info = 0;
#line 180 "cgbtrs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 181 "cgbtrs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 183 "cgbtrs.f"
	*info = -1;
#line 184 "cgbtrs.f"
    } else if (*n < 0) {
#line 185 "cgbtrs.f"
	*info = -2;
#line 186 "cgbtrs.f"
    } else if (*kl < 0) {
#line 187 "cgbtrs.f"
	*info = -3;
#line 188 "cgbtrs.f"
    } else if (*ku < 0) {
#line 189 "cgbtrs.f"
	*info = -4;
#line 190 "cgbtrs.f"
    } else if (*nrhs < 0) {
#line 191 "cgbtrs.f"
	*info = -5;
#line 192 "cgbtrs.f"
    } else if (*ldab < (*kl << 1) + *ku + 1) {
#line 193 "cgbtrs.f"
	*info = -7;
#line 194 "cgbtrs.f"
    } else if (*ldb < max(1,*n)) {
#line 195 "cgbtrs.f"
	*info = -10;
#line 196 "cgbtrs.f"
    }
#line 197 "cgbtrs.f"
    if (*info != 0) {
#line 198 "cgbtrs.f"
	i__1 = -(*info);
#line 198 "cgbtrs.f"
	xerbla_("CGBTRS", &i__1, (ftnlen)6);
#line 199 "cgbtrs.f"
	return 0;
#line 200 "cgbtrs.f"
    }

/*     Quick return if possible */

#line 204 "cgbtrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 204 "cgbtrs.f"
	return 0;
#line 204 "cgbtrs.f"
    }

#line 207 "cgbtrs.f"
    kd = *ku + *kl + 1;
#line 208 "cgbtrs.f"
    lnoti = *kl > 0;

#line 210 "cgbtrs.f"
    if (notran) {

/*        Solve  A*X = B. */

/*        Solve L*X = B, overwriting B with X. */

/*        L is represented as a product of permutations and unit lower */
/*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
/*        where each transformation L(i) is a rank-one modification of */
/*        the identity matrix. */

#line 221 "cgbtrs.f"
	if (lnoti) {
#line 222 "cgbtrs.f"
	    i__1 = *n - 1;
#line 222 "cgbtrs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 223 "cgbtrs.f"
		i__2 = *kl, i__3 = *n - j;
#line 223 "cgbtrs.f"
		lm = min(i__2,i__3);
#line 224 "cgbtrs.f"
		l = ipiv[j];
#line 225 "cgbtrs.f"
		if (l != j) {
#line 225 "cgbtrs.f"
		    cswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 225 "cgbtrs.f"
		}
#line 227 "cgbtrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 227 "cgbtrs.f"
		cgeru_(&lm, nrhs, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
			j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
#line 229 "cgbtrs.f"
/* L10: */
#line 229 "cgbtrs.f"
	    }
#line 230 "cgbtrs.f"
	}

#line 232 "cgbtrs.f"
	i__1 = *nrhs;
#line 232 "cgbtrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U*X = B, overwriting B with X. */

#line 236 "cgbtrs.f"
	    i__2 = *kl + *ku;
#line 236 "cgbtrs.f"
	    ctbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[
		    ab_offset], ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, 
		    (ftnlen)12, (ftnlen)8);
#line 238 "cgbtrs.f"
/* L20: */
#line 238 "cgbtrs.f"
	}

#line 240 "cgbtrs.f"
    } else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*        Solve A**T * X = B. */

#line 244 "cgbtrs.f"
	i__1 = *nrhs;
#line 244 "cgbtrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U**T * X = B, overwriting B with X. */

#line 248 "cgbtrs.f"
	    i__2 = *kl + *ku;
#line 248 "cgbtrs.f"
	    ctbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset],
		     ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9, 
		    (ftnlen)8);
#line 250 "cgbtrs.f"
/* L30: */
#line 250 "cgbtrs.f"
	}

/*        Solve L**T * X = B, overwriting B with X. */

#line 254 "cgbtrs.f"
	if (lnoti) {
#line 255 "cgbtrs.f"
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 256 "cgbtrs.f"
		i__1 = *kl, i__2 = *n - j;
#line 256 "cgbtrs.f"
		lm = min(i__1,i__2);
#line 257 "cgbtrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 257 "cgbtrs.f"
		cgemv_("Transpose", &lm, nrhs, &z__1, &b[j + 1 + b_dim1], ldb,
			 &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1, &b[j + 
			b_dim1], ldb, (ftnlen)9);
#line 259 "cgbtrs.f"
		l = ipiv[j];
#line 260 "cgbtrs.f"
		if (l != j) {
#line 260 "cgbtrs.f"
		    cswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 260 "cgbtrs.f"
		}
#line 262 "cgbtrs.f"
/* L40: */
#line 262 "cgbtrs.f"
	    }
#line 263 "cgbtrs.f"
	}

#line 265 "cgbtrs.f"
    } else {

/*        Solve A**H * X = B. */

#line 269 "cgbtrs.f"
	i__1 = *nrhs;
#line 269 "cgbtrs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Solve U**H * X = B, overwriting B with X. */

#line 273 "cgbtrs.f"
	    i__2 = *kl + *ku;
#line 273 "cgbtrs.f"
	    ctbsv_("Upper", "Conjugate transpose", "Non-unit", n, &i__2, &ab[
		    ab_offset], ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, 
		    (ftnlen)19, (ftnlen)8);
#line 275 "cgbtrs.f"
/* L50: */
#line 275 "cgbtrs.f"
	}

/*        Solve L**H * X = B, overwriting B with X. */

#line 279 "cgbtrs.f"
	if (lnoti) {
#line 280 "cgbtrs.f"
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
#line 281 "cgbtrs.f"
		i__1 = *kl, i__2 = *n - j;
#line 281 "cgbtrs.f"
		lm = min(i__1,i__2);
#line 282 "cgbtrs.f"
		clacgv_(nrhs, &b[j + b_dim1], ldb);
#line 283 "cgbtrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 283 "cgbtrs.f"
		cgemv_("Conjugate transpose", &lm, nrhs, &z__1, &b[j + 1 + 
			b_dim1], ldb, &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1,
			 &b[j + b_dim1], ldb, (ftnlen)19);
#line 286 "cgbtrs.f"
		clacgv_(nrhs, &b[j + b_dim1], ldb);
#line 287 "cgbtrs.f"
		l = ipiv[j];
#line 288 "cgbtrs.f"
		if (l != j) {
#line 288 "cgbtrs.f"
		    cswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 288 "cgbtrs.f"
		}
#line 290 "cgbtrs.f"
/* L60: */
#line 290 "cgbtrs.f"
	    }
#line 291 "cgbtrs.f"
	}
#line 292 "cgbtrs.f"
    }
#line 293 "cgbtrs.f"
    return 0;

/*     End of CGBTRS */

} /* cgbtrs_ */

