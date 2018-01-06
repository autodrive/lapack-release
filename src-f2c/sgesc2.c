#line 1 "sgesc2.f"
/* sgesc2.f -- translated by f2c (version 20100827).
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

#line 1 "sgesc2.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b SGESC2 solves a system of linear equations using the LU factorization with complete pivoting co
mputed by sgetc2. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesc2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesc2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesc2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       REAL               A( LDA, * ), RHS( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGESC2 solves a system of linear equations */
/* > */
/* >           A * X = scale* RHS */
/* > */
/* > with a general N-by-N matrix A using the LU factorization with */
/* > complete pivoting computed by SGETC2. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the  LU part of the factorization of the n-by-n */
/* >          matrix A computed by SGETC2:  A = P * L * U * Q */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHS */
/* > \verbatim */
/* >          RHS is REAL array, dimension (N). */
/* >          On entry, the right hand side vector b. */
/* >          On exit, the solution vector X. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N). */
/* >          The pivot indices; for 1 <= i <= N, row i of the */
/* >          matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] JPIV */
/* > \verbatim */
/* >          JPIV is INTEGER array, dimension (N). */
/* >          The pivot indices; for 1 <= j <= N, column j of the */
/* >          matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >           On exit, SCALE contains the scale factor. SCALE is chosen */
/* >           0 <= SCALE <= 1 to prevent owerflow in the solution. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGEauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int sgesc2_(integer *n, doublereal *a, integer *lda, 
	doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal eps, temp;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slaswp_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *);
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*      Set constant to control owerflow */

#line 156 "sgesc2.f"
    /* Parameter adjustments */
#line 156 "sgesc2.f"
    a_dim1 = *lda;
#line 156 "sgesc2.f"
    a_offset = 1 + a_dim1;
#line 156 "sgesc2.f"
    a -= a_offset;
#line 156 "sgesc2.f"
    --rhs;
#line 156 "sgesc2.f"
    --ipiv;
#line 156 "sgesc2.f"
    --jpiv;
#line 156 "sgesc2.f"

#line 156 "sgesc2.f"
    /* Function Body */
#line 156 "sgesc2.f"
    eps = slamch_("P", (ftnlen)1);
#line 157 "sgesc2.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 158 "sgesc2.f"
    bignum = 1. / smlnum;
#line 159 "sgesc2.f"
    slabad_(&smlnum, &bignum);

/*     Apply permutations IPIV to RHS */

#line 163 "sgesc2.f"
    i__1 = *n - 1;
#line 163 "sgesc2.f"
    slaswp_(&c__1, &rhs[1], lda, &c__1, &i__1, &ipiv[1], &c__1);

/*     Solve for L part */

#line 167 "sgesc2.f"
    i__1 = *n - 1;
#line 167 "sgesc2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 168 "sgesc2.f"
	i__2 = *n;
#line 168 "sgesc2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 169 "sgesc2.f"
	    rhs[j] -= a[j + i__ * a_dim1] * rhs[i__];
#line 170 "sgesc2.f"
/* L10: */
#line 170 "sgesc2.f"
	}
#line 171 "sgesc2.f"
/* L20: */
#line 171 "sgesc2.f"
    }

/*     Solve for U part */

#line 175 "sgesc2.f"
    *scale = 1.;

/*     Check for scaling */

#line 179 "sgesc2.f"
    i__ = isamax_(n, &rhs[1], &c__1);
#line 180 "sgesc2.f"
    if (smlnum * 2. * (d__1 = rhs[i__], abs(d__1)) > (d__2 = a[*n + *n * 
	    a_dim1], abs(d__2))) {
#line 181 "sgesc2.f"
	temp = .5 / (d__1 = rhs[i__], abs(d__1));
#line 182 "sgesc2.f"
	sscal_(n, &temp, &rhs[1], &c__1);
#line 183 "sgesc2.f"
	*scale *= temp;
#line 184 "sgesc2.f"
    }

#line 186 "sgesc2.f"
    for (i__ = *n; i__ >= 1; --i__) {
#line 187 "sgesc2.f"
	temp = 1. / a[i__ + i__ * a_dim1];
#line 188 "sgesc2.f"
	rhs[i__] *= temp;
#line 189 "sgesc2.f"
	i__1 = *n;
#line 189 "sgesc2.f"
	for (j = i__ + 1; j <= i__1; ++j) {
#line 190 "sgesc2.f"
	    rhs[i__] -= rhs[j] * (a[i__ + j * a_dim1] * temp);
#line 191 "sgesc2.f"
/* L30: */
#line 191 "sgesc2.f"
	}
#line 192 "sgesc2.f"
/* L40: */
#line 192 "sgesc2.f"
    }

/*     Apply permutations JPIV to the solution (RHS) */

#line 196 "sgesc2.f"
    i__1 = *n - 1;
#line 196 "sgesc2.f"
    slaswp_(&c__1, &rhs[1], lda, &c__1, &i__1, &jpiv[1], &c_n1);
#line 197 "sgesc2.f"
    return 0;

/*     End of SGESC2 */

} /* sgesc2_ */

