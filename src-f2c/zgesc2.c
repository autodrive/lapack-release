#line 1 "zgesc2.f"
/* zgesc2.f -- translated by f2c (version 20100827).
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

#line 1 "zgesc2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b13 = {1.,0.};
static integer c_n1 = -1;

/* > \brief \b ZGESC2 solves a system of linear equations using the LU factorization with complete pivoting co
mputed by sgetc2. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGESC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesc2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesc2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesc2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), RHS( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGESC2 solves a system of linear equations */
/* > */
/* >           A * X = scale* RHS */
/* > */
/* > with a general N-by-N matrix A using the LU factorization with */
/* > complete pivoting computed by ZGETC2. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the  LU part of the factorization of the n-by-n */
/* >          matrix A computed by ZGETC2:  A = P * L * U * Q */
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
/* >          RHS is COMPLEX*16 array, dimension N. */
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
/* >          SCALE is DOUBLE PRECISION */
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

/* > \ingroup complex16GEauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int zgesc2_(integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *rhs, integer *ipiv, integer *jpiv, doublereal *scale)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal eps;
    static doublecomplex temp;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal bignum;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static doublereal smlnum;
    extern /* Subroutine */ int zlaswp_(integer *, doublecomplex *, integer *,
	     integer *, integer *, integer *, integer *);


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

/*     Set constant to control overflow */

#line 158 "zgesc2.f"
    /* Parameter adjustments */
#line 158 "zgesc2.f"
    a_dim1 = *lda;
#line 158 "zgesc2.f"
    a_offset = 1 + a_dim1;
#line 158 "zgesc2.f"
    a -= a_offset;
#line 158 "zgesc2.f"
    --rhs;
#line 158 "zgesc2.f"
    --ipiv;
#line 158 "zgesc2.f"
    --jpiv;
#line 158 "zgesc2.f"

#line 158 "zgesc2.f"
    /* Function Body */
#line 158 "zgesc2.f"
    eps = dlamch_("P", (ftnlen)1);
#line 159 "zgesc2.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 160 "zgesc2.f"
    bignum = 1. / smlnum;
#line 161 "zgesc2.f"
    dlabad_(&smlnum, &bignum);

/*     Apply permutations IPIV to RHS */

#line 165 "zgesc2.f"
    i__1 = *n - 1;
#line 165 "zgesc2.f"
    zlaswp_(&c__1, &rhs[1], lda, &c__1, &i__1, &ipiv[1], &c__1);

/*     Solve for L part */

#line 169 "zgesc2.f"
    i__1 = *n - 1;
#line 169 "zgesc2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 170 "zgesc2.f"
	i__2 = *n;
#line 170 "zgesc2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 171 "zgesc2.f"
	    i__3 = j;
#line 171 "zgesc2.f"
	    i__4 = j;
#line 171 "zgesc2.f"
	    i__5 = j + i__ * a_dim1;
#line 171 "zgesc2.f"
	    i__6 = i__;
#line 171 "zgesc2.f"
	    z__2.r = a[i__5].r * rhs[i__6].r - a[i__5].i * rhs[i__6].i, 
		    z__2.i = a[i__5].r * rhs[i__6].i + a[i__5].i * rhs[i__6]
		    .r;
#line 171 "zgesc2.f"
	    z__1.r = rhs[i__4].r - z__2.r, z__1.i = rhs[i__4].i - z__2.i;
#line 171 "zgesc2.f"
	    rhs[i__3].r = z__1.r, rhs[i__3].i = z__1.i;
#line 172 "zgesc2.f"
/* L10: */
#line 172 "zgesc2.f"
	}
#line 173 "zgesc2.f"
/* L20: */
#line 173 "zgesc2.f"
    }

/*     Solve for U part */

#line 177 "zgesc2.f"
    *scale = 1.;

/*     Check for scaling */

#line 181 "zgesc2.f"
    i__ = izamax_(n, &rhs[1], &c__1);
#line 182 "zgesc2.f"
    if (smlnum * 2. * z_abs(&rhs[i__]) > z_abs(&a[*n + *n * a_dim1])) {
#line 183 "zgesc2.f"
	d__1 = z_abs(&rhs[i__]);
#line 183 "zgesc2.f"
	z__1.r = .5 / d__1, z__1.i = 0. / d__1;
#line 183 "zgesc2.f"
	temp.r = z__1.r, temp.i = z__1.i;
#line 184 "zgesc2.f"
	zscal_(n, &temp, &rhs[1], &c__1);
#line 185 "zgesc2.f"
	*scale *= temp.r;
#line 186 "zgesc2.f"
    }
#line 187 "zgesc2.f"
    for (i__ = *n; i__ >= 1; --i__) {
#line 188 "zgesc2.f"
	z_div(&z__1, &c_b13, &a[i__ + i__ * a_dim1]);
#line 188 "zgesc2.f"
	temp.r = z__1.r, temp.i = z__1.i;
#line 189 "zgesc2.f"
	i__1 = i__;
#line 189 "zgesc2.f"
	i__2 = i__;
#line 189 "zgesc2.f"
	z__1.r = rhs[i__2].r * temp.r - rhs[i__2].i * temp.i, z__1.i = rhs[
		i__2].r * temp.i + rhs[i__2].i * temp.r;
#line 189 "zgesc2.f"
	rhs[i__1].r = z__1.r, rhs[i__1].i = z__1.i;
#line 190 "zgesc2.f"
	i__1 = *n;
#line 190 "zgesc2.f"
	for (j = i__ + 1; j <= i__1; ++j) {
#line 191 "zgesc2.f"
	    i__2 = i__;
#line 191 "zgesc2.f"
	    i__3 = i__;
#line 191 "zgesc2.f"
	    i__4 = j;
#line 191 "zgesc2.f"
	    i__5 = i__ + j * a_dim1;
#line 191 "zgesc2.f"
	    z__3.r = a[i__5].r * temp.r - a[i__5].i * temp.i, z__3.i = a[i__5]
		    .r * temp.i + a[i__5].i * temp.r;
#line 191 "zgesc2.f"
	    z__2.r = rhs[i__4].r * z__3.r - rhs[i__4].i * z__3.i, z__2.i = 
		    rhs[i__4].r * z__3.i + rhs[i__4].i * z__3.r;
#line 191 "zgesc2.f"
	    z__1.r = rhs[i__3].r - z__2.r, z__1.i = rhs[i__3].i - z__2.i;
#line 191 "zgesc2.f"
	    rhs[i__2].r = z__1.r, rhs[i__2].i = z__1.i;
#line 192 "zgesc2.f"
/* L30: */
#line 192 "zgesc2.f"
	}
#line 193 "zgesc2.f"
/* L40: */
#line 193 "zgesc2.f"
    }

/*     Apply permutations JPIV to the solution (RHS) */

#line 197 "zgesc2.f"
    i__1 = *n - 1;
#line 197 "zgesc2.f"
    zlaswp_(&c__1, &rhs[1], lda, &c__1, &i__1, &jpiv[1], &c_n1);
#line 198 "zgesc2.f"
    return 0;

/*     End of ZGESC2 */

} /* zgesc2_ */

