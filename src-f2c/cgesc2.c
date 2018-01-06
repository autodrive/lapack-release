#line 1 "cgesc2.f"
/* cgesc2.f -- translated by f2c (version 20100827).
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

#line 1 "cgesc2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b13 = {1.,0.};
static integer c_n1 = -1;

/* > \brief \b CGESC2 solves a system of linear equations using the LU factorization with complete pivoting co
mputed by sgetc2. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesc2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesc2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesc2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       COMPLEX            A( LDA, * ), RHS( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESC2 solves a system of linear equations */
/* > */
/* >           A * X = scale* RHS */
/* > */
/* > with a general N-by-N matrix A using the LU factorization with */
/* > complete pivoting computed by CGETC2. */
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
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the  LU part of the factorization of the n-by-n */
/* >          matrix A computed by CGETC2:  A = P * L * U * Q */
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
/* >          RHS is COMPLEX array, dimension N. */
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

/* > \date September 2012 */

/* > \ingroup complexGEauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int cgesc2_(integer *n, doublecomplex *a, integer *lda, 
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), slabad_(doublereal *, doublereal *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int claswp_(integer *, doublecomplex *, integer *,
	     integer *, integer *, integer *, integer *);
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 158 "cgesc2.f"
    /* Parameter adjustments */
#line 158 "cgesc2.f"
    a_dim1 = *lda;
#line 158 "cgesc2.f"
    a_offset = 1 + a_dim1;
#line 158 "cgesc2.f"
    a -= a_offset;
#line 158 "cgesc2.f"
    --rhs;
#line 158 "cgesc2.f"
    --ipiv;
#line 158 "cgesc2.f"
    --jpiv;
#line 158 "cgesc2.f"

#line 158 "cgesc2.f"
    /* Function Body */
#line 158 "cgesc2.f"
    eps = slamch_("P", (ftnlen)1);
#line 159 "cgesc2.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 160 "cgesc2.f"
    bignum = 1. / smlnum;
#line 161 "cgesc2.f"
    slabad_(&smlnum, &bignum);

/*     Apply permutations IPIV to RHS */

#line 165 "cgesc2.f"
    i__1 = *n - 1;
#line 165 "cgesc2.f"
    claswp_(&c__1, &rhs[1], lda, &c__1, &i__1, &ipiv[1], &c__1);

/*     Solve for L part */

#line 169 "cgesc2.f"
    i__1 = *n - 1;
#line 169 "cgesc2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 170 "cgesc2.f"
	i__2 = *n;
#line 170 "cgesc2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 171 "cgesc2.f"
	    i__3 = j;
#line 171 "cgesc2.f"
	    i__4 = j;
#line 171 "cgesc2.f"
	    i__5 = j + i__ * a_dim1;
#line 171 "cgesc2.f"
	    i__6 = i__;
#line 171 "cgesc2.f"
	    z__2.r = a[i__5].r * rhs[i__6].r - a[i__5].i * rhs[i__6].i, 
		    z__2.i = a[i__5].r * rhs[i__6].i + a[i__5].i * rhs[i__6]
		    .r;
#line 171 "cgesc2.f"
	    z__1.r = rhs[i__4].r - z__2.r, z__1.i = rhs[i__4].i - z__2.i;
#line 171 "cgesc2.f"
	    rhs[i__3].r = z__1.r, rhs[i__3].i = z__1.i;
#line 172 "cgesc2.f"
/* L10: */
#line 172 "cgesc2.f"
	}
#line 173 "cgesc2.f"
/* L20: */
#line 173 "cgesc2.f"
    }

/*     Solve for U part */

#line 177 "cgesc2.f"
    *scale = 1.;

/*     Check for scaling */

#line 181 "cgesc2.f"
    i__ = icamax_(n, &rhs[1], &c__1);
#line 182 "cgesc2.f"
    if (smlnum * 2. * z_abs(&rhs[i__]) > z_abs(&a[*n + *n * a_dim1])) {
#line 183 "cgesc2.f"
	d__1 = z_abs(&rhs[i__]);
#line 183 "cgesc2.f"
	z__1.r = .5 / d__1, z__1.i = 0. / d__1;
#line 183 "cgesc2.f"
	temp.r = z__1.r, temp.i = z__1.i;
#line 184 "cgesc2.f"
	cscal_(n, &temp, &rhs[1], &c__1);
#line 185 "cgesc2.f"
	*scale *= temp.r;
#line 186 "cgesc2.f"
    }
#line 187 "cgesc2.f"
    for (i__ = *n; i__ >= 1; --i__) {
#line 188 "cgesc2.f"
	z_div(&z__1, &c_b13, &a[i__ + i__ * a_dim1]);
#line 188 "cgesc2.f"
	temp.r = z__1.r, temp.i = z__1.i;
#line 189 "cgesc2.f"
	i__1 = i__;
#line 189 "cgesc2.f"
	i__2 = i__;
#line 189 "cgesc2.f"
	z__1.r = rhs[i__2].r * temp.r - rhs[i__2].i * temp.i, z__1.i = rhs[
		i__2].r * temp.i + rhs[i__2].i * temp.r;
#line 189 "cgesc2.f"
	rhs[i__1].r = z__1.r, rhs[i__1].i = z__1.i;
#line 190 "cgesc2.f"
	i__1 = *n;
#line 190 "cgesc2.f"
	for (j = i__ + 1; j <= i__1; ++j) {
#line 191 "cgesc2.f"
	    i__2 = i__;
#line 191 "cgesc2.f"
	    i__3 = i__;
#line 191 "cgesc2.f"
	    i__4 = j;
#line 191 "cgesc2.f"
	    i__5 = i__ + j * a_dim1;
#line 191 "cgesc2.f"
	    z__3.r = a[i__5].r * temp.r - a[i__5].i * temp.i, z__3.i = a[i__5]
		    .r * temp.i + a[i__5].i * temp.r;
#line 191 "cgesc2.f"
	    z__2.r = rhs[i__4].r * z__3.r - rhs[i__4].i * z__3.i, z__2.i = 
		    rhs[i__4].r * z__3.i + rhs[i__4].i * z__3.r;
#line 191 "cgesc2.f"
	    z__1.r = rhs[i__3].r - z__2.r, z__1.i = rhs[i__3].i - z__2.i;
#line 191 "cgesc2.f"
	    rhs[i__2].r = z__1.r, rhs[i__2].i = z__1.i;
#line 192 "cgesc2.f"
/* L30: */
#line 192 "cgesc2.f"
	}
#line 193 "cgesc2.f"
/* L40: */
#line 193 "cgesc2.f"
    }

/*     Apply permutations JPIV to the solution (RHS) */

#line 197 "cgesc2.f"
    i__1 = *n - 1;
#line 197 "cgesc2.f"
    claswp_(&c__1, &rhs[1], lda, &c__1, &i__1, &jpiv[1], &c_n1);
#line 198 "cgesc2.f"
    return 0;

/*     End of CGESC2 */

} /* cgesc2_ */

