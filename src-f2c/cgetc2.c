#line 1 "cgetc2.f"
/* cgetc2.f -- translated by f2c (version 20100827).
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

#line 1 "cgetc2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b10 = {-1.,-0.};

/* > \brief \b CGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGETC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetc2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetc2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetc2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGETC2( N, A, LDA, IPIV, JPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETC2 computes an LU factorization, using complete pivoting, of the */
/* > n-by-n matrix A. The factorization has the form A = P * L * U * Q, */
/* > where P and Q are permutation matrices, L is lower triangular with */
/* > unit diagonal elements and U is upper triangular. */
/* > */
/* > This is a level 1 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the n-by-n matrix to be factored. */
/* >          On exit, the factors L and U from the factorization */
/* >          A = P*L*U*Q; the unit diagonal elements of L are not stored. */
/* >          If U(k, k) appears to be less than SMIN, U(k, k) is given the */
/* >          value of SMIN, giving a nonsingular perturbed system. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N). */
/* >          The pivot indices; for 1 <= i <= N, row i of the */
/* >          matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] JPIV */
/* > \verbatim */
/* >          JPIV is INTEGER array, dimension (N). */
/* >          The pivot indices; for 1 <= j <= N, column j of the */
/* >          matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0: successful exit */
/* >           > 0: if INFO = k, U(k, k) is likely to produce overflow if */
/* >                one tries to solve for x in Ax = b. So U is perturbed */
/* >                to avoid the overflow. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup complexGEauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int cgetc2_(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, integer *jpiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ip, jp;
    static doublereal eps;
    static integer ipv, jpv;
    static doublereal smin, xmax;
    extern /* Subroutine */ int cgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), slabad_(doublereal *, 
	    doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

/*     Set constants to control overflow */

#line 151 "cgetc2.f"
    /* Parameter adjustments */
#line 151 "cgetc2.f"
    a_dim1 = *lda;
#line 151 "cgetc2.f"
    a_offset = 1 + a_dim1;
#line 151 "cgetc2.f"
    a -= a_offset;
#line 151 "cgetc2.f"
    --ipiv;
#line 151 "cgetc2.f"
    --jpiv;
#line 151 "cgetc2.f"

#line 151 "cgetc2.f"
    /* Function Body */
#line 151 "cgetc2.f"
    *info = 0;
#line 152 "cgetc2.f"
    eps = slamch_("P", (ftnlen)1);
#line 153 "cgetc2.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 154 "cgetc2.f"
    bignum = 1. / smlnum;
#line 155 "cgetc2.f"
    slabad_(&smlnum, &bignum);

/*     Factorize A using complete pivoting. */
/*     Set pivots less than SMIN to SMIN */

#line 160 "cgetc2.f"
    i__1 = *n - 1;
#line 160 "cgetc2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Find max element in matrix A */

#line 164 "cgetc2.f"
	xmax = 0.;
#line 165 "cgetc2.f"
	i__2 = *n;
#line 165 "cgetc2.f"
	for (ip = i__; ip <= i__2; ++ip) {
#line 166 "cgetc2.f"
	    i__3 = *n;
#line 166 "cgetc2.f"
	    for (jp = i__; jp <= i__3; ++jp) {
#line 167 "cgetc2.f"
		if (z_abs(&a[ip + jp * a_dim1]) >= xmax) {
#line 168 "cgetc2.f"
		    xmax = z_abs(&a[ip + jp * a_dim1]);
#line 169 "cgetc2.f"
		    ipv = ip;
#line 170 "cgetc2.f"
		    jpv = jp;
#line 171 "cgetc2.f"
		}
#line 172 "cgetc2.f"
/* L10: */
#line 172 "cgetc2.f"
	    }
#line 173 "cgetc2.f"
/* L20: */
#line 173 "cgetc2.f"
	}
#line 174 "cgetc2.f"
	if (i__ == 1) {
/* Computing MAX */
#line 174 "cgetc2.f"
	    d__1 = eps * xmax;
#line 174 "cgetc2.f"
	    smin = max(d__1,smlnum);
#line 174 "cgetc2.f"
	}

/*        Swap rows */

#line 179 "cgetc2.f"
	if (ipv != i__) {
#line 179 "cgetc2.f"
	    cswap_(n, &a[ipv + a_dim1], lda, &a[i__ + a_dim1], lda);
#line 179 "cgetc2.f"
	}
#line 181 "cgetc2.f"
	ipiv[i__] = ipv;

/*        Swap columns */

#line 185 "cgetc2.f"
	if (jpv != i__) {
#line 185 "cgetc2.f"
	    cswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
#line 185 "cgetc2.f"
	}
#line 187 "cgetc2.f"
	jpiv[i__] = jpv;

/*        Check for singularity */

#line 191 "cgetc2.f"
	if (z_abs(&a[i__ + i__ * a_dim1]) < smin) {
#line 192 "cgetc2.f"
	    *info = i__;
#line 193 "cgetc2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 193 "cgetc2.f"
	    z__1.r = smin, z__1.i = 0.;
#line 193 "cgetc2.f"
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 194 "cgetc2.f"
	}
#line 195 "cgetc2.f"
	i__2 = *n;
#line 195 "cgetc2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 196 "cgetc2.f"
	    i__3 = j + i__ * a_dim1;
#line 196 "cgetc2.f"
	    z_div(&z__1, &a[j + i__ * a_dim1], &a[i__ + i__ * a_dim1]);
#line 196 "cgetc2.f"
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 197 "cgetc2.f"
/* L30: */
#line 197 "cgetc2.f"
	}
#line 198 "cgetc2.f"
	i__2 = *n - i__;
#line 198 "cgetc2.f"
	i__3 = *n - i__;
#line 198 "cgetc2.f"
	cgeru_(&i__2, &i__3, &c_b10, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[
		i__ + (i__ + 1) * a_dim1], lda, &a[i__ + 1 + (i__ + 1) * 
		a_dim1], lda);
#line 200 "cgetc2.f"
/* L40: */
#line 200 "cgetc2.f"
    }

#line 202 "cgetc2.f"
    if (z_abs(&a[*n + *n * a_dim1]) < smin) {
#line 203 "cgetc2.f"
	*info = *n;
#line 204 "cgetc2.f"
	i__1 = *n + *n * a_dim1;
#line 204 "cgetc2.f"
	z__1.r = smin, z__1.i = 0.;
#line 204 "cgetc2.f"
	a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 205 "cgetc2.f"
    }

/*     Set last pivots to N */

#line 209 "cgetc2.f"
    ipiv[*n] = *n;
#line 210 "cgetc2.f"
    jpiv[*n] = *n;

#line 212 "cgetc2.f"
    return 0;

/*     End of CGETC2 */

} /* cgetc2_ */

