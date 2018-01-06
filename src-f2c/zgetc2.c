#line 1 "zgetc2.f"
/* zgetc2.f -- translated by f2c (version 20100827).
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

#line 1 "zgetc2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b10 = {-1.,-0.};

/* > \brief \b ZGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGETC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetc2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetc2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetc2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGETC2( N, A, LDA, IPIV, JPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), JPIV( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGETC2 computes an LU factorization, using complete pivoting, of the */
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
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
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

/* > \date June 2016 */

/* > \ingroup complex16GEauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int zgetc2_(integer *n, doublecomplex *a, integer *lda, 
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
    extern /* Subroutine */ int zgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), dlabad_(doublereal *, 
	    doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 149 "zgetc2.f"
    /* Parameter adjustments */
#line 149 "zgetc2.f"
    a_dim1 = *lda;
#line 149 "zgetc2.f"
    a_offset = 1 + a_dim1;
#line 149 "zgetc2.f"
    a -= a_offset;
#line 149 "zgetc2.f"
    --ipiv;
#line 149 "zgetc2.f"
    --jpiv;
#line 149 "zgetc2.f"

#line 149 "zgetc2.f"
    /* Function Body */
#line 149 "zgetc2.f"
    *info = 0;

/*     Quick return if possible */

#line 153 "zgetc2.f"
    if (*n == 0) {
#line 153 "zgetc2.f"
	return 0;
#line 153 "zgetc2.f"
    }

/*     Set constants to control overflow */

#line 158 "zgetc2.f"
    eps = dlamch_("P", (ftnlen)1);
#line 159 "zgetc2.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 160 "zgetc2.f"
    bignum = 1. / smlnum;
#line 161 "zgetc2.f"
    dlabad_(&smlnum, &bignum);

/*     Handle the case N=1 by itself */

#line 165 "zgetc2.f"
    if (*n == 1) {
#line 166 "zgetc2.f"
	ipiv[1] = 1;
#line 167 "zgetc2.f"
	jpiv[1] = 1;
#line 168 "zgetc2.f"
	if (z_abs(&a[a_dim1 + 1]) < smlnum) {
#line 169 "zgetc2.f"
	    *info = 1;
#line 170 "zgetc2.f"
	    i__1 = a_dim1 + 1;
#line 170 "zgetc2.f"
	    z__1.r = smlnum, z__1.i = 0.;
#line 170 "zgetc2.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 171 "zgetc2.f"
	}
#line 172 "zgetc2.f"
	return 0;
#line 173 "zgetc2.f"
    }

/*     Factorize A using complete pivoting. */
/*     Set pivots less than SMIN to SMIN */

#line 178 "zgetc2.f"
    i__1 = *n - 1;
#line 178 "zgetc2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Find max element in matrix A */

#line 182 "zgetc2.f"
	xmax = 0.;
#line 183 "zgetc2.f"
	i__2 = *n;
#line 183 "zgetc2.f"
	for (ip = i__; ip <= i__2; ++ip) {
#line 184 "zgetc2.f"
	    i__3 = *n;
#line 184 "zgetc2.f"
	    for (jp = i__; jp <= i__3; ++jp) {
#line 185 "zgetc2.f"
		if (z_abs(&a[ip + jp * a_dim1]) >= xmax) {
#line 186 "zgetc2.f"
		    xmax = z_abs(&a[ip + jp * a_dim1]);
#line 187 "zgetc2.f"
		    ipv = ip;
#line 188 "zgetc2.f"
		    jpv = jp;
#line 189 "zgetc2.f"
		}
#line 190 "zgetc2.f"
/* L10: */
#line 190 "zgetc2.f"
	    }
#line 191 "zgetc2.f"
/* L20: */
#line 191 "zgetc2.f"
	}
#line 192 "zgetc2.f"
	if (i__ == 1) {
/* Computing MAX */
#line 192 "zgetc2.f"
	    d__1 = eps * xmax;
#line 192 "zgetc2.f"
	    smin = max(d__1,smlnum);
#line 192 "zgetc2.f"
	}

/*        Swap rows */

#line 197 "zgetc2.f"
	if (ipv != i__) {
#line 197 "zgetc2.f"
	    zswap_(n, &a[ipv + a_dim1], lda, &a[i__ + a_dim1], lda);
#line 197 "zgetc2.f"
	}
#line 199 "zgetc2.f"
	ipiv[i__] = ipv;

/*        Swap columns */

#line 203 "zgetc2.f"
	if (jpv != i__) {
#line 203 "zgetc2.f"
	    zswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
#line 203 "zgetc2.f"
	}
#line 205 "zgetc2.f"
	jpiv[i__] = jpv;

/*        Check for singularity */

#line 209 "zgetc2.f"
	if (z_abs(&a[i__ + i__ * a_dim1]) < smin) {
#line 210 "zgetc2.f"
	    *info = i__;
#line 211 "zgetc2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 211 "zgetc2.f"
	    z__1.r = smin, z__1.i = 0.;
#line 211 "zgetc2.f"
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 212 "zgetc2.f"
	}
#line 213 "zgetc2.f"
	i__2 = *n;
#line 213 "zgetc2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 214 "zgetc2.f"
	    i__3 = j + i__ * a_dim1;
#line 214 "zgetc2.f"
	    z_div(&z__1, &a[j + i__ * a_dim1], &a[i__ + i__ * a_dim1]);
#line 214 "zgetc2.f"
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 215 "zgetc2.f"
/* L30: */
#line 215 "zgetc2.f"
	}
#line 216 "zgetc2.f"
	i__2 = *n - i__;
#line 216 "zgetc2.f"
	i__3 = *n - i__;
#line 216 "zgetc2.f"
	zgeru_(&i__2, &i__3, &c_b10, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[
		i__ + (i__ + 1) * a_dim1], lda, &a[i__ + 1 + (i__ + 1) * 
		a_dim1], lda);
#line 218 "zgetc2.f"
/* L40: */
#line 218 "zgetc2.f"
    }

#line 220 "zgetc2.f"
    if (z_abs(&a[*n + *n * a_dim1]) < smin) {
#line 221 "zgetc2.f"
	*info = *n;
#line 222 "zgetc2.f"
	i__1 = *n + *n * a_dim1;
#line 222 "zgetc2.f"
	z__1.r = smin, z__1.i = 0.;
#line 222 "zgetc2.f"
	a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 223 "zgetc2.f"
    }

/*     Set last pivots to N */

#line 227 "zgetc2.f"
    ipiv[*n] = *n;
#line 228 "zgetc2.f"
    jpiv[*n] = *n;

#line 230 "zgetc2.f"
    return 0;

/*     End of ZGETC2 */

} /* zgetc2_ */

