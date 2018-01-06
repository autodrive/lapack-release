#line 1 "zgetf2.f"
/* zgetf2.f -- translated by f2c (version 20100827).
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

#line 1 "zgetf2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row
 interchanges (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGETF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGETF2 computes an LU factorization of a general m-by-n matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > The factorization has the form */
/* >    A = P * L * U */
/* > where P is a permutation matrix, L is lower triangular with unit */
/* > diagonal elements (lower trapezoidal if m > n), and U is upper */
/* > triangular (upper trapezoidal if m < n). */
/* > */
/* > This is the right-looking Level 2 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the m by n matrix to be factored. */
/* >          On exit, the factors L and U from the factorization */
/* >          A = P*L*U; the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (min(M,N)) */
/* >          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* >          > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
/* >               has been completed, but the factor U is exactly */
/* >               singular, and division by zero will occur if it is used */
/* >               to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgetf2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, jp;
    static doublereal sfmin;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgeru_(integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zswap_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer izamax_(integer *, doublecomplex *, integer *);


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 150 "zgetf2.f"
    /* Parameter adjustments */
#line 150 "zgetf2.f"
    a_dim1 = *lda;
#line 150 "zgetf2.f"
    a_offset = 1 + a_dim1;
#line 150 "zgetf2.f"
    a -= a_offset;
#line 150 "zgetf2.f"
    --ipiv;
#line 150 "zgetf2.f"

#line 150 "zgetf2.f"
    /* Function Body */
#line 150 "zgetf2.f"
    *info = 0;
#line 151 "zgetf2.f"
    if (*m < 0) {
#line 152 "zgetf2.f"
	*info = -1;
#line 153 "zgetf2.f"
    } else if (*n < 0) {
#line 154 "zgetf2.f"
	*info = -2;
#line 155 "zgetf2.f"
    } else if (*lda < max(1,*m)) {
#line 156 "zgetf2.f"
	*info = -4;
#line 157 "zgetf2.f"
    }
#line 158 "zgetf2.f"
    if (*info != 0) {
#line 159 "zgetf2.f"
	i__1 = -(*info);
#line 159 "zgetf2.f"
	xerbla_("ZGETF2", &i__1, (ftnlen)6);
#line 160 "zgetf2.f"
	return 0;
#line 161 "zgetf2.f"
    }

/*     Quick return if possible */

#line 165 "zgetf2.f"
    if (*m == 0 || *n == 0) {
#line 165 "zgetf2.f"
	return 0;
#line 165 "zgetf2.f"
    }

/*     Compute machine safe minimum */

#line 170 "zgetf2.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 172 "zgetf2.f"
    i__1 = min(*m,*n);
#line 172 "zgetf2.f"
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

#line 176 "zgetf2.f"
	i__2 = *m - j + 1;
#line 176 "zgetf2.f"
	jp = j - 1 + izamax_(&i__2, &a[j + j * a_dim1], &c__1);
#line 177 "zgetf2.f"
	ipiv[j] = jp;
#line 178 "zgetf2.f"
	i__2 = jp + j * a_dim1;
#line 178 "zgetf2.f"
	if (a[i__2].r != 0. || a[i__2].i != 0.) {

/*           Apply the interchange to columns 1:N. */

#line 182 "zgetf2.f"
	    if (jp != j) {
#line 182 "zgetf2.f"
		zswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
#line 182 "zgetf2.f"
	    }

/*           Compute elements J+1:M of J-th column. */

#line 187 "zgetf2.f"
	    if (j < *m) {
#line 188 "zgetf2.f"
		if (z_abs(&a[j + j * a_dim1]) >= sfmin) {
#line 189 "zgetf2.f"
		    i__2 = *m - j;
#line 189 "zgetf2.f"
		    z_div(&z__1, &c_b1, &a[j + j * a_dim1]);
#line 189 "zgetf2.f"
		    zscal_(&i__2, &z__1, &a[j + 1 + j * a_dim1], &c__1);
#line 190 "zgetf2.f"
		} else {
#line 191 "zgetf2.f"
		    i__2 = *m - j;
#line 191 "zgetf2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 192 "zgetf2.f"
			i__3 = j + i__ + j * a_dim1;
#line 192 "zgetf2.f"
			z_div(&z__1, &a[j + i__ + j * a_dim1], &a[j + j * 
				a_dim1]);
#line 192 "zgetf2.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 193 "zgetf2.f"
/* L20: */
#line 193 "zgetf2.f"
		    }
#line 194 "zgetf2.f"
		}
#line 195 "zgetf2.f"
	    }

#line 197 "zgetf2.f"
	} else if (*info == 0) {

#line 199 "zgetf2.f"
	    *info = j;
#line 200 "zgetf2.f"
	}

#line 202 "zgetf2.f"
	if (j < min(*m,*n)) {

/*           Update trailing submatrix. */

#line 206 "zgetf2.f"
	    i__2 = *m - j;
#line 206 "zgetf2.f"
	    i__3 = *n - j;
#line 206 "zgetf2.f"
	    z__1.r = -1., z__1.i = -0.;
#line 206 "zgetf2.f"
	    zgeru_(&i__2, &i__3, &z__1, &a[j + 1 + j * a_dim1], &c__1, &a[j + 
		    (j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda)
		    ;
#line 208 "zgetf2.f"
	}
#line 209 "zgetf2.f"
/* L10: */
#line 209 "zgetf2.f"
    }
#line 210 "zgetf2.f"
    return 0;

/*     End of ZGETF2 */

} /* zgetf2_ */

