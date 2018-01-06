#line 1 "sgetf2.f"
/* sgetf2.f -- translated by f2c (version 20100827).
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

#line 1 "sgetf2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = -1.;

/* > \brief \b SGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row
 interchanges (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGETF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGETF2 computes an LU factorization of a general m-by-n matrix A */
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
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgetf2_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, jp;
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), sscal_(integer *, doublereal *, doublereal *, integer 
	    *);
    static doublereal sfmin;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);


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

#line 149 "sgetf2.f"
    /* Parameter adjustments */
#line 149 "sgetf2.f"
    a_dim1 = *lda;
#line 149 "sgetf2.f"
    a_offset = 1 + a_dim1;
#line 149 "sgetf2.f"
    a -= a_offset;
#line 149 "sgetf2.f"
    --ipiv;
#line 149 "sgetf2.f"

#line 149 "sgetf2.f"
    /* Function Body */
#line 149 "sgetf2.f"
    *info = 0;
#line 150 "sgetf2.f"
    if (*m < 0) {
#line 151 "sgetf2.f"
	*info = -1;
#line 152 "sgetf2.f"
    } else if (*n < 0) {
#line 153 "sgetf2.f"
	*info = -2;
#line 154 "sgetf2.f"
    } else if (*lda < max(1,*m)) {
#line 155 "sgetf2.f"
	*info = -4;
#line 156 "sgetf2.f"
    }
#line 157 "sgetf2.f"
    if (*info != 0) {
#line 158 "sgetf2.f"
	i__1 = -(*info);
#line 158 "sgetf2.f"
	xerbla_("SGETF2", &i__1, (ftnlen)6);
#line 159 "sgetf2.f"
	return 0;
#line 160 "sgetf2.f"
    }

/*     Quick return if possible */

#line 164 "sgetf2.f"
    if (*m == 0 || *n == 0) {
#line 164 "sgetf2.f"
	return 0;
#line 164 "sgetf2.f"
    }

/*     Compute machine safe minimum */

#line 169 "sgetf2.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 171 "sgetf2.f"
    i__1 = min(*m,*n);
#line 171 "sgetf2.f"
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

#line 175 "sgetf2.f"
	i__2 = *m - j + 1;
#line 175 "sgetf2.f"
	jp = j - 1 + isamax_(&i__2, &a[j + j * a_dim1], &c__1);
#line 176 "sgetf2.f"
	ipiv[j] = jp;
#line 177 "sgetf2.f"
	if (a[jp + j * a_dim1] != 0.) {

/*           Apply the interchange to columns 1:N. */

#line 181 "sgetf2.f"
	    if (jp != j) {
#line 181 "sgetf2.f"
		sswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
#line 181 "sgetf2.f"
	    }

/*           Compute elements J+1:M of J-th column. */

#line 186 "sgetf2.f"
	    if (j < *m) {
#line 187 "sgetf2.f"
		if ((d__1 = a[j + j * a_dim1], abs(d__1)) >= sfmin) {
#line 188 "sgetf2.f"
		    i__2 = *m - j;
#line 188 "sgetf2.f"
		    d__1 = 1. / a[j + j * a_dim1];
#line 188 "sgetf2.f"
		    sscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
#line 189 "sgetf2.f"
		} else {
#line 190 "sgetf2.f"
		    i__2 = *m - j;
#line 190 "sgetf2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 191 "sgetf2.f"
			a[j + i__ + j * a_dim1] /= a[j + j * a_dim1];
#line 192 "sgetf2.f"
/* L20: */
#line 192 "sgetf2.f"
		    }
#line 193 "sgetf2.f"
		}
#line 194 "sgetf2.f"
	    }

#line 196 "sgetf2.f"
	} else if (*info == 0) {

#line 198 "sgetf2.f"
	    *info = j;
#line 199 "sgetf2.f"
	}

#line 201 "sgetf2.f"
	if (j < min(*m,*n)) {

/*           Update trailing submatrix. */

#line 205 "sgetf2.f"
	    i__2 = *m - j;
#line 205 "sgetf2.f"
	    i__3 = *n - j;
#line 205 "sgetf2.f"
	    sger_(&i__2, &i__3, &c_b8, &a[j + 1 + j * a_dim1], &c__1, &a[j + (
		    j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda);
#line 207 "sgetf2.f"
	}
#line 208 "sgetf2.f"
/* L10: */
#line 208 "sgetf2.f"
    }
#line 209 "sgetf2.f"
    return 0;

/*     End of SGETF2 */

} /* sgetf2_ */

