#line 1 "zgetrf.f"
/* zgetrf.f -- translated by f2c (version 20100827).
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

#line 1 "zgetrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZGETRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGETRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO ) */

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
/* > ZGETRF computes an LU factorization of a general M-by-N matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > The factorization has the form */
/* >    A = P * L * U */
/* > where P is a permutation matrix, L is lower triangular with unit */
/* > diagonal elements (lower trapezoidal if m > n), and U is upper */
/* > triangular (upper trapezoidal if m < n). */
/* > */
/* > This is the right-looking Level 3 BLAS version of the algorithm. */
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
/* >          On entry, the M-by-N matrix to be factored. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
/* >                has been completed, but the factor U is exactly */
/* >                singular, and division by zero will occur if it is used */
/* >                to solve a system of equations. */
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
/* Subroutine */ int zgetrf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, jb, nb, iinfo;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), ztrsm_(char *, char *, char *, char *,
	     integer *, integer *, doublecomplex *, doublecomplex *, integer *
	    , doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    zgetf2_(integer *, integer *, doublecomplex *, integer *, integer 
	    *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 147 "zgetrf.f"
    /* Parameter adjustments */
#line 147 "zgetrf.f"
    a_dim1 = *lda;
#line 147 "zgetrf.f"
    a_offset = 1 + a_dim1;
#line 147 "zgetrf.f"
    a -= a_offset;
#line 147 "zgetrf.f"
    --ipiv;
#line 147 "zgetrf.f"

#line 147 "zgetrf.f"
    /* Function Body */
#line 147 "zgetrf.f"
    *info = 0;
#line 148 "zgetrf.f"
    if (*m < 0) {
#line 149 "zgetrf.f"
	*info = -1;
#line 150 "zgetrf.f"
    } else if (*n < 0) {
#line 151 "zgetrf.f"
	*info = -2;
#line 152 "zgetrf.f"
    } else if (*lda < max(1,*m)) {
#line 153 "zgetrf.f"
	*info = -4;
#line 154 "zgetrf.f"
    }
#line 155 "zgetrf.f"
    if (*info != 0) {
#line 156 "zgetrf.f"
	i__1 = -(*info);
#line 156 "zgetrf.f"
	xerbla_("ZGETRF", &i__1, (ftnlen)6);
#line 157 "zgetrf.f"
	return 0;
#line 158 "zgetrf.f"
    }

/*     Quick return if possible */

#line 162 "zgetrf.f"
    if (*m == 0 || *n == 0) {
#line 162 "zgetrf.f"
	return 0;
#line 162 "zgetrf.f"
    }

/*     Determine the block size for this environment. */

#line 167 "zgetrf.f"
    nb = ilaenv_(&c__1, "ZGETRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 168 "zgetrf.f"
    if (nb <= 1 || nb >= min(*m,*n)) {

/*        Use unblocked code. */

#line 172 "zgetrf.f"
	zgetf2_(m, n, &a[a_offset], lda, &ipiv[1], info);
#line 173 "zgetrf.f"
    } else {

/*        Use blocked code. */

#line 177 "zgetrf.f"
	i__1 = min(*m,*n);
#line 177 "zgetrf.f"
	i__2 = nb;
#line 177 "zgetrf.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 178 "zgetrf.f"
	    i__3 = min(*m,*n) - j + 1;
#line 178 "zgetrf.f"
	    jb = min(i__3,nb);

/*           Factor diagonal and subdiagonal blocks and test for exact */
/*           singularity. */

#line 183 "zgetrf.f"
	    i__3 = *m - j + 1;
#line 183 "zgetrf.f"
	    zgetf2_(&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

#line 187 "zgetrf.f"
	    if (*info == 0 && iinfo > 0) {
#line 187 "zgetrf.f"
		*info = iinfo + j - 1;
#line 187 "zgetrf.f"
	    }
/* Computing MIN */
#line 189 "zgetrf.f"
	    i__4 = *m, i__5 = j + jb - 1;
#line 189 "zgetrf.f"
	    i__3 = min(i__4,i__5);
#line 189 "zgetrf.f"
	    for (i__ = j; i__ <= i__3; ++i__) {
#line 190 "zgetrf.f"
		ipiv[i__] = j - 1 + ipiv[i__];
#line 191 "zgetrf.f"
/* L10: */
#line 191 "zgetrf.f"
	    }

/*           Apply interchanges to columns 1:J-1. */

#line 195 "zgetrf.f"
	    i__3 = j - 1;
#line 195 "zgetrf.f"
	    i__4 = j + jb - 1;
#line 195 "zgetrf.f"
	    zlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

#line 197 "zgetrf.f"
	    if (j + jb <= *n) {

/*              Apply interchanges to columns J+JB:N. */

#line 201 "zgetrf.f"
		i__3 = *n - j - jb + 1;
#line 201 "zgetrf.f"
		i__4 = j + jb - 1;
#line 201 "zgetrf.f"
		zlaswp_(&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &
			ipiv[1], &c__1);

/*              Compute block row of U. */

#line 206 "zgetrf.f"
		i__3 = *n - j - jb + 1;
#line 206 "zgetrf.f"
		ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &
			c_b1, &a[j + j * a_dim1], lda, &a[j + (j + jb) * 
			a_dim1], lda, (ftnlen)4, (ftnlen)5, (ftnlen)12, (
			ftnlen)4);
#line 209 "zgetrf.f"
		if (j + jb <= *m) {

/*                 Update trailing submatrix. */

#line 213 "zgetrf.f"
		    i__3 = *m - j - jb + 1;
#line 213 "zgetrf.f"
		    i__4 = *n - j - jb + 1;
#line 213 "zgetrf.f"
		    z__1.r = -1., z__1.i = -0.;
#line 213 "zgetrf.f"
		    zgemm_("No transpose", "No transpose", &i__3, &i__4, &jb, 
			    &z__1, &a[j + jb + j * a_dim1], lda, &a[j + (j + 
			    jb) * a_dim1], lda, &c_b1, &a[j + jb + (j + jb) * 
			    a_dim1], lda, (ftnlen)12, (ftnlen)12);
#line 217 "zgetrf.f"
		}
#line 218 "zgetrf.f"
	    }
#line 219 "zgetrf.f"
/* L20: */
#line 219 "zgetrf.f"
	}
#line 220 "zgetrf.f"
    }
#line 221 "zgetrf.f"
    return 0;

/*     End of ZGETRF */

} /* zgetrf_ */

