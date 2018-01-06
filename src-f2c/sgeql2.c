#line 1 "sgeql2.f"
/* sgeql2.f -- translated by f2c (version 20100827).
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

#line 1 "sgeql2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SGEQL2 computes the QL factorization of a general rectangular matrix using an unblocked algorit
hm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEQL2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeql2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeql2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeql2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEQL2( M, N, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQL2 computes a QL factorization of a real m by n matrix A: */
/* > A = Q * L. */
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
/* >          On entry, the m by n matrix A. */
/* >          On exit, if m >= n, the lower triangle of the subarray */
/* >          A(m-n+1:m,1:n) contains the n by n lower triangular matrix L; */
/* >          if m <= n, the elements on and below the (n-m)-th */
/* >          superdiagonal contain the m by n lower trapezoidal matrix L; */
/* >          the remaining elements, with the array TAU, represent the */
/* >          orthogonal matrix Q as a product of elementary reflectors */
/* >          (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(k) . . . H(2) H(1), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in */
/* >  A(1:m-k+i-1,n-k+i), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgeql2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal aii;
    extern /* Subroutine */ int slarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    slarfg_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 158 "sgeql2.f"
    /* Parameter adjustments */
#line 158 "sgeql2.f"
    a_dim1 = *lda;
#line 158 "sgeql2.f"
    a_offset = 1 + a_dim1;
#line 158 "sgeql2.f"
    a -= a_offset;
#line 158 "sgeql2.f"
    --tau;
#line 158 "sgeql2.f"
    --work;
#line 158 "sgeql2.f"

#line 158 "sgeql2.f"
    /* Function Body */
#line 158 "sgeql2.f"
    *info = 0;
#line 159 "sgeql2.f"
    if (*m < 0) {
#line 160 "sgeql2.f"
	*info = -1;
#line 161 "sgeql2.f"
    } else if (*n < 0) {
#line 162 "sgeql2.f"
	*info = -2;
#line 163 "sgeql2.f"
    } else if (*lda < max(1,*m)) {
#line 164 "sgeql2.f"
	*info = -4;
#line 165 "sgeql2.f"
    }
#line 166 "sgeql2.f"
    if (*info != 0) {
#line 167 "sgeql2.f"
	i__1 = -(*info);
#line 167 "sgeql2.f"
	xerbla_("SGEQL2", &i__1, (ftnlen)6);
#line 168 "sgeql2.f"
	return 0;
#line 169 "sgeql2.f"
    }

#line 171 "sgeql2.f"
    k = min(*m,*n);

#line 173 "sgeql2.f"
    for (i__ = k; i__ >= 1; --i__) {

/*        Generate elementary reflector H(i) to annihilate */
/*        A(1:m-k+i-1,n-k+i) */

#line 178 "sgeql2.f"
	i__1 = *m - k + i__;
#line 178 "sgeql2.f"
	slarfg_(&i__1, &a[*m - k + i__ + (*n - k + i__) * a_dim1], &a[(*n - k 
		+ i__) * a_dim1 + 1], &c__1, &tau[i__]);

/*        Apply H(i) to A(1:m-k+i,1:n-k+i-1) from the left */

#line 183 "sgeql2.f"
	aii = a[*m - k + i__ + (*n - k + i__) * a_dim1];
#line 184 "sgeql2.f"
	a[*m - k + i__ + (*n - k + i__) * a_dim1] = 1.;
#line 185 "sgeql2.f"
	i__1 = *m - k + i__;
#line 185 "sgeql2.f"
	i__2 = *n - k + i__ - 1;
#line 185 "sgeql2.f"
	slarf_("Left", &i__1, &i__2, &a[(*n - k + i__) * a_dim1 + 1], &c__1, &
		tau[i__], &a[a_offset], lda, &work[1], (ftnlen)4);
#line 187 "sgeql2.f"
	a[*m - k + i__ + (*n - k + i__) * a_dim1] = aii;
#line 188 "sgeql2.f"
/* L10: */
#line 188 "sgeql2.f"
    }
#line 189 "sgeql2.f"
    return 0;

/*     End of SGEQL2 */

} /* sgeql2_ */

