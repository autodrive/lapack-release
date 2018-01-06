#line 1 "cgelq2.f"
/* cgelq2.f -- translated by f2c (version 20100827).
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

#line 1 "cgelq2.f"
/* > \brief \b CGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorit
hm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGELQ2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelq2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelq2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelq2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELQ2( M, N, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGELQ2 computes an LQ factorization of a complex m by n matrix A: */
/* > A = L * Q. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the m by n matrix A. */
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the m by min(m,n) lower trapezoidal matrix L (L is */
/* >          lower triangular if m <= n); the elements above the diagonal, */
/* >          with the array TAU, represent the unitary matrix Q as a */
/* >          product of elementary reflectors (see Further Details). */
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
/* >          TAU is COMPLEX array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (M) */
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

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(k)**H . . . H(2)**H H(1)**H, where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in */
/* >  A(i,i+1:n), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgelq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k;
    static doublecomplex alpha;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), clarfg_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *), 
	    clacgv_(integer *, doublecomplex *, integer *), xerbla_(char *, 
	    integer *, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 156 "cgelq2.f"
    /* Parameter adjustments */
#line 156 "cgelq2.f"
    a_dim1 = *lda;
#line 156 "cgelq2.f"
    a_offset = 1 + a_dim1;
#line 156 "cgelq2.f"
    a -= a_offset;
#line 156 "cgelq2.f"
    --tau;
#line 156 "cgelq2.f"
    --work;
#line 156 "cgelq2.f"

#line 156 "cgelq2.f"
    /* Function Body */
#line 156 "cgelq2.f"
    *info = 0;
#line 157 "cgelq2.f"
    if (*m < 0) {
#line 158 "cgelq2.f"
	*info = -1;
#line 159 "cgelq2.f"
    } else if (*n < 0) {
#line 160 "cgelq2.f"
	*info = -2;
#line 161 "cgelq2.f"
    } else if (*lda < max(1,*m)) {
#line 162 "cgelq2.f"
	*info = -4;
#line 163 "cgelq2.f"
    }
#line 164 "cgelq2.f"
    if (*info != 0) {
#line 165 "cgelq2.f"
	i__1 = -(*info);
#line 165 "cgelq2.f"
	xerbla_("CGELQ2", &i__1, (ftnlen)6);
#line 166 "cgelq2.f"
	return 0;
#line 167 "cgelq2.f"
    }

#line 169 "cgelq2.f"
    k = min(*m,*n);

#line 171 "cgelq2.f"
    i__1 = k;
#line 171 "cgelq2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i,i+1:n) */

#line 175 "cgelq2.f"
	i__2 = *n - i__ + 1;
#line 175 "cgelq2.f"
	clacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);
#line 176 "cgelq2.f"
	i__2 = i__ + i__ * a_dim1;
#line 176 "cgelq2.f"
	alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 177 "cgelq2.f"
	i__2 = *n - i__ + 1;
/* Computing MIN */
#line 177 "cgelq2.f"
	i__3 = i__ + 1;
#line 177 "cgelq2.f"
	clarfg_(&i__2, &alpha, &a[i__ + min(i__3,*n) * a_dim1], lda, &tau[i__]
		);
#line 179 "cgelq2.f"
	if (i__ < *m) {

/*           Apply H(i) to A(i+1:m,i:n) from the right */

#line 183 "cgelq2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 183 "cgelq2.f"
	    a[i__2].r = 1., a[i__2].i = 0.;
#line 184 "cgelq2.f"
	    i__2 = *m - i__;
#line 184 "cgelq2.f"
	    i__3 = *n - i__ + 1;
#line 184 "cgelq2.f"
	    clarf_("Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[
		    i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1], (ftnlen)
		    5);
#line 186 "cgelq2.f"
	}
#line 187 "cgelq2.f"
	i__2 = i__ + i__ * a_dim1;
#line 187 "cgelq2.f"
	a[i__2].r = alpha.r, a[i__2].i = alpha.i;
#line 188 "cgelq2.f"
	i__2 = *n - i__ + 1;
#line 188 "cgelq2.f"
	clacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);
#line 189 "cgelq2.f"
/* L10: */
#line 189 "cgelq2.f"
    }
#line 190 "cgelq2.f"
    return 0;

/*     End of CGELQ2 */

} /* cgelq2_ */

