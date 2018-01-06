#line 1 "cgeqr2.f"
/* cgeqr2.f -- translated by f2c (version 20100827).
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

#line 1 "cgeqr2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorit
hm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEQR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQR2( M, N, A, LDA, TAU, WORK, INFO ) */

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
/* > CGEQR2 computes a QR factorization of a complex m by n matrix A: */
/* > A = Q * R. */
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
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(m,n) by n upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n); the elements below the diagonal, */
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
/* >          WORK is COMPLEX array, dimension (N) */
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

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgeqr2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k;
    static doublecomplex alpha;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), clarfg_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *), 
	    xerbla_(char *, integer *, ftnlen);


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

#line 156 "cgeqr2.f"
    /* Parameter adjustments */
#line 156 "cgeqr2.f"
    a_dim1 = *lda;
#line 156 "cgeqr2.f"
    a_offset = 1 + a_dim1;
#line 156 "cgeqr2.f"
    a -= a_offset;
#line 156 "cgeqr2.f"
    --tau;
#line 156 "cgeqr2.f"
    --work;
#line 156 "cgeqr2.f"

#line 156 "cgeqr2.f"
    /* Function Body */
#line 156 "cgeqr2.f"
    *info = 0;
#line 157 "cgeqr2.f"
    if (*m < 0) {
#line 158 "cgeqr2.f"
	*info = -1;
#line 159 "cgeqr2.f"
    } else if (*n < 0) {
#line 160 "cgeqr2.f"
	*info = -2;
#line 161 "cgeqr2.f"
    } else if (*lda < max(1,*m)) {
#line 162 "cgeqr2.f"
	*info = -4;
#line 163 "cgeqr2.f"
    }
#line 164 "cgeqr2.f"
    if (*info != 0) {
#line 165 "cgeqr2.f"
	i__1 = -(*info);
#line 165 "cgeqr2.f"
	xerbla_("CGEQR2", &i__1, (ftnlen)6);
#line 166 "cgeqr2.f"
	return 0;
#line 167 "cgeqr2.f"
    }

#line 169 "cgeqr2.f"
    k = min(*m,*n);

#line 171 "cgeqr2.f"
    i__1 = k;
#line 171 "cgeqr2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

#line 175 "cgeqr2.f"
	i__2 = *m - i__ + 1;
/* Computing MIN */
#line 175 "cgeqr2.f"
	i__3 = i__ + 1;
#line 175 "cgeqr2.f"
	clarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * a_dim1]
		, &c__1, &tau[i__]);
#line 177 "cgeqr2.f"
	if (i__ < *n) {

/*           Apply H(i)**H to A(i:m,i+1:n) from the left */

#line 181 "cgeqr2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 181 "cgeqr2.f"
	    alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 182 "cgeqr2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 182 "cgeqr2.f"
	    a[i__2].r = 1., a[i__2].i = 0.;
#line 183 "cgeqr2.f"
	    i__2 = *m - i__ + 1;
#line 183 "cgeqr2.f"
	    i__3 = *n - i__;
#line 183 "cgeqr2.f"
	    d_cnjg(&z__1, &tau[i__]);
#line 183 "cgeqr2.f"
	    clarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &z__1,
		     &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)4);
#line 185 "cgeqr2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 185 "cgeqr2.f"
	    a[i__2].r = alpha.r, a[i__2].i = alpha.i;
#line 186 "cgeqr2.f"
	}
#line 187 "cgeqr2.f"
/* L10: */
#line 187 "cgeqr2.f"
    }
#line 188 "cgeqr2.f"
    return 0;

/*     End of CGEQR2 */

} /* cgeqr2_ */

