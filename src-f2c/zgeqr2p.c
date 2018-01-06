#line 1 "zgeqr2p.f"
/* zgeqr2p.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqr2p.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZGEQR2P computes the QR factorization of a general rectangular matrix with non-negative diagona
l elements using an unblocked algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQR2P + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqr2p
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqr2p
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqr2p
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQR2P( M, N, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQR2P computes a QR factorization of a complex m by n matrix A: */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          TAU is COMPLEX*16 array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
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

/* > \ingroup complex16GEcomputational */

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
/* Subroutine */ int zgeqr2p_(integer *m, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), zlarfgp_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *);


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

#line 156 "zgeqr2p.f"
    /* Parameter adjustments */
#line 156 "zgeqr2p.f"
    a_dim1 = *lda;
#line 156 "zgeqr2p.f"
    a_offset = 1 + a_dim1;
#line 156 "zgeqr2p.f"
    a -= a_offset;
#line 156 "zgeqr2p.f"
    --tau;
#line 156 "zgeqr2p.f"
    --work;
#line 156 "zgeqr2p.f"

#line 156 "zgeqr2p.f"
    /* Function Body */
#line 156 "zgeqr2p.f"
    *info = 0;
#line 157 "zgeqr2p.f"
    if (*m < 0) {
#line 158 "zgeqr2p.f"
	*info = -1;
#line 159 "zgeqr2p.f"
    } else if (*n < 0) {
#line 160 "zgeqr2p.f"
	*info = -2;
#line 161 "zgeqr2p.f"
    } else if (*lda < max(1,*m)) {
#line 162 "zgeqr2p.f"
	*info = -4;
#line 163 "zgeqr2p.f"
    }
#line 164 "zgeqr2p.f"
    if (*info != 0) {
#line 165 "zgeqr2p.f"
	i__1 = -(*info);
#line 165 "zgeqr2p.f"
	xerbla_("ZGEQR2P", &i__1, (ftnlen)7);
#line 166 "zgeqr2p.f"
	return 0;
#line 167 "zgeqr2p.f"
    }

#line 169 "zgeqr2p.f"
    k = min(*m,*n);

#line 171 "zgeqr2p.f"
    i__1 = k;
#line 171 "zgeqr2p.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

#line 175 "zgeqr2p.f"
	i__2 = *m - i__ + 1;
/* Computing MIN */
#line 175 "zgeqr2p.f"
	i__3 = i__ + 1;
#line 175 "zgeqr2p.f"
	zlarfgp_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * 
		a_dim1], &c__1, &tau[i__]);
#line 177 "zgeqr2p.f"
	if (i__ < *n) {

/*           Apply H(i)**H to A(i:m,i+1:n) from the left */

#line 181 "zgeqr2p.f"
	    i__2 = i__ + i__ * a_dim1;
#line 181 "zgeqr2p.f"
	    alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 182 "zgeqr2p.f"
	    i__2 = i__ + i__ * a_dim1;
#line 182 "zgeqr2p.f"
	    a[i__2].r = 1., a[i__2].i = 0.;
#line 183 "zgeqr2p.f"
	    i__2 = *m - i__ + 1;
#line 183 "zgeqr2p.f"
	    i__3 = *n - i__;
#line 183 "zgeqr2p.f"
	    d_cnjg(&z__1, &tau[i__]);
#line 183 "zgeqr2p.f"
	    zlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &z__1,
		     &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)4);
#line 185 "zgeqr2p.f"
	    i__2 = i__ + i__ * a_dim1;
#line 185 "zgeqr2p.f"
	    a[i__2].r = alpha.r, a[i__2].i = alpha.i;
#line 186 "zgeqr2p.f"
	}
#line 187 "zgeqr2p.f"
/* L10: */
#line 187 "zgeqr2p.f"
    }
#line 188 "zgeqr2p.f"
    return 0;

/*     End of ZGEQR2P */

} /* zgeqr2p_ */

