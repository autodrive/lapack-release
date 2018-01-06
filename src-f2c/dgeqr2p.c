#line 1 "dgeqr2p.f"
/* dgeqr2p.f -- translated by f2c (version 20100827).
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

#line 1 "dgeqr2p.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGEQR2P computes the QR factorization of a general rectangular matrix with non-negative diagona
l elements using an unblocked algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQR2P + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2p
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2p
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2p
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQR2P( M, N, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQR2 computes a QR factorization of a real m by n matrix A: */
/* > A = Q * R. The diagonal entries of R are nonnegative. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the m by n matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(m,n) by n upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n). The diagonal entries of R are */
/* >          nonnegative; the elements below the diagonal, */
/* >          with the array TAU, represent the orthogonal matrix Q as a */
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
/* >          TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \date November 2015 */

/* > \ingroup doubleGEcomputational */

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
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > */
/* > See Lapack Working Note 203 for details */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgeqr2p_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k;
    static doublereal aii;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    dlarfgp_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 159 "dgeqr2p.f"
    /* Parameter adjustments */
#line 159 "dgeqr2p.f"
    a_dim1 = *lda;
#line 159 "dgeqr2p.f"
    a_offset = 1 + a_dim1;
#line 159 "dgeqr2p.f"
    a -= a_offset;
#line 159 "dgeqr2p.f"
    --tau;
#line 159 "dgeqr2p.f"
    --work;
#line 159 "dgeqr2p.f"

#line 159 "dgeqr2p.f"
    /* Function Body */
#line 159 "dgeqr2p.f"
    *info = 0;
#line 160 "dgeqr2p.f"
    if (*m < 0) {
#line 161 "dgeqr2p.f"
	*info = -1;
#line 162 "dgeqr2p.f"
    } else if (*n < 0) {
#line 163 "dgeqr2p.f"
	*info = -2;
#line 164 "dgeqr2p.f"
    } else if (*lda < max(1,*m)) {
#line 165 "dgeqr2p.f"
	*info = -4;
#line 166 "dgeqr2p.f"
    }
#line 167 "dgeqr2p.f"
    if (*info != 0) {
#line 168 "dgeqr2p.f"
	i__1 = -(*info);
#line 168 "dgeqr2p.f"
	xerbla_("DGEQR2P", &i__1, (ftnlen)7);
#line 169 "dgeqr2p.f"
	return 0;
#line 170 "dgeqr2p.f"
    }

#line 172 "dgeqr2p.f"
    k = min(*m,*n);

#line 174 "dgeqr2p.f"
    i__1 = k;
#line 174 "dgeqr2p.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

#line 178 "dgeqr2p.f"
	i__2 = *m - i__ + 1;
/* Computing MIN */
#line 178 "dgeqr2p.f"
	i__3 = i__ + 1;
#line 178 "dgeqr2p.f"
	dlarfgp_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * 
		a_dim1], &c__1, &tau[i__]);
#line 180 "dgeqr2p.f"
	if (i__ < *n) {

/*           Apply H(i) to A(i:m,i+1:n) from the left */

#line 184 "dgeqr2p.f"
	    aii = a[i__ + i__ * a_dim1];
#line 185 "dgeqr2p.f"
	    a[i__ + i__ * a_dim1] = 1.;
#line 186 "dgeqr2p.f"
	    i__2 = *m - i__ + 1;
#line 186 "dgeqr2p.f"
	    i__3 = *n - i__;
#line 186 "dgeqr2p.f"
	    dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[
		    i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (
		    ftnlen)4);
#line 188 "dgeqr2p.f"
	    a[i__ + i__ * a_dim1] = aii;
#line 189 "dgeqr2p.f"
	}
#line 190 "dgeqr2p.f"
/* L10: */
#line 190 "dgeqr2p.f"
    }
#line 191 "dgeqr2p.f"
    return 0;

/*     End of DGEQR2P */

} /* dgeqr2p_ */

