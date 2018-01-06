#line 1 "zgerq2.f"
/* zgerq2.f -- translated by f2c (version 20100827).
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

#line 1 "zgerq2.f"
/* > \brief \b ZGERQ2 computes the RQ factorization of a general rectangular matrix using an unblocked algorit
hm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGERQ2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgerq2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgerq2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgerq2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGERQ2( M, N, A, LDA, TAU, WORK, INFO ) */

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
/* > ZGERQ2 computes an RQ factorization of a complex m by n matrix A: */
/* > A = R * Q. */
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
/* >          On exit, if m <= n, the upper triangle of the subarray */
/* >          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R; */
/* >          if m >= n, the elements on and above the (m-n)-th subdiagonal */
/* >          contain the m by n upper trapezoidal matrix R; the remaining */
/* >          elements, with the array TAU, represent the unitary matrix */
/* >          Q as a product of elementary reflectors (see Further */
/* >          Details). */
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
/* >          WORK is COMPLEX*16 array, dimension (M) */
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

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1)**H H(2)**H . . . H(k)**H, where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on */
/* >  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgerq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublecomplex alpha;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), zlarfg_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), zlacgv_(integer *, doublecomplex *, 
	    integer *);


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

#line 158 "zgerq2.f"
    /* Parameter adjustments */
#line 158 "zgerq2.f"
    a_dim1 = *lda;
#line 158 "zgerq2.f"
    a_offset = 1 + a_dim1;
#line 158 "zgerq2.f"
    a -= a_offset;
#line 158 "zgerq2.f"
    --tau;
#line 158 "zgerq2.f"
    --work;
#line 158 "zgerq2.f"

#line 158 "zgerq2.f"
    /* Function Body */
#line 158 "zgerq2.f"
    *info = 0;
#line 159 "zgerq2.f"
    if (*m < 0) {
#line 160 "zgerq2.f"
	*info = -1;
#line 161 "zgerq2.f"
    } else if (*n < 0) {
#line 162 "zgerq2.f"
	*info = -2;
#line 163 "zgerq2.f"
    } else if (*lda < max(1,*m)) {
#line 164 "zgerq2.f"
	*info = -4;
#line 165 "zgerq2.f"
    }
#line 166 "zgerq2.f"
    if (*info != 0) {
#line 167 "zgerq2.f"
	i__1 = -(*info);
#line 167 "zgerq2.f"
	xerbla_("ZGERQ2", &i__1, (ftnlen)6);
#line 168 "zgerq2.f"
	return 0;
#line 169 "zgerq2.f"
    }

#line 171 "zgerq2.f"
    k = min(*m,*n);

#line 173 "zgerq2.f"
    for (i__ = k; i__ >= 1; --i__) {

/*        Generate elementary reflector H(i) to annihilate */
/*        A(m-k+i,1:n-k+i-1) */

#line 178 "zgerq2.f"
	i__1 = *n - k + i__;
#line 178 "zgerq2.f"
	zlacgv_(&i__1, &a[*m - k + i__ + a_dim1], lda);
#line 179 "zgerq2.f"
	i__1 = *m - k + i__ + (*n - k + i__) * a_dim1;
#line 179 "zgerq2.f"
	alpha.r = a[i__1].r, alpha.i = a[i__1].i;
#line 180 "zgerq2.f"
	i__1 = *n - k + i__;
#line 180 "zgerq2.f"
	zlarfg_(&i__1, &alpha, &a[*m - k + i__ + a_dim1], lda, &tau[i__]);

/*        Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right */

#line 184 "zgerq2.f"
	i__1 = *m - k + i__ + (*n - k + i__) * a_dim1;
#line 184 "zgerq2.f"
	a[i__1].r = 1., a[i__1].i = 0.;
#line 185 "zgerq2.f"
	i__1 = *m - k + i__ - 1;
#line 185 "zgerq2.f"
	i__2 = *n - k + i__;
#line 185 "zgerq2.f"
	zlarf_("Right", &i__1, &i__2, &a[*m - k + i__ + a_dim1], lda, &tau[
		i__], &a[a_offset], lda, &work[1], (ftnlen)5);
#line 187 "zgerq2.f"
	i__1 = *m - k + i__ + (*n - k + i__) * a_dim1;
#line 187 "zgerq2.f"
	a[i__1].r = alpha.r, a[i__1].i = alpha.i;
#line 188 "zgerq2.f"
	i__1 = *n - k + i__ - 1;
#line 188 "zgerq2.f"
	zlacgv_(&i__1, &a[*m - k + i__ + a_dim1], lda);
#line 189 "zgerq2.f"
/* L10: */
#line 189 "zgerq2.f"
    }
#line 190 "zgerq2.f"
    return 0;

/*     End of ZGERQ2 */

} /* zgerq2_ */

