#line 1 "slahrd.f"
/* slahrd.f -- translated by f2c (version 20100827).
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

#line 1 "slahrd.f"
/* Table of constant values */

static doublereal c_b4 = -1.;
static doublereal c_b5 = 1.;
static integer c__1 = 1;
static doublereal c_b38 = 0.;

/* > \brief \b SLAHRD reduces the first nb columns of a general rectangular matrix A so that elements below th
e k-th subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformati
on to the unreduced part of A. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            K, LDA, LDT, LDY, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), T( LDT, NB ), TAU( NB ), */
/*      $                   Y( LDY, NB ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAHRD reduces the first NB columns of a real general n-by-(n-k+1) */
/* > matrix A so that elements below the k-th subdiagonal are zero. The */
/* > reduction is performed by an orthogonal similarity transformation */
/* > Q**T * A * Q. The routine returns the matrices V and T which determine */
/* > Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T. */
/* > */
/* > This is an OBSOLETE auxiliary routine. */
/* > This routine will be 'deprecated' in a  future release. */
/* > Please use the new routine SLAHR2 instead. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The offset for the reduction. Elements below the k-th */
/* >          subdiagonal in the first NB columns are reduced to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The number of columns to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N-K+1) */
/* >          On entry, the n-by-(n-k+1) general matrix A. */
/* >          On exit, the elements on and above the k-th subdiagonal in */
/* >          the first NB columns are overwritten with the corresponding */
/* >          elements of the reduced matrix; the elements below the k-th */
/* >          subdiagonal, with the array TAU, represent the matrix Q as a */
/* >          product of elementary reflectors. The other columns of A are */
/* >          unchanged. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (NB) */
/* >          The scalar factors of the elementary reflectors. See Further */
/* >          Details. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,NB) */
/* >          The upper triangular matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* >          Y is REAL array, dimension (LDY,NB) */
/* >          The n-by-nb matrix Y. */
/* > \endverbatim */
/* > */
/* > \param[in] LDY */
/* > \verbatim */
/* >          LDY is INTEGER */
/* >          The leading dimension of the array Y. LDY >= N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of nb elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(nb). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in */
/* >  A(i+k+1:n,i), and tau in TAU(i). */
/* > */
/* >  The elements of the vectors v together form the (n-k+1)-by-nb matrix */
/* >  V which is needed, with T and Y, to apply the transformation to the */
/* >  unreduced part of the matrix, using an update of the form: */
/* >  A := (I - V*T*V**T) * (A - Y*V**T). */
/* > */
/* >  The contents of A on exit are illustrated by the following example */
/* >  with n = 7, k = 3 and nb = 2: */
/* > */
/* >     ( a   h   a   a   a ) */
/* >     ( a   h   a   a   a ) */
/* >     ( a   h   a   a   a ) */
/* >     ( h   h   a   a   a ) */
/* >     ( v1  h   a   a   a ) */
/* >     ( v1  v2  a   a   a ) */
/* >     ( v1  v2  a   a   a ) */
/* > */
/* >  where a denotes an element of the original matrix A, h denotes a */
/* >  modified element of the upper Hessenberg matrix H, and vi denotes an */
/* >  element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slahrd_(integer *n, integer *k, integer *nb, doublereal *
	a, integer *lda, doublereal *tau, doublereal *t, integer *ldt, 
	doublereal *y, integer *ldy)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal ei;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), scopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), saxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *), strmv_(char 
	    *, char *, char *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen, ftnlen, ftnlen), slarfg_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
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

/*     Quick return if possible */

#line 205 "slahrd.f"
    /* Parameter adjustments */
#line 205 "slahrd.f"
    --tau;
#line 205 "slahrd.f"
    a_dim1 = *lda;
#line 205 "slahrd.f"
    a_offset = 1 + a_dim1;
#line 205 "slahrd.f"
    a -= a_offset;
#line 205 "slahrd.f"
    t_dim1 = *ldt;
#line 205 "slahrd.f"
    t_offset = 1 + t_dim1;
#line 205 "slahrd.f"
    t -= t_offset;
#line 205 "slahrd.f"
    y_dim1 = *ldy;
#line 205 "slahrd.f"
    y_offset = 1 + y_dim1;
#line 205 "slahrd.f"
    y -= y_offset;
#line 205 "slahrd.f"

#line 205 "slahrd.f"
    /* Function Body */
#line 205 "slahrd.f"
    if (*n <= 1) {
#line 205 "slahrd.f"
	return 0;
#line 205 "slahrd.f"
    }

#line 208 "slahrd.f"
    i__1 = *nb;
#line 208 "slahrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "slahrd.f"
	if (i__ > 1) {

/*           Update A(1:n,i) */

/*           Compute i-th column of A - Y * V**T */

#line 215 "slahrd.f"
	    i__2 = i__ - 1;
#line 215 "slahrd.f"
	    sgemv_("No transpose", n, &i__2, &c_b4, &y[y_offset], ldy, &a[*k 
		    + i__ - 1 + a_dim1], lda, &c_b5, &a[i__ * a_dim1 + 1], &
		    c__1, (ftnlen)12);

/*           Apply I - V * T**T * V**T to this column (call it b) from the */
/*           left, using the last column of T as workspace */

/*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows) */
/*                    ( V2 )             ( b2 ) */

/*           where V1 is unit lower triangular */

/*           w := V1**T * b1 */

#line 228 "slahrd.f"
	    i__2 = i__ - 1;
#line 228 "slahrd.f"
	    scopy_(&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 
		    1], &c__1);
#line 229 "slahrd.f"
	    i__2 = i__ - 1;
#line 229 "slahrd.f"
	    strmv_("Lower", "Transpose", "Unit", &i__2, &a[*k + 1 + a_dim1], 
		    lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9, (
		    ftnlen)4);

/*           w := w + V2**T *b2 */

#line 234 "slahrd.f"
	    i__2 = *n - *k - i__ + 1;
#line 234 "slahrd.f"
	    i__3 = i__ - 1;
#line 234 "slahrd.f"
	    sgemv_("Transpose", &i__2, &i__3, &c_b5, &a[*k + i__ + a_dim1], 
		    lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b5, &t[*nb * 
		    t_dim1 + 1], &c__1, (ftnlen)9);

/*           w := T**T *w */

#line 239 "slahrd.f"
	    i__2 = i__ - 1;
#line 239 "slahrd.f"
	    strmv_("Upper", "Transpose", "Non-unit", &i__2, &t[t_offset], ldt,
		     &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);

/*           b2 := b2 - V2*w */

#line 244 "slahrd.f"
	    i__2 = *n - *k - i__ + 1;
#line 244 "slahrd.f"
	    i__3 = i__ - 1;
#line 244 "slahrd.f"
	    sgemv_("No transpose", &i__2, &i__3, &c_b4, &a[*k + i__ + a_dim1],
		     lda, &t[*nb * t_dim1 + 1], &c__1, &c_b5, &a[*k + i__ + 
		    i__ * a_dim1], &c__1, (ftnlen)12);

/*           b1 := b1 - V1*w */

#line 249 "slahrd.f"
	    i__2 = i__ - 1;
#line 249 "slahrd.f"
	    strmv_("Lower", "No transpose", "Unit", &i__2, &a[*k + 1 + a_dim1]
		    , lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12,
		     (ftnlen)4);
#line 251 "slahrd.f"
	    i__2 = i__ - 1;
#line 251 "slahrd.f"
	    saxpy_(&i__2, &c_b4, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ 
		    * a_dim1], &c__1);

#line 253 "slahrd.f"
	    a[*k + i__ - 1 + (i__ - 1) * a_dim1] = ei;
#line 254 "slahrd.f"
	}

/*        Generate the elementary reflector H(i) to annihilate */
/*        A(k+i+1:n,i) */

#line 259 "slahrd.f"
	i__2 = *n - *k - i__ + 1;
/* Computing MIN */
#line 259 "slahrd.f"
	i__3 = *k + i__ + 1;
#line 259 "slahrd.f"
	slarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[min(i__3,*n) + i__ * 
		a_dim1], &c__1, &tau[i__]);
#line 261 "slahrd.f"
	ei = a[*k + i__ + i__ * a_dim1];
#line 262 "slahrd.f"
	a[*k + i__ + i__ * a_dim1] = 1.;

/*        Compute  Y(1:n,i) */

#line 266 "slahrd.f"
	i__2 = *n - *k - i__ + 1;
#line 266 "slahrd.f"
	sgemv_("No transpose", n, &i__2, &c_b5, &a[(i__ + 1) * a_dim1 + 1], 
		lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b38, &y[i__ * 
		y_dim1 + 1], &c__1, (ftnlen)12);
#line 268 "slahrd.f"
	i__2 = *n - *k - i__ + 1;
#line 268 "slahrd.f"
	i__3 = i__ - 1;
#line 268 "slahrd.f"
	sgemv_("Transpose", &i__2, &i__3, &c_b5, &a[*k + i__ + a_dim1], lda, &
		a[*k + i__ + i__ * a_dim1], &c__1, &c_b38, &t[i__ * t_dim1 + 
		1], &c__1, (ftnlen)9);
#line 270 "slahrd.f"
	i__2 = i__ - 1;
#line 270 "slahrd.f"
	sgemv_("No transpose", n, &i__2, &c_b4, &y[y_offset], ldy, &t[i__ * 
		t_dim1 + 1], &c__1, &c_b5, &y[i__ * y_dim1 + 1], &c__1, (
		ftnlen)12);
#line 272 "slahrd.f"
	sscal_(n, &tau[i__], &y[i__ * y_dim1 + 1], &c__1);

/*        Compute T(1:i,i) */

#line 276 "slahrd.f"
	i__2 = i__ - 1;
#line 276 "slahrd.f"
	d__1 = -tau[i__];
#line 276 "slahrd.f"
	sscal_(&i__2, &d__1, &t[i__ * t_dim1 + 1], &c__1);
#line 277 "slahrd.f"
	i__2 = i__ - 1;
#line 277 "slahrd.f"
	strmv_("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt, 
		&t[i__ * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8)
		;
#line 279 "slahrd.f"
	t[i__ + i__ * t_dim1] = tau[i__];

#line 281 "slahrd.f"
/* L10: */
#line 281 "slahrd.f"
    }
#line 282 "slahrd.f"
    a[*k + *nb + *nb * a_dim1] = ei;

#line 284 "slahrd.f"
    return 0;

/*     End of SLAHRD */

} /* slahrd_ */

