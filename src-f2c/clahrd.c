#line 1 "clahrd.f"
/* clahrd.f -- translated by f2c (version 20100827).
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

#line 1 "clahrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLAHRD reduces the first nb columns of a general rectangular matrix A so that elements below th
e k-th subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformati
on to the unreduced part of A. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            K, LDA, LDT, LDY, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), T( LDT, NB ), TAU( NB ), */
/*      $                   Y( LDY, NB ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAHRD reduces the first NB columns of a complex general n-by-(n-k+1) */
/* > matrix A so that elements below the k-th subdiagonal are zero. The */
/* > reduction is performed by a unitary similarity transformation */
/* > Q**H * A * Q. The routine returns the matrices V and T which determine */
/* > Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T. */
/* > */
/* > This is an OBSOLETE auxiliary routine. */
/* > This routine will be 'deprecated' in a  future release. */
/* > Please use the new routine CLAHR2 instead. */
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
/* >          A is COMPLEX array, dimension (LDA,N-K+1) */
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
/* >          TAU is COMPLEX array, dimension (NB) */
/* >          The scalar factors of the elementary reflectors. See Further */
/* >          Details. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,NB) */
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
/* >          Y is COMPLEX array, dimension (LDY,NB) */
/* >          The n-by-nb matrix Y. */
/* > \endverbatim */
/* > */
/* > \param[in] LDY */
/* > \verbatim */
/* >          LDY is INTEGER */
/* >          The leading dimension of the array Y. LDY >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

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
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in */
/* >  A(i+k+1:n,i), and tau in TAU(i). */
/* > */
/* >  The elements of the vectors v together form the (n-k+1)-by-nb matrix */
/* >  V which is needed, with T and Y, to apply the transformation to the */
/* >  unreduced part of the matrix, using an update of the form: */
/* >  A := (I - V*T*V**H) * (A - Y*V**H). */
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
/* Subroutine */ int clahrd_(integer *n, integer *k, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t, 
	integer *ldt, doublecomplex *y, integer *ldy)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__;
    static doublecomplex ei;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), caxpy_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ctrmv_(char *, char *, 
	    char *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen), clarfg_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *), 
	    clacgv_(integer *, doublecomplex *, integer *);


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

#line 207 "clahrd.f"
    /* Parameter adjustments */
#line 207 "clahrd.f"
    --tau;
#line 207 "clahrd.f"
    a_dim1 = *lda;
#line 207 "clahrd.f"
    a_offset = 1 + a_dim1;
#line 207 "clahrd.f"
    a -= a_offset;
#line 207 "clahrd.f"
    t_dim1 = *ldt;
#line 207 "clahrd.f"
    t_offset = 1 + t_dim1;
#line 207 "clahrd.f"
    t -= t_offset;
#line 207 "clahrd.f"
    y_dim1 = *ldy;
#line 207 "clahrd.f"
    y_offset = 1 + y_dim1;
#line 207 "clahrd.f"
    y -= y_offset;
#line 207 "clahrd.f"

#line 207 "clahrd.f"
    /* Function Body */
#line 207 "clahrd.f"
    if (*n <= 1) {
#line 207 "clahrd.f"
	return 0;
#line 207 "clahrd.f"
    }

#line 210 "clahrd.f"
    i__1 = *nb;
#line 210 "clahrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "clahrd.f"
	if (i__ > 1) {

/*           Update A(1:n,i) */

/*           Compute i-th column of A - Y * V**H */

#line 217 "clahrd.f"
	    i__2 = i__ - 1;
#line 217 "clahrd.f"
	    clacgv_(&i__2, &a[*k + i__ - 1 + a_dim1], lda);
#line 218 "clahrd.f"
	    i__2 = i__ - 1;
#line 218 "clahrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 218 "clahrd.f"
	    cgemv_("No transpose", n, &i__2, &z__1, &y[y_offset], ldy, &a[*k 
		    + i__ - 1 + a_dim1], lda, &c_b2, &a[i__ * a_dim1 + 1], &
		    c__1, (ftnlen)12);
#line 220 "clahrd.f"
	    i__2 = i__ - 1;
#line 220 "clahrd.f"
	    clacgv_(&i__2, &a[*k + i__ - 1 + a_dim1], lda);

/*           Apply I - V * T**H * V**H to this column (call it b) from the */
/*           left, using the last column of T as workspace */

/*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows) */
/*                    ( V2 )             ( b2 ) */

/*           where V1 is unit lower triangular */

/*           w := V1**H * b1 */

#line 232 "clahrd.f"
	    i__2 = i__ - 1;
#line 232 "clahrd.f"
	    ccopy_(&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 
		    1], &c__1);
#line 233 "clahrd.f"
	    i__2 = i__ - 1;
#line 233 "clahrd.f"
	    ctrmv_("Lower", "Conjugate transpose", "Unit", &i__2, &a[*k + 1 + 
		    a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (
		    ftnlen)19, (ftnlen)4);

/*           w := w + V2**H *b2 */

#line 238 "clahrd.f"
	    i__2 = *n - *k - i__ + 1;
#line 238 "clahrd.f"
	    i__3 = i__ - 1;
#line 238 "clahrd.f"
	    cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i__ + 
		    a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b2, &
		    t[*nb * t_dim1 + 1], &c__1, (ftnlen)19);

/*           w := T**H *w */

#line 244 "clahrd.f"
	    i__2 = i__ - 1;
#line 244 "clahrd.f"
	    ctrmv_("Upper", "Conjugate transpose", "Non-unit", &i__2, &t[
		    t_offset], ldt, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

/*           b2 := b2 - V2*w */

#line 249 "clahrd.f"
	    i__2 = *n - *k - i__ + 1;
#line 249 "clahrd.f"
	    i__3 = i__ - 1;
#line 249 "clahrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 249 "clahrd.f"
	    cgemv_("No transpose", &i__2, &i__3, &z__1, &a[*k + i__ + a_dim1],
		     lda, &t[*nb * t_dim1 + 1], &c__1, &c_b2, &a[*k + i__ + 
		    i__ * a_dim1], &c__1, (ftnlen)12);

/*           b1 := b1 - V1*w */

#line 254 "clahrd.f"
	    i__2 = i__ - 1;
#line 254 "clahrd.f"
	    ctrmv_("Lower", "No transpose", "Unit", &i__2, &a[*k + 1 + a_dim1]
		    , lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12,
		     (ftnlen)4);
#line 256 "clahrd.f"
	    i__2 = i__ - 1;
#line 256 "clahrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 256 "clahrd.f"
	    caxpy_(&i__2, &z__1, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ 
		    * a_dim1], &c__1);

#line 258 "clahrd.f"
	    i__2 = *k + i__ - 1 + (i__ - 1) * a_dim1;
#line 258 "clahrd.f"
	    a[i__2].r = ei.r, a[i__2].i = ei.i;
#line 259 "clahrd.f"
	}

/*        Generate the elementary reflector H(i) to annihilate */
/*        A(k+i+1:n,i) */

#line 264 "clahrd.f"
	i__2 = *k + i__ + i__ * a_dim1;
#line 264 "clahrd.f"
	ei.r = a[i__2].r, ei.i = a[i__2].i;
#line 265 "clahrd.f"
	i__2 = *n - *k - i__ + 1;
/* Computing MIN */
#line 265 "clahrd.f"
	i__3 = *k + i__ + 1;
#line 265 "clahrd.f"
	clarfg_(&i__2, &ei, &a[min(i__3,*n) + i__ * a_dim1], &c__1, &tau[i__])
		;
#line 267 "clahrd.f"
	i__2 = *k + i__ + i__ * a_dim1;
#line 267 "clahrd.f"
	a[i__2].r = 1., a[i__2].i = 0.;

/*        Compute  Y(1:n,i) */

#line 271 "clahrd.f"
	i__2 = *n - *k - i__ + 1;
#line 271 "clahrd.f"
	cgemv_("No transpose", n, &i__2, &c_b2, &a[(i__ + 1) * a_dim1 + 1], 
		lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b1, &y[i__ * 
		y_dim1 + 1], &c__1, (ftnlen)12);
#line 273 "clahrd.f"
	i__2 = *n - *k - i__ + 1;
#line 273 "clahrd.f"
	i__3 = i__ - 1;
#line 273 "clahrd.f"
	cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i__ + 
		a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b1, &t[
		i__ * t_dim1 + 1], &c__1, (ftnlen)19);
#line 276 "clahrd.f"
	i__2 = i__ - 1;
#line 276 "clahrd.f"
	z__1.r = -1., z__1.i = -0.;
#line 276 "clahrd.f"
	cgemv_("No transpose", n, &i__2, &z__1, &y[y_offset], ldy, &t[i__ * 
		t_dim1 + 1], &c__1, &c_b2, &y[i__ * y_dim1 + 1], &c__1, (
		ftnlen)12);
#line 278 "clahrd.f"
	cscal_(n, &tau[i__], &y[i__ * y_dim1 + 1], &c__1);

/*        Compute T(1:i,i) */

#line 282 "clahrd.f"
	i__2 = i__ - 1;
#line 282 "clahrd.f"
	i__3 = i__;
#line 282 "clahrd.f"
	z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
#line 282 "clahrd.f"
	cscal_(&i__2, &z__1, &t[i__ * t_dim1 + 1], &c__1);
#line 283 "clahrd.f"
	i__2 = i__ - 1;
#line 283 "clahrd.f"
	ctrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt, 
		&t[i__ * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8)
		;
#line 285 "clahrd.f"
	i__2 = i__ + i__ * t_dim1;
#line 285 "clahrd.f"
	i__3 = i__;
#line 285 "clahrd.f"
	t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;

#line 287 "clahrd.f"
/* L10: */
#line 287 "clahrd.f"
    }
#line 288 "clahrd.f"
    i__1 = *k + *nb + *nb * a_dim1;
#line 288 "clahrd.f"
    a[i__1].r = ei.r, a[i__1].i = ei.i;

#line 290 "clahrd.f"
    return 0;

/*     End of CLAHRD */

} /* clahrd_ */

