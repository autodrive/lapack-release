#line 1 "clahr2.f"
/* clahr2.f -- translated by f2c (version 20100827).
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

#line 1 "clahr2.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that 
elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to 
apply the transformation to the unreduced part */
/* of A. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAHR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY ) */

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
/* > CLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1) */
/* > matrix A so that elements below the k-th subdiagonal are zero. The */
/* > reduction is performed by an unitary similarity transformation */
/* > Q**H * A * Q. The routine returns the matrices V and T which determine */
/* > Q as a block reflector I - V*T*v**H, and also the matrix Y = A * V * T. */
/* > */
/* > This is an auxiliary routine called by CGEHRD. */
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
/* >          K < N. */
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
/* >          The leading dimension of the array Y. LDY >= N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

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
/* >     ( a   a   a   a   a ) */
/* >     ( a   a   a   a   a ) */
/* >     ( a   a   a   a   a ) */
/* >     ( h   h   a   a   a ) */
/* >     ( v1  h   a   a   a ) */
/* >     ( v1  v2  a   a   a ) */
/* >     ( v1  v2  a   a   a ) */
/* > */
/* >  where a denotes an element of the original matrix A, h denotes a */
/* >  modified element of the upper Hessenberg matrix H, and vi denotes an */
/* >  element of the vector defining H(i). */
/* > */
/* >  This subroutine is a slight modification of LAPACK-3.0's DLAHRD */
/* >  incorporating improvements proposed by Quintana-Orti and Van de */
/* >  Gejin. Note that the entries of A(1:K,2:NB) differ from those */
/* >  returned by the original LAPACK-3.0's DLAHRD routine. (This */
/* >  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.) */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the */
/* >  performance of reduction to Hessenberg form," ACM Transactions on */
/* >  Mathematical Software, 32(2):180-194, June 2006. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clahr2_(integer *n, integer *k, integer *nb, 
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
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), ctrmm_(char *, char *, char *, char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    caxpy_(integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ctrmv_(char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), clarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), clacgv_(integer *, 
	    doublecomplex *, integer *), clacpy_(char *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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

/*     Quick return if possible */

#line 219 "clahr2.f"
    /* Parameter adjustments */
#line 219 "clahr2.f"
    --tau;
#line 219 "clahr2.f"
    a_dim1 = *lda;
#line 219 "clahr2.f"
    a_offset = 1 + a_dim1;
#line 219 "clahr2.f"
    a -= a_offset;
#line 219 "clahr2.f"
    t_dim1 = *ldt;
#line 219 "clahr2.f"
    t_offset = 1 + t_dim1;
#line 219 "clahr2.f"
    t -= t_offset;
#line 219 "clahr2.f"
    y_dim1 = *ldy;
#line 219 "clahr2.f"
    y_offset = 1 + y_dim1;
#line 219 "clahr2.f"
    y -= y_offset;
#line 219 "clahr2.f"

#line 219 "clahr2.f"
    /* Function Body */
#line 219 "clahr2.f"
    if (*n <= 1) {
#line 219 "clahr2.f"
	return 0;
#line 219 "clahr2.f"
    }

#line 222 "clahr2.f"
    i__1 = *nb;
#line 222 "clahr2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 223 "clahr2.f"
	if (i__ > 1) {

/*           Update A(K+1:N,I) */

/*           Update I-th column of A - Y * V**H */

#line 229 "clahr2.f"
	    i__2 = i__ - 1;
#line 229 "clahr2.f"
	    clacgv_(&i__2, &a[*k + i__ - 1 + a_dim1], lda);
#line 230 "clahr2.f"
	    i__2 = *n - *k;
#line 230 "clahr2.f"
	    i__3 = i__ - 1;
#line 230 "clahr2.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "clahr2.f"
	    cgemv_("NO TRANSPOSE", &i__2, &i__3, &z__1, &y[*k + 1 + y_dim1], 
		    ldy, &a[*k + i__ - 1 + a_dim1], lda, &c_b2, &a[*k + 1 + 
		    i__ * a_dim1], &c__1, (ftnlen)12);
#line 232 "clahr2.f"
	    i__2 = i__ - 1;
#line 232 "clahr2.f"
	    clacgv_(&i__2, &a[*k + i__ - 1 + a_dim1], lda);

/*           Apply I - V * T**H * V**H to this column (call it b) from the */
/*           left, using the last column of T as workspace */

/*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows) */
/*                    ( V2 )             ( b2 ) */

/*           where V1 is unit lower triangular */

/*           w := V1**H * b1 */

#line 244 "clahr2.f"
	    i__2 = i__ - 1;
#line 244 "clahr2.f"
	    ccopy_(&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 
		    1], &c__1);
#line 245 "clahr2.f"
	    i__2 = i__ - 1;
#line 245 "clahr2.f"
	    ctrmv_("Lower", "Conjugate transpose", "UNIT", &i__2, &a[*k + 1 + 
		    a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (
		    ftnlen)19, (ftnlen)4);

/*           w := w + V2**H * b2 */

#line 251 "clahr2.f"
	    i__2 = *n - *k - i__ + 1;
#line 251 "clahr2.f"
	    i__3 = i__ - 1;
#line 251 "clahr2.f"
	    cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i__ + 
		    a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b2, &
		    t[*nb * t_dim1 + 1], &c__1, (ftnlen)19);

/*           w := T**H * w */

#line 257 "clahr2.f"
	    i__2 = i__ - 1;
#line 257 "clahr2.f"
	    ctrmv_("Upper", "Conjugate transpose", "NON-UNIT", &i__2, &t[
		    t_offset], ldt, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (
		    ftnlen)19, (ftnlen)8);

/*           b2 := b2 - V2*w */

#line 263 "clahr2.f"
	    i__2 = *n - *k - i__ + 1;
#line 263 "clahr2.f"
	    i__3 = i__ - 1;
#line 263 "clahr2.f"
	    z__1.r = -1., z__1.i = -0.;
#line 263 "clahr2.f"
	    cgemv_("NO TRANSPOSE", &i__2, &i__3, &z__1, &a[*k + i__ + a_dim1],
		     lda, &t[*nb * t_dim1 + 1], &c__1, &c_b2, &a[*k + i__ + 
		    i__ * a_dim1], &c__1, (ftnlen)12);

/*           b1 := b1 - V1*w */

#line 269 "clahr2.f"
	    i__2 = i__ - 1;
#line 269 "clahr2.f"
	    ctrmv_("Lower", "NO TRANSPOSE", "UNIT", &i__2, &a[*k + 1 + a_dim1]
		    , lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12,
		     (ftnlen)4);
#line 272 "clahr2.f"
	    i__2 = i__ - 1;
#line 272 "clahr2.f"
	    z__1.r = -1., z__1.i = -0.;
#line 272 "clahr2.f"
	    caxpy_(&i__2, &z__1, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ 
		    * a_dim1], &c__1);

#line 274 "clahr2.f"
	    i__2 = *k + i__ - 1 + (i__ - 1) * a_dim1;
#line 274 "clahr2.f"
	    a[i__2].r = ei.r, a[i__2].i = ei.i;
#line 275 "clahr2.f"
	}

/*        Generate the elementary reflector H(I) to annihilate */
/*        A(K+I+1:N,I) */

#line 280 "clahr2.f"
	i__2 = *n - *k - i__ + 1;
/* Computing MIN */
#line 280 "clahr2.f"
	i__3 = *k + i__ + 1;
#line 280 "clahr2.f"
	clarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[min(i__3,*n) + i__ * 
		a_dim1], &c__1, &tau[i__]);
#line 282 "clahr2.f"
	i__2 = *k + i__ + i__ * a_dim1;
#line 282 "clahr2.f"
	ei.r = a[i__2].r, ei.i = a[i__2].i;
#line 283 "clahr2.f"
	i__2 = *k + i__ + i__ * a_dim1;
#line 283 "clahr2.f"
	a[i__2].r = 1., a[i__2].i = 0.;

/*        Compute  Y(K+1:N,I) */

#line 287 "clahr2.f"
	i__2 = *n - *k;
#line 287 "clahr2.f"
	i__3 = *n - *k - i__ + 1;
#line 287 "clahr2.f"
	cgemv_("NO TRANSPOSE", &i__2, &i__3, &c_b2, &a[*k + 1 + (i__ + 1) * 
		a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b1, &y[*
		k + 1 + i__ * y_dim1], &c__1, (ftnlen)12);
#line 290 "clahr2.f"
	i__2 = *n - *k - i__ + 1;
#line 290 "clahr2.f"
	i__3 = i__ - 1;
#line 290 "clahr2.f"
	cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i__ + 
		a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b1, &t[
		i__ * t_dim1 + 1], &c__1, (ftnlen)19);
#line 293 "clahr2.f"
	i__2 = *n - *k;
#line 293 "clahr2.f"
	i__3 = i__ - 1;
#line 293 "clahr2.f"
	z__1.r = -1., z__1.i = -0.;
#line 293 "clahr2.f"
	cgemv_("NO TRANSPOSE", &i__2, &i__3, &z__1, &y[*k + 1 + y_dim1], ldy, 
		&t[i__ * t_dim1 + 1], &c__1, &c_b2, &y[*k + 1 + i__ * y_dim1],
		 &c__1, (ftnlen)12);
#line 296 "clahr2.f"
	i__2 = *n - *k;
#line 296 "clahr2.f"
	cscal_(&i__2, &tau[i__], &y[*k + 1 + i__ * y_dim1], &c__1);

/*        Compute T(1:I,I) */

#line 300 "clahr2.f"
	i__2 = i__ - 1;
#line 300 "clahr2.f"
	i__3 = i__;
#line 300 "clahr2.f"
	z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
#line 300 "clahr2.f"
	cscal_(&i__2, &z__1, &t[i__ * t_dim1 + 1], &c__1);
#line 301 "clahr2.f"
	i__2 = i__ - 1;
#line 301 "clahr2.f"
	ctrmv_("Upper", "No Transpose", "NON-UNIT", &i__2, &t[t_offset], ldt, 
		&t[i__ * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8)
		;
#line 304 "clahr2.f"
	i__2 = i__ + i__ * t_dim1;
#line 304 "clahr2.f"
	i__3 = i__;
#line 304 "clahr2.f"
	t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;

#line 306 "clahr2.f"
/* L10: */
#line 306 "clahr2.f"
    }
#line 307 "clahr2.f"
    i__1 = *k + *nb + *nb * a_dim1;
#line 307 "clahr2.f"
    a[i__1].r = ei.r, a[i__1].i = ei.i;

/*     Compute Y(1:K,1:NB) */

#line 311 "clahr2.f"
    clacpy_("ALL", k, nb, &a[(a_dim1 << 1) + 1], lda, &y[y_offset], ldy, (
	    ftnlen)3);
#line 312 "clahr2.f"
    ctrmm_("RIGHT", "Lower", "NO TRANSPOSE", "UNIT", k, nb, &c_b2, &a[*k + 1 
	    + a_dim1], lda, &y[y_offset], ldy, (ftnlen)5, (ftnlen)5, (ftnlen)
	    12, (ftnlen)4);
#line 315 "clahr2.f"
    if (*n > *k + *nb) {
#line 315 "clahr2.f"
	i__1 = *n - *k - *nb;
#line 315 "clahr2.f"
	cgemm_("NO TRANSPOSE", "NO TRANSPOSE", k, nb, &i__1, &c_b2, &a[(*nb + 
		2) * a_dim1 + 1], lda, &a[*k + 1 + *nb + a_dim1], lda, &c_b2, 
		&y[y_offset], ldy, (ftnlen)12, (ftnlen)12);
#line 315 "clahr2.f"
    }
#line 320 "clahr2.f"
    ctrmm_("RIGHT", "Upper", "NO TRANSPOSE", "NON-UNIT", k, nb, &c_b2, &t[
	    t_offset], ldt, &y[y_offset], ldy, (ftnlen)5, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);

#line 324 "clahr2.f"
    return 0;

/*     End of CLAHR2 */

} /* clahr2_ */

