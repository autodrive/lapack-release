#line 1 "clabrd.f"
/* clabrd.f -- translated by f2c (version 20100827).
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

#line 1 "clabrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLABRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clabrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clabrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clabrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, */
/*                          LDY ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, LDX, LDY, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ) */
/*       COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ), */
/*      $                   Y( LDY, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLABRD reduces the first NB rows and columns of a complex general */
/* > m by n matrix A to upper or lower real bidiagonal form by a unitary */
/* > transformation Q**H * A * P, and returns the matrices X and Y which */
/* > are needed to apply the transformation to the unreduced part of A. */
/* > */
/* > If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower */
/* > bidiagonal form. */
/* > */
/* > This is an auxiliary routine called by CGEBRD */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows in the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns in the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The number of leading rows and columns of A to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the m by n general matrix to be reduced. */
/* >          On exit, the first NB rows and columns of the matrix are */
/* >          overwritten; the rest of the array is unchanged. */
/* >          If m >= n, elements on and below the diagonal in the first NB */
/* >            columns, with the array TAUQ, represent the unitary */
/* >            matrix Q as a product of elementary reflectors; and */
/* >            elements above the diagonal in the first NB rows, with the */
/* >            array TAUP, represent the unitary matrix P as a product */
/* >            of elementary reflectors. */
/* >          If m < n, elements below the diagonal in the first NB */
/* >            columns, with the array TAUQ, represent the unitary */
/* >            matrix Q as a product of elementary reflectors, and */
/* >            elements on and above the diagonal in the first NB rows, */
/* >            with the array TAUP, represent the unitary matrix P as */
/* >            a product of elementary reflectors. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (NB) */
/* >          The diagonal elements of the first NB rows and columns of */
/* >          the reduced matrix.  D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (NB) */
/* >          The off-diagonal elements of the first NB rows and columns of */
/* >          the reduced matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* >          TAUQ is COMPLEX array dimension (NB) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the unitary matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* >          TAUP is COMPLEX array, dimension (NB) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the unitary matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (LDX,NB) */
/* >          The m-by-nb matrix X required to update the unreduced part */
/* >          of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X. LDX >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension (LDY,NB) */
/* >          The n-by-nb matrix Y required to update the unreduced part */
/* >          of A. */
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
/* >  The matrices Q and P are represented as products of elementary */
/* >  reflectors: */
/* > */
/* >     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb) */
/* > */
/* >  Each H(i) and G(i) has the form: */
/* > */
/* >     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H */
/* > */
/* >  where tauq and taup are complex scalars, and v and u are complex */
/* >  vectors. */
/* > */
/* >  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in */
/* >  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in */
/* >  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in */
/* >  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in */
/* >  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  The elements of the vectors v and u together form the m-by-nb matrix */
/* >  V and the nb-by-n matrix U**H which are needed, with X and Y, to apply */
/* >  the transformation to the unreduced part of the matrix, using a block */
/* >  update of the form:  A := A - V*Y**H - X*U**H. */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with nb = 2: */
/* > */
/* >  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n): */
/* > */
/* >    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 ) */
/* >    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 ) */
/* >    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  ) */
/* >    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  ) */
/* >    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  ) */
/* >    (  v1  v2  a   a   a  ) */
/* > */
/* >  where a denotes an element of the original matrix which is unchanged, */
/* >  vi denotes an element of the vector defining H(i), and ui an element */
/* >  of the vector defining G(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clabrd_(integer *m, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, 
	doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, integer *
	ldx, doublecomplex *y, integer *ldy)
{
    /* System generated locals */
    integer a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__;
    static doublecomplex alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    clarfg_(integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *), clacgv_(integer *, doublecomplex *, integer *);


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

#line 250 "clabrd.f"
    /* Parameter adjustments */
#line 250 "clabrd.f"
    a_dim1 = *lda;
#line 250 "clabrd.f"
    a_offset = 1 + a_dim1;
#line 250 "clabrd.f"
    a -= a_offset;
#line 250 "clabrd.f"
    --d__;
#line 250 "clabrd.f"
    --e;
#line 250 "clabrd.f"
    --tauq;
#line 250 "clabrd.f"
    --taup;
#line 250 "clabrd.f"
    x_dim1 = *ldx;
#line 250 "clabrd.f"
    x_offset = 1 + x_dim1;
#line 250 "clabrd.f"
    x -= x_offset;
#line 250 "clabrd.f"
    y_dim1 = *ldy;
#line 250 "clabrd.f"
    y_offset = 1 + y_dim1;
#line 250 "clabrd.f"
    y -= y_offset;
#line 250 "clabrd.f"

#line 250 "clabrd.f"
    /* Function Body */
#line 250 "clabrd.f"
    if (*m <= 0 || *n <= 0) {
#line 250 "clabrd.f"
	return 0;
#line 250 "clabrd.f"
    }

#line 253 "clabrd.f"
    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

#line 257 "clabrd.f"
	i__1 = *nb;
#line 257 "clabrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i:m,i) */

#line 261 "clabrd.f"
	    i__2 = i__ - 1;
#line 261 "clabrd.f"
	    clacgv_(&i__2, &y[i__ + y_dim1], ldy);
#line 262 "clabrd.f"
	    i__2 = *m - i__ + 1;
#line 262 "clabrd.f"
	    i__3 = i__ - 1;
#line 262 "clabrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 262 "clabrd.f"
	    cgemv_("No transpose", &i__2, &i__3, &z__1, &a[i__ + a_dim1], lda,
		     &y[i__ + y_dim1], ldy, &c_b2, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 264 "clabrd.f"
	    i__2 = i__ - 1;
#line 264 "clabrd.f"
	    clacgv_(&i__2, &y[i__ + y_dim1], ldy);
#line 265 "clabrd.f"
	    i__2 = *m - i__ + 1;
#line 265 "clabrd.f"
	    i__3 = i__ - 1;
#line 265 "clabrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 265 "clabrd.f"
	    cgemv_("No transpose", &i__2, &i__3, &z__1, &x[i__ + x_dim1], ldx,
		     &a[i__ * a_dim1 + 1], &c__1, &c_b2, &a[i__ + i__ * 
		    a_dim1], &c__1, (ftnlen)12);

/*           Generate reflection Q(i) to annihilate A(i+1:m,i) */

#line 270 "clabrd.f"
	    i__2 = i__ + i__ * a_dim1;
#line 270 "clabrd.f"
	    alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 271 "clabrd.f"
	    i__2 = *m - i__ + 1;
/* Computing MIN */
#line 271 "clabrd.f"
	    i__3 = i__ + 1;
#line 271 "clabrd.f"
	    clarfg_(&i__2, &alpha, &a[min(i__3,*m) + i__ * a_dim1], &c__1, &
		    tauq[i__]);
#line 273 "clabrd.f"
	    i__2 = i__;
#line 273 "clabrd.f"
	    d__[i__2] = alpha.r;
#line 274 "clabrd.f"
	    if (i__ < *n) {
#line 275 "clabrd.f"
		i__2 = i__ + i__ * a_dim1;
#line 275 "clabrd.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute Y(i+1:n,i) */

#line 279 "clabrd.f"
		i__2 = *m - i__ + 1;
#line 279 "clabrd.f"
		i__3 = *n - i__;
#line 279 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[i__ + (
			i__ + 1) * a_dim1], lda, &a[i__ + i__ * a_dim1], &
			c__1, &c_b1, &y[i__ + 1 + i__ * y_dim1], &c__1, (
			ftnlen)19);
#line 282 "clabrd.f"
		i__2 = *m - i__ + 1;
#line 282 "clabrd.f"
		i__3 = i__ - 1;
#line 282 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[i__ + 
			a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &c_b1, &
			y[i__ * y_dim1 + 1], &c__1, (ftnlen)19);
#line 285 "clabrd.f"
		i__2 = *n - i__;
#line 285 "clabrd.f"
		i__3 = i__ - 1;
#line 285 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 285 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &c_b2, &y[
			i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)12);
#line 287 "clabrd.f"
		i__2 = *m - i__ + 1;
#line 287 "clabrd.f"
		i__3 = i__ - 1;
#line 287 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &x[i__ + 
			x_dim1], ldx, &a[i__ + i__ * a_dim1], &c__1, &c_b1, &
			y[i__ * y_dim1 + 1], &c__1, (ftnlen)19);
#line 290 "clabrd.f"
		i__2 = i__ - 1;
#line 290 "clabrd.f"
		i__3 = *n - i__;
#line 290 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 290 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &z__1, &a[(i__ + 
			1) * a_dim1 + 1], lda, &y[i__ * y_dim1 + 1], &c__1, &
			c_b2, &y[i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)19);
#line 293 "clabrd.f"
		i__2 = *n - i__;
#line 293 "clabrd.f"
		cscal_(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);

/*              Update A(i,i+1:n) */

#line 297 "clabrd.f"
		i__2 = *n - i__;
#line 297 "clabrd.f"
		clacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 298 "clabrd.f"
		clacgv_(&i__, &a[i__ + a_dim1], lda);
#line 299 "clabrd.f"
		i__2 = *n - i__;
#line 299 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 299 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__, &z__1, &y[i__ + 1 + 
			y_dim1], ldy, &a[i__ + a_dim1], lda, &c_b2, &a[i__ + (
			i__ + 1) * a_dim1], lda, (ftnlen)12);
#line 301 "clabrd.f"
		clacgv_(&i__, &a[i__ + a_dim1], lda);
#line 302 "clabrd.f"
		i__2 = i__ - 1;
#line 302 "clabrd.f"
		clacgv_(&i__2, &x[i__ + x_dim1], ldx);
#line 303 "clabrd.f"
		i__2 = i__ - 1;
#line 303 "clabrd.f"
		i__3 = *n - i__;
#line 303 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 303 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &z__1, &a[(i__ + 
			1) * a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &c_b2, &
			a[i__ + (i__ + 1) * a_dim1], lda, (ftnlen)19);
#line 306 "clabrd.f"
		i__2 = i__ - 1;
#line 306 "clabrd.f"
		clacgv_(&i__2, &x[i__ + x_dim1], ldx);

/*              Generate reflection P(i) to annihilate A(i,i+2:n) */

#line 310 "clabrd.f"
		i__2 = i__ + (i__ + 1) * a_dim1;
#line 310 "clabrd.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 311 "clabrd.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 311 "clabrd.f"
		i__3 = i__ + 2;
#line 311 "clabrd.f"
		clarfg_(&i__2, &alpha, &a[i__ + min(i__3,*n) * a_dim1], lda, &
			taup[i__]);
#line 313 "clabrd.f"
		i__2 = i__;
#line 313 "clabrd.f"
		e[i__2] = alpha.r;
#line 314 "clabrd.f"
		i__2 = i__ + (i__ + 1) * a_dim1;
#line 314 "clabrd.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute X(i+1:m,i) */

#line 318 "clabrd.f"
		i__2 = *m - i__;
#line 318 "clabrd.f"
		i__3 = *n - i__;
#line 318 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 + (i__ 
			+ 1) * a_dim1], lda, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &c_b1, &x[i__ + 1 + i__ * x_dim1], &c__1, (
			ftnlen)12);
#line 320 "clabrd.f"
		i__2 = *n - i__;
#line 320 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__, &c_b2, &y[i__ + 1 
			+ y_dim1], ldy, &a[i__ + (i__ + 1) * a_dim1], lda, &
			c_b1, &x[i__ * x_dim1 + 1], &c__1, (ftnlen)19);
#line 323 "clabrd.f"
		i__2 = *m - i__;
#line 323 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 323 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__, &z__1, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &c_b2, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 325 "clabrd.f"
		i__2 = i__ - 1;
#line 325 "clabrd.f"
		i__3 = *n - i__;
#line 325 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &a[i__ + (i__ + 1) * a_dim1], lda, &
			c_b1, &x[i__ * x_dim1 + 1], &c__1, (ftnlen)12);
#line 327 "clabrd.f"
		i__2 = *m - i__;
#line 327 "clabrd.f"
		i__3 = i__ - 1;
#line 327 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 327 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &c_b2, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 329 "clabrd.f"
		i__2 = *m - i__;
#line 329 "clabrd.f"
		cscal_(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
#line 330 "clabrd.f"
		i__2 = *n - i__;
#line 330 "clabrd.f"
		clacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 331 "clabrd.f"
	    }
#line 332 "clabrd.f"
/* L10: */
#line 332 "clabrd.f"
	}
#line 333 "clabrd.f"
    } else {

/*        Reduce to lower bidiagonal form */

#line 337 "clabrd.f"
	i__1 = *nb;
#line 337 "clabrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i,i:n) */

#line 341 "clabrd.f"
	    i__2 = *n - i__ + 1;
#line 341 "clabrd.f"
	    clacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);
#line 342 "clabrd.f"
	    i__2 = i__ - 1;
#line 342 "clabrd.f"
	    clacgv_(&i__2, &a[i__ + a_dim1], lda);
#line 343 "clabrd.f"
	    i__2 = *n - i__ + 1;
#line 343 "clabrd.f"
	    i__3 = i__ - 1;
#line 343 "clabrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 343 "clabrd.f"
	    cgemv_("No transpose", &i__2, &i__3, &z__1, &y[i__ + y_dim1], ldy,
		     &a[i__ + a_dim1], lda, &c_b2, &a[i__ + i__ * a_dim1], 
		    lda, (ftnlen)12);
#line 345 "clabrd.f"
	    i__2 = i__ - 1;
#line 345 "clabrd.f"
	    clacgv_(&i__2, &a[i__ + a_dim1], lda);
#line 346 "clabrd.f"
	    i__2 = i__ - 1;
#line 346 "clabrd.f"
	    clacgv_(&i__2, &x[i__ + x_dim1], ldx);
#line 347 "clabrd.f"
	    i__2 = i__ - 1;
#line 347 "clabrd.f"
	    i__3 = *n - i__ + 1;
#line 347 "clabrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 347 "clabrd.f"
	    cgemv_("Conjugate transpose", &i__2, &i__3, &z__1, &a[i__ * 
		    a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &c_b2, &a[i__ + 
		    i__ * a_dim1], lda, (ftnlen)19);
#line 350 "clabrd.f"
	    i__2 = i__ - 1;
#line 350 "clabrd.f"
	    clacgv_(&i__2, &x[i__ + x_dim1], ldx);

/*           Generate reflection P(i) to annihilate A(i,i+1:n) */

#line 354 "clabrd.f"
	    i__2 = i__ + i__ * a_dim1;
#line 354 "clabrd.f"
	    alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 355 "clabrd.f"
	    i__2 = *n - i__ + 1;
/* Computing MIN */
#line 355 "clabrd.f"
	    i__3 = i__ + 1;
#line 355 "clabrd.f"
	    clarfg_(&i__2, &alpha, &a[i__ + min(i__3,*n) * a_dim1], lda, &
		    taup[i__]);
#line 357 "clabrd.f"
	    i__2 = i__;
#line 357 "clabrd.f"
	    d__[i__2] = alpha.r;
#line 358 "clabrd.f"
	    if (i__ < *m) {
#line 359 "clabrd.f"
		i__2 = i__ + i__ * a_dim1;
#line 359 "clabrd.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute X(i+1:m,i) */

#line 363 "clabrd.f"
		i__2 = *m - i__;
#line 363 "clabrd.f"
		i__3 = *n - i__ + 1;
#line 363 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 + i__ *
			 a_dim1], lda, &a[i__ + i__ * a_dim1], lda, &c_b1, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 365 "clabrd.f"
		i__2 = *n - i__ + 1;
#line 365 "clabrd.f"
		i__3 = i__ - 1;
#line 365 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &y[i__ + 
			y_dim1], ldy, &a[i__ + i__ * a_dim1], lda, &c_b1, &x[
			i__ * x_dim1 + 1], &c__1, (ftnlen)19);
#line 368 "clabrd.f"
		i__2 = *m - i__;
#line 368 "clabrd.f"
		i__3 = i__ - 1;
#line 368 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 368 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &c_b2, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 370 "clabrd.f"
		i__2 = i__ - 1;
#line 370 "clabrd.f"
		i__3 = *n - i__ + 1;
#line 370 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &a[i__ * a_dim1 + 
			1], lda, &a[i__ + i__ * a_dim1], lda, &c_b1, &x[i__ * 
			x_dim1 + 1], &c__1, (ftnlen)12);
#line 372 "clabrd.f"
		i__2 = *m - i__;
#line 372 "clabrd.f"
		i__3 = i__ - 1;
#line 372 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 372 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &c_b2, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 374 "clabrd.f"
		i__2 = *m - i__;
#line 374 "clabrd.f"
		cscal_(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
#line 375 "clabrd.f"
		i__2 = *n - i__ + 1;
#line 375 "clabrd.f"
		clacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);

/*              Update A(i+1:m,i) */

#line 379 "clabrd.f"
		i__2 = i__ - 1;
#line 379 "clabrd.f"
		clacgv_(&i__2, &y[i__ + y_dim1], ldy);
#line 380 "clabrd.f"
		i__2 = *m - i__;
#line 380 "clabrd.f"
		i__3 = i__ - 1;
#line 380 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 380 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &a[i__ + 1 + 
			a_dim1], lda, &y[i__ + y_dim1], ldy, &c_b2, &a[i__ + 
			1 + i__ * a_dim1], &c__1, (ftnlen)12);
#line 382 "clabrd.f"
		i__2 = i__ - 1;
#line 382 "clabrd.f"
		clacgv_(&i__2, &y[i__ + y_dim1], ldy);
#line 383 "clabrd.f"
		i__2 = *m - i__;
#line 383 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 383 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__, &z__1, &x[i__ + 1 + 
			x_dim1], ldx, &a[i__ * a_dim1 + 1], &c__1, &c_b2, &a[
			i__ + 1 + i__ * a_dim1], &c__1, (ftnlen)12);

/*              Generate reflection Q(i) to annihilate A(i+2:m,i) */

#line 388 "clabrd.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 388 "clabrd.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 389 "clabrd.f"
		i__2 = *m - i__;
/* Computing MIN */
#line 389 "clabrd.f"
		i__3 = i__ + 2;
#line 389 "clabrd.f"
		clarfg_(&i__2, &alpha, &a[min(i__3,*m) + i__ * a_dim1], &c__1,
			 &tauq[i__]);
#line 391 "clabrd.f"
		i__2 = i__;
#line 391 "clabrd.f"
		e[i__2] = alpha.r;
#line 392 "clabrd.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 392 "clabrd.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute Y(i+1:n,i) */

#line 396 "clabrd.f"
		i__2 = *m - i__;
#line 396 "clabrd.f"
		i__3 = *n - i__;
#line 396 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 
			+ (i__ + 1) * a_dim1], lda, &a[i__ + 1 + i__ * a_dim1]
			, &c__1, &c_b1, &y[i__ + 1 + i__ * y_dim1], &c__1, (
			ftnlen)19);
#line 399 "clabrd.f"
		i__2 = *m - i__;
#line 399 "clabrd.f"
		i__3 = i__ - 1;
#line 399 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 
			+ a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &
			c_b1, &y[i__ * y_dim1 + 1], &c__1, (ftnlen)19);
#line 402 "clabrd.f"
		i__2 = *n - i__;
#line 402 "clabrd.f"
		i__3 = i__ - 1;
#line 402 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 402 "clabrd.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &c_b2, &y[
			i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)12);
#line 404 "clabrd.f"
		i__2 = *m - i__;
#line 404 "clabrd.f"
		cgemv_("Conjugate transpose", &i__2, &i__, &c_b2, &x[i__ + 1 
			+ x_dim1], ldx, &a[i__ + 1 + i__ * a_dim1], &c__1, &
			c_b1, &y[i__ * y_dim1 + 1], &c__1, (ftnlen)19);
#line 407 "clabrd.f"
		i__2 = *n - i__;
#line 407 "clabrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 407 "clabrd.f"
		cgemv_("Conjugate transpose", &i__, &i__2, &z__1, &a[(i__ + 1)
			 * a_dim1 + 1], lda, &y[i__ * y_dim1 + 1], &c__1, &
			c_b2, &y[i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)19);
#line 410 "clabrd.f"
		i__2 = *n - i__;
#line 410 "clabrd.f"
		cscal_(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
#line 411 "clabrd.f"
	    } else {
#line 412 "clabrd.f"
		i__2 = *n - i__ + 1;
#line 412 "clabrd.f"
		clacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);
#line 413 "clabrd.f"
	    }
#line 414 "clabrd.f"
/* L20: */
#line 414 "clabrd.f"
	}
#line 415 "clabrd.f"
    }
#line 416 "clabrd.f"
    return 0;

/*     End of CLABRD */

} /* clabrd_ */

