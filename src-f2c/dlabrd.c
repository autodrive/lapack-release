#line 1 "dlabrd.f"
/* dlabrd.f -- translated by f2c (version 20100827).
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

#line 1 "dlabrd.f"
/* Table of constant values */

static doublereal c_b4 = -1.;
static doublereal c_b5 = 1.;
static integer c__1 = 1;
static doublereal c_b16 = 0.;

/* > \brief \b DLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLABRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlabrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlabrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, */
/*                          LDY ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, LDX, LDY, M, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ), */
/*      $                   TAUQ( * ), X( LDX, * ), Y( LDY, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLABRD reduces the first NB rows and columns of a real general */
/* > m by n matrix A to upper or lower bidiagonal form by an orthogonal */
/* > transformation Q**T * A * P, and returns the matrices X and Y which */
/* > are needed to apply the transformation to the unreduced part of A. */
/* > */
/* > If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower */
/* > bidiagonal form. */
/* > */
/* > This is an auxiliary routine called by DGEBRD */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the m by n general matrix to be reduced. */
/* >          On exit, the first NB rows and columns of the matrix are */
/* >          overwritten; the rest of the array is unchanged. */
/* >          If m >= n, elements on and below the diagonal in the first NB */
/* >            columns, with the array TAUQ, represent the orthogonal */
/* >            matrix Q as a product of elementary reflectors; and */
/* >            elements above the diagonal in the first NB rows, with the */
/* >            array TAUP, represent the orthogonal matrix P as a product */
/* >            of elementary reflectors. */
/* >          If m < n, elements below the diagonal in the first NB */
/* >            columns, with the array TAUQ, represent the orthogonal */
/* >            matrix Q as a product of elementary reflectors, and */
/* >            elements on and above the diagonal in the first NB rows, */
/* >            with the array TAUP, represent the orthogonal matrix P as */
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
/* >          D is DOUBLE PRECISION array, dimension (NB) */
/* >          The diagonal elements of the first NB rows and columns of */
/* >          the reduced matrix.  D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (NB) */
/* >          The off-diagonal elements of the first NB rows and columns of */
/* >          the reduced matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* >          TAUQ is DOUBLE PRECISION array dimension (NB) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* >          TAUP is DOUBLE PRECISION array, dimension (NB) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NB) */
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
/* >          Y is DOUBLE PRECISION array, dimension (LDY,NB) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

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
/* >     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T */
/* > */
/* >  where tauq and taup are real scalars, and v and u are real vectors. */
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
/* >  V and the nb-by-n matrix U**T which are needed, with X and Y, to apply */
/* >  the transformation to the unreduced part of the matrix, using a block */
/* >  update of the form:  A := A - V*Y**T - X*U**T. */
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
/* Subroutine */ int dlabrd_(integer *m, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, 
	doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer 
	*ldy)
{
    /* System generated locals */
    integer a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *);


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

#line 245 "dlabrd.f"
    /* Parameter adjustments */
#line 245 "dlabrd.f"
    a_dim1 = *lda;
#line 245 "dlabrd.f"
    a_offset = 1 + a_dim1;
#line 245 "dlabrd.f"
    a -= a_offset;
#line 245 "dlabrd.f"
    --d__;
#line 245 "dlabrd.f"
    --e;
#line 245 "dlabrd.f"
    --tauq;
#line 245 "dlabrd.f"
    --taup;
#line 245 "dlabrd.f"
    x_dim1 = *ldx;
#line 245 "dlabrd.f"
    x_offset = 1 + x_dim1;
#line 245 "dlabrd.f"
    x -= x_offset;
#line 245 "dlabrd.f"
    y_dim1 = *ldy;
#line 245 "dlabrd.f"
    y_offset = 1 + y_dim1;
#line 245 "dlabrd.f"
    y -= y_offset;
#line 245 "dlabrd.f"

#line 245 "dlabrd.f"
    /* Function Body */
#line 245 "dlabrd.f"
    if (*m <= 0 || *n <= 0) {
#line 245 "dlabrd.f"
	return 0;
#line 245 "dlabrd.f"
    }

#line 248 "dlabrd.f"
    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

#line 252 "dlabrd.f"
	i__1 = *nb;
#line 252 "dlabrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i:m,i) */

#line 256 "dlabrd.f"
	    i__2 = *m - i__ + 1;
#line 256 "dlabrd.f"
	    i__3 = i__ - 1;
#line 256 "dlabrd.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b4, &a[i__ + a_dim1], lda,
		     &y[i__ + y_dim1], ldy, &c_b5, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 258 "dlabrd.f"
	    i__2 = *m - i__ + 1;
#line 258 "dlabrd.f"
	    i__3 = i__ - 1;
#line 258 "dlabrd.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b4, &x[i__ + x_dim1], ldx,
		     &a[i__ * a_dim1 + 1], &c__1, &c_b5, &a[i__ + i__ * 
		    a_dim1], &c__1, (ftnlen)12);

/*           Generate reflection Q(i) to annihilate A(i+1:m,i) */

#line 263 "dlabrd.f"
	    i__2 = *m - i__ + 1;
/* Computing MIN */
#line 263 "dlabrd.f"
	    i__3 = i__ + 1;
#line 263 "dlabrd.f"
	    dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * 
		    a_dim1], &c__1, &tauq[i__]);
#line 265 "dlabrd.f"
	    d__[i__] = a[i__ + i__ * a_dim1];
#line 266 "dlabrd.f"
	    if (i__ < *n) {
#line 267 "dlabrd.f"
		a[i__ + i__ * a_dim1] = 1.;

/*              Compute Y(i+1:n,i) */

#line 271 "dlabrd.f"
		i__2 = *m - i__ + 1;
#line 271 "dlabrd.f"
		i__3 = *n - i__;
#line 271 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b5, &a[i__ + (i__ + 1) * 
			a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &c_b16, &
			y[i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)9);
#line 273 "dlabrd.f"
		i__2 = *m - i__ + 1;
#line 273 "dlabrd.f"
		i__3 = i__ - 1;
#line 273 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b5, &a[i__ + a_dim1], 
			lda, &a[i__ + i__ * a_dim1], &c__1, &c_b16, &y[i__ * 
			y_dim1 + 1], &c__1, (ftnlen)9);
#line 275 "dlabrd.f"
		i__2 = *n - i__;
#line 275 "dlabrd.f"
		i__3 = i__ - 1;
#line 275 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b4, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[
			i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)12);
#line 277 "dlabrd.f"
		i__2 = *m - i__ + 1;
#line 277 "dlabrd.f"
		i__3 = i__ - 1;
#line 277 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b5, &x[i__ + x_dim1], 
			ldx, &a[i__ + i__ * a_dim1], &c__1, &c_b16, &y[i__ * 
			y_dim1 + 1], &c__1, (ftnlen)9);
#line 279 "dlabrd.f"
		i__2 = i__ - 1;
#line 279 "dlabrd.f"
		i__3 = *n - i__;
#line 279 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b4, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &y[i__ * y_dim1 + 1], &c__1, &c_b5, 
			&y[i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)9);
#line 281 "dlabrd.f"
		i__2 = *n - i__;
#line 281 "dlabrd.f"
		dscal_(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);

/*              Update A(i,i+1:n) */

#line 285 "dlabrd.f"
		i__2 = *n - i__;
#line 285 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__, &c_b4, &y[i__ + 1 + 
			y_dim1], ldy, &a[i__ + a_dim1], lda, &c_b5, &a[i__ + (
			i__ + 1) * a_dim1], lda, (ftnlen)12);
#line 287 "dlabrd.f"
		i__2 = i__ - 1;
#line 287 "dlabrd.f"
		i__3 = *n - i__;
#line 287 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b4, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &c_b5, &a[
			i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);

/*              Generate reflection P(i) to annihilate A(i,i+2:n) */

#line 292 "dlabrd.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 292 "dlabrd.f"
		i__3 = i__ + 2;
#line 292 "dlabrd.f"
		dlarfg_(&i__2, &a[i__ + (i__ + 1) * a_dim1], &a[i__ + min(
			i__3,*n) * a_dim1], lda, &taup[i__]);
#line 294 "dlabrd.f"
		e[i__] = a[i__ + (i__ + 1) * a_dim1];
#line 295 "dlabrd.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;

/*              Compute X(i+1:m,i) */

#line 299 "dlabrd.f"
		i__2 = *m - i__;
#line 299 "dlabrd.f"
		i__3 = *n - i__;
#line 299 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + (i__ 
			+ 1) * a_dim1], lda, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &c_b16, &x[i__ + 1 + i__ * x_dim1], &c__1, (
			ftnlen)12);
#line 301 "dlabrd.f"
		i__2 = *n - i__;
#line 301 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__, &c_b5, &y[i__ + 1 + y_dim1], 
			ldy, &a[i__ + (i__ + 1) * a_dim1], lda, &c_b16, &x[
			i__ * x_dim1 + 1], &c__1, (ftnlen)9);
#line 303 "dlabrd.f"
		i__2 = *m - i__;
#line 303 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__, &c_b4, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 305 "dlabrd.f"
		i__2 = i__ - 1;
#line 305 "dlabrd.f"
		i__3 = *n - i__;
#line 305 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &a[i__ + (i__ + 1) * a_dim1], lda, &
			c_b16, &x[i__ * x_dim1 + 1], &c__1, (ftnlen)12);
#line 307 "dlabrd.f"
		i__2 = *m - i__;
#line 307 "dlabrd.f"
		i__3 = i__ - 1;
#line 307 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b4, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 309 "dlabrd.f"
		i__2 = *m - i__;
#line 309 "dlabrd.f"
		dscal_(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
#line 310 "dlabrd.f"
	    }
#line 311 "dlabrd.f"
/* L10: */
#line 311 "dlabrd.f"
	}
#line 312 "dlabrd.f"
    } else {

/*        Reduce to lower bidiagonal form */

#line 316 "dlabrd.f"
	i__1 = *nb;
#line 316 "dlabrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i,i:n) */

#line 320 "dlabrd.f"
	    i__2 = *n - i__ + 1;
#line 320 "dlabrd.f"
	    i__3 = i__ - 1;
#line 320 "dlabrd.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b4, &y[i__ + y_dim1], ldy,
		     &a[i__ + a_dim1], lda, &c_b5, &a[i__ + i__ * a_dim1], 
		    lda, (ftnlen)12);
#line 322 "dlabrd.f"
	    i__2 = i__ - 1;
#line 322 "dlabrd.f"
	    i__3 = *n - i__ + 1;
#line 322 "dlabrd.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b4, &a[i__ * a_dim1 + 1], 
		    lda, &x[i__ + x_dim1], ldx, &c_b5, &a[i__ + i__ * a_dim1],
		     lda, (ftnlen)9);

/*           Generate reflection P(i) to annihilate A(i,i+1:n) */

#line 327 "dlabrd.f"
	    i__2 = *n - i__ + 1;
/* Computing MIN */
#line 327 "dlabrd.f"
	    i__3 = i__ + 1;
#line 327 "dlabrd.f"
	    dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + min(i__3,*n) * 
		    a_dim1], lda, &taup[i__]);
#line 329 "dlabrd.f"
	    d__[i__] = a[i__ + i__ * a_dim1];
#line 330 "dlabrd.f"
	    if (i__ < *m) {
#line 331 "dlabrd.f"
		a[i__ + i__ * a_dim1] = 1.;

/*              Compute X(i+1:m,i) */

#line 335 "dlabrd.f"
		i__2 = *m - i__;
#line 335 "dlabrd.f"
		i__3 = *n - i__ + 1;
#line 335 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + i__ *
			 a_dim1], lda, &a[i__ + i__ * a_dim1], lda, &c_b16, &
			x[i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 337 "dlabrd.f"
		i__2 = *n - i__ + 1;
#line 337 "dlabrd.f"
		i__3 = i__ - 1;
#line 337 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b5, &y[i__ + y_dim1], 
			ldy, &a[i__ + i__ * a_dim1], lda, &c_b16, &x[i__ * 
			x_dim1 + 1], &c__1, (ftnlen)9);
#line 339 "dlabrd.f"
		i__2 = *m - i__;
#line 339 "dlabrd.f"
		i__3 = i__ - 1;
#line 339 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b4, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 341 "dlabrd.f"
		i__2 = i__ - 1;
#line 341 "dlabrd.f"
		i__3 = *n - i__ + 1;
#line 341 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ * a_dim1 + 
			1], lda, &a[i__ + i__ * a_dim1], lda, &c_b16, &x[i__ *
			 x_dim1 + 1], &c__1, (ftnlen)12);
#line 343 "dlabrd.f"
		i__2 = *m - i__;
#line 343 "dlabrd.f"
		i__3 = i__ - 1;
#line 343 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b4, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[
			i__ + 1 + i__ * x_dim1], &c__1, (ftnlen)12);
#line 345 "dlabrd.f"
		i__2 = *m - i__;
#line 345 "dlabrd.f"
		dscal_(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);

/*              Update A(i+1:m,i) */

#line 349 "dlabrd.f"
		i__2 = *m - i__;
#line 349 "dlabrd.f"
		i__3 = i__ - 1;
#line 349 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b4, &a[i__ + 1 + 
			a_dim1], lda, &y[i__ + y_dim1], ldy, &c_b5, &a[i__ + 
			1 + i__ * a_dim1], &c__1, (ftnlen)12);
#line 351 "dlabrd.f"
		i__2 = *m - i__;
#line 351 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__, &c_b4, &x[i__ + 1 + 
			x_dim1], ldx, &a[i__ * a_dim1 + 1], &c__1, &c_b5, &a[
			i__ + 1 + i__ * a_dim1], &c__1, (ftnlen)12);

/*              Generate reflection Q(i) to annihilate A(i+2:m,i) */

#line 356 "dlabrd.f"
		i__2 = *m - i__;
/* Computing MIN */
#line 356 "dlabrd.f"
		i__3 = i__ + 2;
#line 356 "dlabrd.f"
		dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*m) + 
			i__ * a_dim1], &c__1, &tauq[i__]);
#line 358 "dlabrd.f"
		e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 359 "dlabrd.f"
		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Compute Y(i+1:n,i) */

#line 363 "dlabrd.f"
		i__2 = *m - i__;
#line 363 "dlabrd.f"
		i__3 = *n - i__;
#line 363 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + (i__ + 
			1) * a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, 
			&c_b16, &y[i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)9);
#line 365 "dlabrd.f"
		i__2 = *m - i__;
#line 365 "dlabrd.f"
		i__3 = i__ - 1;
#line 365 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + a_dim1],
			 lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &y[
			i__ * y_dim1 + 1], &c__1, (ftnlen)9);
#line 367 "dlabrd.f"
		i__2 = *n - i__;
#line 367 "dlabrd.f"
		i__3 = i__ - 1;
#line 367 "dlabrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b4, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[
			i__ + 1 + i__ * y_dim1], &c__1, (ftnlen)12);
#line 369 "dlabrd.f"
		i__2 = *m - i__;
#line 369 "dlabrd.f"
		dgemv_("Transpose", &i__2, &i__, &c_b5, &x[i__ + 1 + x_dim1], 
			ldx, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &y[
			i__ * y_dim1 + 1], &c__1, (ftnlen)9);
#line 371 "dlabrd.f"
		i__2 = *n - i__;
#line 371 "dlabrd.f"
		dgemv_("Transpose", &i__, &i__2, &c_b4, &a[(i__ + 1) * a_dim1 
			+ 1], lda, &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[i__ 
			+ 1 + i__ * y_dim1], &c__1, (ftnlen)9);
#line 373 "dlabrd.f"
		i__2 = *n - i__;
#line 373 "dlabrd.f"
		dscal_(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
#line 374 "dlabrd.f"
	    }
#line 375 "dlabrd.f"
/* L20: */
#line 375 "dlabrd.f"
	}
#line 376 "dlabrd.f"
    }
#line 377 "dlabrd.f"
    return 0;

/*     End of DLABRD */

} /* dlabrd_ */

