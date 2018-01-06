#line 1 "dlatrd.f"
/* dlatrd.f -- translated by f2c (version 20100827).
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

#line 1 "dlatrd.f"
/* Table of constant values */

static doublereal c_b5 = -1.;
static doublereal c_b6 = 1.;
static integer c__1 = 1;
static doublereal c_b16 = 0.;

/* > \brief \b DLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiago
nal form by an orthogonal similarity transformation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLATRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATRD reduces NB rows and columns of a real symmetric matrix A to */
/* > symmetric tridiagonal form by an orthogonal similarity */
/* > transformation Q**T * A * Q, and returns the matrices V and W which are */
/* > needed to apply the transformation to the unreduced part of A. */
/* > */
/* > If UPLO = 'U', DLATRD reduces the last NB rows and columns of a */
/* > matrix, of which the upper triangle is supplied; */
/* > if UPLO = 'L', DLATRD reduces the first NB rows and columns of a */
/* > matrix, of which the lower triangle is supplied. */
/* > */
/* > This is an auxiliary routine called by DSYTRD. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored: */
/* >          = 'U': Upper triangular */
/* >          = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The number of rows and columns to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          n-by-n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n-by-n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit: */
/* >          if UPLO = 'U', the last NB columns have been reduced to */
/* >            tridiagonal form, with the diagonal elements overwriting */
/* >            the diagonal elements of A; the elements above the diagonal */
/* >            with the array TAU, represent the orthogonal matrix Q as a */
/* >            product of elementary reflectors; */
/* >          if UPLO = 'L', the first NB columns have been reduced to */
/* >            tridiagonal form, with the diagonal elements overwriting */
/* >            the diagonal elements of A; the elements below the diagonal */
/* >            with the array TAU, represent the  orthogonal matrix Q as a */
/* >            product of elementary reflectors. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= (1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal */
/* >          elements of the last NB columns of the reduced matrix; */
/* >          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of */
/* >          the first NB columns of the reduced matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors, stored in */
/* >          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (LDW,NB) */
/* >          The n-by-nb matrix W required to update the unreduced part */
/* >          of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDW */
/* > \verbatim */
/* >          LDW is INTEGER */
/* >          The leading dimension of the array W. LDW >= max(1,N). */
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
/* >  If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(n) H(n-1) . . . H(n-nb+1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i), */
/* >  and tau in TAU(i-1). */
/* > */
/* >  If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(nb). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), */
/* >  and tau in TAU(i). */
/* > */
/* >  The elements of the vectors v together form the n-by-nb matrix V */
/* >  which is needed, with W, to apply the transformation to the unreduced */
/* >  part of the matrix, using a symmetric rank-2k update of the form: */
/* >  A := A - V*W**T - W*V**T. */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with n = 5 and nb = 2: */
/* > */
/* >  if UPLO = 'U':                       if UPLO = 'L': */
/* > */
/* >    (  a   a   a   v4  v5 )              (  d                  ) */
/* >    (      a   a   v4  v5 )              (  1   d              ) */
/* >    (          a   1   v5 )              (  v1  1   a          ) */
/* >    (              d   1  )              (  v1  v2  a   a      ) */
/* >    (                  d  )              (  v1  v2  a   a   a  ) */
/* > */
/* >  where d denotes a diagonal element of the reduced matrix, a denotes */
/* >  an element of the original matrix that is unchanged, and vi denotes */
/* >  an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlatrd_(char *uplo, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, 
	integer *ldw, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, iw;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    dsymv_(char *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen), dlarfg_(integer *, doublereal *, doublereal *, integer *,
	     doublereal *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 239 "dlatrd.f"
    /* Parameter adjustments */
#line 239 "dlatrd.f"
    a_dim1 = *lda;
#line 239 "dlatrd.f"
    a_offset = 1 + a_dim1;
#line 239 "dlatrd.f"
    a -= a_offset;
#line 239 "dlatrd.f"
    --e;
#line 239 "dlatrd.f"
    --tau;
#line 239 "dlatrd.f"
    w_dim1 = *ldw;
#line 239 "dlatrd.f"
    w_offset = 1 + w_dim1;
#line 239 "dlatrd.f"
    w -= w_offset;
#line 239 "dlatrd.f"

#line 239 "dlatrd.f"
    /* Function Body */
#line 239 "dlatrd.f"
    if (*n <= 0) {
#line 239 "dlatrd.f"
	return 0;
#line 239 "dlatrd.f"
    }

#line 242 "dlatrd.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Reduce last NB columns of upper triangle */

#line 246 "dlatrd.f"
	i__1 = *n - *nb + 1;
#line 246 "dlatrd.f"
	for (i__ = *n; i__ >= i__1; --i__) {
#line 247 "dlatrd.f"
	    iw = i__ - *n + *nb;
#line 248 "dlatrd.f"
	    if (i__ < *n) {

/*              Update A(1:i,i) */

#line 252 "dlatrd.f"
		i__2 = *n - i__;
#line 252 "dlatrd.f"
		dgemv_("No transpose", &i__, &i__2, &c_b5, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &w[i__ + (iw + 1) * w_dim1], ldw, &
			c_b6, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 254 "dlatrd.f"
		i__2 = *n - i__;
#line 254 "dlatrd.f"
		dgemv_("No transpose", &i__, &i__2, &c_b5, &w[(iw + 1) * 
			w_dim1 + 1], ldw, &a[i__ + (i__ + 1) * a_dim1], lda, &
			c_b6, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 256 "dlatrd.f"
	    }
#line 257 "dlatrd.f"
	    if (i__ > 1) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(1:i-2,i) */

#line 262 "dlatrd.f"
		i__2 = i__ - 1;
#line 262 "dlatrd.f"
		dlarfg_(&i__2, &a[i__ - 1 + i__ * a_dim1], &a[i__ * a_dim1 + 
			1], &c__1, &tau[i__ - 1]);
#line 263 "dlatrd.f"
		e[i__ - 1] = a[i__ - 1 + i__ * a_dim1];
#line 264 "dlatrd.f"
		a[i__ - 1 + i__ * a_dim1] = 1.;

/*              Compute W(1:i-1,i) */

#line 268 "dlatrd.f"
		i__2 = i__ - 1;
#line 268 "dlatrd.f"
		dsymv_("Upper", &i__2, &c_b6, &a[a_offset], lda, &a[i__ * 
			a_dim1 + 1], &c__1, &c_b16, &w[iw * w_dim1 + 1], &
			c__1, (ftnlen)5);
#line 270 "dlatrd.f"
		if (i__ < *n) {
#line 271 "dlatrd.f"
		    i__2 = i__ - 1;
#line 271 "dlatrd.f"
		    i__3 = *n - i__;
#line 271 "dlatrd.f"
		    dgemv_("Transpose", &i__2, &i__3, &c_b6, &w[(iw + 1) * 
			    w_dim1 + 1], ldw, &a[i__ * a_dim1 + 1], &c__1, &
			    c_b16, &w[i__ + 1 + iw * w_dim1], &c__1, (ftnlen)
			    9);
#line 273 "dlatrd.f"
		    i__2 = i__ - 1;
#line 273 "dlatrd.f"
		    i__3 = *n - i__;
#line 273 "dlatrd.f"
		    dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[(i__ + 1) *
			     a_dim1 + 1], lda, &w[i__ + 1 + iw * w_dim1], &
			    c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1, (ftnlen)
			    12);
#line 276 "dlatrd.f"
		    i__2 = i__ - 1;
#line 276 "dlatrd.f"
		    i__3 = *n - i__;
#line 276 "dlatrd.f"
		    dgemv_("Transpose", &i__2, &i__3, &c_b6, &a[(i__ + 1) * 
			    a_dim1 + 1], lda, &a[i__ * a_dim1 + 1], &c__1, &
			    c_b16, &w[i__ + 1 + iw * w_dim1], &c__1, (ftnlen)
			    9);
#line 278 "dlatrd.f"
		    i__2 = i__ - 1;
#line 278 "dlatrd.f"
		    i__3 = *n - i__;
#line 278 "dlatrd.f"
		    dgemv_("No transpose", &i__2, &i__3, &c_b5, &w[(iw + 1) * 
			    w_dim1 + 1], ldw, &w[i__ + 1 + iw * w_dim1], &
			    c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1, (ftnlen)
			    12);
#line 281 "dlatrd.f"
		}
#line 282 "dlatrd.f"
		i__2 = i__ - 1;
#line 282 "dlatrd.f"
		dscal_(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
#line 283 "dlatrd.f"
		i__2 = i__ - 1;
#line 283 "dlatrd.f"
		alpha = tau[i__ - 1] * -.5 * ddot_(&i__2, &w[iw * w_dim1 + 1],
			 &c__1, &a[i__ * a_dim1 + 1], &c__1);
#line 285 "dlatrd.f"
		i__2 = i__ - 1;
#line 285 "dlatrd.f"
		daxpy_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &w[iw * 
			w_dim1 + 1], &c__1);
#line 286 "dlatrd.f"
	    }

#line 288 "dlatrd.f"
/* L10: */
#line 288 "dlatrd.f"
	}
#line 289 "dlatrd.f"
    } else {

/*        Reduce first NB columns of lower triangle */

#line 293 "dlatrd.f"
	i__1 = *nb;
#line 293 "dlatrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i:n,i) */

#line 297 "dlatrd.f"
	    i__2 = *n - i__ + 1;
#line 297 "dlatrd.f"
	    i__3 = i__ - 1;
#line 297 "dlatrd.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + a_dim1], lda,
		     &w[i__ + w_dim1], ldw, &c_b6, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 299 "dlatrd.f"
	    i__2 = *n - i__ + 1;
#line 299 "dlatrd.f"
	    i__3 = i__ - 1;
#line 299 "dlatrd.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b5, &w[i__ + w_dim1], ldw,
		     &a[i__ + a_dim1], lda, &c_b6, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 301 "dlatrd.f"
	    if (i__ < *n) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(i+2:n,i) */

#line 306 "dlatrd.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 306 "dlatrd.f"
		i__3 = i__ + 2;
#line 306 "dlatrd.f"
		dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + 
			i__ * a_dim1], &c__1, &tau[i__]);
#line 308 "dlatrd.f"
		e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 309 "dlatrd.f"
		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Compute W(i+1:n,i) */

#line 313 "dlatrd.f"
		i__2 = *n - i__;
#line 313 "dlatrd.f"
		dsymv_("Lower", &i__2, &c_b6, &a[i__ + 1 + (i__ + 1) * a_dim1]
			, lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)5);
#line 315 "dlatrd.f"
		i__2 = *n - i__;
#line 315 "dlatrd.f"
		i__3 = i__ - 1;
#line 315 "dlatrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b6, &w[i__ + 1 + w_dim1],
			 ldw, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
			i__ * w_dim1 + 1], &c__1, (ftnlen)9);
#line 317 "dlatrd.f"
		i__2 = *n - i__;
#line 317 "dlatrd.f"
		i__3 = i__ - 1;
#line 317 "dlatrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + 
			a_dim1], lda, &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)12);
#line 319 "dlatrd.f"
		i__2 = *n - i__;
#line 319 "dlatrd.f"
		i__3 = i__ - 1;
#line 319 "dlatrd.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b6, &a[i__ + 1 + a_dim1],
			 lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
			i__ * w_dim1 + 1], &c__1, (ftnlen)9);
#line 321 "dlatrd.f"
		i__2 = *n - i__;
#line 321 "dlatrd.f"
		i__3 = i__ - 1;
#line 321 "dlatrd.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &w[i__ + 1 + 
			w_dim1], ldw, &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)12);
#line 323 "dlatrd.f"
		i__2 = *n - i__;
#line 323 "dlatrd.f"
		dscal_(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
#line 324 "dlatrd.f"
		i__2 = *n - i__;
#line 324 "dlatrd.f"
		alpha = tau[i__] * -.5 * ddot_(&i__2, &w[i__ + 1 + i__ * 
			w_dim1], &c__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 326 "dlatrd.f"
		i__2 = *n - i__;
#line 326 "dlatrd.f"
		daxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &w[
			i__ + 1 + i__ * w_dim1], &c__1);
#line 327 "dlatrd.f"
	    }

#line 329 "dlatrd.f"
/* L20: */
#line 329 "dlatrd.f"
	}
#line 330 "dlatrd.f"
    }

#line 332 "dlatrd.f"
    return 0;

/*     End of DLATRD */

} /* dlatrd_ */

