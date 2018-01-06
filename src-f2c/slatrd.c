#line 1 "slatrd.f"
/* slatrd.f -- translated by f2c (version 20100827).
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

#line 1 "slatrd.f"
/* Table of constant values */

static doublereal c_b5 = -1.;
static doublereal c_b6 = 1.;
static integer c__1 = 1;
static doublereal c_b16 = 0.;

/* > \brief \b SLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiago
nal form by an orthogonal similarity transformation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLATRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), E( * ), TAU( * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLATRD reduces NB rows and columns of a real symmetric matrix A to */
/* > symmetric tridiagonal form by an orthogonal similarity */
/* > transformation Q**T * A * Q, and returns the matrices V and W which are */
/* > needed to apply the transformation to the unreduced part of A. */
/* > */
/* > If UPLO = 'U', SLATRD reduces the last NB rows and columns of a */
/* > matrix, of which the upper triangle is supplied; */
/* > if UPLO = 'L', SLATRD reduces the first NB rows and columns of a */
/* > matrix, of which the lower triangle is supplied. */
/* > */
/* > This is an auxiliary routine called by SSYTRD. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          E is REAL array, dimension (N-1) */
/* >          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal */
/* >          elements of the last NB columns of the reduced matrix; */
/* >          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of */
/* >          the first NB columns of the reduced matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors, stored in */
/* >          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (LDW,NB) */
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

/* > \date September 2012 */

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
/* Subroutine */ int slatrd_(char *uplo, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, 
	integer *ldw, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, iw;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), ssymv_(char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen), 
	    slarfg_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 239 "slatrd.f"
    /* Parameter adjustments */
#line 239 "slatrd.f"
    a_dim1 = *lda;
#line 239 "slatrd.f"
    a_offset = 1 + a_dim1;
#line 239 "slatrd.f"
    a -= a_offset;
#line 239 "slatrd.f"
    --e;
#line 239 "slatrd.f"
    --tau;
#line 239 "slatrd.f"
    w_dim1 = *ldw;
#line 239 "slatrd.f"
    w_offset = 1 + w_dim1;
#line 239 "slatrd.f"
    w -= w_offset;
#line 239 "slatrd.f"

#line 239 "slatrd.f"
    /* Function Body */
#line 239 "slatrd.f"
    if (*n <= 0) {
#line 239 "slatrd.f"
	return 0;
#line 239 "slatrd.f"
    }

#line 242 "slatrd.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Reduce last NB columns of upper triangle */

#line 246 "slatrd.f"
	i__1 = *n - *nb + 1;
#line 246 "slatrd.f"
	for (i__ = *n; i__ >= i__1; --i__) {
#line 247 "slatrd.f"
	    iw = i__ - *n + *nb;
#line 248 "slatrd.f"
	    if (i__ < *n) {

/*              Update A(1:i,i) */

#line 252 "slatrd.f"
		i__2 = *n - i__;
#line 252 "slatrd.f"
		sgemv_("No transpose", &i__, &i__2, &c_b5, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &w[i__ + (iw + 1) * w_dim1], ldw, &
			c_b6, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 254 "slatrd.f"
		i__2 = *n - i__;
#line 254 "slatrd.f"
		sgemv_("No transpose", &i__, &i__2, &c_b5, &w[(iw + 1) * 
			w_dim1 + 1], ldw, &a[i__ + (i__ + 1) * a_dim1], lda, &
			c_b6, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 256 "slatrd.f"
	    }
#line 257 "slatrd.f"
	    if (i__ > 1) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(1:i-2,i) */

#line 262 "slatrd.f"
		i__2 = i__ - 1;
#line 262 "slatrd.f"
		slarfg_(&i__2, &a[i__ - 1 + i__ * a_dim1], &a[i__ * a_dim1 + 
			1], &c__1, &tau[i__ - 1]);
#line 263 "slatrd.f"
		e[i__ - 1] = a[i__ - 1 + i__ * a_dim1];
#line 264 "slatrd.f"
		a[i__ - 1 + i__ * a_dim1] = 1.;

/*              Compute W(1:i-1,i) */

#line 268 "slatrd.f"
		i__2 = i__ - 1;
#line 268 "slatrd.f"
		ssymv_("Upper", &i__2, &c_b6, &a[a_offset], lda, &a[i__ * 
			a_dim1 + 1], &c__1, &c_b16, &w[iw * w_dim1 + 1], &
			c__1, (ftnlen)5);
#line 270 "slatrd.f"
		if (i__ < *n) {
#line 271 "slatrd.f"
		    i__2 = i__ - 1;
#line 271 "slatrd.f"
		    i__3 = *n - i__;
#line 271 "slatrd.f"
		    sgemv_("Transpose", &i__2, &i__3, &c_b6, &w[(iw + 1) * 
			    w_dim1 + 1], ldw, &a[i__ * a_dim1 + 1], &c__1, &
			    c_b16, &w[i__ + 1 + iw * w_dim1], &c__1, (ftnlen)
			    9);
#line 273 "slatrd.f"
		    i__2 = i__ - 1;
#line 273 "slatrd.f"
		    i__3 = *n - i__;
#line 273 "slatrd.f"
		    sgemv_("No transpose", &i__2, &i__3, &c_b5, &a[(i__ + 1) *
			     a_dim1 + 1], lda, &w[i__ + 1 + iw * w_dim1], &
			    c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1, (ftnlen)
			    12);
#line 276 "slatrd.f"
		    i__2 = i__ - 1;
#line 276 "slatrd.f"
		    i__3 = *n - i__;
#line 276 "slatrd.f"
		    sgemv_("Transpose", &i__2, &i__3, &c_b6, &a[(i__ + 1) * 
			    a_dim1 + 1], lda, &a[i__ * a_dim1 + 1], &c__1, &
			    c_b16, &w[i__ + 1 + iw * w_dim1], &c__1, (ftnlen)
			    9);
#line 278 "slatrd.f"
		    i__2 = i__ - 1;
#line 278 "slatrd.f"
		    i__3 = *n - i__;
#line 278 "slatrd.f"
		    sgemv_("No transpose", &i__2, &i__3, &c_b5, &w[(iw + 1) * 
			    w_dim1 + 1], ldw, &w[i__ + 1 + iw * w_dim1], &
			    c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1, (ftnlen)
			    12);
#line 281 "slatrd.f"
		}
#line 282 "slatrd.f"
		i__2 = i__ - 1;
#line 282 "slatrd.f"
		sscal_(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
#line 283 "slatrd.f"
		i__2 = i__ - 1;
#line 283 "slatrd.f"
		alpha = tau[i__ - 1] * -.5 * sdot_(&i__2, &w[iw * w_dim1 + 1],
			 &c__1, &a[i__ * a_dim1 + 1], &c__1);
#line 285 "slatrd.f"
		i__2 = i__ - 1;
#line 285 "slatrd.f"
		saxpy_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &w[iw * 
			w_dim1 + 1], &c__1);
#line 286 "slatrd.f"
	    }

#line 288 "slatrd.f"
/* L10: */
#line 288 "slatrd.f"
	}
#line 289 "slatrd.f"
    } else {

/*        Reduce first NB columns of lower triangle */

#line 293 "slatrd.f"
	i__1 = *nb;
#line 293 "slatrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i:n,i) */

#line 297 "slatrd.f"
	    i__2 = *n - i__ + 1;
#line 297 "slatrd.f"
	    i__3 = i__ - 1;
#line 297 "slatrd.f"
	    sgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + a_dim1], lda,
		     &w[i__ + w_dim1], ldw, &c_b6, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 299 "slatrd.f"
	    i__2 = *n - i__ + 1;
#line 299 "slatrd.f"
	    i__3 = i__ - 1;
#line 299 "slatrd.f"
	    sgemv_("No transpose", &i__2, &i__3, &c_b5, &w[i__ + w_dim1], ldw,
		     &a[i__ + a_dim1], lda, &c_b6, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 301 "slatrd.f"
	    if (i__ < *n) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(i+2:n,i) */

#line 306 "slatrd.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 306 "slatrd.f"
		i__3 = i__ + 2;
#line 306 "slatrd.f"
		slarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + 
			i__ * a_dim1], &c__1, &tau[i__]);
#line 308 "slatrd.f"
		e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 309 "slatrd.f"
		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Compute W(i+1:n,i) */

#line 313 "slatrd.f"
		i__2 = *n - i__;
#line 313 "slatrd.f"
		ssymv_("Lower", &i__2, &c_b6, &a[i__ + 1 + (i__ + 1) * a_dim1]
			, lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)5);
#line 315 "slatrd.f"
		i__2 = *n - i__;
#line 315 "slatrd.f"
		i__3 = i__ - 1;
#line 315 "slatrd.f"
		sgemv_("Transpose", &i__2, &i__3, &c_b6, &w[i__ + 1 + w_dim1],
			 ldw, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
			i__ * w_dim1 + 1], &c__1, (ftnlen)9);
#line 317 "slatrd.f"
		i__2 = *n - i__;
#line 317 "slatrd.f"
		i__3 = i__ - 1;
#line 317 "slatrd.f"
		sgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + 
			a_dim1], lda, &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)12);
#line 319 "slatrd.f"
		i__2 = *n - i__;
#line 319 "slatrd.f"
		i__3 = i__ - 1;
#line 319 "slatrd.f"
		sgemv_("Transpose", &i__2, &i__3, &c_b6, &a[i__ + 1 + a_dim1],
			 lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
			i__ * w_dim1 + 1], &c__1, (ftnlen)9);
#line 321 "slatrd.f"
		i__2 = *n - i__;
#line 321 "slatrd.f"
		i__3 = i__ - 1;
#line 321 "slatrd.f"
		sgemv_("No transpose", &i__2, &i__3, &c_b5, &w[i__ + 1 + 
			w_dim1], ldw, &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)12);
#line 323 "slatrd.f"
		i__2 = *n - i__;
#line 323 "slatrd.f"
		sscal_(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
#line 324 "slatrd.f"
		i__2 = *n - i__;
#line 324 "slatrd.f"
		alpha = tau[i__] * -.5 * sdot_(&i__2, &w[i__ + 1 + i__ * 
			w_dim1], &c__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 326 "slatrd.f"
		i__2 = *n - i__;
#line 326 "slatrd.f"
		saxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &w[
			i__ + 1 + i__ * w_dim1], &c__1);
#line 327 "slatrd.f"
	    }

#line 329 "slatrd.f"
/* L20: */
#line 329 "slatrd.f"
	}
#line 330 "slatrd.f"
    }

#line 332 "slatrd.f"
    return 0;

/*     End of SLATRD */

} /* slatrd_ */

