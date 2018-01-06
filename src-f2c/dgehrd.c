#line 1 "dgehrd.f"
/* dgehrd.f -- translated by f2c (version 20100827).
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

#line 1 "dgehrd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__65 = 65;
static doublereal c_b25 = -1.;
static doublereal c_b26 = 1.;

/* > \brief \b DGEHRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEHRD reduces a real general matrix A to upper Hessenberg form H by */
/* > an orthogonal similarity transformation:  Q**T * A * Q = H . */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* > */
/* >          It is assumed that A is already upper triangular in rows */
/* >          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/* >          set by a previous call to DGEBAL; otherwise they should be */
/* >          set to 1 and N respectively. See Further Details. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the N-by-N general matrix to be reduced. */
/* >          On exit, the upper triangle and the first subdiagonal of A */
/* >          are overwritten with the upper Hessenberg matrix H, and the */
/* >          elements below the first subdiagonal, with the array TAU, */
/* >          represent the orthogonal matrix Q as a product of elementary */
/* >          reflectors. See Further Details. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to */
/* >          zero. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= max(1,N). */
/* >          For optimum performance LWORK >= N*NB, where NB is the */
/* >          optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of (ihi-ilo) elementary */
/* >  reflectors */
/* > */
/* >     Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on */
/* >  exit in A(i+2:ihi,i), and tau in TAU(i). */
/* > */
/* >  The contents of A are illustrated by the following example, with */
/* >  n = 7, ilo = 2 and ihi = 6: */
/* > */
/* >  on entry,                        on exit, */
/* > */
/* >  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a ) */
/* >  (     a   a   a   a   a   a )    (      a   h   h   h   h   a ) */
/* >  (     a   a   a   a   a   a )    (      h   h   h   h   h   h ) */
/* >  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h ) */
/* >  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h ) */
/* >  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h ) */
/* >  (                         a )    (                          a ) */
/* > */
/* >  where a denotes an element of the original matrix A, h denotes a */
/* >  modified element of the upper Hessenberg matrix H, and vi denotes an */
/* >  element of the vector defining H(i). */
/* > */
/* >  This file is a slight modification of LAPACK-3.0's DGEHRD */
/* >  subroutine incorporating improvements proposed by Quintana-Orti and */
/* >  Van de Geijn (2006). (See DLAHR2.) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j;
    static doublereal t[4160]	/* was [65][64] */;
    static integer ib;
    static doublereal ei;
    static integer nb, nh, nx, iws;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dgehd2_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *), dlahr2_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    dlarfb_(char *, char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 216 "dgehrd.f"
    /* Parameter adjustments */
#line 216 "dgehrd.f"
    a_dim1 = *lda;
#line 216 "dgehrd.f"
    a_offset = 1 + a_dim1;
#line 216 "dgehrd.f"
    a -= a_offset;
#line 216 "dgehrd.f"
    --tau;
#line 216 "dgehrd.f"
    --work;
#line 216 "dgehrd.f"

#line 216 "dgehrd.f"
    /* Function Body */
#line 216 "dgehrd.f"
    *info = 0;
/* Computing MIN */
#line 217 "dgehrd.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 217 "dgehrd.f"
    nb = min(i__1,i__2);
#line 218 "dgehrd.f"
    lwkopt = *n * nb;
#line 219 "dgehrd.f"
    work[1] = (doublereal) lwkopt;
#line 220 "dgehrd.f"
    lquery = *lwork == -1;
#line 221 "dgehrd.f"
    if (*n < 0) {
#line 222 "dgehrd.f"
	*info = -1;
#line 223 "dgehrd.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 224 "dgehrd.f"
	*info = -2;
#line 225 "dgehrd.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 226 "dgehrd.f"
	*info = -3;
#line 227 "dgehrd.f"
    } else if (*lda < max(1,*n)) {
#line 228 "dgehrd.f"
	*info = -5;
#line 229 "dgehrd.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 230 "dgehrd.f"
	*info = -8;
#line 231 "dgehrd.f"
    }
#line 232 "dgehrd.f"
    if (*info != 0) {
#line 233 "dgehrd.f"
	i__1 = -(*info);
#line 233 "dgehrd.f"
	xerbla_("DGEHRD", &i__1, (ftnlen)6);
#line 234 "dgehrd.f"
	return 0;
#line 235 "dgehrd.f"
    } else if (lquery) {
#line 236 "dgehrd.f"
	return 0;
#line 237 "dgehrd.f"
    }

/*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */

#line 241 "dgehrd.f"
    i__1 = *ilo - 1;
#line 241 "dgehrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "dgehrd.f"
	tau[i__] = 0.;
#line 243 "dgehrd.f"
/* L10: */
#line 243 "dgehrd.f"
    }
#line 244 "dgehrd.f"
    i__1 = *n - 1;
#line 244 "dgehrd.f"
    for (i__ = max(1,*ihi); i__ <= i__1; ++i__) {
#line 245 "dgehrd.f"
	tau[i__] = 0.;
#line 246 "dgehrd.f"
/* L20: */
#line 246 "dgehrd.f"
    }

/*     Quick return if possible */

#line 250 "dgehrd.f"
    nh = *ihi - *ilo + 1;
#line 251 "dgehrd.f"
    if (nh <= 1) {
#line 252 "dgehrd.f"
	work[1] = 1.;
#line 253 "dgehrd.f"
	return 0;
#line 254 "dgehrd.f"
    }

/*     Determine the block size */

/* Computing MIN */
#line 258 "dgehrd.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 258 "dgehrd.f"
    nb = min(i__1,i__2);
#line 259 "dgehrd.f"
    nbmin = 2;
#line 260 "dgehrd.f"
    iws = 1;
#line 261 "dgehrd.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to cross over from blocked to unblocked code */
/*        (last block is always handled by unblocked code) */

/* Computing MAX */
#line 266 "dgehrd.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "DGEHRD", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 266 "dgehrd.f"
	nx = max(i__1,i__2);
#line 267 "dgehrd.f"
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code */

#line 271 "dgehrd.f"
	    iws = *n * nb;
#line 272 "dgehrd.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code */

/* Computing MAX */
#line 278 "dgehrd.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGEHRD", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 278 "dgehrd.f"
		nbmin = max(i__1,i__2);
#line 280 "dgehrd.f"
		if (*lwork >= *n * nbmin) {
#line 281 "dgehrd.f"
		    nb = *lwork / *n;
#line 282 "dgehrd.f"
		} else {
#line 283 "dgehrd.f"
		    nb = 1;
#line 284 "dgehrd.f"
		}
#line 285 "dgehrd.f"
	    }
#line 286 "dgehrd.f"
	}
#line 287 "dgehrd.f"
    }
#line 288 "dgehrd.f"
    ldwork = *n;

#line 290 "dgehrd.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

#line 294 "dgehrd.f"
	i__ = *ilo;

#line 296 "dgehrd.f"
    } else {

/*        Use blocked code */

#line 300 "dgehrd.f"
	i__1 = *ihi - 1 - nx;
#line 300 "dgehrd.f"
	i__2 = nb;
#line 300 "dgehrd.f"
	for (i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 301 "dgehrd.f"
	    i__3 = nb, i__4 = *ihi - i__;
#line 301 "dgehrd.f"
	    ib = min(i__3,i__4);

/*           Reduce columns i:i+ib-1 to Hessenberg form, returning the */
/*           matrices V and T of the block reflector H = I - V*T*V**T */
/*           which performs the reduction, and also the matrix Y = A*V*T */

#line 307 "dgehrd.f"
	    dlahr2_(ihi, &i__, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], t, &
		    c__65, &work[1], &ldwork);

/*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the */
/*           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set */
/*           to 1 */

#line 314 "dgehrd.f"
	    ei = a[i__ + ib + (i__ + ib - 1) * a_dim1];
#line 315 "dgehrd.f"
	    a[i__ + ib + (i__ + ib - 1) * a_dim1] = 1.;
#line 316 "dgehrd.f"
	    i__3 = *ihi - i__ - ib + 1;
#line 316 "dgehrd.f"
	    dgemm_("No transpose", "Transpose", ihi, &i__3, &ib, &c_b25, &
		    work[1], &ldwork, &a[i__ + ib + i__ * a_dim1], lda, &
		    c_b26, &a[(i__ + ib) * a_dim1 + 1], lda, (ftnlen)12, (
		    ftnlen)9);
#line 320 "dgehrd.f"
	    a[i__ + ib + (i__ + ib - 1) * a_dim1] = ei;

/*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the */
/*           right */

#line 325 "dgehrd.f"
	    i__3 = ib - 1;
#line 325 "dgehrd.f"
	    dtrmm_("Right", "Lower", "Transpose", "Unit", &i__, &i__3, &c_b26,
		     &a[i__ + 1 + i__ * a_dim1], lda, &work[1], &ldwork, (
		    ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 328 "dgehrd.f"
	    i__3 = ib - 2;
#line 328 "dgehrd.f"
	    for (j = 0; j <= i__3; ++j) {
#line 329 "dgehrd.f"
		daxpy_(&i__, &c_b25, &work[ldwork * j + 1], &c__1, &a[(i__ + 
			j + 1) * a_dim1 + 1], &c__1);
#line 331 "dgehrd.f"
/* L30: */
#line 331 "dgehrd.f"
	    }

/*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the */
/*           left */

#line 336 "dgehrd.f"
	    i__3 = *ihi - i__;
#line 336 "dgehrd.f"
	    i__4 = *n - i__ - ib + 1;
#line 336 "dgehrd.f"
	    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
		    i__4, &ib, &a[i__ + 1 + i__ * a_dim1], lda, t, &c__65, &a[
		    i__ + 1 + (i__ + ib) * a_dim1], lda, &work[1], &ldwork, (
		    ftnlen)4, (ftnlen)9, (ftnlen)7, (ftnlen)10);
#line 340 "dgehrd.f"
/* L40: */
#line 340 "dgehrd.f"
	}
#line 341 "dgehrd.f"
    }

/*     Use unblocked code to reduce the rest of the matrix */

#line 345 "dgehrd.f"
    dgehd2_(n, &i__, ihi, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 346 "dgehrd.f"
    work[1] = (doublereal) iws;

#line 348 "dgehrd.f"
    return 0;

/*     End of DGEHRD */

} /* dgehrd_ */

