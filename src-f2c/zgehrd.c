#line 1 "zgehrd.f"
/* zgehrd.f -- translated by f2c (version 20100827).
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

#line 1 "zgehrd.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b ZGEHRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16        A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEHRD reduces a complex general matrix A to upper Hessenberg form H by */
/* > an unitary similarity transformation:  Q**H * A * Q = H . */
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
/* >          set by a previous call to ZGEBAL; otherwise they should be */
/* >          set to 1 and N respectively. See Further Details. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the N-by-N general matrix to be reduced. */
/* >          On exit, the upper triangle and the first subdiagonal of A */
/* >          are overwritten with the upper Hessenberg matrix H, and the */
/* >          elements below the first subdiagonal, with the array TAU, */
/* >          represent the unitary matrix Q as a product of elementary */
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
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to */
/* >          zero. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
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

/* > \ingroup complex16GEcomputational */

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
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
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
/* Subroutine */ int zgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublecomplex t[4160]	/* was [65][64] */;
    static integer ib;
    static doublecomplex ei;
    static integer nb, nh, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), ztrmm_(char *, char *, char *, char *,
	     integer *, integer *, doublecomplex *, doublecomplex *, integer *
	    , doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zgehd2_(integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *), zlahr2_(integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 216 "zgehrd.f"
    /* Parameter adjustments */
#line 216 "zgehrd.f"
    a_dim1 = *lda;
#line 216 "zgehrd.f"
    a_offset = 1 + a_dim1;
#line 216 "zgehrd.f"
    a -= a_offset;
#line 216 "zgehrd.f"
    --tau;
#line 216 "zgehrd.f"
    --work;
#line 216 "zgehrd.f"

#line 216 "zgehrd.f"
    /* Function Body */
#line 216 "zgehrd.f"
    *info = 0;
/* Computing MIN */
#line 217 "zgehrd.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "ZGEHRD", " ", n, ilo, ihi, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 217 "zgehrd.f"
    nb = min(i__1,i__2);
#line 218 "zgehrd.f"
    lwkopt = *n * nb;
#line 219 "zgehrd.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 220 "zgehrd.f"
    lquery = *lwork == -1;
#line 221 "zgehrd.f"
    if (*n < 0) {
#line 222 "zgehrd.f"
	*info = -1;
#line 223 "zgehrd.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 224 "zgehrd.f"
	*info = -2;
#line 225 "zgehrd.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 226 "zgehrd.f"
	*info = -3;
#line 227 "zgehrd.f"
    } else if (*lda < max(1,*n)) {
#line 228 "zgehrd.f"
	*info = -5;
#line 229 "zgehrd.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 230 "zgehrd.f"
	*info = -8;
#line 231 "zgehrd.f"
    }
#line 232 "zgehrd.f"
    if (*info != 0) {
#line 233 "zgehrd.f"
	i__1 = -(*info);
#line 233 "zgehrd.f"
	xerbla_("ZGEHRD", &i__1, (ftnlen)6);
#line 234 "zgehrd.f"
	return 0;
#line 235 "zgehrd.f"
    } else if (lquery) {
#line 236 "zgehrd.f"
	return 0;
#line 237 "zgehrd.f"
    }

/*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */

#line 241 "zgehrd.f"
    i__1 = *ilo - 1;
#line 241 "zgehrd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "zgehrd.f"
	i__2 = i__;
#line 242 "zgehrd.f"
	tau[i__2].r = 0., tau[i__2].i = 0.;
#line 243 "zgehrd.f"
/* L10: */
#line 243 "zgehrd.f"
    }
#line 244 "zgehrd.f"
    i__1 = *n - 1;
#line 244 "zgehrd.f"
    for (i__ = max(1,*ihi); i__ <= i__1; ++i__) {
#line 245 "zgehrd.f"
	i__2 = i__;
#line 245 "zgehrd.f"
	tau[i__2].r = 0., tau[i__2].i = 0.;
#line 246 "zgehrd.f"
/* L20: */
#line 246 "zgehrd.f"
    }

/*     Quick return if possible */

#line 250 "zgehrd.f"
    nh = *ihi - *ilo + 1;
#line 251 "zgehrd.f"
    if (nh <= 1) {
#line 252 "zgehrd.f"
	work[1].r = 1., work[1].i = 0.;
#line 253 "zgehrd.f"
	return 0;
#line 254 "zgehrd.f"
    }

/*     Determine the block size */

/* Computing MIN */
#line 258 "zgehrd.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "ZGEHRD", " ", n, ilo, ihi, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 258 "zgehrd.f"
    nb = min(i__1,i__2);
#line 259 "zgehrd.f"
    nbmin = 2;
#line 260 "zgehrd.f"
    iws = 1;
#line 261 "zgehrd.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to cross over from blocked to unblocked code */
/*        (last block is always handled by unblocked code) */

/* Computing MAX */
#line 266 "zgehrd.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "ZGEHRD", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 266 "zgehrd.f"
	nx = max(i__1,i__2);
#line 267 "zgehrd.f"
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code */

#line 271 "zgehrd.f"
	    iws = *n * nb;
#line 272 "zgehrd.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code */

/* Computing MAX */
#line 278 "zgehrd.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZGEHRD", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 278 "zgehrd.f"
		nbmin = max(i__1,i__2);
#line 280 "zgehrd.f"
		if (*lwork >= *n * nbmin) {
#line 281 "zgehrd.f"
		    nb = *lwork / *n;
#line 282 "zgehrd.f"
		} else {
#line 283 "zgehrd.f"
		    nb = 1;
#line 284 "zgehrd.f"
		}
#line 285 "zgehrd.f"
	    }
#line 286 "zgehrd.f"
	}
#line 287 "zgehrd.f"
    }
#line 288 "zgehrd.f"
    ldwork = *n;

#line 290 "zgehrd.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

#line 294 "zgehrd.f"
	i__ = *ilo;

#line 296 "zgehrd.f"
    } else {

/*        Use blocked code */

#line 300 "zgehrd.f"
	i__1 = *ihi - 1 - nx;
#line 300 "zgehrd.f"
	i__2 = nb;
#line 300 "zgehrd.f"
	for (i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 301 "zgehrd.f"
	    i__3 = nb, i__4 = *ihi - i__;
#line 301 "zgehrd.f"
	    ib = min(i__3,i__4);

/*           Reduce columns i:i+ib-1 to Hessenberg form, returning the */
/*           matrices V and T of the block reflector H = I - V*T*V**H */
/*           which performs the reduction, and also the matrix Y = A*V*T */

#line 307 "zgehrd.f"
	    zlahr2_(ihi, &i__, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], t, &
		    c__65, &work[1], &ldwork);

/*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the */
/*           right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set */
/*           to 1 */

#line 314 "zgehrd.f"
	    i__3 = i__ + ib + (i__ + ib - 1) * a_dim1;
#line 314 "zgehrd.f"
	    ei.r = a[i__3].r, ei.i = a[i__3].i;
#line 315 "zgehrd.f"
	    i__3 = i__ + ib + (i__ + ib - 1) * a_dim1;
#line 315 "zgehrd.f"
	    a[i__3].r = 1., a[i__3].i = 0.;
#line 316 "zgehrd.f"
	    i__3 = *ihi - i__ - ib + 1;
#line 316 "zgehrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 316 "zgehrd.f"
	    zgemm_("No transpose", "Conjugate transpose", ihi, &i__3, &ib, &
		    z__1, &work[1], &ldwork, &a[i__ + ib + i__ * a_dim1], lda,
		     &c_b2, &a[(i__ + ib) * a_dim1 + 1], lda, (ftnlen)12, (
		    ftnlen)19);
#line 320 "zgehrd.f"
	    i__3 = i__ + ib + (i__ + ib - 1) * a_dim1;
#line 320 "zgehrd.f"
	    a[i__3].r = ei.r, a[i__3].i = ei.i;

/*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the */
/*           right */

#line 325 "zgehrd.f"
	    i__3 = ib - 1;
#line 325 "zgehrd.f"
	    ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", &i__, &
		    i__3, &c_b2, &a[i__ + 1 + i__ * a_dim1], lda, &work[1], &
		    ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);
#line 328 "zgehrd.f"
	    i__3 = ib - 2;
#line 328 "zgehrd.f"
	    for (j = 0; j <= i__3; ++j) {
#line 329 "zgehrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 329 "zgehrd.f"
		zaxpy_(&i__, &z__1, &work[ldwork * j + 1], &c__1, &a[(i__ + j 
			+ 1) * a_dim1 + 1], &c__1);
#line 331 "zgehrd.f"
/* L30: */
#line 331 "zgehrd.f"
	    }

/*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the */
/*           left */

#line 336 "zgehrd.f"
	    i__3 = *ihi - i__;
#line 336 "zgehrd.f"
	    i__4 = *n - i__ - ib + 1;
#line 336 "zgehrd.f"
	    zlarfb_("Left", "Conjugate transpose", "Forward", "Columnwise", &
		    i__3, &i__4, &ib, &a[i__ + 1 + i__ * a_dim1], lda, t, &
		    c__65, &a[i__ + 1 + (i__ + ib) * a_dim1], lda, &work[1], &
		    ldwork, (ftnlen)4, (ftnlen)19, (ftnlen)7, (ftnlen)10);
#line 340 "zgehrd.f"
/* L40: */
#line 340 "zgehrd.f"
	}
#line 341 "zgehrd.f"
    }

/*     Use unblocked code to reduce the rest of the matrix */

#line 345 "zgehrd.f"
    zgehd2_(n, &i__, ihi, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 346 "zgehrd.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;

#line 348 "zgehrd.f"
    return 0;

/*     End of ZGEHRD */

} /* zgehrd_ */

