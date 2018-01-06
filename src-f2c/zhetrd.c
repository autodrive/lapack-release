#line 1 "zhetrd.f"
/* zhetrd.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static doublereal c_b23 = 1.;

/* > \brief \b ZHETRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRD reduces a complex Hermitian matrix A to real symmetric */
/* > tridiagonal form T by a unitary similarity transformation: */
/* > Q**H * A * Q = T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          tridiagonal matrix T, and the elements above the first */
/* >          superdiagonal, with the array TAU, represent the unitary */
/* >          matrix Q as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and first subdiagonal of A are over- */
/* >          written by the corresponding elements of the tridiagonal */
/* >          matrix T, and the elements below the first subdiagonal, with */
/* >          the array TAU, represent the unitary matrix Q as a product */
/* >          of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T: */
/* >          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= 1. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16HEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(n-1) . . . H(2) H(1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in */
/* >  A(1:i-1,i+1), and tau in TAU(i). */
/* > */
/* >  If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(n-1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), */
/* >  and tau in TAU(i). */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with n = 5: */
/* > */
/* >  if UPLO = 'U':                       if UPLO = 'L': */
/* > */
/* >    (  d   e   v2  v3  v4 )              (  d                  ) */
/* >    (      d   e   v3  v4 )              (  e   d              ) */
/* >    (          d   e   v4 )              (  v1  e   d          ) */
/* >    (              d   e  )              (  v1  v2  e   d      ) */
/* >    (                  d  )              (  v1  v2  v3  e   d  ) */
/* > */
/* >  where d and e denote diagonal and off-diagonal elements of T, and vi */
/* >  denotes an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhetrd_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, 
	doublecomplex *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, nb, kk, nx, iws;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    static logical upper;
    extern /* Subroutine */ int zhetd2_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     ftnlen), zher2k_(char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlatrd_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen);
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 237 "zhetrd.f"
    /* Parameter adjustments */
#line 237 "zhetrd.f"
    a_dim1 = *lda;
#line 237 "zhetrd.f"
    a_offset = 1 + a_dim1;
#line 237 "zhetrd.f"
    a -= a_offset;
#line 237 "zhetrd.f"
    --d__;
#line 237 "zhetrd.f"
    --e;
#line 237 "zhetrd.f"
    --tau;
#line 237 "zhetrd.f"
    --work;
#line 237 "zhetrd.f"

#line 237 "zhetrd.f"
    /* Function Body */
#line 237 "zhetrd.f"
    *info = 0;
#line 238 "zhetrd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 239 "zhetrd.f"
    lquery = *lwork == -1;
#line 240 "zhetrd.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 241 "zhetrd.f"
	*info = -1;
#line 242 "zhetrd.f"
    } else if (*n < 0) {
#line 243 "zhetrd.f"
	*info = -2;
#line 244 "zhetrd.f"
    } else if (*lda < max(1,*n)) {
#line 245 "zhetrd.f"
	*info = -4;
#line 246 "zhetrd.f"
    } else if (*lwork < 1 && ! lquery) {
#line 247 "zhetrd.f"
	*info = -9;
#line 248 "zhetrd.f"
    }

#line 250 "zhetrd.f"
    if (*info == 0) {

/*        Determine the block size. */

#line 254 "zhetrd.f"
	nb = ilaenv_(&c__1, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 255 "zhetrd.f"
	lwkopt = *n * nb;
#line 256 "zhetrd.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 257 "zhetrd.f"
    }

#line 259 "zhetrd.f"
    if (*info != 0) {
#line 260 "zhetrd.f"
	i__1 = -(*info);
#line 260 "zhetrd.f"
	xerbla_("ZHETRD", &i__1, (ftnlen)6);
#line 261 "zhetrd.f"
	return 0;
#line 262 "zhetrd.f"
    } else if (lquery) {
#line 263 "zhetrd.f"
	return 0;
#line 264 "zhetrd.f"
    }

/*     Quick return if possible */

#line 268 "zhetrd.f"
    if (*n == 0) {
#line 269 "zhetrd.f"
	work[1].r = 1., work[1].i = 0.;
#line 270 "zhetrd.f"
	return 0;
#line 271 "zhetrd.f"
    }

#line 273 "zhetrd.f"
    nx = *n;
#line 274 "zhetrd.f"
    iws = 1;
#line 275 "zhetrd.f"
    if (nb > 1 && nb < *n) {

/*        Determine when to cross over from blocked to unblocked code */
/*        (last block is always handled by unblocked code). */

/* Computing MAX */
#line 280 "zhetrd.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "ZHETRD", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
#line 280 "zhetrd.f"
	nx = max(i__1,i__2);
#line 281 "zhetrd.f"
	if (nx < *n) {

/*           Determine if workspace is large enough for blocked code. */

#line 285 "zhetrd.f"
	    ldwork = *n;
#line 286 "zhetrd.f"
	    iws = ldwork * nb;
#line 287 "zhetrd.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code by setting NX = N. */

/* Computing MAX */
#line 293 "zhetrd.f"
		i__1 = *lwork / ldwork;
#line 293 "zhetrd.f"
		nb = max(i__1,1);
#line 294 "zhetrd.f"
		nbmin = ilaenv_(&c__2, "ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 295 "zhetrd.f"
		if (nb < nbmin) {
#line 295 "zhetrd.f"
		    nx = *n;
#line 295 "zhetrd.f"
		}
#line 297 "zhetrd.f"
	    }
#line 298 "zhetrd.f"
	} else {
#line 299 "zhetrd.f"
	    nx = *n;
#line 300 "zhetrd.f"
	}
#line 301 "zhetrd.f"
    } else {
#line 302 "zhetrd.f"
	nb = 1;
#line 303 "zhetrd.f"
    }

#line 305 "zhetrd.f"
    if (upper) {

/*        Reduce the upper triangle of A. */
/*        Columns 1:kk are handled by the unblocked method. */

#line 310 "zhetrd.f"
	kk = *n - (*n - nx + nb - 1) / nb * nb;
#line 311 "zhetrd.f"
	i__1 = kk + 1;
#line 311 "zhetrd.f"
	i__2 = -nb;
#line 311 "zhetrd.f"
	for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form the */
/*           matrix W which is needed to update the unreduced part of */
/*           the matrix */

#line 317 "zhetrd.f"
	    i__3 = i__ + nb - 1;
#line 317 "zhetrd.f"
	    zlatrd_(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
		    work[1], &ldwork, (ftnlen)1);

/*           Update the unreduced submatrix A(1:i-1,1:i-1), using an */
/*           update of the form:  A := A - V*W**H - W*V**H */

#line 323 "zhetrd.f"
	    i__3 = i__ - 1;
#line 323 "zhetrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 323 "zhetrd.f"
	    zher2k_(uplo, "No transpose", &i__3, &nb, &z__1, &a[i__ * a_dim1 
		    + 1], lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda, (
		    ftnlen)1, (ftnlen)12);

/*           Copy superdiagonal elements back into A, and diagonal */
/*           elements into D */

#line 329 "zhetrd.f"
	    i__3 = i__ + nb - 1;
#line 329 "zhetrd.f"
	    for (j = i__; j <= i__3; ++j) {
#line 330 "zhetrd.f"
		i__4 = j - 1 + j * a_dim1;
#line 330 "zhetrd.f"
		i__5 = j - 1;
#line 330 "zhetrd.f"
		a[i__4].r = e[i__5], a[i__4].i = 0.;
#line 331 "zhetrd.f"
		i__4 = j;
#line 331 "zhetrd.f"
		i__5 = j + j * a_dim1;
#line 331 "zhetrd.f"
		d__[i__4] = a[i__5].r;
#line 332 "zhetrd.f"
/* L10: */
#line 332 "zhetrd.f"
	    }
#line 333 "zhetrd.f"
/* L20: */
#line 333 "zhetrd.f"
	}

/*        Use unblocked code to reduce the last or only block */

#line 337 "zhetrd.f"
	zhetd2_(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo,
		 (ftnlen)1);
#line 338 "zhetrd.f"
    } else {

/*        Reduce the lower triangle of A */

#line 342 "zhetrd.f"
	i__2 = *n - nx;
#line 342 "zhetrd.f"
	i__1 = nb;
#line 342 "zhetrd.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form the */
/*           matrix W which is needed to update the unreduced part of */
/*           the matrix */

#line 348 "zhetrd.f"
	    i__3 = *n - i__ + 1;
#line 348 "zhetrd.f"
	    zlatrd_(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &
		    tau[i__], &work[1], &ldwork, (ftnlen)1);

/*           Update the unreduced submatrix A(i+nb:n,i+nb:n), using */
/*           an update of the form:  A := A - V*W**H - W*V**H */

#line 354 "zhetrd.f"
	    i__3 = *n - i__ - nb + 1;
#line 354 "zhetrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 354 "zhetrd.f"
	    zher2k_(uplo, "No transpose", &i__3, &nb, &z__1, &a[i__ + nb + 
		    i__ * a_dim1], lda, &work[nb + 1], &ldwork, &c_b23, &a[
		    i__ + nb + (i__ + nb) * a_dim1], lda, (ftnlen)1, (ftnlen)
		    12);

/*           Copy subdiagonal elements back into A, and diagonal */
/*           elements into D */

#line 361 "zhetrd.f"
	    i__3 = i__ + nb - 1;
#line 361 "zhetrd.f"
	    for (j = i__; j <= i__3; ++j) {
#line 362 "zhetrd.f"
		i__4 = j + 1 + j * a_dim1;
#line 362 "zhetrd.f"
		i__5 = j;
#line 362 "zhetrd.f"
		a[i__4].r = e[i__5], a[i__4].i = 0.;
#line 363 "zhetrd.f"
		i__4 = j;
#line 363 "zhetrd.f"
		i__5 = j + j * a_dim1;
#line 363 "zhetrd.f"
		d__[i__4] = a[i__5].r;
#line 364 "zhetrd.f"
/* L30: */
#line 364 "zhetrd.f"
	    }
#line 365 "zhetrd.f"
/* L40: */
#line 365 "zhetrd.f"
	}

/*        Use unblocked code to reduce the last or only block */

#line 369 "zhetrd.f"
	i__1 = *n - i__ + 1;
#line 369 "zhetrd.f"
	zhetd2_(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], 
		&tau[i__], &iinfo, (ftnlen)1);
#line 371 "zhetrd.f"
    }

#line 373 "zhetrd.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 374 "zhetrd.f"
    return 0;

/*     End of ZHETRD */

} /* zhetrd_ */

