#line 1 "zgeqrfp.f"
/* zgeqrfp.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqrfp.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b ZGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQRFP computes a QR factorization of a complex M-by-N matrix A: */
/* > A = Q * R. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n); the elements below the diagonal, */
/* >          with the array TAU, represent the unitary matrix Q as a */
/* >          product of min(m,n) elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (min(M,N)) */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
/* >          For optimum performance LWORK >= N*NB, where NB is */
/* >          the optimal blocksize. */
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

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqrfp_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int zlarft_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zgeqr2p_(integer *, integer *, doublecomplex *
	    , integer *, doublecomplex *, doublecomplex *, integer *);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 172 "zgeqrfp.f"
    /* Parameter adjustments */
#line 172 "zgeqrfp.f"
    a_dim1 = *lda;
#line 172 "zgeqrfp.f"
    a_offset = 1 + a_dim1;
#line 172 "zgeqrfp.f"
    a -= a_offset;
#line 172 "zgeqrfp.f"
    --tau;
#line 172 "zgeqrfp.f"
    --work;
#line 172 "zgeqrfp.f"

#line 172 "zgeqrfp.f"
    /* Function Body */
#line 172 "zgeqrfp.f"
    *info = 0;
#line 173 "zgeqrfp.f"
    nb = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 174 "zgeqrfp.f"
    lwkopt = *n * nb;
#line 175 "zgeqrfp.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 176 "zgeqrfp.f"
    lquery = *lwork == -1;
#line 177 "zgeqrfp.f"
    if (*m < 0) {
#line 178 "zgeqrfp.f"
	*info = -1;
#line 179 "zgeqrfp.f"
    } else if (*n < 0) {
#line 180 "zgeqrfp.f"
	*info = -2;
#line 181 "zgeqrfp.f"
    } else if (*lda < max(1,*m)) {
#line 182 "zgeqrfp.f"
	*info = -4;
#line 183 "zgeqrfp.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 184 "zgeqrfp.f"
	*info = -7;
#line 185 "zgeqrfp.f"
    }
#line 186 "zgeqrfp.f"
    if (*info != 0) {
#line 187 "zgeqrfp.f"
	i__1 = -(*info);
#line 187 "zgeqrfp.f"
	xerbla_("ZGEQRFP", &i__1, (ftnlen)7);
#line 188 "zgeqrfp.f"
	return 0;
#line 189 "zgeqrfp.f"
    } else if (lquery) {
#line 190 "zgeqrfp.f"
	return 0;
#line 191 "zgeqrfp.f"
    }

/*     Quick return if possible */

#line 195 "zgeqrfp.f"
    k = min(*m,*n);
#line 196 "zgeqrfp.f"
    if (k == 0) {
#line 197 "zgeqrfp.f"
	work[1].r = 1., work[1].i = 0.;
#line 198 "zgeqrfp.f"
	return 0;
#line 199 "zgeqrfp.f"
    }

#line 201 "zgeqrfp.f"
    nbmin = 2;
#line 202 "zgeqrfp.f"
    nx = 0;
#line 203 "zgeqrfp.f"
    iws = *n;
#line 204 "zgeqrfp.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 208 "zgeqrfp.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 208 "zgeqrfp.f"
	nx = max(i__1,i__2);
#line 209 "zgeqrfp.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 213 "zgeqrfp.f"
	    ldwork = *n;
#line 214 "zgeqrfp.f"
	    iws = ldwork * nb;
#line 215 "zgeqrfp.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 220 "zgeqrfp.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 221 "zgeqrfp.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 221 "zgeqrfp.f"
		nbmin = max(i__1,i__2);
#line 223 "zgeqrfp.f"
	    }
#line 224 "zgeqrfp.f"
	}
#line 225 "zgeqrfp.f"
    }

#line 227 "zgeqrfp.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

#line 231 "zgeqrfp.f"
	i__1 = k - nx;
#line 231 "zgeqrfp.f"
	i__2 = nb;
#line 231 "zgeqrfp.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 232 "zgeqrfp.f"
	    i__3 = k - i__ + 1;
#line 232 "zgeqrfp.f"
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

#line 237 "zgeqrfp.f"
	    i__3 = *m - i__ + 1;
#line 237 "zgeqrfp.f"
	    zgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
#line 239 "zgeqrfp.f"
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 244 "zgeqrfp.f"
		i__3 = *m - i__ + 1;
#line 244 "zgeqrfp.f"
		zlarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**H to A(i:m,i+ib:n) from the left */

#line 249 "zgeqrfp.f"
		i__3 = *m - i__ + 1;
#line 249 "zgeqrfp.f"
		i__4 = *n - i__ - ib + 1;
#line 249 "zgeqrfp.f"
		zlarfb_("Left", "Conjugate transpose", "Forward", "Columnwise"
			, &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &
			work[1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, 
			&work[ib + 1], &ldwork, (ftnlen)4, (ftnlen)19, (
			ftnlen)7, (ftnlen)10);
#line 253 "zgeqrfp.f"
	    }
#line 254 "zgeqrfp.f"
/* L10: */
#line 254 "zgeqrfp.f"
	}
#line 255 "zgeqrfp.f"
    } else {
#line 256 "zgeqrfp.f"
	i__ = 1;
#line 257 "zgeqrfp.f"
    }

/*     Use unblocked code to factor the last or only block. */

#line 261 "zgeqrfp.f"
    if (i__ <= k) {
#line 261 "zgeqrfp.f"
	i__2 = *m - i__ + 1;
#line 261 "zgeqrfp.f"
	i__1 = *n - i__ + 1;
#line 261 "zgeqrfp.f"
	zgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
#line 261 "zgeqrfp.f"
    }

#line 265 "zgeqrfp.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 266 "zgeqrfp.f"
    return 0;

/*     End of ZGEQRFP */

} /* zgeqrfp_ */

