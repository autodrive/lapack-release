#line 1 "dgeqrfp.f"
/* dgeqrfp.f -- translated by f2c (version 20100827).
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

#line 1 "dgeqrfp.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEQRFP computes a QR factorization of a real M-by-N matrix A: */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n); the elements below the diagonal, */
/* >          with the array TAU, represent the orthogonal matrix Q as a */
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
/* >          TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup doubleGEcomputational */

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
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/* >  and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgeqrfp_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarft_(char *, char *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), xerbla_(char *, integer 
	    *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int dgeqr2p_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);


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

#line 172 "dgeqrfp.f"
    /* Parameter adjustments */
#line 172 "dgeqrfp.f"
    a_dim1 = *lda;
#line 172 "dgeqrfp.f"
    a_offset = 1 + a_dim1;
#line 172 "dgeqrfp.f"
    a -= a_offset;
#line 172 "dgeqrfp.f"
    --tau;
#line 172 "dgeqrfp.f"
    --work;
#line 172 "dgeqrfp.f"

#line 172 "dgeqrfp.f"
    /* Function Body */
#line 172 "dgeqrfp.f"
    *info = 0;
#line 173 "dgeqrfp.f"
    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 174 "dgeqrfp.f"
    lwkopt = *n * nb;
#line 175 "dgeqrfp.f"
    work[1] = (doublereal) lwkopt;
#line 176 "dgeqrfp.f"
    lquery = *lwork == -1;
#line 177 "dgeqrfp.f"
    if (*m < 0) {
#line 178 "dgeqrfp.f"
	*info = -1;
#line 179 "dgeqrfp.f"
    } else if (*n < 0) {
#line 180 "dgeqrfp.f"
	*info = -2;
#line 181 "dgeqrfp.f"
    } else if (*lda < max(1,*m)) {
#line 182 "dgeqrfp.f"
	*info = -4;
#line 183 "dgeqrfp.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 184 "dgeqrfp.f"
	*info = -7;
#line 185 "dgeqrfp.f"
    }
#line 186 "dgeqrfp.f"
    if (*info != 0) {
#line 187 "dgeqrfp.f"
	i__1 = -(*info);
#line 187 "dgeqrfp.f"
	xerbla_("DGEQRFP", &i__1, (ftnlen)7);
#line 188 "dgeqrfp.f"
	return 0;
#line 189 "dgeqrfp.f"
    } else if (lquery) {
#line 190 "dgeqrfp.f"
	return 0;
#line 191 "dgeqrfp.f"
    }

/*     Quick return if possible */

#line 195 "dgeqrfp.f"
    k = min(*m,*n);
#line 196 "dgeqrfp.f"
    if (k == 0) {
#line 197 "dgeqrfp.f"
	work[1] = 1.;
#line 198 "dgeqrfp.f"
	return 0;
#line 199 "dgeqrfp.f"
    }

#line 201 "dgeqrfp.f"
    nbmin = 2;
#line 202 "dgeqrfp.f"
    nx = 0;
#line 203 "dgeqrfp.f"
    iws = *n;
#line 204 "dgeqrfp.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 208 "dgeqrfp.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 208 "dgeqrfp.f"
	nx = max(i__1,i__2);
#line 209 "dgeqrfp.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 213 "dgeqrfp.f"
	    ldwork = *n;
#line 214 "dgeqrfp.f"
	    iws = ldwork * nb;
#line 215 "dgeqrfp.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 220 "dgeqrfp.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 221 "dgeqrfp.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 221 "dgeqrfp.f"
		nbmin = max(i__1,i__2);
#line 223 "dgeqrfp.f"
	    }
#line 224 "dgeqrfp.f"
	}
#line 225 "dgeqrfp.f"
    }

#line 227 "dgeqrfp.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

#line 231 "dgeqrfp.f"
	i__1 = k - nx;
#line 231 "dgeqrfp.f"
	i__2 = nb;
#line 231 "dgeqrfp.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 232 "dgeqrfp.f"
	    i__3 = k - i__ + 1;
#line 232 "dgeqrfp.f"
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

#line 237 "dgeqrfp.f"
	    i__3 = *m - i__ + 1;
#line 237 "dgeqrfp.f"
	    dgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
#line 239 "dgeqrfp.f"
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 244 "dgeqrfp.f"
		i__3 = *m - i__ + 1;
#line 244 "dgeqrfp.f"
		dlarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**T to A(i:m,i+ib:n) from the left */

#line 249 "dgeqrfp.f"
		i__3 = *m - i__ + 1;
#line 249 "dgeqrfp.f"
		i__4 = *n - i__ - ib + 1;
#line 249 "dgeqrfp.f"
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib 
			+ 1], &ldwork, (ftnlen)4, (ftnlen)9, (ftnlen)7, (
			ftnlen)10);
#line 253 "dgeqrfp.f"
	    }
#line 254 "dgeqrfp.f"
/* L10: */
#line 254 "dgeqrfp.f"
	}
#line 255 "dgeqrfp.f"
    } else {
#line 256 "dgeqrfp.f"
	i__ = 1;
#line 257 "dgeqrfp.f"
    }

/*     Use unblocked code to factor the last or only block. */

#line 261 "dgeqrfp.f"
    if (i__ <= k) {
#line 261 "dgeqrfp.f"
	i__2 = *m - i__ + 1;
#line 261 "dgeqrfp.f"
	i__1 = *n - i__ + 1;
#line 261 "dgeqrfp.f"
	dgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
#line 261 "dgeqrfp.f"
    }

#line 265 "dgeqrfp.f"
    work[1] = (doublereal) iws;
#line 266 "dgeqrfp.f"
    return 0;

/*     End of DGEQRFP */

} /* dgeqrfp_ */

