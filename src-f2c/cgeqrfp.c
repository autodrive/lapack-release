#line 1 "cgeqrfp.f"
/* cgeqrfp.f -- translated by f2c (version 20100827).
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

#line 1 "cgeqrfp.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b CGEQRFP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrfp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrfp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrfp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQRFP computes a QR factorization of a complex M-by-N matrix A: */
/* > A = Q * R. The diagonal entries of R are real and nonnegative. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
/* >          upper triangular if m >= n). The diagonal entries of R */
/* >          are real and nonnegative; the elements below the diagonal, */
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
/* >          TAU is COMPLEX array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

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
/* > */
/* > See Lapack Working Note 203 for details */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgeqrfp_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    clarft_(char *, char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int cgeqr2p_(integer *, integer *, doublecomplex *
	    , integer *, doublecomplex *, doublecomplex *, integer *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 175 "cgeqrfp.f"
    /* Parameter adjustments */
#line 175 "cgeqrfp.f"
    a_dim1 = *lda;
#line 175 "cgeqrfp.f"
    a_offset = 1 + a_dim1;
#line 175 "cgeqrfp.f"
    a -= a_offset;
#line 175 "cgeqrfp.f"
    --tau;
#line 175 "cgeqrfp.f"
    --work;
#line 175 "cgeqrfp.f"

#line 175 "cgeqrfp.f"
    /* Function Body */
#line 175 "cgeqrfp.f"
    *info = 0;
#line 176 "cgeqrfp.f"
    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 177 "cgeqrfp.f"
    lwkopt = *n * nb;
#line 178 "cgeqrfp.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 179 "cgeqrfp.f"
    lquery = *lwork == -1;
#line 180 "cgeqrfp.f"
    if (*m < 0) {
#line 181 "cgeqrfp.f"
	*info = -1;
#line 182 "cgeqrfp.f"
    } else if (*n < 0) {
#line 183 "cgeqrfp.f"
	*info = -2;
#line 184 "cgeqrfp.f"
    } else if (*lda < max(1,*m)) {
#line 185 "cgeqrfp.f"
	*info = -4;
#line 186 "cgeqrfp.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 187 "cgeqrfp.f"
	*info = -7;
#line 188 "cgeqrfp.f"
    }
#line 189 "cgeqrfp.f"
    if (*info != 0) {
#line 190 "cgeqrfp.f"
	i__1 = -(*info);
#line 190 "cgeqrfp.f"
	xerbla_("CGEQRFP", &i__1, (ftnlen)7);
#line 191 "cgeqrfp.f"
	return 0;
#line 192 "cgeqrfp.f"
    } else if (lquery) {
#line 193 "cgeqrfp.f"
	return 0;
#line 194 "cgeqrfp.f"
    }

/*     Quick return if possible */

#line 198 "cgeqrfp.f"
    k = min(*m,*n);
#line 199 "cgeqrfp.f"
    if (k == 0) {
#line 200 "cgeqrfp.f"
	work[1].r = 1., work[1].i = 0.;
#line 201 "cgeqrfp.f"
	return 0;
#line 202 "cgeqrfp.f"
    }

#line 204 "cgeqrfp.f"
    nbmin = 2;
#line 205 "cgeqrfp.f"
    nx = 0;
#line 206 "cgeqrfp.f"
    iws = *n;
#line 207 "cgeqrfp.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 211 "cgeqrfp.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "CGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 211 "cgeqrfp.f"
	nx = max(i__1,i__2);
#line 212 "cgeqrfp.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 216 "cgeqrfp.f"
	    ldwork = *n;
#line 217 "cgeqrfp.f"
	    iws = ldwork * nb;
#line 218 "cgeqrfp.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 223 "cgeqrfp.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 224 "cgeqrfp.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "CGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 224 "cgeqrfp.f"
		nbmin = max(i__1,i__2);
#line 226 "cgeqrfp.f"
	    }
#line 227 "cgeqrfp.f"
	}
#line 228 "cgeqrfp.f"
    }

#line 230 "cgeqrfp.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

#line 234 "cgeqrfp.f"
	i__1 = k - nx;
#line 234 "cgeqrfp.f"
	i__2 = nb;
#line 234 "cgeqrfp.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 235 "cgeqrfp.f"
	    i__3 = k - i__ + 1;
#line 235 "cgeqrfp.f"
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

#line 240 "cgeqrfp.f"
	    i__3 = *m - i__ + 1;
#line 240 "cgeqrfp.f"
	    cgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
#line 242 "cgeqrfp.f"
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 247 "cgeqrfp.f"
		i__3 = *m - i__ + 1;
#line 247 "cgeqrfp.f"
		clarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**H to A(i:m,i+ib:n) from the left */

#line 252 "cgeqrfp.f"
		i__3 = *m - i__ + 1;
#line 252 "cgeqrfp.f"
		i__4 = *n - i__ - ib + 1;
#line 252 "cgeqrfp.f"
		clarfb_("Left", "Conjugate transpose", "Forward", "Columnwise"
			, &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &
			work[1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, 
			&work[ib + 1], &ldwork, (ftnlen)4, (ftnlen)19, (
			ftnlen)7, (ftnlen)10);
#line 256 "cgeqrfp.f"
	    }
#line 257 "cgeqrfp.f"
/* L10: */
#line 257 "cgeqrfp.f"
	}
#line 258 "cgeqrfp.f"
    } else {
#line 259 "cgeqrfp.f"
	i__ = 1;
#line 260 "cgeqrfp.f"
    }

/*     Use unblocked code to factor the last or only block. */

#line 264 "cgeqrfp.f"
    if (i__ <= k) {
#line 264 "cgeqrfp.f"
	i__2 = *m - i__ + 1;
#line 264 "cgeqrfp.f"
	i__1 = *n - i__ + 1;
#line 264 "cgeqrfp.f"
	cgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		1], &iinfo);
#line 264 "cgeqrfp.f"
    }

#line 268 "cgeqrfp.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 269 "cgeqrfp.f"
    return 0;

/*     End of CGEQRFP */

} /* cgeqrfp_ */

