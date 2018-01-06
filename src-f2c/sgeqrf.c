#line 1 "sgeqrf.f"
/* sgeqrf.f -- translated by f2c (version 20100827).
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

#line 1 "sgeqrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SGEQRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEQRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQRF computes a QR factorization of a real M-by-N matrix A: */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          TAU is REAL array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup realGEcomputational */

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
/* Subroutine */ int sgeqrf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), slarfb_(char *,
	     char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;


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

#line 172 "sgeqrf.f"
    /* Parameter adjustments */
#line 172 "sgeqrf.f"
    a_dim1 = *lda;
#line 172 "sgeqrf.f"
    a_offset = 1 + a_dim1;
#line 172 "sgeqrf.f"
    a -= a_offset;
#line 172 "sgeqrf.f"
    --tau;
#line 172 "sgeqrf.f"
    --work;
#line 172 "sgeqrf.f"

#line 172 "sgeqrf.f"
    /* Function Body */
#line 172 "sgeqrf.f"
    *info = 0;
#line 173 "sgeqrf.f"
    nb = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 174 "sgeqrf.f"
    lwkopt = *n * nb;
#line 175 "sgeqrf.f"
    work[1] = (doublereal) lwkopt;
#line 176 "sgeqrf.f"
    lquery = *lwork == -1;
#line 177 "sgeqrf.f"
    if (*m < 0) {
#line 178 "sgeqrf.f"
	*info = -1;
#line 179 "sgeqrf.f"
    } else if (*n < 0) {
#line 180 "sgeqrf.f"
	*info = -2;
#line 181 "sgeqrf.f"
    } else if (*lda < max(1,*m)) {
#line 182 "sgeqrf.f"
	*info = -4;
#line 183 "sgeqrf.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 184 "sgeqrf.f"
	*info = -7;
#line 185 "sgeqrf.f"
    }
#line 186 "sgeqrf.f"
    if (*info != 0) {
#line 187 "sgeqrf.f"
	i__1 = -(*info);
#line 187 "sgeqrf.f"
	xerbla_("SGEQRF", &i__1, (ftnlen)6);
#line 188 "sgeqrf.f"
	return 0;
#line 189 "sgeqrf.f"
    } else if (lquery) {
#line 190 "sgeqrf.f"
	return 0;
#line 191 "sgeqrf.f"
    }

/*     Quick return if possible */

#line 195 "sgeqrf.f"
    k = min(*m,*n);
#line 196 "sgeqrf.f"
    if (k == 0) {
#line 197 "sgeqrf.f"
	work[1] = 1.;
#line 198 "sgeqrf.f"
	return 0;
#line 199 "sgeqrf.f"
    }

#line 201 "sgeqrf.f"
    nbmin = 2;
#line 202 "sgeqrf.f"
    nx = 0;
#line 203 "sgeqrf.f"
    iws = *n;
#line 204 "sgeqrf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 208 "sgeqrf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 208 "sgeqrf.f"
	nx = max(i__1,i__2);
#line 209 "sgeqrf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 213 "sgeqrf.f"
	    ldwork = *n;
#line 214 "sgeqrf.f"
	    iws = ldwork * nb;
#line 215 "sgeqrf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 220 "sgeqrf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 221 "sgeqrf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 221 "sgeqrf.f"
		nbmin = max(i__1,i__2);
#line 223 "sgeqrf.f"
	    }
#line 224 "sgeqrf.f"
	}
#line 225 "sgeqrf.f"
    }

#line 227 "sgeqrf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

#line 231 "sgeqrf.f"
	i__1 = k - nx;
#line 231 "sgeqrf.f"
	i__2 = nb;
#line 231 "sgeqrf.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 232 "sgeqrf.f"
	    i__3 = k - i__ + 1;
#line 232 "sgeqrf.f"
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block */
/*           A(i:m,i:i+ib-1) */

#line 237 "sgeqrf.f"
	    i__3 = *m - i__ + 1;
#line 237 "sgeqrf.f"
	    sgeqr2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
#line 239 "sgeqrf.f"
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 244 "sgeqrf.f"
		i__3 = *m - i__ + 1;
#line 244 "sgeqrf.f"
		slarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H**T to A(i:m,i+ib:n) from the left */

#line 249 "sgeqrf.f"
		i__3 = *m - i__ + 1;
#line 249 "sgeqrf.f"
		i__4 = *n - i__ - ib + 1;
#line 249 "sgeqrf.f"
		slarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib 
			+ 1], &ldwork, (ftnlen)4, (ftnlen)9, (ftnlen)7, (
			ftnlen)10);
#line 253 "sgeqrf.f"
	    }
#line 254 "sgeqrf.f"
/* L10: */
#line 254 "sgeqrf.f"
	}
#line 255 "sgeqrf.f"
    } else {
#line 256 "sgeqrf.f"
	i__ = 1;
#line 257 "sgeqrf.f"
    }

/*     Use unblocked code to factor the last or only block. */

#line 261 "sgeqrf.f"
    if (i__ <= k) {
#line 261 "sgeqrf.f"
	i__2 = *m - i__ + 1;
#line 261 "sgeqrf.f"
	i__1 = *n - i__ + 1;
#line 261 "sgeqrf.f"
	sgeqr2_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
#line 261 "sgeqrf.f"
    }

#line 265 "sgeqrf.f"
    work[1] = (doublereal) iws;
#line 266 "sgeqrf.f"
    return 0;

/*     End of SGEQRF */

} /* sgeqrf_ */

