#line 1 "cgelqf.f"
/* cgelqf.f -- translated by f2c (version 20100827).
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

#line 1 "cgelqf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b CGELQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGELQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > CGELQF computes an LQ factorization of a complex M-by-N matrix A: */
/* > A = L * Q. */
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
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the m-by-min(m,n) lower trapezoidal matrix L (L is */
/* >          lower triangular if m <= n); the elements above the diagonal, */
/* >          with the array TAU, represent the unitary matrix Q as a */
/* >          product of elementary reflectors (see Further Details). */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= M*NB, where NB is the */
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

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(k)**H . . . H(2)**H H(1)**H, where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in */
/* >  A(i,i+1:n), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgelqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int cgelq2_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *), clarfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), clarft_(char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
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

#line 171 "cgelqf.f"
    /* Parameter adjustments */
#line 171 "cgelqf.f"
    a_dim1 = *lda;
#line 171 "cgelqf.f"
    a_offset = 1 + a_dim1;
#line 171 "cgelqf.f"
    a -= a_offset;
#line 171 "cgelqf.f"
    --tau;
#line 171 "cgelqf.f"
    --work;
#line 171 "cgelqf.f"

#line 171 "cgelqf.f"
    /* Function Body */
#line 171 "cgelqf.f"
    *info = 0;
#line 172 "cgelqf.f"
    nb = ilaenv_(&c__1, "CGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 173 "cgelqf.f"
    lwkopt = *m * nb;
#line 174 "cgelqf.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 175 "cgelqf.f"
    lquery = *lwork == -1;
#line 176 "cgelqf.f"
    if (*m < 0) {
#line 177 "cgelqf.f"
	*info = -1;
#line 178 "cgelqf.f"
    } else if (*n < 0) {
#line 179 "cgelqf.f"
	*info = -2;
#line 180 "cgelqf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "cgelqf.f"
	*info = -4;
#line 182 "cgelqf.f"
    } else if (*lwork < max(1,*m) && ! lquery) {
#line 183 "cgelqf.f"
	*info = -7;
#line 184 "cgelqf.f"
    }
#line 185 "cgelqf.f"
    if (*info != 0) {
#line 186 "cgelqf.f"
	i__1 = -(*info);
#line 186 "cgelqf.f"
	xerbla_("CGELQF", &i__1, (ftnlen)6);
#line 187 "cgelqf.f"
	return 0;
#line 188 "cgelqf.f"
    } else if (lquery) {
#line 189 "cgelqf.f"
	return 0;
#line 190 "cgelqf.f"
    }

/*     Quick return if possible */

#line 194 "cgelqf.f"
    k = min(*m,*n);
#line 195 "cgelqf.f"
    if (k == 0) {
#line 196 "cgelqf.f"
	work[1].r = 1., work[1].i = 0.;
#line 197 "cgelqf.f"
	return 0;
#line 198 "cgelqf.f"
    }

#line 200 "cgelqf.f"
    nbmin = 2;
#line 201 "cgelqf.f"
    nx = 0;
#line 202 "cgelqf.f"
    iws = *m;
#line 203 "cgelqf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 207 "cgelqf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "CGELQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 207 "cgelqf.f"
	nx = max(i__1,i__2);
#line 208 "cgelqf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 212 "cgelqf.f"
	    ldwork = *m;
#line 213 "cgelqf.f"
	    iws = ldwork * nb;
#line 214 "cgelqf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 219 "cgelqf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 220 "cgelqf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "CGELQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 220 "cgelqf.f"
		nbmin = max(i__1,i__2);
#line 222 "cgelqf.f"
	    }
#line 223 "cgelqf.f"
	}
#line 224 "cgelqf.f"
    }

#line 226 "cgelqf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

#line 230 "cgelqf.f"
	i__1 = k - nx;
#line 230 "cgelqf.f"
	i__2 = nb;
#line 230 "cgelqf.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 231 "cgelqf.f"
	    i__3 = k - i__ + 1;
#line 231 "cgelqf.f"
	    ib = min(i__3,nb);

/*           Compute the LQ factorization of the current block */
/*           A(i:i+ib-1,i:n) */

#line 236 "cgelqf.f"
	    i__3 = *n - i__ + 1;
#line 236 "cgelqf.f"
	    cgelq2_(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
#line 238 "cgelqf.f"
	    if (i__ + ib <= *m) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 243 "cgelqf.f"
		i__3 = *n - i__ + 1;
#line 243 "cgelqf.f"
		clarft_("Forward", "Rowwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)7);

/*              Apply H to A(i+ib:m,i:n) from the right */

#line 248 "cgelqf.f"
		i__3 = *m - i__ - ib + 1;
#line 248 "cgelqf.f"
		i__4 = *n - i__ + 1;
#line 248 "cgelqf.f"
		clarfb_("Right", "No transpose", "Forward", "Rowwise", &i__3, 
			&i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork, (ftnlen)5, (ftnlen)12, (ftnlen)7, (
			ftnlen)7);
#line 252 "cgelqf.f"
	    }
#line 253 "cgelqf.f"
/* L10: */
#line 253 "cgelqf.f"
	}
#line 254 "cgelqf.f"
    } else {
#line 255 "cgelqf.f"
	i__ = 1;
#line 256 "cgelqf.f"
    }

/*     Use unblocked code to factor the last or only block. */

#line 260 "cgelqf.f"
    if (i__ <= k) {
#line 260 "cgelqf.f"
	i__2 = *m - i__ + 1;
#line 260 "cgelqf.f"
	i__1 = *n - i__ + 1;
#line 260 "cgelqf.f"
	cgelq2_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
#line 260 "cgelqf.f"
    }

#line 264 "cgelqf.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 265 "cgelqf.f"
    return 0;

/*     End of CGELQF */

} /* cgelqf_ */

