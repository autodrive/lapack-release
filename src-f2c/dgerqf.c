#line 1 "dgerqf.f"
/* dgerqf.f -- translated by f2c (version 20100827).
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

#line 1 "dgerqf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DGERQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGERQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > DGERQF computes an RQ factorization of a real M-by-N matrix A: */
/* > A = R * Q. */
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
/* >          On exit, */
/* >          if m <= n, the upper triangle of the subarray */
/* >          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R; */
/* >          if m >= n, the elements on and above the (m-n)-th subdiagonal */
/* >          contain the M-by-N upper trapezoidal matrix R; */
/* >          the remaining elements, with the array TAU, represent the */
/* >          orthogonal matrix Q as a product of min(m,n) elementary */
/* >          reflectors (see Further Details). */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= M*NB, where NB is */
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
/* >  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in */
/* >  A(m-k+i,1:n-k+i-1), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgerqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dgerq2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dlarfb_(char *,
	     char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), dlarft_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
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

#line 174 "dgerqf.f"
    /* Parameter adjustments */
#line 174 "dgerqf.f"
    a_dim1 = *lda;
#line 174 "dgerqf.f"
    a_offset = 1 + a_dim1;
#line 174 "dgerqf.f"
    a -= a_offset;
#line 174 "dgerqf.f"
    --tau;
#line 174 "dgerqf.f"
    --work;
#line 174 "dgerqf.f"

#line 174 "dgerqf.f"
    /* Function Body */
#line 174 "dgerqf.f"
    *info = 0;
#line 175 "dgerqf.f"
    lquery = *lwork == -1;
#line 176 "dgerqf.f"
    if (*m < 0) {
#line 177 "dgerqf.f"
	*info = -1;
#line 178 "dgerqf.f"
    } else if (*n < 0) {
#line 179 "dgerqf.f"
	*info = -2;
#line 180 "dgerqf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "dgerqf.f"
	*info = -4;
#line 182 "dgerqf.f"
    }

#line 184 "dgerqf.f"
    if (*info == 0) {
#line 185 "dgerqf.f"
	k = min(*m,*n);
#line 186 "dgerqf.f"
	if (k == 0) {
#line 187 "dgerqf.f"
	    lwkopt = 1;
#line 188 "dgerqf.f"
	} else {
#line 189 "dgerqf.f"
	    nb = ilaenv_(&c__1, "DGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 190 "dgerqf.f"
	    lwkopt = *m * nb;
#line 191 "dgerqf.f"
	}
#line 192 "dgerqf.f"
	work[1] = (doublereal) lwkopt;

#line 194 "dgerqf.f"
	if (*lwork < max(1,*m) && ! lquery) {
#line 195 "dgerqf.f"
	    *info = -7;
#line 196 "dgerqf.f"
	}
#line 197 "dgerqf.f"
    }

#line 199 "dgerqf.f"
    if (*info != 0) {
#line 200 "dgerqf.f"
	i__1 = -(*info);
#line 200 "dgerqf.f"
	xerbla_("DGERQF", &i__1, (ftnlen)6);
#line 201 "dgerqf.f"
	return 0;
#line 202 "dgerqf.f"
    } else if (lquery) {
#line 203 "dgerqf.f"
	return 0;
#line 204 "dgerqf.f"
    }

/*     Quick return if possible */

#line 208 "dgerqf.f"
    if (k == 0) {
#line 209 "dgerqf.f"
	return 0;
#line 210 "dgerqf.f"
    }

#line 212 "dgerqf.f"
    nbmin = 2;
#line 213 "dgerqf.f"
    nx = 1;
#line 214 "dgerqf.f"
    iws = *m;
#line 215 "dgerqf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 219 "dgerqf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGERQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 219 "dgerqf.f"
	nx = max(i__1,i__2);
#line 220 "dgerqf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 224 "dgerqf.f"
	    ldwork = *m;
#line 225 "dgerqf.f"
	    iws = ldwork * nb;
#line 226 "dgerqf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 231 "dgerqf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 232 "dgerqf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGERQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 232 "dgerqf.f"
		nbmin = max(i__1,i__2);
#line 234 "dgerqf.f"
	    }
#line 235 "dgerqf.f"
	}
#line 236 "dgerqf.f"
    }

#line 238 "dgerqf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially. */
/*        The last kk rows are handled by the block method. */

#line 243 "dgerqf.f"
	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
#line 244 "dgerqf.f"
	i__1 = k, i__2 = ki + nb;
#line 244 "dgerqf.f"
	kk = min(i__1,i__2);

#line 246 "dgerqf.f"
	i__1 = k - kk + 1;
#line 246 "dgerqf.f"
	i__2 = -nb;
#line 246 "dgerqf.f"
	for (i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {
/* Computing MIN */
#line 247 "dgerqf.f"
	    i__3 = k - i__ + 1;
#line 247 "dgerqf.f"
	    ib = min(i__3,nb);

/*           Compute the RQ factorization of the current block */
/*           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1) */

#line 252 "dgerqf.f"
	    i__3 = *n - k + i__ + ib - 1;
#line 252 "dgerqf.f"
	    dgerq2_(&ib, &i__3, &a[*m - k + i__ + a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
#line 254 "dgerqf.f"
	    if (*m - k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 259 "dgerqf.f"
		i__3 = *n - k + i__ + ib - 1;
#line 259 "dgerqf.f"
		dlarft_("Backward", "Rowwise", &i__3, &ib, &a[*m - k + i__ + 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)8,
			 (ftnlen)7);

/*              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */

#line 264 "dgerqf.f"
		i__3 = *m - k + i__ - 1;
#line 264 "dgerqf.f"
		i__4 = *n - k + i__ + ib - 1;
#line 264 "dgerqf.f"
		dlarfb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
			 &i__4, &ib, &a[*m - k + i__ + a_dim1], lda, &work[1],
			 &ldwork, &a[a_offset], lda, &work[ib + 1], &ldwork, (
			ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)7);
#line 268 "dgerqf.f"
	    }
#line 269 "dgerqf.f"
/* L10: */
#line 269 "dgerqf.f"
	}
#line 270 "dgerqf.f"
	mu = *m - k + i__ + nb - 1;
#line 271 "dgerqf.f"
	nu = *n - k + i__ + nb - 1;
#line 272 "dgerqf.f"
    } else {
#line 273 "dgerqf.f"
	mu = *m;
#line 274 "dgerqf.f"
	nu = *n;
#line 275 "dgerqf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 279 "dgerqf.f"
    if (mu > 0 && nu > 0) {
#line 279 "dgerqf.f"
	dgerq2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 279 "dgerqf.f"
    }

#line 282 "dgerqf.f"
    work[1] = (doublereal) iws;
#line 283 "dgerqf.f"
    return 0;

/*     End of DGERQF */

} /* dgerqf_ */

