#line 1 "cgerqf.f"
/* cgerqf.f -- translated by f2c (version 20100827).
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

#line 1 "cgerqf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b CGERQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGERQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgerqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgerqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgerqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > CGERQF computes an RQ factorization of a complex M-by-N matrix A: */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          if m <= n, the upper triangle of the subarray */
/* >          A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R; */
/* >          if m >= n, the elements on and above the (m-n)-th subdiagonal */
/* >          contain the M-by-N upper trapezoidal matrix R; */
/* >          the remaining elements, with the array TAU, represent the */
/* >          unitary matrix Q as a product of min(m,n) elementary */
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

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1)**H H(2)**H . . . H(k)**H, where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on */
/* >  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgerqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int cgerq2_(integer *, integer *, doublecomplex *,
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

#line 174 "cgerqf.f"
    /* Parameter adjustments */
#line 174 "cgerqf.f"
    a_dim1 = *lda;
#line 174 "cgerqf.f"
    a_offset = 1 + a_dim1;
#line 174 "cgerqf.f"
    a -= a_offset;
#line 174 "cgerqf.f"
    --tau;
#line 174 "cgerqf.f"
    --work;
#line 174 "cgerqf.f"

#line 174 "cgerqf.f"
    /* Function Body */
#line 174 "cgerqf.f"
    *info = 0;
#line 175 "cgerqf.f"
    lquery = *lwork == -1;
#line 176 "cgerqf.f"
    if (*m < 0) {
#line 177 "cgerqf.f"
	*info = -1;
#line 178 "cgerqf.f"
    } else if (*n < 0) {
#line 179 "cgerqf.f"
	*info = -2;
#line 180 "cgerqf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "cgerqf.f"
	*info = -4;
#line 182 "cgerqf.f"
    }

#line 184 "cgerqf.f"
    if (*info == 0) {
#line 185 "cgerqf.f"
	k = min(*m,*n);
#line 186 "cgerqf.f"
	if (k == 0) {
#line 187 "cgerqf.f"
	    lwkopt = 1;
#line 188 "cgerqf.f"
	} else {
#line 189 "cgerqf.f"
	    nb = ilaenv_(&c__1, "CGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 190 "cgerqf.f"
	    lwkopt = *m * nb;
#line 191 "cgerqf.f"
	}
#line 192 "cgerqf.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 194 "cgerqf.f"
	if (*lwork < max(1,*m) && ! lquery) {
#line 195 "cgerqf.f"
	    *info = -7;
#line 196 "cgerqf.f"
	}
#line 197 "cgerqf.f"
    }

#line 199 "cgerqf.f"
    if (*info != 0) {
#line 200 "cgerqf.f"
	i__1 = -(*info);
#line 200 "cgerqf.f"
	xerbla_("CGERQF", &i__1, (ftnlen)6);
#line 201 "cgerqf.f"
	return 0;
#line 202 "cgerqf.f"
    } else if (lquery) {
#line 203 "cgerqf.f"
	return 0;
#line 204 "cgerqf.f"
    }

/*     Quick return if possible */

#line 208 "cgerqf.f"
    if (k == 0) {
#line 209 "cgerqf.f"
	return 0;
#line 210 "cgerqf.f"
    }

#line 212 "cgerqf.f"
    nbmin = 2;
#line 213 "cgerqf.f"
    nx = 1;
#line 214 "cgerqf.f"
    iws = *m;
#line 215 "cgerqf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 219 "cgerqf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "CGERQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 219 "cgerqf.f"
	nx = max(i__1,i__2);
#line 220 "cgerqf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 224 "cgerqf.f"
	    ldwork = *m;
#line 225 "cgerqf.f"
	    iws = ldwork * nb;
#line 226 "cgerqf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 231 "cgerqf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 232 "cgerqf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "CGERQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 232 "cgerqf.f"
		nbmin = max(i__1,i__2);
#line 234 "cgerqf.f"
	    }
#line 235 "cgerqf.f"
	}
#line 236 "cgerqf.f"
    }

#line 238 "cgerqf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially. */
/*        The last kk rows are handled by the block method. */

#line 243 "cgerqf.f"
	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
#line 244 "cgerqf.f"
	i__1 = k, i__2 = ki + nb;
#line 244 "cgerqf.f"
	kk = min(i__1,i__2);

#line 246 "cgerqf.f"
	i__1 = k - kk + 1;
#line 246 "cgerqf.f"
	i__2 = -nb;
#line 246 "cgerqf.f"
	for (i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {
/* Computing MIN */
#line 247 "cgerqf.f"
	    i__3 = k - i__ + 1;
#line 247 "cgerqf.f"
	    ib = min(i__3,nb);

/*           Compute the RQ factorization of the current block */
/*           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1) */

#line 252 "cgerqf.f"
	    i__3 = *n - k + i__ + ib - 1;
#line 252 "cgerqf.f"
	    cgerq2_(&ib, &i__3, &a[*m - k + i__ + a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
#line 254 "cgerqf.f"
	    if (*m - k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 259 "cgerqf.f"
		i__3 = *n - k + i__ + ib - 1;
#line 259 "cgerqf.f"
		clarft_("Backward", "Rowwise", &i__3, &ib, &a[*m - k + i__ + 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)8,
			 (ftnlen)7);

/*              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */

#line 264 "cgerqf.f"
		i__3 = *m - k + i__ - 1;
#line 264 "cgerqf.f"
		i__4 = *n - k + i__ + ib - 1;
#line 264 "cgerqf.f"
		clarfb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
			 &i__4, &ib, &a[*m - k + i__ + a_dim1], lda, &work[1],
			 &ldwork, &a[a_offset], lda, &work[ib + 1], &ldwork, (
			ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)7);
#line 268 "cgerqf.f"
	    }
#line 269 "cgerqf.f"
/* L10: */
#line 269 "cgerqf.f"
	}
#line 270 "cgerqf.f"
	mu = *m - k + i__ + nb - 1;
#line 271 "cgerqf.f"
	nu = *n - k + i__ + nb - 1;
#line 272 "cgerqf.f"
    } else {
#line 273 "cgerqf.f"
	mu = *m;
#line 274 "cgerqf.f"
	nu = *n;
#line 275 "cgerqf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 279 "cgerqf.f"
    if (mu > 0 && nu > 0) {
#line 279 "cgerqf.f"
	cgerq2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 279 "cgerqf.f"
    }

#line 282 "cgerqf.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 283 "cgerqf.f"
    return 0;

/*     End of CGERQF */

} /* cgerqf_ */

