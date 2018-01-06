#line 1 "sgerqf.f"
/* sgerqf.f -- translated by f2c (version 20100827).
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

#line 1 "sgerqf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SGERQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGERQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgerqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgerqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgerqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > SGERQF computes an RQ factorization of a real M-by-N matrix A: */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in */
/* >  A(m-k+i,1:n-k+i-1), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgerqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sgerq2_(integer *, integer *, doublereal *, 
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

#line 174 "sgerqf.f"
    /* Parameter adjustments */
#line 174 "sgerqf.f"
    a_dim1 = *lda;
#line 174 "sgerqf.f"
    a_offset = 1 + a_dim1;
#line 174 "sgerqf.f"
    a -= a_offset;
#line 174 "sgerqf.f"
    --tau;
#line 174 "sgerqf.f"
    --work;
#line 174 "sgerqf.f"

#line 174 "sgerqf.f"
    /* Function Body */
#line 174 "sgerqf.f"
    *info = 0;
#line 175 "sgerqf.f"
    lquery = *lwork == -1;
#line 176 "sgerqf.f"
    if (*m < 0) {
#line 177 "sgerqf.f"
	*info = -1;
#line 178 "sgerqf.f"
    } else if (*n < 0) {
#line 179 "sgerqf.f"
	*info = -2;
#line 180 "sgerqf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "sgerqf.f"
	*info = -4;
#line 182 "sgerqf.f"
    } else if (*lwork < max(1,*m) && ! lquery) {
#line 183 "sgerqf.f"
	*info = -7;
#line 184 "sgerqf.f"
    }

#line 186 "sgerqf.f"
    if (*info == 0) {
#line 187 "sgerqf.f"
	k = min(*m,*n);
#line 188 "sgerqf.f"
	if (k == 0) {
#line 189 "sgerqf.f"
	    lwkopt = 1;
#line 190 "sgerqf.f"
	} else {
#line 191 "sgerqf.f"
	    nb = ilaenv_(&c__1, "SGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 192 "sgerqf.f"
	    lwkopt = *m * nb;
#line 193 "sgerqf.f"
	    work[1] = (doublereal) lwkopt;
#line 194 "sgerqf.f"
	}
#line 195 "sgerqf.f"
	work[1] = (doublereal) lwkopt;

#line 197 "sgerqf.f"
	if (*lwork < max(1,*m) && ! lquery) {
#line 198 "sgerqf.f"
	    *info = -7;
#line 199 "sgerqf.f"
	}
#line 200 "sgerqf.f"
    }

#line 202 "sgerqf.f"
    if (*info != 0) {
#line 203 "sgerqf.f"
	i__1 = -(*info);
#line 203 "sgerqf.f"
	xerbla_("SGERQF", &i__1, (ftnlen)6);
#line 204 "sgerqf.f"
	return 0;
#line 205 "sgerqf.f"
    } else if (lquery) {
#line 206 "sgerqf.f"
	return 0;
#line 207 "sgerqf.f"
    }

/*     Quick return if possible */

#line 211 "sgerqf.f"
    if (k == 0) {
#line 212 "sgerqf.f"
	return 0;
#line 213 "sgerqf.f"
    }

#line 215 "sgerqf.f"
    nbmin = 2;
#line 216 "sgerqf.f"
    nx = 1;
#line 217 "sgerqf.f"
    iws = *m;
#line 218 "sgerqf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 222 "sgerqf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SGERQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 222 "sgerqf.f"
	nx = max(i__1,i__2);
#line 223 "sgerqf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 227 "sgerqf.f"
	    ldwork = *m;
#line 228 "sgerqf.f"
	    iws = ldwork * nb;
#line 229 "sgerqf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 234 "sgerqf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 235 "sgerqf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGERQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 235 "sgerqf.f"
		nbmin = max(i__1,i__2);
#line 237 "sgerqf.f"
	    }
#line 238 "sgerqf.f"
	}
#line 239 "sgerqf.f"
    }

#line 241 "sgerqf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially. */
/*        The last kk rows are handled by the block method. */

#line 246 "sgerqf.f"
	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
#line 247 "sgerqf.f"
	i__1 = k, i__2 = ki + nb;
#line 247 "sgerqf.f"
	kk = min(i__1,i__2);

#line 249 "sgerqf.f"
	i__1 = k - kk + 1;
#line 249 "sgerqf.f"
	i__2 = -nb;
#line 249 "sgerqf.f"
	for (i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {
/* Computing MIN */
#line 250 "sgerqf.f"
	    i__3 = k - i__ + 1;
#line 250 "sgerqf.f"
	    ib = min(i__3,nb);

/*           Compute the RQ factorization of the current block */
/*           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1) */

#line 255 "sgerqf.f"
	    i__3 = *n - k + i__ + ib - 1;
#line 255 "sgerqf.f"
	    sgerq2_(&ib, &i__3, &a[*m - k + i__ + a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);
#line 257 "sgerqf.f"
	    if (*m - k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 262 "sgerqf.f"
		i__3 = *n - k + i__ + ib - 1;
#line 262 "sgerqf.f"
		slarft_("Backward", "Rowwise", &i__3, &ib, &a[*m - k + i__ + 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)8,
			 (ftnlen)7);

/*              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */

#line 267 "sgerqf.f"
		i__3 = *m - k + i__ - 1;
#line 267 "sgerqf.f"
		i__4 = *n - k + i__ + ib - 1;
#line 267 "sgerqf.f"
		slarfb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
			 &i__4, &ib, &a[*m - k + i__ + a_dim1], lda, &work[1],
			 &ldwork, &a[a_offset], lda, &work[ib + 1], &ldwork, (
			ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)7);
#line 271 "sgerqf.f"
	    }
#line 272 "sgerqf.f"
/* L10: */
#line 272 "sgerqf.f"
	}
#line 273 "sgerqf.f"
	mu = *m - k + i__ + nb - 1;
#line 274 "sgerqf.f"
	nu = *n - k + i__ + nb - 1;
#line 275 "sgerqf.f"
    } else {
#line 276 "sgerqf.f"
	mu = *m;
#line 277 "sgerqf.f"
	nu = *n;
#line 278 "sgerqf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 282 "sgerqf.f"
    if (mu > 0 && nu > 0) {
#line 282 "sgerqf.f"
	sgerq2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 282 "sgerqf.f"
    }

#line 285 "sgerqf.f"
    work[1] = (doublereal) iws;
#line 286 "sgerqf.f"
    return 0;

/*     End of SGERQF */

} /* sgerqf_ */

