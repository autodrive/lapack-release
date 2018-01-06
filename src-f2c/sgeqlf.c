#line 1 "sgeqlf.f"
/* sgeqlf.f -- translated by f2c (version 20100827).
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

#line 1 "sgeqlf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SGEQLF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEQLF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqlf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqlf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqlf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > SGEQLF computes a QL factorization of a real M-by-N matrix A: */
/* > A = Q * L. */
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
/* >          if m >= n, the lower triangle of the subarray */
/* >          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L; */
/* >          if m <= n, the elements on and below the (n-m)-th */
/* >          superdiagonal contain the M-by-N lower trapezoidal matrix L; */
/* >          the remaining elements, with the array TAU, represent the */
/* >          orthogonal matrix Q as a product of elementary reflectors */
/* >          (see Further Details). */
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

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(k) . . . H(2) H(1), where k = min(m,n). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in */
/* >  A(1:m-k+i-1,n-k+i), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgeqlf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sgeql2_(integer *, integer *, doublereal *, 
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

#line 174 "sgeqlf.f"
    /* Parameter adjustments */
#line 174 "sgeqlf.f"
    a_dim1 = *lda;
#line 174 "sgeqlf.f"
    a_offset = 1 + a_dim1;
#line 174 "sgeqlf.f"
    a -= a_offset;
#line 174 "sgeqlf.f"
    --tau;
#line 174 "sgeqlf.f"
    --work;
#line 174 "sgeqlf.f"

#line 174 "sgeqlf.f"
    /* Function Body */
#line 174 "sgeqlf.f"
    *info = 0;
#line 175 "sgeqlf.f"
    lquery = *lwork == -1;
#line 176 "sgeqlf.f"
    if (*m < 0) {
#line 177 "sgeqlf.f"
	*info = -1;
#line 178 "sgeqlf.f"
    } else if (*n < 0) {
#line 179 "sgeqlf.f"
	*info = -2;
#line 180 "sgeqlf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "sgeqlf.f"
	*info = -4;
#line 182 "sgeqlf.f"
    }

#line 184 "sgeqlf.f"
    if (*info == 0) {
#line 185 "sgeqlf.f"
	k = min(*m,*n);
#line 186 "sgeqlf.f"
	if (k == 0) {
#line 187 "sgeqlf.f"
	    lwkopt = 1;
#line 188 "sgeqlf.f"
	} else {
#line 189 "sgeqlf.f"
	    nb = ilaenv_(&c__1, "SGEQLF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 190 "sgeqlf.f"
	    lwkopt = *n * nb;
#line 191 "sgeqlf.f"
	}
#line 192 "sgeqlf.f"
	work[1] = (doublereal) lwkopt;

#line 194 "sgeqlf.f"
	if (*lwork < max(1,*n) && ! lquery) {
#line 195 "sgeqlf.f"
	    *info = -7;
#line 196 "sgeqlf.f"
	}
#line 197 "sgeqlf.f"
    }

#line 199 "sgeqlf.f"
    if (*info != 0) {
#line 200 "sgeqlf.f"
	i__1 = -(*info);
#line 200 "sgeqlf.f"
	xerbla_("SGEQLF", &i__1, (ftnlen)6);
#line 201 "sgeqlf.f"
	return 0;
#line 202 "sgeqlf.f"
    } else if (lquery) {
#line 203 "sgeqlf.f"
	return 0;
#line 204 "sgeqlf.f"
    }

/*     Quick return if possible */

#line 208 "sgeqlf.f"
    if (k == 0) {
#line 209 "sgeqlf.f"
	return 0;
#line 210 "sgeqlf.f"
    }

#line 212 "sgeqlf.f"
    nbmin = 2;
#line 213 "sgeqlf.f"
    nx = 1;
#line 214 "sgeqlf.f"
    iws = *n;
#line 215 "sgeqlf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 219 "sgeqlf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SGEQLF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 219 "sgeqlf.f"
	nx = max(i__1,i__2);
#line 220 "sgeqlf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 224 "sgeqlf.f"
	    ldwork = *n;
#line 225 "sgeqlf.f"
	    iws = ldwork * nb;
#line 226 "sgeqlf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 231 "sgeqlf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 232 "sgeqlf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGEQLF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 232 "sgeqlf.f"
		nbmin = max(i__1,i__2);
#line 234 "sgeqlf.f"
	    }
#line 235 "sgeqlf.f"
	}
#line 236 "sgeqlf.f"
    }

#line 238 "sgeqlf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially. */
/*        The last kk columns are handled by the block method. */

#line 243 "sgeqlf.f"
	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
#line 244 "sgeqlf.f"
	i__1 = k, i__2 = ki + nb;
#line 244 "sgeqlf.f"
	kk = min(i__1,i__2);

#line 246 "sgeqlf.f"
	i__1 = k - kk + 1;
#line 246 "sgeqlf.f"
	i__2 = -nb;
#line 246 "sgeqlf.f"
	for (i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {
/* Computing MIN */
#line 247 "sgeqlf.f"
	    i__3 = k - i__ + 1;
#line 247 "sgeqlf.f"
	    ib = min(i__3,nb);

/*           Compute the QL factorization of the current block */
/*           A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1) */

#line 252 "sgeqlf.f"
	    i__3 = *m - k + i__ + ib - 1;
#line 252 "sgeqlf.f"
	    sgeql2_(&i__3, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda, &tau[
		    i__], &work[1], &iinfo);
#line 254 "sgeqlf.f"
	    if (*n - k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 259 "sgeqlf.f"
		i__3 = *m - k + i__ + ib - 1;
#line 259 "sgeqlf.f"
		slarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - k + 
			i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork,
			 (ftnlen)8, (ftnlen)10);

/*              Apply H**T to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

#line 264 "sgeqlf.f"
		i__3 = *m - k + i__ + ib - 1;
#line 264 "sgeqlf.f"
		i__4 = *n - k + i__ - 1;
#line 264 "sgeqlf.f"
		slarfb_("Left", "Transpose", "Backward", "Columnwise", &i__3, 
			&i__4, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda, &
			work[1], &ldwork, &a[a_offset], lda, &work[ib + 1], &
			ldwork, (ftnlen)4, (ftnlen)9, (ftnlen)8, (ftnlen)10);
#line 268 "sgeqlf.f"
	    }
#line 269 "sgeqlf.f"
/* L10: */
#line 269 "sgeqlf.f"
	}
#line 270 "sgeqlf.f"
	mu = *m - k + i__ + nb - 1;
#line 271 "sgeqlf.f"
	nu = *n - k + i__ + nb - 1;
#line 272 "sgeqlf.f"
    } else {
#line 273 "sgeqlf.f"
	mu = *m;
#line 274 "sgeqlf.f"
	nu = *n;
#line 275 "sgeqlf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 279 "sgeqlf.f"
    if (mu > 0 && nu > 0) {
#line 279 "sgeqlf.f"
	sgeql2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 279 "sgeqlf.f"
    }

#line 282 "sgeqlf.f"
    work[1] = (doublereal) iws;
#line 283 "sgeqlf.f"
    return 0;

/*     End of SGEQLF */

} /* sgeqlf_ */

