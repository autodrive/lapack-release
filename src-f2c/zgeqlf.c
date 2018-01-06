#line 1 "zgeqlf.f"
/* zgeqlf.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqlf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b ZGEQLF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQLF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqlf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqlf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqlf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > ZGEQLF computes a QL factorization of a complex M-by-N matrix A: */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          if m >= n, the lower triangle of the subarray */
/* >          A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L; */
/* >          if m <= n, the elements on and below the (n-m)-th */
/* >          superdiagonal contain the M-by-N lower trapezoidal matrix L; */
/* >          the remaining elements, with the array TAU, represent the */
/* >          unitary matrix Q as a product of elementary reflectors */
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

/* > \date December 2016 */

/* > \ingroup complex16GEcomputational */

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
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in */
/* >  A(1:m-k+i-1,n-k+i), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqlf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int zgeql2_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *), xerbla_(
	    char *, integer *, ftnlen);
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

#line 174 "zgeqlf.f"
    /* Parameter adjustments */
#line 174 "zgeqlf.f"
    a_dim1 = *lda;
#line 174 "zgeqlf.f"
    a_offset = 1 + a_dim1;
#line 174 "zgeqlf.f"
    a -= a_offset;
#line 174 "zgeqlf.f"
    --tau;
#line 174 "zgeqlf.f"
    --work;
#line 174 "zgeqlf.f"

#line 174 "zgeqlf.f"
    /* Function Body */
#line 174 "zgeqlf.f"
    *info = 0;
#line 175 "zgeqlf.f"
    lquery = *lwork == -1;
#line 176 "zgeqlf.f"
    if (*m < 0) {
#line 177 "zgeqlf.f"
	*info = -1;
#line 178 "zgeqlf.f"
    } else if (*n < 0) {
#line 179 "zgeqlf.f"
	*info = -2;
#line 180 "zgeqlf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "zgeqlf.f"
	*info = -4;
#line 182 "zgeqlf.f"
    }

#line 184 "zgeqlf.f"
    if (*info == 0) {
#line 185 "zgeqlf.f"
	k = min(*m,*n);
#line 186 "zgeqlf.f"
	if (k == 0) {
#line 187 "zgeqlf.f"
	    lwkopt = 1;
#line 188 "zgeqlf.f"
	} else {
#line 189 "zgeqlf.f"
	    nb = ilaenv_(&c__1, "ZGEQLF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 190 "zgeqlf.f"
	    lwkopt = *n * nb;
#line 191 "zgeqlf.f"
	}
#line 192 "zgeqlf.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 194 "zgeqlf.f"
	if (*lwork < max(1,*n) && ! lquery) {
#line 195 "zgeqlf.f"
	    *info = -7;
#line 196 "zgeqlf.f"
	}
#line 197 "zgeqlf.f"
    }

#line 199 "zgeqlf.f"
    if (*info != 0) {
#line 200 "zgeqlf.f"
	i__1 = -(*info);
#line 200 "zgeqlf.f"
	xerbla_("ZGEQLF", &i__1, (ftnlen)6);
#line 201 "zgeqlf.f"
	return 0;
#line 202 "zgeqlf.f"
    } else if (lquery) {
#line 203 "zgeqlf.f"
	return 0;
#line 204 "zgeqlf.f"
    }

/*     Quick return if possible */

#line 208 "zgeqlf.f"
    if (k == 0) {
#line 209 "zgeqlf.f"
	return 0;
#line 210 "zgeqlf.f"
    }

#line 212 "zgeqlf.f"
    nbmin = 2;
#line 213 "zgeqlf.f"
    nx = 1;
#line 214 "zgeqlf.f"
    iws = *n;
#line 215 "zgeqlf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 219 "zgeqlf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZGEQLF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 219 "zgeqlf.f"
	nx = max(i__1,i__2);
#line 220 "zgeqlf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 224 "zgeqlf.f"
	    ldwork = *n;
#line 225 "zgeqlf.f"
	    iws = ldwork * nb;
#line 226 "zgeqlf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 231 "zgeqlf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 232 "zgeqlf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZGEQLF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 232 "zgeqlf.f"
		nbmin = max(i__1,i__2);
#line 234 "zgeqlf.f"
	    }
#line 235 "zgeqlf.f"
	}
#line 236 "zgeqlf.f"
    }

#line 238 "zgeqlf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially. */
/*        The last kk columns are handled by the block method. */

#line 243 "zgeqlf.f"
	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
#line 244 "zgeqlf.f"
	i__1 = k, i__2 = ki + nb;
#line 244 "zgeqlf.f"
	kk = min(i__1,i__2);

#line 246 "zgeqlf.f"
	i__1 = k - kk + 1;
#line 246 "zgeqlf.f"
	i__2 = -nb;
#line 246 "zgeqlf.f"
	for (i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {
/* Computing MIN */
#line 247 "zgeqlf.f"
	    i__3 = k - i__ + 1;
#line 247 "zgeqlf.f"
	    ib = min(i__3,nb);

/*           Compute the QL factorization of the current block */
/*           A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1) */

#line 252 "zgeqlf.f"
	    i__3 = *m - k + i__ + ib - 1;
#line 252 "zgeqlf.f"
	    zgeql2_(&i__3, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda, &tau[
		    i__], &work[1], &iinfo);
#line 254 "zgeqlf.f"
	    if (*n - k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 259 "zgeqlf.f"
		i__3 = *m - k + i__ + ib - 1;
#line 259 "zgeqlf.f"
		zlarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - k + 
			i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork,
			 (ftnlen)8, (ftnlen)10);

/*              Apply H**H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

#line 264 "zgeqlf.f"
		i__3 = *m - k + i__ + ib - 1;
#line 264 "zgeqlf.f"
		i__4 = *n - k + i__ - 1;
#line 264 "zgeqlf.f"
		zlarfb_("Left", "Conjugate transpose", "Backward", "Columnwi"\
			"se", &i__3, &i__4, &ib, &a[(*n - k + i__) * a_dim1 + 
			1], lda, &work[1], &ldwork, &a[a_offset], lda, &work[
			ib + 1], &ldwork, (ftnlen)4, (ftnlen)19, (ftnlen)8, (
			ftnlen)10);
#line 268 "zgeqlf.f"
	    }
#line 269 "zgeqlf.f"
/* L10: */
#line 269 "zgeqlf.f"
	}
#line 270 "zgeqlf.f"
	mu = *m - k + i__ + nb - 1;
#line 271 "zgeqlf.f"
	nu = *n - k + i__ + nb - 1;
#line 272 "zgeqlf.f"
    } else {
#line 273 "zgeqlf.f"
	mu = *m;
#line 274 "zgeqlf.f"
	nu = *n;
#line 275 "zgeqlf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 279 "zgeqlf.f"
    if (mu > 0 && nu > 0) {
#line 279 "zgeqlf.f"
	zgeql2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
#line 279 "zgeqlf.f"
    }

#line 282 "zgeqlf.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 283 "zgeqlf.f"
    return 0;

/*     End of ZGEQLF */

} /* zgeqlf_ */

