#line 1 "zunglq.f"
/* zunglq.f -- translated by f2c (version 20100827).
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

#line 1 "zunglq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b ZUNGLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunglq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunglq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunglq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNGLQ generates an M-by-N complex matrix Q with orthonormal rows, */
/* > which is defined as the first M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* >       Q  =  H(k)**H . . . H(2)**H H(1)**H */
/* > */
/* > as returned by ZGELQF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Q. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines the */
/* >          matrix Q. M >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by ZGELQF in the first k rows of its array argument A. */
/* >          On exit, the M-by-N matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The first dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGELQF. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,M). */
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
/* >          = 0:  successful exit; */
/* >          < 0:  if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zunglq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int zungl2_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
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
    static logical lquery;
    static integer lwkopt;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
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

#line 167 "zunglq.f"
    /* Parameter adjustments */
#line 167 "zunglq.f"
    a_dim1 = *lda;
#line 167 "zunglq.f"
    a_offset = 1 + a_dim1;
#line 167 "zunglq.f"
    a -= a_offset;
#line 167 "zunglq.f"
    --tau;
#line 167 "zunglq.f"
    --work;
#line 167 "zunglq.f"

#line 167 "zunglq.f"
    /* Function Body */
#line 167 "zunglq.f"
    *info = 0;
#line 168 "zunglq.f"
    nb = ilaenv_(&c__1, "ZUNGLQ", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
#line 169 "zunglq.f"
    lwkopt = max(1,*m) * nb;
#line 170 "zunglq.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 171 "zunglq.f"
    lquery = *lwork == -1;
#line 172 "zunglq.f"
    if (*m < 0) {
#line 173 "zunglq.f"
	*info = -1;
#line 174 "zunglq.f"
    } else if (*n < *m) {
#line 175 "zunglq.f"
	*info = -2;
#line 176 "zunglq.f"
    } else if (*k < 0 || *k > *m) {
#line 177 "zunglq.f"
	*info = -3;
#line 178 "zunglq.f"
    } else if (*lda < max(1,*m)) {
#line 179 "zunglq.f"
	*info = -5;
#line 180 "zunglq.f"
    } else if (*lwork < max(1,*m) && ! lquery) {
#line 181 "zunglq.f"
	*info = -8;
#line 182 "zunglq.f"
    }
#line 183 "zunglq.f"
    if (*info != 0) {
#line 184 "zunglq.f"
	i__1 = -(*info);
#line 184 "zunglq.f"
	xerbla_("ZUNGLQ", &i__1, (ftnlen)6);
#line 185 "zunglq.f"
	return 0;
#line 186 "zunglq.f"
    } else if (lquery) {
#line 187 "zunglq.f"
	return 0;
#line 188 "zunglq.f"
    }

/*     Quick return if possible */

#line 192 "zunglq.f"
    if (*m <= 0) {
#line 193 "zunglq.f"
	work[1].r = 1., work[1].i = 0.;
#line 194 "zunglq.f"
	return 0;
#line 195 "zunglq.f"
    }

#line 197 "zunglq.f"
    nbmin = 2;
#line 198 "zunglq.f"
    nx = 0;
#line 199 "zunglq.f"
    iws = *m;
#line 200 "zunglq.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 204 "zunglq.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZUNGLQ", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 204 "zunglq.f"
	nx = max(i__1,i__2);
#line 205 "zunglq.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 209 "zunglq.f"
	    ldwork = *m;
#line 210 "zunglq.f"
	    iws = ldwork * nb;
#line 211 "zunglq.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 216 "zunglq.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 217 "zunglq.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNGLQ", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 217 "zunglq.f"
		nbmin = max(i__1,i__2);
#line 218 "zunglq.f"
	    }
#line 219 "zunglq.f"
	}
#line 220 "zunglq.f"
    }

#line 222 "zunglq.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block. */
/*        The first kk rows are handled by the block method. */

#line 227 "zunglq.f"
	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
#line 228 "zunglq.f"
	i__1 = *k, i__2 = ki + nb;
#line 228 "zunglq.f"
	kk = min(i__1,i__2);

/*        Set A(kk+1:m,1:kk) to zero. */

#line 232 "zunglq.f"
	i__1 = kk;
#line 232 "zunglq.f"
	for (j = 1; j <= i__1; ++j) {
#line 233 "zunglq.f"
	    i__2 = *m;
#line 233 "zunglq.f"
	    for (i__ = kk + 1; i__ <= i__2; ++i__) {
#line 234 "zunglq.f"
		i__3 = i__ + j * a_dim1;
#line 234 "zunglq.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 235 "zunglq.f"
/* L10: */
#line 235 "zunglq.f"
	    }
#line 236 "zunglq.f"
/* L20: */
#line 236 "zunglq.f"
	}
#line 237 "zunglq.f"
    } else {
#line 238 "zunglq.f"
	kk = 0;
#line 239 "zunglq.f"
    }

/*     Use unblocked code for the last or only block. */

#line 243 "zunglq.f"
    if (kk < *m) {
#line 243 "zunglq.f"
	i__1 = *m - kk;
#line 243 "zunglq.f"
	i__2 = *n - kk;
#line 243 "zunglq.f"
	i__3 = *k - kk;
#line 243 "zunglq.f"
	zungl2_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
#line 243 "zunglq.f"
    }

#line 247 "zunglq.f"
    if (kk > 0) {

/*        Use blocked code */

#line 251 "zunglq.f"
	i__1 = -nb;
#line 251 "zunglq.f"
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 252 "zunglq.f"
	    i__2 = nb, i__3 = *k - i__ + 1;
#line 252 "zunglq.f"
	    ib = min(i__2,i__3);
#line 253 "zunglq.f"
	    if (i__ + ib <= *m) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 258 "zunglq.f"
		i__2 = *n - i__ + 1;
#line 258 "zunglq.f"
		zlarft_("Forward", "Rowwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)7);

/*              Apply H**H to A(i+ib:m,i:n) from the right */

#line 263 "zunglq.f"
		i__2 = *m - i__ - ib + 1;
#line 263 "zunglq.f"
		i__3 = *n - i__ + 1;
#line 263 "zunglq.f"
		zlarfb_("Right", "Conjugate transpose", "Forward", "Rowwise", 
			&i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[
			1], &ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[
			ib + 1], &ldwork, (ftnlen)5, (ftnlen)19, (ftnlen)7, (
			ftnlen)7);
#line 267 "zunglq.f"
	    }

/*           Apply H**H to columns i:n of current block */

#line 271 "zunglq.f"
	    i__2 = *n - i__ + 1;
#line 271 "zunglq.f"
	    zungl2_(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

/*           Set columns 1:i-1 of current block to zero */

#line 276 "zunglq.f"
	    i__2 = i__ - 1;
#line 276 "zunglq.f"
	    for (j = 1; j <= i__2; ++j) {
#line 277 "zunglq.f"
		i__3 = i__ + ib - 1;
#line 277 "zunglq.f"
		for (l = i__; l <= i__3; ++l) {
#line 278 "zunglq.f"
		    i__4 = l + j * a_dim1;
#line 278 "zunglq.f"
		    a[i__4].r = 0., a[i__4].i = 0.;
#line 279 "zunglq.f"
/* L30: */
#line 279 "zunglq.f"
		}
#line 280 "zunglq.f"
/* L40: */
#line 280 "zunglq.f"
	    }
#line 281 "zunglq.f"
/* L50: */
#line 281 "zunglq.f"
	}
#line 282 "zunglq.f"
    }

#line 284 "zunglq.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 285 "zunglq.f"
    return 0;

/*     End of ZUNGLQ */

} /* zunglq_ */

