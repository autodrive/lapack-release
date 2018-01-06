#line 1 "sorglq.f"
/* sorglq.f -- translated by f2c (version 20100827).
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

#line 1 "sorglq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SORGLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorglq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorglq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorglq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORGLQ generates an M-by-N real matrix Q with orthonormal rows, */
/* > which is defined as the first M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGELQF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by SGELQF in the first k rows of its array argument A. */
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
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SGELQF. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorglq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sorgl2_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    slarfb_(char *, char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
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

#line 167 "sorglq.f"
    /* Parameter adjustments */
#line 167 "sorglq.f"
    a_dim1 = *lda;
#line 167 "sorglq.f"
    a_offset = 1 + a_dim1;
#line 167 "sorglq.f"
    a -= a_offset;
#line 167 "sorglq.f"
    --tau;
#line 167 "sorglq.f"
    --work;
#line 167 "sorglq.f"

#line 167 "sorglq.f"
    /* Function Body */
#line 167 "sorglq.f"
    *info = 0;
#line 168 "sorglq.f"
    nb = ilaenv_(&c__1, "SORGLQ", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
#line 169 "sorglq.f"
    lwkopt = max(1,*m) * nb;
#line 170 "sorglq.f"
    work[1] = (doublereal) lwkopt;
#line 171 "sorglq.f"
    lquery = *lwork == -1;
#line 172 "sorglq.f"
    if (*m < 0) {
#line 173 "sorglq.f"
	*info = -1;
#line 174 "sorglq.f"
    } else if (*n < *m) {
#line 175 "sorglq.f"
	*info = -2;
#line 176 "sorglq.f"
    } else if (*k < 0 || *k > *m) {
#line 177 "sorglq.f"
	*info = -3;
#line 178 "sorglq.f"
    } else if (*lda < max(1,*m)) {
#line 179 "sorglq.f"
	*info = -5;
#line 180 "sorglq.f"
    } else if (*lwork < max(1,*m) && ! lquery) {
#line 181 "sorglq.f"
	*info = -8;
#line 182 "sorglq.f"
    }
#line 183 "sorglq.f"
    if (*info != 0) {
#line 184 "sorglq.f"
	i__1 = -(*info);
#line 184 "sorglq.f"
	xerbla_("SORGLQ", &i__1, (ftnlen)6);
#line 185 "sorglq.f"
	return 0;
#line 186 "sorglq.f"
    } else if (lquery) {
#line 187 "sorglq.f"
	return 0;
#line 188 "sorglq.f"
    }

/*     Quick return if possible */

#line 192 "sorglq.f"
    if (*m <= 0) {
#line 193 "sorglq.f"
	work[1] = 1.;
#line 194 "sorglq.f"
	return 0;
#line 195 "sorglq.f"
    }

#line 197 "sorglq.f"
    nbmin = 2;
#line 198 "sorglq.f"
    nx = 0;
#line 199 "sorglq.f"
    iws = *m;
#line 200 "sorglq.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 204 "sorglq.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SORGLQ", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 204 "sorglq.f"
	nx = max(i__1,i__2);
#line 205 "sorglq.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 209 "sorglq.f"
	    ldwork = *m;
#line 210 "sorglq.f"
	    iws = ldwork * nb;
#line 211 "sorglq.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 216 "sorglq.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 217 "sorglq.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SORGLQ", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 217 "sorglq.f"
		nbmin = max(i__1,i__2);
#line 218 "sorglq.f"
	    }
#line 219 "sorglq.f"
	}
#line 220 "sorglq.f"
    }

#line 222 "sorglq.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block. */
/*        The first kk rows are handled by the block method. */

#line 227 "sorglq.f"
	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
#line 228 "sorglq.f"
	i__1 = *k, i__2 = ki + nb;
#line 228 "sorglq.f"
	kk = min(i__1,i__2);

/*        Set A(kk+1:m,1:kk) to zero. */

#line 232 "sorglq.f"
	i__1 = kk;
#line 232 "sorglq.f"
	for (j = 1; j <= i__1; ++j) {
#line 233 "sorglq.f"
	    i__2 = *m;
#line 233 "sorglq.f"
	    for (i__ = kk + 1; i__ <= i__2; ++i__) {
#line 234 "sorglq.f"
		a[i__ + j * a_dim1] = 0.;
#line 235 "sorglq.f"
/* L10: */
#line 235 "sorglq.f"
	    }
#line 236 "sorglq.f"
/* L20: */
#line 236 "sorglq.f"
	}
#line 237 "sorglq.f"
    } else {
#line 238 "sorglq.f"
	kk = 0;
#line 239 "sorglq.f"
    }

/*     Use unblocked code for the last or only block. */

#line 243 "sorglq.f"
    if (kk < *m) {
#line 243 "sorglq.f"
	i__1 = *m - kk;
#line 243 "sorglq.f"
	i__2 = *n - kk;
#line 243 "sorglq.f"
	i__3 = *k - kk;
#line 243 "sorglq.f"
	sorgl2_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
#line 243 "sorglq.f"
    }

#line 247 "sorglq.f"
    if (kk > 0) {

/*        Use blocked code */

#line 251 "sorglq.f"
	i__1 = -nb;
#line 251 "sorglq.f"
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 252 "sorglq.f"
	    i__2 = nb, i__3 = *k - i__ + 1;
#line 252 "sorglq.f"
	    ib = min(i__2,i__3);
#line 253 "sorglq.f"
	    if (i__ + ib <= *m) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 258 "sorglq.f"
		i__2 = *n - i__ + 1;
#line 258 "sorglq.f"
		slarft_("Forward", "Rowwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)7);

/*              Apply H**T to A(i+ib:m,i:n) from the right */

#line 263 "sorglq.f"
		i__2 = *m - i__ - ib + 1;
#line 263 "sorglq.f"
		i__3 = *n - i__ + 1;
#line 263 "sorglq.f"
		slarfb_("Right", "Transpose", "Forward", "Rowwise", &i__2, &
			i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork, (ftnlen)5, (ftnlen)9, (ftnlen)7, (ftnlen)
			7);
#line 267 "sorglq.f"
	    }

/*           Apply H**T to columns i:n of current block */

#line 271 "sorglq.f"
	    i__2 = *n - i__ + 1;
#line 271 "sorglq.f"
	    sorgl2_(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

/*           Set columns 1:i-1 of current block to zero */

#line 276 "sorglq.f"
	    i__2 = i__ - 1;
#line 276 "sorglq.f"
	    for (j = 1; j <= i__2; ++j) {
#line 277 "sorglq.f"
		i__3 = i__ + ib - 1;
#line 277 "sorglq.f"
		for (l = i__; l <= i__3; ++l) {
#line 278 "sorglq.f"
		    a[l + j * a_dim1] = 0.;
#line 279 "sorglq.f"
/* L30: */
#line 279 "sorglq.f"
		}
#line 280 "sorglq.f"
/* L40: */
#line 280 "sorglq.f"
	    }
#line 281 "sorglq.f"
/* L50: */
#line 281 "sorglq.f"
	}
#line 282 "sorglq.f"
    }

#line 284 "sorglq.f"
    work[1] = (doublereal) iws;
#line 285 "sorglq.f"
    return 0;

/*     End of SORGLQ */

} /* sorglq_ */

