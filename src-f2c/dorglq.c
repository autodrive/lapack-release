#line 1 "dorglq.f"
/* dorglq.f -- translated by f2c (version 20100827).
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

#line 1 "dorglq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DORGLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorglq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorglq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorglq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGLQ generates an M-by-N real matrix Q with orthonormal rows, */
/* > which is defined as the first M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by DGELQF. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by DGELQF in the first k rows of its array argument A. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGELQF. */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorglq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dorgl2_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dlarfb_(char *, char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), dlarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
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

#line 167 "dorglq.f"
    /* Parameter adjustments */
#line 167 "dorglq.f"
    a_dim1 = *lda;
#line 167 "dorglq.f"
    a_offset = 1 + a_dim1;
#line 167 "dorglq.f"
    a -= a_offset;
#line 167 "dorglq.f"
    --tau;
#line 167 "dorglq.f"
    --work;
#line 167 "dorglq.f"

#line 167 "dorglq.f"
    /* Function Body */
#line 167 "dorglq.f"
    *info = 0;
#line 168 "dorglq.f"
    nb = ilaenv_(&c__1, "DORGLQ", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
#line 169 "dorglq.f"
    lwkopt = max(1,*m) * nb;
#line 170 "dorglq.f"
    work[1] = (doublereal) lwkopt;
#line 171 "dorglq.f"
    lquery = *lwork == -1;
#line 172 "dorglq.f"
    if (*m < 0) {
#line 173 "dorglq.f"
	*info = -1;
#line 174 "dorglq.f"
    } else if (*n < *m) {
#line 175 "dorglq.f"
	*info = -2;
#line 176 "dorglq.f"
    } else if (*k < 0 || *k > *m) {
#line 177 "dorglq.f"
	*info = -3;
#line 178 "dorglq.f"
    } else if (*lda < max(1,*m)) {
#line 179 "dorglq.f"
	*info = -5;
#line 180 "dorglq.f"
    } else if (*lwork < max(1,*m) && ! lquery) {
#line 181 "dorglq.f"
	*info = -8;
#line 182 "dorglq.f"
    }
#line 183 "dorglq.f"
    if (*info != 0) {
#line 184 "dorglq.f"
	i__1 = -(*info);
#line 184 "dorglq.f"
	xerbla_("DORGLQ", &i__1, (ftnlen)6);
#line 185 "dorglq.f"
	return 0;
#line 186 "dorglq.f"
    } else if (lquery) {
#line 187 "dorglq.f"
	return 0;
#line 188 "dorglq.f"
    }

/*     Quick return if possible */

#line 192 "dorglq.f"
    if (*m <= 0) {
#line 193 "dorglq.f"
	work[1] = 1.;
#line 194 "dorglq.f"
	return 0;
#line 195 "dorglq.f"
    }

#line 197 "dorglq.f"
    nbmin = 2;
#line 198 "dorglq.f"
    nx = 0;
#line 199 "dorglq.f"
    iws = *m;
#line 200 "dorglq.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 204 "dorglq.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGLQ", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 204 "dorglq.f"
	nx = max(i__1,i__2);
#line 205 "dorglq.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 209 "dorglq.f"
	    ldwork = *m;
#line 210 "dorglq.f"
	    iws = ldwork * nb;
#line 211 "dorglq.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 216 "dorglq.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 217 "dorglq.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGLQ", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 217 "dorglq.f"
		nbmin = max(i__1,i__2);
#line 218 "dorglq.f"
	    }
#line 219 "dorglq.f"
	}
#line 220 "dorglq.f"
    }

#line 222 "dorglq.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block. */
/*        The first kk rows are handled by the block method. */

#line 227 "dorglq.f"
	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
#line 228 "dorglq.f"
	i__1 = *k, i__2 = ki + nb;
#line 228 "dorglq.f"
	kk = min(i__1,i__2);

/*        Set A(kk+1:m,1:kk) to zero. */

#line 232 "dorglq.f"
	i__1 = kk;
#line 232 "dorglq.f"
	for (j = 1; j <= i__1; ++j) {
#line 233 "dorglq.f"
	    i__2 = *m;
#line 233 "dorglq.f"
	    for (i__ = kk + 1; i__ <= i__2; ++i__) {
#line 234 "dorglq.f"
		a[i__ + j * a_dim1] = 0.;
#line 235 "dorglq.f"
/* L10: */
#line 235 "dorglq.f"
	    }
#line 236 "dorglq.f"
/* L20: */
#line 236 "dorglq.f"
	}
#line 237 "dorglq.f"
    } else {
#line 238 "dorglq.f"
	kk = 0;
#line 239 "dorglq.f"
    }

/*     Use unblocked code for the last or only block. */

#line 243 "dorglq.f"
    if (kk < *m) {
#line 243 "dorglq.f"
	i__1 = *m - kk;
#line 243 "dorglq.f"
	i__2 = *n - kk;
#line 243 "dorglq.f"
	i__3 = *k - kk;
#line 243 "dorglq.f"
	dorgl2_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
#line 243 "dorglq.f"
    }

#line 247 "dorglq.f"
    if (kk > 0) {

/*        Use blocked code */

#line 251 "dorglq.f"
	i__1 = -nb;
#line 251 "dorglq.f"
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 252 "dorglq.f"
	    i__2 = nb, i__3 = *k - i__ + 1;
#line 252 "dorglq.f"
	    ib = min(i__2,i__3);
#line 253 "dorglq.f"
	    if (i__ + ib <= *m) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 258 "dorglq.f"
		i__2 = *n - i__ + 1;
#line 258 "dorglq.f"
		dlarft_("Forward", "Rowwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)7);

/*              Apply H**T to A(i+ib:m,i:n) from the right */

#line 263 "dorglq.f"
		i__2 = *m - i__ - ib + 1;
#line 263 "dorglq.f"
		i__3 = *n - i__ + 1;
#line 263 "dorglq.f"
		dlarfb_("Right", "Transpose", "Forward", "Rowwise", &i__2, &
			i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork, (ftnlen)5, (ftnlen)9, (ftnlen)7, (ftnlen)
			7);
#line 267 "dorglq.f"
	    }

/*           Apply H**T to columns i:n of current block */

#line 271 "dorglq.f"
	    i__2 = *n - i__ + 1;
#line 271 "dorglq.f"
	    dorgl2_(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

/*           Set columns 1:i-1 of current block to zero */

#line 276 "dorglq.f"
	    i__2 = i__ - 1;
#line 276 "dorglq.f"
	    for (j = 1; j <= i__2; ++j) {
#line 277 "dorglq.f"
		i__3 = i__ + ib - 1;
#line 277 "dorglq.f"
		for (l = i__; l <= i__3; ++l) {
#line 278 "dorglq.f"
		    a[l + j * a_dim1] = 0.;
#line 279 "dorglq.f"
/* L30: */
#line 279 "dorglq.f"
		}
#line 280 "dorglq.f"
/* L40: */
#line 280 "dorglq.f"
	    }
#line 281 "dorglq.f"
/* L50: */
#line 281 "dorglq.f"
	}
#line 282 "dorglq.f"
    }

#line 284 "dorglq.f"
    work[1] = (doublereal) iws;
#line 285 "dorglq.f"
    return 0;

/*     End of DORGLQ */

} /* dorglq_ */

