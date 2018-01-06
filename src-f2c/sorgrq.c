#line 1 "sorgrq.f"
/* sorgrq.f -- translated by f2c (version 20100827).
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

#line 1 "sorgrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SORGRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > SORGRQ generates an M-by-N real matrix Q with orthonormal rows, */
/* > which is defined as the last M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* >       Q  =  H(1) H(2) . . . H(k) */
/* > */
/* > as returned by SGERQF. */
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
/* >          On entry, the (m-k+i)-th row must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by SGERQF in the last k rows of its array argument */
/* >          A. */
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
/* >          reflector H(i), as returned by SGERQF. */
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
/* Subroutine */ int sorgrq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ib, nb, ii, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sorgr2_(integer *, integer *, integer *, 
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

#line 168 "sorgrq.f"
    /* Parameter adjustments */
#line 168 "sorgrq.f"
    a_dim1 = *lda;
#line 168 "sorgrq.f"
    a_offset = 1 + a_dim1;
#line 168 "sorgrq.f"
    a -= a_offset;
#line 168 "sorgrq.f"
    --tau;
#line 168 "sorgrq.f"
    --work;
#line 168 "sorgrq.f"

#line 168 "sorgrq.f"
    /* Function Body */
#line 168 "sorgrq.f"
    *info = 0;
#line 169 "sorgrq.f"
    lquery = *lwork == -1;
#line 170 "sorgrq.f"
    if (*m < 0) {
#line 171 "sorgrq.f"
	*info = -1;
#line 172 "sorgrq.f"
    } else if (*n < *m) {
#line 173 "sorgrq.f"
	*info = -2;
#line 174 "sorgrq.f"
    } else if (*k < 0 || *k > *m) {
#line 175 "sorgrq.f"
	*info = -3;
#line 176 "sorgrq.f"
    } else if (*lda < max(1,*m)) {
#line 177 "sorgrq.f"
	*info = -5;
#line 178 "sorgrq.f"
    }

#line 180 "sorgrq.f"
    if (*info == 0) {
#line 181 "sorgrq.f"
	if (*m <= 0) {
#line 182 "sorgrq.f"
	    lwkopt = 1;
#line 183 "sorgrq.f"
	} else {
#line 184 "sorgrq.f"
	    nb = ilaenv_(&c__1, "SORGRQ", " ", m, n, k, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 185 "sorgrq.f"
	    lwkopt = *m * nb;
#line 186 "sorgrq.f"
	}
#line 187 "sorgrq.f"
	work[1] = (doublereal) lwkopt;

#line 189 "sorgrq.f"
	if (*lwork < max(1,*m) && ! lquery) {
#line 190 "sorgrq.f"
	    *info = -8;
#line 191 "sorgrq.f"
	}
#line 192 "sorgrq.f"
    }

#line 194 "sorgrq.f"
    if (*info != 0) {
#line 195 "sorgrq.f"
	i__1 = -(*info);
#line 195 "sorgrq.f"
	xerbla_("SORGRQ", &i__1, (ftnlen)6);
#line 196 "sorgrq.f"
	return 0;
#line 197 "sorgrq.f"
    } else if (lquery) {
#line 198 "sorgrq.f"
	return 0;
#line 199 "sorgrq.f"
    }

/*     Quick return if possible */

#line 203 "sorgrq.f"
    if (*m <= 0) {
#line 204 "sorgrq.f"
	return 0;
#line 205 "sorgrq.f"
    }

#line 207 "sorgrq.f"
    nbmin = 2;
#line 208 "sorgrq.f"
    nx = 0;
#line 209 "sorgrq.f"
    iws = *m;
#line 210 "sorgrq.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 214 "sorgrq.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SORGRQ", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 214 "sorgrq.f"
	nx = max(i__1,i__2);
#line 215 "sorgrq.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 219 "sorgrq.f"
	    ldwork = *m;
#line 220 "sorgrq.f"
	    iws = ldwork * nb;
#line 221 "sorgrq.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 226 "sorgrq.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 227 "sorgrq.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SORGRQ", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 227 "sorgrq.f"
		nbmin = max(i__1,i__2);
#line 228 "sorgrq.f"
	    }
#line 229 "sorgrq.f"
	}
#line 230 "sorgrq.f"
    }

#line 232 "sorgrq.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block. */
/*        The last kk rows are handled by the block method. */

/* Computing MIN */
#line 237 "sorgrq.f"
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
#line 237 "sorgrq.f"
	kk = min(i__1,i__2);

/*        Set A(1:m-kk,n-kk+1:n) to zero. */

#line 241 "sorgrq.f"
	i__1 = *n;
#line 241 "sorgrq.f"
	for (j = *n - kk + 1; j <= i__1; ++j) {
#line 242 "sorgrq.f"
	    i__2 = *m - kk;
#line 242 "sorgrq.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 243 "sorgrq.f"
		a[i__ + j * a_dim1] = 0.;
#line 244 "sorgrq.f"
/* L10: */
#line 244 "sorgrq.f"
	    }
#line 245 "sorgrq.f"
/* L20: */
#line 245 "sorgrq.f"
	}
#line 246 "sorgrq.f"
    } else {
#line 247 "sorgrq.f"
	kk = 0;
#line 248 "sorgrq.f"
    }

/*     Use unblocked code for the first or only block. */

#line 252 "sorgrq.f"
    i__1 = *m - kk;
#line 252 "sorgrq.f"
    i__2 = *n - kk;
#line 252 "sorgrq.f"
    i__3 = *k - kk;
#line 252 "sorgrq.f"
    sorgr2_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

#line 254 "sorgrq.f"
    if (kk > 0) {

/*        Use blocked code */

#line 258 "sorgrq.f"
	i__1 = *k;
#line 258 "sorgrq.f"
	i__2 = nb;
#line 258 "sorgrq.f"
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
#line 259 "sorgrq.f"
	    i__3 = nb, i__4 = *k - i__ + 1;
#line 259 "sorgrq.f"
	    ib = min(i__3,i__4);
#line 260 "sorgrq.f"
	    ii = *m - *k + i__;
#line 261 "sorgrq.f"
	    if (ii > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 266 "sorgrq.f"
		i__3 = *n - *k + i__ + ib - 1;
#line 266 "sorgrq.f"
		slarft_("Backward", "Rowwise", &i__3, &ib, &a[ii + a_dim1], 
			lda, &tau[i__], &work[1], &ldwork, (ftnlen)8, (ftnlen)
			7);

/*              Apply H**T to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */

#line 271 "sorgrq.f"
		i__3 = ii - 1;
#line 271 "sorgrq.f"
		i__4 = *n - *k + i__ + ib - 1;
#line 271 "sorgrq.f"
		slarfb_("Right", "Transpose", "Backward", "Rowwise", &i__3, &
			i__4, &ib, &a[ii + a_dim1], lda, &work[1], &ldwork, &
			a[a_offset], lda, &work[ib + 1], &ldwork, (ftnlen)5, (
			ftnlen)9, (ftnlen)8, (ftnlen)7);
#line 274 "sorgrq.f"
	    }

/*           Apply H**T to columns 1:n-k+i+ib-1 of current block */

#line 278 "sorgrq.f"
	    i__3 = *n - *k + i__ + ib - 1;
#line 278 "sorgrq.f"
	    sorgr2_(&ib, &i__3, &ib, &a[ii + a_dim1], lda, &tau[i__], &work[1]
		    , &iinfo);

/*           Set columns n-k+i+ib:n of current block to zero */

#line 283 "sorgrq.f"
	    i__3 = *n;
#line 283 "sorgrq.f"
	    for (l = *n - *k + i__ + ib; l <= i__3; ++l) {
#line 284 "sorgrq.f"
		i__4 = ii + ib - 1;
#line 284 "sorgrq.f"
		for (j = ii; j <= i__4; ++j) {
#line 285 "sorgrq.f"
		    a[j + l * a_dim1] = 0.;
#line 286 "sorgrq.f"
/* L30: */
#line 286 "sorgrq.f"
		}
#line 287 "sorgrq.f"
/* L40: */
#line 287 "sorgrq.f"
	    }
#line 288 "sorgrq.f"
/* L50: */
#line 288 "sorgrq.f"
	}
#line 289 "sorgrq.f"
    }

#line 291 "sorgrq.f"
    work[1] = (doublereal) iws;
#line 292 "sorgrq.f"
    return 0;

/*     End of SORGRQ */

} /* sorgrq_ */

