#line 1 "cungrq.f"
/* cungrq.f -- translated by f2c (version 20100827).
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

#line 1 "cungrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b CUNGRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNGRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNGRQ generates an M-by-N complex matrix Q with orthonormal rows, */
/* > which is defined as the last M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* >       Q  =  H(1)**H H(2)**H . . . H(k)**H */
/* > */
/* > as returned by CGERQF. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the (m-k+i)-th row must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by CGERQF in the last k rows of its array argument */
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
/* >          TAU is COMPLEX array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CGERQF. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cungrq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, l, ib, nb, ii, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int cungr2_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), clarfb_(char *, char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , ftnlen, ftnlen, ftnlen, ftnlen), clarft_(char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

#line 168 "cungrq.f"
    /* Parameter adjustments */
#line 168 "cungrq.f"
    a_dim1 = *lda;
#line 168 "cungrq.f"
    a_offset = 1 + a_dim1;
#line 168 "cungrq.f"
    a -= a_offset;
#line 168 "cungrq.f"
    --tau;
#line 168 "cungrq.f"
    --work;
#line 168 "cungrq.f"

#line 168 "cungrq.f"
    /* Function Body */
#line 168 "cungrq.f"
    *info = 0;
#line 169 "cungrq.f"
    lquery = *lwork == -1;
#line 170 "cungrq.f"
    if (*m < 0) {
#line 171 "cungrq.f"
	*info = -1;
#line 172 "cungrq.f"
    } else if (*n < *m) {
#line 173 "cungrq.f"
	*info = -2;
#line 174 "cungrq.f"
    } else if (*k < 0 || *k > *m) {
#line 175 "cungrq.f"
	*info = -3;
#line 176 "cungrq.f"
    } else if (*lda < max(1,*m)) {
#line 177 "cungrq.f"
	*info = -5;
#line 178 "cungrq.f"
    }

#line 180 "cungrq.f"
    if (*info == 0) {
#line 181 "cungrq.f"
	if (*m <= 0) {
#line 182 "cungrq.f"
	    lwkopt = 1;
#line 183 "cungrq.f"
	} else {
#line 184 "cungrq.f"
	    nb = ilaenv_(&c__1, "CUNGRQ", " ", m, n, k, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 185 "cungrq.f"
	    lwkopt = *m * nb;
#line 186 "cungrq.f"
	}
#line 187 "cungrq.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 189 "cungrq.f"
	if (*lwork < max(1,*m) && ! lquery) {
#line 190 "cungrq.f"
	    *info = -8;
#line 191 "cungrq.f"
	}
#line 192 "cungrq.f"
    }

#line 194 "cungrq.f"
    if (*info != 0) {
#line 195 "cungrq.f"
	i__1 = -(*info);
#line 195 "cungrq.f"
	xerbla_("CUNGRQ", &i__1, (ftnlen)6);
#line 196 "cungrq.f"
	return 0;
#line 197 "cungrq.f"
    } else if (lquery) {
#line 198 "cungrq.f"
	return 0;
#line 199 "cungrq.f"
    }

/*     Quick return if possible */

#line 203 "cungrq.f"
    if (*m <= 0) {
#line 204 "cungrq.f"
	return 0;
#line 205 "cungrq.f"
    }

#line 207 "cungrq.f"
    nbmin = 2;
#line 208 "cungrq.f"
    nx = 0;
#line 209 "cungrq.f"
    iws = *m;
#line 210 "cungrq.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 214 "cungrq.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "CUNGRQ", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 214 "cungrq.f"
	nx = max(i__1,i__2);
#line 215 "cungrq.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 219 "cungrq.f"
	    ldwork = *m;
#line 220 "cungrq.f"
	    iws = ldwork * nb;
#line 221 "cungrq.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 226 "cungrq.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 227 "cungrq.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "CUNGRQ", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 227 "cungrq.f"
		nbmin = max(i__1,i__2);
#line 228 "cungrq.f"
	    }
#line 229 "cungrq.f"
	}
#line 230 "cungrq.f"
    }

#line 232 "cungrq.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block. */
/*        The last kk rows are handled by the block method. */

/* Computing MIN */
#line 237 "cungrq.f"
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
#line 237 "cungrq.f"
	kk = min(i__1,i__2);

/*        Set A(1:m-kk,n-kk+1:n) to zero. */

#line 241 "cungrq.f"
	i__1 = *n;
#line 241 "cungrq.f"
	for (j = *n - kk + 1; j <= i__1; ++j) {
#line 242 "cungrq.f"
	    i__2 = *m - kk;
#line 242 "cungrq.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 243 "cungrq.f"
		i__3 = i__ + j * a_dim1;
#line 243 "cungrq.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 244 "cungrq.f"
/* L10: */
#line 244 "cungrq.f"
	    }
#line 245 "cungrq.f"
/* L20: */
#line 245 "cungrq.f"
	}
#line 246 "cungrq.f"
    } else {
#line 247 "cungrq.f"
	kk = 0;
#line 248 "cungrq.f"
    }

/*     Use unblocked code for the first or only block. */

#line 252 "cungrq.f"
    i__1 = *m - kk;
#line 252 "cungrq.f"
    i__2 = *n - kk;
#line 252 "cungrq.f"
    i__3 = *k - kk;
#line 252 "cungrq.f"
    cungr2_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

#line 254 "cungrq.f"
    if (kk > 0) {

/*        Use blocked code */

#line 258 "cungrq.f"
	i__1 = *k;
#line 258 "cungrq.f"
	i__2 = nb;
#line 258 "cungrq.f"
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
#line 259 "cungrq.f"
	    i__3 = nb, i__4 = *k - i__ + 1;
#line 259 "cungrq.f"
	    ib = min(i__3,i__4);
#line 260 "cungrq.f"
	    ii = *m - *k + i__;
#line 261 "cungrq.f"
	    if (ii > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 266 "cungrq.f"
		i__3 = *n - *k + i__ + ib - 1;
#line 266 "cungrq.f"
		clarft_("Backward", "Rowwise", &i__3, &ib, &a[ii + a_dim1], 
			lda, &tau[i__], &work[1], &ldwork, (ftnlen)8, (ftnlen)
			7);

/*              Apply H**H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */

#line 271 "cungrq.f"
		i__3 = ii - 1;
#line 271 "cungrq.f"
		i__4 = *n - *k + i__ + ib - 1;
#line 271 "cungrq.f"
		clarfb_("Right", "Conjugate transpose", "Backward", "Rowwise",
			 &i__3, &i__4, &ib, &a[ii + a_dim1], lda, &work[1], &
			ldwork, &a[a_offset], lda, &work[ib + 1], &ldwork, (
			ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)7);
#line 275 "cungrq.f"
	    }

/*           Apply H**H to columns 1:n-k+i+ib-1 of current block */

#line 279 "cungrq.f"
	    i__3 = *n - *k + i__ + ib - 1;
#line 279 "cungrq.f"
	    cungr2_(&ib, &i__3, &ib, &a[ii + a_dim1], lda, &tau[i__], &work[1]
		    , &iinfo);

/*           Set columns n-k+i+ib:n of current block to zero */

#line 284 "cungrq.f"
	    i__3 = *n;
#line 284 "cungrq.f"
	    for (l = *n - *k + i__ + ib; l <= i__3; ++l) {
#line 285 "cungrq.f"
		i__4 = ii + ib - 1;
#line 285 "cungrq.f"
		for (j = ii; j <= i__4; ++j) {
#line 286 "cungrq.f"
		    i__5 = j + l * a_dim1;
#line 286 "cungrq.f"
		    a[i__5].r = 0., a[i__5].i = 0.;
#line 287 "cungrq.f"
/* L30: */
#line 287 "cungrq.f"
		}
#line 288 "cungrq.f"
/* L40: */
#line 288 "cungrq.f"
	    }
#line 289 "cungrq.f"
/* L50: */
#line 289 "cungrq.f"
	}
#line 290 "cungrq.f"
    }

#line 292 "cungrq.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 293 "cungrq.f"
    return 0;

/*     End of CUNGRQ */

} /* cungrq_ */

