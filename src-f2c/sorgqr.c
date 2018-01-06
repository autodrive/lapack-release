#line 1 "sorgqr.f"
/* sorgqr.f -- translated by f2c (version 20100827).
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

#line 1 "sorgqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SORGQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > SORGQR generates an M-by-N real matrix Q with orthonormal columns, */
/* > which is defined as the first N columns of a product of K elementary */
/* > reflectors of order M */
/* > */
/* >       Q  =  H(1) H(2) . . . H(k) */
/* > */
/* > as returned by SGEQRF. */
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
/* >          The number of columns of the matrix Q. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines the */
/* >          matrix Q. N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the i-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by SGEQRF in the first k columns of its array */
/* >          argument A. */
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
/* >          reflector H(i), as returned by SGEQRF. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N). */
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
/* >          < 0:  if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorgqr_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sorg2r_(integer *, integer *, integer *, 
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

#line 168 "sorgqr.f"
    /* Parameter adjustments */
#line 168 "sorgqr.f"
    a_dim1 = *lda;
#line 168 "sorgqr.f"
    a_offset = 1 + a_dim1;
#line 168 "sorgqr.f"
    a -= a_offset;
#line 168 "sorgqr.f"
    --tau;
#line 168 "sorgqr.f"
    --work;
#line 168 "sorgqr.f"

#line 168 "sorgqr.f"
    /* Function Body */
#line 168 "sorgqr.f"
    *info = 0;
#line 169 "sorgqr.f"
    nb = ilaenv_(&c__1, "SORGQR", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
#line 170 "sorgqr.f"
    lwkopt = max(1,*n) * nb;
#line 171 "sorgqr.f"
    work[1] = (doublereal) lwkopt;
#line 172 "sorgqr.f"
    lquery = *lwork == -1;
#line 173 "sorgqr.f"
    if (*m < 0) {
#line 174 "sorgqr.f"
	*info = -1;
#line 175 "sorgqr.f"
    } else if (*n < 0 || *n > *m) {
#line 176 "sorgqr.f"
	*info = -2;
#line 177 "sorgqr.f"
    } else if (*k < 0 || *k > *n) {
#line 178 "sorgqr.f"
	*info = -3;
#line 179 "sorgqr.f"
    } else if (*lda < max(1,*m)) {
#line 180 "sorgqr.f"
	*info = -5;
#line 181 "sorgqr.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 182 "sorgqr.f"
	*info = -8;
#line 183 "sorgqr.f"
    }
#line 184 "sorgqr.f"
    if (*info != 0) {
#line 185 "sorgqr.f"
	i__1 = -(*info);
#line 185 "sorgqr.f"
	xerbla_("SORGQR", &i__1, (ftnlen)6);
#line 186 "sorgqr.f"
	return 0;
#line 187 "sorgqr.f"
    } else if (lquery) {
#line 188 "sorgqr.f"
	return 0;
#line 189 "sorgqr.f"
    }

/*     Quick return if possible */

#line 193 "sorgqr.f"
    if (*n <= 0) {
#line 194 "sorgqr.f"
	work[1] = 1.;
#line 195 "sorgqr.f"
	return 0;
#line 196 "sorgqr.f"
    }

#line 198 "sorgqr.f"
    nbmin = 2;
#line 199 "sorgqr.f"
    nx = 0;
#line 200 "sorgqr.f"
    iws = *n;
#line 201 "sorgqr.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 205 "sorgqr.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SORGQR", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 205 "sorgqr.f"
	nx = max(i__1,i__2);
#line 206 "sorgqr.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 210 "sorgqr.f"
	    ldwork = *n;
#line 211 "sorgqr.f"
	    iws = ldwork * nb;
#line 212 "sorgqr.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 217 "sorgqr.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 218 "sorgqr.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SORGQR", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 218 "sorgqr.f"
		nbmin = max(i__1,i__2);
#line 219 "sorgqr.f"
	    }
#line 220 "sorgqr.f"
	}
#line 221 "sorgqr.f"
    }

#line 223 "sorgqr.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block. */
/*        The first kk columns are handled by the block method. */

#line 228 "sorgqr.f"
	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
#line 229 "sorgqr.f"
	i__1 = *k, i__2 = ki + nb;
#line 229 "sorgqr.f"
	kk = min(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

#line 233 "sorgqr.f"
	i__1 = *n;
#line 233 "sorgqr.f"
	for (j = kk + 1; j <= i__1; ++j) {
#line 234 "sorgqr.f"
	    i__2 = kk;
#line 234 "sorgqr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 235 "sorgqr.f"
		a[i__ + j * a_dim1] = 0.;
#line 236 "sorgqr.f"
/* L10: */
#line 236 "sorgqr.f"
	    }
#line 237 "sorgqr.f"
/* L20: */
#line 237 "sorgqr.f"
	}
#line 238 "sorgqr.f"
    } else {
#line 239 "sorgqr.f"
	kk = 0;
#line 240 "sorgqr.f"
    }

/*     Use unblocked code for the last or only block. */

#line 244 "sorgqr.f"
    if (kk < *n) {
#line 244 "sorgqr.f"
	i__1 = *m - kk;
#line 244 "sorgqr.f"
	i__2 = *n - kk;
#line 244 "sorgqr.f"
	i__3 = *k - kk;
#line 244 "sorgqr.f"
	sorg2r_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
#line 244 "sorgqr.f"
    }

#line 248 "sorgqr.f"
    if (kk > 0) {

/*        Use blocked code */

#line 252 "sorgqr.f"
	i__1 = -nb;
#line 252 "sorgqr.f"
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 253 "sorgqr.f"
	    i__2 = nb, i__3 = *k - i__ + 1;
#line 253 "sorgqr.f"
	    ib = min(i__2,i__3);
#line 254 "sorgqr.f"
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 259 "sorgqr.f"
		i__2 = *m - i__ + 1;
#line 259 "sorgqr.f"
		slarft_("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)10);

/*              Apply H to A(i:m,i+ib:n) from the left */

#line 264 "sorgqr.f"
		i__2 = *m - i__ + 1;
#line 264 "sorgqr.f"
		i__3 = *n - i__ - ib + 1;
#line 264 "sorgqr.f"
		slarfb_("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[
			1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &
			work[ib + 1], &ldwork, (ftnlen)4, (ftnlen)12, (ftnlen)
			7, (ftnlen)10);
#line 268 "sorgqr.f"
	    }

/*           Apply H to rows i:m of current block */

#line 272 "sorgqr.f"
	    i__2 = *m - i__ + 1;
#line 272 "sorgqr.f"
	    sorg2r_(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

/*           Set rows 1:i-1 of current block to zero */

#line 277 "sorgqr.f"
	    i__2 = i__ + ib - 1;
#line 277 "sorgqr.f"
	    for (j = i__; j <= i__2; ++j) {
#line 278 "sorgqr.f"
		i__3 = i__ - 1;
#line 278 "sorgqr.f"
		for (l = 1; l <= i__3; ++l) {
#line 279 "sorgqr.f"
		    a[l + j * a_dim1] = 0.;
#line 280 "sorgqr.f"
/* L30: */
#line 280 "sorgqr.f"
		}
#line 281 "sorgqr.f"
/* L40: */
#line 281 "sorgqr.f"
	    }
#line 282 "sorgqr.f"
/* L50: */
#line 282 "sorgqr.f"
	}
#line 283 "sorgqr.f"
    }

#line 285 "sorgqr.f"
    work[1] = (doublereal) iws;
#line 286 "sorgqr.f"
    return 0;

/*     End of SORGQR */

} /* sorgqr_ */

