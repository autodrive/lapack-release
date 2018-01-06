#line 1 "dgelqf.f"
/* dgelqf.f -- translated by f2c (version 20100827).
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

#line 1 "dgelqf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DGELQF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGELQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelqf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelqf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelqf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGELQF computes an LQ factorization of a real M-by-N matrix A: */
/* > A = L * Q. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and below the diagonal of the array */
/* >          contain the m-by-min(m,n) lower trapezoidal matrix L (L is */
/* >          lower triangular if m <= n); the elements above the diagonal, */
/* >          with the array TAU, represent the orthogonal matrix Q as a */
/* >          product of elementary reflectors (see Further Details). */
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
/* >          TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,M). */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleGEcomputational */

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
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n), */
/* >  and tau in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgelqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dgelq2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dlarfb_(char *,
	     char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), dlarft_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
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

#line 171 "dgelqf.f"
    /* Parameter adjustments */
#line 171 "dgelqf.f"
    a_dim1 = *lda;
#line 171 "dgelqf.f"
    a_offset = 1 + a_dim1;
#line 171 "dgelqf.f"
    a -= a_offset;
#line 171 "dgelqf.f"
    --tau;
#line 171 "dgelqf.f"
    --work;
#line 171 "dgelqf.f"

#line 171 "dgelqf.f"
    /* Function Body */
#line 171 "dgelqf.f"
    *info = 0;
#line 172 "dgelqf.f"
    nb = ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
#line 173 "dgelqf.f"
    lwkopt = *m * nb;
#line 174 "dgelqf.f"
    work[1] = (doublereal) lwkopt;
#line 175 "dgelqf.f"
    lquery = *lwork == -1;
#line 176 "dgelqf.f"
    if (*m < 0) {
#line 177 "dgelqf.f"
	*info = -1;
#line 178 "dgelqf.f"
    } else if (*n < 0) {
#line 179 "dgelqf.f"
	*info = -2;
#line 180 "dgelqf.f"
    } else if (*lda < max(1,*m)) {
#line 181 "dgelqf.f"
	*info = -4;
#line 182 "dgelqf.f"
    } else if (*lwork < max(1,*m) && ! lquery) {
#line 183 "dgelqf.f"
	*info = -7;
#line 184 "dgelqf.f"
    }
#line 185 "dgelqf.f"
    if (*info != 0) {
#line 186 "dgelqf.f"
	i__1 = -(*info);
#line 186 "dgelqf.f"
	xerbla_("DGELQF", &i__1, (ftnlen)6);
#line 187 "dgelqf.f"
	return 0;
#line 188 "dgelqf.f"
    } else if (lquery) {
#line 189 "dgelqf.f"
	return 0;
#line 190 "dgelqf.f"
    }

/*     Quick return if possible */

#line 194 "dgelqf.f"
    k = min(*m,*n);
#line 195 "dgelqf.f"
    if (k == 0) {
#line 196 "dgelqf.f"
	work[1] = 1.;
#line 197 "dgelqf.f"
	return 0;
#line 198 "dgelqf.f"
    }

#line 200 "dgelqf.f"
    nbmin = 2;
#line 201 "dgelqf.f"
    nx = 0;
#line 202 "dgelqf.f"
    iws = *m;
#line 203 "dgelqf.f"
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 207 "dgelqf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGELQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 207 "dgelqf.f"
	nx = max(i__1,i__2);
#line 208 "dgelqf.f"
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

#line 212 "dgelqf.f"
	    ldwork = *m;
#line 213 "dgelqf.f"
	    iws = ldwork * nb;
#line 214 "dgelqf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 219 "dgelqf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 220 "dgelqf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGELQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 220 "dgelqf.f"
		nbmin = max(i__1,i__2);
#line 222 "dgelqf.f"
	    }
#line 223 "dgelqf.f"
	}
#line 224 "dgelqf.f"
    }

#line 226 "dgelqf.f"
    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

#line 230 "dgelqf.f"
	i__1 = k - nx;
#line 230 "dgelqf.f"
	i__2 = nb;
#line 230 "dgelqf.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 231 "dgelqf.f"
	    i__3 = k - i__ + 1;
#line 231 "dgelqf.f"
	    ib = min(i__3,nb);

/*           Compute the LQ factorization of the current block */
/*           A(i:i+ib-1,i:n) */

#line 236 "dgelqf.f"
	    i__3 = *n - i__ + 1;
#line 236 "dgelqf.f"
	    dgelq2_(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
#line 238 "dgelqf.f"
	    if (i__ + ib <= *m) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

#line 243 "dgelqf.f"
		i__3 = *n - i__ + 1;
#line 243 "dgelqf.f"
		dlarft_("Forward", "Rowwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
			 (ftnlen)7);

/*              Apply H to A(i+ib:m,i:n) from the right */

#line 248 "dgelqf.f"
		i__3 = *m - i__ - ib + 1;
#line 248 "dgelqf.f"
		i__4 = *n - i__ + 1;
#line 248 "dgelqf.f"
		dlarfb_("Right", "No transpose", "Forward", "Rowwise", &i__3, 
			&i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork, (ftnlen)5, (ftnlen)12, (ftnlen)7, (
			ftnlen)7);
#line 252 "dgelqf.f"
	    }
#line 253 "dgelqf.f"
/* L10: */
#line 253 "dgelqf.f"
	}
#line 254 "dgelqf.f"
    } else {
#line 255 "dgelqf.f"
	i__ = 1;
#line 256 "dgelqf.f"
    }

/*     Use unblocked code to factor the last or only block. */

#line 260 "dgelqf.f"
    if (i__ <= k) {
#line 260 "dgelqf.f"
	i__2 = *m - i__ + 1;
#line 260 "dgelqf.f"
	i__1 = *n - i__ + 1;
#line 260 "dgelqf.f"
	dgelq2_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
#line 260 "dgelqf.f"
    }

#line 264 "dgelqf.f"
    work[1] = (doublereal) iws;
#line 265 "dgelqf.f"
    return 0;

/*     End of DGELQF */

} /* dgelqf_ */

