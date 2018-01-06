#line 1 "dorgql.f"
/* dorgql.f -- translated by f2c (version 20100827).
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

#line 1 "dorgql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DORGQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > DORGQL generates an M-by-N real matrix Q with orthonormal columns, */
/* > which is defined as the last N columns of a product of K elementary */
/* > reflectors of order M */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by DGEQLF. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the (n-k+i)-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by DGEQLF in the last k columns of its array */
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
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGEQLF. */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorgql_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ib, nb, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dorg2l_(integer *, integer *, integer *, 
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

#line 168 "dorgql.f"
    /* Parameter adjustments */
#line 168 "dorgql.f"
    a_dim1 = *lda;
#line 168 "dorgql.f"
    a_offset = 1 + a_dim1;
#line 168 "dorgql.f"
    a -= a_offset;
#line 168 "dorgql.f"
    --tau;
#line 168 "dorgql.f"
    --work;
#line 168 "dorgql.f"

#line 168 "dorgql.f"
    /* Function Body */
#line 168 "dorgql.f"
    *info = 0;
#line 169 "dorgql.f"
    lquery = *lwork == -1;
#line 170 "dorgql.f"
    if (*m < 0) {
#line 171 "dorgql.f"
	*info = -1;
#line 172 "dorgql.f"
    } else if (*n < 0 || *n > *m) {
#line 173 "dorgql.f"
	*info = -2;
#line 174 "dorgql.f"
    } else if (*k < 0 || *k > *n) {
#line 175 "dorgql.f"
	*info = -3;
#line 176 "dorgql.f"
    } else if (*lda < max(1,*m)) {
#line 177 "dorgql.f"
	*info = -5;
#line 178 "dorgql.f"
    }

#line 180 "dorgql.f"
    if (*info == 0) {
#line 181 "dorgql.f"
	if (*n == 0) {
#line 182 "dorgql.f"
	    lwkopt = 1;
#line 183 "dorgql.f"
	} else {
#line 184 "dorgql.f"
	    nb = ilaenv_(&c__1, "DORGQL", " ", m, n, k, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 185 "dorgql.f"
	    lwkopt = *n * nb;
#line 186 "dorgql.f"
	}
#line 187 "dorgql.f"
	work[1] = (doublereal) lwkopt;

#line 189 "dorgql.f"
	if (*lwork < max(1,*n) && ! lquery) {
#line 190 "dorgql.f"
	    *info = -8;
#line 191 "dorgql.f"
	}
#line 192 "dorgql.f"
    }

#line 194 "dorgql.f"
    if (*info != 0) {
#line 195 "dorgql.f"
	i__1 = -(*info);
#line 195 "dorgql.f"
	xerbla_("DORGQL", &i__1, (ftnlen)6);
#line 196 "dorgql.f"
	return 0;
#line 197 "dorgql.f"
    } else if (lquery) {
#line 198 "dorgql.f"
	return 0;
#line 199 "dorgql.f"
    }

/*     Quick return if possible */

#line 203 "dorgql.f"
    if (*n <= 0) {
#line 204 "dorgql.f"
	return 0;
#line 205 "dorgql.f"
    }

#line 207 "dorgql.f"
    nbmin = 2;
#line 208 "dorgql.f"
    nx = 0;
#line 209 "dorgql.f"
    iws = *n;
#line 210 "dorgql.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 214 "dorgql.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQL", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 214 "dorgql.f"
	nx = max(i__1,i__2);
#line 215 "dorgql.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 219 "dorgql.f"
	    ldwork = *n;
#line 220 "dorgql.f"
	    iws = ldwork * nb;
#line 221 "dorgql.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 226 "dorgql.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 227 "dorgql.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQL", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 227 "dorgql.f"
		nbmin = max(i__1,i__2);
#line 228 "dorgql.f"
	    }
#line 229 "dorgql.f"
	}
#line 230 "dorgql.f"
    }

#line 232 "dorgql.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block. */
/*        The last kk columns are handled by the block method. */

/* Computing MIN */
#line 237 "dorgql.f"
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
#line 237 "dorgql.f"
	kk = min(i__1,i__2);

/*        Set A(m-kk+1:m,1:n-kk) to zero. */

#line 241 "dorgql.f"
	i__1 = *n - kk;
#line 241 "dorgql.f"
	for (j = 1; j <= i__1; ++j) {
#line 242 "dorgql.f"
	    i__2 = *m;
#line 242 "dorgql.f"
	    for (i__ = *m - kk + 1; i__ <= i__2; ++i__) {
#line 243 "dorgql.f"
		a[i__ + j * a_dim1] = 0.;
#line 244 "dorgql.f"
/* L10: */
#line 244 "dorgql.f"
	    }
#line 245 "dorgql.f"
/* L20: */
#line 245 "dorgql.f"
	}
#line 246 "dorgql.f"
    } else {
#line 247 "dorgql.f"
	kk = 0;
#line 248 "dorgql.f"
    }

/*     Use unblocked code for the first or only block. */

#line 252 "dorgql.f"
    i__1 = *m - kk;
#line 252 "dorgql.f"
    i__2 = *n - kk;
#line 252 "dorgql.f"
    i__3 = *k - kk;
#line 252 "dorgql.f"
    dorg2l_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

#line 254 "dorgql.f"
    if (kk > 0) {

/*        Use blocked code */

#line 258 "dorgql.f"
	i__1 = *k;
#line 258 "dorgql.f"
	i__2 = nb;
#line 258 "dorgql.f"
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
#line 259 "dorgql.f"
	    i__3 = nb, i__4 = *k - i__ + 1;
#line 259 "dorgql.f"
	    ib = min(i__3,i__4);
#line 260 "dorgql.f"
	    if (*n - *k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 265 "dorgql.f"
		i__3 = *m - *k + i__ + ib - 1;
#line 265 "dorgql.f"
		dlarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - *k + 
			i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork,
			 (ftnlen)8, (ftnlen)10);

/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

#line 270 "dorgql.f"
		i__3 = *m - *k + i__ + ib - 1;
#line 270 "dorgql.f"
		i__4 = *n - *k + i__ - 1;
#line 270 "dorgql.f"
		dlarfb_("Left", "No transpose", "Backward", "Columnwise", &
			i__3, &i__4, &ib, &a[(*n - *k + i__) * a_dim1 + 1], 
			lda, &work[1], &ldwork, &a[a_offset], lda, &work[ib + 
			1], &ldwork, (ftnlen)4, (ftnlen)12, (ftnlen)8, (
			ftnlen)10);
#line 274 "dorgql.f"
	    }

/*           Apply H to rows 1:m-k+i+ib-1 of current block */

#line 278 "dorgql.f"
	    i__3 = *m - *k + i__ + ib - 1;
#line 278 "dorgql.f"
	    dorg2l_(&i__3, &ib, &ib, &a[(*n - *k + i__) * a_dim1 + 1], lda, &
		    tau[i__], &work[1], &iinfo);

/*           Set rows m-k+i+ib:m of current block to zero */

#line 283 "dorgql.f"
	    i__3 = *n - *k + i__ + ib - 1;
#line 283 "dorgql.f"
	    for (j = *n - *k + i__; j <= i__3; ++j) {
#line 284 "dorgql.f"
		i__4 = *m;
#line 284 "dorgql.f"
		for (l = *m - *k + i__ + ib; l <= i__4; ++l) {
#line 285 "dorgql.f"
		    a[l + j * a_dim1] = 0.;
#line 286 "dorgql.f"
/* L30: */
#line 286 "dorgql.f"
		}
#line 287 "dorgql.f"
/* L40: */
#line 287 "dorgql.f"
	    }
#line 288 "dorgql.f"
/* L50: */
#line 288 "dorgql.f"
	}
#line 289 "dorgql.f"
    }

#line 291 "dorgql.f"
    work[1] = (doublereal) iws;
#line 292 "dorgql.f"
    return 0;

/*     End of DORGQL */

} /* dorgql_ */

