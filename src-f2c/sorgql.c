#line 1 "sorgql.f"
/* sorgql.f -- translated by f2c (version 20100827).
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

#line 1 "sorgql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b SORGQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > SORGQL generates an M-by-N real matrix Q with orthonormal columns, */
/* > which is defined as the last N columns of a product of K elementary */
/* > reflectors of order M */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGEQLF. */
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
/* >          On entry, the (n-k+i)-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by SGEQLF in the last k columns of its array */
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
/* >          reflector H(i), as returned by SGEQLF. */
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

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorgql_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ib, nb, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int sorg2l_(integer *, integer *, integer *, 
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

#line 168 "sorgql.f"
    /* Parameter adjustments */
#line 168 "sorgql.f"
    a_dim1 = *lda;
#line 168 "sorgql.f"
    a_offset = 1 + a_dim1;
#line 168 "sorgql.f"
    a -= a_offset;
#line 168 "sorgql.f"
    --tau;
#line 168 "sorgql.f"
    --work;
#line 168 "sorgql.f"

#line 168 "sorgql.f"
    /* Function Body */
#line 168 "sorgql.f"
    *info = 0;
#line 169 "sorgql.f"
    lquery = *lwork == -1;
#line 170 "sorgql.f"
    if (*m < 0) {
#line 171 "sorgql.f"
	*info = -1;
#line 172 "sorgql.f"
    } else if (*n < 0 || *n > *m) {
#line 173 "sorgql.f"
	*info = -2;
#line 174 "sorgql.f"
    } else if (*k < 0 || *k > *n) {
#line 175 "sorgql.f"
	*info = -3;
#line 176 "sorgql.f"
    } else if (*lda < max(1,*m)) {
#line 177 "sorgql.f"
	*info = -5;
#line 178 "sorgql.f"
    }

#line 180 "sorgql.f"
    if (*info == 0) {
#line 181 "sorgql.f"
	if (*n == 0) {
#line 182 "sorgql.f"
	    lwkopt = 1;
#line 183 "sorgql.f"
	} else {
#line 184 "sorgql.f"
	    nb = ilaenv_(&c__1, "SORGQL", " ", m, n, k, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 185 "sorgql.f"
	    lwkopt = *n * nb;
#line 186 "sorgql.f"
	}
#line 187 "sorgql.f"
	work[1] = (doublereal) lwkopt;

#line 189 "sorgql.f"
	if (*lwork < max(1,*n) && ! lquery) {
#line 190 "sorgql.f"
	    *info = -8;
#line 191 "sorgql.f"
	}
#line 192 "sorgql.f"
    }

#line 194 "sorgql.f"
    if (*info != 0) {
#line 195 "sorgql.f"
	i__1 = -(*info);
#line 195 "sorgql.f"
	xerbla_("SORGQL", &i__1, (ftnlen)6);
#line 196 "sorgql.f"
	return 0;
#line 197 "sorgql.f"
    } else if (lquery) {
#line 198 "sorgql.f"
	return 0;
#line 199 "sorgql.f"
    }

/*     Quick return if possible */

#line 203 "sorgql.f"
    if (*n <= 0) {
#line 204 "sorgql.f"
	return 0;
#line 205 "sorgql.f"
    }

#line 207 "sorgql.f"
    nbmin = 2;
#line 208 "sorgql.f"
    nx = 0;
#line 209 "sorgql.f"
    iws = *n;
#line 210 "sorgql.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 214 "sorgql.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SORGQL", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 214 "sorgql.f"
	nx = max(i__1,i__2);
#line 215 "sorgql.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 219 "sorgql.f"
	    ldwork = *n;
#line 220 "sorgql.f"
	    iws = ldwork * nb;
#line 221 "sorgql.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 226 "sorgql.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 227 "sorgql.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SORGQL", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 227 "sorgql.f"
		nbmin = max(i__1,i__2);
#line 228 "sorgql.f"
	    }
#line 229 "sorgql.f"
	}
#line 230 "sorgql.f"
    }

#line 232 "sorgql.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block. */
/*        The last kk columns are handled by the block method. */

/* Computing MIN */
#line 237 "sorgql.f"
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
#line 237 "sorgql.f"
	kk = min(i__1,i__2);

/*        Set A(m-kk+1:m,1:n-kk) to zero. */

#line 241 "sorgql.f"
	i__1 = *n - kk;
#line 241 "sorgql.f"
	for (j = 1; j <= i__1; ++j) {
#line 242 "sorgql.f"
	    i__2 = *m;
#line 242 "sorgql.f"
	    for (i__ = *m - kk + 1; i__ <= i__2; ++i__) {
#line 243 "sorgql.f"
		a[i__ + j * a_dim1] = 0.;
#line 244 "sorgql.f"
/* L10: */
#line 244 "sorgql.f"
	    }
#line 245 "sorgql.f"
/* L20: */
#line 245 "sorgql.f"
	}
#line 246 "sorgql.f"
    } else {
#line 247 "sorgql.f"
	kk = 0;
#line 248 "sorgql.f"
    }

/*     Use unblocked code for the first or only block. */

#line 252 "sorgql.f"
    i__1 = *m - kk;
#line 252 "sorgql.f"
    i__2 = *n - kk;
#line 252 "sorgql.f"
    i__3 = *k - kk;
#line 252 "sorgql.f"
    sorg2l_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

#line 254 "sorgql.f"
    if (kk > 0) {

/*        Use blocked code */

#line 258 "sorgql.f"
	i__1 = *k;
#line 258 "sorgql.f"
	i__2 = nb;
#line 258 "sorgql.f"
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
#line 259 "sorgql.f"
	    i__3 = nb, i__4 = *k - i__ + 1;
#line 259 "sorgql.f"
	    ib = min(i__3,i__4);
#line 260 "sorgql.f"
	    if (*n - *k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 265 "sorgql.f"
		i__3 = *m - *k + i__ + ib - 1;
#line 265 "sorgql.f"
		slarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - *k + 
			i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork,
			 (ftnlen)8, (ftnlen)10);

/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

#line 270 "sorgql.f"
		i__3 = *m - *k + i__ + ib - 1;
#line 270 "sorgql.f"
		i__4 = *n - *k + i__ - 1;
#line 270 "sorgql.f"
		slarfb_("Left", "No transpose", "Backward", "Columnwise", &
			i__3, &i__4, &ib, &a[(*n - *k + i__) * a_dim1 + 1], 
			lda, &work[1], &ldwork, &a[a_offset], lda, &work[ib + 
			1], &ldwork, (ftnlen)4, (ftnlen)12, (ftnlen)8, (
			ftnlen)10);
#line 274 "sorgql.f"
	    }

/*           Apply H to rows 1:m-k+i+ib-1 of current block */

#line 278 "sorgql.f"
	    i__3 = *m - *k + i__ + ib - 1;
#line 278 "sorgql.f"
	    sorg2l_(&i__3, &ib, &ib, &a[(*n - *k + i__) * a_dim1 + 1], lda, &
		    tau[i__], &work[1], &iinfo);

/*           Set rows m-k+i+ib:m of current block to zero */

#line 283 "sorgql.f"
	    i__3 = *n - *k + i__ + ib - 1;
#line 283 "sorgql.f"
	    for (j = *n - *k + i__; j <= i__3; ++j) {
#line 284 "sorgql.f"
		i__4 = *m;
#line 284 "sorgql.f"
		for (l = *m - *k + i__ + ib; l <= i__4; ++l) {
#line 285 "sorgql.f"
		    a[l + j * a_dim1] = 0.;
#line 286 "sorgql.f"
/* L30: */
#line 286 "sorgql.f"
		}
#line 287 "sorgql.f"
/* L40: */
#line 287 "sorgql.f"
	    }
#line 288 "sorgql.f"
/* L50: */
#line 288 "sorgql.f"
	}
#line 289 "sorgql.f"
    }

#line 291 "sorgql.f"
    work[1] = (doublereal) iws;
#line 292 "sorgql.f"
    return 0;

/*     End of SORGQL */

} /* sorgql_ */

