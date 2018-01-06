#line 1 "zungql.f"
/* zungql.f -- translated by f2c (version 20100827).
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

#line 1 "zungql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b ZUNGQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

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
/* > ZUNGQL generates an M-by-N complex matrix Q with orthonormal columns, */
/* > which is defined as the last N columns of a product of K elementary */
/* > reflectors of order M */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by ZGEQLF. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the (n-k+i)-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by ZGEQLF in the last k columns of its array */
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
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGEQLF. */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zungql_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, l, ib, nb, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int zung2l_(integer *, integer *, integer *, 
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

#line 168 "zungql.f"
    /* Parameter adjustments */
#line 168 "zungql.f"
    a_dim1 = *lda;
#line 168 "zungql.f"
    a_offset = 1 + a_dim1;
#line 168 "zungql.f"
    a -= a_offset;
#line 168 "zungql.f"
    --tau;
#line 168 "zungql.f"
    --work;
#line 168 "zungql.f"

#line 168 "zungql.f"
    /* Function Body */
#line 168 "zungql.f"
    *info = 0;
#line 169 "zungql.f"
    lquery = *lwork == -1;
#line 170 "zungql.f"
    if (*m < 0) {
#line 171 "zungql.f"
	*info = -1;
#line 172 "zungql.f"
    } else if (*n < 0 || *n > *m) {
#line 173 "zungql.f"
	*info = -2;
#line 174 "zungql.f"
    } else if (*k < 0 || *k > *n) {
#line 175 "zungql.f"
	*info = -3;
#line 176 "zungql.f"
    } else if (*lda < max(1,*m)) {
#line 177 "zungql.f"
	*info = -5;
#line 178 "zungql.f"
    }

#line 180 "zungql.f"
    if (*info == 0) {
#line 181 "zungql.f"
	if (*n == 0) {
#line 182 "zungql.f"
	    lwkopt = 1;
#line 183 "zungql.f"
	} else {
#line 184 "zungql.f"
	    nb = ilaenv_(&c__1, "ZUNGQL", " ", m, n, k, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 185 "zungql.f"
	    lwkopt = *n * nb;
#line 186 "zungql.f"
	}
#line 187 "zungql.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 189 "zungql.f"
	if (*lwork < max(1,*n) && ! lquery) {
#line 190 "zungql.f"
	    *info = -8;
#line 191 "zungql.f"
	}
#line 192 "zungql.f"
    }

#line 194 "zungql.f"
    if (*info != 0) {
#line 195 "zungql.f"
	i__1 = -(*info);
#line 195 "zungql.f"
	xerbla_("ZUNGQL", &i__1, (ftnlen)6);
#line 196 "zungql.f"
	return 0;
#line 197 "zungql.f"
    } else if (lquery) {
#line 198 "zungql.f"
	return 0;
#line 199 "zungql.f"
    }

/*     Quick return if possible */

#line 203 "zungql.f"
    if (*n <= 0) {
#line 204 "zungql.f"
	return 0;
#line 205 "zungql.f"
    }

#line 207 "zungql.f"
    nbmin = 2;
#line 208 "zungql.f"
    nx = 0;
#line 209 "zungql.f"
    iws = *n;
#line 210 "zungql.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 214 "zungql.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZUNGQL", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 214 "zungql.f"
	nx = max(i__1,i__2);
#line 215 "zungql.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

#line 219 "zungql.f"
	    ldwork = *n;
#line 220 "zungql.f"
	    iws = ldwork * nb;
#line 221 "zungql.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 226 "zungql.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 227 "zungql.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNGQL", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
#line 227 "zungql.f"
		nbmin = max(i__1,i__2);
#line 228 "zungql.f"
	    }
#line 229 "zungql.f"
	}
#line 230 "zungql.f"
    }

#line 232 "zungql.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block. */
/*        The last kk columns are handled by the block method. */

/* Computing MIN */
#line 237 "zungql.f"
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
#line 237 "zungql.f"
	kk = min(i__1,i__2);

/*        Set A(m-kk+1:m,1:n-kk) to zero. */

#line 241 "zungql.f"
	i__1 = *n - kk;
#line 241 "zungql.f"
	for (j = 1; j <= i__1; ++j) {
#line 242 "zungql.f"
	    i__2 = *m;
#line 242 "zungql.f"
	    for (i__ = *m - kk + 1; i__ <= i__2; ++i__) {
#line 243 "zungql.f"
		i__3 = i__ + j * a_dim1;
#line 243 "zungql.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 244 "zungql.f"
/* L10: */
#line 244 "zungql.f"
	    }
#line 245 "zungql.f"
/* L20: */
#line 245 "zungql.f"
	}
#line 246 "zungql.f"
    } else {
#line 247 "zungql.f"
	kk = 0;
#line 248 "zungql.f"
    }

/*     Use unblocked code for the first or only block. */

#line 252 "zungql.f"
    i__1 = *m - kk;
#line 252 "zungql.f"
    i__2 = *n - kk;
#line 252 "zungql.f"
    i__3 = *k - kk;
#line 252 "zungql.f"
    zung2l_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

#line 254 "zungql.f"
    if (kk > 0) {

/*        Use blocked code */

#line 258 "zungql.f"
	i__1 = *k;
#line 258 "zungql.f"
	i__2 = nb;
#line 258 "zungql.f"
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
#line 259 "zungql.f"
	    i__3 = nb, i__4 = *k - i__ + 1;
#line 259 "zungql.f"
	    ib = min(i__3,i__4);
#line 260 "zungql.f"
	    if (*n - *k + i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 265 "zungql.f"
		i__3 = *m - *k + i__ + ib - 1;
#line 265 "zungql.f"
		zlarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - *k + 
			i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork,
			 (ftnlen)8, (ftnlen)10);

/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

#line 270 "zungql.f"
		i__3 = *m - *k + i__ + ib - 1;
#line 270 "zungql.f"
		i__4 = *n - *k + i__ - 1;
#line 270 "zungql.f"
		zlarfb_("Left", "No transpose", "Backward", "Columnwise", &
			i__3, &i__4, &ib, &a[(*n - *k + i__) * a_dim1 + 1], 
			lda, &work[1], &ldwork, &a[a_offset], lda, &work[ib + 
			1], &ldwork, (ftnlen)4, (ftnlen)12, (ftnlen)8, (
			ftnlen)10);
#line 274 "zungql.f"
	    }

/*           Apply H to rows 1:m-k+i+ib-1 of current block */

#line 278 "zungql.f"
	    i__3 = *m - *k + i__ + ib - 1;
#line 278 "zungql.f"
	    zung2l_(&i__3, &ib, &ib, &a[(*n - *k + i__) * a_dim1 + 1], lda, &
		    tau[i__], &work[1], &iinfo);

/*           Set rows m-k+i+ib:m of current block to zero */

#line 283 "zungql.f"
	    i__3 = *n - *k + i__ + ib - 1;
#line 283 "zungql.f"
	    for (j = *n - *k + i__; j <= i__3; ++j) {
#line 284 "zungql.f"
		i__4 = *m;
#line 284 "zungql.f"
		for (l = *m - *k + i__ + ib; l <= i__4; ++l) {
#line 285 "zungql.f"
		    i__5 = l + j * a_dim1;
#line 285 "zungql.f"
		    a[i__5].r = 0., a[i__5].i = 0.;
#line 286 "zungql.f"
/* L30: */
#line 286 "zungql.f"
		}
#line 287 "zungql.f"
/* L40: */
#line 287 "zungql.f"
	    }
#line 288 "zungql.f"
/* L50: */
#line 288 "zungql.f"
	}
#line 289 "zungql.f"
    }

#line 291 "zungql.f"
    work[1].r = (doublereal) iws, work[1].i = 0.;
#line 292 "zungql.f"
    return 0;

/*     End of ZUNGQL */

} /* zungql_ */

