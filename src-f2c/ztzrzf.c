#line 1 "ztzrzf.f"
/* ztzrzf.f -- translated by f2c (version 20100827).
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

#line 1 "ztzrzf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b ZTZRZF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTZRZF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztzrzf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztzrzf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztzrzf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTZRZF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A */
/* > to upper triangular form by means of unitary transformations. */
/* > */
/* > The upper trapezoidal matrix A is factored as */
/* > */
/* >    A = ( R  0 ) * Z, */
/* > */
/* > where Z is an N-by-N unitary matrix and R is an M-by-M upper */
/* > triangular matrix. */
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
/* >          The number of columns of the matrix A.  N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the leading M-by-N upper trapezoidal part of the */
/* >          array A must contain the matrix to be factorized. */
/* >          On exit, the leading M-by-M upper triangular part of A */
/* >          contains the upper triangular matrix R, and elements M+1 to */
/* >          N of the first M rows of A, with the array TAU, represent the */
/* >          unitary matrix Z as a product of M elementary reflectors. */
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
/* >          TAU is COMPLEX*16 array, dimension (M) */
/* >          The scalar factors of the elementary reflectors. */
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
/* >          The dimension of the array WORK.  LWORK >= max(1,M). */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The N-by-N matrix Z can be computed by */
/* > */
/* >     Z =  Z(1)*Z(2)* ... *Z(M) */
/* > */
/* >  where each N-by-N Z(k) is given by */
/* > */
/* >     Z(k) = I - tau(k)*v(k)*v(k)**H */
/* > */
/* >  with v(k) is the kth row vector of the M-by-N matrix */
/* > */
/* >     V = ( I   A(:,M+1:N) ) */
/* > */
/* >  I is the M-by-M identity matrix, A(:,M+1:N) */
/* >  is the output stored in A on exit from DTZRZF, */
/* >  and tau(k) is the kth element of the array TAU. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztzrzf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, m1, ib, nb, ki, kk, mu, nx, iws, nbmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkmin, ldwork;
    extern /* Subroutine */ int zlarzb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zlarzt_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zlatrz_(integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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

#line 191 "ztzrzf.f"
    /* Parameter adjustments */
#line 191 "ztzrzf.f"
    a_dim1 = *lda;
#line 191 "ztzrzf.f"
    a_offset = 1 + a_dim1;
#line 191 "ztzrzf.f"
    a -= a_offset;
#line 191 "ztzrzf.f"
    --tau;
#line 191 "ztzrzf.f"
    --work;
#line 191 "ztzrzf.f"

#line 191 "ztzrzf.f"
    /* Function Body */
#line 191 "ztzrzf.f"
    *info = 0;
#line 192 "ztzrzf.f"
    lquery = *lwork == -1;
#line 193 "ztzrzf.f"
    if (*m < 0) {
#line 194 "ztzrzf.f"
	*info = -1;
#line 195 "ztzrzf.f"
    } else if (*n < *m) {
#line 196 "ztzrzf.f"
	*info = -2;
#line 197 "ztzrzf.f"
    } else if (*lda < max(1,*m)) {
#line 198 "ztzrzf.f"
	*info = -4;
#line 199 "ztzrzf.f"
    }

#line 201 "ztzrzf.f"
    if (*info == 0) {
#line 202 "ztzrzf.f"
	if (*m == 0 || *m == *n) {
#line 203 "ztzrzf.f"
	    lwkopt = 1;
#line 204 "ztzrzf.f"
	    lwkmin = 1;
#line 205 "ztzrzf.f"
	} else {

/*           Determine the block size. */

#line 209 "ztzrzf.f"
	    nb = ilaenv_(&c__1, "ZGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 210 "ztzrzf.f"
	    lwkopt = *m * nb;
#line 211 "ztzrzf.f"
	    lwkmin = max(1,*m);
#line 212 "ztzrzf.f"
	}
#line 213 "ztzrzf.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 215 "ztzrzf.f"
	if (*lwork < lwkmin && ! lquery) {
#line 216 "ztzrzf.f"
	    *info = -7;
#line 217 "ztzrzf.f"
	}
#line 218 "ztzrzf.f"
    }

#line 220 "ztzrzf.f"
    if (*info != 0) {
#line 221 "ztzrzf.f"
	i__1 = -(*info);
#line 221 "ztzrzf.f"
	xerbla_("ZTZRZF", &i__1, (ftnlen)6);
#line 222 "ztzrzf.f"
	return 0;
#line 223 "ztzrzf.f"
    } else if (lquery) {
#line 224 "ztzrzf.f"
	return 0;
#line 225 "ztzrzf.f"
    }

/*     Quick return if possible */

#line 229 "ztzrzf.f"
    if (*m == 0) {
#line 230 "ztzrzf.f"
	return 0;
#line 231 "ztzrzf.f"
    } else if (*m == *n) {
#line 232 "ztzrzf.f"
	i__1 = *n;
#line 232 "ztzrzf.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "ztzrzf.f"
	    i__2 = i__;
#line 233 "ztzrzf.f"
	    tau[i__2].r = 0., tau[i__2].i = 0.;
#line 234 "ztzrzf.f"
/* L10: */
#line 234 "ztzrzf.f"
	}
#line 235 "ztzrzf.f"
	return 0;
#line 236 "ztzrzf.f"
    }

#line 238 "ztzrzf.f"
    nbmin = 2;
#line 239 "ztzrzf.f"
    nx = 1;
#line 240 "ztzrzf.f"
    iws = *m;
#line 241 "ztzrzf.f"
    if (nb > 1 && nb < *m) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 245 "ztzrzf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "ZGERQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 245 "ztzrzf.f"
	nx = max(i__1,i__2);
#line 246 "ztzrzf.f"
	if (nx < *m) {

/*           Determine if workspace is large enough for blocked code. */

#line 250 "ztzrzf.f"
	    ldwork = *m;
#line 251 "ztzrzf.f"
	    iws = ldwork * nb;
#line 252 "ztzrzf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 257 "ztzrzf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 258 "ztzrzf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZGERQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 258 "ztzrzf.f"
		nbmin = max(i__1,i__2);
#line 260 "ztzrzf.f"
	    }
#line 261 "ztzrzf.f"
	}
#line 262 "ztzrzf.f"
    }

#line 264 "ztzrzf.f"
    if (nb >= nbmin && nb < *m && nx < *m) {

/*        Use blocked code initially. */
/*        The last kk rows are handled by the block method. */

/* Computing MIN */
#line 269 "ztzrzf.f"
	i__1 = *m + 1;
#line 269 "ztzrzf.f"
	m1 = min(i__1,*n);
#line 270 "ztzrzf.f"
	ki = (*m - nx - 1) / nb * nb;
/* Computing MIN */
#line 271 "ztzrzf.f"
	i__1 = *m, i__2 = ki + nb;
#line 271 "ztzrzf.f"
	kk = min(i__1,i__2);

#line 273 "ztzrzf.f"
	i__1 = *m - kk + 1;
#line 273 "ztzrzf.f"
	i__2 = -nb;
#line 273 "ztzrzf.f"
	for (i__ = *m - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; 
		i__ += i__2) {
/* Computing MIN */
#line 274 "ztzrzf.f"
	    i__3 = *m - i__ + 1;
#line 274 "ztzrzf.f"
	    ib = min(i__3,nb);

/*           Compute the TZ factorization of the current block */
/*           A(i:i+ib-1,i:n) */

#line 279 "ztzrzf.f"
	    i__3 = *n - i__ + 1;
#line 279 "ztzrzf.f"
	    i__4 = *n - *m;
#line 279 "ztzrzf.f"
	    zlatrz_(&ib, &i__3, &i__4, &a[i__ + i__ * a_dim1], lda, &tau[i__],
		     &work[1]);
#line 281 "ztzrzf.f"
	    if (i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 286 "ztzrzf.f"
		i__3 = *n - *m;
#line 286 "ztzrzf.f"
		zlarzt_("Backward", "Rowwise", &i__3, &ib, &a[i__ + m1 * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)8,
			 (ftnlen)7);

/*              Apply H to A(1:i-1,i:n) from the right */

#line 291 "ztzrzf.f"
		i__3 = i__ - 1;
#line 291 "ztzrzf.f"
		i__4 = *n - i__ + 1;
#line 291 "ztzrzf.f"
		i__5 = *n - *m;
#line 291 "ztzrzf.f"
		zlarzb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
			 &i__4, &ib, &i__5, &a[i__ + m1 * a_dim1], lda, &work[
			1], &ldwork, &a[i__ * a_dim1 + 1], lda, &work[ib + 1],
			 &ldwork, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)7)
			;
#line 295 "ztzrzf.f"
	    }
#line 296 "ztzrzf.f"
/* L20: */
#line 296 "ztzrzf.f"
	}
#line 297 "ztzrzf.f"
	mu = i__ + nb - 1;
#line 298 "ztzrzf.f"
    } else {
#line 299 "ztzrzf.f"
	mu = *m;
#line 300 "ztzrzf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 304 "ztzrzf.f"
    if (mu > 0) {
#line 304 "ztzrzf.f"
	i__2 = *n - *m;
#line 304 "ztzrzf.f"
	zlatrz_(&mu, n, &i__2, &a[a_offset], lda, &tau[1], &work[1]);
#line 304 "ztzrzf.f"
    }

#line 307 "ztzrzf.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 309 "ztzrzf.f"
    return 0;

/*     End of ZTZRZF */

} /* ztzrzf_ */

