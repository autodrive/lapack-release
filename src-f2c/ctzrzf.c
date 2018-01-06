#line 1 "ctzrzf.f"
/* ctzrzf.f -- translated by f2c (version 20100827).
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

#line 1 "ctzrzf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b CTZRZF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTZRZF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctzrzf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctzrzf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctzrzf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTZRZF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          TAU is COMPLEX array, dimension (M) */
/* >          The scalar factors of the elementary reflectors. */
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

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int ctzrzf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, m1, ib, nb, ki, kk, mu, nx, iws, nbmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), clarzb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int clarzt_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), clatrz_(integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *);
    static integer lwkmin, ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.1) -- */
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

#line 191 "ctzrzf.f"
    /* Parameter adjustments */
#line 191 "ctzrzf.f"
    a_dim1 = *lda;
#line 191 "ctzrzf.f"
    a_offset = 1 + a_dim1;
#line 191 "ctzrzf.f"
    a -= a_offset;
#line 191 "ctzrzf.f"
    --tau;
#line 191 "ctzrzf.f"
    --work;
#line 191 "ctzrzf.f"

#line 191 "ctzrzf.f"
    /* Function Body */
#line 191 "ctzrzf.f"
    *info = 0;
#line 192 "ctzrzf.f"
    lquery = *lwork == -1;
#line 193 "ctzrzf.f"
    if (*m < 0) {
#line 194 "ctzrzf.f"
	*info = -1;
#line 195 "ctzrzf.f"
    } else if (*n < *m) {
#line 196 "ctzrzf.f"
	*info = -2;
#line 197 "ctzrzf.f"
    } else if (*lda < max(1,*m)) {
#line 198 "ctzrzf.f"
	*info = -4;
#line 199 "ctzrzf.f"
    }

#line 201 "ctzrzf.f"
    if (*info == 0) {
#line 202 "ctzrzf.f"
	if (*m == 0 || *m == *n) {
#line 203 "ctzrzf.f"
	    lwkopt = 1;
#line 204 "ctzrzf.f"
	    lwkmin = 1;
#line 205 "ctzrzf.f"
	} else {

/*           Determine the block size. */

#line 209 "ctzrzf.f"
	    nb = ilaenv_(&c__1, "CGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 210 "ctzrzf.f"
	    lwkopt = *m * nb;
#line 211 "ctzrzf.f"
	    lwkmin = max(1,*m);
#line 212 "ctzrzf.f"
	}
#line 213 "ctzrzf.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 215 "ctzrzf.f"
	if (*lwork < lwkmin && ! lquery) {
#line 216 "ctzrzf.f"
	    *info = -7;
#line 217 "ctzrzf.f"
	}
#line 218 "ctzrzf.f"
    }

#line 220 "ctzrzf.f"
    if (*info != 0) {
#line 221 "ctzrzf.f"
	i__1 = -(*info);
#line 221 "ctzrzf.f"
	xerbla_("CTZRZF", &i__1, (ftnlen)6);
#line 222 "ctzrzf.f"
	return 0;
#line 223 "ctzrzf.f"
    } else if (lquery) {
#line 224 "ctzrzf.f"
	return 0;
#line 225 "ctzrzf.f"
    }

/*     Quick return if possible */

#line 229 "ctzrzf.f"
    if (*m == 0) {
#line 230 "ctzrzf.f"
	return 0;
#line 231 "ctzrzf.f"
    } else if (*m == *n) {
#line 232 "ctzrzf.f"
	i__1 = *n;
#line 232 "ctzrzf.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "ctzrzf.f"
	    i__2 = i__;
#line 233 "ctzrzf.f"
	    tau[i__2].r = 0., tau[i__2].i = 0.;
#line 234 "ctzrzf.f"
/* L10: */
#line 234 "ctzrzf.f"
	}
#line 235 "ctzrzf.f"
	return 0;
#line 236 "ctzrzf.f"
    }

#line 238 "ctzrzf.f"
    nbmin = 2;
#line 239 "ctzrzf.f"
    nx = 1;
#line 240 "ctzrzf.f"
    iws = *m;
#line 241 "ctzrzf.f"
    if (nb > 1 && nb < *m) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 245 "ctzrzf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "CGERQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 245 "ctzrzf.f"
	nx = max(i__1,i__2);
#line 246 "ctzrzf.f"
	if (nx < *m) {

/*           Determine if workspace is large enough for blocked code. */

#line 250 "ctzrzf.f"
	    ldwork = *m;
#line 251 "ctzrzf.f"
	    iws = ldwork * nb;
#line 252 "ctzrzf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 257 "ctzrzf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 258 "ctzrzf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "CGERQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 258 "ctzrzf.f"
		nbmin = max(i__1,i__2);
#line 260 "ctzrzf.f"
	    }
#line 261 "ctzrzf.f"
	}
#line 262 "ctzrzf.f"
    }

#line 264 "ctzrzf.f"
    if (nb >= nbmin && nb < *m && nx < *m) {

/*        Use blocked code initially. */
/*        The last kk rows are handled by the block method. */

/* Computing MIN */
#line 269 "ctzrzf.f"
	i__1 = *m + 1;
#line 269 "ctzrzf.f"
	m1 = min(i__1,*n);
#line 270 "ctzrzf.f"
	ki = (*m - nx - 1) / nb * nb;
/* Computing MIN */
#line 271 "ctzrzf.f"
	i__1 = *m, i__2 = ki + nb;
#line 271 "ctzrzf.f"
	kk = min(i__1,i__2);

#line 273 "ctzrzf.f"
	i__1 = *m - kk + 1;
#line 273 "ctzrzf.f"
	i__2 = -nb;
#line 273 "ctzrzf.f"
	for (i__ = *m - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; 
		i__ += i__2) {
/* Computing MIN */
#line 274 "ctzrzf.f"
	    i__3 = *m - i__ + 1;
#line 274 "ctzrzf.f"
	    ib = min(i__3,nb);

/*           Compute the TZ factorization of the current block */
/*           A(i:i+ib-1,i:n) */

#line 279 "ctzrzf.f"
	    i__3 = *n - i__ + 1;
#line 279 "ctzrzf.f"
	    i__4 = *n - *m;
#line 279 "ctzrzf.f"
	    clatrz_(&ib, &i__3, &i__4, &a[i__ + i__ * a_dim1], lda, &tau[i__],
		     &work[1]);
#line 281 "ctzrzf.f"
	    if (i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 286 "ctzrzf.f"
		i__3 = *n - *m;
#line 286 "ctzrzf.f"
		clarzt_("Backward", "Rowwise", &i__3, &ib, &a[i__ + m1 * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)8,
			 (ftnlen)7);

/*              Apply H to A(1:i-1,i:n) from the right */

#line 291 "ctzrzf.f"
		i__3 = i__ - 1;
#line 291 "ctzrzf.f"
		i__4 = *n - i__ + 1;
#line 291 "ctzrzf.f"
		i__5 = *n - *m;
#line 291 "ctzrzf.f"
		clarzb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
			 &i__4, &ib, &i__5, &a[i__ + m1 * a_dim1], lda, &work[
			1], &ldwork, &a[i__ * a_dim1 + 1], lda, &work[ib + 1],
			 &ldwork, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)7)
			;
#line 295 "ctzrzf.f"
	    }
#line 296 "ctzrzf.f"
/* L20: */
#line 296 "ctzrzf.f"
	}
#line 297 "ctzrzf.f"
	mu = i__ + nb - 1;
#line 298 "ctzrzf.f"
    } else {
#line 299 "ctzrzf.f"
	mu = *m;
#line 300 "ctzrzf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 304 "ctzrzf.f"
    if (mu > 0) {
#line 304 "ctzrzf.f"
	i__2 = *n - *m;
#line 304 "ctzrzf.f"
	clatrz_(&mu, n, &i__2, &a[a_offset], lda, &tau[1], &work[1]);
#line 304 "ctzrzf.f"
    }

#line 307 "ctzrzf.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 309 "ctzrzf.f"
    return 0;

/*     End of CTZRZF */

} /* ctzrzf_ */

