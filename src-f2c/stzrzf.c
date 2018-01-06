#line 1 "stzrzf.f"
/* stzrzf.f -- translated by f2c (version 20100827).
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

#line 1 "stzrzf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b STZRZF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STZRZF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stzrzf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stzrzf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stzrzf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A */
/* > to upper triangular form by means of orthogonal transformations. */
/* > */
/* > The upper trapezoidal matrix A is factored as */
/* > */
/* >    A = ( R  0 ) * Z, */
/* > */
/* > where Z is an N-by-N orthogonal matrix and R is an M-by-M upper */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the leading M-by-N upper trapezoidal part of the */
/* >          array A must contain the matrix to be factorized. */
/* >          On exit, the leading M-by-M upper triangular part of A */
/* >          contains the upper triangular matrix R, and elements M+1 to */
/* >          N of the first M rows of A, with the array TAU, represent the */
/* >          orthogonal matrix Z as a product of M elementary reflectors. */
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
/* >          TAU is REAL array, dimension (M) */
/* >          The scalar factors of the elementary reflectors. */
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

/* > \ingroup realOTHERcomputational */

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
/* >     Z(k) = I - tau(k)*v(k)*v(k)**T */
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
/* Subroutine */ int stzrzf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, m1, ib, nb, ki, kk, mu, nx, iws, nbmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarzb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer lwkmin, ldwork;
    extern /* Subroutine */ int slarzt_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int slatrz_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);


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

#line 191 "stzrzf.f"
    /* Parameter adjustments */
#line 191 "stzrzf.f"
    a_dim1 = *lda;
#line 191 "stzrzf.f"
    a_offset = 1 + a_dim1;
#line 191 "stzrzf.f"
    a -= a_offset;
#line 191 "stzrzf.f"
    --tau;
#line 191 "stzrzf.f"
    --work;
#line 191 "stzrzf.f"

#line 191 "stzrzf.f"
    /* Function Body */
#line 191 "stzrzf.f"
    *info = 0;
#line 192 "stzrzf.f"
    lquery = *lwork == -1;
#line 193 "stzrzf.f"
    if (*m < 0) {
#line 194 "stzrzf.f"
	*info = -1;
#line 195 "stzrzf.f"
    } else if (*n < *m) {
#line 196 "stzrzf.f"
	*info = -2;
#line 197 "stzrzf.f"
    } else if (*lda < max(1,*m)) {
#line 198 "stzrzf.f"
	*info = -4;
#line 199 "stzrzf.f"
    }

#line 201 "stzrzf.f"
    if (*info == 0) {
#line 202 "stzrzf.f"
	if (*m == 0 || *m == *n) {
#line 203 "stzrzf.f"
	    lwkopt = 1;
#line 204 "stzrzf.f"
	    lwkmin = 1;
#line 205 "stzrzf.f"
	} else {

/*           Determine the block size. */

#line 209 "stzrzf.f"
	    nb = ilaenv_(&c__1, "SGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 210 "stzrzf.f"
	    lwkopt = *m * nb;
#line 211 "stzrzf.f"
	    lwkmin = max(1,*m);
#line 212 "stzrzf.f"
	}
#line 213 "stzrzf.f"
	work[1] = (doublereal) lwkopt;

#line 215 "stzrzf.f"
	if (*lwork < lwkmin && ! lquery) {
#line 216 "stzrzf.f"
	    *info = -7;
#line 217 "stzrzf.f"
	}
#line 218 "stzrzf.f"
    }

#line 220 "stzrzf.f"
    if (*info != 0) {
#line 221 "stzrzf.f"
	i__1 = -(*info);
#line 221 "stzrzf.f"
	xerbla_("STZRZF", &i__1, (ftnlen)6);
#line 222 "stzrzf.f"
	return 0;
#line 223 "stzrzf.f"
    } else if (lquery) {
#line 224 "stzrzf.f"
	return 0;
#line 225 "stzrzf.f"
    }

/*     Quick return if possible */

#line 229 "stzrzf.f"
    if (*m == 0) {
#line 230 "stzrzf.f"
	return 0;
#line 231 "stzrzf.f"
    } else if (*m == *n) {
#line 232 "stzrzf.f"
	i__1 = *n;
#line 232 "stzrzf.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "stzrzf.f"
	    tau[i__] = 0.;
#line 234 "stzrzf.f"
/* L10: */
#line 234 "stzrzf.f"
	}
#line 235 "stzrzf.f"
	return 0;
#line 236 "stzrzf.f"
    }

#line 238 "stzrzf.f"
    nbmin = 2;
#line 239 "stzrzf.f"
    nx = 1;
#line 240 "stzrzf.f"
    iws = *m;
#line 241 "stzrzf.f"
    if (nb > 1 && nb < *m) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 245 "stzrzf.f"
	i__1 = 0, i__2 = ilaenv_(&c__3, "SGERQF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 245 "stzrzf.f"
	nx = max(i__1,i__2);
#line 246 "stzrzf.f"
	if (nx < *m) {

/*           Determine if workspace is large enough for blocked code. */

#line 250 "stzrzf.f"
	    ldwork = *m;
#line 251 "stzrzf.f"
	    iws = ldwork * nb;
#line 252 "stzrzf.f"
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 257 "stzrzf.f"
		nb = *lwork / ldwork;
/* Computing MAX */
#line 258 "stzrzf.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGERQF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 258 "stzrzf.f"
		nbmin = max(i__1,i__2);
#line 260 "stzrzf.f"
	    }
#line 261 "stzrzf.f"
	}
#line 262 "stzrzf.f"
    }

#line 264 "stzrzf.f"
    if (nb >= nbmin && nb < *m && nx < *m) {

/*        Use blocked code initially. */
/*        The last kk rows are handled by the block method. */

/* Computing MIN */
#line 269 "stzrzf.f"
	i__1 = *m + 1;
#line 269 "stzrzf.f"
	m1 = min(i__1,*n);
#line 270 "stzrzf.f"
	ki = (*m - nx - 1) / nb * nb;
/* Computing MIN */
#line 271 "stzrzf.f"
	i__1 = *m, i__2 = ki + nb;
#line 271 "stzrzf.f"
	kk = min(i__1,i__2);

#line 273 "stzrzf.f"
	i__1 = *m - kk + 1;
#line 273 "stzrzf.f"
	i__2 = -nb;
#line 273 "stzrzf.f"
	for (i__ = *m - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; 
		i__ += i__2) {
/* Computing MIN */
#line 274 "stzrzf.f"
	    i__3 = *m - i__ + 1;
#line 274 "stzrzf.f"
	    ib = min(i__3,nb);

/*           Compute the TZ factorization of the current block */
/*           A(i:i+ib-1,i:n) */

#line 279 "stzrzf.f"
	    i__3 = *n - i__ + 1;
#line 279 "stzrzf.f"
	    i__4 = *n - *m;
#line 279 "stzrzf.f"
	    slatrz_(&ib, &i__3, &i__4, &a[i__ + i__ * a_dim1], lda, &tau[i__],
		     &work[1]);
#line 281 "stzrzf.f"
	    if (i__ > 1) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i+ib-1) . . . H(i+1) H(i) */

#line 286 "stzrzf.f"
		i__3 = *n - *m;
#line 286 "stzrzf.f"
		slarzt_("Backward", "Rowwise", &i__3, &ib, &a[i__ + m1 * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)8,
			 (ftnlen)7);

/*              Apply H to A(1:i-1,i:n) from the right */

#line 291 "stzrzf.f"
		i__3 = i__ - 1;
#line 291 "stzrzf.f"
		i__4 = *n - i__ + 1;
#line 291 "stzrzf.f"
		i__5 = *n - *m;
#line 291 "stzrzf.f"
		slarzb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
			 &i__4, &ib, &i__5, &a[i__ + m1 * a_dim1], lda, &work[
			1], &ldwork, &a[i__ * a_dim1 + 1], lda, &work[ib + 1],
			 &ldwork, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)7)
			;
#line 295 "stzrzf.f"
	    }
#line 296 "stzrzf.f"
/* L20: */
#line 296 "stzrzf.f"
	}
#line 297 "stzrzf.f"
	mu = i__ + nb - 1;
#line 298 "stzrzf.f"
    } else {
#line 299 "stzrzf.f"
	mu = *m;
#line 300 "stzrzf.f"
    }

/*     Use unblocked code to factor the last or only block */

#line 304 "stzrzf.f"
    if (mu > 0) {
#line 304 "stzrzf.f"
	i__2 = *n - *m;
#line 304 "stzrzf.f"
	slatrz_(&mu, n, &i__2, &a[a_offset], lda, &tau[1], &work[1]);
#line 304 "stzrzf.f"
    }

#line 307 "stzrzf.f"
    work[1] = (doublereal) lwkopt;

#line 309 "stzrzf.f"
    return 0;

/*     End of STZRZF */

} /* stzrzf_ */

