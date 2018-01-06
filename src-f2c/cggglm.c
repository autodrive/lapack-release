#line 1 "cggglm.f"
/* cggglm.f -- translated by f2c (version 20100827).
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

#line 1 "cggglm.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CGGGLM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGGLM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggglm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggglm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggglm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), */
/*      $                   X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGGLM solves a general Gauss-Markov linear model (GLM) problem: */
/* > */
/* >         minimize || y ||_2   subject to   d = A*x + B*y */
/* >             x */
/* > */
/* > where A is an N-by-M matrix, B is an N-by-P matrix, and d is a */
/* > given N-vector. It is assumed that M <= N <= M+P, and */
/* > */
/* >            rank(A) = M    and    rank( A B ) = N. */
/* > */
/* > Under these assumptions, the constrained equation is always */
/* > consistent, and there is a unique solution x and a minimal 2-norm */
/* > solution y, which is obtained using a generalized QR factorization */
/* > of the matrices (A, B) given by */
/* > */
/* >    A = Q*(R),   B = Q*T*Z. */
/* >          (0) */
/* > */
/* > In particular, if matrix B is square nonsingular, then the problem */
/* > GLM is equivalent to the following weighted linear least squares */
/* > problem */
/* > */
/* >              minimize || inv(B)*(d-A*x) ||_2 */
/* >                  x */
/* > */
/* > where inv(B) denotes the inverse of B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of rows of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of columns of the matrix A.  0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >          The number of columns of the matrix B.  P >= N-M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,M) */
/* >          On entry, the N-by-M matrix A. */
/* >          On exit, the upper triangular part of the array A contains */
/* >          the M-by-M upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,P) */
/* >          On entry, the N-by-P matrix B. */
/* >          On exit, if N <= P, the upper triangle of the subarray */
/* >          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T; */
/* >          if N > P, the elements on and above the (N-P)th subdiagonal */
/* >          contain the N-by-P upper trapezoidal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          On entry, D is the left hand side of the GLM equation. */
/* >          On exit, D is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (M) */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension (P) */
/* > */
/* >          On exit, X and Y are the solutions of the GLM problem. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N+M+P). */
/* >          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB, */
/* >          where NB is an upper bound for the optimal blocksizes for */
/* >          CGEQRF, CGERQF, CUNMQR and CUNMRQ. */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          = 1:  the upper triangular factor R associated with A in the */
/* >                generalized QR factorization of the pair (A, B) is */
/* >                singular, so that rank(A) < M; the least squares */
/* >                solution could not be computed. */
/* >          = 2:  the bottom (N-M) by (N-M) part of the upper trapezoidal */
/* >                factor T associated with B in the generalized QR */
/* >                factorization of the pair (A, B) is singular, so that */
/* >                rank( A B ) < N; the least squares solution could not */
/* >                be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complexOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int cggglm_(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex 
	*work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, nb, np, nb1, nb2, nb3, nb4, lopt;
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), cggqrf_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkmin;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), cunmrq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int ctrtrs_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  =================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 228 "cggglm.f"
    /* Parameter adjustments */
#line 228 "cggglm.f"
    a_dim1 = *lda;
#line 228 "cggglm.f"
    a_offset = 1 + a_dim1;
#line 228 "cggglm.f"
    a -= a_offset;
#line 228 "cggglm.f"
    b_dim1 = *ldb;
#line 228 "cggglm.f"
    b_offset = 1 + b_dim1;
#line 228 "cggglm.f"
    b -= b_offset;
#line 228 "cggglm.f"
    --d__;
#line 228 "cggglm.f"
    --x;
#line 228 "cggglm.f"
    --y;
#line 228 "cggglm.f"
    --work;
#line 228 "cggglm.f"

#line 228 "cggglm.f"
    /* Function Body */
#line 228 "cggglm.f"
    *info = 0;
#line 229 "cggglm.f"
    np = min(*n,*p);
#line 230 "cggglm.f"
    lquery = *lwork == -1;
#line 231 "cggglm.f"
    if (*n < 0) {
#line 232 "cggglm.f"
	*info = -1;
#line 233 "cggglm.f"
    } else if (*m < 0 || *m > *n) {
#line 234 "cggglm.f"
	*info = -2;
#line 235 "cggglm.f"
    } else if (*p < 0 || *p < *n - *m) {
#line 236 "cggglm.f"
	*info = -3;
#line 237 "cggglm.f"
    } else if (*lda < max(1,*n)) {
#line 238 "cggglm.f"
	*info = -5;
#line 239 "cggglm.f"
    } else if (*ldb < max(1,*n)) {
#line 240 "cggglm.f"
	*info = -7;
#line 241 "cggglm.f"
    }

/*     Calculate workspace */

#line 245 "cggglm.f"
    if (*info == 0) {
#line 246 "cggglm.f"
	if (*n == 0) {
#line 247 "cggglm.f"
	    lwkmin = 1;
#line 248 "cggglm.f"
	    lwkopt = 1;
#line 249 "cggglm.f"
	} else {
#line 250 "cggglm.f"
	    nb1 = ilaenv_(&c__1, "CGEQRF", " ", n, m, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 251 "cggglm.f"
	    nb2 = ilaenv_(&c__1, "CGERQF", " ", n, m, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 252 "cggglm.f"
	    nb3 = ilaenv_(&c__1, "CUNMQR", " ", n, m, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 253 "cggglm.f"
	    nb4 = ilaenv_(&c__1, "CUNMRQ", " ", n, m, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
/* Computing MAX */
#line 254 "cggglm.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 254 "cggglm.f"
	    nb = max(i__1,nb4);
#line 255 "cggglm.f"
	    lwkmin = *m + *n + *p;
#line 256 "cggglm.f"
	    lwkopt = *m + np + max(*n,*p) * nb;
#line 257 "cggglm.f"
	}
#line 258 "cggglm.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 260 "cggglm.f"
	if (*lwork < lwkmin && ! lquery) {
#line 261 "cggglm.f"
	    *info = -12;
#line 262 "cggglm.f"
	}
#line 263 "cggglm.f"
    }

#line 265 "cggglm.f"
    if (*info != 0) {
#line 266 "cggglm.f"
	i__1 = -(*info);
#line 266 "cggglm.f"
	xerbla_("CGGGLM", &i__1, (ftnlen)6);
#line 267 "cggglm.f"
	return 0;
#line 268 "cggglm.f"
    } else if (lquery) {
#line 269 "cggglm.f"
	return 0;
#line 270 "cggglm.f"
    }

/*     Quick return if possible */

#line 274 "cggglm.f"
    if (*n == 0) {
#line 274 "cggglm.f"
	return 0;
#line 274 "cggglm.f"
    }

/*     Compute the GQR factorization of matrices A and B: */

/*          Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M */
/*                   (  0  ) N-M                 (  0    T22 ) N-M */
/*                      M                         M+P-N  N-M */

/*     where R11 and T22 are upper triangular, and Q and Z are */
/*     unitary. */

#line 286 "cggglm.f"
    i__1 = *lwork - *m - np;
#line 286 "cggglm.f"
    cggqrf_(n, m, p, &a[a_offset], lda, &work[1], &b[b_offset], ldb, &work[*m 
	    + 1], &work[*m + np + 1], &i__1, info);
#line 288 "cggglm.f"
    i__1 = *m + np + 1;
#line 288 "cggglm.f"
    lopt = (integer) work[i__1].r;

/*     Update left-hand-side vector d = Q**H*d = ( d1 ) M */
/*                                               ( d2 ) N-M */

#line 293 "cggglm.f"
    i__1 = max(1,*n);
#line 293 "cggglm.f"
    i__2 = *lwork - *m - np;
#line 293 "cggglm.f"
    cunmqr_("Left", "Conjugate transpose", n, &c__1, m, &a[a_offset], lda, &
	    work[1], &d__[1], &i__1, &work[*m + np + 1], &i__2, info, (ftnlen)
	    4, (ftnlen)19);
/* Computing MAX */
#line 295 "cggglm.f"
    i__3 = *m + np + 1;
#line 295 "cggglm.f"
    i__1 = lopt, i__2 = (integer) work[i__3].r;
#line 295 "cggglm.f"
    lopt = max(i__1,i__2);

/*     Solve T22*y2 = d2 for y2 */

#line 299 "cggglm.f"
    if (*n > *m) {
#line 300 "cggglm.f"
	i__1 = *n - *m;
#line 300 "cggglm.f"
	i__2 = *n - *m;
#line 300 "cggglm.f"
	ctrtrs_("Upper", "No transpose", "Non unit", &i__1, &c__1, &b[*m + 1 
		+ (*m + *p - *n + 1) * b_dim1], ldb, &d__[*m + 1], &i__2, 
		info, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 303 "cggglm.f"
	if (*info > 0) {
#line 304 "cggglm.f"
	    *info = 1;
#line 305 "cggglm.f"
	    return 0;
#line 306 "cggglm.f"
	}

#line 308 "cggglm.f"
	i__1 = *n - *m;
#line 308 "cggglm.f"
	ccopy_(&i__1, &d__[*m + 1], &c__1, &y[*m + *p - *n + 1], &c__1);
#line 309 "cggglm.f"
    }

/*     Set y1 = 0 */

#line 313 "cggglm.f"
    i__1 = *m + *p - *n;
#line 313 "cggglm.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 314 "cggglm.f"
	i__2 = i__;
#line 314 "cggglm.f"
	y[i__2].r = 0., y[i__2].i = 0.;
#line 315 "cggglm.f"
/* L10: */
#line 315 "cggglm.f"
    }

/*     Update d1 = d1 - T12*y2 */

#line 319 "cggglm.f"
    i__1 = *n - *m;
#line 319 "cggglm.f"
    z__1.r = -1., z__1.i = -0.;
#line 319 "cggglm.f"
    cgemv_("No transpose", m, &i__1, &z__1, &b[(*m + *p - *n + 1) * b_dim1 + 
	    1], ldb, &y[*m + *p - *n + 1], &c__1, &c_b2, &d__[1], &c__1, (
	    ftnlen)12);

/*     Solve triangular system: R11*x = d1 */

#line 324 "cggglm.f"
    if (*m > 0) {
#line 325 "cggglm.f"
	ctrtrs_("Upper", "No Transpose", "Non unit", m, &c__1, &a[a_offset], 
		lda, &d__[1], m, info, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 328 "cggglm.f"
	if (*info > 0) {
#line 329 "cggglm.f"
	    *info = 2;
#line 330 "cggglm.f"
	    return 0;
#line 331 "cggglm.f"
	}

/*        Copy D to X */

#line 335 "cggglm.f"
	ccopy_(m, &d__[1], &c__1, &x[1], &c__1);
#line 336 "cggglm.f"
    }

/*     Backward transformation y = Z**H *y */

/* Computing MAX */
#line 340 "cggglm.f"
    i__1 = 1, i__2 = *n - *p + 1;
#line 340 "cggglm.f"
    i__3 = max(1,*p);
#line 340 "cggglm.f"
    i__4 = *lwork - *m - np;
#line 340 "cggglm.f"
    cunmrq_("Left", "Conjugate transpose", p, &c__1, &np, &b[max(i__1,i__2) + 
	    b_dim1], ldb, &work[*m + 1], &y[1], &i__3, &work[*m + np + 1], &
	    i__4, info, (ftnlen)4, (ftnlen)19);
/* Computing MAX */
#line 343 "cggglm.f"
    i__4 = *m + np + 1;
#line 343 "cggglm.f"
    i__2 = lopt, i__3 = (integer) work[i__4].r;
#line 343 "cggglm.f"
    i__1 = *m + np + max(i__2,i__3);
#line 343 "cggglm.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;

#line 345 "cggglm.f"
    return 0;

/*     End of CGGGLM */

} /* cggglm_ */

