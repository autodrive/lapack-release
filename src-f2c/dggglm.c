#line 1 "dggglm.f"
/* dggglm.f -- translated by f2c (version 20100827).
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

#line 1 "dggglm.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b32 = -1.;
static doublereal c_b34 = 1.;

/* > \brief \b DGGGLM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGGLM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggglm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggglm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggglm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), */
/*      $                   X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGGLM solves a general Gauss-Markov linear model (GLM) problem: */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,M) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,P) */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, D is the left hand side of the GLM equation. */
/* >          On exit, D is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (M) */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension (P) */
/* > */
/* >          On exit, X and Y are the solutions of the GLM problem. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N+M+P). */
/* >          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB, */
/* >          where NB is an upper bound for the optimal blocksizes for */
/* >          DGEQRF, SGERQF, DORMQR and SORMRQ. */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int dggglm_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *d__, 
	doublereal *x, doublereal *y, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, nb, np, nb1, nb2, nb3, nb4, lopt;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dggqrf_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkmin;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dormrq_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int dtrtrs_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 227 "dggglm.f"
    /* Parameter adjustments */
#line 227 "dggglm.f"
    a_dim1 = *lda;
#line 227 "dggglm.f"
    a_offset = 1 + a_dim1;
#line 227 "dggglm.f"
    a -= a_offset;
#line 227 "dggglm.f"
    b_dim1 = *ldb;
#line 227 "dggglm.f"
    b_offset = 1 + b_dim1;
#line 227 "dggglm.f"
    b -= b_offset;
#line 227 "dggglm.f"
    --d__;
#line 227 "dggglm.f"
    --x;
#line 227 "dggglm.f"
    --y;
#line 227 "dggglm.f"
    --work;
#line 227 "dggglm.f"

#line 227 "dggglm.f"
    /* Function Body */
#line 227 "dggglm.f"
    *info = 0;
#line 228 "dggglm.f"
    np = min(*n,*p);
#line 229 "dggglm.f"
    lquery = *lwork == -1;
#line 230 "dggglm.f"
    if (*n < 0) {
#line 231 "dggglm.f"
	*info = -1;
#line 232 "dggglm.f"
    } else if (*m < 0 || *m > *n) {
#line 233 "dggglm.f"
	*info = -2;
#line 234 "dggglm.f"
    } else if (*p < 0 || *p < *n - *m) {
#line 235 "dggglm.f"
	*info = -3;
#line 236 "dggglm.f"
    } else if (*lda < max(1,*n)) {
#line 237 "dggglm.f"
	*info = -5;
#line 238 "dggglm.f"
    } else if (*ldb < max(1,*n)) {
#line 239 "dggglm.f"
	*info = -7;
#line 240 "dggglm.f"
    }

/*     Calculate workspace */

#line 244 "dggglm.f"
    if (*info == 0) {
#line 245 "dggglm.f"
	if (*n == 0) {
#line 246 "dggglm.f"
	    lwkmin = 1;
#line 247 "dggglm.f"
	    lwkopt = 1;
#line 248 "dggglm.f"
	} else {
#line 249 "dggglm.f"
	    nb1 = ilaenv_(&c__1, "DGEQRF", " ", n, m, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 250 "dggglm.f"
	    nb2 = ilaenv_(&c__1, "DGERQF", " ", n, m, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 251 "dggglm.f"
	    nb3 = ilaenv_(&c__1, "DORMQR", " ", n, m, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 252 "dggglm.f"
	    nb4 = ilaenv_(&c__1, "DORMRQ", " ", n, m, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
/* Computing MAX */
#line 253 "dggglm.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 253 "dggglm.f"
	    nb = max(i__1,nb4);
#line 254 "dggglm.f"
	    lwkmin = *m + *n + *p;
#line 255 "dggglm.f"
	    lwkopt = *m + np + max(*n,*p) * nb;
#line 256 "dggglm.f"
	}
#line 257 "dggglm.f"
	work[1] = (doublereal) lwkopt;

#line 259 "dggglm.f"
	if (*lwork < lwkmin && ! lquery) {
#line 260 "dggglm.f"
	    *info = -12;
#line 261 "dggglm.f"
	}
#line 262 "dggglm.f"
    }

#line 264 "dggglm.f"
    if (*info != 0) {
#line 265 "dggglm.f"
	i__1 = -(*info);
#line 265 "dggglm.f"
	xerbla_("DGGGLM", &i__1, (ftnlen)6);
#line 266 "dggglm.f"
	return 0;
#line 267 "dggglm.f"
    } else if (lquery) {
#line 268 "dggglm.f"
	return 0;
#line 269 "dggglm.f"
    }

/*     Quick return if possible */

#line 273 "dggglm.f"
    if (*n == 0) {
#line 273 "dggglm.f"
	return 0;
#line 273 "dggglm.f"
    }

/*     Compute the GQR factorization of matrices A and B: */

/*          Q**T*A = ( R11 ) M,    Q**T*B*Z**T = ( T11   T12 ) M */
/*                   (  0  ) N-M                 (  0    T22 ) N-M */
/*                      M                         M+P-N  N-M */

/*     where R11 and T22 are upper triangular, and Q and Z are */
/*     orthogonal. */

#line 285 "dggglm.f"
    i__1 = *lwork - *m - np;
#line 285 "dggglm.f"
    dggqrf_(n, m, p, &a[a_offset], lda, &work[1], &b[b_offset], ldb, &work[*m 
	    + 1], &work[*m + np + 1], &i__1, info);
#line 287 "dggglm.f"
    lopt = (integer) work[*m + np + 1];

/*     Update left-hand-side vector d = Q**T*d = ( d1 ) M */
/*                                               ( d2 ) N-M */

#line 292 "dggglm.f"
    i__1 = max(1,*n);
#line 292 "dggglm.f"
    i__2 = *lwork - *m - np;
#line 292 "dggglm.f"
    dormqr_("Left", "Transpose", n, &c__1, m, &a[a_offset], lda, &work[1], &
	    d__[1], &i__1, &work[*m + np + 1], &i__2, info, (ftnlen)4, (
	    ftnlen)9);
/* Computing MAX */
#line 294 "dggglm.f"
    i__1 = lopt, i__2 = (integer) work[*m + np + 1];
#line 294 "dggglm.f"
    lopt = max(i__1,i__2);

/*     Solve T22*y2 = d2 for y2 */

#line 298 "dggglm.f"
    if (*n > *m) {
#line 299 "dggglm.f"
	i__1 = *n - *m;
#line 299 "dggglm.f"
	i__2 = *n - *m;
#line 299 "dggglm.f"
	dtrtrs_("Upper", "No transpose", "Non unit", &i__1, &c__1, &b[*m + 1 
		+ (*m + *p - *n + 1) * b_dim1], ldb, &d__[*m + 1], &i__2, 
		info, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 302 "dggglm.f"
	if (*info > 0) {
#line 303 "dggglm.f"
	    *info = 1;
#line 304 "dggglm.f"
	    return 0;
#line 305 "dggglm.f"
	}

#line 307 "dggglm.f"
	i__1 = *n - *m;
#line 307 "dggglm.f"
	dcopy_(&i__1, &d__[*m + 1], &c__1, &y[*m + *p - *n + 1], &c__1);
#line 308 "dggglm.f"
    }

/*     Set y1 = 0 */

#line 312 "dggglm.f"
    i__1 = *m + *p - *n;
#line 312 "dggglm.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "dggglm.f"
	y[i__] = 0.;
#line 314 "dggglm.f"
/* L10: */
#line 314 "dggglm.f"
    }

/*     Update d1 = d1 - T12*y2 */

#line 318 "dggglm.f"
    i__1 = *n - *m;
#line 318 "dggglm.f"
    dgemv_("No transpose", m, &i__1, &c_b32, &b[(*m + *p - *n + 1) * b_dim1 + 
	    1], ldb, &y[*m + *p - *n + 1], &c__1, &c_b34, &d__[1], &c__1, (
	    ftnlen)12);

/*     Solve triangular system: R11*x = d1 */

#line 323 "dggglm.f"
    if (*m > 0) {
#line 324 "dggglm.f"
	dtrtrs_("Upper", "No Transpose", "Non unit", m, &c__1, &a[a_offset], 
		lda, &d__[1], m, info, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 327 "dggglm.f"
	if (*info > 0) {
#line 328 "dggglm.f"
	    *info = 2;
#line 329 "dggglm.f"
	    return 0;
#line 330 "dggglm.f"
	}

/*        Copy D to X */

#line 334 "dggglm.f"
	dcopy_(m, &d__[1], &c__1, &x[1], &c__1);
#line 335 "dggglm.f"
    }

/*     Backward transformation y = Z**T *y */

/* Computing MAX */
#line 339 "dggglm.f"
    i__1 = 1, i__2 = *n - *p + 1;
#line 339 "dggglm.f"
    i__3 = max(1,*p);
#line 339 "dggglm.f"
    i__4 = *lwork - *m - np;
#line 339 "dggglm.f"
    dormrq_("Left", "Transpose", p, &c__1, &np, &b[max(i__1,i__2) + b_dim1], 
	    ldb, &work[*m + 1], &y[1], &i__3, &work[*m + np + 1], &i__4, info,
	     (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 342 "dggglm.f"
    i__1 = lopt, i__2 = (integer) work[*m + np + 1];
#line 342 "dggglm.f"
    work[1] = (doublereal) (*m + np + max(i__1,i__2));

#line 344 "dggglm.f"
    return 0;

/*     End of DGGGLM */

} /* dggglm_ */

