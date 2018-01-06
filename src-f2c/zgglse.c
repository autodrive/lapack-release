#line 1 "zgglse.f"
/* zgglse.f -- translated by f2c (version 20100827).
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

#line 1 "zgglse.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> ZGGLSE solves overdetermined or underdetermined systems for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGLSE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgglse.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgglse.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgglse.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( * ), D( * ), */
/*      $                   WORK( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGLSE solves the linear equality-constrained least squares (LSE) */
/* > problem: */
/* > */
/* >         minimize || c - A*x ||_2   subject to   B*x = d */
/* > */
/* > where A is an M-by-N matrix, B is a P-by-N matrix, c is a given */
/* > M-vector, and d is a given P-vector. It is assumed that */
/* > P <= N <= M+P, and */
/* > */
/* >          rank(B) = P and  rank( (A) ) = N. */
/* >                               ( (B) ) */
/* > */
/* > These conditions ensure that the LSE problem has a unique solution, */
/* > which is obtained using a generalized RQ factorization of the */
/* > matrices (B, A) given by */
/* > */
/* >    B = (0 R)*Q,   A = Z*T*Q. */
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
/* >          The number of columns of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >          The number of rows of the matrix B. 0 <= P <= N <= M+P. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the elements on and above the diagonal of the array */
/* >          contain the min(M,N)-by-N upper trapezoidal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On entry, the P-by-N matrix B. */
/* >          On exit, the upper triangle of the subarray B(1:P,N-P+1:N) */
/* >          contains the P-by-P upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (M) */
/* >          On entry, C contains the right hand side vector for the */
/* >          least squares part of the LSE problem. */
/* >          On exit, the residual sum of squares for the solution */
/* >          is given by the sum of squares of elements N-P+1 to M of */
/* >          vector C. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is COMPLEX*16 array, dimension (P) */
/* >          On entry, D contains the right hand side vector for the */
/* >          constrained equation. */
/* >          On exit, D is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (N) */
/* >          On exit, X is the solution of the LSE problem. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,M+N+P). */
/* >          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB, */
/* >          where NB is an upper bound for the optimal blocksizes for */
/* >          ZGEQRF, CGERQF, ZUNMQR and CUNMRQ. */
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
/* >          = 1:  the upper triangular factor R associated with B in the */
/* >                generalized RQ factorization of the pair (B, A) is */
/* >                singular, so that rank(B) < P; the least squares */
/* >                solution could not be computed. */
/* >          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor */
/* >                T associated with A in the generalized RQ factorization */
/* >                of the pair (B, A) is singular, so that */
/* >                rank( (A) ) < N; the least squares solution could not */
/* >                    ( (B) ) */
/* >                be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHERsolve */

/*  ===================================================================== */
/* Subroutine */ int zgglse_(integer *m, integer *n, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, doublecomplex *d__, doublecomplex *x, 
	doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer nb, mn, nr, nb1, nb2, nb3, nb4, lopt;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ztrmv_(char *, char *, 
	    char *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zggrqf_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;
    static integer lwkmin, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmrq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), ztrtrs_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 222 "zgglse.f"
    /* Parameter adjustments */
#line 222 "zgglse.f"
    a_dim1 = *lda;
#line 222 "zgglse.f"
    a_offset = 1 + a_dim1;
#line 222 "zgglse.f"
    a -= a_offset;
#line 222 "zgglse.f"
    b_dim1 = *ldb;
#line 222 "zgglse.f"
    b_offset = 1 + b_dim1;
#line 222 "zgglse.f"
    b -= b_offset;
#line 222 "zgglse.f"
    --c__;
#line 222 "zgglse.f"
    --d__;
#line 222 "zgglse.f"
    --x;
#line 222 "zgglse.f"
    --work;
#line 222 "zgglse.f"

#line 222 "zgglse.f"
    /* Function Body */
#line 222 "zgglse.f"
    *info = 0;
#line 223 "zgglse.f"
    mn = min(*m,*n);
#line 224 "zgglse.f"
    lquery = *lwork == -1;
#line 225 "zgglse.f"
    if (*m < 0) {
#line 226 "zgglse.f"
	*info = -1;
#line 227 "zgglse.f"
    } else if (*n < 0) {
#line 228 "zgglse.f"
	*info = -2;
#line 229 "zgglse.f"
    } else if (*p < 0 || *p > *n || *p < *n - *m) {
#line 230 "zgglse.f"
	*info = -3;
#line 231 "zgglse.f"
    } else if (*lda < max(1,*m)) {
#line 232 "zgglse.f"
	*info = -5;
#line 233 "zgglse.f"
    } else if (*ldb < max(1,*p)) {
#line 234 "zgglse.f"
	*info = -7;
#line 235 "zgglse.f"
    }

/*     Calculate workspace */

#line 239 "zgglse.f"
    if (*info == 0) {
#line 240 "zgglse.f"
	if (*n == 0) {
#line 241 "zgglse.f"
	    lwkmin = 1;
#line 242 "zgglse.f"
	    lwkopt = 1;
#line 243 "zgglse.f"
	} else {
#line 244 "zgglse.f"
	    nb1 = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 245 "zgglse.f"
	    nb2 = ilaenv_(&c__1, "ZGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 246 "zgglse.f"
	    nb3 = ilaenv_(&c__1, "ZUNMQR", " ", m, n, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 247 "zgglse.f"
	    nb4 = ilaenv_(&c__1, "ZUNMRQ", " ", m, n, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
/* Computing MAX */
#line 248 "zgglse.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 248 "zgglse.f"
	    nb = max(i__1,nb4);
#line 249 "zgglse.f"
	    lwkmin = *m + *n + *p;
#line 250 "zgglse.f"
	    lwkopt = *p + mn + max(*m,*n) * nb;
#line 251 "zgglse.f"
	}
#line 252 "zgglse.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 254 "zgglse.f"
	if (*lwork < lwkmin && ! lquery) {
#line 255 "zgglse.f"
	    *info = -12;
#line 256 "zgglse.f"
	}
#line 257 "zgglse.f"
    }

#line 259 "zgglse.f"
    if (*info != 0) {
#line 260 "zgglse.f"
	i__1 = -(*info);
#line 260 "zgglse.f"
	xerbla_("ZGGLSE", &i__1, (ftnlen)6);
#line 261 "zgglse.f"
	return 0;
#line 262 "zgglse.f"
    } else if (lquery) {
#line 263 "zgglse.f"
	return 0;
#line 264 "zgglse.f"
    }

/*     Quick return if possible */

#line 268 "zgglse.f"
    if (*n == 0) {
#line 268 "zgglse.f"
	return 0;
#line 268 "zgglse.f"
    }

/*     Compute the GRQ factorization of matrices B and A: */

/*            B*Q**H = (  0  T12 ) P   Z**H*A*Q**H = ( R11 R12 ) N-P */
/*                        N-P  P                     (  0  R22 ) M+P-N */
/*                                                      N-P  P */

/*     where T12 and R11 are upper triangular, and Q and Z are */
/*     unitary. */

#line 280 "zgglse.f"
    i__1 = *lwork - *p - mn;
#line 280 "zgglse.f"
    zggrqf_(p, m, n, &b[b_offset], ldb, &work[1], &a[a_offset], lda, &work[*p 
	    + 1], &work[*p + mn + 1], &i__1, info);
#line 282 "zgglse.f"
    i__1 = *p + mn + 1;
#line 282 "zgglse.f"
    lopt = (integer) work[i__1].r;

/*     Update c = Z**H *c = ( c1 ) N-P */
/*                       ( c2 ) M+P-N */

#line 287 "zgglse.f"
    i__1 = max(1,*m);
#line 287 "zgglse.f"
    i__2 = *lwork - *p - mn;
#line 287 "zgglse.f"
    zunmqr_("Left", "Conjugate Transpose", m, &c__1, &mn, &a[a_offset], lda, &
	    work[*p + 1], &c__[1], &i__1, &work[*p + mn + 1], &i__2, info, (
	    ftnlen)4, (ftnlen)19);
/* Computing MAX */
#line 290 "zgglse.f"
    i__3 = *p + mn + 1;
#line 290 "zgglse.f"
    i__1 = lopt, i__2 = (integer) work[i__3].r;
#line 290 "zgglse.f"
    lopt = max(i__1,i__2);

/*     Solve T12*x2 = d for x2 */

#line 294 "zgglse.f"
    if (*p > 0) {
#line 295 "zgglse.f"
	ztrtrs_("Upper", "No transpose", "Non-unit", p, &c__1, &b[(*n - *p + 
		1) * b_dim1 + 1], ldb, &d__[1], p, info, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);

#line 298 "zgglse.f"
	if (*info > 0) {
#line 299 "zgglse.f"
	    *info = 1;
#line 300 "zgglse.f"
	    return 0;
#line 301 "zgglse.f"
	}

/*        Put the solution in X */

#line 305 "zgglse.f"
	zcopy_(p, &d__[1], &c__1, &x[*n - *p + 1], &c__1);

/*        Update c1 */

#line 309 "zgglse.f"
	i__1 = *n - *p;
#line 309 "zgglse.f"
	z__1.r = -1., z__1.i = -0.;
#line 309 "zgglse.f"
	zgemv_("No transpose", &i__1, p, &z__1, &a[(*n - *p + 1) * a_dim1 + 1]
		, lda, &d__[1], &c__1, &c_b1, &c__[1], &c__1, (ftnlen)12);
#line 311 "zgglse.f"
    }

/*     Solve R11*x1 = c1 for x1 */

#line 315 "zgglse.f"
    if (*n > *p) {
#line 316 "zgglse.f"
	i__1 = *n - *p;
#line 316 "zgglse.f"
	i__2 = *n - *p;
#line 316 "zgglse.f"
	ztrtrs_("Upper", "No transpose", "Non-unit", &i__1, &c__1, &a[
		a_offset], lda, &c__[1], &i__2, info, (ftnlen)5, (ftnlen)12, (
		ftnlen)8);

#line 319 "zgglse.f"
	if (*info > 0) {
#line 320 "zgglse.f"
	    *info = 2;
#line 321 "zgglse.f"
	    return 0;
#line 322 "zgglse.f"
	}

/*        Put the solutions in X */

#line 326 "zgglse.f"
	i__1 = *n - *p;
#line 326 "zgglse.f"
	zcopy_(&i__1, &c__[1], &c__1, &x[1], &c__1);
#line 327 "zgglse.f"
    }

/*     Compute the residual vector: */

#line 331 "zgglse.f"
    if (*m < *n) {
#line 332 "zgglse.f"
	nr = *m + *p - *n;
#line 333 "zgglse.f"
	if (nr > 0) {
#line 333 "zgglse.f"
	    i__1 = *n - *m;
#line 333 "zgglse.f"
	    z__1.r = -1., z__1.i = -0.;
#line 333 "zgglse.f"
	    zgemv_("No transpose", &nr, &i__1, &z__1, &a[*n - *p + 1 + (*m + 
		    1) * a_dim1], lda, &d__[nr + 1], &c__1, &c_b1, &c__[*n - *
		    p + 1], &c__1, (ftnlen)12);
#line 333 "zgglse.f"
	}
#line 336 "zgglse.f"
    } else {
#line 337 "zgglse.f"
	nr = *p;
#line 338 "zgglse.f"
    }
#line 339 "zgglse.f"
    if (nr > 0) {
#line 340 "zgglse.f"
	ztrmv_("Upper", "No transpose", "Non unit", &nr, &a[*n - *p + 1 + (*n 
		- *p + 1) * a_dim1], lda, &d__[1], &c__1, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);
#line 342 "zgglse.f"
	z__1.r = -1., z__1.i = -0.;
#line 342 "zgglse.f"
	zaxpy_(&nr, &z__1, &d__[1], &c__1, &c__[*n - *p + 1], &c__1);
#line 343 "zgglse.f"
    }

/*     Backward transformation x = Q**H*x */

#line 347 "zgglse.f"
    i__1 = *lwork - *p - mn;
#line 347 "zgglse.f"
    zunmrq_("Left", "Conjugate Transpose", n, &c__1, p, &b[b_offset], ldb, &
	    work[1], &x[1], n, &work[*p + mn + 1], &i__1, info, (ftnlen)4, (
	    ftnlen)19);
/* Computing MAX */
#line 349 "zgglse.f"
    i__4 = *p + mn + 1;
#line 349 "zgglse.f"
    i__2 = lopt, i__3 = (integer) work[i__4].r;
#line 349 "zgglse.f"
    i__1 = *p + mn + max(i__2,i__3);
#line 349 "zgglse.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;

#line 351 "zgglse.f"
    return 0;

/*     End of ZGGLSE */

} /* zgglse_ */

