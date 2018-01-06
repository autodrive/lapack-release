#line 1 "sgglse.f"
/* sgglse.f -- translated by f2c (version 20100827).
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

#line 1 "sgglse.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b31 = -1.;
static doublereal c_b33 = 1.;

/* > \brief <b> SGGLSE solves overdetermined or underdetermined systems for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGLSE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgglse.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgglse.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgglse.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), C( * ), D( * ), */
/*      $                   WORK( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGLSE solves the linear equality-constrained least squares (LSE) */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          B is REAL array, dimension (LDB,N) */
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
/* >          C is REAL array, dimension (M) */
/* >          On entry, C contains the right hand side vector for the */
/* >          least squares part of the LSE problem. */
/* >          On exit, the residual sum of squares for the solution */
/* >          is given by the sum of squares of elements N-P+1 to M of */
/* >          vector C. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (P) */
/* >          On entry, D contains the right hand side vector for the */
/* >          constrained equation. */
/* >          On exit, D is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (N) */
/* >          On exit, X is the solution of the LSE problem. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,M+N+P). */
/* >          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB, */
/* >          where NB is an upper bound for the optimal blocksizes for */
/* >          SGEQRF, SGERQF, SORMQR and SORMRQ. */
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

/* > \date December 2016 */

/* > \ingroup realOTHERsolve */

/*  ===================================================================== */
/* Subroutine */ int sgglse_(integer *m, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	doublereal *d__, doublereal *x, doublereal *work, integer *lwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer nb, mn, nr, nb1, nb2, nb3, nb4, lopt;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), scopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), saxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    , strmv_(char *, char *, char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), xerbla_(char 
	    *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sggrqf_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer lwkmin, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    sormrq_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), strtrs_(char 
	    *, char *, char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 222 "sgglse.f"
    /* Parameter adjustments */
#line 222 "sgglse.f"
    a_dim1 = *lda;
#line 222 "sgglse.f"
    a_offset = 1 + a_dim1;
#line 222 "sgglse.f"
    a -= a_offset;
#line 222 "sgglse.f"
    b_dim1 = *ldb;
#line 222 "sgglse.f"
    b_offset = 1 + b_dim1;
#line 222 "sgglse.f"
    b -= b_offset;
#line 222 "sgglse.f"
    --c__;
#line 222 "sgglse.f"
    --d__;
#line 222 "sgglse.f"
    --x;
#line 222 "sgglse.f"
    --work;
#line 222 "sgglse.f"

#line 222 "sgglse.f"
    /* Function Body */
#line 222 "sgglse.f"
    *info = 0;
#line 223 "sgglse.f"
    mn = min(*m,*n);
#line 224 "sgglse.f"
    lquery = *lwork == -1;
#line 225 "sgglse.f"
    if (*m < 0) {
#line 226 "sgglse.f"
	*info = -1;
#line 227 "sgglse.f"
    } else if (*n < 0) {
#line 228 "sgglse.f"
	*info = -2;
#line 229 "sgglse.f"
    } else if (*p < 0 || *p > *n || *p < *n - *m) {
#line 230 "sgglse.f"
	*info = -3;
#line 231 "sgglse.f"
    } else if (*lda < max(1,*m)) {
#line 232 "sgglse.f"
	*info = -5;
#line 233 "sgglse.f"
    } else if (*ldb < max(1,*p)) {
#line 234 "sgglse.f"
	*info = -7;
#line 235 "sgglse.f"
    }

/*     Calculate workspace */

#line 239 "sgglse.f"
    if (*info == 0) {
#line 240 "sgglse.f"
	if (*n == 0) {
#line 241 "sgglse.f"
	    lwkmin = 1;
#line 242 "sgglse.f"
	    lwkopt = 1;
#line 243 "sgglse.f"
	} else {
#line 244 "sgglse.f"
	    nb1 = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 245 "sgglse.f"
	    nb2 = ilaenv_(&c__1, "SGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 246 "sgglse.f"
	    nb3 = ilaenv_(&c__1, "SORMQR", " ", m, n, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 247 "sgglse.f"
	    nb4 = ilaenv_(&c__1, "SORMRQ", " ", m, n, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
/* Computing MAX */
#line 248 "sgglse.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 248 "sgglse.f"
	    nb = max(i__1,nb4);
#line 249 "sgglse.f"
	    lwkmin = *m + *n + *p;
#line 250 "sgglse.f"
	    lwkopt = *p + mn + max(*m,*n) * nb;
#line 251 "sgglse.f"
	}
#line 252 "sgglse.f"
	work[1] = (doublereal) lwkopt;

#line 254 "sgglse.f"
	if (*lwork < lwkmin && ! lquery) {
#line 255 "sgglse.f"
	    *info = -12;
#line 256 "sgglse.f"
	}
#line 257 "sgglse.f"
    }

#line 259 "sgglse.f"
    if (*info != 0) {
#line 260 "sgglse.f"
	i__1 = -(*info);
#line 260 "sgglse.f"
	xerbla_("SGGLSE", &i__1, (ftnlen)6);
#line 261 "sgglse.f"
	return 0;
#line 262 "sgglse.f"
    } else if (lquery) {
#line 263 "sgglse.f"
	return 0;
#line 264 "sgglse.f"
    }

/*     Quick return if possible */

#line 268 "sgglse.f"
    if (*n == 0) {
#line 268 "sgglse.f"
	return 0;
#line 268 "sgglse.f"
    }

/*     Compute the GRQ factorization of matrices B and A: */

/*            B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P */
/*                        N-P  P                     (  0  R22 ) M+P-N */
/*                                                      N-P  P */

/*     where T12 and R11 are upper triangular, and Q and Z are */
/*     orthogonal. */

#line 280 "sgglse.f"
    i__1 = *lwork - *p - mn;
#line 280 "sgglse.f"
    sggrqf_(p, m, n, &b[b_offset], ldb, &work[1], &a[a_offset], lda, &work[*p 
	    + 1], &work[*p + mn + 1], &i__1, info);
#line 282 "sgglse.f"
    lopt = (integer) work[*p + mn + 1];

/*     Update c = Z**T *c = ( c1 ) N-P */
/*                          ( c2 ) M+P-N */

#line 287 "sgglse.f"
    i__1 = max(1,*m);
#line 287 "sgglse.f"
    i__2 = *lwork - *p - mn;
#line 287 "sgglse.f"
    sormqr_("Left", "Transpose", m, &c__1, &mn, &a[a_offset], lda, &work[*p + 
	    1], &c__[1], &i__1, &work[*p + mn + 1], &i__2, info, (ftnlen)4, (
	    ftnlen)9);
/* Computing MAX */
#line 289 "sgglse.f"
    i__1 = lopt, i__2 = (integer) work[*p + mn + 1];
#line 289 "sgglse.f"
    lopt = max(i__1,i__2);

/*     Solve T12*x2 = d for x2 */

#line 293 "sgglse.f"
    if (*p > 0) {
#line 294 "sgglse.f"
	strtrs_("Upper", "No transpose", "Non-unit", p, &c__1, &b[(*n - *p + 
		1) * b_dim1 + 1], ldb, &d__[1], p, info, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);

#line 297 "sgglse.f"
	if (*info > 0) {
#line 298 "sgglse.f"
	    *info = 1;
#line 299 "sgglse.f"
	    return 0;
#line 300 "sgglse.f"
	}

/*        Put the solution in X */

#line 304 "sgglse.f"
	scopy_(p, &d__[1], &c__1, &x[*n - *p + 1], &c__1);

/*        Update c1 */

#line 308 "sgglse.f"
	i__1 = *n - *p;
#line 308 "sgglse.f"
	sgemv_("No transpose", &i__1, p, &c_b31, &a[(*n - *p + 1) * a_dim1 + 
		1], lda, &d__[1], &c__1, &c_b33, &c__[1], &c__1, (ftnlen)12);
#line 310 "sgglse.f"
    }

/*     Solve R11*x1 = c1 for x1 */

#line 314 "sgglse.f"
    if (*n > *p) {
#line 315 "sgglse.f"
	i__1 = *n - *p;
#line 315 "sgglse.f"
	i__2 = *n - *p;
#line 315 "sgglse.f"
	strtrs_("Upper", "No transpose", "Non-unit", &i__1, &c__1, &a[
		a_offset], lda, &c__[1], &i__2, info, (ftnlen)5, (ftnlen)12, (
		ftnlen)8);

#line 318 "sgglse.f"
	if (*info > 0) {
#line 319 "sgglse.f"
	    *info = 2;
#line 320 "sgglse.f"
	    return 0;
#line 321 "sgglse.f"
	}

/*        Put the solutions in X */

#line 325 "sgglse.f"
	i__1 = *n - *p;
#line 325 "sgglse.f"
	scopy_(&i__1, &c__[1], &c__1, &x[1], &c__1);
#line 326 "sgglse.f"
    }

/*     Compute the residual vector: */

#line 330 "sgglse.f"
    if (*m < *n) {
#line 331 "sgglse.f"
	nr = *m + *p - *n;
#line 332 "sgglse.f"
	if (nr > 0) {
#line 332 "sgglse.f"
	    i__1 = *n - *m;
#line 332 "sgglse.f"
	    sgemv_("No transpose", &nr, &i__1, &c_b31, &a[*n - *p + 1 + (*m + 
		    1) * a_dim1], lda, &d__[nr + 1], &c__1, &c_b33, &c__[*n - 
		    *p + 1], &c__1, (ftnlen)12);
#line 332 "sgglse.f"
	}
#line 335 "sgglse.f"
    } else {
#line 336 "sgglse.f"
	nr = *p;
#line 337 "sgglse.f"
    }
#line 338 "sgglse.f"
    if (nr > 0) {
#line 339 "sgglse.f"
	strmv_("Upper", "No transpose", "Non unit", &nr, &a[*n - *p + 1 + (*n 
		- *p + 1) * a_dim1], lda, &d__[1], &c__1, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);
#line 341 "sgglse.f"
	saxpy_(&nr, &c_b31, &d__[1], &c__1, &c__[*n - *p + 1], &c__1);
#line 342 "sgglse.f"
    }

/*     Backward transformation x = Q**T*x */

#line 346 "sgglse.f"
    i__1 = *lwork - *p - mn;
#line 346 "sgglse.f"
    sormrq_("Left", "Transpose", n, &c__1, p, &b[b_offset], ldb, &work[1], &x[
	    1], n, &work[*p + mn + 1], &i__1, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 348 "sgglse.f"
    i__1 = lopt, i__2 = (integer) work[*p + mn + 1];
#line 348 "sgglse.f"
    work[1] = (doublereal) (*p + mn + max(i__1,i__2));

#line 350 "sgglse.f"
    return 0;

/*     End of SGGLSE */

} /* sgglse_ */

