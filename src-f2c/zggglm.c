#line 1 "zggglm.f"
/* zggglm.f -- translated by f2c (version 20100827).
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

#line 1 "zggglm.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief <b> ZGGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGGLM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggglm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggglm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggglm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), */
/*      $                   X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGGLM solves a general Gauss-Markov linear model (GLM) problem: */
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
/* >          A is COMPLEX*16 array, dimension (LDA,M) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,P) */
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
/* >          D is COMPLEX*16 array, dimension (N) */
/* >          On entry, D is the left hand side of the GLM equation. */
/* >          On exit, D is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (M) */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array, dimension (P) */
/* > */
/* >          On exit, X and Y are the solutions of the GLM problem. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N+M+P). */
/* >          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB, */
/* >          where NB is an upper bound for the optimal blocksizes for */
/* >          ZGEQRF, ZGERQF, ZUNMQR and ZUNMRQ. */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHEReigen */

/*  ===================================================================== */
/* Subroutine */ int zggglm_(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex 
	*work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, nb, np, nb1, nb2, nb3, nb4, lopt;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zggqrf_(integer *, integer *, integer *, 
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

#line 228 "zggglm.f"
    /* Parameter adjustments */
#line 228 "zggglm.f"
    a_dim1 = *lda;
#line 228 "zggglm.f"
    a_offset = 1 + a_dim1;
#line 228 "zggglm.f"
    a -= a_offset;
#line 228 "zggglm.f"
    b_dim1 = *ldb;
#line 228 "zggglm.f"
    b_offset = 1 + b_dim1;
#line 228 "zggglm.f"
    b -= b_offset;
#line 228 "zggglm.f"
    --d__;
#line 228 "zggglm.f"
    --x;
#line 228 "zggglm.f"
    --y;
#line 228 "zggglm.f"
    --work;
#line 228 "zggglm.f"

#line 228 "zggglm.f"
    /* Function Body */
#line 228 "zggglm.f"
    *info = 0;
#line 229 "zggglm.f"
    np = min(*n,*p);
#line 230 "zggglm.f"
    lquery = *lwork == -1;
#line 231 "zggglm.f"
    if (*n < 0) {
#line 232 "zggglm.f"
	*info = -1;
#line 233 "zggglm.f"
    } else if (*m < 0 || *m > *n) {
#line 234 "zggglm.f"
	*info = -2;
#line 235 "zggglm.f"
    } else if (*p < 0 || *p < *n - *m) {
#line 236 "zggglm.f"
	*info = -3;
#line 237 "zggglm.f"
    } else if (*lda < max(1,*n)) {
#line 238 "zggglm.f"
	*info = -5;
#line 239 "zggglm.f"
    } else if (*ldb < max(1,*n)) {
#line 240 "zggglm.f"
	*info = -7;
#line 241 "zggglm.f"
    }

/*     Calculate workspace */

#line 245 "zggglm.f"
    if (*info == 0) {
#line 246 "zggglm.f"
	if (*n == 0) {
#line 247 "zggglm.f"
	    lwkmin = 1;
#line 248 "zggglm.f"
	    lwkopt = 1;
#line 249 "zggglm.f"
	} else {
#line 250 "zggglm.f"
	    nb1 = ilaenv_(&c__1, "ZGEQRF", " ", n, m, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 251 "zggglm.f"
	    nb2 = ilaenv_(&c__1, "ZGERQF", " ", n, m, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 252 "zggglm.f"
	    nb3 = ilaenv_(&c__1, "ZUNMQR", " ", n, m, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 253 "zggglm.f"
	    nb4 = ilaenv_(&c__1, "ZUNMRQ", " ", n, m, p, &c_n1, (ftnlen)6, (
		    ftnlen)1);
/* Computing MAX */
#line 254 "zggglm.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 254 "zggglm.f"
	    nb = max(i__1,nb4);
#line 255 "zggglm.f"
	    lwkmin = *m + *n + *p;
#line 256 "zggglm.f"
	    lwkopt = *m + np + max(*n,*p) * nb;
#line 257 "zggglm.f"
	}
#line 258 "zggglm.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 260 "zggglm.f"
	if (*lwork < lwkmin && ! lquery) {
#line 261 "zggglm.f"
	    *info = -12;
#line 262 "zggglm.f"
	}
#line 263 "zggglm.f"
    }

#line 265 "zggglm.f"
    if (*info != 0) {
#line 266 "zggglm.f"
	i__1 = -(*info);
#line 266 "zggglm.f"
	xerbla_("ZGGGLM", &i__1, (ftnlen)6);
#line 267 "zggglm.f"
	return 0;
#line 268 "zggglm.f"
    } else if (lquery) {
#line 269 "zggglm.f"
	return 0;
#line 270 "zggglm.f"
    }

/*     Quick return if possible */

#line 274 "zggglm.f"
    if (*n == 0) {
#line 274 "zggglm.f"
	return 0;
#line 274 "zggglm.f"
    }

/*     Compute the GQR factorization of matrices A and B: */

/*          Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M */
/*                   (  0  ) N-M                 (  0    T22 ) N-M */
/*                      M                         M+P-N  N-M */

/*     where R11 and T22 are upper triangular, and Q and Z are */
/*     unitary. */

#line 286 "zggglm.f"
    i__1 = *lwork - *m - np;
#line 286 "zggglm.f"
    zggqrf_(n, m, p, &a[a_offset], lda, &work[1], &b[b_offset], ldb, &work[*m 
	    + 1], &work[*m + np + 1], &i__1, info);
#line 288 "zggglm.f"
    i__1 = *m + np + 1;
#line 288 "zggglm.f"
    lopt = (integer) work[i__1].r;

/*     Update left-hand-side vector d = Q**H*d = ( d1 ) M */
/*                                               ( d2 ) N-M */

#line 293 "zggglm.f"
    i__1 = max(1,*n);
#line 293 "zggglm.f"
    i__2 = *lwork - *m - np;
#line 293 "zggglm.f"
    zunmqr_("Left", "Conjugate transpose", n, &c__1, m, &a[a_offset], lda, &
	    work[1], &d__[1], &i__1, &work[*m + np + 1], &i__2, info, (ftnlen)
	    4, (ftnlen)19);
/* Computing MAX */
#line 295 "zggglm.f"
    i__3 = *m + np + 1;
#line 295 "zggglm.f"
    i__1 = lopt, i__2 = (integer) work[i__3].r;
#line 295 "zggglm.f"
    lopt = max(i__1,i__2);

/*     Solve T22*y2 = d2 for y2 */

#line 299 "zggglm.f"
    if (*n > *m) {
#line 300 "zggglm.f"
	i__1 = *n - *m;
#line 300 "zggglm.f"
	i__2 = *n - *m;
#line 300 "zggglm.f"
	ztrtrs_("Upper", "No transpose", "Non unit", &i__1, &c__1, &b[*m + 1 
		+ (*m + *p - *n + 1) * b_dim1], ldb, &d__[*m + 1], &i__2, 
		info, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 303 "zggglm.f"
	if (*info > 0) {
#line 304 "zggglm.f"
	    *info = 1;
#line 305 "zggglm.f"
	    return 0;
#line 306 "zggglm.f"
	}

#line 308 "zggglm.f"
	i__1 = *n - *m;
#line 308 "zggglm.f"
	zcopy_(&i__1, &d__[*m + 1], &c__1, &y[*m + *p - *n + 1], &c__1);
#line 309 "zggglm.f"
    }

/*     Set y1 = 0 */

#line 313 "zggglm.f"
    i__1 = *m + *p - *n;
#line 313 "zggglm.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 314 "zggglm.f"
	i__2 = i__;
#line 314 "zggglm.f"
	y[i__2].r = 0., y[i__2].i = 0.;
#line 315 "zggglm.f"
/* L10: */
#line 315 "zggglm.f"
    }

/*     Update d1 = d1 - T12*y2 */

#line 319 "zggglm.f"
    i__1 = *n - *m;
#line 319 "zggglm.f"
    z__1.r = -1., z__1.i = -0.;
#line 319 "zggglm.f"
    zgemv_("No transpose", m, &i__1, &z__1, &b[(*m + *p - *n + 1) * b_dim1 + 
	    1], ldb, &y[*m + *p - *n + 1], &c__1, &c_b2, &d__[1], &c__1, (
	    ftnlen)12);

/*     Solve triangular system: R11*x = d1 */

#line 324 "zggglm.f"
    if (*m > 0) {
#line 325 "zggglm.f"
	ztrtrs_("Upper", "No Transpose", "Non unit", m, &c__1, &a[a_offset], 
		lda, &d__[1], m, info, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 328 "zggglm.f"
	if (*info > 0) {
#line 329 "zggglm.f"
	    *info = 2;
#line 330 "zggglm.f"
	    return 0;
#line 331 "zggglm.f"
	}

/*        Copy D to X */

#line 335 "zggglm.f"
	zcopy_(m, &d__[1], &c__1, &x[1], &c__1);
#line 336 "zggglm.f"
    }

/*     Backward transformation y = Z**H *y */

/* Computing MAX */
#line 340 "zggglm.f"
    i__1 = 1, i__2 = *n - *p + 1;
#line 340 "zggglm.f"
    i__3 = max(1,*p);
#line 340 "zggglm.f"
    i__4 = *lwork - *m - np;
#line 340 "zggglm.f"
    zunmrq_("Left", "Conjugate transpose", p, &c__1, &np, &b[max(i__1,i__2) + 
	    b_dim1], ldb, &work[*m + 1], &y[1], &i__3, &work[*m + np + 1], &
	    i__4, info, (ftnlen)4, (ftnlen)19);
/* Computing MAX */
#line 343 "zggglm.f"
    i__4 = *m + np + 1;
#line 343 "zggglm.f"
    i__2 = lopt, i__3 = (integer) work[i__4].r;
#line 343 "zggglm.f"
    i__1 = *m + np + max(i__2,i__3);
#line 343 "zggglm.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;

#line 345 "zggglm.f"
    return 0;

/*     End of ZGGGLM */

} /* zggglm_ */

