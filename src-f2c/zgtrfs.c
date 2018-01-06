#line 1 "zgtrfs.f"
/* zgtrfs.f -- translated by f2c (version 20100827).
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

#line 1 "zgtrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = -1.;
static doublereal c_b19 = 1.;
static doublecomplex c_b26 = {1.,0.};

/* > \brief \b ZGTRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgtrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgtrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgtrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, */
/*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX*16         B( LDB, * ), D( * ), DF( * ), DL( * ), */
/*      $                   DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is tridiagonal, and provides */
/* > error bounds and backward error estimates for the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX*16 array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DLF */
/* > \verbatim */
/* >          DLF is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A as computed by ZGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* >          DF is COMPLEX*16 array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DUF */
/* > \verbatim */
/* >          DUF is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX*16 array, dimension (N-2) */
/* >          The (n-2) elements of the second superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          The right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by ZGTTRS. */
/* >          On exit, the improved solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* >          FERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The estimated forward error bound for each solution vector */
/* >          X(j) (the j-th column of the solution matrix X). */
/* >          If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* >          is an estimated upper bound for the magnitude of the largest */
/* >          element in (X(j) - XTRUE) divided by the magnitude of the */
/* >          largest element in X(j).  The estimate is as reliable as */
/* >          the estimate for RCOND, and is almost always a slight */
/* >          overestimate of the true error. */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in */
/* >          any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  ITMAX is the maximum number of steps of iterative refinement. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GTcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgtrfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3], count;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zlacn2_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlagtm_(
	    char *, integer *, integer *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublereal *, doublecomplex *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
    static doublereal lstres;
    extern /* Subroutine */ int zgttrs_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, doublecomplex *, integer *, integer *, ftnlen);


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 273 "zgtrfs.f"
    /* Parameter adjustments */
#line 273 "zgtrfs.f"
    --dl;
#line 273 "zgtrfs.f"
    --d__;
#line 273 "zgtrfs.f"
    --du;
#line 273 "zgtrfs.f"
    --dlf;
#line 273 "zgtrfs.f"
    --df;
#line 273 "zgtrfs.f"
    --duf;
#line 273 "zgtrfs.f"
    --du2;
#line 273 "zgtrfs.f"
    --ipiv;
#line 273 "zgtrfs.f"
    b_dim1 = *ldb;
#line 273 "zgtrfs.f"
    b_offset = 1 + b_dim1;
#line 273 "zgtrfs.f"
    b -= b_offset;
#line 273 "zgtrfs.f"
    x_dim1 = *ldx;
#line 273 "zgtrfs.f"
    x_offset = 1 + x_dim1;
#line 273 "zgtrfs.f"
    x -= x_offset;
#line 273 "zgtrfs.f"
    --ferr;
#line 273 "zgtrfs.f"
    --berr;
#line 273 "zgtrfs.f"
    --work;
#line 273 "zgtrfs.f"
    --rwork;
#line 273 "zgtrfs.f"

#line 273 "zgtrfs.f"
    /* Function Body */
#line 273 "zgtrfs.f"
    *info = 0;
#line 274 "zgtrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 275 "zgtrfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 277 "zgtrfs.f"
	*info = -1;
#line 278 "zgtrfs.f"
    } else if (*n < 0) {
#line 279 "zgtrfs.f"
	*info = -2;
#line 280 "zgtrfs.f"
    } else if (*nrhs < 0) {
#line 281 "zgtrfs.f"
	*info = -3;
#line 282 "zgtrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 283 "zgtrfs.f"
	*info = -13;
#line 284 "zgtrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 285 "zgtrfs.f"
	*info = -15;
#line 286 "zgtrfs.f"
    }
#line 287 "zgtrfs.f"
    if (*info != 0) {
#line 288 "zgtrfs.f"
	i__1 = -(*info);
#line 288 "zgtrfs.f"
	xerbla_("ZGTRFS", &i__1, (ftnlen)6);
#line 289 "zgtrfs.f"
	return 0;
#line 290 "zgtrfs.f"
    }

/*     Quick return if possible */

#line 294 "zgtrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 295 "zgtrfs.f"
	i__1 = *nrhs;
#line 295 "zgtrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 296 "zgtrfs.f"
	    ferr[j] = 0.;
#line 297 "zgtrfs.f"
	    berr[j] = 0.;
#line 298 "zgtrfs.f"
/* L10: */
#line 298 "zgtrfs.f"
	}
#line 299 "zgtrfs.f"
	return 0;
#line 300 "zgtrfs.f"
    }

#line 302 "zgtrfs.f"
    if (notran) {
#line 303 "zgtrfs.f"
	*(unsigned char *)transn = 'N';
#line 304 "zgtrfs.f"
	*(unsigned char *)transt = 'C';
#line 305 "zgtrfs.f"
    } else {
#line 306 "zgtrfs.f"
	*(unsigned char *)transn = 'C';
#line 307 "zgtrfs.f"
	*(unsigned char *)transt = 'N';
#line 308 "zgtrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 312 "zgtrfs.f"
    nz = 4;
#line 313 "zgtrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 314 "zgtrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 315 "zgtrfs.f"
    safe1 = nz * safmin;
#line 316 "zgtrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 320 "zgtrfs.f"
    i__1 = *nrhs;
#line 320 "zgtrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 322 "zgtrfs.f"
	count = 1;
#line 323 "zgtrfs.f"
	lstres = 3.;
#line 324 "zgtrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 331 "zgtrfs.f"
	zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 332 "zgtrfs.f"
	zlagtm_(trans, n, &c__1, &c_b18, &dl[1], &d__[1], &du[1], &x[j * 
		x_dim1 + 1], ldx, &c_b19, &work[1], n, (ftnlen)1);

/*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward */
/*        error bound. */

#line 338 "zgtrfs.f"
	if (notran) {
#line 339 "zgtrfs.f"
	    if (*n == 1) {
#line 340 "zgtrfs.f"
		i__2 = j * b_dim1 + 1;
#line 340 "zgtrfs.f"
		i__3 = j * x_dim1 + 1;
#line 340 "zgtrfs.f"
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6)));
#line 342 "zgtrfs.f"
	    } else {
#line 343 "zgtrfs.f"
		i__2 = j * b_dim1 + 1;
#line 343 "zgtrfs.f"
		i__3 = j * x_dim1 + 1;
#line 343 "zgtrfs.f"
		i__4 = j * x_dim1 + 2;
#line 343 "zgtrfs.f"
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6))) + ((d__7 = du[1].r, abs(
			d__7)) + (d__8 = d_imag(&du[1]), abs(d__8))) * ((d__9 
			= x[i__4].r, abs(d__9)) + (d__10 = d_imag(&x[j * 
			x_dim1 + 2]), abs(d__10)));
#line 346 "zgtrfs.f"
		i__2 = *n - 1;
#line 346 "zgtrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 347 "zgtrfs.f"
		    i__3 = i__ + j * b_dim1;
#line 347 "zgtrfs.f"
		    i__4 = i__ - 1;
#line 347 "zgtrfs.f"
		    i__5 = i__ - 1 + j * x_dim1;
#line 347 "zgtrfs.f"
		    i__6 = i__;
#line 347 "zgtrfs.f"
		    i__7 = i__ + j * x_dim1;
#line 347 "zgtrfs.f"
		    i__8 = i__;
#line 347 "zgtrfs.f"
		    i__9 = i__ + 1 + j * x_dim1;
#line 347 "zgtrfs.f"
		    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&b[i__ + j * b_dim1]), abs(d__2)) + ((d__3 
			    = dl[i__4].r, abs(d__3)) + (d__4 = d_imag(&dl[i__ 
			    - 1]), abs(d__4))) * ((d__5 = x[i__5].r, abs(d__5)
			    ) + (d__6 = d_imag(&x[i__ - 1 + j * x_dim1]), abs(
			    d__6))) + ((d__7 = d__[i__6].r, abs(d__7)) + (
			    d__8 = d_imag(&d__[i__]), abs(d__8))) * ((d__9 = 
			    x[i__7].r, abs(d__9)) + (d__10 = d_imag(&x[i__ + 
			    j * x_dim1]), abs(d__10))) + ((d__11 = du[i__8].r,
			     abs(d__11)) + (d__12 = d_imag(&du[i__]), abs(
			    d__12))) * ((d__13 = x[i__9].r, abs(d__13)) + (
			    d__14 = d_imag(&x[i__ + 1 + j * x_dim1]), abs(
			    d__14)));
#line 351 "zgtrfs.f"
/* L30: */
#line 351 "zgtrfs.f"
		}
#line 352 "zgtrfs.f"
		i__2 = *n + j * b_dim1;
#line 352 "zgtrfs.f"
		i__3 = *n - 1;
#line 352 "zgtrfs.f"
		i__4 = *n - 1 + j * x_dim1;
#line 352 "zgtrfs.f"
		i__5 = *n;
#line 352 "zgtrfs.f"
		i__6 = *n + j * x_dim1;
#line 352 "zgtrfs.f"
		rwork[*n] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			*n + j * b_dim1]), abs(d__2)) + ((d__3 = dl[i__3].r, 
			abs(d__3)) + (d__4 = d_imag(&dl[*n - 1]), abs(d__4))) 
			* ((d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[*
			n - 1 + j * x_dim1]), abs(d__6))) + ((d__7 = d__[i__5]
			.r, abs(d__7)) + (d__8 = d_imag(&d__[*n]), abs(d__8)))
			 * ((d__9 = x[i__6].r, abs(d__9)) + (d__10 = d_imag(&
			x[*n + j * x_dim1]), abs(d__10)));
#line 355 "zgtrfs.f"
	    }
#line 356 "zgtrfs.f"
	} else {
#line 357 "zgtrfs.f"
	    if (*n == 1) {
#line 358 "zgtrfs.f"
		i__2 = j * b_dim1 + 1;
#line 358 "zgtrfs.f"
		i__3 = j * x_dim1 + 1;
#line 358 "zgtrfs.f"
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6)));
#line 360 "zgtrfs.f"
	    } else {
#line 361 "zgtrfs.f"
		i__2 = j * b_dim1 + 1;
#line 361 "zgtrfs.f"
		i__3 = j * x_dim1 + 1;
#line 361 "zgtrfs.f"
		i__4 = j * x_dim1 + 2;
#line 361 "zgtrfs.f"
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6))) + ((d__7 = dl[1].r, abs(
			d__7)) + (d__8 = d_imag(&dl[1]), abs(d__8))) * ((d__9 
			= x[i__4].r, abs(d__9)) + (d__10 = d_imag(&x[j * 
			x_dim1 + 2]), abs(d__10)));
#line 364 "zgtrfs.f"
		i__2 = *n - 1;
#line 364 "zgtrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 365 "zgtrfs.f"
		    i__3 = i__ + j * b_dim1;
#line 365 "zgtrfs.f"
		    i__4 = i__ - 1;
#line 365 "zgtrfs.f"
		    i__5 = i__ - 1 + j * x_dim1;
#line 365 "zgtrfs.f"
		    i__6 = i__;
#line 365 "zgtrfs.f"
		    i__7 = i__ + j * x_dim1;
#line 365 "zgtrfs.f"
		    i__8 = i__;
#line 365 "zgtrfs.f"
		    i__9 = i__ + 1 + j * x_dim1;
#line 365 "zgtrfs.f"
		    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&b[i__ + j * b_dim1]), abs(d__2)) + ((d__3 
			    = du[i__4].r, abs(d__3)) + (d__4 = d_imag(&du[i__ 
			    - 1]), abs(d__4))) * ((d__5 = x[i__5].r, abs(d__5)
			    ) + (d__6 = d_imag(&x[i__ - 1 + j * x_dim1]), abs(
			    d__6))) + ((d__7 = d__[i__6].r, abs(d__7)) + (
			    d__8 = d_imag(&d__[i__]), abs(d__8))) * ((d__9 = 
			    x[i__7].r, abs(d__9)) + (d__10 = d_imag(&x[i__ + 
			    j * x_dim1]), abs(d__10))) + ((d__11 = dl[i__8].r,
			     abs(d__11)) + (d__12 = d_imag(&dl[i__]), abs(
			    d__12))) * ((d__13 = x[i__9].r, abs(d__13)) + (
			    d__14 = d_imag(&x[i__ + 1 + j * x_dim1]), abs(
			    d__14)));
#line 369 "zgtrfs.f"
/* L40: */
#line 369 "zgtrfs.f"
		}
#line 370 "zgtrfs.f"
		i__2 = *n + j * b_dim1;
#line 370 "zgtrfs.f"
		i__3 = *n - 1;
#line 370 "zgtrfs.f"
		i__4 = *n - 1 + j * x_dim1;
#line 370 "zgtrfs.f"
		i__5 = *n;
#line 370 "zgtrfs.f"
		i__6 = *n + j * x_dim1;
#line 370 "zgtrfs.f"
		rwork[*n] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			*n + j * b_dim1]), abs(d__2)) + ((d__3 = du[i__3].r, 
			abs(d__3)) + (d__4 = d_imag(&du[*n - 1]), abs(d__4))) 
			* ((d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[*
			n - 1 + j * x_dim1]), abs(d__6))) + ((d__7 = d__[i__5]
			.r, abs(d__7)) + (d__8 = d_imag(&d__[*n]), abs(d__8)))
			 * ((d__9 = x[i__6].r, abs(d__9)) + (d__10 = d_imag(&
			x[*n + j * x_dim1]), abs(d__10)));
#line 373 "zgtrfs.f"
	    }
#line 374 "zgtrfs.f"
	}

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 385 "zgtrfs.f"
	s = 0.;
#line 386 "zgtrfs.f"
	i__2 = *n;
#line 386 "zgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 387 "zgtrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 388 "zgtrfs.f"
		i__3 = i__;
#line 388 "zgtrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 388 "zgtrfs.f"
		s = max(d__3,d__4);
#line 389 "zgtrfs.f"
	    } else {
/* Computing MAX */
#line 390 "zgtrfs.f"
		i__3 = i__;
#line 390 "zgtrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 390 "zgtrfs.f"
		s = max(d__3,d__4);
#line 392 "zgtrfs.f"
	    }
#line 393 "zgtrfs.f"
/* L50: */
#line 393 "zgtrfs.f"
	}
#line 394 "zgtrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 402 "zgtrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 407 "zgtrfs.f"
	    zgttrs_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[
		    1], &work[1], n, info, (ftnlen)1);
#line 409 "zgtrfs.f"
	    zaxpy_(n, &c_b26, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 410 "zgtrfs.f"
	    lstres = berr[j];
#line 411 "zgtrfs.f"
	    ++count;
#line 412 "zgtrfs.f"
	    goto L20;
#line 413 "zgtrfs.f"
	}

/*        Bound error from formula */

/*        norm(X - XTRUE) / norm(X) .le. FERR = */
/*        norm( abs(inv(op(A)))* */
/*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X) */

/*        where */
/*          norm(Z) is the magnitude of the largest component of Z */
/*          inv(op(A)) is the inverse of op(A) */
/*          abs(Z) is the componentwise absolute value of the matrix or */
/*             vector Z */
/*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
/*          EPS is machine epsilon */

/*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B)) */
/*        is incremented by SAFE1 if the i-th component of */
/*        abs(op(A))*abs(X) + abs(B) is less than SAFE2. */

/*        Use ZLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 437 "zgtrfs.f"
	i__2 = *n;
#line 437 "zgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 438 "zgtrfs.f"
	    if (rwork[i__] > safe2) {
#line 439 "zgtrfs.f"
		i__3 = i__;
#line 439 "zgtrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 440 "zgtrfs.f"
	    } else {
#line 441 "zgtrfs.f"
		i__3 = i__;
#line 441 "zgtrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 443 "zgtrfs.f"
	    }
#line 444 "zgtrfs.f"
/* L60: */
#line 444 "zgtrfs.f"
	}

#line 446 "zgtrfs.f"
	kase = 0;
#line 447 "zgtrfs.f"
L70:
#line 448 "zgtrfs.f"
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 449 "zgtrfs.f"
	if (kase != 0) {
#line 450 "zgtrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 454 "zgtrfs.f"
		zgttrs_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[1], n, info, (ftnlen)1);
#line 456 "zgtrfs.f"
		i__2 = *n;
#line 456 "zgtrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 457 "zgtrfs.f"
		    i__3 = i__;
#line 457 "zgtrfs.f"
		    i__4 = i__;
#line 457 "zgtrfs.f"
		    i__5 = i__;
#line 457 "zgtrfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 457 "zgtrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 458 "zgtrfs.f"
/* L80: */
#line 458 "zgtrfs.f"
		}
#line 459 "zgtrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 463 "zgtrfs.f"
		i__2 = *n;
#line 463 "zgtrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 464 "zgtrfs.f"
		    i__3 = i__;
#line 464 "zgtrfs.f"
		    i__4 = i__;
#line 464 "zgtrfs.f"
		    i__5 = i__;
#line 464 "zgtrfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 464 "zgtrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 465 "zgtrfs.f"
/* L90: */
#line 465 "zgtrfs.f"
		}
#line 466 "zgtrfs.f"
		zgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[1], n, info, (ftnlen)1);
#line 468 "zgtrfs.f"
	    }
#line 469 "zgtrfs.f"
	    goto L70;
#line 470 "zgtrfs.f"
	}

/*        Normalize error. */

#line 474 "zgtrfs.f"
	lstres = 0.;
#line 475 "zgtrfs.f"
	i__2 = *n;
#line 475 "zgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 476 "zgtrfs.f"
	    i__3 = i__ + j * x_dim1;
#line 476 "zgtrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 476 "zgtrfs.f"
	    lstres = max(d__3,d__4);
#line 477 "zgtrfs.f"
/* L100: */
#line 477 "zgtrfs.f"
	}
#line 478 "zgtrfs.f"
	if (lstres != 0.) {
#line 478 "zgtrfs.f"
	    ferr[j] /= lstres;
#line 478 "zgtrfs.f"
	}

#line 481 "zgtrfs.f"
/* L110: */
#line 481 "zgtrfs.f"
    }

#line 483 "zgtrfs.f"
    return 0;

/*     End of ZGTRFS */

} /* zgtrfs_ */

