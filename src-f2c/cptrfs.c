#line 1 "cptrfs.f"
/* cptrfs.f -- translated by f2c (version 20100827).
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

#line 1 "cptrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b16 = {1.,0.};

/* > \brief \b CPTRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, */
/*                          FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), D( * ), DF( * ), FERR( * ), */
/*      $                   RWORK( * ) */
/*       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian positive definite */
/* > and tridiagonal, and provides error bounds and backward error */
/* > estimates for the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the superdiagonal or the subdiagonal of the */
/* >          tridiagonal matrix A is stored and the form of the */
/* >          factorization: */
/* >          = 'U':  E is the superdiagonal of A, and A = U**H*D*U; */
/* >          = 'L':  E is the subdiagonal of A, and A = L*D*L**H. */
/* >          (The two forms are equivalent if A is real.) */
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
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n real diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N-1) */
/* >          The (n-1) off-diagonal elements of the tridiagonal matrix A */
/* >          (see UPLO). */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* >          DF is REAL array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from */
/* >          the factorization computed by CPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] EF */
/* > \verbatim */
/* >          EF is COMPLEX array, dimension (N-1) */
/* >          The (n-1) off-diagonal elements of the unit bidiagonal */
/* >          factor U or L from the factorization computed by CPTTRF */
/* >          (see UPLO). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by CPTTRS. */
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
/* >          FERR is REAL array, dimension (NRHS) */
/* >          The forward error bound for each solution vector */
/* >          X(j) (the j-th column of the solution matrix X). */
/* >          If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* >          is an estimated upper bound for the magnitude of the largest */
/* >          element in (X(j) - XTRUE) divided by the magnitude of the */
/* >          largest element in X(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is REAL array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in */
/* >          any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
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

/* > \ingroup complexPTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cptrfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static doublecomplex bi, cx, dx, ex;
    static integer ix, nz;
    static doublereal eps, safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer count;
    static logical upper;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal lstres;
    extern /* Subroutine */ int cpttrs_(char *, integer *, integer *, 
	    doublereal *, doublecomplex *, doublecomplex *, integer *, 
	    integer *, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 244 "cptrfs.f"
    /* Parameter adjustments */
#line 244 "cptrfs.f"
    --d__;
#line 244 "cptrfs.f"
    --e;
#line 244 "cptrfs.f"
    --df;
#line 244 "cptrfs.f"
    --ef;
#line 244 "cptrfs.f"
    b_dim1 = *ldb;
#line 244 "cptrfs.f"
    b_offset = 1 + b_dim1;
#line 244 "cptrfs.f"
    b -= b_offset;
#line 244 "cptrfs.f"
    x_dim1 = *ldx;
#line 244 "cptrfs.f"
    x_offset = 1 + x_dim1;
#line 244 "cptrfs.f"
    x -= x_offset;
#line 244 "cptrfs.f"
    --ferr;
#line 244 "cptrfs.f"
    --berr;
#line 244 "cptrfs.f"
    --work;
#line 244 "cptrfs.f"
    --rwork;
#line 244 "cptrfs.f"

#line 244 "cptrfs.f"
    /* Function Body */
#line 244 "cptrfs.f"
    *info = 0;
#line 245 "cptrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 246 "cptrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 247 "cptrfs.f"
	*info = -1;
#line 248 "cptrfs.f"
    } else if (*n < 0) {
#line 249 "cptrfs.f"
	*info = -2;
#line 250 "cptrfs.f"
    } else if (*nrhs < 0) {
#line 251 "cptrfs.f"
	*info = -3;
#line 252 "cptrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 253 "cptrfs.f"
	*info = -9;
#line 254 "cptrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 255 "cptrfs.f"
	*info = -11;
#line 256 "cptrfs.f"
    }
#line 257 "cptrfs.f"
    if (*info != 0) {
#line 258 "cptrfs.f"
	i__1 = -(*info);
#line 258 "cptrfs.f"
	xerbla_("CPTRFS", &i__1, (ftnlen)6);
#line 259 "cptrfs.f"
	return 0;
#line 260 "cptrfs.f"
    }

/*     Quick return if possible */

#line 264 "cptrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 265 "cptrfs.f"
	i__1 = *nrhs;
#line 265 "cptrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 266 "cptrfs.f"
	    ferr[j] = 0.;
#line 267 "cptrfs.f"
	    berr[j] = 0.;
#line 268 "cptrfs.f"
/* L10: */
#line 268 "cptrfs.f"
	}
#line 269 "cptrfs.f"
	return 0;
#line 270 "cptrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 274 "cptrfs.f"
    nz = 4;
#line 275 "cptrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 276 "cptrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 277 "cptrfs.f"
    safe1 = nz * safmin;
#line 278 "cptrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 282 "cptrfs.f"
    i__1 = *nrhs;
#line 282 "cptrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 284 "cptrfs.f"
	count = 1;
#line 285 "cptrfs.f"
	lstres = 3.;
#line 286 "cptrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X.  Also compute */
/*        abs(A)*abs(x) + abs(b) for use in the backward error bound. */

#line 293 "cptrfs.f"
	if (upper) {
#line 294 "cptrfs.f"
	    if (*n == 1) {
#line 295 "cptrfs.f"
		i__2 = j * b_dim1 + 1;
#line 295 "cptrfs.f"
		bi.r = b[i__2].r, bi.i = b[i__2].i;
#line 296 "cptrfs.f"
		i__2 = j * x_dim1 + 1;
#line 296 "cptrfs.f"
		z__1.r = d__[1] * x[i__2].r, z__1.i = d__[1] * x[i__2].i;
#line 296 "cptrfs.f"
		dx.r = z__1.r, dx.i = z__1.i;
#line 297 "cptrfs.f"
		z__1.r = bi.r - dx.r, z__1.i = bi.i - dx.i;
#line 297 "cptrfs.f"
		work[1].r = z__1.r, work[1].i = z__1.i;
#line 298 "cptrfs.f"
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4)));
#line 299 "cptrfs.f"
	    } else {
#line 300 "cptrfs.f"
		i__2 = j * b_dim1 + 1;
#line 300 "cptrfs.f"
		bi.r = b[i__2].r, bi.i = b[i__2].i;
#line 301 "cptrfs.f"
		i__2 = j * x_dim1 + 1;
#line 301 "cptrfs.f"
		z__1.r = d__[1] * x[i__2].r, z__1.i = d__[1] * x[i__2].i;
#line 301 "cptrfs.f"
		dx.r = z__1.r, dx.i = z__1.i;
#line 302 "cptrfs.f"
		i__2 = j * x_dim1 + 2;
#line 302 "cptrfs.f"
		z__1.r = e[1].r * x[i__2].r - e[1].i * x[i__2].i, z__1.i = e[
			1].r * x[i__2].i + e[1].i * x[i__2].r;
#line 302 "cptrfs.f"
		ex.r = z__1.r, ex.i = z__1.i;
#line 303 "cptrfs.f"
		z__2.r = bi.r - dx.r, z__2.i = bi.i - dx.i;
#line 303 "cptrfs.f"
		z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
#line 303 "cptrfs.f"
		work[1].r = z__1.r, work[1].i = z__1.i;
#line 304 "cptrfs.f"
		i__2 = j * x_dim1 + 2;
#line 304 "cptrfs.f"
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4))) + ((d__5 = e[1].r, abs(d__5))
			 + (d__6 = d_imag(&e[1]), abs(d__6))) * ((d__7 = x[
			i__2].r, abs(d__7)) + (d__8 = d_imag(&x[j * x_dim1 + 
			2]), abs(d__8)));
#line 306 "cptrfs.f"
		i__2 = *n - 1;
#line 306 "cptrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 307 "cptrfs.f"
		    i__3 = i__ + j * b_dim1;
#line 307 "cptrfs.f"
		    bi.r = b[i__3].r, bi.i = b[i__3].i;
#line 308 "cptrfs.f"
		    d_cnjg(&z__2, &e[i__ - 1]);
#line 308 "cptrfs.f"
		    i__3 = i__ - 1 + j * x_dim1;
#line 308 "cptrfs.f"
		    z__1.r = z__2.r * x[i__3].r - z__2.i * x[i__3].i, z__1.i =
			     z__2.r * x[i__3].i + z__2.i * x[i__3].r;
#line 308 "cptrfs.f"
		    cx.r = z__1.r, cx.i = z__1.i;
#line 309 "cptrfs.f"
		    i__3 = i__;
#line 309 "cptrfs.f"
		    i__4 = i__ + j * x_dim1;
#line 309 "cptrfs.f"
		    z__1.r = d__[i__3] * x[i__4].r, z__1.i = d__[i__3] * x[
			    i__4].i;
#line 309 "cptrfs.f"
		    dx.r = z__1.r, dx.i = z__1.i;
#line 310 "cptrfs.f"
		    i__3 = i__;
#line 310 "cptrfs.f"
		    i__4 = i__ + 1 + j * x_dim1;
#line 310 "cptrfs.f"
		    z__1.r = e[i__3].r * x[i__4].r - e[i__3].i * x[i__4].i, 
			    z__1.i = e[i__3].r * x[i__4].i + e[i__3].i * x[
			    i__4].r;
#line 310 "cptrfs.f"
		    ex.r = z__1.r, ex.i = z__1.i;
#line 311 "cptrfs.f"
		    i__3 = i__;
#line 311 "cptrfs.f"
		    z__3.r = bi.r - cx.r, z__3.i = bi.i - cx.i;
#line 311 "cptrfs.f"
		    z__2.r = z__3.r - dx.r, z__2.i = z__3.i - dx.i;
#line 311 "cptrfs.f"
		    z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
#line 311 "cptrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 312 "cptrfs.f"
		    i__3 = i__ - 1;
#line 312 "cptrfs.f"
		    i__4 = i__ - 1 + j * x_dim1;
#line 312 "cptrfs.f"
		    i__5 = i__;
#line 312 "cptrfs.f"
		    i__6 = i__ + 1 + j * x_dim1;
#line 312 "cptrfs.f"
		    rwork[i__] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&
			    bi), abs(d__2)) + ((d__3 = e[i__3].r, abs(d__3)) 
			    + (d__4 = d_imag(&e[i__ - 1]), abs(d__4))) * ((
			    d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[
			    i__ - 1 + j * x_dim1]), abs(d__6))) + ((d__7 = 
			    dx.r, abs(d__7)) + (d__8 = d_imag(&dx), abs(d__8))
			    ) + ((d__9 = e[i__5].r, abs(d__9)) + (d__10 = 
			    d_imag(&e[i__]), abs(d__10))) * ((d__11 = x[i__6]
			    .r, abs(d__11)) + (d__12 = d_imag(&x[i__ + 1 + j *
			     x_dim1]), abs(d__12)));
#line 316 "cptrfs.f"
/* L30: */
#line 316 "cptrfs.f"
		}
#line 317 "cptrfs.f"
		i__2 = *n + j * b_dim1;
#line 317 "cptrfs.f"
		bi.r = b[i__2].r, bi.i = b[i__2].i;
#line 318 "cptrfs.f"
		d_cnjg(&z__2, &e[*n - 1]);
#line 318 "cptrfs.f"
		i__2 = *n - 1 + j * x_dim1;
#line 318 "cptrfs.f"
		z__1.r = z__2.r * x[i__2].r - z__2.i * x[i__2].i, z__1.i = 
			z__2.r * x[i__2].i + z__2.i * x[i__2].r;
#line 318 "cptrfs.f"
		cx.r = z__1.r, cx.i = z__1.i;
#line 319 "cptrfs.f"
		i__2 = *n;
#line 319 "cptrfs.f"
		i__3 = *n + j * x_dim1;
#line 319 "cptrfs.f"
		z__1.r = d__[i__2] * x[i__3].r, z__1.i = d__[i__2] * x[i__3]
			.i;
#line 319 "cptrfs.f"
		dx.r = z__1.r, dx.i = z__1.i;
#line 320 "cptrfs.f"
		i__2 = *n;
#line 320 "cptrfs.f"
		z__2.r = bi.r - cx.r, z__2.i = bi.i - cx.i;
#line 320 "cptrfs.f"
		z__1.r = z__2.r - dx.r, z__1.i = z__2.i - dx.i;
#line 320 "cptrfs.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 321 "cptrfs.f"
		i__2 = *n - 1;
#line 321 "cptrfs.f"
		i__3 = *n - 1 + j * x_dim1;
#line 321 "cptrfs.f"
		rwork[*n] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = e[i__2].r, abs(d__3)) + (d__4 = 
			d_imag(&e[*n - 1]), abs(d__4))) * ((d__5 = x[i__3].r, 
			abs(d__5)) + (d__6 = d_imag(&x[*n - 1 + j * x_dim1]), 
			abs(d__6))) + ((d__7 = dx.r, abs(d__7)) + (d__8 = 
			d_imag(&dx), abs(d__8)));
#line 323 "cptrfs.f"
	    }
#line 324 "cptrfs.f"
	} else {
#line 325 "cptrfs.f"
	    if (*n == 1) {
#line 326 "cptrfs.f"
		i__2 = j * b_dim1 + 1;
#line 326 "cptrfs.f"
		bi.r = b[i__2].r, bi.i = b[i__2].i;
#line 327 "cptrfs.f"
		i__2 = j * x_dim1 + 1;
#line 327 "cptrfs.f"
		z__1.r = d__[1] * x[i__2].r, z__1.i = d__[1] * x[i__2].i;
#line 327 "cptrfs.f"
		dx.r = z__1.r, dx.i = z__1.i;
#line 328 "cptrfs.f"
		z__1.r = bi.r - dx.r, z__1.i = bi.i - dx.i;
#line 328 "cptrfs.f"
		work[1].r = z__1.r, work[1].i = z__1.i;
#line 329 "cptrfs.f"
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4)));
#line 330 "cptrfs.f"
	    } else {
#line 331 "cptrfs.f"
		i__2 = j * b_dim1 + 1;
#line 331 "cptrfs.f"
		bi.r = b[i__2].r, bi.i = b[i__2].i;
#line 332 "cptrfs.f"
		i__2 = j * x_dim1 + 1;
#line 332 "cptrfs.f"
		z__1.r = d__[1] * x[i__2].r, z__1.i = d__[1] * x[i__2].i;
#line 332 "cptrfs.f"
		dx.r = z__1.r, dx.i = z__1.i;
#line 333 "cptrfs.f"
		d_cnjg(&z__2, &e[1]);
#line 333 "cptrfs.f"
		i__2 = j * x_dim1 + 2;
#line 333 "cptrfs.f"
		z__1.r = z__2.r * x[i__2].r - z__2.i * x[i__2].i, z__1.i = 
			z__2.r * x[i__2].i + z__2.i * x[i__2].r;
#line 333 "cptrfs.f"
		ex.r = z__1.r, ex.i = z__1.i;
#line 334 "cptrfs.f"
		z__2.r = bi.r - dx.r, z__2.i = bi.i - dx.i;
#line 334 "cptrfs.f"
		z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
#line 334 "cptrfs.f"
		work[1].r = z__1.r, work[1].i = z__1.i;
#line 335 "cptrfs.f"
		i__2 = j * x_dim1 + 2;
#line 335 "cptrfs.f"
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4))) + ((d__5 = e[1].r, abs(d__5))
			 + (d__6 = d_imag(&e[1]), abs(d__6))) * ((d__7 = x[
			i__2].r, abs(d__7)) + (d__8 = d_imag(&x[j * x_dim1 + 
			2]), abs(d__8)));
#line 337 "cptrfs.f"
		i__2 = *n - 1;
#line 337 "cptrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 338 "cptrfs.f"
		    i__3 = i__ + j * b_dim1;
#line 338 "cptrfs.f"
		    bi.r = b[i__3].r, bi.i = b[i__3].i;
#line 339 "cptrfs.f"
		    i__3 = i__ - 1;
#line 339 "cptrfs.f"
		    i__4 = i__ - 1 + j * x_dim1;
#line 339 "cptrfs.f"
		    z__1.r = e[i__3].r * x[i__4].r - e[i__3].i * x[i__4].i, 
			    z__1.i = e[i__3].r * x[i__4].i + e[i__3].i * x[
			    i__4].r;
#line 339 "cptrfs.f"
		    cx.r = z__1.r, cx.i = z__1.i;
#line 340 "cptrfs.f"
		    i__3 = i__;
#line 340 "cptrfs.f"
		    i__4 = i__ + j * x_dim1;
#line 340 "cptrfs.f"
		    z__1.r = d__[i__3] * x[i__4].r, z__1.i = d__[i__3] * x[
			    i__4].i;
#line 340 "cptrfs.f"
		    dx.r = z__1.r, dx.i = z__1.i;
#line 341 "cptrfs.f"
		    d_cnjg(&z__2, &e[i__]);
#line 341 "cptrfs.f"
		    i__3 = i__ + 1 + j * x_dim1;
#line 341 "cptrfs.f"
		    z__1.r = z__2.r * x[i__3].r - z__2.i * x[i__3].i, z__1.i =
			     z__2.r * x[i__3].i + z__2.i * x[i__3].r;
#line 341 "cptrfs.f"
		    ex.r = z__1.r, ex.i = z__1.i;
#line 342 "cptrfs.f"
		    i__3 = i__;
#line 342 "cptrfs.f"
		    z__3.r = bi.r - cx.r, z__3.i = bi.i - cx.i;
#line 342 "cptrfs.f"
		    z__2.r = z__3.r - dx.r, z__2.i = z__3.i - dx.i;
#line 342 "cptrfs.f"
		    z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
#line 342 "cptrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 343 "cptrfs.f"
		    i__3 = i__ - 1;
#line 343 "cptrfs.f"
		    i__4 = i__ - 1 + j * x_dim1;
#line 343 "cptrfs.f"
		    i__5 = i__;
#line 343 "cptrfs.f"
		    i__6 = i__ + 1 + j * x_dim1;
#line 343 "cptrfs.f"
		    rwork[i__] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&
			    bi), abs(d__2)) + ((d__3 = e[i__3].r, abs(d__3)) 
			    + (d__4 = d_imag(&e[i__ - 1]), abs(d__4))) * ((
			    d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[
			    i__ - 1 + j * x_dim1]), abs(d__6))) + ((d__7 = 
			    dx.r, abs(d__7)) + (d__8 = d_imag(&dx), abs(d__8))
			    ) + ((d__9 = e[i__5].r, abs(d__9)) + (d__10 = 
			    d_imag(&e[i__]), abs(d__10))) * ((d__11 = x[i__6]
			    .r, abs(d__11)) + (d__12 = d_imag(&x[i__ + 1 + j *
			     x_dim1]), abs(d__12)));
#line 347 "cptrfs.f"
/* L40: */
#line 347 "cptrfs.f"
		}
#line 348 "cptrfs.f"
		i__2 = *n + j * b_dim1;
#line 348 "cptrfs.f"
		bi.r = b[i__2].r, bi.i = b[i__2].i;
#line 349 "cptrfs.f"
		i__2 = *n - 1;
#line 349 "cptrfs.f"
		i__3 = *n - 1 + j * x_dim1;
#line 349 "cptrfs.f"
		z__1.r = e[i__2].r * x[i__3].r - e[i__2].i * x[i__3].i, 
			z__1.i = e[i__2].r * x[i__3].i + e[i__2].i * x[i__3]
			.r;
#line 349 "cptrfs.f"
		cx.r = z__1.r, cx.i = z__1.i;
#line 350 "cptrfs.f"
		i__2 = *n;
#line 350 "cptrfs.f"
		i__3 = *n + j * x_dim1;
#line 350 "cptrfs.f"
		z__1.r = d__[i__2] * x[i__3].r, z__1.i = d__[i__2] * x[i__3]
			.i;
#line 350 "cptrfs.f"
		dx.r = z__1.r, dx.i = z__1.i;
#line 351 "cptrfs.f"
		i__2 = *n;
#line 351 "cptrfs.f"
		z__2.r = bi.r - cx.r, z__2.i = bi.i - cx.i;
#line 351 "cptrfs.f"
		z__1.r = z__2.r - dx.r, z__1.i = z__2.i - dx.i;
#line 351 "cptrfs.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 352 "cptrfs.f"
		i__2 = *n - 1;
#line 352 "cptrfs.f"
		i__3 = *n - 1 + j * x_dim1;
#line 352 "cptrfs.f"
		rwork[*n] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = e[i__2].r, abs(d__3)) + (d__4 = 
			d_imag(&e[*n - 1]), abs(d__4))) * ((d__5 = x[i__3].r, 
			abs(d__5)) + (d__6 = d_imag(&x[*n - 1 + j * x_dim1]), 
			abs(d__6))) + ((d__7 = dx.r, abs(d__7)) + (d__8 = 
			d_imag(&dx), abs(d__8)));
#line 354 "cptrfs.f"
	    }
#line 355 "cptrfs.f"
	}

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 366 "cptrfs.f"
	s = 0.;
#line 367 "cptrfs.f"
	i__2 = *n;
#line 367 "cptrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 368 "cptrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 369 "cptrfs.f"
		i__3 = i__;
#line 369 "cptrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 369 "cptrfs.f"
		s = max(d__3,d__4);
#line 370 "cptrfs.f"
	    } else {
/* Computing MAX */
#line 371 "cptrfs.f"
		i__3 = i__;
#line 371 "cptrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 371 "cptrfs.f"
		s = max(d__3,d__4);
#line 373 "cptrfs.f"
	    }
#line 374 "cptrfs.f"
/* L50: */
#line 374 "cptrfs.f"
	}
#line 375 "cptrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 383 "cptrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 388 "cptrfs.f"
	    cpttrs_(uplo, n, &c__1, &df[1], &ef[1], &work[1], n, info, (
		    ftnlen)1);
#line 389 "cptrfs.f"
	    caxpy_(n, &c_b16, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 390 "cptrfs.f"
	    lstres = berr[j];
#line 391 "cptrfs.f"
	    ++count;
#line 392 "cptrfs.f"
	    goto L20;
#line 393 "cptrfs.f"
	}

/*        Bound error from formula */

/*        norm(X - XTRUE) / norm(X) .le. FERR = */
/*        norm( abs(inv(A))* */
/*           ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X) */

/*        where */
/*          norm(Z) is the magnitude of the largest component of Z */
/*          inv(A) is the inverse of A */
/*          abs(Z) is the componentwise absolute value of the matrix or */
/*             vector Z */
/*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
/*          EPS is machine epsilon */

/*        The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B)) */
/*        is incremented by SAFE1 if the i-th component of */
/*        abs(A)*abs(X) + abs(B) is less than SAFE2. */

#line 413 "cptrfs.f"
	i__2 = *n;
#line 413 "cptrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 414 "cptrfs.f"
	    if (rwork[i__] > safe2) {
#line 415 "cptrfs.f"
		i__3 = i__;
#line 415 "cptrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 416 "cptrfs.f"
	    } else {
#line 417 "cptrfs.f"
		i__3 = i__;
#line 417 "cptrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 419 "cptrfs.f"
	    }
#line 420 "cptrfs.f"
/* L60: */
#line 420 "cptrfs.f"
	}
#line 421 "cptrfs.f"
	ix = isamax_(n, &rwork[1], &c__1);
#line 422 "cptrfs.f"
	ferr[j] = rwork[ix];

/*        Estimate the norm of inv(A). */

/*        Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */

/*           m(i,j) =  abs(A(i,j)), i = j, */
/*           m(i,j) = -abs(A(i,j)), i .ne. j, */

/*        and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H. */

/*        Solve M(L) * x = e. */

#line 435 "cptrfs.f"
	rwork[1] = 1.;
#line 436 "cptrfs.f"
	i__2 = *n;
#line 436 "cptrfs.f"
	for (i__ = 2; i__ <= i__2; ++i__) {
#line 437 "cptrfs.f"
	    rwork[i__] = rwork[i__ - 1] * z_abs(&ef[i__ - 1]) + 1.;
#line 438 "cptrfs.f"
/* L70: */
#line 438 "cptrfs.f"
	}

/*        Solve D * M(L)**H * x = b. */

#line 442 "cptrfs.f"
	rwork[*n] /= df[*n];
#line 443 "cptrfs.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {
#line 444 "cptrfs.f"
	    rwork[i__] = rwork[i__] / df[i__] + rwork[i__ + 1] * z_abs(&ef[
		    i__]);
#line 446 "cptrfs.f"
/* L80: */
#line 446 "cptrfs.f"
	}

/*        Compute norm(inv(A)) = max(x(i)), 1<=i<=n. */

#line 450 "cptrfs.f"
	ix = isamax_(n, &rwork[1], &c__1);
#line 451 "cptrfs.f"
	ferr[j] *= (d__1 = rwork[ix], abs(d__1));

/*        Normalize error. */

#line 455 "cptrfs.f"
	lstres = 0.;
#line 456 "cptrfs.f"
	i__2 = *n;
#line 456 "cptrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 457 "cptrfs.f"
	    d__1 = lstres, d__2 = z_abs(&x[i__ + j * x_dim1]);
#line 457 "cptrfs.f"
	    lstres = max(d__1,d__2);
#line 458 "cptrfs.f"
/* L90: */
#line 458 "cptrfs.f"
	}
#line 459 "cptrfs.f"
	if (lstres != 0.) {
#line 459 "cptrfs.f"
	    ferr[j] /= lstres;
#line 459 "cptrfs.f"
	}

#line 462 "cptrfs.f"
/* L100: */
#line 462 "cptrfs.f"
    }

#line 464 "cptrfs.f"
    return 0;

/*     End of CPTRFS */

} /* cptrfs_ */

