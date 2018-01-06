#line 1 "zlatbs.f"
/* zlatbs.f -- translated by f2c (version 20100827).
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

#line 1 "zlatbs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b ZLATBS solves a triangular banded system of equations. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLATBS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatbs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatbs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatbs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, */
/*                          SCALE, CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   CNORM( * ) */
/*       COMPLEX*16         AB( LDAB, * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLATBS solves one of the triangular systems */
/* > */
/* >    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular band matrix.  Here A**T denotes the transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine ZTBSV is called.  If the matrix A */
/* > is singular (A(j,j) = 0 for some j), then s is set to 0 and a */
/* > non-trivial solution to A*x = 0 is returned. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower triangular. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the operation applied to A. */
/* >          = 'N':  Solve A * x = s*b     (No transpose) */
/* >          = 'T':  Solve A**T * x = s*b  (Transpose) */
/* >          = 'C':  Solve A**H * x = s*b  (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A is unit triangular. */
/* >          = 'N':  Non-unit triangular */
/* >          = 'U':  Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] NORMIN */
/* > \verbatim */
/* >          NORMIN is CHARACTER*1 */
/* >          Specifies whether CNORM has been set or not. */
/* >          = 'Y':  CNORM contains the column norms on entry */
/* >          = 'N':  CNORM is not set on entry.  On exit, the norms will */
/* >                  be computed and stored in CNORM. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of subdiagonals or superdiagonals in the */
/* >          triangular matrix A.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          The upper or lower triangular band matrix A, stored in the */
/* >          first KD+1 rows of the array. The j-th column of A is stored */
/* >          in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b. */
/* >          If SCALE = 0, the matrix A is singular or badly scaled, and */
/* >          the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* >          CNORM is DOUBLE PRECISION array, dimension (N) */
/* > */
/* >          If NORMIN = 'Y', CNORM is an input argument and CNORM(j) */
/* >          contains the norm of the off-diagonal part of the j-th column */
/* >          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal */
/* >          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j) */
/* >          must be greater than or equal to the 1-norm. */
/* > */
/* >          If NORMIN = 'N', CNORM is an output argument and CNORM(j) */
/* >          returns the 1-norm of the offdiagonal part of the j-th column */
/* >          of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, ZTBSV */
/* >  is called, otherwise, specific code is used which checks for possible */
/* >  overflow or divide-by-zero at every operation. */
/* > */
/* >  A columnwise scheme is used for solving A*x = b.  The basic algorithm */
/* >  if A is lower triangular is */
/* > */
/* >       x[1:n] := b[1:n] */
/* >       for j = 1, ..., n */
/* >            x(j) := x(j) / A(j,j) */
/* >            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j] */
/* >       end */
/* > */
/* >  Define bounds on the components of x after j iterations of the loop: */
/* >     M(j) = bound on x[1:j] */
/* >     G(j) = bound on x[j+1:n] */
/* >  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}. */
/* > */
/* >  Then for iteration j+1 we have */
/* >     M(j+1) <= G(j) / | A(j+1,j+1) | */
/* >     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] | */
/* >            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | ) */
/* > */
/* >  where CNORM(j+1) is greater than or equal to the infinity-norm of */
/* >  column j+1 of A, not counting the diagonal.  Hence */
/* > */
/* >     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | ) */
/* >                  1<=i<=j */
/* >  and */
/* > */
/* >     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) */
/* >                                   1<=i< j */
/* > */
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTBSV if the */
/* >  reciprocal of the largest M(j), j=1,..,n, is larger than */
/* >  max(underflow, 1/overflow). */
/* > */
/* >  The bound on x(j) is also used to determine when a step in the */
/* >  columnwise method can be performed without fear of overflow.  If */
/* >  the computed bound is greater than a large constant, x is scaled to */
/* >  prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
/* >  1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */
/* > */
/* >  Similarly, a row-wise scheme is used to solve A**T *x = b  or */
/* >  A**H *x = b.  The basic algorithm for A upper triangular is */
/* > */
/* >       for j = 1, ..., n */
/* >            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j) */
/* >       end */
/* > */
/* >  We simultaneously compute two bounds */
/* >       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j */
/* >       M(j) = bound on x(i), 1<=i<=j */
/* > */
/* >  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we */
/* >  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. */
/* >  Then the bound on x(j) is */
/* > */
/* >       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) | */
/* > */
/* >            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| ) */
/* >                      1<=i<=j */
/* > */
/* >  and we can safely call ZTBSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab, 
	doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info,
	 ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen 
	normin_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax;
    static doublecomplex tjjs;
    static doublereal xmax, grow;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer maind;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal tscal;
    static doublecomplex uscal;
    static integer jlast;
    static doublecomplex csumj;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Double Complex */ VOID zdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int ztbsv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal bignum;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static logical notran;
    static integer jfirst;
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);
    static doublereal smlnum;
    static logical nounit;


/*  -- LAPACK auxiliary routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 299 "zlatbs.f"
    /* Parameter adjustments */
#line 299 "zlatbs.f"
    ab_dim1 = *ldab;
#line 299 "zlatbs.f"
    ab_offset = 1 + ab_dim1;
#line 299 "zlatbs.f"
    ab -= ab_offset;
#line 299 "zlatbs.f"
    --x;
#line 299 "zlatbs.f"
    --cnorm;
#line 299 "zlatbs.f"

#line 299 "zlatbs.f"
    /* Function Body */
#line 299 "zlatbs.f"
    *info = 0;
#line 300 "zlatbs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 301 "zlatbs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 302 "zlatbs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 306 "zlatbs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 307 "zlatbs.f"
	*info = -1;
#line 308 "zlatbs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 310 "zlatbs.f"
	*info = -2;
#line 311 "zlatbs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 312 "zlatbs.f"
	*info = -3;
#line 313 "zlatbs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 315 "zlatbs.f"
	*info = -4;
#line 316 "zlatbs.f"
    } else if (*n < 0) {
#line 317 "zlatbs.f"
	*info = -5;
#line 318 "zlatbs.f"
    } else if (*kd < 0) {
#line 319 "zlatbs.f"
	*info = -6;
#line 320 "zlatbs.f"
    } else if (*ldab < *kd + 1) {
#line 321 "zlatbs.f"
	*info = -8;
#line 322 "zlatbs.f"
    }
#line 323 "zlatbs.f"
    if (*info != 0) {
#line 324 "zlatbs.f"
	i__1 = -(*info);
#line 324 "zlatbs.f"
	xerbla_("ZLATBS", &i__1, (ftnlen)6);
#line 325 "zlatbs.f"
	return 0;
#line 326 "zlatbs.f"
    }

/*     Quick return if possible */

#line 330 "zlatbs.f"
    if (*n == 0) {
#line 330 "zlatbs.f"
	return 0;
#line 330 "zlatbs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 335 "zlatbs.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 336 "zlatbs.f"
    bignum = 1. / smlnum;
#line 337 "zlatbs.f"
    dlabad_(&smlnum, &bignum);
#line 338 "zlatbs.f"
    smlnum /= dlamch_("Precision", (ftnlen)9);
#line 339 "zlatbs.f"
    bignum = 1. / smlnum;
#line 340 "zlatbs.f"
    *scale = 1.;

#line 342 "zlatbs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 346 "zlatbs.f"
	if (upper) {

/*           A is upper triangular. */

#line 350 "zlatbs.f"
	    i__1 = *n;
#line 350 "zlatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 351 "zlatbs.f"
		i__2 = *kd, i__3 = j - 1;
#line 351 "zlatbs.f"
		jlen = min(i__2,i__3);
#line 352 "zlatbs.f"
		cnorm[j] = dzasum_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], &
			c__1);
#line 353 "zlatbs.f"
/* L10: */
#line 353 "zlatbs.f"
	    }
#line 354 "zlatbs.f"
	} else {

/*           A is lower triangular. */

#line 358 "zlatbs.f"
	    i__1 = *n;
#line 358 "zlatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 359 "zlatbs.f"
		i__2 = *kd, i__3 = *n - j;
#line 359 "zlatbs.f"
		jlen = min(i__2,i__3);
#line 360 "zlatbs.f"
		if (jlen > 0) {
#line 361 "zlatbs.f"
		    cnorm[j] = dzasum_(&jlen, &ab[j * ab_dim1 + 2], &c__1);
#line 362 "zlatbs.f"
		} else {
#line 363 "zlatbs.f"
		    cnorm[j] = 0.;
#line 364 "zlatbs.f"
		}
#line 365 "zlatbs.f"
/* L20: */
#line 365 "zlatbs.f"
	    }
#line 366 "zlatbs.f"
	}
#line 367 "zlatbs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM/2. */

#line 372 "zlatbs.f"
    imax = idamax_(n, &cnorm[1], &c__1);
#line 373 "zlatbs.f"
    tmax = cnorm[imax];
#line 374 "zlatbs.f"
    if (tmax <= bignum * .5) {
#line 375 "zlatbs.f"
	tscal = 1.;
#line 376 "zlatbs.f"
    } else {
#line 377 "zlatbs.f"
	tscal = .5 / (smlnum * tmax);
#line 378 "zlatbs.f"
	dscal_(n, &tscal, &cnorm[1], &c__1);
#line 379 "zlatbs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine ZTBSV can be used. */

#line 384 "zlatbs.f"
    xmax = 0.;
#line 385 "zlatbs.f"
    i__1 = *n;
#line 385 "zlatbs.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 386 "zlatbs.f"
	i__2 = j;
#line 386 "zlatbs.f"
	d__3 = xmax, d__4 = (d__1 = x[i__2].r / 2., abs(d__1)) + (d__2 = 
		d_imag(&x[j]) / 2., abs(d__2));
#line 386 "zlatbs.f"
	xmax = max(d__3,d__4);
#line 387 "zlatbs.f"
/* L30: */
#line 387 "zlatbs.f"
    }
#line 388 "zlatbs.f"
    xbnd = xmax;
#line 389 "zlatbs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 393 "zlatbs.f"
	if (upper) {
#line 394 "zlatbs.f"
	    jfirst = *n;
#line 395 "zlatbs.f"
	    jlast = 1;
#line 396 "zlatbs.f"
	    jinc = -1;
#line 397 "zlatbs.f"
	    maind = *kd + 1;
#line 398 "zlatbs.f"
	} else {
#line 399 "zlatbs.f"
	    jfirst = 1;
#line 400 "zlatbs.f"
	    jlast = *n;
#line 401 "zlatbs.f"
	    jinc = 1;
#line 402 "zlatbs.f"
	    maind = 1;
#line 403 "zlatbs.f"
	}

#line 405 "zlatbs.f"
	if (tscal != 1.) {
#line 406 "zlatbs.f"
	    grow = 0.;
#line 407 "zlatbs.f"
	    goto L60;
#line 408 "zlatbs.f"
	}

#line 410 "zlatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 417 "zlatbs.f"
	    grow = .5 / max(xbnd,smlnum);
#line 418 "zlatbs.f"
	    xbnd = grow;
#line 419 "zlatbs.f"
	    i__1 = jlast;
#line 419 "zlatbs.f"
	    i__2 = jinc;
#line 419 "zlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 423 "zlatbs.f"
		if (grow <= smlnum) {
#line 423 "zlatbs.f"
		    goto L60;
#line 423 "zlatbs.f"
		}

#line 426 "zlatbs.f"
		i__3 = maind + j * ab_dim1;
#line 426 "zlatbs.f"
		tjjs.r = ab[i__3].r, tjjs.i = ab[i__3].i;
#line 427 "zlatbs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 429 "zlatbs.f"
		if (tjj >= smlnum) {

/*                 M(j) = G(j-1) / abs(A(j,j)) */

/* Computing MIN */
#line 433 "zlatbs.f"
		    d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 433 "zlatbs.f"
		    xbnd = min(d__1,d__2);
#line 434 "zlatbs.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 438 "zlatbs.f"
		    xbnd = 0.;
#line 439 "zlatbs.f"
		}

#line 441 "zlatbs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 445 "zlatbs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 446 "zlatbs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 450 "zlatbs.f"
		    grow = 0.;
#line 451 "zlatbs.f"
		}
#line 452 "zlatbs.f"
/* L40: */
#line 452 "zlatbs.f"
	    }
#line 453 "zlatbs.f"
	    grow = xbnd;
#line 454 "zlatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 460 "zlatbs.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 460 "zlatbs.f"
	    grow = min(d__1,d__2);
#line 461 "zlatbs.f"
	    i__2 = jlast;
#line 461 "zlatbs.f"
	    i__1 = jinc;
#line 461 "zlatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 465 "zlatbs.f"
		if (grow <= smlnum) {
#line 465 "zlatbs.f"
		    goto L60;
#line 465 "zlatbs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 470 "zlatbs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 471 "zlatbs.f"
/* L50: */
#line 471 "zlatbs.f"
	    }
#line 472 "zlatbs.f"
	}
#line 473 "zlatbs.f"
L60:

#line 475 "zlatbs.f"
	;
#line 475 "zlatbs.f"
    } else {

/*        Compute the growth in A**T * x = b  or  A**H * x = b. */

#line 479 "zlatbs.f"
	if (upper) {
#line 480 "zlatbs.f"
	    jfirst = 1;
#line 481 "zlatbs.f"
	    jlast = *n;
#line 482 "zlatbs.f"
	    jinc = 1;
#line 483 "zlatbs.f"
	    maind = *kd + 1;
#line 484 "zlatbs.f"
	} else {
#line 485 "zlatbs.f"
	    jfirst = *n;
#line 486 "zlatbs.f"
	    jlast = 1;
#line 487 "zlatbs.f"
	    jinc = -1;
#line 488 "zlatbs.f"
	    maind = 1;
#line 489 "zlatbs.f"
	}

#line 491 "zlatbs.f"
	if (tscal != 1.) {
#line 492 "zlatbs.f"
	    grow = 0.;
#line 493 "zlatbs.f"
	    goto L90;
#line 494 "zlatbs.f"
	}

#line 496 "zlatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 503 "zlatbs.f"
	    grow = .5 / max(xbnd,smlnum);
#line 504 "zlatbs.f"
	    xbnd = grow;
#line 505 "zlatbs.f"
	    i__1 = jlast;
#line 505 "zlatbs.f"
	    i__2 = jinc;
#line 505 "zlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 509 "zlatbs.f"
		if (grow <= smlnum) {
#line 509 "zlatbs.f"
		    goto L90;
#line 509 "zlatbs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 514 "zlatbs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 515 "zlatbs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 515 "zlatbs.f"
		grow = min(d__1,d__2);

#line 517 "zlatbs.f"
		i__3 = maind + j * ab_dim1;
#line 517 "zlatbs.f"
		tjjs.r = ab[i__3].r, tjjs.i = ab[i__3].i;
#line 518 "zlatbs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 520 "zlatbs.f"
		if (tjj >= smlnum) {

/*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 524 "zlatbs.f"
		    if (xj > tjj) {
#line 524 "zlatbs.f"
			xbnd *= tjj / xj;
#line 524 "zlatbs.f"
		    }
#line 526 "zlatbs.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 530 "zlatbs.f"
		    xbnd = 0.;
#line 531 "zlatbs.f"
		}
#line 532 "zlatbs.f"
/* L70: */
#line 532 "zlatbs.f"
	    }
#line 533 "zlatbs.f"
	    grow = min(grow,xbnd);
#line 534 "zlatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 540 "zlatbs.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 540 "zlatbs.f"
	    grow = min(d__1,d__2);
#line 541 "zlatbs.f"
	    i__2 = jlast;
#line 541 "zlatbs.f"
	    i__1 = jinc;
#line 541 "zlatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 545 "zlatbs.f"
		if (grow <= smlnum) {
#line 545 "zlatbs.f"
		    goto L90;
#line 545 "zlatbs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 550 "zlatbs.f"
		xj = cnorm[j] + 1.;
#line 551 "zlatbs.f"
		grow /= xj;
#line 552 "zlatbs.f"
/* L80: */
#line 552 "zlatbs.f"
	    }
#line 553 "zlatbs.f"
	}
#line 554 "zlatbs.f"
L90:
#line 555 "zlatbs.f"
	;
#line 555 "zlatbs.f"
    }

#line 557 "zlatbs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 562 "zlatbs.f"
	ztbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &x[1], &c__1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 563 "zlatbs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 567 "zlatbs.f"
	if (xmax > bignum * .5) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 572 "zlatbs.f"
	    *scale = bignum * .5 / xmax;
#line 573 "zlatbs.f"
	    zdscal_(n, scale, &x[1], &c__1);
#line 574 "zlatbs.f"
	    xmax = bignum;
#line 575 "zlatbs.f"
	} else {
#line 576 "zlatbs.f"
	    xmax *= 2.;
#line 577 "zlatbs.f"
	}

#line 579 "zlatbs.f"
	if (notran) {

/*           Solve A * x = b */

#line 583 "zlatbs.f"
	    i__1 = jlast;
#line 583 "zlatbs.f"
	    i__2 = jinc;
#line 583 "zlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 587 "zlatbs.f"
		i__3 = j;
#line 587 "zlatbs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 588 "zlatbs.f"
		if (nounit) {
#line 589 "zlatbs.f"
		    i__3 = maind + j * ab_dim1;
#line 589 "zlatbs.f"
		    z__1.r = tscal * ab[i__3].r, z__1.i = tscal * ab[i__3].i;
#line 589 "zlatbs.f"
		    tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 590 "zlatbs.f"
		} else {
#line 591 "zlatbs.f"
		    tjjs.r = tscal, tjjs.i = 0.;
#line 592 "zlatbs.f"
		    if (tscal == 1.) {
#line 592 "zlatbs.f"
			goto L110;
#line 592 "zlatbs.f"
		    }
#line 594 "zlatbs.f"
		}
#line 595 "zlatbs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));
#line 596 "zlatbs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 600 "zlatbs.f"
		    if (tjj < 1.) {
#line 601 "zlatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 605 "zlatbs.f"
			    rec = 1. / xj;
#line 606 "zlatbs.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 607 "zlatbs.f"
			    *scale *= rec;
#line 608 "zlatbs.f"
			    xmax *= rec;
#line 609 "zlatbs.f"
			}
#line 610 "zlatbs.f"
		    }
#line 611 "zlatbs.f"
		    i__3 = j;
#line 611 "zlatbs.f"
		    zladiv_(&z__1, &x[j], &tjjs);
#line 611 "zlatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 612 "zlatbs.f"
		    i__3 = j;
#line 612 "zlatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 613 "zlatbs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 617 "zlatbs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 622 "zlatbs.f"
			rec = tjj * bignum / xj;
#line 623 "zlatbs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 628 "zlatbs.f"
			    rec /= cnorm[j];
#line 629 "zlatbs.f"
			}
#line 630 "zlatbs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 631 "zlatbs.f"
			*scale *= rec;
#line 632 "zlatbs.f"
			xmax *= rec;
#line 633 "zlatbs.f"
		    }
#line 634 "zlatbs.f"
		    i__3 = j;
#line 634 "zlatbs.f"
		    zladiv_(&z__1, &x[j], &tjjs);
#line 634 "zlatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 635 "zlatbs.f"
		    i__3 = j;
#line 635 "zlatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 636 "zlatbs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 641 "zlatbs.f"
		    i__3 = *n;
#line 641 "zlatbs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 642 "zlatbs.f"
			i__4 = i__;
#line 642 "zlatbs.f"
			x[i__4].r = 0., x[i__4].i = 0.;
#line 643 "zlatbs.f"
/* L100: */
#line 643 "zlatbs.f"
		    }
#line 644 "zlatbs.f"
		    i__3 = j;
#line 644 "zlatbs.f"
		    x[i__3].r = 1., x[i__3].i = 0.;
#line 645 "zlatbs.f"
		    xj = 1.;
#line 646 "zlatbs.f"
		    *scale = 0.;
#line 647 "zlatbs.f"
		    xmax = 0.;
#line 648 "zlatbs.f"
		}
#line 649 "zlatbs.f"
L110:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 654 "zlatbs.f"
		if (xj > 1.) {
#line 655 "zlatbs.f"
		    rec = 1. / xj;
#line 656 "zlatbs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 660 "zlatbs.f"
			rec *= .5;
#line 661 "zlatbs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 662 "zlatbs.f"
			*scale *= rec;
#line 663 "zlatbs.f"
		    }
#line 664 "zlatbs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 668 "zlatbs.f"
		    zdscal_(n, &c_b36, &x[1], &c__1);
#line 669 "zlatbs.f"
		    *scale *= .5;
#line 670 "zlatbs.f"
		}

#line 672 "zlatbs.f"
		if (upper) {
#line 673 "zlatbs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) - */
/*                                             x(j)* A(max(1,j-kd):j-1,j) */

/* Computing MIN */
#line 679 "zlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 679 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 680 "zlatbs.f"
			i__3 = j;
#line 680 "zlatbs.f"
			z__2.r = -x[i__3].r, z__2.i = -x[i__3].i;
#line 680 "zlatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 680 "zlatbs.f"
			zaxpy_(&jlen, &z__1, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 682 "zlatbs.f"
			i__3 = j - 1;
#line 682 "zlatbs.f"
			i__ = izamax_(&i__3, &x[1], &c__1);
#line 683 "zlatbs.f"
			i__3 = i__;
#line 683 "zlatbs.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 684 "zlatbs.f"
		    }
#line 685 "zlatbs.f"
		} else if (j < *n) {

/*                 Compute the update */
/*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) - */
/*                                          x(j) * A(j+1:min(j+kd,n),j) */

/* Computing MIN */
#line 691 "zlatbs.f"
		    i__3 = *kd, i__4 = *n - j;
#line 691 "zlatbs.f"
		    jlen = min(i__3,i__4);
#line 692 "zlatbs.f"
		    if (jlen > 0) {
#line 692 "zlatbs.f"
			i__3 = j;
#line 692 "zlatbs.f"
			z__2.r = -x[i__3].r, z__2.i = -x[i__3].i;
#line 692 "zlatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 692 "zlatbs.f"
			zaxpy_(&jlen, &z__1, &ab[j * ab_dim1 + 2], &c__1, &x[
				j + 1], &c__1);
#line 692 "zlatbs.f"
		    }
#line 695 "zlatbs.f"
		    i__3 = *n - j;
#line 695 "zlatbs.f"
		    i__ = j + izamax_(&i__3, &x[j + 1], &c__1);
#line 696 "zlatbs.f"
		    i__3 = i__;
#line 696 "zlatbs.f"
		    xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
			    i__]), abs(d__2));
#line 697 "zlatbs.f"
		}
#line 698 "zlatbs.f"
/* L120: */
#line 698 "zlatbs.f"
	    }

#line 700 "zlatbs.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Solve A**T * x = b */

#line 704 "zlatbs.f"
	    i__2 = jlast;
#line 704 "zlatbs.f"
	    i__1 = jinc;
#line 704 "zlatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 709 "zlatbs.f"
		i__3 = j;
#line 709 "zlatbs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 710 "zlatbs.f"
		uscal.r = tscal, uscal.i = 0.;
#line 711 "zlatbs.f"
		rec = 1. / max(xmax,1.);
#line 712 "zlatbs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 716 "zlatbs.f"
		    rec *= .5;
#line 717 "zlatbs.f"
		    if (nounit) {
#line 718 "zlatbs.f"
			i__3 = maind + j * ab_dim1;
#line 718 "zlatbs.f"
			z__1.r = tscal * ab[i__3].r, z__1.i = tscal * ab[i__3]
				.i;
#line 718 "zlatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 719 "zlatbs.f"
		    } else {
#line 720 "zlatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 721 "zlatbs.f"
		    }
#line 722 "zlatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 723 "zlatbs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 727 "zlatbs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 727 "zlatbs.f"
			rec = min(d__1,d__2);
#line 728 "zlatbs.f"
			zladiv_(&z__1, &uscal, &tjjs);
#line 728 "zlatbs.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 729 "zlatbs.f"
		    }
#line 730 "zlatbs.f"
		    if (rec < 1.) {
#line 731 "zlatbs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 732 "zlatbs.f"
			*scale *= rec;
#line 733 "zlatbs.f"
			xmax *= rec;
#line 734 "zlatbs.f"
		    }
#line 735 "zlatbs.f"
		}

#line 737 "zlatbs.f"
		csumj.r = 0., csumj.i = 0.;
#line 738 "zlatbs.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call ZDOTU to perform the dot product. */

#line 743 "zlatbs.f"
		    if (upper) {
/* Computing MIN */
#line 744 "zlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 744 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 745 "zlatbs.f"
			zdotu_(&z__1, &jlen, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 745 "zlatbs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 747 "zlatbs.f"
		    } else {
/* Computing MIN */
#line 748 "zlatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 748 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 749 "zlatbs.f"
			if (jlen > 1) {
#line 749 "zlatbs.f"
			    zdotu_(&z__1, &jlen, &ab[j * ab_dim1 + 2], &c__1, 
				    &x[j + 1], &c__1);
#line 749 "zlatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 749 "zlatbs.f"
			}
#line 752 "zlatbs.f"
		    }
#line 753 "zlatbs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 757 "zlatbs.f"
		    if (upper) {
/* Computing MIN */
#line 758 "zlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 758 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 759 "zlatbs.f"
			i__3 = jlen;
#line 759 "zlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 760 "zlatbs.f"
			    i__4 = *kd + i__ - jlen + j * ab_dim1;
#line 760 "zlatbs.f"
			    z__3.r = ab[i__4].r * uscal.r - ab[i__4].i * 
				    uscal.i, z__3.i = ab[i__4].r * uscal.i + 
				    ab[i__4].i * uscal.r;
#line 760 "zlatbs.f"
			    i__5 = j - jlen - 1 + i__;
#line 760 "zlatbs.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 760 "zlatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 760 "zlatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 762 "zlatbs.f"
/* L130: */
#line 762 "zlatbs.f"
			}
#line 763 "zlatbs.f"
		    } else {
/* Computing MIN */
#line 764 "zlatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 764 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 765 "zlatbs.f"
			i__3 = jlen;
#line 765 "zlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 766 "zlatbs.f"
			    i__4 = i__ + 1 + j * ab_dim1;
#line 766 "zlatbs.f"
			    z__3.r = ab[i__4].r * uscal.r - ab[i__4].i * 
				    uscal.i, z__3.i = ab[i__4].r * uscal.i + 
				    ab[i__4].i * uscal.r;
#line 766 "zlatbs.f"
			    i__5 = j + i__;
#line 766 "zlatbs.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 766 "zlatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 766 "zlatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 767 "zlatbs.f"
/* L140: */
#line 767 "zlatbs.f"
			}
#line 768 "zlatbs.f"
		    }
#line 769 "zlatbs.f"
		}

#line 771 "zlatbs.f"
		z__1.r = tscal, z__1.i = 0.;
#line 771 "zlatbs.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 776 "zlatbs.f"
		    i__3 = j;
#line 776 "zlatbs.f"
		    i__4 = j;
#line 776 "zlatbs.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 776 "zlatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 777 "zlatbs.f"
		    i__3 = j;
#line 777 "zlatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 778 "zlatbs.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 782 "zlatbs.f"
			i__3 = maind + j * ab_dim1;
#line 782 "zlatbs.f"
			z__1.r = tscal * ab[i__3].r, z__1.i = tscal * ab[i__3]
				.i;
#line 782 "zlatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 783 "zlatbs.f"
		    } else {
#line 784 "zlatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 785 "zlatbs.f"
			if (tscal == 1.) {
#line 785 "zlatbs.f"
			    goto L160;
#line 785 "zlatbs.f"
			}
#line 787 "zlatbs.f"
		    }
#line 788 "zlatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 789 "zlatbs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 793 "zlatbs.f"
			if (tjj < 1.) {
#line 794 "zlatbs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 798 "zlatbs.f"
				rec = 1. / xj;
#line 799 "zlatbs.f"
				zdscal_(n, &rec, &x[1], &c__1);
#line 800 "zlatbs.f"
				*scale *= rec;
#line 801 "zlatbs.f"
				xmax *= rec;
#line 802 "zlatbs.f"
			    }
#line 803 "zlatbs.f"
			}
#line 804 "zlatbs.f"
			i__3 = j;
#line 804 "zlatbs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 804 "zlatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 805 "zlatbs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 809 "zlatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 813 "zlatbs.f"
			    rec = tjj * bignum / xj;
#line 814 "zlatbs.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 815 "zlatbs.f"
			    *scale *= rec;
#line 816 "zlatbs.f"
			    xmax *= rec;
#line 817 "zlatbs.f"
			}
#line 818 "zlatbs.f"
			i__3 = j;
#line 818 "zlatbs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 818 "zlatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 819 "zlatbs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**T *x = 0. */

#line 824 "zlatbs.f"
			i__3 = *n;
#line 824 "zlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 825 "zlatbs.f"
			    i__4 = i__;
#line 825 "zlatbs.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 826 "zlatbs.f"
/* L150: */
#line 826 "zlatbs.f"
			}
#line 827 "zlatbs.f"
			i__3 = j;
#line 827 "zlatbs.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 828 "zlatbs.f"
			*scale = 0.;
#line 829 "zlatbs.f"
			xmax = 0.;
#line 830 "zlatbs.f"
		    }
#line 831 "zlatbs.f"
L160:
#line 832 "zlatbs.f"
		    ;
#line 832 "zlatbs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 837 "zlatbs.f"
		    i__3 = j;
#line 837 "zlatbs.f"
		    zladiv_(&z__2, &x[j], &tjjs);
#line 837 "zlatbs.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 837 "zlatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 838 "zlatbs.f"
		}
/* Computing MAX */
#line 839 "zlatbs.f"
		i__3 = j;
#line 839 "zlatbs.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 839 "zlatbs.f"
		xmax = max(d__3,d__4);
#line 840 "zlatbs.f"
/* L170: */
#line 840 "zlatbs.f"
	    }

#line 842 "zlatbs.f"
	} else {

/*           Solve A**H * x = b */

#line 846 "zlatbs.f"
	    i__1 = jlast;
#line 846 "zlatbs.f"
	    i__2 = jinc;
#line 846 "zlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 851 "zlatbs.f"
		i__3 = j;
#line 851 "zlatbs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 852 "zlatbs.f"
		uscal.r = tscal, uscal.i = 0.;
#line 853 "zlatbs.f"
		rec = 1. / max(xmax,1.);
#line 854 "zlatbs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 858 "zlatbs.f"
		    rec *= .5;
#line 859 "zlatbs.f"
		    if (nounit) {
#line 860 "zlatbs.f"
			d_cnjg(&z__2, &ab[maind + j * ab_dim1]);
#line 860 "zlatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 860 "zlatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 861 "zlatbs.f"
		    } else {
#line 862 "zlatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 863 "zlatbs.f"
		    }
#line 864 "zlatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 865 "zlatbs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 869 "zlatbs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 869 "zlatbs.f"
			rec = min(d__1,d__2);
#line 870 "zlatbs.f"
			zladiv_(&z__1, &uscal, &tjjs);
#line 870 "zlatbs.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 871 "zlatbs.f"
		    }
#line 872 "zlatbs.f"
		    if (rec < 1.) {
#line 873 "zlatbs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 874 "zlatbs.f"
			*scale *= rec;
#line 875 "zlatbs.f"
			xmax *= rec;
#line 876 "zlatbs.f"
		    }
#line 877 "zlatbs.f"
		}

#line 879 "zlatbs.f"
		csumj.r = 0., csumj.i = 0.;
#line 880 "zlatbs.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call ZDOTC to perform the dot product. */

#line 885 "zlatbs.f"
		    if (upper) {
/* Computing MIN */
#line 886 "zlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 886 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 887 "zlatbs.f"
			zdotc_(&z__1, &jlen, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 887 "zlatbs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 889 "zlatbs.f"
		    } else {
/* Computing MIN */
#line 890 "zlatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 890 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 891 "zlatbs.f"
			if (jlen > 1) {
#line 891 "zlatbs.f"
			    zdotc_(&z__1, &jlen, &ab[j * ab_dim1 + 2], &c__1, 
				    &x[j + 1], &c__1);
#line 891 "zlatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 891 "zlatbs.f"
			}
#line 894 "zlatbs.f"
		    }
#line 895 "zlatbs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 899 "zlatbs.f"
		    if (upper) {
/* Computing MIN */
#line 900 "zlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 900 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 901 "zlatbs.f"
			i__3 = jlen;
#line 901 "zlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 902 "zlatbs.f"
			    d_cnjg(&z__4, &ab[*kd + i__ - jlen + j * ab_dim1])
				    ;
#line 902 "zlatbs.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 902 "zlatbs.f"
			    i__4 = j - jlen - 1 + i__;
#line 902 "zlatbs.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 902 "zlatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 902 "zlatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 904 "zlatbs.f"
/* L180: */
#line 904 "zlatbs.f"
			}
#line 905 "zlatbs.f"
		    } else {
/* Computing MIN */
#line 906 "zlatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 906 "zlatbs.f"
			jlen = min(i__3,i__4);
#line 907 "zlatbs.f"
			i__3 = jlen;
#line 907 "zlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 908 "zlatbs.f"
			    d_cnjg(&z__4, &ab[i__ + 1 + j * ab_dim1]);
#line 908 "zlatbs.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 908 "zlatbs.f"
			    i__4 = j + i__;
#line 908 "zlatbs.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 908 "zlatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 908 "zlatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 910 "zlatbs.f"
/* L190: */
#line 910 "zlatbs.f"
			}
#line 911 "zlatbs.f"
		    }
#line 912 "zlatbs.f"
		}

#line 914 "zlatbs.f"
		z__1.r = tscal, z__1.i = 0.;
#line 914 "zlatbs.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 919 "zlatbs.f"
		    i__3 = j;
#line 919 "zlatbs.f"
		    i__4 = j;
#line 919 "zlatbs.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 919 "zlatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 920 "zlatbs.f"
		    i__3 = j;
#line 920 "zlatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 921 "zlatbs.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 925 "zlatbs.f"
			d_cnjg(&z__2, &ab[maind + j * ab_dim1]);
#line 925 "zlatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 925 "zlatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 926 "zlatbs.f"
		    } else {
#line 927 "zlatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 928 "zlatbs.f"
			if (tscal == 1.) {
#line 928 "zlatbs.f"
			    goto L210;
#line 928 "zlatbs.f"
			}
#line 930 "zlatbs.f"
		    }
#line 931 "zlatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 932 "zlatbs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 936 "zlatbs.f"
			if (tjj < 1.) {
#line 937 "zlatbs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 941 "zlatbs.f"
				rec = 1. / xj;
#line 942 "zlatbs.f"
				zdscal_(n, &rec, &x[1], &c__1);
#line 943 "zlatbs.f"
				*scale *= rec;
#line 944 "zlatbs.f"
				xmax *= rec;
#line 945 "zlatbs.f"
			    }
#line 946 "zlatbs.f"
			}
#line 947 "zlatbs.f"
			i__3 = j;
#line 947 "zlatbs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 947 "zlatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 948 "zlatbs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 952 "zlatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 956 "zlatbs.f"
			    rec = tjj * bignum / xj;
#line 957 "zlatbs.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 958 "zlatbs.f"
			    *scale *= rec;
#line 959 "zlatbs.f"
			    xmax *= rec;
#line 960 "zlatbs.f"
			}
#line 961 "zlatbs.f"
			i__3 = j;
#line 961 "zlatbs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 961 "zlatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 962 "zlatbs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**H *x = 0. */

#line 967 "zlatbs.f"
			i__3 = *n;
#line 967 "zlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 968 "zlatbs.f"
			    i__4 = i__;
#line 968 "zlatbs.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 969 "zlatbs.f"
/* L200: */
#line 969 "zlatbs.f"
			}
#line 970 "zlatbs.f"
			i__3 = j;
#line 970 "zlatbs.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 971 "zlatbs.f"
			*scale = 0.;
#line 972 "zlatbs.f"
			xmax = 0.;
#line 973 "zlatbs.f"
		    }
#line 974 "zlatbs.f"
L210:
#line 975 "zlatbs.f"
		    ;
#line 975 "zlatbs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 980 "zlatbs.f"
		    i__3 = j;
#line 980 "zlatbs.f"
		    zladiv_(&z__2, &x[j], &tjjs);
#line 980 "zlatbs.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 980 "zlatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 981 "zlatbs.f"
		}
/* Computing MAX */
#line 982 "zlatbs.f"
		i__3 = j;
#line 982 "zlatbs.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 982 "zlatbs.f"
		xmax = max(d__3,d__4);
#line 983 "zlatbs.f"
/* L220: */
#line 983 "zlatbs.f"
	    }
#line 984 "zlatbs.f"
	}
#line 985 "zlatbs.f"
	*scale /= tscal;
#line 986 "zlatbs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 990 "zlatbs.f"
    if (tscal != 1.) {
#line 991 "zlatbs.f"
	d__1 = 1. / tscal;
#line 991 "zlatbs.f"
	dscal_(n, &d__1, &cnorm[1], &c__1);
#line 992 "zlatbs.f"
    }

#line 994 "zlatbs.f"
    return 0;

/*     End of ZLATBS */

} /* zlatbs_ */

