#line 1 "dlatrs.f"
/* dlatrs.f -- translated by f2c (version 20100827).
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

#line 1 "dlatrs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b DLATRS solves a triangular system of equations with the scale factor set to prevent overflow. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLATRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), CNORM( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATRS solves one of the triangular systems */
/* > */
/* >    A *x = s*b  or  A**T *x = s*b */
/* > */
/* > with scaling to prevent overflow.  Here A is an upper or lower */
/* > triangular matrix, A**T denotes the transpose of A, x and b are */
/* > n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A */
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
/* >          = 'N':  Solve A * x = s*b  (No transpose) */
/* >          = 'T':  Solve A**T* x = s*b  (Transpose) */
/* >          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose) */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The triangular matrix A.  If UPLO = 'U', the leading n by n */
/* >          upper triangular part of the array A contains the upper */
/* >          triangular matrix, and the strictly lower triangular part of */
/* >          A is not referenced.  If UPLO = 'L', the leading n by n lower */
/* >          triangular part of the array A contains the lower triangular */
/* >          matrix, and the strictly upper triangular part of A is not */
/* >          referenced.  If DIAG = 'U', the diagonal elements of A are */
/* >          also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max (1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b  or  A**T* x = s*b. */
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

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, DTRSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the */
/* >  reciprocal of the largest M(j), j=1,..,n, is larger than */
/* >  max(underflow, 1/overflow). */
/* > */
/* >  The bound on x(j) is also used to determine when a step in the */
/* >  columnwise method can be performed without fear of overflow.  If */
/* >  the computed bound is greater than a large constant, x is scaled to */
/* >  prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
/* >  1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */
/* > */
/* >  Similarly, a row-wise scheme is used to solve A**T*x = b.  The basic */
/* >  algorithm for A upper triangular is */
/* > */
/* >       for j = 1, ..., n */
/* >            x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j) */
/* >       end */
/* > */
/* >  We simultaneously compute two bounds */
/* >       G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j */
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
/* >  and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *a, integer *lda, doublereal *x, 
	doublereal *scale, doublereal *cnorm, integer *info, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax, tjjs, xmax, grow, sumj;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal tscal, uscal;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer jlast;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notran;
    static integer jfirst;
    static doublereal smlnum;
    static logical nounit;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
/*     .. Executable Statements .. */

#line 281 "dlatrs.f"
    /* Parameter adjustments */
#line 281 "dlatrs.f"
    a_dim1 = *lda;
#line 281 "dlatrs.f"
    a_offset = 1 + a_dim1;
#line 281 "dlatrs.f"
    a -= a_offset;
#line 281 "dlatrs.f"
    --x;
#line 281 "dlatrs.f"
    --cnorm;
#line 281 "dlatrs.f"

#line 281 "dlatrs.f"
    /* Function Body */
#line 281 "dlatrs.f"
    *info = 0;
#line 282 "dlatrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 283 "dlatrs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 284 "dlatrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 288 "dlatrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 289 "dlatrs.f"
	*info = -1;
#line 290 "dlatrs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 292 "dlatrs.f"
	*info = -2;
#line 293 "dlatrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 294 "dlatrs.f"
	*info = -3;
#line 295 "dlatrs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 297 "dlatrs.f"
	*info = -4;
#line 298 "dlatrs.f"
    } else if (*n < 0) {
#line 299 "dlatrs.f"
	*info = -5;
#line 300 "dlatrs.f"
    } else if (*lda < max(1,*n)) {
#line 301 "dlatrs.f"
	*info = -7;
#line 302 "dlatrs.f"
    }
#line 303 "dlatrs.f"
    if (*info != 0) {
#line 304 "dlatrs.f"
	i__1 = -(*info);
#line 304 "dlatrs.f"
	xerbla_("DLATRS", &i__1, (ftnlen)6);
#line 305 "dlatrs.f"
	return 0;
#line 306 "dlatrs.f"
    }

/*     Quick return if possible */

#line 310 "dlatrs.f"
    if (*n == 0) {
#line 310 "dlatrs.f"
	return 0;
#line 310 "dlatrs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 315 "dlatrs.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 316 "dlatrs.f"
    bignum = 1. / smlnum;
#line 317 "dlatrs.f"
    *scale = 1.;

#line 319 "dlatrs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 323 "dlatrs.f"
	if (upper) {

/*           A is upper triangular. */

#line 327 "dlatrs.f"
	    i__1 = *n;
#line 327 "dlatrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 328 "dlatrs.f"
		i__2 = j - 1;
#line 328 "dlatrs.f"
		cnorm[j] = dasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
#line 329 "dlatrs.f"
/* L10: */
#line 329 "dlatrs.f"
	    }
#line 330 "dlatrs.f"
	} else {

/*           A is lower triangular. */

#line 334 "dlatrs.f"
	    i__1 = *n - 1;
#line 334 "dlatrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 335 "dlatrs.f"
		i__2 = *n - j;
#line 335 "dlatrs.f"
		cnorm[j] = dasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
#line 336 "dlatrs.f"
/* L20: */
#line 336 "dlatrs.f"
	    }
#line 337 "dlatrs.f"
	    cnorm[*n] = 0.;
#line 338 "dlatrs.f"
	}
#line 339 "dlatrs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM. */

#line 344 "dlatrs.f"
    imax = idamax_(n, &cnorm[1], &c__1);
#line 345 "dlatrs.f"
    tmax = cnorm[imax];
#line 346 "dlatrs.f"
    if (tmax <= bignum) {
#line 347 "dlatrs.f"
	tscal = 1.;
#line 348 "dlatrs.f"
    } else {
#line 349 "dlatrs.f"
	tscal = 1. / (smlnum * tmax);
#line 350 "dlatrs.f"
	dscal_(n, &tscal, &cnorm[1], &c__1);
#line 351 "dlatrs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine DTRSV can be used. */

#line 356 "dlatrs.f"
    j = idamax_(n, &x[1], &c__1);
#line 357 "dlatrs.f"
    xmax = (d__1 = x[j], abs(d__1));
#line 358 "dlatrs.f"
    xbnd = xmax;
#line 359 "dlatrs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 363 "dlatrs.f"
	if (upper) {
#line 364 "dlatrs.f"
	    jfirst = *n;
#line 365 "dlatrs.f"
	    jlast = 1;
#line 366 "dlatrs.f"
	    jinc = -1;
#line 367 "dlatrs.f"
	} else {
#line 368 "dlatrs.f"
	    jfirst = 1;
#line 369 "dlatrs.f"
	    jlast = *n;
#line 370 "dlatrs.f"
	    jinc = 1;
#line 371 "dlatrs.f"
	}

#line 373 "dlatrs.f"
	if (tscal != 1.) {
#line 374 "dlatrs.f"
	    grow = 0.;
#line 375 "dlatrs.f"
	    goto L50;
#line 376 "dlatrs.f"
	}

#line 378 "dlatrs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 385 "dlatrs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 386 "dlatrs.f"
	    xbnd = grow;
#line 387 "dlatrs.f"
	    i__1 = jlast;
#line 387 "dlatrs.f"
	    i__2 = jinc;
#line 387 "dlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 391 "dlatrs.f"
		if (grow <= smlnum) {
#line 391 "dlatrs.f"
		    goto L50;
#line 391 "dlatrs.f"
		}

/*              M(j) = G(j-1) / abs(A(j,j)) */

#line 396 "dlatrs.f"
		tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
/* Computing MIN */
#line 397 "dlatrs.f"
		d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 397 "dlatrs.f"
		xbnd = min(d__1,d__2);
#line 398 "dlatrs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 402 "dlatrs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 403 "dlatrs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 407 "dlatrs.f"
		    grow = 0.;
#line 408 "dlatrs.f"
		}
#line 409 "dlatrs.f"
/* L30: */
#line 409 "dlatrs.f"
	    }
#line 410 "dlatrs.f"
	    grow = xbnd;
#line 411 "dlatrs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 417 "dlatrs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 417 "dlatrs.f"
	    grow = min(d__1,d__2);
#line 418 "dlatrs.f"
	    i__2 = jlast;
#line 418 "dlatrs.f"
	    i__1 = jinc;
#line 418 "dlatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 422 "dlatrs.f"
		if (grow <= smlnum) {
#line 422 "dlatrs.f"
		    goto L50;
#line 422 "dlatrs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 427 "dlatrs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 428 "dlatrs.f"
/* L40: */
#line 428 "dlatrs.f"
	    }
#line 429 "dlatrs.f"
	}
#line 430 "dlatrs.f"
L50:

#line 432 "dlatrs.f"
	;
#line 432 "dlatrs.f"
    } else {

/*        Compute the growth in A**T * x = b. */

#line 436 "dlatrs.f"
	if (upper) {
#line 437 "dlatrs.f"
	    jfirst = 1;
#line 438 "dlatrs.f"
	    jlast = *n;
#line 439 "dlatrs.f"
	    jinc = 1;
#line 440 "dlatrs.f"
	} else {
#line 441 "dlatrs.f"
	    jfirst = *n;
#line 442 "dlatrs.f"
	    jlast = 1;
#line 443 "dlatrs.f"
	    jinc = -1;
#line 444 "dlatrs.f"
	}

#line 446 "dlatrs.f"
	if (tscal != 1.) {
#line 447 "dlatrs.f"
	    grow = 0.;
#line 448 "dlatrs.f"
	    goto L80;
#line 449 "dlatrs.f"
	}

#line 451 "dlatrs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 458 "dlatrs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 459 "dlatrs.f"
	    xbnd = grow;
#line 460 "dlatrs.f"
	    i__1 = jlast;
#line 460 "dlatrs.f"
	    i__2 = jinc;
#line 460 "dlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 464 "dlatrs.f"
		if (grow <= smlnum) {
#line 464 "dlatrs.f"
		    goto L80;
#line 464 "dlatrs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 469 "dlatrs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 470 "dlatrs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 470 "dlatrs.f"
		grow = min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 474 "dlatrs.f"
		tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 475 "dlatrs.f"
		if (xj > tjj) {
#line 475 "dlatrs.f"
		    xbnd *= tjj / xj;
#line 475 "dlatrs.f"
		}
#line 477 "dlatrs.f"
/* L60: */
#line 477 "dlatrs.f"
	    }
#line 478 "dlatrs.f"
	    grow = min(grow,xbnd);
#line 479 "dlatrs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 485 "dlatrs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 485 "dlatrs.f"
	    grow = min(d__1,d__2);
#line 486 "dlatrs.f"
	    i__2 = jlast;
#line 486 "dlatrs.f"
	    i__1 = jinc;
#line 486 "dlatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 490 "dlatrs.f"
		if (grow <= smlnum) {
#line 490 "dlatrs.f"
		    goto L80;
#line 490 "dlatrs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 495 "dlatrs.f"
		xj = cnorm[j] + 1.;
#line 496 "dlatrs.f"
		grow /= xj;
#line 497 "dlatrs.f"
/* L70: */
#line 497 "dlatrs.f"
	    }
#line 498 "dlatrs.f"
	}
#line 499 "dlatrs.f"
L80:
#line 500 "dlatrs.f"
	;
#line 500 "dlatrs.f"
    }

#line 502 "dlatrs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 507 "dlatrs.f"
	dtrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, (ftnlen)
		1, (ftnlen)1, (ftnlen)1);
#line 508 "dlatrs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 512 "dlatrs.f"
	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 517 "dlatrs.f"
	    *scale = bignum / xmax;
#line 518 "dlatrs.f"
	    dscal_(n, scale, &x[1], &c__1);
#line 519 "dlatrs.f"
	    xmax = bignum;
#line 520 "dlatrs.f"
	}

#line 522 "dlatrs.f"
	if (notran) {

/*           Solve A * x = b */

#line 526 "dlatrs.f"
	    i__1 = jlast;
#line 526 "dlatrs.f"
	    i__2 = jinc;
#line 526 "dlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 530 "dlatrs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 531 "dlatrs.f"
		if (nounit) {
#line 532 "dlatrs.f"
		    tjjs = a[j + j * a_dim1] * tscal;
#line 533 "dlatrs.f"
		} else {
#line 534 "dlatrs.f"
		    tjjs = tscal;
#line 535 "dlatrs.f"
		    if (tscal == 1.) {
#line 535 "dlatrs.f"
			goto L100;
#line 535 "dlatrs.f"
		    }
#line 537 "dlatrs.f"
		}
#line 538 "dlatrs.f"
		tjj = abs(tjjs);
#line 539 "dlatrs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 543 "dlatrs.f"
		    if (tjj < 1.) {
#line 544 "dlatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 548 "dlatrs.f"
			    rec = 1. / xj;
#line 549 "dlatrs.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 550 "dlatrs.f"
			    *scale *= rec;
#line 551 "dlatrs.f"
			    xmax *= rec;
#line 552 "dlatrs.f"
			}
#line 553 "dlatrs.f"
		    }
#line 554 "dlatrs.f"
		    x[j] /= tjjs;
#line 555 "dlatrs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 556 "dlatrs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 560 "dlatrs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 565 "dlatrs.f"
			rec = tjj * bignum / xj;
#line 566 "dlatrs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 571 "dlatrs.f"
			    rec /= cnorm[j];
#line 572 "dlatrs.f"
			}
#line 573 "dlatrs.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 574 "dlatrs.f"
			*scale *= rec;
#line 575 "dlatrs.f"
			xmax *= rec;
#line 576 "dlatrs.f"
		    }
#line 577 "dlatrs.f"
		    x[j] /= tjjs;
#line 578 "dlatrs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 579 "dlatrs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 584 "dlatrs.f"
		    i__3 = *n;
#line 584 "dlatrs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 585 "dlatrs.f"
			x[i__] = 0.;
#line 586 "dlatrs.f"
/* L90: */
#line 586 "dlatrs.f"
		    }
#line 587 "dlatrs.f"
		    x[j] = 1.;
#line 588 "dlatrs.f"
		    xj = 1.;
#line 589 "dlatrs.f"
		    *scale = 0.;
#line 590 "dlatrs.f"
		    xmax = 0.;
#line 591 "dlatrs.f"
		}
#line 592 "dlatrs.f"
L100:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 597 "dlatrs.f"
		if (xj > 1.) {
#line 598 "dlatrs.f"
		    rec = 1. / xj;
#line 599 "dlatrs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 603 "dlatrs.f"
			rec *= .5;
#line 604 "dlatrs.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 605 "dlatrs.f"
			*scale *= rec;
#line 606 "dlatrs.f"
		    }
#line 607 "dlatrs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 611 "dlatrs.f"
		    dscal_(n, &c_b36, &x[1], &c__1);
#line 612 "dlatrs.f"
		    *scale *= .5;
#line 613 "dlatrs.f"
		}

#line 615 "dlatrs.f"
		if (upper) {
#line 616 "dlatrs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 621 "dlatrs.f"
			i__3 = j - 1;
#line 621 "dlatrs.f"
			d__1 = -x[j] * tscal;
#line 621 "dlatrs.f"
			daxpy_(&i__3, &d__1, &a[j * a_dim1 + 1], &c__1, &x[1],
				 &c__1);
#line 623 "dlatrs.f"
			i__3 = j - 1;
#line 623 "dlatrs.f"
			i__ = idamax_(&i__3, &x[1], &c__1);
#line 624 "dlatrs.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 625 "dlatrs.f"
		    }
#line 626 "dlatrs.f"
		} else {
#line 627 "dlatrs.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 632 "dlatrs.f"
			i__3 = *n - j;
#line 632 "dlatrs.f"
			d__1 = -x[j] * tscal;
#line 632 "dlatrs.f"
			daxpy_(&i__3, &d__1, &a[j + 1 + j * a_dim1], &c__1, &
				x[j + 1], &c__1);
#line 634 "dlatrs.f"
			i__3 = *n - j;
#line 634 "dlatrs.f"
			i__ = j + idamax_(&i__3, &x[j + 1], &c__1);
#line 635 "dlatrs.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 636 "dlatrs.f"
		    }
#line 637 "dlatrs.f"
		}
#line 638 "dlatrs.f"
/* L110: */
#line 638 "dlatrs.f"
	    }

#line 640 "dlatrs.f"
	} else {

/*           Solve A**T * x = b */

#line 644 "dlatrs.f"
	    i__2 = jlast;
#line 644 "dlatrs.f"
	    i__1 = jinc;
#line 644 "dlatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 649 "dlatrs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 650 "dlatrs.f"
		uscal = tscal;
#line 651 "dlatrs.f"
		rec = 1. / max(xmax,1.);
#line 652 "dlatrs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 656 "dlatrs.f"
		    rec *= .5;
#line 657 "dlatrs.f"
		    if (nounit) {
#line 658 "dlatrs.f"
			tjjs = a[j + j * a_dim1] * tscal;
#line 659 "dlatrs.f"
		    } else {
#line 660 "dlatrs.f"
			tjjs = tscal;
#line 661 "dlatrs.f"
		    }
#line 662 "dlatrs.f"
		    tjj = abs(tjjs);
#line 663 "dlatrs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 667 "dlatrs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 667 "dlatrs.f"
			rec = min(d__1,d__2);
#line 668 "dlatrs.f"
			uscal /= tjjs;
#line 669 "dlatrs.f"
		    }
#line 670 "dlatrs.f"
		    if (rec < 1.) {
#line 671 "dlatrs.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 672 "dlatrs.f"
			*scale *= rec;
#line 673 "dlatrs.f"
			xmax *= rec;
#line 674 "dlatrs.f"
		    }
#line 675 "dlatrs.f"
		}

#line 677 "dlatrs.f"
		sumj = 0.;
#line 678 "dlatrs.f"
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call DDOT to perform the dot product. */

#line 683 "dlatrs.f"
		    if (upper) {
#line 684 "dlatrs.f"
			i__3 = j - 1;
#line 684 "dlatrs.f"
			sumj = ddot_(&i__3, &a[j * a_dim1 + 1], &c__1, &x[1], 
				&c__1);
#line 685 "dlatrs.f"
		    } else if (j < *n) {
#line 686 "dlatrs.f"
			i__3 = *n - j;
#line 686 "dlatrs.f"
			sumj = ddot_(&i__3, &a[j + 1 + j * a_dim1], &c__1, &x[
				j + 1], &c__1);
#line 687 "dlatrs.f"
		    }
#line 688 "dlatrs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 692 "dlatrs.f"
		    if (upper) {
#line 693 "dlatrs.f"
			i__3 = j - 1;
#line 693 "dlatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 694 "dlatrs.f"
			    sumj += a[i__ + j * a_dim1] * uscal * x[i__];
#line 695 "dlatrs.f"
/* L120: */
#line 695 "dlatrs.f"
			}
#line 696 "dlatrs.f"
		    } else if (j < *n) {
#line 697 "dlatrs.f"
			i__3 = *n;
#line 697 "dlatrs.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 698 "dlatrs.f"
			    sumj += a[i__ + j * a_dim1] * uscal * x[i__];
#line 699 "dlatrs.f"
/* L130: */
#line 699 "dlatrs.f"
			}
#line 700 "dlatrs.f"
		    }
#line 701 "dlatrs.f"
		}

#line 703 "dlatrs.f"
		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 708 "dlatrs.f"
		    x[j] -= sumj;
#line 709 "dlatrs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 710 "dlatrs.f"
		    if (nounit) {
#line 711 "dlatrs.f"
			tjjs = a[j + j * a_dim1] * tscal;
#line 712 "dlatrs.f"
		    } else {
#line 713 "dlatrs.f"
			tjjs = tscal;
#line 714 "dlatrs.f"
			if (tscal == 1.) {
#line 714 "dlatrs.f"
			    goto L150;
#line 714 "dlatrs.f"
			}
#line 716 "dlatrs.f"
		    }

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 720 "dlatrs.f"
		    tjj = abs(tjjs);
#line 721 "dlatrs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 725 "dlatrs.f"
			if (tjj < 1.) {
#line 726 "dlatrs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 730 "dlatrs.f"
				rec = 1. / xj;
#line 731 "dlatrs.f"
				dscal_(n, &rec, &x[1], &c__1);
#line 732 "dlatrs.f"
				*scale *= rec;
#line 733 "dlatrs.f"
				xmax *= rec;
#line 734 "dlatrs.f"
			    }
#line 735 "dlatrs.f"
			}
#line 736 "dlatrs.f"
			x[j] /= tjjs;
#line 737 "dlatrs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 741 "dlatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 745 "dlatrs.f"
			    rec = tjj * bignum / xj;
#line 746 "dlatrs.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 747 "dlatrs.f"
			    *scale *= rec;
#line 748 "dlatrs.f"
			    xmax *= rec;
#line 749 "dlatrs.f"
			}
#line 750 "dlatrs.f"
			x[j] /= tjjs;
#line 751 "dlatrs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0, and compute a solution to A**T*x = 0. */

#line 756 "dlatrs.f"
			i__3 = *n;
#line 756 "dlatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 757 "dlatrs.f"
			    x[i__] = 0.;
#line 758 "dlatrs.f"
/* L140: */
#line 758 "dlatrs.f"
			}
#line 759 "dlatrs.f"
			x[j] = 1.;
#line 760 "dlatrs.f"
			*scale = 0.;
#line 761 "dlatrs.f"
			xmax = 0.;
#line 762 "dlatrs.f"
		    }
#line 763 "dlatrs.f"
L150:
#line 764 "dlatrs.f"
		    ;
#line 764 "dlatrs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 769 "dlatrs.f"
		    x[j] = x[j] / tjjs - sumj;
#line 770 "dlatrs.f"
		}
/* Computing MAX */
#line 771 "dlatrs.f"
		d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
#line 771 "dlatrs.f"
		xmax = max(d__2,d__3);
#line 772 "dlatrs.f"
/* L160: */
#line 772 "dlatrs.f"
	    }
#line 773 "dlatrs.f"
	}
#line 774 "dlatrs.f"
	*scale /= tscal;
#line 775 "dlatrs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 779 "dlatrs.f"
    if (tscal != 1.) {
#line 780 "dlatrs.f"
	d__1 = 1. / tscal;
#line 780 "dlatrs.f"
	dscal_(n, &d__1, &cnorm[1], &c__1);
#line 781 "dlatrs.f"
    }

#line 783 "dlatrs.f"
    return 0;

/*     End of DLATRS */

} /* dlatrs_ */

