#line 1 "dlatbs.f"
/* dlatbs.f -- translated by f2c (version 20100827).
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

#line 1 "dlatbs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b DLATBS solves a triangular banded system of equations. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLATBS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatbs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatbs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatbs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, */
/*                          SCALE, CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), CNORM( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATBS solves one of the triangular systems */
/* > */
/* >    A *x = s*b  or  A**T*x = s*b */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular band matrix.  Here A**T denotes the transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine DTBSV is called.  If the matrix A */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of subdiagonals or superdiagonals in the */
/* >          triangular matrix A.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, DTBSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTBSV if the */
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
/* >  and we can safely call DTBSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublereal *ab, integer *ldab, 
	doublereal *x, doublereal *scale, doublereal *cnorm, integer *info, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax, tjjs, xmax, grow, sumj;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer maind;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal tscal, uscal;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer jlast;
    extern /* Subroutine */ int dtbsv_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notran;
    static integer jfirst;
    static doublereal smlnum;
    static logical nounit;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. Executable Statements .. */

#line 285 "dlatbs.f"
    /* Parameter adjustments */
#line 285 "dlatbs.f"
    ab_dim1 = *ldab;
#line 285 "dlatbs.f"
    ab_offset = 1 + ab_dim1;
#line 285 "dlatbs.f"
    ab -= ab_offset;
#line 285 "dlatbs.f"
    --x;
#line 285 "dlatbs.f"
    --cnorm;
#line 285 "dlatbs.f"

#line 285 "dlatbs.f"
    /* Function Body */
#line 285 "dlatbs.f"
    *info = 0;
#line 286 "dlatbs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 287 "dlatbs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 288 "dlatbs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 292 "dlatbs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 293 "dlatbs.f"
	*info = -1;
#line 294 "dlatbs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 296 "dlatbs.f"
	*info = -2;
#line 297 "dlatbs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 298 "dlatbs.f"
	*info = -3;
#line 299 "dlatbs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 301 "dlatbs.f"
	*info = -4;
#line 302 "dlatbs.f"
    } else if (*n < 0) {
#line 303 "dlatbs.f"
	*info = -5;
#line 304 "dlatbs.f"
    } else if (*kd < 0) {
#line 305 "dlatbs.f"
	*info = -6;
#line 306 "dlatbs.f"
    } else if (*ldab < *kd + 1) {
#line 307 "dlatbs.f"
	*info = -8;
#line 308 "dlatbs.f"
    }
#line 309 "dlatbs.f"
    if (*info != 0) {
#line 310 "dlatbs.f"
	i__1 = -(*info);
#line 310 "dlatbs.f"
	xerbla_("DLATBS", &i__1, (ftnlen)6);
#line 311 "dlatbs.f"
	return 0;
#line 312 "dlatbs.f"
    }

/*     Quick return if possible */

#line 316 "dlatbs.f"
    if (*n == 0) {
#line 316 "dlatbs.f"
	return 0;
#line 316 "dlatbs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 321 "dlatbs.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 322 "dlatbs.f"
    bignum = 1. / smlnum;
#line 323 "dlatbs.f"
    *scale = 1.;

#line 325 "dlatbs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 329 "dlatbs.f"
	if (upper) {

/*           A is upper triangular. */

#line 333 "dlatbs.f"
	    i__1 = *n;
#line 333 "dlatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 334 "dlatbs.f"
		i__2 = *kd, i__3 = j - 1;
#line 334 "dlatbs.f"
		jlen = min(i__2,i__3);
#line 335 "dlatbs.f"
		cnorm[j] = dasum_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], &
			c__1);
#line 336 "dlatbs.f"
/* L10: */
#line 336 "dlatbs.f"
	    }
#line 337 "dlatbs.f"
	} else {

/*           A is lower triangular. */

#line 341 "dlatbs.f"
	    i__1 = *n;
#line 341 "dlatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 342 "dlatbs.f"
		i__2 = *kd, i__3 = *n - j;
#line 342 "dlatbs.f"
		jlen = min(i__2,i__3);
#line 343 "dlatbs.f"
		if (jlen > 0) {
#line 344 "dlatbs.f"
		    cnorm[j] = dasum_(&jlen, &ab[j * ab_dim1 + 2], &c__1);
#line 345 "dlatbs.f"
		} else {
#line 346 "dlatbs.f"
		    cnorm[j] = 0.;
#line 347 "dlatbs.f"
		}
#line 348 "dlatbs.f"
/* L20: */
#line 348 "dlatbs.f"
	    }
#line 349 "dlatbs.f"
	}
#line 350 "dlatbs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM. */

#line 355 "dlatbs.f"
    imax = idamax_(n, &cnorm[1], &c__1);
#line 356 "dlatbs.f"
    tmax = cnorm[imax];
#line 357 "dlatbs.f"
    if (tmax <= bignum) {
#line 358 "dlatbs.f"
	tscal = 1.;
#line 359 "dlatbs.f"
    } else {
#line 360 "dlatbs.f"
	tscal = 1. / (smlnum * tmax);
#line 361 "dlatbs.f"
	dscal_(n, &tscal, &cnorm[1], &c__1);
#line 362 "dlatbs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine DTBSV can be used. */

#line 367 "dlatbs.f"
    j = idamax_(n, &x[1], &c__1);
#line 368 "dlatbs.f"
    xmax = (d__1 = x[j], abs(d__1));
#line 369 "dlatbs.f"
    xbnd = xmax;
#line 370 "dlatbs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 374 "dlatbs.f"
	if (upper) {
#line 375 "dlatbs.f"
	    jfirst = *n;
#line 376 "dlatbs.f"
	    jlast = 1;
#line 377 "dlatbs.f"
	    jinc = -1;
#line 378 "dlatbs.f"
	    maind = *kd + 1;
#line 379 "dlatbs.f"
	} else {
#line 380 "dlatbs.f"
	    jfirst = 1;
#line 381 "dlatbs.f"
	    jlast = *n;
#line 382 "dlatbs.f"
	    jinc = 1;
#line 383 "dlatbs.f"
	    maind = 1;
#line 384 "dlatbs.f"
	}

#line 386 "dlatbs.f"
	if (tscal != 1.) {
#line 387 "dlatbs.f"
	    grow = 0.;
#line 388 "dlatbs.f"
	    goto L50;
#line 389 "dlatbs.f"
	}

#line 391 "dlatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 398 "dlatbs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 399 "dlatbs.f"
	    xbnd = grow;
#line 400 "dlatbs.f"
	    i__1 = jlast;
#line 400 "dlatbs.f"
	    i__2 = jinc;
#line 400 "dlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 404 "dlatbs.f"
		if (grow <= smlnum) {
#line 404 "dlatbs.f"
		    goto L50;
#line 404 "dlatbs.f"
		}

/*              M(j) = G(j-1) / abs(A(j,j)) */

#line 409 "dlatbs.f"
		tjj = (d__1 = ab[maind + j * ab_dim1], abs(d__1));
/* Computing MIN */
#line 410 "dlatbs.f"
		d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 410 "dlatbs.f"
		xbnd = min(d__1,d__2);
#line 411 "dlatbs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 415 "dlatbs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 416 "dlatbs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 420 "dlatbs.f"
		    grow = 0.;
#line 421 "dlatbs.f"
		}
#line 422 "dlatbs.f"
/* L30: */
#line 422 "dlatbs.f"
	    }
#line 423 "dlatbs.f"
	    grow = xbnd;
#line 424 "dlatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 430 "dlatbs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 430 "dlatbs.f"
	    grow = min(d__1,d__2);
#line 431 "dlatbs.f"
	    i__2 = jlast;
#line 431 "dlatbs.f"
	    i__1 = jinc;
#line 431 "dlatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 435 "dlatbs.f"
		if (grow <= smlnum) {
#line 435 "dlatbs.f"
		    goto L50;
#line 435 "dlatbs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 440 "dlatbs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 441 "dlatbs.f"
/* L40: */
#line 441 "dlatbs.f"
	    }
#line 442 "dlatbs.f"
	}
#line 443 "dlatbs.f"
L50:

#line 445 "dlatbs.f"
	;
#line 445 "dlatbs.f"
    } else {

/*        Compute the growth in A**T * x = b. */

#line 449 "dlatbs.f"
	if (upper) {
#line 450 "dlatbs.f"
	    jfirst = 1;
#line 451 "dlatbs.f"
	    jlast = *n;
#line 452 "dlatbs.f"
	    jinc = 1;
#line 453 "dlatbs.f"
	    maind = *kd + 1;
#line 454 "dlatbs.f"
	} else {
#line 455 "dlatbs.f"
	    jfirst = *n;
#line 456 "dlatbs.f"
	    jlast = 1;
#line 457 "dlatbs.f"
	    jinc = -1;
#line 458 "dlatbs.f"
	    maind = 1;
#line 459 "dlatbs.f"
	}

#line 461 "dlatbs.f"
	if (tscal != 1.) {
#line 462 "dlatbs.f"
	    grow = 0.;
#line 463 "dlatbs.f"
	    goto L80;
#line 464 "dlatbs.f"
	}

#line 466 "dlatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 473 "dlatbs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 474 "dlatbs.f"
	    xbnd = grow;
#line 475 "dlatbs.f"
	    i__1 = jlast;
#line 475 "dlatbs.f"
	    i__2 = jinc;
#line 475 "dlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 479 "dlatbs.f"
		if (grow <= smlnum) {
#line 479 "dlatbs.f"
		    goto L80;
#line 479 "dlatbs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 484 "dlatbs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 485 "dlatbs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 485 "dlatbs.f"
		grow = min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 489 "dlatbs.f"
		tjj = (d__1 = ab[maind + j * ab_dim1], abs(d__1));
#line 490 "dlatbs.f"
		if (xj > tjj) {
#line 490 "dlatbs.f"
		    xbnd *= tjj / xj;
#line 490 "dlatbs.f"
		}
#line 492 "dlatbs.f"
/* L60: */
#line 492 "dlatbs.f"
	    }
#line 493 "dlatbs.f"
	    grow = min(grow,xbnd);
#line 494 "dlatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 500 "dlatbs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 500 "dlatbs.f"
	    grow = min(d__1,d__2);
#line 501 "dlatbs.f"
	    i__2 = jlast;
#line 501 "dlatbs.f"
	    i__1 = jinc;
#line 501 "dlatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 505 "dlatbs.f"
		if (grow <= smlnum) {
#line 505 "dlatbs.f"
		    goto L80;
#line 505 "dlatbs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 510 "dlatbs.f"
		xj = cnorm[j] + 1.;
#line 511 "dlatbs.f"
		grow /= xj;
#line 512 "dlatbs.f"
/* L70: */
#line 512 "dlatbs.f"
	    }
#line 513 "dlatbs.f"
	}
#line 514 "dlatbs.f"
L80:
#line 515 "dlatbs.f"
	;
#line 515 "dlatbs.f"
    }

#line 517 "dlatbs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 522 "dlatbs.f"
	dtbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &x[1], &c__1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 523 "dlatbs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 527 "dlatbs.f"
	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 532 "dlatbs.f"
	    *scale = bignum / xmax;
#line 533 "dlatbs.f"
	    dscal_(n, scale, &x[1], &c__1);
#line 534 "dlatbs.f"
	    xmax = bignum;
#line 535 "dlatbs.f"
	}

#line 537 "dlatbs.f"
	if (notran) {

/*           Solve A * x = b */

#line 541 "dlatbs.f"
	    i__1 = jlast;
#line 541 "dlatbs.f"
	    i__2 = jinc;
#line 541 "dlatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 545 "dlatbs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 546 "dlatbs.f"
		if (nounit) {
#line 547 "dlatbs.f"
		    tjjs = ab[maind + j * ab_dim1] * tscal;
#line 548 "dlatbs.f"
		} else {
#line 549 "dlatbs.f"
		    tjjs = tscal;
#line 550 "dlatbs.f"
		    if (tscal == 1.) {
#line 550 "dlatbs.f"
			goto L100;
#line 550 "dlatbs.f"
		    }
#line 552 "dlatbs.f"
		}
#line 553 "dlatbs.f"
		tjj = abs(tjjs);
#line 554 "dlatbs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 558 "dlatbs.f"
		    if (tjj < 1.) {
#line 559 "dlatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 563 "dlatbs.f"
			    rec = 1. / xj;
#line 564 "dlatbs.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 565 "dlatbs.f"
			    *scale *= rec;
#line 566 "dlatbs.f"
			    xmax *= rec;
#line 567 "dlatbs.f"
			}
#line 568 "dlatbs.f"
		    }
#line 569 "dlatbs.f"
		    x[j] /= tjjs;
#line 570 "dlatbs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 571 "dlatbs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 575 "dlatbs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 580 "dlatbs.f"
			rec = tjj * bignum / xj;
#line 581 "dlatbs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 586 "dlatbs.f"
			    rec /= cnorm[j];
#line 587 "dlatbs.f"
			}
#line 588 "dlatbs.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 589 "dlatbs.f"
			*scale *= rec;
#line 590 "dlatbs.f"
			xmax *= rec;
#line 591 "dlatbs.f"
		    }
#line 592 "dlatbs.f"
		    x[j] /= tjjs;
#line 593 "dlatbs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 594 "dlatbs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 599 "dlatbs.f"
		    i__3 = *n;
#line 599 "dlatbs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 600 "dlatbs.f"
			x[i__] = 0.;
#line 601 "dlatbs.f"
/* L90: */
#line 601 "dlatbs.f"
		    }
#line 602 "dlatbs.f"
		    x[j] = 1.;
#line 603 "dlatbs.f"
		    xj = 1.;
#line 604 "dlatbs.f"
		    *scale = 0.;
#line 605 "dlatbs.f"
		    xmax = 0.;
#line 606 "dlatbs.f"
		}
#line 607 "dlatbs.f"
L100:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 612 "dlatbs.f"
		if (xj > 1.) {
#line 613 "dlatbs.f"
		    rec = 1. / xj;
#line 614 "dlatbs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 618 "dlatbs.f"
			rec *= .5;
#line 619 "dlatbs.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 620 "dlatbs.f"
			*scale *= rec;
#line 621 "dlatbs.f"
		    }
#line 622 "dlatbs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 626 "dlatbs.f"
		    dscal_(n, &c_b36, &x[1], &c__1);
#line 627 "dlatbs.f"
		    *scale *= .5;
#line 628 "dlatbs.f"
		}

#line 630 "dlatbs.f"
		if (upper) {
#line 631 "dlatbs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) - */
/*                                             x(j)* A(max(1,j-kd):j-1,j) */

/* Computing MIN */
#line 637 "dlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 637 "dlatbs.f"
			jlen = min(i__3,i__4);
#line 638 "dlatbs.f"
			d__1 = -x[j] * tscal;
#line 638 "dlatbs.f"
			daxpy_(&jlen, &d__1, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 640 "dlatbs.f"
			i__3 = j - 1;
#line 640 "dlatbs.f"
			i__ = idamax_(&i__3, &x[1], &c__1);
#line 641 "dlatbs.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 642 "dlatbs.f"
		    }
#line 643 "dlatbs.f"
		} else if (j < *n) {

/*                 Compute the update */
/*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) - */
/*                                          x(j) * A(j+1:min(j+kd,n),j) */

/* Computing MIN */
#line 649 "dlatbs.f"
		    i__3 = *kd, i__4 = *n - j;
#line 649 "dlatbs.f"
		    jlen = min(i__3,i__4);
#line 650 "dlatbs.f"
		    if (jlen > 0) {
#line 650 "dlatbs.f"
			d__1 = -x[j] * tscal;
#line 650 "dlatbs.f"
			daxpy_(&jlen, &d__1, &ab[j * ab_dim1 + 2], &c__1, &x[
				j + 1], &c__1);
#line 650 "dlatbs.f"
		    }
#line 653 "dlatbs.f"
		    i__3 = *n - j;
#line 653 "dlatbs.f"
		    i__ = j + idamax_(&i__3, &x[j + 1], &c__1);
#line 654 "dlatbs.f"
		    xmax = (d__1 = x[i__], abs(d__1));
#line 655 "dlatbs.f"
		}
#line 656 "dlatbs.f"
/* L110: */
#line 656 "dlatbs.f"
	    }

#line 658 "dlatbs.f"
	} else {

/*           Solve A**T * x = b */

#line 662 "dlatbs.f"
	    i__2 = jlast;
#line 662 "dlatbs.f"
	    i__1 = jinc;
#line 662 "dlatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 667 "dlatbs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 668 "dlatbs.f"
		uscal = tscal;
#line 669 "dlatbs.f"
		rec = 1. / max(xmax,1.);
#line 670 "dlatbs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 674 "dlatbs.f"
		    rec *= .5;
#line 675 "dlatbs.f"
		    if (nounit) {
#line 676 "dlatbs.f"
			tjjs = ab[maind + j * ab_dim1] * tscal;
#line 677 "dlatbs.f"
		    } else {
#line 678 "dlatbs.f"
			tjjs = tscal;
#line 679 "dlatbs.f"
		    }
#line 680 "dlatbs.f"
		    tjj = abs(tjjs);
#line 681 "dlatbs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 685 "dlatbs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 685 "dlatbs.f"
			rec = min(d__1,d__2);
#line 686 "dlatbs.f"
			uscal /= tjjs;
#line 687 "dlatbs.f"
		    }
#line 688 "dlatbs.f"
		    if (rec < 1.) {
#line 689 "dlatbs.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 690 "dlatbs.f"
			*scale *= rec;
#line 691 "dlatbs.f"
			xmax *= rec;
#line 692 "dlatbs.f"
		    }
#line 693 "dlatbs.f"
		}

#line 695 "dlatbs.f"
		sumj = 0.;
#line 696 "dlatbs.f"
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call DDOT to perform the dot product. */

#line 701 "dlatbs.f"
		    if (upper) {
/* Computing MIN */
#line 702 "dlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 702 "dlatbs.f"
			jlen = min(i__3,i__4);
#line 703 "dlatbs.f"
			sumj = ddot_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1],
				 &c__1, &x[j - jlen], &c__1);
#line 705 "dlatbs.f"
		    } else {
/* Computing MIN */
#line 706 "dlatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 706 "dlatbs.f"
			jlen = min(i__3,i__4);
#line 707 "dlatbs.f"
			if (jlen > 0) {
#line 707 "dlatbs.f"
			    sumj = ddot_(&jlen, &ab[j * ab_dim1 + 2], &c__1, &
				    x[j + 1], &c__1);
#line 707 "dlatbs.f"
			}
#line 709 "dlatbs.f"
		    }
#line 710 "dlatbs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 714 "dlatbs.f"
		    if (upper) {
/* Computing MIN */
#line 715 "dlatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 715 "dlatbs.f"
			jlen = min(i__3,i__4);
#line 716 "dlatbs.f"
			i__3 = jlen;
#line 716 "dlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 717 "dlatbs.f"
			    sumj += ab[*kd + i__ - jlen + j * ab_dim1] * 
				    uscal * x[j - jlen - 1 + i__];
#line 719 "dlatbs.f"
/* L120: */
#line 719 "dlatbs.f"
			}
#line 720 "dlatbs.f"
		    } else {
/* Computing MIN */
#line 721 "dlatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 721 "dlatbs.f"
			jlen = min(i__3,i__4);
#line 722 "dlatbs.f"
			i__3 = jlen;
#line 722 "dlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 723 "dlatbs.f"
			    sumj += ab[i__ + 1 + j * ab_dim1] * uscal * x[j + 
				    i__];
#line 724 "dlatbs.f"
/* L130: */
#line 724 "dlatbs.f"
			}
#line 725 "dlatbs.f"
		    }
#line 726 "dlatbs.f"
		}

#line 728 "dlatbs.f"
		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 733 "dlatbs.f"
		    x[j] -= sumj;
#line 734 "dlatbs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 735 "dlatbs.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 739 "dlatbs.f"
			tjjs = ab[maind + j * ab_dim1] * tscal;
#line 740 "dlatbs.f"
		    } else {
#line 741 "dlatbs.f"
			tjjs = tscal;
#line 742 "dlatbs.f"
			if (tscal == 1.) {
#line 742 "dlatbs.f"
			    goto L150;
#line 742 "dlatbs.f"
			}
#line 744 "dlatbs.f"
		    }
#line 745 "dlatbs.f"
		    tjj = abs(tjjs);
#line 746 "dlatbs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 750 "dlatbs.f"
			if (tjj < 1.) {
#line 751 "dlatbs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 755 "dlatbs.f"
				rec = 1. / xj;
#line 756 "dlatbs.f"
				dscal_(n, &rec, &x[1], &c__1);
#line 757 "dlatbs.f"
				*scale *= rec;
#line 758 "dlatbs.f"
				xmax *= rec;
#line 759 "dlatbs.f"
			    }
#line 760 "dlatbs.f"
			}
#line 761 "dlatbs.f"
			x[j] /= tjjs;
#line 762 "dlatbs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 766 "dlatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 770 "dlatbs.f"
			    rec = tjj * bignum / xj;
#line 771 "dlatbs.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 772 "dlatbs.f"
			    *scale *= rec;
#line 773 "dlatbs.f"
			    xmax *= rec;
#line 774 "dlatbs.f"
			}
#line 775 "dlatbs.f"
			x[j] /= tjjs;
#line 776 "dlatbs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0, and compute a solution to A**T*x = 0. */

#line 781 "dlatbs.f"
			i__3 = *n;
#line 781 "dlatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 782 "dlatbs.f"
			    x[i__] = 0.;
#line 783 "dlatbs.f"
/* L140: */
#line 783 "dlatbs.f"
			}
#line 784 "dlatbs.f"
			x[j] = 1.;
#line 785 "dlatbs.f"
			*scale = 0.;
#line 786 "dlatbs.f"
			xmax = 0.;
#line 787 "dlatbs.f"
		    }
#line 788 "dlatbs.f"
L150:
#line 789 "dlatbs.f"
		    ;
#line 789 "dlatbs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - sumj if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 794 "dlatbs.f"
		    x[j] = x[j] / tjjs - sumj;
#line 795 "dlatbs.f"
		}
/* Computing MAX */
#line 796 "dlatbs.f"
		d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
#line 796 "dlatbs.f"
		xmax = max(d__2,d__3);
#line 797 "dlatbs.f"
/* L160: */
#line 797 "dlatbs.f"
	    }
#line 798 "dlatbs.f"
	}
#line 799 "dlatbs.f"
	*scale /= tscal;
#line 800 "dlatbs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 804 "dlatbs.f"
    if (tscal != 1.) {
#line 805 "dlatbs.f"
	d__1 = 1. / tscal;
#line 805 "dlatbs.f"
	dscal_(n, &d__1, &cnorm[1], &c__1);
#line 806 "dlatbs.f"
    }

#line 808 "dlatbs.f"
    return 0;

/*     End of DLATBS */

} /* dlatbs_ */

