#line 1 "zlatrs.f"
/* zlatrs.f -- translated by f2c (version 20100827).
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

#line 1 "zlatrs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b ZLATRS solves a triangular system of equations with the scale factor set to prevent overflow. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLATRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   CNORM( * ) */
/*       COMPLEX*16         A( LDA, * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLATRS solves one of the triangular systems */
/* > */
/* >    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow.  Here A is an upper or lower */
/* > triangular matrix, A**T denotes the transpose of A, A**H denotes the */
/* > conjugate transpose of A, x and b are n-element vectors, and s is a */
/* > scaling factor, usually less than or equal to 1, chosen so that the */
/* > components of x will be less than the overflow threshold.  If the */
/* > unscaled problem will not cause overflow, the Level 2 BLAS routine */
/* > ZTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j), */
/* > then s is set to 0 and a non-trivial solution to A*x = 0 is returned. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, ZTRSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTRSV if the */
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
/* >  and we can safely call ZTRSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, 
	doublereal *scale, doublereal *cnorm, integer *info, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc;
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax;
    static doublecomplex tjjs;
    static doublereal xmax, grow;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
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
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztrsv_(
	    char *, char *, char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), dlabad_(
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 295 "zlatrs.f"
    /* Parameter adjustments */
#line 295 "zlatrs.f"
    a_dim1 = *lda;
#line 295 "zlatrs.f"
    a_offset = 1 + a_dim1;
#line 295 "zlatrs.f"
    a -= a_offset;
#line 295 "zlatrs.f"
    --x;
#line 295 "zlatrs.f"
    --cnorm;
#line 295 "zlatrs.f"

#line 295 "zlatrs.f"
    /* Function Body */
#line 295 "zlatrs.f"
    *info = 0;
#line 296 "zlatrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 297 "zlatrs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 298 "zlatrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 302 "zlatrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 303 "zlatrs.f"
	*info = -1;
#line 304 "zlatrs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 306 "zlatrs.f"
	*info = -2;
#line 307 "zlatrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 308 "zlatrs.f"
	*info = -3;
#line 309 "zlatrs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 311 "zlatrs.f"
	*info = -4;
#line 312 "zlatrs.f"
    } else if (*n < 0) {
#line 313 "zlatrs.f"
	*info = -5;
#line 314 "zlatrs.f"
    } else if (*lda < max(1,*n)) {
#line 315 "zlatrs.f"
	*info = -7;
#line 316 "zlatrs.f"
    }
#line 317 "zlatrs.f"
    if (*info != 0) {
#line 318 "zlatrs.f"
	i__1 = -(*info);
#line 318 "zlatrs.f"
	xerbla_("ZLATRS", &i__1, (ftnlen)6);
#line 319 "zlatrs.f"
	return 0;
#line 320 "zlatrs.f"
    }

/*     Quick return if possible */

#line 324 "zlatrs.f"
    if (*n == 0) {
#line 324 "zlatrs.f"
	return 0;
#line 324 "zlatrs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 329 "zlatrs.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 330 "zlatrs.f"
    bignum = 1. / smlnum;
#line 331 "zlatrs.f"
    dlabad_(&smlnum, &bignum);
#line 332 "zlatrs.f"
    smlnum /= dlamch_("Precision", (ftnlen)9);
#line 333 "zlatrs.f"
    bignum = 1. / smlnum;
#line 334 "zlatrs.f"
    *scale = 1.;

#line 336 "zlatrs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 340 "zlatrs.f"
	if (upper) {

/*           A is upper triangular. */

#line 344 "zlatrs.f"
	    i__1 = *n;
#line 344 "zlatrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 345 "zlatrs.f"
		i__2 = j - 1;
#line 345 "zlatrs.f"
		cnorm[j] = dzasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
#line 346 "zlatrs.f"
/* L10: */
#line 346 "zlatrs.f"
	    }
#line 347 "zlatrs.f"
	} else {

/*           A is lower triangular. */

#line 351 "zlatrs.f"
	    i__1 = *n - 1;
#line 351 "zlatrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 352 "zlatrs.f"
		i__2 = *n - j;
#line 352 "zlatrs.f"
		cnorm[j] = dzasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
#line 353 "zlatrs.f"
/* L20: */
#line 353 "zlatrs.f"
	    }
#line 354 "zlatrs.f"
	    cnorm[*n] = 0.;
#line 355 "zlatrs.f"
	}
#line 356 "zlatrs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM/2. */

#line 361 "zlatrs.f"
    imax = idamax_(n, &cnorm[1], &c__1);
#line 362 "zlatrs.f"
    tmax = cnorm[imax];
#line 363 "zlatrs.f"
    if (tmax <= bignum * .5) {
#line 364 "zlatrs.f"
	tscal = 1.;
#line 365 "zlatrs.f"
    } else {
#line 366 "zlatrs.f"
	tscal = .5 / (smlnum * tmax);
#line 367 "zlatrs.f"
	dscal_(n, &tscal, &cnorm[1], &c__1);
#line 368 "zlatrs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine ZTRSV can be used. */

#line 373 "zlatrs.f"
    xmax = 0.;
#line 374 "zlatrs.f"
    i__1 = *n;
#line 374 "zlatrs.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 375 "zlatrs.f"
	i__2 = j;
#line 375 "zlatrs.f"
	d__3 = xmax, d__4 = (d__1 = x[i__2].r / 2., abs(d__1)) + (d__2 = 
		d_imag(&x[j]) / 2., abs(d__2));
#line 375 "zlatrs.f"
	xmax = max(d__3,d__4);
#line 376 "zlatrs.f"
/* L30: */
#line 376 "zlatrs.f"
    }
#line 377 "zlatrs.f"
    xbnd = xmax;

#line 379 "zlatrs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 383 "zlatrs.f"
	if (upper) {
#line 384 "zlatrs.f"
	    jfirst = *n;
#line 385 "zlatrs.f"
	    jlast = 1;
#line 386 "zlatrs.f"
	    jinc = -1;
#line 387 "zlatrs.f"
	} else {
#line 388 "zlatrs.f"
	    jfirst = 1;
#line 389 "zlatrs.f"
	    jlast = *n;
#line 390 "zlatrs.f"
	    jinc = 1;
#line 391 "zlatrs.f"
	}

#line 393 "zlatrs.f"
	if (tscal != 1.) {
#line 394 "zlatrs.f"
	    grow = 0.;
#line 395 "zlatrs.f"
	    goto L60;
#line 396 "zlatrs.f"
	}

#line 398 "zlatrs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 405 "zlatrs.f"
	    grow = .5 / max(xbnd,smlnum);
#line 406 "zlatrs.f"
	    xbnd = grow;
#line 407 "zlatrs.f"
	    i__1 = jlast;
#line 407 "zlatrs.f"
	    i__2 = jinc;
#line 407 "zlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 411 "zlatrs.f"
		if (grow <= smlnum) {
#line 411 "zlatrs.f"
		    goto L60;
#line 411 "zlatrs.f"
		}

#line 414 "zlatrs.f"
		i__3 = j + j * a_dim1;
#line 414 "zlatrs.f"
		tjjs.r = a[i__3].r, tjjs.i = a[i__3].i;
#line 415 "zlatrs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 417 "zlatrs.f"
		if (tjj >= smlnum) {

/*                 M(j) = G(j-1) / abs(A(j,j)) */

/* Computing MIN */
#line 421 "zlatrs.f"
		    d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 421 "zlatrs.f"
		    xbnd = min(d__1,d__2);
#line 422 "zlatrs.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 426 "zlatrs.f"
		    xbnd = 0.;
#line 427 "zlatrs.f"
		}

#line 429 "zlatrs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 433 "zlatrs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 434 "zlatrs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 438 "zlatrs.f"
		    grow = 0.;
#line 439 "zlatrs.f"
		}
#line 440 "zlatrs.f"
/* L40: */
#line 440 "zlatrs.f"
	    }
#line 441 "zlatrs.f"
	    grow = xbnd;
#line 442 "zlatrs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 448 "zlatrs.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 448 "zlatrs.f"
	    grow = min(d__1,d__2);
#line 449 "zlatrs.f"
	    i__2 = jlast;
#line 449 "zlatrs.f"
	    i__1 = jinc;
#line 449 "zlatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 453 "zlatrs.f"
		if (grow <= smlnum) {
#line 453 "zlatrs.f"
		    goto L60;
#line 453 "zlatrs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 458 "zlatrs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 459 "zlatrs.f"
/* L50: */
#line 459 "zlatrs.f"
	    }
#line 460 "zlatrs.f"
	}
#line 461 "zlatrs.f"
L60:

#line 463 "zlatrs.f"
	;
#line 463 "zlatrs.f"
    } else {

/*        Compute the growth in A**T * x = b  or  A**H * x = b. */

#line 467 "zlatrs.f"
	if (upper) {
#line 468 "zlatrs.f"
	    jfirst = 1;
#line 469 "zlatrs.f"
	    jlast = *n;
#line 470 "zlatrs.f"
	    jinc = 1;
#line 471 "zlatrs.f"
	} else {
#line 472 "zlatrs.f"
	    jfirst = *n;
#line 473 "zlatrs.f"
	    jlast = 1;
#line 474 "zlatrs.f"
	    jinc = -1;
#line 475 "zlatrs.f"
	}

#line 477 "zlatrs.f"
	if (tscal != 1.) {
#line 478 "zlatrs.f"
	    grow = 0.;
#line 479 "zlatrs.f"
	    goto L90;
#line 480 "zlatrs.f"
	}

#line 482 "zlatrs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 489 "zlatrs.f"
	    grow = .5 / max(xbnd,smlnum);
#line 490 "zlatrs.f"
	    xbnd = grow;
#line 491 "zlatrs.f"
	    i__1 = jlast;
#line 491 "zlatrs.f"
	    i__2 = jinc;
#line 491 "zlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 495 "zlatrs.f"
		if (grow <= smlnum) {
#line 495 "zlatrs.f"
		    goto L90;
#line 495 "zlatrs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 500 "zlatrs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 501 "zlatrs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 501 "zlatrs.f"
		grow = min(d__1,d__2);

#line 503 "zlatrs.f"
		i__3 = j + j * a_dim1;
#line 503 "zlatrs.f"
		tjjs.r = a[i__3].r, tjjs.i = a[i__3].i;
#line 504 "zlatrs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 506 "zlatrs.f"
		if (tjj >= smlnum) {

/*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 510 "zlatrs.f"
		    if (xj > tjj) {
#line 510 "zlatrs.f"
			xbnd *= tjj / xj;
#line 510 "zlatrs.f"
		    }
#line 512 "zlatrs.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 516 "zlatrs.f"
		    xbnd = 0.;
#line 517 "zlatrs.f"
		}
#line 518 "zlatrs.f"
/* L70: */
#line 518 "zlatrs.f"
	    }
#line 519 "zlatrs.f"
	    grow = min(grow,xbnd);
#line 520 "zlatrs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 526 "zlatrs.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 526 "zlatrs.f"
	    grow = min(d__1,d__2);
#line 527 "zlatrs.f"
	    i__2 = jlast;
#line 527 "zlatrs.f"
	    i__1 = jinc;
#line 527 "zlatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 531 "zlatrs.f"
		if (grow <= smlnum) {
#line 531 "zlatrs.f"
		    goto L90;
#line 531 "zlatrs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 536 "zlatrs.f"
		xj = cnorm[j] + 1.;
#line 537 "zlatrs.f"
		grow /= xj;
#line 538 "zlatrs.f"
/* L80: */
#line 538 "zlatrs.f"
	    }
#line 539 "zlatrs.f"
	}
#line 540 "zlatrs.f"
L90:
#line 541 "zlatrs.f"
	;
#line 541 "zlatrs.f"
    }

#line 543 "zlatrs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 548 "zlatrs.f"
	ztrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, (ftnlen)
		1, (ftnlen)1, (ftnlen)1);
#line 549 "zlatrs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 553 "zlatrs.f"
	if (xmax > bignum * .5) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 558 "zlatrs.f"
	    *scale = bignum * .5 / xmax;
#line 559 "zlatrs.f"
	    zdscal_(n, scale, &x[1], &c__1);
#line 560 "zlatrs.f"
	    xmax = bignum;
#line 561 "zlatrs.f"
	} else {
#line 562 "zlatrs.f"
	    xmax *= 2.;
#line 563 "zlatrs.f"
	}

#line 565 "zlatrs.f"
	if (notran) {

/*           Solve A * x = b */

#line 569 "zlatrs.f"
	    i__1 = jlast;
#line 569 "zlatrs.f"
	    i__2 = jinc;
#line 569 "zlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 573 "zlatrs.f"
		i__3 = j;
#line 573 "zlatrs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 574 "zlatrs.f"
		if (nounit) {
#line 575 "zlatrs.f"
		    i__3 = j + j * a_dim1;
#line 575 "zlatrs.f"
		    z__1.r = tscal * a[i__3].r, z__1.i = tscal * a[i__3].i;
#line 575 "zlatrs.f"
		    tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 576 "zlatrs.f"
		} else {
#line 577 "zlatrs.f"
		    tjjs.r = tscal, tjjs.i = 0.;
#line 578 "zlatrs.f"
		    if (tscal == 1.) {
#line 578 "zlatrs.f"
			goto L110;
#line 578 "zlatrs.f"
		    }
#line 580 "zlatrs.f"
		}
#line 581 "zlatrs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));
#line 582 "zlatrs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 586 "zlatrs.f"
		    if (tjj < 1.) {
#line 587 "zlatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 591 "zlatrs.f"
			    rec = 1. / xj;
#line 592 "zlatrs.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 593 "zlatrs.f"
			    *scale *= rec;
#line 594 "zlatrs.f"
			    xmax *= rec;
#line 595 "zlatrs.f"
			}
#line 596 "zlatrs.f"
		    }
#line 597 "zlatrs.f"
		    i__3 = j;
#line 597 "zlatrs.f"
		    zladiv_(&z__1, &x[j], &tjjs);
#line 597 "zlatrs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 598 "zlatrs.f"
		    i__3 = j;
#line 598 "zlatrs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 599 "zlatrs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 603 "zlatrs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 608 "zlatrs.f"
			rec = tjj * bignum / xj;
#line 609 "zlatrs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 614 "zlatrs.f"
			    rec /= cnorm[j];
#line 615 "zlatrs.f"
			}
#line 616 "zlatrs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 617 "zlatrs.f"
			*scale *= rec;
#line 618 "zlatrs.f"
			xmax *= rec;
#line 619 "zlatrs.f"
		    }
#line 620 "zlatrs.f"
		    i__3 = j;
#line 620 "zlatrs.f"
		    zladiv_(&z__1, &x[j], &tjjs);
#line 620 "zlatrs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 621 "zlatrs.f"
		    i__3 = j;
#line 621 "zlatrs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 622 "zlatrs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 627 "zlatrs.f"
		    i__3 = *n;
#line 627 "zlatrs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 628 "zlatrs.f"
			i__4 = i__;
#line 628 "zlatrs.f"
			x[i__4].r = 0., x[i__4].i = 0.;
#line 629 "zlatrs.f"
/* L100: */
#line 629 "zlatrs.f"
		    }
#line 630 "zlatrs.f"
		    i__3 = j;
#line 630 "zlatrs.f"
		    x[i__3].r = 1., x[i__3].i = 0.;
#line 631 "zlatrs.f"
		    xj = 1.;
#line 632 "zlatrs.f"
		    *scale = 0.;
#line 633 "zlatrs.f"
		    xmax = 0.;
#line 634 "zlatrs.f"
		}
#line 635 "zlatrs.f"
L110:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 640 "zlatrs.f"
		if (xj > 1.) {
#line 641 "zlatrs.f"
		    rec = 1. / xj;
#line 642 "zlatrs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 646 "zlatrs.f"
			rec *= .5;
#line 647 "zlatrs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 648 "zlatrs.f"
			*scale *= rec;
#line 649 "zlatrs.f"
		    }
#line 650 "zlatrs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 654 "zlatrs.f"
		    zdscal_(n, &c_b36, &x[1], &c__1);
#line 655 "zlatrs.f"
		    *scale *= .5;
#line 656 "zlatrs.f"
		}

#line 658 "zlatrs.f"
		if (upper) {
#line 659 "zlatrs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 664 "zlatrs.f"
			i__3 = j - 1;
#line 664 "zlatrs.f"
			i__4 = j;
#line 664 "zlatrs.f"
			z__2.r = -x[i__4].r, z__2.i = -x[i__4].i;
#line 664 "zlatrs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 664 "zlatrs.f"
			zaxpy_(&i__3, &z__1, &a[j * a_dim1 + 1], &c__1, &x[1],
				 &c__1);
#line 666 "zlatrs.f"
			i__3 = j - 1;
#line 666 "zlatrs.f"
			i__ = izamax_(&i__3, &x[1], &c__1);
#line 667 "zlatrs.f"
			i__3 = i__;
#line 667 "zlatrs.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 668 "zlatrs.f"
		    }
#line 669 "zlatrs.f"
		} else {
#line 670 "zlatrs.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 675 "zlatrs.f"
			i__3 = *n - j;
#line 675 "zlatrs.f"
			i__4 = j;
#line 675 "zlatrs.f"
			z__2.r = -x[i__4].r, z__2.i = -x[i__4].i;
#line 675 "zlatrs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 675 "zlatrs.f"
			zaxpy_(&i__3, &z__1, &a[j + 1 + j * a_dim1], &c__1, &
				x[j + 1], &c__1);
#line 677 "zlatrs.f"
			i__3 = *n - j;
#line 677 "zlatrs.f"
			i__ = j + izamax_(&i__3, &x[j + 1], &c__1);
#line 678 "zlatrs.f"
			i__3 = i__;
#line 678 "zlatrs.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 679 "zlatrs.f"
		    }
#line 680 "zlatrs.f"
		}
#line 681 "zlatrs.f"
/* L120: */
#line 681 "zlatrs.f"
	    }

#line 683 "zlatrs.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Solve A**T * x = b */

#line 687 "zlatrs.f"
	    i__2 = jlast;
#line 687 "zlatrs.f"
	    i__1 = jinc;
#line 687 "zlatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 692 "zlatrs.f"
		i__3 = j;
#line 692 "zlatrs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 693 "zlatrs.f"
		uscal.r = tscal, uscal.i = 0.;
#line 694 "zlatrs.f"
		rec = 1. / max(xmax,1.);
#line 695 "zlatrs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 699 "zlatrs.f"
		    rec *= .5;
#line 700 "zlatrs.f"
		    if (nounit) {
#line 701 "zlatrs.f"
			i__3 = j + j * a_dim1;
#line 701 "zlatrs.f"
			z__1.r = tscal * a[i__3].r, z__1.i = tscal * a[i__3]
				.i;
#line 701 "zlatrs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 702 "zlatrs.f"
		    } else {
#line 703 "zlatrs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 704 "zlatrs.f"
		    }
#line 705 "zlatrs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 706 "zlatrs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 710 "zlatrs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 710 "zlatrs.f"
			rec = min(d__1,d__2);
#line 711 "zlatrs.f"
			zladiv_(&z__1, &uscal, &tjjs);
#line 711 "zlatrs.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 712 "zlatrs.f"
		    }
#line 713 "zlatrs.f"
		    if (rec < 1.) {
#line 714 "zlatrs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 715 "zlatrs.f"
			*scale *= rec;
#line 716 "zlatrs.f"
			xmax *= rec;
#line 717 "zlatrs.f"
		    }
#line 718 "zlatrs.f"
		}

#line 720 "zlatrs.f"
		csumj.r = 0., csumj.i = 0.;
#line 721 "zlatrs.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call ZDOTU to perform the dot product. */

#line 726 "zlatrs.f"
		    if (upper) {
#line 727 "zlatrs.f"
			i__3 = j - 1;
#line 727 "zlatrs.f"
			zdotu_(&z__1, &i__3, &a[j * a_dim1 + 1], &c__1, &x[1],
				 &c__1);
#line 727 "zlatrs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 728 "zlatrs.f"
		    } else if (j < *n) {
#line 729 "zlatrs.f"
			i__3 = *n - j;
#line 729 "zlatrs.f"
			zdotu_(&z__1, &i__3, &a[j + 1 + j * a_dim1], &c__1, &
				x[j + 1], &c__1);
#line 729 "zlatrs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 730 "zlatrs.f"
		    }
#line 731 "zlatrs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 735 "zlatrs.f"
		    if (upper) {
#line 736 "zlatrs.f"
			i__3 = j - 1;
#line 736 "zlatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 737 "zlatrs.f"
			    i__4 = i__ + j * a_dim1;
#line 737 "zlatrs.f"
			    z__3.r = a[i__4].r * uscal.r - a[i__4].i * 
				    uscal.i, z__3.i = a[i__4].r * uscal.i + a[
				    i__4].i * uscal.r;
#line 737 "zlatrs.f"
			    i__5 = i__;
#line 737 "zlatrs.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 737 "zlatrs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 737 "zlatrs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 738 "zlatrs.f"
/* L130: */
#line 738 "zlatrs.f"
			}
#line 739 "zlatrs.f"
		    } else if (j < *n) {
#line 740 "zlatrs.f"
			i__3 = *n;
#line 740 "zlatrs.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 741 "zlatrs.f"
			    i__4 = i__ + j * a_dim1;
#line 741 "zlatrs.f"
			    z__3.r = a[i__4].r * uscal.r - a[i__4].i * 
				    uscal.i, z__3.i = a[i__4].r * uscal.i + a[
				    i__4].i * uscal.r;
#line 741 "zlatrs.f"
			    i__5 = i__;
#line 741 "zlatrs.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 741 "zlatrs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 741 "zlatrs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 742 "zlatrs.f"
/* L140: */
#line 742 "zlatrs.f"
			}
#line 743 "zlatrs.f"
		    }
#line 744 "zlatrs.f"
		}

#line 746 "zlatrs.f"
		z__1.r = tscal, z__1.i = 0.;
#line 746 "zlatrs.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 751 "zlatrs.f"
		    i__3 = j;
#line 751 "zlatrs.f"
		    i__4 = j;
#line 751 "zlatrs.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 751 "zlatrs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 752 "zlatrs.f"
		    i__3 = j;
#line 752 "zlatrs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 753 "zlatrs.f"
		    if (nounit) {
#line 754 "zlatrs.f"
			i__3 = j + j * a_dim1;
#line 754 "zlatrs.f"
			z__1.r = tscal * a[i__3].r, z__1.i = tscal * a[i__3]
				.i;
#line 754 "zlatrs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 755 "zlatrs.f"
		    } else {
#line 756 "zlatrs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 757 "zlatrs.f"
			if (tscal == 1.) {
#line 757 "zlatrs.f"
			    goto L160;
#line 757 "zlatrs.f"
			}
#line 759 "zlatrs.f"
		    }

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 763 "zlatrs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 764 "zlatrs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 768 "zlatrs.f"
			if (tjj < 1.) {
#line 769 "zlatrs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 773 "zlatrs.f"
				rec = 1. / xj;
#line 774 "zlatrs.f"
				zdscal_(n, &rec, &x[1], &c__1);
#line 775 "zlatrs.f"
				*scale *= rec;
#line 776 "zlatrs.f"
				xmax *= rec;
#line 777 "zlatrs.f"
			    }
#line 778 "zlatrs.f"
			}
#line 779 "zlatrs.f"
			i__3 = j;
#line 779 "zlatrs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 779 "zlatrs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 780 "zlatrs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 784 "zlatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 788 "zlatrs.f"
			    rec = tjj * bignum / xj;
#line 789 "zlatrs.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 790 "zlatrs.f"
			    *scale *= rec;
#line 791 "zlatrs.f"
			    xmax *= rec;
#line 792 "zlatrs.f"
			}
#line 793 "zlatrs.f"
			i__3 = j;
#line 793 "zlatrs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 793 "zlatrs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 794 "zlatrs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**T *x = 0. */

#line 799 "zlatrs.f"
			i__3 = *n;
#line 799 "zlatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 800 "zlatrs.f"
			    i__4 = i__;
#line 800 "zlatrs.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 801 "zlatrs.f"
/* L150: */
#line 801 "zlatrs.f"
			}
#line 802 "zlatrs.f"
			i__3 = j;
#line 802 "zlatrs.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 803 "zlatrs.f"
			*scale = 0.;
#line 804 "zlatrs.f"
			xmax = 0.;
#line 805 "zlatrs.f"
		    }
#line 806 "zlatrs.f"
L160:
#line 807 "zlatrs.f"
		    ;
#line 807 "zlatrs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 812 "zlatrs.f"
		    i__3 = j;
#line 812 "zlatrs.f"
		    zladiv_(&z__2, &x[j], &tjjs);
#line 812 "zlatrs.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 812 "zlatrs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 813 "zlatrs.f"
		}
/* Computing MAX */
#line 814 "zlatrs.f"
		i__3 = j;
#line 814 "zlatrs.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 814 "zlatrs.f"
		xmax = max(d__3,d__4);
#line 815 "zlatrs.f"
/* L170: */
#line 815 "zlatrs.f"
	    }

#line 817 "zlatrs.f"
	} else {

/*           Solve A**H * x = b */

#line 821 "zlatrs.f"
	    i__1 = jlast;
#line 821 "zlatrs.f"
	    i__2 = jinc;
#line 821 "zlatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 826 "zlatrs.f"
		i__3 = j;
#line 826 "zlatrs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 827 "zlatrs.f"
		uscal.r = tscal, uscal.i = 0.;
#line 828 "zlatrs.f"
		rec = 1. / max(xmax,1.);
#line 829 "zlatrs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 833 "zlatrs.f"
		    rec *= .5;
#line 834 "zlatrs.f"
		    if (nounit) {
#line 835 "zlatrs.f"
			d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 835 "zlatrs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 835 "zlatrs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 836 "zlatrs.f"
		    } else {
#line 837 "zlatrs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 838 "zlatrs.f"
		    }
#line 839 "zlatrs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 840 "zlatrs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 844 "zlatrs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 844 "zlatrs.f"
			rec = min(d__1,d__2);
#line 845 "zlatrs.f"
			zladiv_(&z__1, &uscal, &tjjs);
#line 845 "zlatrs.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 846 "zlatrs.f"
		    }
#line 847 "zlatrs.f"
		    if (rec < 1.) {
#line 848 "zlatrs.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 849 "zlatrs.f"
			*scale *= rec;
#line 850 "zlatrs.f"
			xmax *= rec;
#line 851 "zlatrs.f"
		    }
#line 852 "zlatrs.f"
		}

#line 854 "zlatrs.f"
		csumj.r = 0., csumj.i = 0.;
#line 855 "zlatrs.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call ZDOTC to perform the dot product. */

#line 860 "zlatrs.f"
		    if (upper) {
#line 861 "zlatrs.f"
			i__3 = j - 1;
#line 861 "zlatrs.f"
			zdotc_(&z__1, &i__3, &a[j * a_dim1 + 1], &c__1, &x[1],
				 &c__1);
#line 861 "zlatrs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 862 "zlatrs.f"
		    } else if (j < *n) {
#line 863 "zlatrs.f"
			i__3 = *n - j;
#line 863 "zlatrs.f"
			zdotc_(&z__1, &i__3, &a[j + 1 + j * a_dim1], &c__1, &
				x[j + 1], &c__1);
#line 863 "zlatrs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 864 "zlatrs.f"
		    }
#line 865 "zlatrs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 869 "zlatrs.f"
		    if (upper) {
#line 870 "zlatrs.f"
			i__3 = j - 1;
#line 870 "zlatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 871 "zlatrs.f"
			    d_cnjg(&z__4, &a[i__ + j * a_dim1]);
#line 871 "zlatrs.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 871 "zlatrs.f"
			    i__4 = i__;
#line 871 "zlatrs.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 871 "zlatrs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 871 "zlatrs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 873 "zlatrs.f"
/* L180: */
#line 873 "zlatrs.f"
			}
#line 874 "zlatrs.f"
		    } else if (j < *n) {
#line 875 "zlatrs.f"
			i__3 = *n;
#line 875 "zlatrs.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 876 "zlatrs.f"
			    d_cnjg(&z__4, &a[i__ + j * a_dim1]);
#line 876 "zlatrs.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 876 "zlatrs.f"
			    i__4 = i__;
#line 876 "zlatrs.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 876 "zlatrs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 876 "zlatrs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 878 "zlatrs.f"
/* L190: */
#line 878 "zlatrs.f"
			}
#line 879 "zlatrs.f"
		    }
#line 880 "zlatrs.f"
		}

#line 882 "zlatrs.f"
		z__1.r = tscal, z__1.i = 0.;
#line 882 "zlatrs.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 887 "zlatrs.f"
		    i__3 = j;
#line 887 "zlatrs.f"
		    i__4 = j;
#line 887 "zlatrs.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 887 "zlatrs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 888 "zlatrs.f"
		    i__3 = j;
#line 888 "zlatrs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 889 "zlatrs.f"
		    if (nounit) {
#line 890 "zlatrs.f"
			d_cnjg(&z__2, &a[j + j * a_dim1]);
#line 890 "zlatrs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 890 "zlatrs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 891 "zlatrs.f"
		    } else {
#line 892 "zlatrs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 893 "zlatrs.f"
			if (tscal == 1.) {
#line 893 "zlatrs.f"
			    goto L210;
#line 893 "zlatrs.f"
			}
#line 895 "zlatrs.f"
		    }

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 899 "zlatrs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 900 "zlatrs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 904 "zlatrs.f"
			if (tjj < 1.) {
#line 905 "zlatrs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 909 "zlatrs.f"
				rec = 1. / xj;
#line 910 "zlatrs.f"
				zdscal_(n, &rec, &x[1], &c__1);
#line 911 "zlatrs.f"
				*scale *= rec;
#line 912 "zlatrs.f"
				xmax *= rec;
#line 913 "zlatrs.f"
			    }
#line 914 "zlatrs.f"
			}
#line 915 "zlatrs.f"
			i__3 = j;
#line 915 "zlatrs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 915 "zlatrs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 916 "zlatrs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 920 "zlatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 924 "zlatrs.f"
			    rec = tjj * bignum / xj;
#line 925 "zlatrs.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 926 "zlatrs.f"
			    *scale *= rec;
#line 927 "zlatrs.f"
			    xmax *= rec;
#line 928 "zlatrs.f"
			}
#line 929 "zlatrs.f"
			i__3 = j;
#line 929 "zlatrs.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 929 "zlatrs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 930 "zlatrs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**H *x = 0. */

#line 935 "zlatrs.f"
			i__3 = *n;
#line 935 "zlatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 936 "zlatrs.f"
			    i__4 = i__;
#line 936 "zlatrs.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 937 "zlatrs.f"
/* L200: */
#line 937 "zlatrs.f"
			}
#line 938 "zlatrs.f"
			i__3 = j;
#line 938 "zlatrs.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 939 "zlatrs.f"
			*scale = 0.;
#line 940 "zlatrs.f"
			xmax = 0.;
#line 941 "zlatrs.f"
		    }
#line 942 "zlatrs.f"
L210:
#line 943 "zlatrs.f"
		    ;
#line 943 "zlatrs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 948 "zlatrs.f"
		    i__3 = j;
#line 948 "zlatrs.f"
		    zladiv_(&z__2, &x[j], &tjjs);
#line 948 "zlatrs.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 948 "zlatrs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 949 "zlatrs.f"
		}
/* Computing MAX */
#line 950 "zlatrs.f"
		i__3 = j;
#line 950 "zlatrs.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 950 "zlatrs.f"
		xmax = max(d__3,d__4);
#line 951 "zlatrs.f"
/* L220: */
#line 951 "zlatrs.f"
	    }
#line 952 "zlatrs.f"
	}
#line 953 "zlatrs.f"
	*scale /= tscal;
#line 954 "zlatrs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 958 "zlatrs.f"
    if (tscal != 1.) {
#line 959 "zlatrs.f"
	d__1 = 1. / tscal;
#line 959 "zlatrs.f"
	dscal_(n, &d__1, &cnorm[1], &c__1);
#line 960 "zlatrs.f"
    }

#line 962 "zlatrs.f"
    return 0;

/*     End of ZLATRS */

} /* zlatrs_ */

