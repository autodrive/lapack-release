#line 1 "clatps.f"
/* clatps.f -- translated by f2c (version 20100827).
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

#line 1 "clatps.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b CLATPS solves a triangular system of equations with the matrix held in packed storage. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLATPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               CNORM( * ) */
/*       COMPLEX            AP( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLATPS solves one of the triangular systems */
/* > */
/* >    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular matrix stored in packed form.  Here A**T denotes the */
/* > transpose of A, A**H denotes the conjugate transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine CTPSV is called. If the matrix A */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b. */
/* >          If SCALE = 0, the matrix A is singular or badly scaled, and */
/* >          the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* >          CNORM is REAL array, dimension (N) */
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

/* > \ingroup complexOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, CTPSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTPSV if the */
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
/* >  and we can safely call CTPSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *
	scale, doublereal *cnorm, integer *info, ftnlen uplo_len, ftnlen 
	trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ip;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax;
    static doublecomplex tjjs;
    static doublereal xmax, grow;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal tscal;
    static doublecomplex uscal;
    static integer jlast;
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublecomplex csumj;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ctpsv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen), slabad_(doublereal *, doublereal *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */ VOID cladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 287 "clatps.f"
    /* Parameter adjustments */
#line 287 "clatps.f"
    --cnorm;
#line 287 "clatps.f"
    --x;
#line 287 "clatps.f"
    --ap;
#line 287 "clatps.f"

#line 287 "clatps.f"
    /* Function Body */
#line 287 "clatps.f"
    *info = 0;
#line 288 "clatps.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 289 "clatps.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 290 "clatps.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 294 "clatps.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 295 "clatps.f"
	*info = -1;
#line 296 "clatps.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 298 "clatps.f"
	*info = -2;
#line 299 "clatps.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 300 "clatps.f"
	*info = -3;
#line 301 "clatps.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 303 "clatps.f"
	*info = -4;
#line 304 "clatps.f"
    } else if (*n < 0) {
#line 305 "clatps.f"
	*info = -5;
#line 306 "clatps.f"
    }
#line 307 "clatps.f"
    if (*info != 0) {
#line 308 "clatps.f"
	i__1 = -(*info);
#line 308 "clatps.f"
	xerbla_("CLATPS", &i__1, (ftnlen)6);
#line 309 "clatps.f"
	return 0;
#line 310 "clatps.f"
    }

/*     Quick return if possible */

#line 314 "clatps.f"
    if (*n == 0) {
#line 314 "clatps.f"
	return 0;
#line 314 "clatps.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 319 "clatps.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 320 "clatps.f"
    bignum = 1. / smlnum;
#line 321 "clatps.f"
    slabad_(&smlnum, &bignum);
#line 322 "clatps.f"
    smlnum /= slamch_("Precision", (ftnlen)9);
#line 323 "clatps.f"
    bignum = 1. / smlnum;
#line 324 "clatps.f"
    *scale = 1.;

#line 326 "clatps.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 330 "clatps.f"
	if (upper) {

/*           A is upper triangular. */

#line 334 "clatps.f"
	    ip = 1;
#line 335 "clatps.f"
	    i__1 = *n;
#line 335 "clatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 336 "clatps.f"
		i__2 = j - 1;
#line 336 "clatps.f"
		cnorm[j] = scasum_(&i__2, &ap[ip], &c__1);
#line 337 "clatps.f"
		ip += j;
#line 338 "clatps.f"
/* L10: */
#line 338 "clatps.f"
	    }
#line 339 "clatps.f"
	} else {

/*           A is lower triangular. */

#line 343 "clatps.f"
	    ip = 1;
#line 344 "clatps.f"
	    i__1 = *n - 1;
#line 344 "clatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 345 "clatps.f"
		i__2 = *n - j;
#line 345 "clatps.f"
		cnorm[j] = scasum_(&i__2, &ap[ip + 1], &c__1);
#line 346 "clatps.f"
		ip = ip + *n - j + 1;
#line 347 "clatps.f"
/* L20: */
#line 347 "clatps.f"
	    }
#line 348 "clatps.f"
	    cnorm[*n] = 0.;
#line 349 "clatps.f"
	}
#line 350 "clatps.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM/2. */

#line 355 "clatps.f"
    imax = isamax_(n, &cnorm[1], &c__1);
#line 356 "clatps.f"
    tmax = cnorm[imax];
#line 357 "clatps.f"
    if (tmax <= bignum * .5) {
#line 358 "clatps.f"
	tscal = 1.;
#line 359 "clatps.f"
    } else {
#line 360 "clatps.f"
	tscal = .5 / (smlnum * tmax);
#line 361 "clatps.f"
	sscal_(n, &tscal, &cnorm[1], &c__1);
#line 362 "clatps.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine CTPSV can be used. */

#line 367 "clatps.f"
    xmax = 0.;
#line 368 "clatps.f"
    i__1 = *n;
#line 368 "clatps.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 369 "clatps.f"
	i__2 = j;
#line 369 "clatps.f"
	d__3 = xmax, d__4 = (d__1 = x[i__2].r / 2., abs(d__1)) + (d__2 = 
		d_imag(&x[j]) / 2., abs(d__2));
#line 369 "clatps.f"
	xmax = max(d__3,d__4);
#line 370 "clatps.f"
/* L30: */
#line 370 "clatps.f"
    }
#line 371 "clatps.f"
    xbnd = xmax;
#line 372 "clatps.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 376 "clatps.f"
	if (upper) {
#line 377 "clatps.f"
	    jfirst = *n;
#line 378 "clatps.f"
	    jlast = 1;
#line 379 "clatps.f"
	    jinc = -1;
#line 380 "clatps.f"
	} else {
#line 381 "clatps.f"
	    jfirst = 1;
#line 382 "clatps.f"
	    jlast = *n;
#line 383 "clatps.f"
	    jinc = 1;
#line 384 "clatps.f"
	}

#line 386 "clatps.f"
	if (tscal != 1.) {
#line 387 "clatps.f"
	    grow = 0.;
#line 388 "clatps.f"
	    goto L60;
#line 389 "clatps.f"
	}

#line 391 "clatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 398 "clatps.f"
	    grow = .5 / max(xbnd,smlnum);
#line 399 "clatps.f"
	    xbnd = grow;
#line 400 "clatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 401 "clatps.f"
	    jlen = *n;
#line 402 "clatps.f"
	    i__1 = jlast;
#line 402 "clatps.f"
	    i__2 = jinc;
#line 402 "clatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 406 "clatps.f"
		if (grow <= smlnum) {
#line 406 "clatps.f"
		    goto L60;
#line 406 "clatps.f"
		}

#line 409 "clatps.f"
		i__3 = ip;
#line 409 "clatps.f"
		tjjs.r = ap[i__3].r, tjjs.i = ap[i__3].i;
#line 410 "clatps.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 412 "clatps.f"
		if (tjj >= smlnum) {

/*                 M(j) = G(j-1) / abs(A(j,j)) */

/* Computing MIN */
#line 416 "clatps.f"
		    d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 416 "clatps.f"
		    xbnd = min(d__1,d__2);
#line 417 "clatps.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 421 "clatps.f"
		    xbnd = 0.;
#line 422 "clatps.f"
		}

#line 424 "clatps.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 428 "clatps.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 429 "clatps.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 433 "clatps.f"
		    grow = 0.;
#line 434 "clatps.f"
		}
#line 435 "clatps.f"
		ip += jinc * jlen;
#line 436 "clatps.f"
		--jlen;
#line 437 "clatps.f"
/* L40: */
#line 437 "clatps.f"
	    }
#line 438 "clatps.f"
	    grow = xbnd;
#line 439 "clatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 445 "clatps.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 445 "clatps.f"
	    grow = min(d__1,d__2);
#line 446 "clatps.f"
	    i__2 = jlast;
#line 446 "clatps.f"
	    i__1 = jinc;
#line 446 "clatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 450 "clatps.f"
		if (grow <= smlnum) {
#line 450 "clatps.f"
		    goto L60;
#line 450 "clatps.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 455 "clatps.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 456 "clatps.f"
/* L50: */
#line 456 "clatps.f"
	    }
#line 457 "clatps.f"
	}
#line 458 "clatps.f"
L60:

#line 460 "clatps.f"
	;
#line 460 "clatps.f"
    } else {

/*        Compute the growth in A**T * x = b  or  A**H * x = b. */

#line 464 "clatps.f"
	if (upper) {
#line 465 "clatps.f"
	    jfirst = 1;
#line 466 "clatps.f"
	    jlast = *n;
#line 467 "clatps.f"
	    jinc = 1;
#line 468 "clatps.f"
	} else {
#line 469 "clatps.f"
	    jfirst = *n;
#line 470 "clatps.f"
	    jlast = 1;
#line 471 "clatps.f"
	    jinc = -1;
#line 472 "clatps.f"
	}

#line 474 "clatps.f"
	if (tscal != 1.) {
#line 475 "clatps.f"
	    grow = 0.;
#line 476 "clatps.f"
	    goto L90;
#line 477 "clatps.f"
	}

#line 479 "clatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 486 "clatps.f"
	    grow = .5 / max(xbnd,smlnum);
#line 487 "clatps.f"
	    xbnd = grow;
#line 488 "clatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 489 "clatps.f"
	    jlen = 1;
#line 490 "clatps.f"
	    i__1 = jlast;
#line 490 "clatps.f"
	    i__2 = jinc;
#line 490 "clatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 494 "clatps.f"
		if (grow <= smlnum) {
#line 494 "clatps.f"
		    goto L90;
#line 494 "clatps.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 499 "clatps.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 500 "clatps.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 500 "clatps.f"
		grow = min(d__1,d__2);

#line 502 "clatps.f"
		i__3 = ip;
#line 502 "clatps.f"
		tjjs.r = ap[i__3].r, tjjs.i = ap[i__3].i;
#line 503 "clatps.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 505 "clatps.f"
		if (tjj >= smlnum) {

/*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 509 "clatps.f"
		    if (xj > tjj) {
#line 509 "clatps.f"
			xbnd *= tjj / xj;
#line 509 "clatps.f"
		    }
#line 511 "clatps.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 515 "clatps.f"
		    xbnd = 0.;
#line 516 "clatps.f"
		}
#line 517 "clatps.f"
		++jlen;
#line 518 "clatps.f"
		ip += jinc * jlen;
#line 519 "clatps.f"
/* L70: */
#line 519 "clatps.f"
	    }
#line 520 "clatps.f"
	    grow = min(grow,xbnd);
#line 521 "clatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 527 "clatps.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 527 "clatps.f"
	    grow = min(d__1,d__2);
#line 528 "clatps.f"
	    i__2 = jlast;
#line 528 "clatps.f"
	    i__1 = jinc;
#line 528 "clatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 532 "clatps.f"
		if (grow <= smlnum) {
#line 532 "clatps.f"
		    goto L90;
#line 532 "clatps.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 537 "clatps.f"
		xj = cnorm[j] + 1.;
#line 538 "clatps.f"
		grow /= xj;
#line 539 "clatps.f"
/* L80: */
#line 539 "clatps.f"
	    }
#line 540 "clatps.f"
	}
#line 541 "clatps.f"
L90:
#line 542 "clatps.f"
	;
#line 542 "clatps.f"
    }

#line 544 "clatps.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 549 "clatps.f"
	ctpsv_(uplo, trans, diag, n, &ap[1], &x[1], &c__1, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
#line 550 "clatps.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 554 "clatps.f"
	if (xmax > bignum * .5) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 559 "clatps.f"
	    *scale = bignum * .5 / xmax;
#line 560 "clatps.f"
	    csscal_(n, scale, &x[1], &c__1);
#line 561 "clatps.f"
	    xmax = bignum;
#line 562 "clatps.f"
	} else {
#line 563 "clatps.f"
	    xmax *= 2.;
#line 564 "clatps.f"
	}

#line 566 "clatps.f"
	if (notran) {

/*           Solve A * x = b */

#line 570 "clatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 571 "clatps.f"
	    i__1 = jlast;
#line 571 "clatps.f"
	    i__2 = jinc;
#line 571 "clatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 575 "clatps.f"
		i__3 = j;
#line 575 "clatps.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 576 "clatps.f"
		if (nounit) {
#line 577 "clatps.f"
		    i__3 = ip;
#line 577 "clatps.f"
		    z__1.r = tscal * ap[i__3].r, z__1.i = tscal * ap[i__3].i;
#line 577 "clatps.f"
		    tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 578 "clatps.f"
		} else {
#line 579 "clatps.f"
		    tjjs.r = tscal, tjjs.i = 0.;
#line 580 "clatps.f"
		    if (tscal == 1.) {
#line 580 "clatps.f"
			goto L105;
#line 580 "clatps.f"
		    }
#line 582 "clatps.f"
		}
#line 583 "clatps.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));
#line 584 "clatps.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 588 "clatps.f"
		    if (tjj < 1.) {
#line 589 "clatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 593 "clatps.f"
			    rec = 1. / xj;
#line 594 "clatps.f"
			    csscal_(n, &rec, &x[1], &c__1);
#line 595 "clatps.f"
			    *scale *= rec;
#line 596 "clatps.f"
			    xmax *= rec;
#line 597 "clatps.f"
			}
#line 598 "clatps.f"
		    }
#line 599 "clatps.f"
		    i__3 = j;
#line 599 "clatps.f"
		    cladiv_(&z__1, &x[j], &tjjs);
#line 599 "clatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 600 "clatps.f"
		    i__3 = j;
#line 600 "clatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 601 "clatps.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 605 "clatps.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 610 "clatps.f"
			rec = tjj * bignum / xj;
#line 611 "clatps.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 616 "clatps.f"
			    rec /= cnorm[j];
#line 617 "clatps.f"
			}
#line 618 "clatps.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 619 "clatps.f"
			*scale *= rec;
#line 620 "clatps.f"
			xmax *= rec;
#line 621 "clatps.f"
		    }
#line 622 "clatps.f"
		    i__3 = j;
#line 622 "clatps.f"
		    cladiv_(&z__1, &x[j], &tjjs);
#line 622 "clatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 623 "clatps.f"
		    i__3 = j;
#line 623 "clatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 624 "clatps.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 629 "clatps.f"
		    i__3 = *n;
#line 629 "clatps.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 630 "clatps.f"
			i__4 = i__;
#line 630 "clatps.f"
			x[i__4].r = 0., x[i__4].i = 0.;
#line 631 "clatps.f"
/* L100: */
#line 631 "clatps.f"
		    }
#line 632 "clatps.f"
		    i__3 = j;
#line 632 "clatps.f"
		    x[i__3].r = 1., x[i__3].i = 0.;
#line 633 "clatps.f"
		    xj = 1.;
#line 634 "clatps.f"
		    *scale = 0.;
#line 635 "clatps.f"
		    xmax = 0.;
#line 636 "clatps.f"
		}
#line 637 "clatps.f"
L105:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 642 "clatps.f"
		if (xj > 1.) {
#line 643 "clatps.f"
		    rec = 1. / xj;
#line 644 "clatps.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 648 "clatps.f"
			rec *= .5;
#line 649 "clatps.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 650 "clatps.f"
			*scale *= rec;
#line 651 "clatps.f"
		    }
#line 652 "clatps.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 656 "clatps.f"
		    csscal_(n, &c_b36, &x[1], &c__1);
#line 657 "clatps.f"
		    *scale *= .5;
#line 658 "clatps.f"
		}

#line 660 "clatps.f"
		if (upper) {
#line 661 "clatps.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 666 "clatps.f"
			i__3 = j - 1;
#line 666 "clatps.f"
			i__4 = j;
#line 666 "clatps.f"
			z__2.r = -x[i__4].r, z__2.i = -x[i__4].i;
#line 666 "clatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 666 "clatps.f"
			caxpy_(&i__3, &z__1, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 668 "clatps.f"
			i__3 = j - 1;
#line 668 "clatps.f"
			i__ = icamax_(&i__3, &x[1], &c__1);
#line 669 "clatps.f"
			i__3 = i__;
#line 669 "clatps.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 670 "clatps.f"
		    }
#line 671 "clatps.f"
		    ip -= j;
#line 672 "clatps.f"
		} else {
#line 673 "clatps.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 678 "clatps.f"
			i__3 = *n - j;
#line 678 "clatps.f"
			i__4 = j;
#line 678 "clatps.f"
			z__2.r = -x[i__4].r, z__2.i = -x[i__4].i;
#line 678 "clatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 678 "clatps.f"
			caxpy_(&i__3, &z__1, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 680 "clatps.f"
			i__3 = *n - j;
#line 680 "clatps.f"
			i__ = j + icamax_(&i__3, &x[j + 1], &c__1);
#line 681 "clatps.f"
			i__3 = i__;
#line 681 "clatps.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 682 "clatps.f"
		    }
#line 683 "clatps.f"
		    ip = ip + *n - j + 1;
#line 684 "clatps.f"
		}
#line 685 "clatps.f"
/* L110: */
#line 685 "clatps.f"
	    }

#line 687 "clatps.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Solve A**T * x = b */

#line 691 "clatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 692 "clatps.f"
	    jlen = 1;
#line 693 "clatps.f"
	    i__2 = jlast;
#line 693 "clatps.f"
	    i__1 = jinc;
#line 693 "clatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 698 "clatps.f"
		i__3 = j;
#line 698 "clatps.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 699 "clatps.f"
		uscal.r = tscal, uscal.i = 0.;
#line 700 "clatps.f"
		rec = 1. / max(xmax,1.);
#line 701 "clatps.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 705 "clatps.f"
		    rec *= .5;
#line 706 "clatps.f"
		    if (nounit) {
#line 707 "clatps.f"
			i__3 = ip;
#line 707 "clatps.f"
			z__1.r = tscal * ap[i__3].r, z__1.i = tscal * ap[i__3]
				.i;
#line 707 "clatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 708 "clatps.f"
		    } else {
#line 709 "clatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 710 "clatps.f"
		    }
#line 711 "clatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 712 "clatps.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 716 "clatps.f"
			d__1 = 1., d__2 = rec * tjj;
#line 716 "clatps.f"
			rec = min(d__1,d__2);
#line 717 "clatps.f"
			cladiv_(&z__1, &uscal, &tjjs);
#line 717 "clatps.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 718 "clatps.f"
		    }
#line 719 "clatps.f"
		    if (rec < 1.) {
#line 720 "clatps.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 721 "clatps.f"
			*scale *= rec;
#line 722 "clatps.f"
			xmax *= rec;
#line 723 "clatps.f"
		    }
#line 724 "clatps.f"
		}

#line 726 "clatps.f"
		csumj.r = 0., csumj.i = 0.;
#line 727 "clatps.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call CDOTU to perform the dot product. */

#line 732 "clatps.f"
		    if (upper) {
#line 733 "clatps.f"
			i__3 = j - 1;
#line 733 "clatps.f"
			cdotu_(&z__1, &i__3, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 733 "clatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 734 "clatps.f"
		    } else if (j < *n) {
#line 735 "clatps.f"
			i__3 = *n - j;
#line 735 "clatps.f"
			cdotu_(&z__1, &i__3, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 735 "clatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 736 "clatps.f"
		    }
#line 737 "clatps.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 741 "clatps.f"
		    if (upper) {
#line 742 "clatps.f"
			i__3 = j - 1;
#line 742 "clatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 743 "clatps.f"
			    i__4 = ip - j + i__;
#line 743 "clatps.f"
			    z__3.r = ap[i__4].r * uscal.r - ap[i__4].i * 
				    uscal.i, z__3.i = ap[i__4].r * uscal.i + 
				    ap[i__4].i * uscal.r;
#line 743 "clatps.f"
			    i__5 = i__;
#line 743 "clatps.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 743 "clatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 743 "clatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 744 "clatps.f"
/* L120: */
#line 744 "clatps.f"
			}
#line 745 "clatps.f"
		    } else if (j < *n) {
#line 746 "clatps.f"
			i__3 = *n - j;
#line 746 "clatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 747 "clatps.f"
			    i__4 = ip + i__;
#line 747 "clatps.f"
			    z__3.r = ap[i__4].r * uscal.r - ap[i__4].i * 
				    uscal.i, z__3.i = ap[i__4].r * uscal.i + 
				    ap[i__4].i * uscal.r;
#line 747 "clatps.f"
			    i__5 = j + i__;
#line 747 "clatps.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 747 "clatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 747 "clatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 748 "clatps.f"
/* L130: */
#line 748 "clatps.f"
			}
#line 749 "clatps.f"
		    }
#line 750 "clatps.f"
		}

#line 752 "clatps.f"
		z__1.r = tscal, z__1.i = 0.;
#line 752 "clatps.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 757 "clatps.f"
		    i__3 = j;
#line 757 "clatps.f"
		    i__4 = j;
#line 757 "clatps.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 757 "clatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 758 "clatps.f"
		    i__3 = j;
#line 758 "clatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 759 "clatps.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 763 "clatps.f"
			i__3 = ip;
#line 763 "clatps.f"
			z__1.r = tscal * ap[i__3].r, z__1.i = tscal * ap[i__3]
				.i;
#line 763 "clatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 764 "clatps.f"
		    } else {
#line 765 "clatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 766 "clatps.f"
			if (tscal == 1.) {
#line 766 "clatps.f"
			    goto L145;
#line 766 "clatps.f"
			}
#line 768 "clatps.f"
		    }
#line 769 "clatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 770 "clatps.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 774 "clatps.f"
			if (tjj < 1.) {
#line 775 "clatps.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 779 "clatps.f"
				rec = 1. / xj;
#line 780 "clatps.f"
				csscal_(n, &rec, &x[1], &c__1);
#line 781 "clatps.f"
				*scale *= rec;
#line 782 "clatps.f"
				xmax *= rec;
#line 783 "clatps.f"
			    }
#line 784 "clatps.f"
			}
#line 785 "clatps.f"
			i__3 = j;
#line 785 "clatps.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 785 "clatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 786 "clatps.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 790 "clatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 794 "clatps.f"
			    rec = tjj * bignum / xj;
#line 795 "clatps.f"
			    csscal_(n, &rec, &x[1], &c__1);
#line 796 "clatps.f"
			    *scale *= rec;
#line 797 "clatps.f"
			    xmax *= rec;
#line 798 "clatps.f"
			}
#line 799 "clatps.f"
			i__3 = j;
#line 799 "clatps.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 799 "clatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 800 "clatps.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**T *x = 0. */

#line 805 "clatps.f"
			i__3 = *n;
#line 805 "clatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 806 "clatps.f"
			    i__4 = i__;
#line 806 "clatps.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 807 "clatps.f"
/* L140: */
#line 807 "clatps.f"
			}
#line 808 "clatps.f"
			i__3 = j;
#line 808 "clatps.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 809 "clatps.f"
			*scale = 0.;
#line 810 "clatps.f"
			xmax = 0.;
#line 811 "clatps.f"
		    }
#line 812 "clatps.f"
L145:
#line 813 "clatps.f"
		    ;
#line 813 "clatps.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 818 "clatps.f"
		    i__3 = j;
#line 818 "clatps.f"
		    cladiv_(&z__2, &x[j], &tjjs);
#line 818 "clatps.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 818 "clatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 819 "clatps.f"
		}
/* Computing MAX */
#line 820 "clatps.f"
		i__3 = j;
#line 820 "clatps.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 820 "clatps.f"
		xmax = max(d__3,d__4);
#line 821 "clatps.f"
		++jlen;
#line 822 "clatps.f"
		ip += jinc * jlen;
#line 823 "clatps.f"
/* L150: */
#line 823 "clatps.f"
	    }

#line 825 "clatps.f"
	} else {

/*           Solve A**H * x = b */

#line 829 "clatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 830 "clatps.f"
	    jlen = 1;
#line 831 "clatps.f"
	    i__1 = jlast;
#line 831 "clatps.f"
	    i__2 = jinc;
#line 831 "clatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 836 "clatps.f"
		i__3 = j;
#line 836 "clatps.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 837 "clatps.f"
		uscal.r = tscal, uscal.i = 0.;
#line 838 "clatps.f"
		rec = 1. / max(xmax,1.);
#line 839 "clatps.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 843 "clatps.f"
		    rec *= .5;
#line 844 "clatps.f"
		    if (nounit) {
#line 845 "clatps.f"
			d_cnjg(&z__2, &ap[ip]);
#line 845 "clatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 845 "clatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 846 "clatps.f"
		    } else {
#line 847 "clatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 848 "clatps.f"
		    }
#line 849 "clatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 850 "clatps.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 854 "clatps.f"
			d__1 = 1., d__2 = rec * tjj;
#line 854 "clatps.f"
			rec = min(d__1,d__2);
#line 855 "clatps.f"
			cladiv_(&z__1, &uscal, &tjjs);
#line 855 "clatps.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 856 "clatps.f"
		    }
#line 857 "clatps.f"
		    if (rec < 1.) {
#line 858 "clatps.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 859 "clatps.f"
			*scale *= rec;
#line 860 "clatps.f"
			xmax *= rec;
#line 861 "clatps.f"
		    }
#line 862 "clatps.f"
		}

#line 864 "clatps.f"
		csumj.r = 0., csumj.i = 0.;
#line 865 "clatps.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call CDOTC to perform the dot product. */

#line 870 "clatps.f"
		    if (upper) {
#line 871 "clatps.f"
			i__3 = j - 1;
#line 871 "clatps.f"
			cdotc_(&z__1, &i__3, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 871 "clatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 872 "clatps.f"
		    } else if (j < *n) {
#line 873 "clatps.f"
			i__3 = *n - j;
#line 873 "clatps.f"
			cdotc_(&z__1, &i__3, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 873 "clatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 874 "clatps.f"
		    }
#line 875 "clatps.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 879 "clatps.f"
		    if (upper) {
#line 880 "clatps.f"
			i__3 = j - 1;
#line 880 "clatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 881 "clatps.f"
			    d_cnjg(&z__4, &ap[ip - j + i__]);
#line 881 "clatps.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 881 "clatps.f"
			    i__4 = i__;
#line 881 "clatps.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 881 "clatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 881 "clatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 883 "clatps.f"
/* L160: */
#line 883 "clatps.f"
			}
#line 884 "clatps.f"
		    } else if (j < *n) {
#line 885 "clatps.f"
			i__3 = *n - j;
#line 885 "clatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 886 "clatps.f"
			    d_cnjg(&z__4, &ap[ip + i__]);
#line 886 "clatps.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 886 "clatps.f"
			    i__4 = j + i__;
#line 886 "clatps.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 886 "clatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 886 "clatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 888 "clatps.f"
/* L170: */
#line 888 "clatps.f"
			}
#line 889 "clatps.f"
		    }
#line 890 "clatps.f"
		}

#line 892 "clatps.f"
		z__1.r = tscal, z__1.i = 0.;
#line 892 "clatps.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 897 "clatps.f"
		    i__3 = j;
#line 897 "clatps.f"
		    i__4 = j;
#line 897 "clatps.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 897 "clatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 898 "clatps.f"
		    i__3 = j;
#line 898 "clatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 899 "clatps.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 903 "clatps.f"
			d_cnjg(&z__2, &ap[ip]);
#line 903 "clatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 903 "clatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 904 "clatps.f"
		    } else {
#line 905 "clatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 906 "clatps.f"
			if (tscal == 1.) {
#line 906 "clatps.f"
			    goto L185;
#line 906 "clatps.f"
			}
#line 908 "clatps.f"
		    }
#line 909 "clatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 910 "clatps.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 914 "clatps.f"
			if (tjj < 1.) {
#line 915 "clatps.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 919 "clatps.f"
				rec = 1. / xj;
#line 920 "clatps.f"
				csscal_(n, &rec, &x[1], &c__1);
#line 921 "clatps.f"
				*scale *= rec;
#line 922 "clatps.f"
				xmax *= rec;
#line 923 "clatps.f"
			    }
#line 924 "clatps.f"
			}
#line 925 "clatps.f"
			i__3 = j;
#line 925 "clatps.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 925 "clatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 926 "clatps.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 930 "clatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 934 "clatps.f"
			    rec = tjj * bignum / xj;
#line 935 "clatps.f"
			    csscal_(n, &rec, &x[1], &c__1);
#line 936 "clatps.f"
			    *scale *= rec;
#line 937 "clatps.f"
			    xmax *= rec;
#line 938 "clatps.f"
			}
#line 939 "clatps.f"
			i__3 = j;
#line 939 "clatps.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 939 "clatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 940 "clatps.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**H *x = 0. */

#line 945 "clatps.f"
			i__3 = *n;
#line 945 "clatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 946 "clatps.f"
			    i__4 = i__;
#line 946 "clatps.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 947 "clatps.f"
/* L180: */
#line 947 "clatps.f"
			}
#line 948 "clatps.f"
			i__3 = j;
#line 948 "clatps.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 949 "clatps.f"
			*scale = 0.;
#line 950 "clatps.f"
			xmax = 0.;
#line 951 "clatps.f"
		    }
#line 952 "clatps.f"
L185:
#line 953 "clatps.f"
		    ;
#line 953 "clatps.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 958 "clatps.f"
		    i__3 = j;
#line 958 "clatps.f"
		    cladiv_(&z__2, &x[j], &tjjs);
#line 958 "clatps.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 958 "clatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 959 "clatps.f"
		}
/* Computing MAX */
#line 960 "clatps.f"
		i__3 = j;
#line 960 "clatps.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 960 "clatps.f"
		xmax = max(d__3,d__4);
#line 961 "clatps.f"
		++jlen;
#line 962 "clatps.f"
		ip += jinc * jlen;
#line 963 "clatps.f"
/* L190: */
#line 963 "clatps.f"
	    }
#line 964 "clatps.f"
	}
#line 965 "clatps.f"
	*scale /= tscal;
#line 966 "clatps.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 970 "clatps.f"
    if (tscal != 1.) {
#line 971 "clatps.f"
	d__1 = 1. / tscal;
#line 971 "clatps.f"
	sscal_(n, &d__1, &cnorm[1], &c__1);
#line 972 "clatps.f"
    }

#line 974 "clatps.f"
    return 0;

/*     End of CLATPS */

} /* clatps_ */

