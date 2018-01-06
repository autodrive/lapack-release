#line 1 "zlals0.f"
/* zlals0.f -- translated by f2c (version 20100827).
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

#line 1 "zlals0.f"
/* Table of constant values */

static doublereal c_b5 = -1.;
static integer c__1 = 1;
static doublereal c_b13 = 1.;
static doublereal c_b15 = 0.;
static integer c__0 = 0;

/* > \brief \b ZLALS0 applies back multiplying factors in solving the least squares problem using divide and c
onquer SVD approach. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLALS0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlals0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlals0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlals0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, */
/*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, */
/*                          POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, */
/*      $                   LDGNUM, NL, NR, NRHS, SQRE */
/*       DOUBLE PRECISION   C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), PERM( * ) */
/*       DOUBLE PRECISION   DIFL( * ), DIFR( LDGNUM, * ), */
/*      $                   GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), */
/*      $                   RWORK( * ), Z( * ) */
/*       COMPLEX*16         B( LDB, * ), BX( LDBX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLALS0 applies back the multiplying factors of either the left or the */
/* > right singular vector matrix of a diagonal matrix appended by a row */
/* > to the right hand side matrix B in solving the least squares problem */
/* > using the divide-and-conquer SVD approach. */
/* > */
/* > For the left singular vector matrix, three types of orthogonal */
/* > matrices are involved: */
/* > */
/* > (1L) Givens rotations: the number of such rotations is GIVPTR; the */
/* >      pairs of columns/rows they were applied to are stored in GIVCOL; */
/* >      and the C- and S-values of these rotations are stored in GIVNUM. */
/* > */
/* > (2L) Permutation. The (NL+1)-st row of B is to be moved to the first */
/* >      row, and for J=2:N, PERM(J)-th row of B is to be moved to the */
/* >      J-th row. */
/* > */
/* > (3L) The left singular vector matrix of the remaining matrix. */
/* > */
/* > For the right singular vector matrix, four types of orthogonal */
/* > matrices are involved: */
/* > */
/* > (1R) The right singular vector matrix of the remaining matrix. */
/* > */
/* > (2R) If SQRE = 1, one extra Givens rotation to generate the right */
/* >      null space. */
/* > */
/* > (3R) The inverse transformation of (2L). */
/* > */
/* > (4R) The inverse transformation of (1L). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >         Specifies whether singular vectors are to be computed in */
/* >         factored form: */
/* >         = 0: Left singular vector matrix. */
/* >         = 1: Right singular vector matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] NL */
/* > \verbatim */
/* >          NL is INTEGER */
/* >         The row dimension of the upper block. NL >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NR */
/* > \verbatim */
/* >          NR is INTEGER */
/* >         The row dimension of the lower block. NR >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* >          SQRE is INTEGER */
/* >         = 0: the lower block is an NR-by-NR square matrix. */
/* >         = 1: the lower block is an NR-by-(NR+1) rectangular matrix. */
/* > */
/* >         The bidiagonal matrix has row dimension N = NL + NR + 1, */
/* >         and column dimension M = N + SQRE. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >         The number of columns of B and BX. NRHS must be at least 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension ( LDB, NRHS ) */
/* >         On input, B contains the right hand sides of the least */
/* >         squares problem in rows 1 through M. On output, B contains */
/* >         the solution X in rows 1 through N. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >         The leading dimension of B. LDB must be at least */
/* >         max(1,MAX( M, N ) ). */
/* > \endverbatim */
/* > */
/* > \param[out] BX */
/* > \verbatim */
/* >          BX is COMPLEX*16 array, dimension ( LDBX, NRHS ) */
/* > \endverbatim */
/* > */
/* > \param[in] LDBX */
/* > \verbatim */
/* >          LDBX is INTEGER */
/* >         The leading dimension of BX. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension ( N ) */
/* >         The permutations (from deflation and sorting) applied */
/* >         to the two blocks. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER */
/* >         The number of Givens rotations which took place in this */
/* >         subproblem. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 ) */
/* >         Each pair of numbers indicates a pair of rows/columns */
/* >         involved in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* >          LDGCOL is INTEGER */
/* >         The leading dimension of GIVCOL, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ) */
/* >         Each number indicates the C or S value used in the */
/* >         corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGNUM */
/* > \verbatim */
/* >          LDGNUM is INTEGER */
/* >         The leading dimension of arrays DIFR, POLES and */
/* >         GIVNUM, must be at least K. */
/* > \endverbatim */
/* > */
/* > \param[in] POLES */
/* > \verbatim */
/* >          POLES is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ) */
/* >         On entry, POLES(1:K, 1) contains the new singular */
/* >         values obtained from solving the secular equation, and */
/* >         POLES(1:K, 2) is an array containing the poles in the secular */
/* >         equation. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFL */
/* > \verbatim */
/* >          DIFL is DOUBLE PRECISION array, dimension ( K ). */
/* >         On entry, DIFL(I) is the distance between I-th updated */
/* >         (undeflated) singular value and the I-th (undeflated) old */
/* >         singular value. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFR */
/* > \verbatim */
/* >          DIFR is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ). */
/* >         On entry, DIFR(I, 1) contains the distances between I-th */
/* >         updated (undeflated) singular value and the I+1-th */
/* >         (undeflated) old singular value. And DIFR(I, 2) is the */
/* >         normalizing factor for the I-th right singular vector. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( K ) */
/* >         Contain the components of the deflation-adjusted updating row */
/* >         vector. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         Contains the dimension of the non-deflated matrix, */
/* >         This is the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >         C contains garbage if SQRE =0 and the C-value of a Givens */
/* >         rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* >         S contains garbage if SQRE =0 and the S-value of a Givens */
/* >         rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension */
/* >         ( K*(1+NRHS) + 2*NRHS ) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int zlals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublecomplex *b, integer *ldb, 
	doublecomplex *bx, integer *ldbx, integer *perm, integer *givptr, 
	integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum,
	 doublereal *poles, doublereal *difl, doublereal *difr, doublereal *
	z__, integer *k, doublereal *c__, doublereal *s, doublereal *rwork, 
	integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, difr_dim1, difr_offset, givnum_dim1, 
	    givnum_offset, poles_dim1, poles_offset, b_dim1, b_offset, 
	    bx_dim1, bx_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, m, n;
    static doublereal dj;
    static integer nlp1, jcol;
    static doublereal temp;
    static integer jrow;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal diflj, difrj, dsigj;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), zdrot_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal dsigjp;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), zlascl_(char *, integer *, integer *,
	     doublereal *, doublereal *, integer *, integer *, doublecomplex *
	    , integer *, integer *, ftnlen), zlacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 316 "zlals0.f"
    /* Parameter adjustments */
#line 316 "zlals0.f"
    b_dim1 = *ldb;
#line 316 "zlals0.f"
    b_offset = 1 + b_dim1;
#line 316 "zlals0.f"
    b -= b_offset;
#line 316 "zlals0.f"
    bx_dim1 = *ldbx;
#line 316 "zlals0.f"
    bx_offset = 1 + bx_dim1;
#line 316 "zlals0.f"
    bx -= bx_offset;
#line 316 "zlals0.f"
    --perm;
#line 316 "zlals0.f"
    givcol_dim1 = *ldgcol;
#line 316 "zlals0.f"
    givcol_offset = 1 + givcol_dim1;
#line 316 "zlals0.f"
    givcol -= givcol_offset;
#line 316 "zlals0.f"
    difr_dim1 = *ldgnum;
#line 316 "zlals0.f"
    difr_offset = 1 + difr_dim1;
#line 316 "zlals0.f"
    difr -= difr_offset;
#line 316 "zlals0.f"
    poles_dim1 = *ldgnum;
#line 316 "zlals0.f"
    poles_offset = 1 + poles_dim1;
#line 316 "zlals0.f"
    poles -= poles_offset;
#line 316 "zlals0.f"
    givnum_dim1 = *ldgnum;
#line 316 "zlals0.f"
    givnum_offset = 1 + givnum_dim1;
#line 316 "zlals0.f"
    givnum -= givnum_offset;
#line 316 "zlals0.f"
    --difl;
#line 316 "zlals0.f"
    --z__;
#line 316 "zlals0.f"
    --rwork;
#line 316 "zlals0.f"

#line 316 "zlals0.f"
    /* Function Body */
#line 316 "zlals0.f"
    *info = 0;

#line 318 "zlals0.f"
    if (*icompq < 0 || *icompq > 1) {
#line 319 "zlals0.f"
	*info = -1;
#line 320 "zlals0.f"
    } else if (*nl < 1) {
#line 321 "zlals0.f"
	*info = -2;
#line 322 "zlals0.f"
    } else if (*nr < 1) {
#line 323 "zlals0.f"
	*info = -3;
#line 324 "zlals0.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 325 "zlals0.f"
	*info = -4;
#line 326 "zlals0.f"
    }

#line 328 "zlals0.f"
    n = *nl + *nr + 1;

#line 330 "zlals0.f"
    if (*nrhs < 1) {
#line 331 "zlals0.f"
	*info = -5;
#line 332 "zlals0.f"
    } else if (*ldb < n) {
#line 333 "zlals0.f"
	*info = -7;
#line 334 "zlals0.f"
    } else if (*ldbx < n) {
#line 335 "zlals0.f"
	*info = -9;
#line 336 "zlals0.f"
    } else if (*givptr < 0) {
#line 337 "zlals0.f"
	*info = -11;
#line 338 "zlals0.f"
    } else if (*ldgcol < n) {
#line 339 "zlals0.f"
	*info = -13;
#line 340 "zlals0.f"
    } else if (*ldgnum < n) {
#line 341 "zlals0.f"
	*info = -15;
#line 342 "zlals0.f"
    } else if (*k < 1) {
#line 343 "zlals0.f"
	*info = -20;
#line 344 "zlals0.f"
    }
#line 345 "zlals0.f"
    if (*info != 0) {
#line 346 "zlals0.f"
	i__1 = -(*info);
#line 346 "zlals0.f"
	xerbla_("ZLALS0", &i__1, (ftnlen)6);
#line 347 "zlals0.f"
	return 0;
#line 348 "zlals0.f"
    }

#line 350 "zlals0.f"
    m = n + *sqre;
#line 351 "zlals0.f"
    nlp1 = *nl + 1;

#line 353 "zlals0.f"
    if (*icompq == 0) {

/*        Apply back orthogonal transformations from the left. */

/*        Step (1L): apply back the Givens rotations performed. */

#line 359 "zlals0.f"
	i__1 = *givptr;
#line 359 "zlals0.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 360 "zlals0.f"
	    zdrot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb, &
		    b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + 
		    (givnum_dim1 << 1)], &givnum[i__ + givnum_dim1]);
#line 363 "zlals0.f"
/* L10: */
#line 363 "zlals0.f"
	}

/*        Step (2L): permute rows of B. */

#line 367 "zlals0.f"
	zcopy_(nrhs, &b[nlp1 + b_dim1], ldb, &bx[bx_dim1 + 1], ldbx);
#line 368 "zlals0.f"
	i__1 = n;
#line 368 "zlals0.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 369 "zlals0.f"
	    zcopy_(nrhs, &b[perm[i__] + b_dim1], ldb, &bx[i__ + bx_dim1], 
		    ldbx);
#line 370 "zlals0.f"
/* L20: */
#line 370 "zlals0.f"
	}

/*        Step (3L): apply the inverse of the left singular vector */
/*        matrix to BX. */

#line 375 "zlals0.f"
	if (*k == 1) {
#line 376 "zlals0.f"
	    zcopy_(nrhs, &bx[bx_offset], ldbx, &b[b_offset], ldb);
#line 377 "zlals0.f"
	    if (z__[1] < 0.) {
#line 378 "zlals0.f"
		zdscal_(nrhs, &c_b5, &b[b_offset], ldb);
#line 379 "zlals0.f"
	    }
#line 380 "zlals0.f"
	} else {
#line 381 "zlals0.f"
	    i__1 = *k;
#line 381 "zlals0.f"
	    for (j = 1; j <= i__1; ++j) {
#line 382 "zlals0.f"
		diflj = difl[j];
#line 383 "zlals0.f"
		dj = poles[j + poles_dim1];
#line 384 "zlals0.f"
		dsigj = -poles[j + (poles_dim1 << 1)];
#line 385 "zlals0.f"
		if (j < *k) {
#line 386 "zlals0.f"
		    difrj = -difr[j + difr_dim1];
#line 387 "zlals0.f"
		    dsigjp = -poles[j + 1 + (poles_dim1 << 1)];
#line 388 "zlals0.f"
		}
#line 389 "zlals0.f"
		if (z__[j] == 0. || poles[j + (poles_dim1 << 1)] == 0.) {
#line 391 "zlals0.f"
		    rwork[j] = 0.;
#line 392 "zlals0.f"
		} else {
#line 393 "zlals0.f"
		    rwork[j] = -poles[j + (poles_dim1 << 1)] * z__[j] / diflj 
			    / (poles[j + (poles_dim1 << 1)] + dj);
#line 395 "zlals0.f"
		}
#line 396 "zlals0.f"
		i__2 = j - 1;
#line 396 "zlals0.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 397 "zlals0.f"
		    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 
			    0.) {
#line 399 "zlals0.f"
			rwork[i__] = 0.;
#line 400 "zlals0.f"
		    } else {
#line 401 "zlals0.f"
			rwork[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__]
				 / (dlamc3_(&poles[i__ + (poles_dim1 << 1)], &
				dsigj) - diflj) / (poles[i__ + (poles_dim1 << 
				1)] + dj);
#line 404 "zlals0.f"
		    }
#line 405 "zlals0.f"
/* L30: */
#line 405 "zlals0.f"
		}
#line 406 "zlals0.f"
		i__2 = *k;
#line 406 "zlals0.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 407 "zlals0.f"
		    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 
			    0.) {
#line 409 "zlals0.f"
			rwork[i__] = 0.;
#line 410 "zlals0.f"
		    } else {
#line 411 "zlals0.f"
			rwork[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__]
				 / (dlamc3_(&poles[i__ + (poles_dim1 << 1)], &
				dsigjp) + difrj) / (poles[i__ + (poles_dim1 <<
				 1)] + dj);
#line 414 "zlals0.f"
		    }
#line 415 "zlals0.f"
/* L40: */
#line 415 "zlals0.f"
		}
#line 416 "zlals0.f"
		rwork[1] = -1.;
#line 417 "zlals0.f"
		temp = dnrm2_(k, &rwork[1], &c__1);

/*              Since B and BX are complex, the following call to DGEMV */
/*              is performed in two steps (real and imaginary parts). */

/*              CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO, */
/*    $                     B( J, 1 ), LDB ) */

#line 425 "zlals0.f"
		i__ = *k + (*nrhs << 1);
#line 426 "zlals0.f"
		i__2 = *nrhs;
#line 426 "zlals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 427 "zlals0.f"
		    i__3 = *k;
#line 427 "zlals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 428 "zlals0.f"
			++i__;
#line 429 "zlals0.f"
			i__4 = jrow + jcol * bx_dim1;
#line 429 "zlals0.f"
			rwork[i__] = bx[i__4].r;
#line 430 "zlals0.f"
/* L50: */
#line 430 "zlals0.f"
		    }
#line 431 "zlals0.f"
/* L60: */
#line 431 "zlals0.f"
		}
#line 432 "zlals0.f"
		dgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1], &c__1, (
			ftnlen)1);
#line 434 "zlals0.f"
		i__ = *k + (*nrhs << 1);
#line 435 "zlals0.f"
		i__2 = *nrhs;
#line 435 "zlals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 436 "zlals0.f"
		    i__3 = *k;
#line 436 "zlals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 437 "zlals0.f"
			++i__;
#line 438 "zlals0.f"
			rwork[i__] = d_imag(&bx[jrow + jcol * bx_dim1]);
#line 439 "zlals0.f"
/* L70: */
#line 439 "zlals0.f"
		    }
#line 440 "zlals0.f"
/* L80: */
#line 440 "zlals0.f"
		}
#line 441 "zlals0.f"
		dgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1 + *nrhs], &
			c__1, (ftnlen)1);
#line 443 "zlals0.f"
		i__2 = *nrhs;
#line 443 "zlals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 444 "zlals0.f"
		    i__3 = j + jcol * b_dim1;
#line 444 "zlals0.f"
		    i__4 = jcol + *k;
#line 444 "zlals0.f"
		    i__5 = jcol + *k + *nrhs;
#line 444 "zlals0.f"
		    z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 444 "zlals0.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 446 "zlals0.f"
/* L90: */
#line 446 "zlals0.f"
		}
#line 447 "zlals0.f"
		zlascl_("G", &c__0, &c__0, &temp, &c_b13, &c__1, nrhs, &b[j + 
			b_dim1], ldb, info, (ftnlen)1);
#line 449 "zlals0.f"
/* L100: */
#line 449 "zlals0.f"
	    }
#line 450 "zlals0.f"
	}

/*        Move the deflated rows of BX to B also. */

#line 454 "zlals0.f"
	if (*k < max(m,n)) {
#line 454 "zlals0.f"
	    i__1 = n - *k;
#line 454 "zlals0.f"
	    zlacpy_("A", &i__1, nrhs, &bx[*k + 1 + bx_dim1], ldbx, &b[*k + 1 
		    + b_dim1], ldb, (ftnlen)1);
#line 454 "zlals0.f"
	}
#line 457 "zlals0.f"
    } else {

/*        Apply back the right orthogonal transformations. */

/*        Step (1R): apply back the new right singular vector matrix */
/*        to B. */

#line 464 "zlals0.f"
	if (*k == 1) {
#line 465 "zlals0.f"
	    zcopy_(nrhs, &b[b_offset], ldb, &bx[bx_offset], ldbx);
#line 466 "zlals0.f"
	} else {
#line 467 "zlals0.f"
	    i__1 = *k;
#line 467 "zlals0.f"
	    for (j = 1; j <= i__1; ++j) {
#line 468 "zlals0.f"
		dsigj = poles[j + (poles_dim1 << 1)];
#line 469 "zlals0.f"
		if (z__[j] == 0.) {
#line 470 "zlals0.f"
		    rwork[j] = 0.;
#line 471 "zlals0.f"
		} else {
#line 472 "zlals0.f"
		    rwork[j] = -z__[j] / difl[j] / (dsigj + poles[j + 
			    poles_dim1]) / difr[j + (difr_dim1 << 1)];
#line 474 "zlals0.f"
		}
#line 475 "zlals0.f"
		i__2 = j - 1;
#line 475 "zlals0.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 476 "zlals0.f"
		    if (z__[j] == 0.) {
#line 477 "zlals0.f"
			rwork[i__] = 0.;
#line 478 "zlals0.f"
		    } else {
#line 479 "zlals0.f"
			d__1 = -poles[i__ + 1 + (poles_dim1 << 1)];
#line 479 "zlals0.f"
			rwork[i__] = z__[j] / (dlamc3_(&dsigj, &d__1) - difr[
				i__ + difr_dim1]) / (dsigj + poles[i__ + 
				poles_dim1]) / difr[i__ + (difr_dim1 << 1)];
#line 482 "zlals0.f"
		    }
#line 483 "zlals0.f"
/* L110: */
#line 483 "zlals0.f"
		}
#line 484 "zlals0.f"
		i__2 = *k;
#line 484 "zlals0.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 485 "zlals0.f"
		    if (z__[j] == 0.) {
#line 486 "zlals0.f"
			rwork[i__] = 0.;
#line 487 "zlals0.f"
		    } else {
#line 488 "zlals0.f"
			d__1 = -poles[i__ + (poles_dim1 << 1)];
#line 488 "zlals0.f"
			rwork[i__] = z__[j] / (dlamc3_(&dsigj, &d__1) - difl[
				i__]) / (dsigj + poles[i__ + poles_dim1]) / 
				difr[i__ + (difr_dim1 << 1)];
#line 491 "zlals0.f"
		    }
#line 492 "zlals0.f"
/* L120: */
#line 492 "zlals0.f"
		}

/*              Since B and BX are complex, the following call to DGEMV */
/*              is performed in two steps (real and imaginary parts). */

/*              CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO, */
/*    $                     BX( J, 1 ), LDBX ) */

#line 500 "zlals0.f"
		i__ = *k + (*nrhs << 1);
#line 501 "zlals0.f"
		i__2 = *nrhs;
#line 501 "zlals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 502 "zlals0.f"
		    i__3 = *k;
#line 502 "zlals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 503 "zlals0.f"
			++i__;
#line 504 "zlals0.f"
			i__4 = jrow + jcol * b_dim1;
#line 504 "zlals0.f"
			rwork[i__] = b[i__4].r;
#line 505 "zlals0.f"
/* L130: */
#line 505 "zlals0.f"
		    }
#line 506 "zlals0.f"
/* L140: */
#line 506 "zlals0.f"
		}
#line 507 "zlals0.f"
		dgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1], &c__1, (
			ftnlen)1);
#line 509 "zlals0.f"
		i__ = *k + (*nrhs << 1);
#line 510 "zlals0.f"
		i__2 = *nrhs;
#line 510 "zlals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 511 "zlals0.f"
		    i__3 = *k;
#line 511 "zlals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 512 "zlals0.f"
			++i__;
#line 513 "zlals0.f"
			rwork[i__] = d_imag(&b[jrow + jcol * b_dim1]);
#line 514 "zlals0.f"
/* L150: */
#line 514 "zlals0.f"
		    }
#line 515 "zlals0.f"
/* L160: */
#line 515 "zlals0.f"
		}
#line 516 "zlals0.f"
		dgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1 + *nrhs], &
			c__1, (ftnlen)1);
#line 518 "zlals0.f"
		i__2 = *nrhs;
#line 518 "zlals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 519 "zlals0.f"
		    i__3 = j + jcol * bx_dim1;
#line 519 "zlals0.f"
		    i__4 = jcol + *k;
#line 519 "zlals0.f"
		    i__5 = jcol + *k + *nrhs;
#line 519 "zlals0.f"
		    z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 519 "zlals0.f"
		    bx[i__3].r = z__1.r, bx[i__3].i = z__1.i;
#line 521 "zlals0.f"
/* L170: */
#line 521 "zlals0.f"
		}
#line 522 "zlals0.f"
/* L180: */
#line 522 "zlals0.f"
	    }
#line 523 "zlals0.f"
	}

/*        Step (2R): if SQRE = 1, apply back the rotation that is */
/*        related to the right null space of the subproblem. */

#line 528 "zlals0.f"
	if (*sqre == 1) {
#line 529 "zlals0.f"
	    zcopy_(nrhs, &b[m + b_dim1], ldb, &bx[m + bx_dim1], ldbx);
#line 530 "zlals0.f"
	    zdrot_(nrhs, &bx[bx_dim1 + 1], ldbx, &bx[m + bx_dim1], ldbx, c__, 
		    s);
#line 531 "zlals0.f"
	}
#line 532 "zlals0.f"
	if (*k < max(m,n)) {
#line 532 "zlals0.f"
	    i__1 = n - *k;
#line 532 "zlals0.f"
	    zlacpy_("A", &i__1, nrhs, &b[*k + 1 + b_dim1], ldb, &bx[*k + 1 + 
		    bx_dim1], ldbx, (ftnlen)1);
#line 532 "zlals0.f"
	}

/*        Step (3R): permute rows of B. */

#line 538 "zlals0.f"
	zcopy_(nrhs, &bx[bx_dim1 + 1], ldbx, &b[nlp1 + b_dim1], ldb);
#line 539 "zlals0.f"
	if (*sqre == 1) {
#line 540 "zlals0.f"
	    zcopy_(nrhs, &bx[m + bx_dim1], ldbx, &b[m + b_dim1], ldb);
#line 541 "zlals0.f"
	}
#line 542 "zlals0.f"
	i__1 = n;
#line 542 "zlals0.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 543 "zlals0.f"
	    zcopy_(nrhs, &bx[i__ + bx_dim1], ldbx, &b[perm[i__] + b_dim1], 
		    ldb);
#line 544 "zlals0.f"
/* L190: */
#line 544 "zlals0.f"
	}

/*        Step (4R): apply back the Givens rotations performed. */

#line 548 "zlals0.f"
	for (i__ = *givptr; i__ >= 1; --i__) {
#line 549 "zlals0.f"
	    d__1 = -givnum[i__ + givnum_dim1];
#line 549 "zlals0.f"
	    zdrot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb, &
		    b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + 
		    (givnum_dim1 << 1)], &d__1);
#line 552 "zlals0.f"
/* L200: */
#line 552 "zlals0.f"
	}
#line 553 "zlals0.f"
    }

#line 555 "zlals0.f"
    return 0;

/*     End of ZLALS0 */

} /* zlals0_ */

