#line 1 "slals0.f"
/* slals0.f -- translated by f2c (version 20100827).
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

#line 1 "slals0.f"
/* Table of constant values */

static doublereal c_b5 = -1.;
static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b13 = 0.;
static integer c__0 = 0;

/* > \brief \b SLALS0 applies back multiplying factors in solving the least squares problem using divide and c
onquer SVD approach. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLALS0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slals0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slals0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slals0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, */
/*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, */
/*                          POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, */
/*      $                   LDGNUM, NL, NR, NRHS, SQRE */
/*       REAL               C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), PERM( * ) */
/*       REAL               B( LDB, * ), BX( LDBX, * ), DIFL( * ), */
/*      $                   DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), */
/*      $                   POLES( LDGNUM, * ), WORK( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALS0 applies back the multiplying factors of either the left or the */
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
/* >          B is REAL array, dimension ( LDB, NRHS ) */
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
/* >          BX is REAL array, dimension ( LDBX, NRHS ) */
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
/* >          GIVNUM is REAL array, dimension ( LDGNUM, 2 ) */
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
/* >          POLES is REAL array, dimension ( LDGNUM, 2 ) */
/* >         On entry, POLES(1:K, 1) contains the new singular */
/* >         values obtained from solving the secular equation, and */
/* >         POLES(1:K, 2) is an array containing the poles in the secular */
/* >         equation. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFL */
/* > \verbatim */
/* >          DIFL is REAL array, dimension ( K ). */
/* >         On entry, DIFL(I) is the distance between I-th updated */
/* >         (undeflated) singular value and the I-th (undeflated) old */
/* >         singular value. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFR */
/* > \verbatim */
/* >          DIFR is REAL array, dimension ( LDGNUM, 2 ). */
/* >         On entry, DIFR(I, 1) contains the distances between I-th */
/* >         updated (undeflated) singular value and the I+1-th */
/* >         (undeflated) old singular value. And DIFR(I, 2) is the */
/* >         normalizing factor for the I-th right singular vector. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension ( K ) */
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
/* >          C is REAL */
/* >         C contains garbage if SQRE =0 and the C-value of a Givens */
/* >         rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL */
/* >         S contains garbage if SQRE =0 and the S-value of a Givens */
/* >         rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension ( K ) */
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

/* > \date November 2015 */

/* > \ingroup realOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int slals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal 
	*bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *
	poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *
	k, doublereal *c__, doublereal *s, doublereal *work, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, b_dim1, b_offset, bx_dim1, bx_offset, 
	    difr_dim1, difr_offset, givnum_dim1, givnum_offset, poles_dim1, 
	    poles_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, m, n;
    static doublereal dj;
    static integer nlp1;
    static doublereal temp;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal diflj, difrj, dsigj;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), scopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern doublereal slamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal dsigjp;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), slacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 313 "slals0.f"
    /* Parameter adjustments */
#line 313 "slals0.f"
    b_dim1 = *ldb;
#line 313 "slals0.f"
    b_offset = 1 + b_dim1;
#line 313 "slals0.f"
    b -= b_offset;
#line 313 "slals0.f"
    bx_dim1 = *ldbx;
#line 313 "slals0.f"
    bx_offset = 1 + bx_dim1;
#line 313 "slals0.f"
    bx -= bx_offset;
#line 313 "slals0.f"
    --perm;
#line 313 "slals0.f"
    givcol_dim1 = *ldgcol;
#line 313 "slals0.f"
    givcol_offset = 1 + givcol_dim1;
#line 313 "slals0.f"
    givcol -= givcol_offset;
#line 313 "slals0.f"
    difr_dim1 = *ldgnum;
#line 313 "slals0.f"
    difr_offset = 1 + difr_dim1;
#line 313 "slals0.f"
    difr -= difr_offset;
#line 313 "slals0.f"
    poles_dim1 = *ldgnum;
#line 313 "slals0.f"
    poles_offset = 1 + poles_dim1;
#line 313 "slals0.f"
    poles -= poles_offset;
#line 313 "slals0.f"
    givnum_dim1 = *ldgnum;
#line 313 "slals0.f"
    givnum_offset = 1 + givnum_dim1;
#line 313 "slals0.f"
    givnum -= givnum_offset;
#line 313 "slals0.f"
    --difl;
#line 313 "slals0.f"
    --z__;
#line 313 "slals0.f"
    --work;
#line 313 "slals0.f"

#line 313 "slals0.f"
    /* Function Body */
#line 313 "slals0.f"
    *info = 0;
#line 314 "slals0.f"
    n = *nl + *nr + 1;

#line 316 "slals0.f"
    if (*icompq < 0 || *icompq > 1) {
#line 317 "slals0.f"
	*info = -1;
#line 318 "slals0.f"
    } else if (*nl < 1) {
#line 319 "slals0.f"
	*info = -2;
#line 320 "slals0.f"
    } else if (*nr < 1) {
#line 321 "slals0.f"
	*info = -3;
#line 322 "slals0.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 323 "slals0.f"
	*info = -4;
#line 324 "slals0.f"
    } else if (*nrhs < 1) {
#line 325 "slals0.f"
	*info = -5;
#line 326 "slals0.f"
    } else if (*ldb < n) {
#line 327 "slals0.f"
	*info = -7;
#line 328 "slals0.f"
    } else if (*ldbx < n) {
#line 329 "slals0.f"
	*info = -9;
#line 330 "slals0.f"
    } else if (*givptr < 0) {
#line 331 "slals0.f"
	*info = -11;
#line 332 "slals0.f"
    } else if (*ldgcol < n) {
#line 333 "slals0.f"
	*info = -13;
#line 334 "slals0.f"
    } else if (*ldgnum < n) {
#line 335 "slals0.f"
	*info = -15;
#line 336 "slals0.f"
    } else if (*k < 1) {
#line 337 "slals0.f"
	*info = -20;
#line 338 "slals0.f"
    }
#line 339 "slals0.f"
    if (*info != 0) {
#line 340 "slals0.f"
	i__1 = -(*info);
#line 340 "slals0.f"
	xerbla_("SLALS0", &i__1, (ftnlen)6);
#line 341 "slals0.f"
	return 0;
#line 342 "slals0.f"
    }

#line 344 "slals0.f"
    m = n + *sqre;
#line 345 "slals0.f"
    nlp1 = *nl + 1;

#line 347 "slals0.f"
    if (*icompq == 0) {

/*        Apply back orthogonal transformations from the left. */

/*        Step (1L): apply back the Givens rotations performed. */

#line 353 "slals0.f"
	i__1 = *givptr;
#line 353 "slals0.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 354 "slals0.f"
	    srot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb, &
		    b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + 
		    (givnum_dim1 << 1)], &givnum[i__ + givnum_dim1]);
#line 357 "slals0.f"
/* L10: */
#line 357 "slals0.f"
	}

/*        Step (2L): permute rows of B. */

#line 361 "slals0.f"
	scopy_(nrhs, &b[nlp1 + b_dim1], ldb, &bx[bx_dim1 + 1], ldbx);
#line 362 "slals0.f"
	i__1 = n;
#line 362 "slals0.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 363 "slals0.f"
	    scopy_(nrhs, &b[perm[i__] + b_dim1], ldb, &bx[i__ + bx_dim1], 
		    ldbx);
#line 364 "slals0.f"
/* L20: */
#line 364 "slals0.f"
	}

/*        Step (3L): apply the inverse of the left singular vector */
/*        matrix to BX. */

#line 369 "slals0.f"
	if (*k == 1) {
#line 370 "slals0.f"
	    scopy_(nrhs, &bx[bx_offset], ldbx, &b[b_offset], ldb);
#line 371 "slals0.f"
	    if (z__[1] < 0.) {
#line 372 "slals0.f"
		sscal_(nrhs, &c_b5, &b[b_offset], ldb);
#line 373 "slals0.f"
	    }
#line 374 "slals0.f"
	} else {
#line 375 "slals0.f"
	    i__1 = *k;
#line 375 "slals0.f"
	    for (j = 1; j <= i__1; ++j) {
#line 376 "slals0.f"
		diflj = difl[j];
#line 377 "slals0.f"
		dj = poles[j + poles_dim1];
#line 378 "slals0.f"
		dsigj = -poles[j + (poles_dim1 << 1)];
#line 379 "slals0.f"
		if (j < *k) {
#line 380 "slals0.f"
		    difrj = -difr[j + difr_dim1];
#line 381 "slals0.f"
		    dsigjp = -poles[j + 1 + (poles_dim1 << 1)];
#line 382 "slals0.f"
		}
#line 383 "slals0.f"
		if (z__[j] == 0. || poles[j + (poles_dim1 << 1)] == 0.) {
#line 385 "slals0.f"
		    work[j] = 0.;
#line 386 "slals0.f"
		} else {
#line 387 "slals0.f"
		    work[j] = -poles[j + (poles_dim1 << 1)] * z__[j] / diflj /
			     (poles[j + (poles_dim1 << 1)] + dj);
#line 389 "slals0.f"
		}
#line 390 "slals0.f"
		i__2 = j - 1;
#line 390 "slals0.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 391 "slals0.f"
		    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 
			    0.) {
#line 393 "slals0.f"
			work[i__] = 0.;
#line 394 "slals0.f"
		    } else {
#line 395 "slals0.f"
			work[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__] 
				/ (slamc3_(&poles[i__ + (poles_dim1 << 1)], &
				dsigj) - diflj) / (poles[i__ + (poles_dim1 << 
				1)] + dj);
#line 398 "slals0.f"
		    }
#line 399 "slals0.f"
/* L30: */
#line 399 "slals0.f"
		}
#line 400 "slals0.f"
		i__2 = *k;
#line 400 "slals0.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 401 "slals0.f"
		    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 
			    0.) {
#line 403 "slals0.f"
			work[i__] = 0.;
#line 404 "slals0.f"
		    } else {
#line 405 "slals0.f"
			work[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__] 
				/ (slamc3_(&poles[i__ + (poles_dim1 << 1)], &
				dsigjp) + difrj) / (poles[i__ + (poles_dim1 <<
				 1)] + dj);
#line 408 "slals0.f"
		    }
#line 409 "slals0.f"
/* L40: */
#line 409 "slals0.f"
		}
#line 410 "slals0.f"
		work[1] = -1.;
#line 411 "slals0.f"
		temp = snrm2_(k, &work[1], &c__1);
#line 412 "slals0.f"
		sgemv_("T", k, nrhs, &c_b11, &bx[bx_offset], ldbx, &work[1], &
			c__1, &c_b13, &b[j + b_dim1], ldb, (ftnlen)1);
#line 414 "slals0.f"
		slascl_("G", &c__0, &c__0, &temp, &c_b11, &c__1, nrhs, &b[j + 
			b_dim1], ldb, info, (ftnlen)1);
#line 416 "slals0.f"
/* L50: */
#line 416 "slals0.f"
	    }
#line 417 "slals0.f"
	}

/*        Move the deflated rows of BX to B also. */

#line 421 "slals0.f"
	if (*k < max(m,n)) {
#line 421 "slals0.f"
	    i__1 = n - *k;
#line 421 "slals0.f"
	    slacpy_("A", &i__1, nrhs, &bx[*k + 1 + bx_dim1], ldbx, &b[*k + 1 
		    + b_dim1], ldb, (ftnlen)1);
#line 421 "slals0.f"
	}
#line 424 "slals0.f"
    } else {

/*        Apply back the right orthogonal transformations. */

/*        Step (1R): apply back the new right singular vector matrix */
/*        to B. */

#line 431 "slals0.f"
	if (*k == 1) {
#line 432 "slals0.f"
	    scopy_(nrhs, &b[b_offset], ldb, &bx[bx_offset], ldbx);
#line 433 "slals0.f"
	} else {
#line 434 "slals0.f"
	    i__1 = *k;
#line 434 "slals0.f"
	    for (j = 1; j <= i__1; ++j) {
#line 435 "slals0.f"
		dsigj = poles[j + (poles_dim1 << 1)];
#line 436 "slals0.f"
		if (z__[j] == 0.) {
#line 437 "slals0.f"
		    work[j] = 0.;
#line 438 "slals0.f"
		} else {
#line 439 "slals0.f"
		    work[j] = -z__[j] / difl[j] / (dsigj + poles[j + 
			    poles_dim1]) / difr[j + (difr_dim1 << 1)];
#line 441 "slals0.f"
		}
#line 442 "slals0.f"
		i__2 = j - 1;
#line 442 "slals0.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 443 "slals0.f"
		    if (z__[j] == 0.) {
#line 444 "slals0.f"
			work[i__] = 0.;
#line 445 "slals0.f"
		    } else {
#line 446 "slals0.f"
			d__1 = -poles[i__ + 1 + (poles_dim1 << 1)];
#line 446 "slals0.f"
			work[i__] = z__[j] / (slamc3_(&dsigj, &d__1) - difr[
				i__ + difr_dim1]) / (dsigj + poles[i__ + 
				poles_dim1]) / difr[i__ + (difr_dim1 << 1)];
#line 449 "slals0.f"
		    }
#line 450 "slals0.f"
/* L60: */
#line 450 "slals0.f"
		}
#line 451 "slals0.f"
		i__2 = *k;
#line 451 "slals0.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 452 "slals0.f"
		    if (z__[j] == 0.) {
#line 453 "slals0.f"
			work[i__] = 0.;
#line 454 "slals0.f"
		    } else {
#line 455 "slals0.f"
			d__1 = -poles[i__ + (poles_dim1 << 1)];
#line 455 "slals0.f"
			work[i__] = z__[j] / (slamc3_(&dsigj, &d__1) - difl[
				i__]) / (dsigj + poles[i__ + poles_dim1]) / 
				difr[i__ + (difr_dim1 << 1)];
#line 458 "slals0.f"
		    }
#line 459 "slals0.f"
/* L70: */
#line 459 "slals0.f"
		}
#line 460 "slals0.f"
		sgemv_("T", k, nrhs, &c_b11, &b[b_offset], ldb, &work[1], &
			c__1, &c_b13, &bx[j + bx_dim1], ldbx, (ftnlen)1);
#line 462 "slals0.f"
/* L80: */
#line 462 "slals0.f"
	    }
#line 463 "slals0.f"
	}

/*        Step (2R): if SQRE = 1, apply back the rotation that is */
/*        related to the right null space of the subproblem. */

#line 468 "slals0.f"
	if (*sqre == 1) {
#line 469 "slals0.f"
	    scopy_(nrhs, &b[m + b_dim1], ldb, &bx[m + bx_dim1], ldbx);
#line 470 "slals0.f"
	    srot_(nrhs, &bx[bx_dim1 + 1], ldbx, &bx[m + bx_dim1], ldbx, c__, 
		    s);
#line 471 "slals0.f"
	}
#line 472 "slals0.f"
	if (*k < max(m,n)) {
#line 472 "slals0.f"
	    i__1 = n - *k;
#line 472 "slals0.f"
	    slacpy_("A", &i__1, nrhs, &b[*k + 1 + b_dim1], ldb, &bx[*k + 1 + 
		    bx_dim1], ldbx, (ftnlen)1);
#line 472 "slals0.f"
	}

/*        Step (3R): permute rows of B. */

#line 478 "slals0.f"
	scopy_(nrhs, &bx[bx_dim1 + 1], ldbx, &b[nlp1 + b_dim1], ldb);
#line 479 "slals0.f"
	if (*sqre == 1) {
#line 480 "slals0.f"
	    scopy_(nrhs, &bx[m + bx_dim1], ldbx, &b[m + b_dim1], ldb);
#line 481 "slals0.f"
	}
#line 482 "slals0.f"
	i__1 = n;
#line 482 "slals0.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 483 "slals0.f"
	    scopy_(nrhs, &bx[i__ + bx_dim1], ldbx, &b[perm[i__] + b_dim1], 
		    ldb);
#line 484 "slals0.f"
/* L90: */
#line 484 "slals0.f"
	}

/*        Step (4R): apply back the Givens rotations performed. */

#line 488 "slals0.f"
	for (i__ = *givptr; i__ >= 1; --i__) {
#line 489 "slals0.f"
	    d__1 = -givnum[i__ + givnum_dim1];
#line 489 "slals0.f"
	    srot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb, &
		    b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + 
		    (givnum_dim1 << 1)], &d__1);
#line 492 "slals0.f"
/* L100: */
#line 492 "slals0.f"
	}
#line 493 "slals0.f"
    }

#line 495 "slals0.f"
    return 0;

/*     End of SLALS0 */

} /* slals0_ */

