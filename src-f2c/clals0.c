#line 1 "clals0.f"
/* clals0.f -- translated by f2c (version 20100827).
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

#line 1 "clals0.f"
/* Table of constant values */

static doublereal c_b5 = -1.;
static integer c__1 = 1;
static doublereal c_b13 = 1.;
static doublereal c_b15 = 0.;
static integer c__0 = 0;

/* > \brief \b CLALS0 applies back multiplying factors in solving the least squares problem using divide and c
onquer SVD approach. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLALS0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clals0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clals0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clals0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, */
/*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, */
/*                          POLES, DIFL, DIFR, Z, K, C, S, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, */
/*      $                   LDGNUM, NL, NR, NRHS, SQRE */
/*       REAL               C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), PERM( * ) */
/*       REAL               DIFL( * ), DIFR( LDGNUM, * ), */
/*      $                   GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), */
/*      $                   RWORK( * ), Z( * ) */
/*       COMPLEX            B( LDB, * ), BX( LDBX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLALS0 applies back the multiplying factors of either the left or the */
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
/* >          B is COMPLEX array, dimension ( LDB, NRHS ) */
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
/* >          BX is COMPLEX array, dimension ( LDBX, NRHS ) */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int clals0_(integer *icompq, integer *nl, integer *nr, 
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal diflj, difrj, dsigj;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), csrot_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *);
    extern doublereal slamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), clacpy_(char *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal dsigjp;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 316 "clals0.f"
    /* Parameter adjustments */
#line 316 "clals0.f"
    b_dim1 = *ldb;
#line 316 "clals0.f"
    b_offset = 1 + b_dim1;
#line 316 "clals0.f"
    b -= b_offset;
#line 316 "clals0.f"
    bx_dim1 = *ldbx;
#line 316 "clals0.f"
    bx_offset = 1 + bx_dim1;
#line 316 "clals0.f"
    bx -= bx_offset;
#line 316 "clals0.f"
    --perm;
#line 316 "clals0.f"
    givcol_dim1 = *ldgcol;
#line 316 "clals0.f"
    givcol_offset = 1 + givcol_dim1;
#line 316 "clals0.f"
    givcol -= givcol_offset;
#line 316 "clals0.f"
    difr_dim1 = *ldgnum;
#line 316 "clals0.f"
    difr_offset = 1 + difr_dim1;
#line 316 "clals0.f"
    difr -= difr_offset;
#line 316 "clals0.f"
    poles_dim1 = *ldgnum;
#line 316 "clals0.f"
    poles_offset = 1 + poles_dim1;
#line 316 "clals0.f"
    poles -= poles_offset;
#line 316 "clals0.f"
    givnum_dim1 = *ldgnum;
#line 316 "clals0.f"
    givnum_offset = 1 + givnum_dim1;
#line 316 "clals0.f"
    givnum -= givnum_offset;
#line 316 "clals0.f"
    --difl;
#line 316 "clals0.f"
    --z__;
#line 316 "clals0.f"
    --rwork;
#line 316 "clals0.f"

#line 316 "clals0.f"
    /* Function Body */
#line 316 "clals0.f"
    *info = 0;
#line 317 "clals0.f"
    n = *nl + *nr + 1;

#line 319 "clals0.f"
    if (*icompq < 0 || *icompq > 1) {
#line 320 "clals0.f"
	*info = -1;
#line 321 "clals0.f"
    } else if (*nl < 1) {
#line 322 "clals0.f"
	*info = -2;
#line 323 "clals0.f"
    } else if (*nr < 1) {
#line 324 "clals0.f"
	*info = -3;
#line 325 "clals0.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 326 "clals0.f"
	*info = -4;
#line 327 "clals0.f"
    } else if (*nrhs < 1) {
#line 328 "clals0.f"
	*info = -5;
#line 329 "clals0.f"
    } else if (*ldb < n) {
#line 330 "clals0.f"
	*info = -7;
#line 331 "clals0.f"
    } else if (*ldbx < n) {
#line 332 "clals0.f"
	*info = -9;
#line 333 "clals0.f"
    } else if (*givptr < 0) {
#line 334 "clals0.f"
	*info = -11;
#line 335 "clals0.f"
    } else if (*ldgcol < n) {
#line 336 "clals0.f"
	*info = -13;
#line 337 "clals0.f"
    } else if (*ldgnum < n) {
#line 338 "clals0.f"
	*info = -15;
#line 339 "clals0.f"
    } else if (*k < 1) {
#line 340 "clals0.f"
	*info = -20;
#line 341 "clals0.f"
    }
#line 342 "clals0.f"
    if (*info != 0) {
#line 343 "clals0.f"
	i__1 = -(*info);
#line 343 "clals0.f"
	xerbla_("CLALS0", &i__1, (ftnlen)6);
#line 344 "clals0.f"
	return 0;
#line 345 "clals0.f"
    }

#line 347 "clals0.f"
    m = n + *sqre;
#line 348 "clals0.f"
    nlp1 = *nl + 1;

#line 350 "clals0.f"
    if (*icompq == 0) {

/*        Apply back orthogonal transformations from the left. */

/*        Step (1L): apply back the Givens rotations performed. */

#line 356 "clals0.f"
	i__1 = *givptr;
#line 356 "clals0.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 357 "clals0.f"
	    csrot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb, &
		    b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + 
		    (givnum_dim1 << 1)], &givnum[i__ + givnum_dim1]);
#line 360 "clals0.f"
/* L10: */
#line 360 "clals0.f"
	}

/*        Step (2L): permute rows of B. */

#line 364 "clals0.f"
	ccopy_(nrhs, &b[nlp1 + b_dim1], ldb, &bx[bx_dim1 + 1], ldbx);
#line 365 "clals0.f"
	i__1 = n;
#line 365 "clals0.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 366 "clals0.f"
	    ccopy_(nrhs, &b[perm[i__] + b_dim1], ldb, &bx[i__ + bx_dim1], 
		    ldbx);
#line 367 "clals0.f"
/* L20: */
#line 367 "clals0.f"
	}

/*        Step (3L): apply the inverse of the left singular vector */
/*        matrix to BX. */

#line 372 "clals0.f"
	if (*k == 1) {
#line 373 "clals0.f"
	    ccopy_(nrhs, &bx[bx_offset], ldbx, &b[b_offset], ldb);
#line 374 "clals0.f"
	    if (z__[1] < 0.) {
#line 375 "clals0.f"
		csscal_(nrhs, &c_b5, &b[b_offset], ldb);
#line 376 "clals0.f"
	    }
#line 377 "clals0.f"
	} else {
#line 378 "clals0.f"
	    i__1 = *k;
#line 378 "clals0.f"
	    for (j = 1; j <= i__1; ++j) {
#line 379 "clals0.f"
		diflj = difl[j];
#line 380 "clals0.f"
		dj = poles[j + poles_dim1];
#line 381 "clals0.f"
		dsigj = -poles[j + (poles_dim1 << 1)];
#line 382 "clals0.f"
		if (j < *k) {
#line 383 "clals0.f"
		    difrj = -difr[j + difr_dim1];
#line 384 "clals0.f"
		    dsigjp = -poles[j + 1 + (poles_dim1 << 1)];
#line 385 "clals0.f"
		}
#line 386 "clals0.f"
		if (z__[j] == 0. || poles[j + (poles_dim1 << 1)] == 0.) {
#line 388 "clals0.f"
		    rwork[j] = 0.;
#line 389 "clals0.f"
		} else {
#line 390 "clals0.f"
		    rwork[j] = -poles[j + (poles_dim1 << 1)] * z__[j] / diflj 
			    / (poles[j + (poles_dim1 << 1)] + dj);
#line 392 "clals0.f"
		}
#line 393 "clals0.f"
		i__2 = j - 1;
#line 393 "clals0.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 394 "clals0.f"
		    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 
			    0.) {
#line 396 "clals0.f"
			rwork[i__] = 0.;
#line 397 "clals0.f"
		    } else {
#line 398 "clals0.f"
			rwork[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__]
				 / (slamc3_(&poles[i__ + (poles_dim1 << 1)], &
				dsigj) - diflj) / (poles[i__ + (poles_dim1 << 
				1)] + dj);
#line 401 "clals0.f"
		    }
#line 402 "clals0.f"
/* L30: */
#line 402 "clals0.f"
		}
#line 403 "clals0.f"
		i__2 = *k;
#line 403 "clals0.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 404 "clals0.f"
		    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 
			    0.) {
#line 406 "clals0.f"
			rwork[i__] = 0.;
#line 407 "clals0.f"
		    } else {
#line 408 "clals0.f"
			rwork[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__]
				 / (slamc3_(&poles[i__ + (poles_dim1 << 1)], &
				dsigjp) + difrj) / (poles[i__ + (poles_dim1 <<
				 1)] + dj);
#line 411 "clals0.f"
		    }
#line 412 "clals0.f"
/* L40: */
#line 412 "clals0.f"
		}
#line 413 "clals0.f"
		rwork[1] = -1.;
#line 414 "clals0.f"
		temp = snrm2_(k, &rwork[1], &c__1);

/*              Since B and BX are complex, the following call to SGEMV */
/*              is performed in two steps (real and imaginary parts). */

/*              CALL SGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO, */
/*    $                     B( J, 1 ), LDB ) */

#line 422 "clals0.f"
		i__ = *k + (*nrhs << 1);
#line 423 "clals0.f"
		i__2 = *nrhs;
#line 423 "clals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 424 "clals0.f"
		    i__3 = *k;
#line 424 "clals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 425 "clals0.f"
			++i__;
#line 426 "clals0.f"
			i__4 = jrow + jcol * bx_dim1;
#line 426 "clals0.f"
			rwork[i__] = bx[i__4].r;
#line 427 "clals0.f"
/* L50: */
#line 427 "clals0.f"
		    }
#line 428 "clals0.f"
/* L60: */
#line 428 "clals0.f"
		}
#line 429 "clals0.f"
		sgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1], &c__1, (
			ftnlen)1);
#line 431 "clals0.f"
		i__ = *k + (*nrhs << 1);
#line 432 "clals0.f"
		i__2 = *nrhs;
#line 432 "clals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 433 "clals0.f"
		    i__3 = *k;
#line 433 "clals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 434 "clals0.f"
			++i__;
#line 435 "clals0.f"
			rwork[i__] = d_imag(&bx[jrow + jcol * bx_dim1]);
#line 436 "clals0.f"
/* L70: */
#line 436 "clals0.f"
		    }
#line 437 "clals0.f"
/* L80: */
#line 437 "clals0.f"
		}
#line 438 "clals0.f"
		sgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1 + *nrhs], &
			c__1, (ftnlen)1);
#line 440 "clals0.f"
		i__2 = *nrhs;
#line 440 "clals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 441 "clals0.f"
		    i__3 = j + jcol * b_dim1;
#line 441 "clals0.f"
		    i__4 = jcol + *k;
#line 441 "clals0.f"
		    i__5 = jcol + *k + *nrhs;
#line 441 "clals0.f"
		    z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 441 "clals0.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 443 "clals0.f"
/* L90: */
#line 443 "clals0.f"
		}
#line 444 "clals0.f"
		clascl_("G", &c__0, &c__0, &temp, &c_b13, &c__1, nrhs, &b[j + 
			b_dim1], ldb, info, (ftnlen)1);
#line 446 "clals0.f"
/* L100: */
#line 446 "clals0.f"
	    }
#line 447 "clals0.f"
	}

/*        Move the deflated rows of BX to B also. */

#line 451 "clals0.f"
	if (*k < max(m,n)) {
#line 451 "clals0.f"
	    i__1 = n - *k;
#line 451 "clals0.f"
	    clacpy_("A", &i__1, nrhs, &bx[*k + 1 + bx_dim1], ldbx, &b[*k + 1 
		    + b_dim1], ldb, (ftnlen)1);
#line 451 "clals0.f"
	}
#line 454 "clals0.f"
    } else {

/*        Apply back the right orthogonal transformations. */

/*        Step (1R): apply back the new right singular vector matrix */
/*        to B. */

#line 461 "clals0.f"
	if (*k == 1) {
#line 462 "clals0.f"
	    ccopy_(nrhs, &b[b_offset], ldb, &bx[bx_offset], ldbx);
#line 463 "clals0.f"
	} else {
#line 464 "clals0.f"
	    i__1 = *k;
#line 464 "clals0.f"
	    for (j = 1; j <= i__1; ++j) {
#line 465 "clals0.f"
		dsigj = poles[j + (poles_dim1 << 1)];
#line 466 "clals0.f"
		if (z__[j] == 0.) {
#line 467 "clals0.f"
		    rwork[j] = 0.;
#line 468 "clals0.f"
		} else {
#line 469 "clals0.f"
		    rwork[j] = -z__[j] / difl[j] / (dsigj + poles[j + 
			    poles_dim1]) / difr[j + (difr_dim1 << 1)];
#line 471 "clals0.f"
		}
#line 472 "clals0.f"
		i__2 = j - 1;
#line 472 "clals0.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 473 "clals0.f"
		    if (z__[j] == 0.) {
#line 474 "clals0.f"
			rwork[i__] = 0.;
#line 475 "clals0.f"
		    } else {
#line 476 "clals0.f"
			d__1 = -poles[i__ + 1 + (poles_dim1 << 1)];
#line 476 "clals0.f"
			rwork[i__] = z__[j] / (slamc3_(&dsigj, &d__1) - difr[
				i__ + difr_dim1]) / (dsigj + poles[i__ + 
				poles_dim1]) / difr[i__ + (difr_dim1 << 1)];
#line 479 "clals0.f"
		    }
#line 480 "clals0.f"
/* L110: */
#line 480 "clals0.f"
		}
#line 481 "clals0.f"
		i__2 = *k;
#line 481 "clals0.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 482 "clals0.f"
		    if (z__[j] == 0.) {
#line 483 "clals0.f"
			rwork[i__] = 0.;
#line 484 "clals0.f"
		    } else {
#line 485 "clals0.f"
			d__1 = -poles[i__ + (poles_dim1 << 1)];
#line 485 "clals0.f"
			rwork[i__] = z__[j] / (slamc3_(&dsigj, &d__1) - difl[
				i__]) / (dsigj + poles[i__ + poles_dim1]) / 
				difr[i__ + (difr_dim1 << 1)];
#line 488 "clals0.f"
		    }
#line 489 "clals0.f"
/* L120: */
#line 489 "clals0.f"
		}

/*              Since B and BX are complex, the following call to SGEMV */
/*              is performed in two steps (real and imaginary parts). */

/*              CALL SGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO, */
/*    $                     BX( J, 1 ), LDBX ) */

#line 497 "clals0.f"
		i__ = *k + (*nrhs << 1);
#line 498 "clals0.f"
		i__2 = *nrhs;
#line 498 "clals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 499 "clals0.f"
		    i__3 = *k;
#line 499 "clals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 500 "clals0.f"
			++i__;
#line 501 "clals0.f"
			i__4 = jrow + jcol * b_dim1;
#line 501 "clals0.f"
			rwork[i__] = b[i__4].r;
#line 502 "clals0.f"
/* L130: */
#line 502 "clals0.f"
		    }
#line 503 "clals0.f"
/* L140: */
#line 503 "clals0.f"
		}
#line 504 "clals0.f"
		sgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1], &c__1, (
			ftnlen)1);
#line 506 "clals0.f"
		i__ = *k + (*nrhs << 1);
#line 507 "clals0.f"
		i__2 = *nrhs;
#line 507 "clals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 508 "clals0.f"
		    i__3 = *k;
#line 508 "clals0.f"
		    for (jrow = 1; jrow <= i__3; ++jrow) {
#line 509 "clals0.f"
			++i__;
#line 510 "clals0.f"
			rwork[i__] = d_imag(&b[jrow + jcol * b_dim1]);
#line 511 "clals0.f"
/* L150: */
#line 511 "clals0.f"
		    }
#line 512 "clals0.f"
/* L160: */
#line 512 "clals0.f"
		}
#line 513 "clals0.f"
		sgemv_("T", k, nrhs, &c_b13, &rwork[*k + 1 + (*nrhs << 1)], k,
			 &rwork[1], &c__1, &c_b15, &rwork[*k + 1 + *nrhs], &
			c__1, (ftnlen)1);
#line 515 "clals0.f"
		i__2 = *nrhs;
#line 515 "clals0.f"
		for (jcol = 1; jcol <= i__2; ++jcol) {
#line 516 "clals0.f"
		    i__3 = j + jcol * bx_dim1;
#line 516 "clals0.f"
		    i__4 = jcol + *k;
#line 516 "clals0.f"
		    i__5 = jcol + *k + *nrhs;
#line 516 "clals0.f"
		    z__1.r = rwork[i__4], z__1.i = rwork[i__5];
#line 516 "clals0.f"
		    bx[i__3].r = z__1.r, bx[i__3].i = z__1.i;
#line 518 "clals0.f"
/* L170: */
#line 518 "clals0.f"
		}
#line 519 "clals0.f"
/* L180: */
#line 519 "clals0.f"
	    }
#line 520 "clals0.f"
	}

/*        Step (2R): if SQRE = 1, apply back the rotation that is */
/*        related to the right null space of the subproblem. */

#line 525 "clals0.f"
	if (*sqre == 1) {
#line 526 "clals0.f"
	    ccopy_(nrhs, &b[m + b_dim1], ldb, &bx[m + bx_dim1], ldbx);
#line 527 "clals0.f"
	    csrot_(nrhs, &bx[bx_dim1 + 1], ldbx, &bx[m + bx_dim1], ldbx, c__, 
		    s);
#line 528 "clals0.f"
	}
#line 529 "clals0.f"
	if (*k < max(m,n)) {
#line 529 "clals0.f"
	    i__1 = n - *k;
#line 529 "clals0.f"
	    clacpy_("A", &i__1, nrhs, &b[*k + 1 + b_dim1], ldb, &bx[*k + 1 + 
		    bx_dim1], ldbx, (ftnlen)1);
#line 529 "clals0.f"
	}

/*        Step (3R): permute rows of B. */

#line 535 "clals0.f"
	ccopy_(nrhs, &bx[bx_dim1 + 1], ldbx, &b[nlp1 + b_dim1], ldb);
#line 536 "clals0.f"
	if (*sqre == 1) {
#line 537 "clals0.f"
	    ccopy_(nrhs, &bx[m + bx_dim1], ldbx, &b[m + b_dim1], ldb);
#line 538 "clals0.f"
	}
#line 539 "clals0.f"
	i__1 = n;
#line 539 "clals0.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 540 "clals0.f"
	    ccopy_(nrhs, &bx[i__ + bx_dim1], ldbx, &b[perm[i__] + b_dim1], 
		    ldb);
#line 541 "clals0.f"
/* L190: */
#line 541 "clals0.f"
	}

/*        Step (4R): apply back the Givens rotations performed. */

#line 545 "clals0.f"
	for (i__ = *givptr; i__ >= 1; --i__) {
#line 546 "clals0.f"
	    d__1 = -givnum[i__ + givnum_dim1];
#line 546 "clals0.f"
	    csrot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb, &
		    b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + 
		    (givnum_dim1 << 1)], &d__1);
#line 549 "clals0.f"
/* L200: */
#line 549 "clals0.f"
	}
#line 550 "clals0.f"
    }

#line 552 "clals0.f"
    return 0;

/*     End of CLALS0 */

} /* clals0_ */

