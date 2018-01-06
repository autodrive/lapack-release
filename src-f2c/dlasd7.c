#line 1 "dlasd7.f"
/* dlasd7.f -- translated by f2c (version 20100827).
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

#line 1 "dlasd7.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLASD7 merges the two sets of singular values together into a single sorted set. Then it tries 
to deflate the size of the problem. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASD7 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd7.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd7.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd7.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL, */
/*                          VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ, */
/*                          PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, */
/*                          C, S, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, */
/*      $                   NR, SQRE */
/*       DOUBLE PRECISION   ALPHA, BETA, C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ), */
/*      $                   IDXQ( * ), PERM( * ) */
/*       DOUBLE PRECISION   D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ), */
/*      $                   VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ), */
/*      $                   ZW( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASD7 merges the two sets of singular values together into a single */
/* > sorted set. Then it tries to deflate the size of the problem. There */
/* > are two ways in which deflation can occur:  when two or more singular */
/* > values are close together or if there is a tiny entry in the Z */
/* > vector. For each such occurrence the order of the related */
/* > secular equation problem is reduced by one. */
/* > */
/* > DLASD7 is called from DLASD6. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >          Specifies whether singular vectors are to be computed */
/* >          in compact form, as follows: */
/* >          = 0: Compute singular values only. */
/* >          = 1: Compute singular vectors of upper */
/* >               bidiagonal matrix in compact form. */
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
/* >         The bidiagonal matrix has */
/* >         N = NL + NR + 1 rows and */
/* >         M = N + SQRE >= N columns. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         Contains the dimension of the non-deflated matrix, this is */
/* >         the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension ( N ) */
/* >         On entry D contains the singular values of the two submatrices */
/* >         to be combined. On exit D contains the trailing (N-K) updated */
/* >         singular values (those which were deflated) sorted into */
/* >         increasing order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( M ) */
/* >         On exit Z contains the updating row vector in the secular */
/* >         equation. */
/* > \endverbatim */
/* > */
/* > \param[out] ZW */
/* > \verbatim */
/* >          ZW is DOUBLE PRECISION array, dimension ( M ) */
/* >         Workspace for Z. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VF */
/* > \verbatim */
/* >          VF is DOUBLE PRECISION array, dimension ( M ) */
/* >         On entry, VF(1:NL+1) contains the first components of all */
/* >         right singular vectors of the upper block; and VF(NL+2:M) */
/* >         contains the first components of all right singular vectors */
/* >         of the lower block. On exit, VF contains the first components */
/* >         of all right singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] VFW */
/* > \verbatim */
/* >          VFW is DOUBLE PRECISION array, dimension ( M ) */
/* >         Workspace for VF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION array, dimension ( M ) */
/* >         On entry, VL(1:NL+1) contains the  last components of all */
/* >         right singular vectors of the upper block; and VL(NL+2:M) */
/* >         contains the last components of all right singular vectors */
/* >         of the lower block. On exit, VL contains the last components */
/* >         of all right singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] VLW */
/* > \verbatim */
/* >          VLW is DOUBLE PRECISION array, dimension ( M ) */
/* >         Workspace for VL. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION */
/* >         Contains the diagonal element associated with the added row. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION */
/* >         Contains the off-diagonal element associated with the added */
/* >         row. */
/* > \endverbatim */
/* > */
/* > \param[out] DSIGMA */
/* > \verbatim */
/* >          DSIGMA is DOUBLE PRECISION array, dimension ( N ) */
/* >         Contains a copy of the diagonal elements (K-1 singular values */
/* >         and one zero) in the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] IDX */
/* > \verbatim */
/* >          IDX is INTEGER array, dimension ( N ) */
/* >         This will contain the permutation used to sort the contents of */
/* >         D into ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXP */
/* > \verbatim */
/* >          IDXP is INTEGER array, dimension ( N ) */
/* >         This will contain the permutation used to place deflated */
/* >         values of D at the end of the array. On output IDXP(2:K) */
/* >         points to the nondeflated D-values and IDXP(K+1:N) */
/* >         points to the deflated singular values. */
/* > \endverbatim */
/* > */
/* > \param[in] IDXQ */
/* > \verbatim */
/* >          IDXQ is INTEGER array, dimension ( N ) */
/* >         This contains the permutation which separately sorts the two */
/* >         sub-problems in D into ascending order.  Note that entries in */
/* >         the first half of this permutation must first be moved one */
/* >         position backward; and entries in the second half */
/* >         must first have NL+1 added to their values. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension ( N ) */
/* >         The permutations (from deflation and sorting) to be applied */
/* >         to each singular block. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER */
/* >         The number of Givens rotations which took place in this */
/* >         subproblem. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 ) */
/* >         Each pair of numbers indicates a pair of columns to take place */
/* >         in a Givens rotation. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* >          LDGCOL is INTEGER */
/* >         The leading dimension of GIVCOL, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ) */
/* >         Each number indicates the C or S value to be used in the */
/* >         corresponding Givens rotation. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGNUM */
/* > \verbatim */
/* >          LDGNUM is INTEGER */
/* >         The leading dimension of GIVNUM, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >         C contains garbage if SQRE =0 and the C-value of a Givens */
/* >         rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* >         S contains garbage if SQRE =0 and the S-value of a Givens */
/* >         rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >         = 0:  successful exit. */
/* >         < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasd7_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *k, doublereal *d__, doublereal *z__, 
	doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl, 
	doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *
	dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, 
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	 integer *ldgnum, doublereal *c__, doublereal *s, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, m, n, k2;
    static doublereal z1;
    static integer jp;
    static doublereal eps, tau, tol;
    static integer nlp1, nlp2, idxi, idxj;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer idxjp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jprev;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal hlftol;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 328 "dlasd7.f"
    /* Parameter adjustments */
#line 328 "dlasd7.f"
    --d__;
#line 328 "dlasd7.f"
    --z__;
#line 328 "dlasd7.f"
    --zw;
#line 328 "dlasd7.f"
    --vf;
#line 328 "dlasd7.f"
    --vfw;
#line 328 "dlasd7.f"
    --vl;
#line 328 "dlasd7.f"
    --vlw;
#line 328 "dlasd7.f"
    --dsigma;
#line 328 "dlasd7.f"
    --idx;
#line 328 "dlasd7.f"
    --idxp;
#line 328 "dlasd7.f"
    --idxq;
#line 328 "dlasd7.f"
    --perm;
#line 328 "dlasd7.f"
    givcol_dim1 = *ldgcol;
#line 328 "dlasd7.f"
    givcol_offset = 1 + givcol_dim1;
#line 328 "dlasd7.f"
    givcol -= givcol_offset;
#line 328 "dlasd7.f"
    givnum_dim1 = *ldgnum;
#line 328 "dlasd7.f"
    givnum_offset = 1 + givnum_dim1;
#line 328 "dlasd7.f"
    givnum -= givnum_offset;
#line 328 "dlasd7.f"

#line 328 "dlasd7.f"
    /* Function Body */
#line 328 "dlasd7.f"
    *info = 0;
#line 329 "dlasd7.f"
    n = *nl + *nr + 1;
#line 330 "dlasd7.f"
    m = n + *sqre;

#line 332 "dlasd7.f"
    if (*icompq < 0 || *icompq > 1) {
#line 333 "dlasd7.f"
	*info = -1;
#line 334 "dlasd7.f"
    } else if (*nl < 1) {
#line 335 "dlasd7.f"
	*info = -2;
#line 336 "dlasd7.f"
    } else if (*nr < 1) {
#line 337 "dlasd7.f"
	*info = -3;
#line 338 "dlasd7.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 339 "dlasd7.f"
	*info = -4;
#line 340 "dlasd7.f"
    } else if (*ldgcol < n) {
#line 341 "dlasd7.f"
	*info = -22;
#line 342 "dlasd7.f"
    } else if (*ldgnum < n) {
#line 343 "dlasd7.f"
	*info = -24;
#line 344 "dlasd7.f"
    }
#line 345 "dlasd7.f"
    if (*info != 0) {
#line 346 "dlasd7.f"
	i__1 = -(*info);
#line 346 "dlasd7.f"
	xerbla_("DLASD7", &i__1, (ftnlen)6);
#line 347 "dlasd7.f"
	return 0;
#line 348 "dlasd7.f"
    }

#line 350 "dlasd7.f"
    nlp1 = *nl + 1;
#line 351 "dlasd7.f"
    nlp2 = *nl + 2;
#line 352 "dlasd7.f"
    if (*icompq == 1) {
#line 353 "dlasd7.f"
	*givptr = 0;
#line 354 "dlasd7.f"
    }

/*     Generate the first part of the vector Z and move the singular */
/*     values in the first part of D one position backward. */

#line 359 "dlasd7.f"
    z1 = *alpha * vl[nlp1];
#line 360 "dlasd7.f"
    vl[nlp1] = 0.;
#line 361 "dlasd7.f"
    tau = vf[nlp1];
#line 362 "dlasd7.f"
    for (i__ = *nl; i__ >= 1; --i__) {
#line 363 "dlasd7.f"
	z__[i__ + 1] = *alpha * vl[i__];
#line 364 "dlasd7.f"
	vl[i__] = 0.;
#line 365 "dlasd7.f"
	vf[i__ + 1] = vf[i__];
#line 366 "dlasd7.f"
	d__[i__ + 1] = d__[i__];
#line 367 "dlasd7.f"
	idxq[i__ + 1] = idxq[i__] + 1;
#line 368 "dlasd7.f"
/* L10: */
#line 368 "dlasd7.f"
    }
#line 369 "dlasd7.f"
    vf[1] = tau;

/*     Generate the second part of the vector Z. */

#line 373 "dlasd7.f"
    i__1 = m;
#line 373 "dlasd7.f"
    for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 374 "dlasd7.f"
	z__[i__] = *beta * vf[i__];
#line 375 "dlasd7.f"
	vf[i__] = 0.;
#line 376 "dlasd7.f"
/* L20: */
#line 376 "dlasd7.f"
    }

/*     Sort the singular values into increasing order */

#line 380 "dlasd7.f"
    i__1 = n;
#line 380 "dlasd7.f"
    for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 381 "dlasd7.f"
	idxq[i__] += nlp1;
#line 382 "dlasd7.f"
/* L30: */
#line 382 "dlasd7.f"
    }

/*     DSIGMA, IDXC, IDXC, and ZW are used as storage space. */

#line 386 "dlasd7.f"
    i__1 = n;
#line 386 "dlasd7.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 387 "dlasd7.f"
	dsigma[i__] = d__[idxq[i__]];
#line 388 "dlasd7.f"
	zw[i__] = z__[idxq[i__]];
#line 389 "dlasd7.f"
	vfw[i__] = vf[idxq[i__]];
#line 390 "dlasd7.f"
	vlw[i__] = vl[idxq[i__]];
#line 391 "dlasd7.f"
/* L40: */
#line 391 "dlasd7.f"
    }

#line 393 "dlasd7.f"
    dlamrg_(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

#line 395 "dlasd7.f"
    i__1 = n;
#line 395 "dlasd7.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 396 "dlasd7.f"
	idxi = idx[i__] + 1;
#line 397 "dlasd7.f"
	d__[i__] = dsigma[idxi];
#line 398 "dlasd7.f"
	z__[i__] = zw[idxi];
#line 399 "dlasd7.f"
	vf[i__] = vfw[idxi];
#line 400 "dlasd7.f"
	vl[i__] = vlw[idxi];
#line 401 "dlasd7.f"
/* L50: */
#line 401 "dlasd7.f"
    }

/*     Calculate the allowable deflation tolerence */

#line 405 "dlasd7.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
/* Computing MAX */
#line 406 "dlasd7.f"
    d__1 = abs(*alpha), d__2 = abs(*beta);
#line 406 "dlasd7.f"
    tol = max(d__1,d__2);
/* Computing MAX */
#line 407 "dlasd7.f"
    d__2 = (d__1 = d__[n], abs(d__1));
#line 407 "dlasd7.f"
    tol = eps * 64. * max(d__2,tol);

/*     There are 2 kinds of deflation -- first a value in the z-vector */
/*     is small, second two (or more) singular values are very close */
/*     together (their difference is small). */

/*     If the value in the z-vector is small, we simply permute the */
/*     array so that the corresponding singular value is moved to the */
/*     end. */

/*     If two values in the D-vector are close, we perform a two-sided */
/*     rotation designed to make one of the corresponding z-vector */
/*     entries zero, and then permute the array so that the deflated */
/*     singular value is moved to the end. */

/*     If there are multiple singular values then the problem deflates. */
/*     Here the number of equal singular values are found.  As each equal */
/*     singular value is found, an elementary reflector is computed to */
/*     rotate the corresponding singular subspace so that the */
/*     corresponding components of Z are zero in this new basis. */

#line 428 "dlasd7.f"
    *k = 1;
#line 429 "dlasd7.f"
    k2 = n + 1;
#line 430 "dlasd7.f"
    i__1 = n;
#line 430 "dlasd7.f"
    for (j = 2; j <= i__1; ++j) {
#line 431 "dlasd7.f"
	if ((d__1 = z__[j], abs(d__1)) <= tol) {

/*           Deflate due to small z component. */

#line 435 "dlasd7.f"
	    --k2;
#line 436 "dlasd7.f"
	    idxp[k2] = j;
#line 437 "dlasd7.f"
	    if (j == n) {
#line 437 "dlasd7.f"
		goto L100;
#line 437 "dlasd7.f"
	    }
#line 439 "dlasd7.f"
	} else {
#line 440 "dlasd7.f"
	    jprev = j;
#line 441 "dlasd7.f"
	    goto L70;
#line 442 "dlasd7.f"
	}
#line 443 "dlasd7.f"
/* L60: */
#line 443 "dlasd7.f"
    }
#line 444 "dlasd7.f"
L70:
#line 445 "dlasd7.f"
    j = jprev;
#line 446 "dlasd7.f"
L80:
#line 447 "dlasd7.f"
    ++j;
#line 448 "dlasd7.f"
    if (j > n) {
#line 448 "dlasd7.f"
	goto L90;
#line 448 "dlasd7.f"
    }
#line 450 "dlasd7.f"
    if ((d__1 = z__[j], abs(d__1)) <= tol) {

/*        Deflate due to small z component. */

#line 454 "dlasd7.f"
	--k2;
#line 455 "dlasd7.f"
	idxp[k2] = j;
#line 456 "dlasd7.f"
    } else {

/*        Check if singular values are close enough to allow deflation. */

#line 460 "dlasd7.f"
	if ((d__1 = d__[j] - d__[jprev], abs(d__1)) <= tol) {

/*           Deflation is possible. */

#line 464 "dlasd7.f"
	    *s = z__[jprev];
#line 465 "dlasd7.f"
	    *c__ = z__[j];

/*           Find sqrt(a**2+b**2) without overflow or */
/*           destructive underflow. */

#line 470 "dlasd7.f"
	    tau = dlapy2_(c__, s);
#line 471 "dlasd7.f"
	    z__[j] = tau;
#line 472 "dlasd7.f"
	    z__[jprev] = 0.;
#line 473 "dlasd7.f"
	    *c__ /= tau;
#line 474 "dlasd7.f"
	    *s = -(*s) / tau;

/*           Record the appropriate Givens rotation */

#line 478 "dlasd7.f"
	    if (*icompq == 1) {
#line 479 "dlasd7.f"
		++(*givptr);
#line 480 "dlasd7.f"
		idxjp = idxq[idx[jprev] + 1];
#line 481 "dlasd7.f"
		idxj = idxq[idx[j] + 1];
#line 482 "dlasd7.f"
		if (idxjp <= nlp1) {
#line 483 "dlasd7.f"
		    --idxjp;
#line 484 "dlasd7.f"
		}
#line 485 "dlasd7.f"
		if (idxj <= nlp1) {
#line 486 "dlasd7.f"
		    --idxj;
#line 487 "dlasd7.f"
		}
#line 488 "dlasd7.f"
		givcol[*givptr + (givcol_dim1 << 1)] = idxjp;
#line 489 "dlasd7.f"
		givcol[*givptr + givcol_dim1] = idxj;
#line 490 "dlasd7.f"
		givnum[*givptr + (givnum_dim1 << 1)] = *c__;
#line 491 "dlasd7.f"
		givnum[*givptr + givnum_dim1] = *s;
#line 492 "dlasd7.f"
	    }
#line 493 "dlasd7.f"
	    drot_(&c__1, &vf[jprev], &c__1, &vf[j], &c__1, c__, s);
#line 494 "dlasd7.f"
	    drot_(&c__1, &vl[jprev], &c__1, &vl[j], &c__1, c__, s);
#line 495 "dlasd7.f"
	    --k2;
#line 496 "dlasd7.f"
	    idxp[k2] = jprev;
#line 497 "dlasd7.f"
	    jprev = j;
#line 498 "dlasd7.f"
	} else {
#line 499 "dlasd7.f"
	    ++(*k);
#line 500 "dlasd7.f"
	    zw[*k] = z__[jprev];
#line 501 "dlasd7.f"
	    dsigma[*k] = d__[jprev];
#line 502 "dlasd7.f"
	    idxp[*k] = jprev;
#line 503 "dlasd7.f"
	    jprev = j;
#line 504 "dlasd7.f"
	}
#line 505 "dlasd7.f"
    }
#line 506 "dlasd7.f"
    goto L80;
#line 507 "dlasd7.f"
L90:

/*     Record the last singular value. */

#line 511 "dlasd7.f"
    ++(*k);
#line 512 "dlasd7.f"
    zw[*k] = z__[jprev];
#line 513 "dlasd7.f"
    dsigma[*k] = d__[jprev];
#line 514 "dlasd7.f"
    idxp[*k] = jprev;

#line 516 "dlasd7.f"
L100:

/*     Sort the singular values into DSIGMA. The singular values which */
/*     were not deflated go into the first K slots of DSIGMA, except */
/*     that DSIGMA(1) is treated separately. */

#line 522 "dlasd7.f"
    i__1 = n;
#line 522 "dlasd7.f"
    for (j = 2; j <= i__1; ++j) {
#line 523 "dlasd7.f"
	jp = idxp[j];
#line 524 "dlasd7.f"
	dsigma[j] = d__[jp];
#line 525 "dlasd7.f"
	vfw[j] = vf[jp];
#line 526 "dlasd7.f"
	vlw[j] = vl[jp];
#line 527 "dlasd7.f"
/* L110: */
#line 527 "dlasd7.f"
    }
#line 528 "dlasd7.f"
    if (*icompq == 1) {
#line 529 "dlasd7.f"
	i__1 = n;
#line 529 "dlasd7.f"
	for (j = 2; j <= i__1; ++j) {
#line 530 "dlasd7.f"
	    jp = idxp[j];
#line 531 "dlasd7.f"
	    perm[j] = idxq[idx[jp] + 1];
#line 532 "dlasd7.f"
	    if (perm[j] <= nlp1) {
#line 533 "dlasd7.f"
		--perm[j];
#line 534 "dlasd7.f"
	    }
#line 535 "dlasd7.f"
/* L120: */
#line 535 "dlasd7.f"
	}
#line 536 "dlasd7.f"
    }

/*     The deflated singular values go back into the last N - K slots of */
/*     D. */

#line 541 "dlasd7.f"
    i__1 = n - *k;
#line 541 "dlasd7.f"
    dcopy_(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);

/*     Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and */
/*     VL(M). */

#line 546 "dlasd7.f"
    dsigma[1] = 0.;
#line 547 "dlasd7.f"
    hlftol = tol / 2.;
#line 548 "dlasd7.f"
    if (abs(dsigma[2]) <= hlftol) {
#line 548 "dlasd7.f"
	dsigma[2] = hlftol;
#line 548 "dlasd7.f"
    }
#line 550 "dlasd7.f"
    if (m > n) {
#line 551 "dlasd7.f"
	z__[1] = dlapy2_(&z1, &z__[m]);
#line 552 "dlasd7.f"
	if (z__[1] <= tol) {
#line 553 "dlasd7.f"
	    *c__ = 1.;
#line 554 "dlasd7.f"
	    *s = 0.;
#line 555 "dlasd7.f"
	    z__[1] = tol;
#line 556 "dlasd7.f"
	} else {
#line 557 "dlasd7.f"
	    *c__ = z1 / z__[1];
#line 558 "dlasd7.f"
	    *s = -z__[m] / z__[1];
#line 559 "dlasd7.f"
	}
#line 560 "dlasd7.f"
	drot_(&c__1, &vf[m], &c__1, &vf[1], &c__1, c__, s);
#line 561 "dlasd7.f"
	drot_(&c__1, &vl[m], &c__1, &vl[1], &c__1, c__, s);
#line 562 "dlasd7.f"
    } else {
#line 563 "dlasd7.f"
	if (abs(z1) <= tol) {
#line 564 "dlasd7.f"
	    z__[1] = tol;
#line 565 "dlasd7.f"
	} else {
#line 566 "dlasd7.f"
	    z__[1] = z1;
#line 567 "dlasd7.f"
	}
#line 568 "dlasd7.f"
    }

/*     Restore Z, VF, and VL. */

#line 572 "dlasd7.f"
    i__1 = *k - 1;
#line 572 "dlasd7.f"
    dcopy_(&i__1, &zw[2], &c__1, &z__[2], &c__1);
#line 573 "dlasd7.f"
    i__1 = n - 1;
#line 573 "dlasd7.f"
    dcopy_(&i__1, &vfw[2], &c__1, &vf[2], &c__1);
#line 574 "dlasd7.f"
    i__1 = n - 1;
#line 574 "dlasd7.f"
    dcopy_(&i__1, &vlw[2], &c__1, &vl[2], &c__1);

#line 576 "dlasd7.f"
    return 0;

/*     End of DLASD7 */

} /* dlasd7_ */

