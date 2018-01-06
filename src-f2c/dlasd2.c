#line 1 "dlasd2.f"
/* dlasd2.f -- translated by f2c (version 20100827).
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

#line 1 "dlasd2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b30 = 0.;

/* > \brief \b DLASD2 merges the two sets of singular values together into a single sorted set. Used by sbdsdc
. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT, */
/*                          LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX, */
/*                          IDXC, IDXQ, COLTYP, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            COLTYP( * ), IDX( * ), IDXC( * ), IDXP( * ), */
/*      $                   IDXQ( * ) */
/*       DOUBLE PRECISION   D( * ), DSIGMA( * ), U( LDU, * ), */
/*      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), */
/*      $                   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASD2 merges the two sets of singular values together into a single */
/* > sorted set.  Then it tries to deflate the size of the problem. */
/* > There are two ways in which deflation can occur:  when two or more */
/* > singular values are close together or if there is a tiny entry in the */
/* > Z vector.  For each such occurrence the order of the related secular */
/* > equation problem is reduced by one. */
/* > */
/* > DLASD2 is called from DLASD1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NL */
/* > \verbatim */
/* >          NL is INTEGER */
/* >         The row dimension of the upper block.  NL >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NR */
/* > \verbatim */
/* >          NR is INTEGER */
/* >         The row dimension of the lower block.  NR >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* >          SQRE is INTEGER */
/* >         = 0: the lower block is an NR-by-NR square matrix. */
/* >         = 1: the lower block is an NR-by-(NR+1) rectangular matrix. */
/* > */
/* >         The bidiagonal matrix has N = NL + NR + 1 rows and */
/* >         M = N + SQRE >= N columns. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         Contains the dimension of the non-deflated matrix, */
/* >         This is the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension(N) */
/* >         On entry D contains the singular values of the two submatrices */
/* >         to be combined.  On exit D contains the trailing (N-K) updated */
/* >         singular values (those which were deflated) sorted into */
/* >         increasing order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension(N) */
/* >         On exit Z contains the updating row vector in the secular */
/* >         equation. */
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
/* > \param[in,out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension(LDU,N) */
/* >         On entry U contains the left singular vectors of two */
/* >         submatrices in the two square blocks with corners at (1,1), */
/* >         (NL, NL), and (NL+2, NL+2), (N,N). */
/* >         On exit U contains the trailing (N-K) updated left singular */
/* >         vectors (those which were deflated) in its last N-K columns. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >         The leading dimension of the array U.  LDU >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension(LDVT,M) */
/* >         On entry VT**T contains the right singular vectors of two */
/* >         submatrices in the two square blocks with corners at (1,1), */
/* >         (NL+1, NL+1), and (NL+2, NL+2), (M,M). */
/* >         On exit VT**T contains the trailing (N-K) updated right singular */
/* >         vectors (those which were deflated) in its last N-K columns. */
/* >         In case SQRE =1, the last row of VT spans the right null */
/* >         space. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >         The leading dimension of the array VT.  LDVT >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] DSIGMA */
/* > \verbatim */
/* >          DSIGMA is DOUBLE PRECISION array, dimension (N) */
/* >         Contains a copy of the diagonal elements (K-1 singular values */
/* >         and one zero) in the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* >          U2 is DOUBLE PRECISION array, dimension(LDU2,N) */
/* >         Contains a copy of the first K-1 left singular vectors which */
/* >         will be used by DLASD3 in a matrix multiply (DGEMM) to solve */
/* >         for the new left singular vectors. U2 is arranged into four */
/* >         blocks. The first block contains a column with 1 at NL+1 and */
/* >         zero everywhere else; the second block contains non-zero */
/* >         entries only at and above NL; the third contains non-zero */
/* >         entries only below NL+1; and the fourth is dense. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* >          LDU2 is INTEGER */
/* >         The leading dimension of the array U2.  LDU2 >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VT2 */
/* > \verbatim */
/* >          VT2 is DOUBLE PRECISION array, dimension(LDVT2,N) */
/* >         VT2**T contains a copy of the first K right singular vectors */
/* >         which will be used by DLASD3 in a matrix multiply (DGEMM) to */
/* >         solve for the new right singular vectors. VT2 is arranged into */
/* >         three blocks. The first block contains a row that corresponds */
/* >         to the special 0 diagonal element in SIGMA; the second block */
/* >         contains non-zeros only at and before NL +1; the third block */
/* >         contains non-zeros only at and after  NL +2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT2 */
/* > \verbatim */
/* >          LDVT2 is INTEGER */
/* >         The leading dimension of the array VT2.  LDVT2 >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXP */
/* > \verbatim */
/* >          IDXP is INTEGER array dimension(N) */
/* >         This will contain the permutation used to place deflated */
/* >         values of D at the end of the array. On output IDXP(2:K) */
/* >         points to the nondeflated D-values and IDXP(K+1:N) */
/* >         points to the deflated singular values. */
/* > \endverbatim */
/* > */
/* > \param[out] IDX */
/* > \verbatim */
/* >          IDX is INTEGER array dimension(N) */
/* >         This will contain the permutation used to sort the contents of */
/* >         D into ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXC */
/* > \verbatim */
/* >          IDXC is INTEGER array dimension(N) */
/* >         This will contain the permutation used to arrange the columns */
/* >         of the deflated U matrix into three groups:  the first group */
/* >         contains non-zero entries only at and above NL, the second */
/* >         contains non-zero entries only below NL+2, and the third is */
/* >         dense. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IDXQ */
/* > \verbatim */
/* >          IDXQ is INTEGER array dimension(N) */
/* >         This contains the permutation which separately sorts the two */
/* >         sub-problems in D into ascending order.  Note that entries in */
/* >         the first hlaf of this permutation must first be moved one */
/* >         position backward; and entries in the second half */
/* >         must first have NL+1 added to their values. */
/* > \endverbatim */
/* > */
/* > \param[out] COLTYP */
/* > \verbatim */
/* >          COLTYP is INTEGER array dimension(N) */
/* >         As workspace, this will contain a label which will indicate */
/* >         which of the following types a column in the U2 matrix or a */
/* >         row in the VT2 matrix is: */
/* >         1 : non-zero in the upper half only */
/* >         2 : non-zero in the lower half only */
/* >         3 : dense */
/* >         4 : deflated */
/* > */
/* >         On exit, it is an array of dimension 4, with COLTYP(I) being */
/* >         the dimension of the I-th type columns. */
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

/* > \ingroup auxOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasd2_(integer *nl, integer *nr, integer *sqre, integer 
	*k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *
	beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, 
	doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, 
	integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
	idxq, integer *coltyp, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset, 
	    vt2_dim1, vt2_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j, m, n;
    static doublereal s;
    static integer k2;
    static doublereal z1;
    static integer ct, jp;
    static doublereal eps, tau, tol;
    static integer psm[4], nlp1, nlp2, idxi, idxj;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer ctot[4], idxjp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jprev;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal hlftol;


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
/*     .. Local Arrays .. */
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

/*     Test the input parameters. */

#line 318 "dlasd2.f"
    /* Parameter adjustments */
#line 318 "dlasd2.f"
    --d__;
#line 318 "dlasd2.f"
    --z__;
#line 318 "dlasd2.f"
    u_dim1 = *ldu;
#line 318 "dlasd2.f"
    u_offset = 1 + u_dim1;
#line 318 "dlasd2.f"
    u -= u_offset;
#line 318 "dlasd2.f"
    vt_dim1 = *ldvt;
#line 318 "dlasd2.f"
    vt_offset = 1 + vt_dim1;
#line 318 "dlasd2.f"
    vt -= vt_offset;
#line 318 "dlasd2.f"
    --dsigma;
#line 318 "dlasd2.f"
    u2_dim1 = *ldu2;
#line 318 "dlasd2.f"
    u2_offset = 1 + u2_dim1;
#line 318 "dlasd2.f"
    u2 -= u2_offset;
#line 318 "dlasd2.f"
    vt2_dim1 = *ldvt2;
#line 318 "dlasd2.f"
    vt2_offset = 1 + vt2_dim1;
#line 318 "dlasd2.f"
    vt2 -= vt2_offset;
#line 318 "dlasd2.f"
    --idxp;
#line 318 "dlasd2.f"
    --idx;
#line 318 "dlasd2.f"
    --idxc;
#line 318 "dlasd2.f"
    --idxq;
#line 318 "dlasd2.f"
    --coltyp;
#line 318 "dlasd2.f"

#line 318 "dlasd2.f"
    /* Function Body */
#line 318 "dlasd2.f"
    *info = 0;

#line 320 "dlasd2.f"
    if (*nl < 1) {
#line 321 "dlasd2.f"
	*info = -1;
#line 322 "dlasd2.f"
    } else if (*nr < 1) {
#line 323 "dlasd2.f"
	*info = -2;
#line 324 "dlasd2.f"
    } else if (*sqre != 1 && *sqre != 0) {
#line 325 "dlasd2.f"
	*info = -3;
#line 326 "dlasd2.f"
    }

#line 328 "dlasd2.f"
    n = *nl + *nr + 1;
#line 329 "dlasd2.f"
    m = n + *sqre;

#line 331 "dlasd2.f"
    if (*ldu < n) {
#line 332 "dlasd2.f"
	*info = -10;
#line 333 "dlasd2.f"
    } else if (*ldvt < m) {
#line 334 "dlasd2.f"
	*info = -12;
#line 335 "dlasd2.f"
    } else if (*ldu2 < n) {
#line 336 "dlasd2.f"
	*info = -15;
#line 337 "dlasd2.f"
    } else if (*ldvt2 < m) {
#line 338 "dlasd2.f"
	*info = -17;
#line 339 "dlasd2.f"
    }
#line 340 "dlasd2.f"
    if (*info != 0) {
#line 341 "dlasd2.f"
	i__1 = -(*info);
#line 341 "dlasd2.f"
	xerbla_("DLASD2", &i__1, (ftnlen)6);
#line 342 "dlasd2.f"
	return 0;
#line 343 "dlasd2.f"
    }

#line 345 "dlasd2.f"
    nlp1 = *nl + 1;
#line 346 "dlasd2.f"
    nlp2 = *nl + 2;

/*     Generate the first part of the vector Z; and move the singular */
/*     values in the first part of D one position backward. */

#line 351 "dlasd2.f"
    z1 = *alpha * vt[nlp1 + nlp1 * vt_dim1];
#line 352 "dlasd2.f"
    z__[1] = z1;
#line 353 "dlasd2.f"
    for (i__ = *nl; i__ >= 1; --i__) {
#line 354 "dlasd2.f"
	z__[i__ + 1] = *alpha * vt[i__ + nlp1 * vt_dim1];
#line 355 "dlasd2.f"
	d__[i__ + 1] = d__[i__];
#line 356 "dlasd2.f"
	idxq[i__ + 1] = idxq[i__] + 1;
#line 357 "dlasd2.f"
/* L10: */
#line 357 "dlasd2.f"
    }

/*     Generate the second part of the vector Z. */

#line 361 "dlasd2.f"
    i__1 = m;
#line 361 "dlasd2.f"
    for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 362 "dlasd2.f"
	z__[i__] = *beta * vt[i__ + nlp2 * vt_dim1];
#line 363 "dlasd2.f"
/* L20: */
#line 363 "dlasd2.f"
    }

/*     Initialize some reference arrays. */

#line 367 "dlasd2.f"
    i__1 = nlp1;
#line 367 "dlasd2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 368 "dlasd2.f"
	coltyp[i__] = 1;
#line 369 "dlasd2.f"
/* L30: */
#line 369 "dlasd2.f"
    }
#line 370 "dlasd2.f"
    i__1 = n;
#line 370 "dlasd2.f"
    for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 371 "dlasd2.f"
	coltyp[i__] = 2;
#line 372 "dlasd2.f"
/* L40: */
#line 372 "dlasd2.f"
    }

/*     Sort the singular values into increasing order */

#line 376 "dlasd2.f"
    i__1 = n;
#line 376 "dlasd2.f"
    for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 377 "dlasd2.f"
	idxq[i__] += nlp1;
#line 378 "dlasd2.f"
/* L50: */
#line 378 "dlasd2.f"
    }

/*     DSIGMA, IDXC, IDXC, and the first column of U2 */
/*     are used as storage space. */

#line 383 "dlasd2.f"
    i__1 = n;
#line 383 "dlasd2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 384 "dlasd2.f"
	dsigma[i__] = d__[idxq[i__]];
#line 385 "dlasd2.f"
	u2[i__ + u2_dim1] = z__[idxq[i__]];
#line 386 "dlasd2.f"
	idxc[i__] = coltyp[idxq[i__]];
#line 387 "dlasd2.f"
/* L60: */
#line 387 "dlasd2.f"
    }

#line 389 "dlasd2.f"
    dlamrg_(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

#line 391 "dlasd2.f"
    i__1 = n;
#line 391 "dlasd2.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 392 "dlasd2.f"
	idxi = idx[i__] + 1;
#line 393 "dlasd2.f"
	d__[i__] = dsigma[idxi];
#line 394 "dlasd2.f"
	z__[i__] = u2[idxi + u2_dim1];
#line 395 "dlasd2.f"
	coltyp[i__] = idxc[idxi];
#line 396 "dlasd2.f"
/* L70: */
#line 396 "dlasd2.f"
    }

/*     Calculate the allowable deflation tolerance */

#line 400 "dlasd2.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
/* Computing MAX */
#line 401 "dlasd2.f"
    d__1 = abs(*alpha), d__2 = abs(*beta);
#line 401 "dlasd2.f"
    tol = max(d__1,d__2);
/* Computing MAX */
#line 402 "dlasd2.f"
    d__2 = (d__1 = d__[n], abs(d__1));
#line 402 "dlasd2.f"
    tol = eps * 8. * max(d__2,tol);

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

#line 423 "dlasd2.f"
    *k = 1;
#line 424 "dlasd2.f"
    k2 = n + 1;
#line 425 "dlasd2.f"
    i__1 = n;
#line 425 "dlasd2.f"
    for (j = 2; j <= i__1; ++j) {
#line 426 "dlasd2.f"
	if ((d__1 = z__[j], abs(d__1)) <= tol) {

/*           Deflate due to small z component. */

#line 430 "dlasd2.f"
	    --k2;
#line 431 "dlasd2.f"
	    idxp[k2] = j;
#line 432 "dlasd2.f"
	    coltyp[j] = 4;
#line 433 "dlasd2.f"
	    if (j == n) {
#line 433 "dlasd2.f"
		goto L120;
#line 433 "dlasd2.f"
	    }
#line 435 "dlasd2.f"
	} else {
#line 436 "dlasd2.f"
	    jprev = j;
#line 437 "dlasd2.f"
	    goto L90;
#line 438 "dlasd2.f"
	}
#line 439 "dlasd2.f"
/* L80: */
#line 439 "dlasd2.f"
    }
#line 440 "dlasd2.f"
L90:
#line 441 "dlasd2.f"
    j = jprev;
#line 442 "dlasd2.f"
L100:
#line 443 "dlasd2.f"
    ++j;
#line 444 "dlasd2.f"
    if (j > n) {
#line 444 "dlasd2.f"
	goto L110;
#line 444 "dlasd2.f"
    }
#line 446 "dlasd2.f"
    if ((d__1 = z__[j], abs(d__1)) <= tol) {

/*        Deflate due to small z component. */

#line 450 "dlasd2.f"
	--k2;
#line 451 "dlasd2.f"
	idxp[k2] = j;
#line 452 "dlasd2.f"
	coltyp[j] = 4;
#line 453 "dlasd2.f"
    } else {

/*        Check if singular values are close enough to allow deflation. */

#line 457 "dlasd2.f"
	if ((d__1 = d__[j] - d__[jprev], abs(d__1)) <= tol) {

/*           Deflation is possible. */

#line 461 "dlasd2.f"
	    s = z__[jprev];
#line 462 "dlasd2.f"
	    c__ = z__[j];

/*           Find sqrt(a**2+b**2) without overflow or */
/*           destructive underflow. */

#line 467 "dlasd2.f"
	    tau = dlapy2_(&c__, &s);
#line 468 "dlasd2.f"
	    c__ /= tau;
#line 469 "dlasd2.f"
	    s = -s / tau;
#line 470 "dlasd2.f"
	    z__[j] = tau;
#line 471 "dlasd2.f"
	    z__[jprev] = 0.;

/*           Apply back the Givens rotation to the left and right */
/*           singular vector matrices. */

#line 476 "dlasd2.f"
	    idxjp = idxq[idx[jprev] + 1];
#line 477 "dlasd2.f"
	    idxj = idxq[idx[j] + 1];
#line 478 "dlasd2.f"
	    if (idxjp <= nlp1) {
#line 479 "dlasd2.f"
		--idxjp;
#line 480 "dlasd2.f"
	    }
#line 481 "dlasd2.f"
	    if (idxj <= nlp1) {
#line 482 "dlasd2.f"
		--idxj;
#line 483 "dlasd2.f"
	    }
#line 484 "dlasd2.f"
	    drot_(&n, &u[idxjp * u_dim1 + 1], &c__1, &u[idxj * u_dim1 + 1], &
		    c__1, &c__, &s);
#line 485 "dlasd2.f"
	    drot_(&m, &vt[idxjp + vt_dim1], ldvt, &vt[idxj + vt_dim1], ldvt, &
		    c__, &s);
#line 487 "dlasd2.f"
	    if (coltyp[j] != coltyp[jprev]) {
#line 488 "dlasd2.f"
		coltyp[j] = 3;
#line 489 "dlasd2.f"
	    }
#line 490 "dlasd2.f"
	    coltyp[jprev] = 4;
#line 491 "dlasd2.f"
	    --k2;
#line 492 "dlasd2.f"
	    idxp[k2] = jprev;
#line 493 "dlasd2.f"
	    jprev = j;
#line 494 "dlasd2.f"
	} else {
#line 495 "dlasd2.f"
	    ++(*k);
#line 496 "dlasd2.f"
	    u2[*k + u2_dim1] = z__[jprev];
#line 497 "dlasd2.f"
	    dsigma[*k] = d__[jprev];
#line 498 "dlasd2.f"
	    idxp[*k] = jprev;
#line 499 "dlasd2.f"
	    jprev = j;
#line 500 "dlasd2.f"
	}
#line 501 "dlasd2.f"
    }
#line 502 "dlasd2.f"
    goto L100;
#line 503 "dlasd2.f"
L110:

/*     Record the last singular value. */

#line 507 "dlasd2.f"
    ++(*k);
#line 508 "dlasd2.f"
    u2[*k + u2_dim1] = z__[jprev];
#line 509 "dlasd2.f"
    dsigma[*k] = d__[jprev];
#line 510 "dlasd2.f"
    idxp[*k] = jprev;

#line 512 "dlasd2.f"
L120:

/*     Count up the total number of the various types of columns, then */
/*     form a permutation which positions the four column types into */
/*     four groups of uniform structure (although one or more of these */
/*     groups may be empty). */

#line 519 "dlasd2.f"
    for (j = 1; j <= 4; ++j) {
#line 520 "dlasd2.f"
	ctot[j - 1] = 0;
#line 521 "dlasd2.f"
/* L130: */
#line 521 "dlasd2.f"
    }
#line 522 "dlasd2.f"
    i__1 = n;
#line 522 "dlasd2.f"
    for (j = 2; j <= i__1; ++j) {
#line 523 "dlasd2.f"
	ct = coltyp[j];
#line 524 "dlasd2.f"
	++ctot[ct - 1];
#line 525 "dlasd2.f"
/* L140: */
#line 525 "dlasd2.f"
    }

/*     PSM(*) = Position in SubMatrix (of types 1 through 4) */

#line 529 "dlasd2.f"
    psm[0] = 2;
#line 530 "dlasd2.f"
    psm[1] = ctot[0] + 2;
#line 531 "dlasd2.f"
    psm[2] = psm[1] + ctot[1];
#line 532 "dlasd2.f"
    psm[3] = psm[2] + ctot[2];

/*     Fill out the IDXC array so that the permutation which it induces */
/*     will place all type-1 columns first, all type-2 columns next, */
/*     then all type-3's, and finally all type-4's, starting from the */
/*     second column. This applies similarly to the rows of VT. */

#line 539 "dlasd2.f"
    i__1 = n;
#line 539 "dlasd2.f"
    for (j = 2; j <= i__1; ++j) {
#line 540 "dlasd2.f"
	jp = idxp[j];
#line 541 "dlasd2.f"
	ct = coltyp[jp];
#line 542 "dlasd2.f"
	idxc[psm[ct - 1]] = j;
#line 543 "dlasd2.f"
	++psm[ct - 1];
#line 544 "dlasd2.f"
/* L150: */
#line 544 "dlasd2.f"
    }

/*     Sort the singular values and corresponding singular vectors into */
/*     DSIGMA, U2, and VT2 respectively.  The singular values/vectors */
/*     which were not deflated go into the first K slots of DSIGMA, U2, */
/*     and VT2 respectively, while those which were deflated go into the */
/*     last N - K slots, except that the first column/row will be treated */
/*     separately. */

#line 553 "dlasd2.f"
    i__1 = n;
#line 553 "dlasd2.f"
    for (j = 2; j <= i__1; ++j) {
#line 554 "dlasd2.f"
	jp = idxp[j];
#line 555 "dlasd2.f"
	dsigma[j] = d__[jp];
#line 556 "dlasd2.f"
	idxj = idxq[idx[idxp[idxc[j]]] + 1];
#line 557 "dlasd2.f"
	if (idxj <= nlp1) {
#line 558 "dlasd2.f"
	    --idxj;
#line 559 "dlasd2.f"
	}
#line 560 "dlasd2.f"
	dcopy_(&n, &u[idxj * u_dim1 + 1], &c__1, &u2[j * u2_dim1 + 1], &c__1);
#line 561 "dlasd2.f"
	dcopy_(&m, &vt[idxj + vt_dim1], ldvt, &vt2[j + vt2_dim1], ldvt2);
#line 562 "dlasd2.f"
/* L160: */
#line 562 "dlasd2.f"
    }

/*     Determine DSIGMA(1), DSIGMA(2) and Z(1) */

#line 566 "dlasd2.f"
    dsigma[1] = 0.;
#line 567 "dlasd2.f"
    hlftol = tol / 2.;
#line 568 "dlasd2.f"
    if (abs(dsigma[2]) <= hlftol) {
#line 568 "dlasd2.f"
	dsigma[2] = hlftol;
#line 568 "dlasd2.f"
    }
#line 570 "dlasd2.f"
    if (m > n) {
#line 571 "dlasd2.f"
	z__[1] = dlapy2_(&z1, &z__[m]);
#line 572 "dlasd2.f"
	if (z__[1] <= tol) {
#line 573 "dlasd2.f"
	    c__ = 1.;
#line 574 "dlasd2.f"
	    s = 0.;
#line 575 "dlasd2.f"
	    z__[1] = tol;
#line 576 "dlasd2.f"
	} else {
#line 577 "dlasd2.f"
	    c__ = z1 / z__[1];
#line 578 "dlasd2.f"
	    s = z__[m] / z__[1];
#line 579 "dlasd2.f"
	}
#line 580 "dlasd2.f"
    } else {
#line 581 "dlasd2.f"
	if (abs(z1) <= tol) {
#line 582 "dlasd2.f"
	    z__[1] = tol;
#line 583 "dlasd2.f"
	} else {
#line 584 "dlasd2.f"
	    z__[1] = z1;
#line 585 "dlasd2.f"
	}
#line 586 "dlasd2.f"
    }

/*     Move the rest of the updating row to Z. */

#line 590 "dlasd2.f"
    i__1 = *k - 1;
#line 590 "dlasd2.f"
    dcopy_(&i__1, &u2[u2_dim1 + 2], &c__1, &z__[2], &c__1);

/*     Determine the first column of U2, the first row of VT2 and the */
/*     last row of VT. */

#line 595 "dlasd2.f"
    dlaset_("A", &n, &c__1, &c_b30, &c_b30, &u2[u2_offset], ldu2, (ftnlen)1);
#line 596 "dlasd2.f"
    u2[nlp1 + u2_dim1] = 1.;
#line 597 "dlasd2.f"
    if (m > n) {
#line 598 "dlasd2.f"
	i__1 = nlp1;
#line 598 "dlasd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 599 "dlasd2.f"
	    vt[m + i__ * vt_dim1] = -s * vt[nlp1 + i__ * vt_dim1];
#line 600 "dlasd2.f"
	    vt2[i__ * vt2_dim1 + 1] = c__ * vt[nlp1 + i__ * vt_dim1];
#line 601 "dlasd2.f"
/* L170: */
#line 601 "dlasd2.f"
	}
#line 602 "dlasd2.f"
	i__1 = m;
#line 602 "dlasd2.f"
	for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 603 "dlasd2.f"
	    vt2[i__ * vt2_dim1 + 1] = s * vt[m + i__ * vt_dim1];
#line 604 "dlasd2.f"
	    vt[m + i__ * vt_dim1] = c__ * vt[m + i__ * vt_dim1];
#line 605 "dlasd2.f"
/* L180: */
#line 605 "dlasd2.f"
	}
#line 606 "dlasd2.f"
    } else {
#line 607 "dlasd2.f"
	dcopy_(&m, &vt[nlp1 + vt_dim1], ldvt, &vt2[vt2_dim1 + 1], ldvt2);
#line 608 "dlasd2.f"
    }
#line 609 "dlasd2.f"
    if (m > n) {
#line 610 "dlasd2.f"
	dcopy_(&m, &vt[m + vt_dim1], ldvt, &vt2[m + vt2_dim1], ldvt2);
#line 611 "dlasd2.f"
    }

/*     The deflated singular values and their corresponding vectors go */
/*     into the back of D, U, and V respectively. */

#line 616 "dlasd2.f"
    if (n > *k) {
#line 617 "dlasd2.f"
	i__1 = n - *k;
#line 617 "dlasd2.f"
	dcopy_(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
#line 618 "dlasd2.f"
	i__1 = n - *k;
#line 618 "dlasd2.f"
	dlacpy_("A", &n, &i__1, &u2[(*k + 1) * u2_dim1 + 1], ldu2, &u[(*k + 1)
		 * u_dim1 + 1], ldu, (ftnlen)1);
#line 620 "dlasd2.f"
	i__1 = n - *k;
#line 620 "dlasd2.f"
	dlacpy_("A", &i__1, &m, &vt2[*k + 1 + vt2_dim1], ldvt2, &vt[*k + 1 + 
		vt_dim1], ldvt, (ftnlen)1);
#line 622 "dlasd2.f"
    }

/*     Copy CTOT into COLTYP for referencing in DLASD3. */

#line 626 "dlasd2.f"
    for (j = 1; j <= 4; ++j) {
#line 627 "dlasd2.f"
	coltyp[j] = ctot[j - 1];
#line 628 "dlasd2.f"
/* L190: */
#line 628 "dlasd2.f"
    }

#line 630 "dlasd2.f"
    return 0;

/*     End of DLASD2 */

} /* dlasd2_ */

