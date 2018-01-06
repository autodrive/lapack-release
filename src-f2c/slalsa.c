#line 1 "slalsa.f"
/* slalsa.f -- translated by f2c (version 20100827).
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

#line 1 "slalsa.f"
/* Table of constant values */

static doublereal c_b7 = 1.;
static doublereal c_b8 = 0.;
static integer c__2 = 2;

/* > \brief \b SLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLALSA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slalsa.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slalsa.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slalsa.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, */
/*                          LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, */
/*                          GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, */
/*      $                   SMLSIZ */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), */
/*      $                   K( * ), PERM( LDGCOL, * ) */
/*       REAL               B( LDB, * ), BX( LDBX, * ), C( * ), */
/*      $                   DIFL( LDU, * ), DIFR( LDU, * ), */
/*      $                   GIVNUM( LDU, * ), POLES( LDU, * ), S( * ), */
/*      $                   U( LDU, * ), VT( LDU, * ), WORK( * ), */
/*      $                   Z( LDU, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALSA is an itermediate step in solving the least squares problem */
/* > by computing the SVD of the coefficient matrix in compact form (The */
/* > singular vectors are computed as products of simple orthorgonal */
/* > matrices.). */
/* > */
/* > If ICOMPQ = 0, SLALSA applies the inverse of the left singular vector */
/* > matrix of an upper bidiagonal matrix to the right hand side; and if */
/* > ICOMPQ = 1, SLALSA applies the right singular vector matrix to the */
/* > right hand side. The singular vector matrices were generated in */
/* > compact form by SLALSA. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >         Specifies whether the left or the right singular vector */
/* >         matrix is involved. */
/* >         = 0: Left singular vector matrix */
/* >         = 1: Right singular vector matrix */
/* > \endverbatim */
/* > */
/* > \param[in] SMLSIZ */
/* > \verbatim */
/* >          SMLSIZ is INTEGER */
/* >         The maximum size of the subproblems at the bottom of the */
/* >         computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The row and column dimensions of the upper bidiagonal matrix. */
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
/* >         squares problem in rows 1 through M. */
/* >         On output, B contains the solution X in rows 1 through N. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >         The leading dimension of B in the calling subprogram. */
/* >         LDB must be at least max(1,MAX( M, N ) ). */
/* > \endverbatim */
/* > */
/* > \param[out] BX */
/* > \verbatim */
/* >          BX is REAL array, dimension ( LDBX, NRHS ) */
/* >         On exit, the result of applying the left or right singular */
/* >         vector matrix to B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBX */
/* > \verbatim */
/* >          LDBX is INTEGER */
/* >         The leading dimension of BX. */
/* > \endverbatim */
/* > */
/* > \param[in] U */
/* > \verbatim */
/* >          U is REAL array, dimension ( LDU, SMLSIZ ). */
/* >         On entry, U contains the left singular vector matrices of all */
/* >         subproblems at the bottom level. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER, LDU = > N. */
/* >         The leading dimension of arrays U, VT, DIFL, DIFR, */
/* >         POLES, GIVNUM, and Z. */
/* > \endverbatim */
/* > */
/* > \param[in] VT */
/* > \verbatim */
/* >          VT is REAL array, dimension ( LDU, SMLSIZ+1 ). */
/* >         On entry, VT**T contains the right singular vector matrices of */
/* >         all subproblems at the bottom level. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER array, dimension ( N ). */
/* > \endverbatim */
/* > */
/* > \param[in] DIFL */
/* > \verbatim */
/* >          DIFL is REAL array, dimension ( LDU, NLVL ). */
/* >         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFR */
/* > \verbatim */
/* >          DIFR is REAL array, dimension ( LDU, 2 * NLVL ). */
/* >         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record */
/* >         distances between singular values on the I-th level and */
/* >         singular values on the (I -1)-th level, and DIFR(*, 2 * I) */
/* >         record the normalizing factors of the right singular vectors */
/* >         matrices of subproblems on I-th level. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension ( LDU, NLVL ). */
/* >         On entry, Z(1, I) contains the components of the deflation- */
/* >         adjusted updating row vector for subproblems on the I-th */
/* >         level. */
/* > \endverbatim */
/* > */
/* > \param[in] POLES */
/* > \verbatim */
/* >          POLES is REAL array, dimension ( LDU, 2 * NLVL ). */
/* >         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old */
/* >         singular values involved in the secular equations on the I-th */
/* >         level. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER array, dimension ( N ). */
/* >         On entry, GIVPTR( I ) records the number of Givens */
/* >         rotations performed on the I-th problem on the computation */
/* >         tree. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, dimension ( LDGCOL, 2 * NLVL ). */
/* >         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the */
/* >         locations of Givens rotations performed on the I-th level on */
/* >         the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* >          LDGCOL is INTEGER, LDGCOL = > N. */
/* >         The leading dimension of arrays GIVCOL and PERM. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension ( LDGCOL, NLVL ). */
/* >         On entry, PERM(*, I) records permutations done on the I-th */
/* >         level of the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is REAL array, dimension ( LDU, 2 * NLVL ). */
/* >         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S- */
/* >         values of Givens rotations performed on the I-th level on the */
/* >         computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension ( N ). */
/* >         On entry, if the I-th subproblem is not square, */
/* >         C( I ) contains the C-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL array, dimension ( N ). */
/* >         On entry, if the I-th subproblem is not square, */
/* >         S( I ) contains the S-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array. */
/* >         The dimension must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array. */
/* >         The dimension must be at least 3 * N */
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

/* > \ingroup realOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int slalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *
	ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k, 
	doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, b_dim1, 
	    b_offset, bx_dim1, bx_offset, difl_dim1, difl_offset, difr_dim1, 
	    difr_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset,
	     u_dim1, u_offset, vt_dim1, vt_offset, z_dim1, z_offset, i__1, 
	    i__2;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, i1, ic, lf, nd, ll, nl, nr, im1, nlf, nrf, lvl, 
	    ndb1, nlp1, lvl2, nrp1, nlvl, sqre, inode, ndiml;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer ndimr;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slals0_(integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *), xerbla_(char *, integer *, ftnlen), slasdt_(integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    );


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 309 "slalsa.f"
    /* Parameter adjustments */
#line 309 "slalsa.f"
    b_dim1 = *ldb;
#line 309 "slalsa.f"
    b_offset = 1 + b_dim1;
#line 309 "slalsa.f"
    b -= b_offset;
#line 309 "slalsa.f"
    bx_dim1 = *ldbx;
#line 309 "slalsa.f"
    bx_offset = 1 + bx_dim1;
#line 309 "slalsa.f"
    bx -= bx_offset;
#line 309 "slalsa.f"
    givnum_dim1 = *ldu;
#line 309 "slalsa.f"
    givnum_offset = 1 + givnum_dim1;
#line 309 "slalsa.f"
    givnum -= givnum_offset;
#line 309 "slalsa.f"
    poles_dim1 = *ldu;
#line 309 "slalsa.f"
    poles_offset = 1 + poles_dim1;
#line 309 "slalsa.f"
    poles -= poles_offset;
#line 309 "slalsa.f"
    z_dim1 = *ldu;
#line 309 "slalsa.f"
    z_offset = 1 + z_dim1;
#line 309 "slalsa.f"
    z__ -= z_offset;
#line 309 "slalsa.f"
    difr_dim1 = *ldu;
#line 309 "slalsa.f"
    difr_offset = 1 + difr_dim1;
#line 309 "slalsa.f"
    difr -= difr_offset;
#line 309 "slalsa.f"
    difl_dim1 = *ldu;
#line 309 "slalsa.f"
    difl_offset = 1 + difl_dim1;
#line 309 "slalsa.f"
    difl -= difl_offset;
#line 309 "slalsa.f"
    vt_dim1 = *ldu;
#line 309 "slalsa.f"
    vt_offset = 1 + vt_dim1;
#line 309 "slalsa.f"
    vt -= vt_offset;
#line 309 "slalsa.f"
    u_dim1 = *ldu;
#line 309 "slalsa.f"
    u_offset = 1 + u_dim1;
#line 309 "slalsa.f"
    u -= u_offset;
#line 309 "slalsa.f"
    --k;
#line 309 "slalsa.f"
    --givptr;
#line 309 "slalsa.f"
    perm_dim1 = *ldgcol;
#line 309 "slalsa.f"
    perm_offset = 1 + perm_dim1;
#line 309 "slalsa.f"
    perm -= perm_offset;
#line 309 "slalsa.f"
    givcol_dim1 = *ldgcol;
#line 309 "slalsa.f"
    givcol_offset = 1 + givcol_dim1;
#line 309 "slalsa.f"
    givcol -= givcol_offset;
#line 309 "slalsa.f"
    --c__;
#line 309 "slalsa.f"
    --s;
#line 309 "slalsa.f"
    --work;
#line 309 "slalsa.f"
    --iwork;
#line 309 "slalsa.f"

#line 309 "slalsa.f"
    /* Function Body */
#line 309 "slalsa.f"
    *info = 0;

#line 311 "slalsa.f"
    if (*icompq < 0 || *icompq > 1) {
#line 312 "slalsa.f"
	*info = -1;
#line 313 "slalsa.f"
    } else if (*smlsiz < 3) {
#line 314 "slalsa.f"
	*info = -2;
#line 315 "slalsa.f"
    } else if (*n < *smlsiz) {
#line 316 "slalsa.f"
	*info = -3;
#line 317 "slalsa.f"
    } else if (*nrhs < 1) {
#line 318 "slalsa.f"
	*info = -4;
#line 319 "slalsa.f"
    } else if (*ldb < *n) {
#line 320 "slalsa.f"
	*info = -6;
#line 321 "slalsa.f"
    } else if (*ldbx < *n) {
#line 322 "slalsa.f"
	*info = -8;
#line 323 "slalsa.f"
    } else if (*ldu < *n) {
#line 324 "slalsa.f"
	*info = -10;
#line 325 "slalsa.f"
    } else if (*ldgcol < *n) {
#line 326 "slalsa.f"
	*info = -19;
#line 327 "slalsa.f"
    }
#line 328 "slalsa.f"
    if (*info != 0) {
#line 329 "slalsa.f"
	i__1 = -(*info);
#line 329 "slalsa.f"
	xerbla_("SLALSA", &i__1, (ftnlen)6);
#line 330 "slalsa.f"
	return 0;
#line 331 "slalsa.f"
    }

/*     Book-keeping and  setting up the computation tree. */

#line 335 "slalsa.f"
    inode = 1;
#line 336 "slalsa.f"
    ndiml = inode + *n;
#line 337 "slalsa.f"
    ndimr = ndiml + *n;

#line 339 "slalsa.f"
    slasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     The following code applies back the left singular vector factors. */
/*     For applying back the right singular vector factors, go to 50. */

#line 345 "slalsa.f"
    if (*icompq == 1) {
#line 346 "slalsa.f"
	goto L50;
#line 347 "slalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by SLASDQ. The corresponding left and right singular vector */
/*     matrices are in explicit form. First apply back the left */
/*     singular vector matrices. */

#line 354 "slalsa.f"
    ndb1 = (nd + 1) / 2;
#line 355 "slalsa.f"
    i__1 = nd;
#line 355 "slalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*        IC : center row of each node */
/*        NL : number of rows of left  subproblem */
/*        NR : number of rows of right subproblem */
/*        NLF: starting row of the left   subproblem */
/*        NRF: starting row of the right  subproblem */

#line 363 "slalsa.f"
	i1 = i__ - 1;
#line 364 "slalsa.f"
	ic = iwork[inode + i1];
#line 365 "slalsa.f"
	nl = iwork[ndiml + i1];
#line 366 "slalsa.f"
	nr = iwork[ndimr + i1];
#line 367 "slalsa.f"
	nlf = ic - nl;
#line 368 "slalsa.f"
	nrf = ic + 1;
#line 369 "slalsa.f"
	sgemm_("T", "N", &nl, nrhs, &nl, &c_b7, &u[nlf + u_dim1], ldu, &b[nlf 
		+ b_dim1], ldb, &c_b8, &bx[nlf + bx_dim1], ldbx, (ftnlen)1, (
		ftnlen)1);
#line 371 "slalsa.f"
	sgemm_("T", "N", &nr, nrhs, &nr, &c_b7, &u[nrf + u_dim1], ldu, &b[nrf 
		+ b_dim1], ldb, &c_b8, &bx[nrf + bx_dim1], ldbx, (ftnlen)1, (
		ftnlen)1);
#line 373 "slalsa.f"
/* L10: */
#line 373 "slalsa.f"
    }

/*     Next copy the rows of B that correspond to unchanged rows */
/*     in the bidiagonal matrix to BX. */

#line 378 "slalsa.f"
    i__1 = nd;
#line 378 "slalsa.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 379 "slalsa.f"
	ic = iwork[inode + i__ - 1];
#line 380 "slalsa.f"
	scopy_(nrhs, &b[ic + b_dim1], ldb, &bx[ic + bx_dim1], ldbx);
#line 381 "slalsa.f"
/* L20: */
#line 381 "slalsa.f"
    }

/*     Finally go through the left singular vector matrices of all */
/*     the other subproblems bottom-up on the tree. */

#line 386 "slalsa.f"
    j = pow_ii(&c__2, &nlvl);
#line 387 "slalsa.f"
    sqre = 0;

#line 389 "slalsa.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {
#line 390 "slalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        find the first node LF and last node LL on */
/*        the current level LVL */

#line 395 "slalsa.f"
	if (lvl == 1) {
#line 396 "slalsa.f"
	    lf = 1;
#line 397 "slalsa.f"
	    ll = 1;
#line 398 "slalsa.f"
	} else {
#line 399 "slalsa.f"
	    i__1 = lvl - 1;
#line 399 "slalsa.f"
	    lf = pow_ii(&c__2, &i__1);
#line 400 "slalsa.f"
	    ll = (lf << 1) - 1;
#line 401 "slalsa.f"
	}
#line 402 "slalsa.f"
	i__1 = ll;
#line 402 "slalsa.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 403 "slalsa.f"
	    im1 = i__ - 1;
#line 404 "slalsa.f"
	    ic = iwork[inode + im1];
#line 405 "slalsa.f"
	    nl = iwork[ndiml + im1];
#line 406 "slalsa.f"
	    nr = iwork[ndimr + im1];
#line 407 "slalsa.f"
	    nlf = ic - nl;
#line 408 "slalsa.f"
	    nrf = ic + 1;
#line 409 "slalsa.f"
	    --j;
#line 410 "slalsa.f"
	    slals0_(icompq, &nl, &nr, &sqre, nrhs, &bx[nlf + bx_dim1], ldbx, &
		    b[nlf + b_dim1], ldb, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &work[1], info);
#line 417 "slalsa.f"
/* L30: */
#line 417 "slalsa.f"
	}
#line 418 "slalsa.f"
/* L40: */
#line 418 "slalsa.f"
    }
#line 419 "slalsa.f"
    goto L90;

/*     ICOMPQ = 1: applying back the right singular vector factors. */

#line 423 "slalsa.f"
L50:

/*     First now go through the right singular vector matrices of all */
/*     the tree nodes top-down. */

#line 428 "slalsa.f"
    j = 0;
#line 429 "slalsa.f"
    i__1 = nlvl;
#line 429 "slalsa.f"
    for (lvl = 1; lvl <= i__1; ++lvl) {
#line 430 "slalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        Find the first node LF and last node LL on */
/*        the current level LVL. */

#line 435 "slalsa.f"
	if (lvl == 1) {
#line 436 "slalsa.f"
	    lf = 1;
#line 437 "slalsa.f"
	    ll = 1;
#line 438 "slalsa.f"
	} else {
#line 439 "slalsa.f"
	    i__2 = lvl - 1;
#line 439 "slalsa.f"
	    lf = pow_ii(&c__2, &i__2);
#line 440 "slalsa.f"
	    ll = (lf << 1) - 1;
#line 441 "slalsa.f"
	}
#line 442 "slalsa.f"
	i__2 = lf;
#line 442 "slalsa.f"
	for (i__ = ll; i__ >= i__2; --i__) {
#line 443 "slalsa.f"
	    im1 = i__ - 1;
#line 444 "slalsa.f"
	    ic = iwork[inode + im1];
#line 445 "slalsa.f"
	    nl = iwork[ndiml + im1];
#line 446 "slalsa.f"
	    nr = iwork[ndimr + im1];
#line 447 "slalsa.f"
	    nlf = ic - nl;
#line 448 "slalsa.f"
	    nrf = ic + 1;
#line 449 "slalsa.f"
	    if (i__ == ll) {
#line 450 "slalsa.f"
		sqre = 0;
#line 451 "slalsa.f"
	    } else {
#line 452 "slalsa.f"
		sqre = 1;
#line 453 "slalsa.f"
	    }
#line 454 "slalsa.f"
	    ++j;
#line 455 "slalsa.f"
	    slals0_(icompq, &nl, &nr, &sqre, nrhs, &b[nlf + b_dim1], ldb, &bx[
		    nlf + bx_dim1], ldbx, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &work[1], info);
#line 462 "slalsa.f"
/* L60: */
#line 462 "slalsa.f"
	}
#line 463 "slalsa.f"
/* L70: */
#line 463 "slalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by SLASDQ. The corresponding right singular vector */
/*     matrices are in explicit form. Apply them back. */

#line 469 "slalsa.f"
    ndb1 = (nd + 1) / 2;
#line 470 "slalsa.f"
    i__1 = nd;
#line 470 "slalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {
#line 471 "slalsa.f"
	i1 = i__ - 1;
#line 472 "slalsa.f"
	ic = iwork[inode + i1];
#line 473 "slalsa.f"
	nl = iwork[ndiml + i1];
#line 474 "slalsa.f"
	nr = iwork[ndimr + i1];
#line 475 "slalsa.f"
	nlp1 = nl + 1;
#line 476 "slalsa.f"
	if (i__ == nd) {
#line 477 "slalsa.f"
	    nrp1 = nr;
#line 478 "slalsa.f"
	} else {
#line 479 "slalsa.f"
	    nrp1 = nr + 1;
#line 480 "slalsa.f"
	}
#line 481 "slalsa.f"
	nlf = ic - nl;
#line 482 "slalsa.f"
	nrf = ic + 1;
#line 483 "slalsa.f"
	sgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b7, &vt[nlf + vt_dim1], ldu, &
		b[nlf + b_dim1], ldb, &c_b8, &bx[nlf + bx_dim1], ldbx, (
		ftnlen)1, (ftnlen)1);
#line 485 "slalsa.f"
	sgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b7, &vt[nrf + vt_dim1], ldu, &
		b[nrf + b_dim1], ldb, &c_b8, &bx[nrf + bx_dim1], ldbx, (
		ftnlen)1, (ftnlen)1);
#line 487 "slalsa.f"
/* L80: */
#line 487 "slalsa.f"
    }

#line 489 "slalsa.f"
L90:

#line 491 "slalsa.f"
    return 0;

/*     End of SLALSA */

} /* slalsa_ */

