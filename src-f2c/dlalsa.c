#line 1 "dlalsa.f"
/* dlalsa.f -- translated by f2c (version 20100827).
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

#line 1 "dlalsa.f"
/* Table of constant values */

static doublereal c_b7 = 1.;
static doublereal c_b8 = 0.;
static integer c__2 = 2;

/* > \brief \b DLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLALSA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlalsa.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlalsa.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlalsa.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, */
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
/*       DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), C( * ), */
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
/* > DLALSA is an itermediate step in solving the least squares problem */
/* > by computing the SVD of the coefficient matrix in compact form (The */
/* > singular vectors are computed as products of simple orthorgonal */
/* > matrices.). */
/* > */
/* > If ICOMPQ = 0, DLALSA applies the inverse of the left singular vector */
/* > matrix of an upper bidiagonal matrix to the right hand side; and if */
/* > ICOMPQ = 1, DLALSA applies the right singular vector matrix to the */
/* > right hand side. The singular vector matrices were generated in */
/* > compact form by DLALSA. */
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
/* >          B is DOUBLE PRECISION array, dimension ( LDB, NRHS ) */
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
/* >          BX is DOUBLE PRECISION array, dimension ( LDBX, NRHS ) */
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
/* >          U is DOUBLE PRECISION array, dimension ( LDU, SMLSIZ ). */
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
/* >          VT is DOUBLE PRECISION array, dimension ( LDU, SMLSIZ+1 ). */
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
/* >          DIFL is DOUBLE PRECISION array, dimension ( LDU, NLVL ). */
/* >         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFR */
/* > \verbatim */
/* >          DIFR is DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ). */
/* >         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record */
/* >         distances between singular values on the I-th level and */
/* >         singular values on the (I -1)-th level, and DIFR(*, 2 * I) */
/* >         record the normalizing factors of the right singular vectors */
/* >         matrices of subproblems on I-th level. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( LDU, NLVL ). */
/* >         On entry, Z(1, I) contains the components of the deflation- */
/* >         adjusted updating row vector for subproblems on the I-th */
/* >         level. */
/* > \endverbatim */
/* > */
/* > \param[in] POLES */
/* > \verbatim */
/* >          POLES is DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ). */
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
/* >          GIVNUM is DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ). */
/* >         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S- */
/* >         values of Givens rotations performed on the I-th level on the */
/* >         computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension ( N ). */
/* >         On entry, if the I-th subproblem is not square, */
/* >         C( I ) contains the C-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension ( N ). */
/* >         On entry, if the I-th subproblem is not square, */
/* >         S( I ) contains the S-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (3*N) */
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

/* > \date June 2017 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int dlalsa_(integer *icompq, integer *smlsiz, integer *n, 
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
	    ndb1, nlp1, lvl2, nrp1, nlvl, sqre;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer inode, ndiml, ndimr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlals0_(integer *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *), dlasdt_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 307 "dlalsa.f"
    /* Parameter adjustments */
#line 307 "dlalsa.f"
    b_dim1 = *ldb;
#line 307 "dlalsa.f"
    b_offset = 1 + b_dim1;
#line 307 "dlalsa.f"
    b -= b_offset;
#line 307 "dlalsa.f"
    bx_dim1 = *ldbx;
#line 307 "dlalsa.f"
    bx_offset = 1 + bx_dim1;
#line 307 "dlalsa.f"
    bx -= bx_offset;
#line 307 "dlalsa.f"
    givnum_dim1 = *ldu;
#line 307 "dlalsa.f"
    givnum_offset = 1 + givnum_dim1;
#line 307 "dlalsa.f"
    givnum -= givnum_offset;
#line 307 "dlalsa.f"
    poles_dim1 = *ldu;
#line 307 "dlalsa.f"
    poles_offset = 1 + poles_dim1;
#line 307 "dlalsa.f"
    poles -= poles_offset;
#line 307 "dlalsa.f"
    z_dim1 = *ldu;
#line 307 "dlalsa.f"
    z_offset = 1 + z_dim1;
#line 307 "dlalsa.f"
    z__ -= z_offset;
#line 307 "dlalsa.f"
    difr_dim1 = *ldu;
#line 307 "dlalsa.f"
    difr_offset = 1 + difr_dim1;
#line 307 "dlalsa.f"
    difr -= difr_offset;
#line 307 "dlalsa.f"
    difl_dim1 = *ldu;
#line 307 "dlalsa.f"
    difl_offset = 1 + difl_dim1;
#line 307 "dlalsa.f"
    difl -= difl_offset;
#line 307 "dlalsa.f"
    vt_dim1 = *ldu;
#line 307 "dlalsa.f"
    vt_offset = 1 + vt_dim1;
#line 307 "dlalsa.f"
    vt -= vt_offset;
#line 307 "dlalsa.f"
    u_dim1 = *ldu;
#line 307 "dlalsa.f"
    u_offset = 1 + u_dim1;
#line 307 "dlalsa.f"
    u -= u_offset;
#line 307 "dlalsa.f"
    --k;
#line 307 "dlalsa.f"
    --givptr;
#line 307 "dlalsa.f"
    perm_dim1 = *ldgcol;
#line 307 "dlalsa.f"
    perm_offset = 1 + perm_dim1;
#line 307 "dlalsa.f"
    perm -= perm_offset;
#line 307 "dlalsa.f"
    givcol_dim1 = *ldgcol;
#line 307 "dlalsa.f"
    givcol_offset = 1 + givcol_dim1;
#line 307 "dlalsa.f"
    givcol -= givcol_offset;
#line 307 "dlalsa.f"
    --c__;
#line 307 "dlalsa.f"
    --s;
#line 307 "dlalsa.f"
    --work;
#line 307 "dlalsa.f"
    --iwork;
#line 307 "dlalsa.f"

#line 307 "dlalsa.f"
    /* Function Body */
#line 307 "dlalsa.f"
    *info = 0;

#line 309 "dlalsa.f"
    if (*icompq < 0 || *icompq > 1) {
#line 310 "dlalsa.f"
	*info = -1;
#line 311 "dlalsa.f"
    } else if (*smlsiz < 3) {
#line 312 "dlalsa.f"
	*info = -2;
#line 313 "dlalsa.f"
    } else if (*n < *smlsiz) {
#line 314 "dlalsa.f"
	*info = -3;
#line 315 "dlalsa.f"
    } else if (*nrhs < 1) {
#line 316 "dlalsa.f"
	*info = -4;
#line 317 "dlalsa.f"
    } else if (*ldb < *n) {
#line 318 "dlalsa.f"
	*info = -6;
#line 319 "dlalsa.f"
    } else if (*ldbx < *n) {
#line 320 "dlalsa.f"
	*info = -8;
#line 321 "dlalsa.f"
    } else if (*ldu < *n) {
#line 322 "dlalsa.f"
	*info = -10;
#line 323 "dlalsa.f"
    } else if (*ldgcol < *n) {
#line 324 "dlalsa.f"
	*info = -19;
#line 325 "dlalsa.f"
    }
#line 326 "dlalsa.f"
    if (*info != 0) {
#line 327 "dlalsa.f"
	i__1 = -(*info);
#line 327 "dlalsa.f"
	xerbla_("DLALSA", &i__1, (ftnlen)6);
#line 328 "dlalsa.f"
	return 0;
#line 329 "dlalsa.f"
    }

/*     Book-keeping and  setting up the computation tree. */

#line 333 "dlalsa.f"
    inode = 1;
#line 334 "dlalsa.f"
    ndiml = inode + *n;
#line 335 "dlalsa.f"
    ndimr = ndiml + *n;

#line 337 "dlalsa.f"
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     The following code applies back the left singular vector factors. */
/*     For applying back the right singular vector factors, go to 50. */

#line 343 "dlalsa.f"
    if (*icompq == 1) {
#line 344 "dlalsa.f"
	goto L50;
#line 345 "dlalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by DLASDQ. The corresponding left and right singular vector */
/*     matrices are in explicit form. First apply back the left */
/*     singular vector matrices. */

#line 352 "dlalsa.f"
    ndb1 = (nd + 1) / 2;
#line 353 "dlalsa.f"
    i__1 = nd;
#line 353 "dlalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*        IC : center row of each node */
/*        NL : number of rows of left  subproblem */
/*        NR : number of rows of right subproblem */
/*        NLF: starting row of the left   subproblem */
/*        NRF: starting row of the right  subproblem */

#line 361 "dlalsa.f"
	i1 = i__ - 1;
#line 362 "dlalsa.f"
	ic = iwork[inode + i1];
#line 363 "dlalsa.f"
	nl = iwork[ndiml + i1];
#line 364 "dlalsa.f"
	nr = iwork[ndimr + i1];
#line 365 "dlalsa.f"
	nlf = ic - nl;
#line 366 "dlalsa.f"
	nrf = ic + 1;
#line 367 "dlalsa.f"
	dgemm_("T", "N", &nl, nrhs, &nl, &c_b7, &u[nlf + u_dim1], ldu, &b[nlf 
		+ b_dim1], ldb, &c_b8, &bx[nlf + bx_dim1], ldbx, (ftnlen)1, (
		ftnlen)1);
#line 369 "dlalsa.f"
	dgemm_("T", "N", &nr, nrhs, &nr, &c_b7, &u[nrf + u_dim1], ldu, &b[nrf 
		+ b_dim1], ldb, &c_b8, &bx[nrf + bx_dim1], ldbx, (ftnlen)1, (
		ftnlen)1);
#line 371 "dlalsa.f"
/* L10: */
#line 371 "dlalsa.f"
    }

/*     Next copy the rows of B that correspond to unchanged rows */
/*     in the bidiagonal matrix to BX. */

#line 376 "dlalsa.f"
    i__1 = nd;
#line 376 "dlalsa.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "dlalsa.f"
	ic = iwork[inode + i__ - 1];
#line 378 "dlalsa.f"
	dcopy_(nrhs, &b[ic + b_dim1], ldb, &bx[ic + bx_dim1], ldbx);
#line 379 "dlalsa.f"
/* L20: */
#line 379 "dlalsa.f"
    }

/*     Finally go through the left singular vector matrices of all */
/*     the other subproblems bottom-up on the tree. */

#line 384 "dlalsa.f"
    j = pow_ii(&c__2, &nlvl);
#line 385 "dlalsa.f"
    sqre = 0;

#line 387 "dlalsa.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {
#line 388 "dlalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        find the first node LF and last node LL on */
/*        the current level LVL */

#line 393 "dlalsa.f"
	if (lvl == 1) {
#line 394 "dlalsa.f"
	    lf = 1;
#line 395 "dlalsa.f"
	    ll = 1;
#line 396 "dlalsa.f"
	} else {
#line 397 "dlalsa.f"
	    i__1 = lvl - 1;
#line 397 "dlalsa.f"
	    lf = pow_ii(&c__2, &i__1);
#line 398 "dlalsa.f"
	    ll = (lf << 1) - 1;
#line 399 "dlalsa.f"
	}
#line 400 "dlalsa.f"
	i__1 = ll;
#line 400 "dlalsa.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 401 "dlalsa.f"
	    im1 = i__ - 1;
#line 402 "dlalsa.f"
	    ic = iwork[inode + im1];
#line 403 "dlalsa.f"
	    nl = iwork[ndiml + im1];
#line 404 "dlalsa.f"
	    nr = iwork[ndimr + im1];
#line 405 "dlalsa.f"
	    nlf = ic - nl;
#line 406 "dlalsa.f"
	    nrf = ic + 1;
#line 407 "dlalsa.f"
	    --j;
#line 408 "dlalsa.f"
	    dlals0_(icompq, &nl, &nr, &sqre, nrhs, &bx[nlf + bx_dim1], ldbx, &
		    b[nlf + b_dim1], ldb, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &work[1], info);
#line 415 "dlalsa.f"
/* L30: */
#line 415 "dlalsa.f"
	}
#line 416 "dlalsa.f"
/* L40: */
#line 416 "dlalsa.f"
    }
#line 417 "dlalsa.f"
    goto L90;

/*     ICOMPQ = 1: applying back the right singular vector factors. */

#line 421 "dlalsa.f"
L50:

/*     First now go through the right singular vector matrices of all */
/*     the tree nodes top-down. */

#line 426 "dlalsa.f"
    j = 0;
#line 427 "dlalsa.f"
    i__1 = nlvl;
#line 427 "dlalsa.f"
    for (lvl = 1; lvl <= i__1; ++lvl) {
#line 428 "dlalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        Find the first node LF and last node LL on */
/*        the current level LVL. */

#line 433 "dlalsa.f"
	if (lvl == 1) {
#line 434 "dlalsa.f"
	    lf = 1;
#line 435 "dlalsa.f"
	    ll = 1;
#line 436 "dlalsa.f"
	} else {
#line 437 "dlalsa.f"
	    i__2 = lvl - 1;
#line 437 "dlalsa.f"
	    lf = pow_ii(&c__2, &i__2);
#line 438 "dlalsa.f"
	    ll = (lf << 1) - 1;
#line 439 "dlalsa.f"
	}
#line 440 "dlalsa.f"
	i__2 = lf;
#line 440 "dlalsa.f"
	for (i__ = ll; i__ >= i__2; --i__) {
#line 441 "dlalsa.f"
	    im1 = i__ - 1;
#line 442 "dlalsa.f"
	    ic = iwork[inode + im1];
#line 443 "dlalsa.f"
	    nl = iwork[ndiml + im1];
#line 444 "dlalsa.f"
	    nr = iwork[ndimr + im1];
#line 445 "dlalsa.f"
	    nlf = ic - nl;
#line 446 "dlalsa.f"
	    nrf = ic + 1;
#line 447 "dlalsa.f"
	    if (i__ == ll) {
#line 448 "dlalsa.f"
		sqre = 0;
#line 449 "dlalsa.f"
	    } else {
#line 450 "dlalsa.f"
		sqre = 1;
#line 451 "dlalsa.f"
	    }
#line 452 "dlalsa.f"
	    ++j;
#line 453 "dlalsa.f"
	    dlals0_(icompq, &nl, &nr, &sqre, nrhs, &b[nlf + b_dim1], ldb, &bx[
		    nlf + bx_dim1], ldbx, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &work[1], info);
#line 460 "dlalsa.f"
/* L60: */
#line 460 "dlalsa.f"
	}
#line 461 "dlalsa.f"
/* L70: */
#line 461 "dlalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by DLASDQ. The corresponding right singular vector */
/*     matrices are in explicit form. Apply them back. */

#line 467 "dlalsa.f"
    ndb1 = (nd + 1) / 2;
#line 468 "dlalsa.f"
    i__1 = nd;
#line 468 "dlalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {
#line 469 "dlalsa.f"
	i1 = i__ - 1;
#line 470 "dlalsa.f"
	ic = iwork[inode + i1];
#line 471 "dlalsa.f"
	nl = iwork[ndiml + i1];
#line 472 "dlalsa.f"
	nr = iwork[ndimr + i1];
#line 473 "dlalsa.f"
	nlp1 = nl + 1;
#line 474 "dlalsa.f"
	if (i__ == nd) {
#line 475 "dlalsa.f"
	    nrp1 = nr;
#line 476 "dlalsa.f"
	} else {
#line 477 "dlalsa.f"
	    nrp1 = nr + 1;
#line 478 "dlalsa.f"
	}
#line 479 "dlalsa.f"
	nlf = ic - nl;
#line 480 "dlalsa.f"
	nrf = ic + 1;
#line 481 "dlalsa.f"
	dgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b7, &vt[nlf + vt_dim1], ldu, &
		b[nlf + b_dim1], ldb, &c_b8, &bx[nlf + bx_dim1], ldbx, (
		ftnlen)1, (ftnlen)1);
#line 483 "dlalsa.f"
	dgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b7, &vt[nrf + vt_dim1], ldu, &
		b[nrf + b_dim1], ldb, &c_b8, &bx[nrf + bx_dim1], ldbx, (
		ftnlen)1, (ftnlen)1);
#line 485 "dlalsa.f"
/* L80: */
#line 485 "dlalsa.f"
    }

#line 487 "dlalsa.f"
L90:

#line 489 "dlalsa.f"
    return 0;

/*     End of DLALSA */

} /* dlalsa_ */

