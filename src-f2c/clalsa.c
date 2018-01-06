#line 1 "clalsa.f"
/* clalsa.f -- translated by f2c (version 20100827).
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

#line 1 "clalsa.f"
/* Table of constant values */

static doublereal c_b9 = 1.;
static doublereal c_b10 = 0.;
static integer c__2 = 2;

/* > \brief \b CLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLALSA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsa.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsa.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsa.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, */
/*                          LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, */
/*                          GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, */
/*      $                   SMLSIZ */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), */
/*      $                   K( * ), PERM( LDGCOL, * ) */
/*       REAL               C( * ), DIFL( LDU, * ), DIFR( LDU, * ), */
/*      $                   GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ), */
/*      $                   S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * ) */
/*       COMPLEX            B( LDB, * ), BX( LDBX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLALSA is an itermediate step in solving the least squares problem */
/* > by computing the SVD of the coefficient matrix in compact form (The */
/* > singular vectors are computed as products of simple orthorgonal */
/* > matrices.). */
/* > */
/* > If ICOMPQ = 0, CLALSA applies the inverse of the left singular vector */
/* > matrix of an upper bidiagonal matrix to the right hand side; and if */
/* > ICOMPQ = 1, CLALSA applies the right singular vector matrix to the */
/* > right hand side. The singular vector matrices were generated in */
/* > compact form by CLALSA. */
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
/* >          B is COMPLEX array, dimension ( LDB, NRHS ) */
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
/* >          BX is COMPLEX array, dimension ( LDBX, NRHS ) */
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
/* >         On entry, VT**H contains the right singular vector matrices of */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension at least */
/* >         MAX( (SMLSZ+1)*NRHS*3, N*(1+NRHS) + 2*NRHS ). */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int clalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx, 
	integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *
	k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	rwork, integer *iwork, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, 
	    difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, 
	    z_dim1, z_offset, b_dim1, b_offset, bx_dim1, bx_offset, i__1, 
	    i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, i1, ic, lf, nd, ll, nl, nr, im1, nlf, nrf, lvl, 
	    ndb1, nlp1, lvl2, nrp1, jcol, nlvl, sqre, jrow, jimag, jreal, 
	    inode, ndiml;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer ndimr;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), clals0_(integer *, integer *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *), xerbla_(char *, integer *, ftnlen), 
	    slasdt_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 309 "clalsa.f"
    /* Parameter adjustments */
#line 309 "clalsa.f"
    b_dim1 = *ldb;
#line 309 "clalsa.f"
    b_offset = 1 + b_dim1;
#line 309 "clalsa.f"
    b -= b_offset;
#line 309 "clalsa.f"
    bx_dim1 = *ldbx;
#line 309 "clalsa.f"
    bx_offset = 1 + bx_dim1;
#line 309 "clalsa.f"
    bx -= bx_offset;
#line 309 "clalsa.f"
    givnum_dim1 = *ldu;
#line 309 "clalsa.f"
    givnum_offset = 1 + givnum_dim1;
#line 309 "clalsa.f"
    givnum -= givnum_offset;
#line 309 "clalsa.f"
    poles_dim1 = *ldu;
#line 309 "clalsa.f"
    poles_offset = 1 + poles_dim1;
#line 309 "clalsa.f"
    poles -= poles_offset;
#line 309 "clalsa.f"
    z_dim1 = *ldu;
#line 309 "clalsa.f"
    z_offset = 1 + z_dim1;
#line 309 "clalsa.f"
    z__ -= z_offset;
#line 309 "clalsa.f"
    difr_dim1 = *ldu;
#line 309 "clalsa.f"
    difr_offset = 1 + difr_dim1;
#line 309 "clalsa.f"
    difr -= difr_offset;
#line 309 "clalsa.f"
    difl_dim1 = *ldu;
#line 309 "clalsa.f"
    difl_offset = 1 + difl_dim1;
#line 309 "clalsa.f"
    difl -= difl_offset;
#line 309 "clalsa.f"
    vt_dim1 = *ldu;
#line 309 "clalsa.f"
    vt_offset = 1 + vt_dim1;
#line 309 "clalsa.f"
    vt -= vt_offset;
#line 309 "clalsa.f"
    u_dim1 = *ldu;
#line 309 "clalsa.f"
    u_offset = 1 + u_dim1;
#line 309 "clalsa.f"
    u -= u_offset;
#line 309 "clalsa.f"
    --k;
#line 309 "clalsa.f"
    --givptr;
#line 309 "clalsa.f"
    perm_dim1 = *ldgcol;
#line 309 "clalsa.f"
    perm_offset = 1 + perm_dim1;
#line 309 "clalsa.f"
    perm -= perm_offset;
#line 309 "clalsa.f"
    givcol_dim1 = *ldgcol;
#line 309 "clalsa.f"
    givcol_offset = 1 + givcol_dim1;
#line 309 "clalsa.f"
    givcol -= givcol_offset;
#line 309 "clalsa.f"
    --c__;
#line 309 "clalsa.f"
    --s;
#line 309 "clalsa.f"
    --rwork;
#line 309 "clalsa.f"
    --iwork;
#line 309 "clalsa.f"

#line 309 "clalsa.f"
    /* Function Body */
#line 309 "clalsa.f"
    *info = 0;

#line 311 "clalsa.f"
    if (*icompq < 0 || *icompq > 1) {
#line 312 "clalsa.f"
	*info = -1;
#line 313 "clalsa.f"
    } else if (*smlsiz < 3) {
#line 314 "clalsa.f"
	*info = -2;
#line 315 "clalsa.f"
    } else if (*n < *smlsiz) {
#line 316 "clalsa.f"
	*info = -3;
#line 317 "clalsa.f"
    } else if (*nrhs < 1) {
#line 318 "clalsa.f"
	*info = -4;
#line 319 "clalsa.f"
    } else if (*ldb < *n) {
#line 320 "clalsa.f"
	*info = -6;
#line 321 "clalsa.f"
    } else if (*ldbx < *n) {
#line 322 "clalsa.f"
	*info = -8;
#line 323 "clalsa.f"
    } else if (*ldu < *n) {
#line 324 "clalsa.f"
	*info = -10;
#line 325 "clalsa.f"
    } else if (*ldgcol < *n) {
#line 326 "clalsa.f"
	*info = -19;
#line 327 "clalsa.f"
    }
#line 328 "clalsa.f"
    if (*info != 0) {
#line 329 "clalsa.f"
	i__1 = -(*info);
#line 329 "clalsa.f"
	xerbla_("CLALSA", &i__1, (ftnlen)6);
#line 330 "clalsa.f"
	return 0;
#line 331 "clalsa.f"
    }

/*     Book-keeping and  setting up the computation tree. */

#line 335 "clalsa.f"
    inode = 1;
#line 336 "clalsa.f"
    ndiml = inode + *n;
#line 337 "clalsa.f"
    ndimr = ndiml + *n;

#line 339 "clalsa.f"
    slasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     The following code applies back the left singular vector factors. */
/*     For applying back the right singular vector factors, go to 170. */

#line 345 "clalsa.f"
    if (*icompq == 1) {
#line 346 "clalsa.f"
	goto L170;
#line 347 "clalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by SLASDQ. The corresponding left and right singular vector */
/*     matrices are in explicit form. First apply back the left */
/*     singular vector matrices. */

#line 354 "clalsa.f"
    ndb1 = (nd + 1) / 2;
#line 355 "clalsa.f"
    i__1 = nd;
#line 355 "clalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*        IC : center row of each node */
/*        NL : number of rows of left  subproblem */
/*        NR : number of rows of right subproblem */
/*        NLF: starting row of the left   subproblem */
/*        NRF: starting row of the right  subproblem */

#line 363 "clalsa.f"
	i1 = i__ - 1;
#line 364 "clalsa.f"
	ic = iwork[inode + i1];
#line 365 "clalsa.f"
	nl = iwork[ndiml + i1];
#line 366 "clalsa.f"
	nr = iwork[ndimr + i1];
#line 367 "clalsa.f"
	nlf = ic - nl;
#line 368 "clalsa.f"
	nrf = ic + 1;

/*        Since B and BX are complex, the following call to SGEMM */
/*        is performed in two steps (real and imaginary parts). */

/*        CALL SGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, */
/*     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX ) */

#line 376 "clalsa.f"
	j = nl * *nrhs << 1;
#line 377 "clalsa.f"
	i__2 = *nrhs;
#line 377 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 378 "clalsa.f"
	    i__3 = nlf + nl - 1;
#line 378 "clalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 379 "clalsa.f"
		++j;
#line 380 "clalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 380 "clalsa.f"
		rwork[j] = b[i__4].r;
#line 381 "clalsa.f"
/* L10: */
#line 381 "clalsa.f"
	    }
#line 382 "clalsa.f"
/* L20: */
#line 382 "clalsa.f"
	}
#line 383 "clalsa.f"
	sgemm_("T", "N", &nl, nrhs, &nl, &c_b9, &u[nlf + u_dim1], ldu, &rwork[
		(nl * *nrhs << 1) + 1], &nl, &c_b10, &rwork[1], &nl, (ftnlen)
		1, (ftnlen)1);
#line 385 "clalsa.f"
	j = nl * *nrhs << 1;
#line 386 "clalsa.f"
	i__2 = *nrhs;
#line 386 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 387 "clalsa.f"
	    i__3 = nlf + nl - 1;
#line 387 "clalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 388 "clalsa.f"
		++j;
#line 389 "clalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 390 "clalsa.f"
/* L30: */
#line 390 "clalsa.f"
	    }
#line 391 "clalsa.f"
/* L40: */
#line 391 "clalsa.f"
	}
#line 392 "clalsa.f"
	sgemm_("T", "N", &nl, nrhs, &nl, &c_b9, &u[nlf + u_dim1], ldu, &rwork[
		(nl * *nrhs << 1) + 1], &nl, &c_b10, &rwork[nl * *nrhs + 1], &
		nl, (ftnlen)1, (ftnlen)1);
#line 395 "clalsa.f"
	jreal = 0;
#line 396 "clalsa.f"
	jimag = nl * *nrhs;
#line 397 "clalsa.f"
	i__2 = *nrhs;
#line 397 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 398 "clalsa.f"
	    i__3 = nlf + nl - 1;
#line 398 "clalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 399 "clalsa.f"
		++jreal;
#line 400 "clalsa.f"
		++jimag;
#line 401 "clalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 401 "clalsa.f"
		i__5 = jreal;
#line 401 "clalsa.f"
		i__6 = jimag;
#line 401 "clalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 401 "clalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 403 "clalsa.f"
/* L50: */
#line 403 "clalsa.f"
	    }
#line 404 "clalsa.f"
/* L60: */
#line 404 "clalsa.f"
	}

/*        Since B and BX are complex, the following call to SGEMM */
/*        is performed in two steps (real and imaginary parts). */

/*        CALL SGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, */
/*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX ) */

#line 412 "clalsa.f"
	j = nr * *nrhs << 1;
#line 413 "clalsa.f"
	i__2 = *nrhs;
#line 413 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 414 "clalsa.f"
	    i__3 = nrf + nr - 1;
#line 414 "clalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 415 "clalsa.f"
		++j;
#line 416 "clalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 416 "clalsa.f"
		rwork[j] = b[i__4].r;
#line 417 "clalsa.f"
/* L70: */
#line 417 "clalsa.f"
	    }
#line 418 "clalsa.f"
/* L80: */
#line 418 "clalsa.f"
	}
#line 419 "clalsa.f"
	sgemm_("T", "N", &nr, nrhs, &nr, &c_b9, &u[nrf + u_dim1], ldu, &rwork[
		(nr * *nrhs << 1) + 1], &nr, &c_b10, &rwork[1], &nr, (ftnlen)
		1, (ftnlen)1);
#line 421 "clalsa.f"
	j = nr * *nrhs << 1;
#line 422 "clalsa.f"
	i__2 = *nrhs;
#line 422 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 423 "clalsa.f"
	    i__3 = nrf + nr - 1;
#line 423 "clalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 424 "clalsa.f"
		++j;
#line 425 "clalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 426 "clalsa.f"
/* L90: */
#line 426 "clalsa.f"
	    }
#line 427 "clalsa.f"
/* L100: */
#line 427 "clalsa.f"
	}
#line 428 "clalsa.f"
	sgemm_("T", "N", &nr, nrhs, &nr, &c_b9, &u[nrf + u_dim1], ldu, &rwork[
		(nr * *nrhs << 1) + 1], &nr, &c_b10, &rwork[nr * *nrhs + 1], &
		nr, (ftnlen)1, (ftnlen)1);
#line 431 "clalsa.f"
	jreal = 0;
#line 432 "clalsa.f"
	jimag = nr * *nrhs;
#line 433 "clalsa.f"
	i__2 = *nrhs;
#line 433 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 434 "clalsa.f"
	    i__3 = nrf + nr - 1;
#line 434 "clalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 435 "clalsa.f"
		++jreal;
#line 436 "clalsa.f"
		++jimag;
#line 437 "clalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 437 "clalsa.f"
		i__5 = jreal;
#line 437 "clalsa.f"
		i__6 = jimag;
#line 437 "clalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 437 "clalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 439 "clalsa.f"
/* L110: */
#line 439 "clalsa.f"
	    }
#line 440 "clalsa.f"
/* L120: */
#line 440 "clalsa.f"
	}

#line 442 "clalsa.f"
/* L130: */
#line 442 "clalsa.f"
    }

/*     Next copy the rows of B that correspond to unchanged rows */
/*     in the bidiagonal matrix to BX. */

#line 447 "clalsa.f"
    i__1 = nd;
#line 447 "clalsa.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 448 "clalsa.f"
	ic = iwork[inode + i__ - 1];
#line 449 "clalsa.f"
	ccopy_(nrhs, &b[ic + b_dim1], ldb, &bx[ic + bx_dim1], ldbx);
#line 450 "clalsa.f"
/* L140: */
#line 450 "clalsa.f"
    }

/*     Finally go through the left singular vector matrices of all */
/*     the other subproblems bottom-up on the tree. */

#line 455 "clalsa.f"
    j = pow_ii(&c__2, &nlvl);
#line 456 "clalsa.f"
    sqre = 0;

#line 458 "clalsa.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {
#line 459 "clalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        find the first node LF and last node LL on */
/*        the current level LVL */

#line 464 "clalsa.f"
	if (lvl == 1) {
#line 465 "clalsa.f"
	    lf = 1;
#line 466 "clalsa.f"
	    ll = 1;
#line 467 "clalsa.f"
	} else {
#line 468 "clalsa.f"
	    i__1 = lvl - 1;
#line 468 "clalsa.f"
	    lf = pow_ii(&c__2, &i__1);
#line 469 "clalsa.f"
	    ll = (lf << 1) - 1;
#line 470 "clalsa.f"
	}
#line 471 "clalsa.f"
	i__1 = ll;
#line 471 "clalsa.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 472 "clalsa.f"
	    im1 = i__ - 1;
#line 473 "clalsa.f"
	    ic = iwork[inode + im1];
#line 474 "clalsa.f"
	    nl = iwork[ndiml + im1];
#line 475 "clalsa.f"
	    nr = iwork[ndimr + im1];
#line 476 "clalsa.f"
	    nlf = ic - nl;
#line 477 "clalsa.f"
	    nrf = ic + 1;
#line 478 "clalsa.f"
	    --j;
#line 479 "clalsa.f"
	    clals0_(icompq, &nl, &nr, &sqre, nrhs, &bx[nlf + bx_dim1], ldbx, &
		    b[nlf + b_dim1], ldb, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &rwork[1], info);
#line 486 "clalsa.f"
/* L150: */
#line 486 "clalsa.f"
	}
#line 487 "clalsa.f"
/* L160: */
#line 487 "clalsa.f"
    }
#line 488 "clalsa.f"
    goto L330;

/*     ICOMPQ = 1: applying back the right singular vector factors. */

#line 492 "clalsa.f"
L170:

/*     First now go through the right singular vector matrices of all */
/*     the tree nodes top-down. */

#line 497 "clalsa.f"
    j = 0;
#line 498 "clalsa.f"
    i__1 = nlvl;
#line 498 "clalsa.f"
    for (lvl = 1; lvl <= i__1; ++lvl) {
#line 499 "clalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        Find the first node LF and last node LL on */
/*        the current level LVL. */

#line 504 "clalsa.f"
	if (lvl == 1) {
#line 505 "clalsa.f"
	    lf = 1;
#line 506 "clalsa.f"
	    ll = 1;
#line 507 "clalsa.f"
	} else {
#line 508 "clalsa.f"
	    i__2 = lvl - 1;
#line 508 "clalsa.f"
	    lf = pow_ii(&c__2, &i__2);
#line 509 "clalsa.f"
	    ll = (lf << 1) - 1;
#line 510 "clalsa.f"
	}
#line 511 "clalsa.f"
	i__2 = lf;
#line 511 "clalsa.f"
	for (i__ = ll; i__ >= i__2; --i__) {
#line 512 "clalsa.f"
	    im1 = i__ - 1;
#line 513 "clalsa.f"
	    ic = iwork[inode + im1];
#line 514 "clalsa.f"
	    nl = iwork[ndiml + im1];
#line 515 "clalsa.f"
	    nr = iwork[ndimr + im1];
#line 516 "clalsa.f"
	    nlf = ic - nl;
#line 517 "clalsa.f"
	    nrf = ic + 1;
#line 518 "clalsa.f"
	    if (i__ == ll) {
#line 519 "clalsa.f"
		sqre = 0;
#line 520 "clalsa.f"
	    } else {
#line 521 "clalsa.f"
		sqre = 1;
#line 522 "clalsa.f"
	    }
#line 523 "clalsa.f"
	    ++j;
#line 524 "clalsa.f"
	    clals0_(icompq, &nl, &nr, &sqre, nrhs, &b[nlf + b_dim1], ldb, &bx[
		    nlf + bx_dim1], ldbx, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &rwork[1], info);
#line 531 "clalsa.f"
/* L180: */
#line 531 "clalsa.f"
	}
#line 532 "clalsa.f"
/* L190: */
#line 532 "clalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by SLASDQ. The corresponding right singular vector */
/*     matrices are in explicit form. Apply them back. */

#line 538 "clalsa.f"
    ndb1 = (nd + 1) / 2;
#line 539 "clalsa.f"
    i__1 = nd;
#line 539 "clalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {
#line 540 "clalsa.f"
	i1 = i__ - 1;
#line 541 "clalsa.f"
	ic = iwork[inode + i1];
#line 542 "clalsa.f"
	nl = iwork[ndiml + i1];
#line 543 "clalsa.f"
	nr = iwork[ndimr + i1];
#line 544 "clalsa.f"
	nlp1 = nl + 1;
#line 545 "clalsa.f"
	if (i__ == nd) {
#line 546 "clalsa.f"
	    nrp1 = nr;
#line 547 "clalsa.f"
	} else {
#line 548 "clalsa.f"
	    nrp1 = nr + 1;
#line 549 "clalsa.f"
	}
#line 550 "clalsa.f"
	nlf = ic - nl;
#line 551 "clalsa.f"
	nrf = ic + 1;

/*        Since B and BX are complex, the following call to SGEMM is */
/*        performed in two steps (real and imaginary parts). */

/*        CALL SGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, */
/*    $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX ) */

#line 559 "clalsa.f"
	j = nlp1 * *nrhs << 1;
#line 560 "clalsa.f"
	i__2 = *nrhs;
#line 560 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 561 "clalsa.f"
	    i__3 = nlf + nlp1 - 1;
#line 561 "clalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 562 "clalsa.f"
		++j;
#line 563 "clalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 563 "clalsa.f"
		rwork[j] = b[i__4].r;
#line 564 "clalsa.f"
/* L200: */
#line 564 "clalsa.f"
	    }
#line 565 "clalsa.f"
/* L210: */
#line 565 "clalsa.f"
	}
#line 566 "clalsa.f"
	sgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b9, &vt[nlf + vt_dim1], ldu, &
		rwork[(nlp1 * *nrhs << 1) + 1], &nlp1, &c_b10, &rwork[1], &
		nlp1, (ftnlen)1, (ftnlen)1);
#line 569 "clalsa.f"
	j = nlp1 * *nrhs << 1;
#line 570 "clalsa.f"
	i__2 = *nrhs;
#line 570 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 571 "clalsa.f"
	    i__3 = nlf + nlp1 - 1;
#line 571 "clalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 572 "clalsa.f"
		++j;
#line 573 "clalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 574 "clalsa.f"
/* L220: */
#line 574 "clalsa.f"
	    }
#line 575 "clalsa.f"
/* L230: */
#line 575 "clalsa.f"
	}
#line 576 "clalsa.f"
	sgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b9, &vt[nlf + vt_dim1], ldu, &
		rwork[(nlp1 * *nrhs << 1) + 1], &nlp1, &c_b10, &rwork[nlp1 * *
		nrhs + 1], &nlp1, (ftnlen)1, (ftnlen)1);
#line 579 "clalsa.f"
	jreal = 0;
#line 580 "clalsa.f"
	jimag = nlp1 * *nrhs;
#line 581 "clalsa.f"
	i__2 = *nrhs;
#line 581 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 582 "clalsa.f"
	    i__3 = nlf + nlp1 - 1;
#line 582 "clalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 583 "clalsa.f"
		++jreal;
#line 584 "clalsa.f"
		++jimag;
#line 585 "clalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 585 "clalsa.f"
		i__5 = jreal;
#line 585 "clalsa.f"
		i__6 = jimag;
#line 585 "clalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 585 "clalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 587 "clalsa.f"
/* L240: */
#line 587 "clalsa.f"
	    }
#line 588 "clalsa.f"
/* L250: */
#line 588 "clalsa.f"
	}

/*        Since B and BX are complex, the following call to SGEMM is */
/*        performed in two steps (real and imaginary parts). */

/*        CALL SGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, */
/*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX ) */

#line 596 "clalsa.f"
	j = nrp1 * *nrhs << 1;
#line 597 "clalsa.f"
	i__2 = *nrhs;
#line 597 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 598 "clalsa.f"
	    i__3 = nrf + nrp1 - 1;
#line 598 "clalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 599 "clalsa.f"
		++j;
#line 600 "clalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 600 "clalsa.f"
		rwork[j] = b[i__4].r;
#line 601 "clalsa.f"
/* L260: */
#line 601 "clalsa.f"
	    }
#line 602 "clalsa.f"
/* L270: */
#line 602 "clalsa.f"
	}
#line 603 "clalsa.f"
	sgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b9, &vt[nrf + vt_dim1], ldu, &
		rwork[(nrp1 * *nrhs << 1) + 1], &nrp1, &c_b10, &rwork[1], &
		nrp1, (ftnlen)1, (ftnlen)1);
#line 606 "clalsa.f"
	j = nrp1 * *nrhs << 1;
#line 607 "clalsa.f"
	i__2 = *nrhs;
#line 607 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 608 "clalsa.f"
	    i__3 = nrf + nrp1 - 1;
#line 608 "clalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 609 "clalsa.f"
		++j;
#line 610 "clalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 611 "clalsa.f"
/* L280: */
#line 611 "clalsa.f"
	    }
#line 612 "clalsa.f"
/* L290: */
#line 612 "clalsa.f"
	}
#line 613 "clalsa.f"
	sgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b9, &vt[nrf + vt_dim1], ldu, &
		rwork[(nrp1 * *nrhs << 1) + 1], &nrp1, &c_b10, &rwork[nrp1 * *
		nrhs + 1], &nrp1, (ftnlen)1, (ftnlen)1);
#line 616 "clalsa.f"
	jreal = 0;
#line 617 "clalsa.f"
	jimag = nrp1 * *nrhs;
#line 618 "clalsa.f"
	i__2 = *nrhs;
#line 618 "clalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 619 "clalsa.f"
	    i__3 = nrf + nrp1 - 1;
#line 619 "clalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 620 "clalsa.f"
		++jreal;
#line 621 "clalsa.f"
		++jimag;
#line 622 "clalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 622 "clalsa.f"
		i__5 = jreal;
#line 622 "clalsa.f"
		i__6 = jimag;
#line 622 "clalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 622 "clalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 624 "clalsa.f"
/* L300: */
#line 624 "clalsa.f"
	    }
#line 625 "clalsa.f"
/* L310: */
#line 625 "clalsa.f"
	}

#line 627 "clalsa.f"
/* L320: */
#line 627 "clalsa.f"
    }

#line 629 "clalsa.f"
L330:

#line 631 "clalsa.f"
    return 0;

/*     End of CLALSA */

} /* clalsa_ */

