#line 1 "zlalsa.f"
/* zlalsa.f -- translated by f2c (version 20100827).
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

#line 1 "zlalsa.f"
/* Table of constant values */

static doublereal c_b9 = 1.;
static doublereal c_b10 = 0.;
static integer c__2 = 2;

/* > \brief \b ZLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLALSA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlalsa.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlalsa.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlalsa.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, */
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
/*       DOUBLE PRECISION   C( * ), DIFL( LDU, * ), DIFR( LDU, * ), */
/*      $                   GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ), */
/*      $                   S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * ) */
/*       COMPLEX*16         B( LDB, * ), BX( LDBX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLALSA is an itermediate step in solving the least squares problem */
/* > by computing the SVD of the coefficient matrix in compact form (The */
/* > singular vectors are computed as products of simple orthorgonal */
/* > matrices.). */
/* > */
/* > If ICOMPQ = 0, ZLALSA applies the inverse of the left singular vector */
/* > matrix of an upper bidiagonal matrix to the right hand side; and if */
/* > ICOMPQ = 1, ZLALSA applies the right singular vector matrix to the */
/* > right hand side. The singular vector matrices were generated in */
/* > compact form by ZLALSA. */
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
/* >          B is COMPLEX*16 array, dimension ( LDB, NRHS ) */
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
/* >          BX is COMPLEX*16 array, dimension ( LDBX, NRHS ) */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension at least */
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

/* > \ingroup complex16OTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int zlalsa_(integer *icompq, integer *smlsiz, integer *n, 
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
	    ndb1, nlp1, lvl2, nrp1, jcol, nlvl, sqre, jrow, jimag;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer jreal, inode, ndiml, ndimr;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlals0_(integer *, integer *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *), dlasdt_(integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);


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

#line 309 "zlalsa.f"
    /* Parameter adjustments */
#line 309 "zlalsa.f"
    b_dim1 = *ldb;
#line 309 "zlalsa.f"
    b_offset = 1 + b_dim1;
#line 309 "zlalsa.f"
    b -= b_offset;
#line 309 "zlalsa.f"
    bx_dim1 = *ldbx;
#line 309 "zlalsa.f"
    bx_offset = 1 + bx_dim1;
#line 309 "zlalsa.f"
    bx -= bx_offset;
#line 309 "zlalsa.f"
    givnum_dim1 = *ldu;
#line 309 "zlalsa.f"
    givnum_offset = 1 + givnum_dim1;
#line 309 "zlalsa.f"
    givnum -= givnum_offset;
#line 309 "zlalsa.f"
    poles_dim1 = *ldu;
#line 309 "zlalsa.f"
    poles_offset = 1 + poles_dim1;
#line 309 "zlalsa.f"
    poles -= poles_offset;
#line 309 "zlalsa.f"
    z_dim1 = *ldu;
#line 309 "zlalsa.f"
    z_offset = 1 + z_dim1;
#line 309 "zlalsa.f"
    z__ -= z_offset;
#line 309 "zlalsa.f"
    difr_dim1 = *ldu;
#line 309 "zlalsa.f"
    difr_offset = 1 + difr_dim1;
#line 309 "zlalsa.f"
    difr -= difr_offset;
#line 309 "zlalsa.f"
    difl_dim1 = *ldu;
#line 309 "zlalsa.f"
    difl_offset = 1 + difl_dim1;
#line 309 "zlalsa.f"
    difl -= difl_offset;
#line 309 "zlalsa.f"
    vt_dim1 = *ldu;
#line 309 "zlalsa.f"
    vt_offset = 1 + vt_dim1;
#line 309 "zlalsa.f"
    vt -= vt_offset;
#line 309 "zlalsa.f"
    u_dim1 = *ldu;
#line 309 "zlalsa.f"
    u_offset = 1 + u_dim1;
#line 309 "zlalsa.f"
    u -= u_offset;
#line 309 "zlalsa.f"
    --k;
#line 309 "zlalsa.f"
    --givptr;
#line 309 "zlalsa.f"
    perm_dim1 = *ldgcol;
#line 309 "zlalsa.f"
    perm_offset = 1 + perm_dim1;
#line 309 "zlalsa.f"
    perm -= perm_offset;
#line 309 "zlalsa.f"
    givcol_dim1 = *ldgcol;
#line 309 "zlalsa.f"
    givcol_offset = 1 + givcol_dim1;
#line 309 "zlalsa.f"
    givcol -= givcol_offset;
#line 309 "zlalsa.f"
    --c__;
#line 309 "zlalsa.f"
    --s;
#line 309 "zlalsa.f"
    --rwork;
#line 309 "zlalsa.f"
    --iwork;
#line 309 "zlalsa.f"

#line 309 "zlalsa.f"
    /* Function Body */
#line 309 "zlalsa.f"
    *info = 0;

#line 311 "zlalsa.f"
    if (*icompq < 0 || *icompq > 1) {
#line 312 "zlalsa.f"
	*info = -1;
#line 313 "zlalsa.f"
    } else if (*smlsiz < 3) {
#line 314 "zlalsa.f"
	*info = -2;
#line 315 "zlalsa.f"
    } else if (*n < *smlsiz) {
#line 316 "zlalsa.f"
	*info = -3;
#line 317 "zlalsa.f"
    } else if (*nrhs < 1) {
#line 318 "zlalsa.f"
	*info = -4;
#line 319 "zlalsa.f"
    } else if (*ldb < *n) {
#line 320 "zlalsa.f"
	*info = -6;
#line 321 "zlalsa.f"
    } else if (*ldbx < *n) {
#line 322 "zlalsa.f"
	*info = -8;
#line 323 "zlalsa.f"
    } else if (*ldu < *n) {
#line 324 "zlalsa.f"
	*info = -10;
#line 325 "zlalsa.f"
    } else if (*ldgcol < *n) {
#line 326 "zlalsa.f"
	*info = -19;
#line 327 "zlalsa.f"
    }
#line 328 "zlalsa.f"
    if (*info != 0) {
#line 329 "zlalsa.f"
	i__1 = -(*info);
#line 329 "zlalsa.f"
	xerbla_("ZLALSA", &i__1, (ftnlen)6);
#line 330 "zlalsa.f"
	return 0;
#line 331 "zlalsa.f"
    }

/*     Book-keeping and  setting up the computation tree. */

#line 335 "zlalsa.f"
    inode = 1;
#line 336 "zlalsa.f"
    ndiml = inode + *n;
#line 337 "zlalsa.f"
    ndimr = ndiml + *n;

#line 339 "zlalsa.f"
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     The following code applies back the left singular vector factors. */
/*     For applying back the right singular vector factors, go to 170. */

#line 345 "zlalsa.f"
    if (*icompq == 1) {
#line 346 "zlalsa.f"
	goto L170;
#line 347 "zlalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by DLASDQ. The corresponding left and right singular vector */
/*     matrices are in explicit form. First apply back the left */
/*     singular vector matrices. */

#line 354 "zlalsa.f"
    ndb1 = (nd + 1) / 2;
#line 355 "zlalsa.f"
    i__1 = nd;
#line 355 "zlalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*        IC : center row of each node */
/*        NL : number of rows of left  subproblem */
/*        NR : number of rows of right subproblem */
/*        NLF: starting row of the left   subproblem */
/*        NRF: starting row of the right  subproblem */

#line 363 "zlalsa.f"
	i1 = i__ - 1;
#line 364 "zlalsa.f"
	ic = iwork[inode + i1];
#line 365 "zlalsa.f"
	nl = iwork[ndiml + i1];
#line 366 "zlalsa.f"
	nr = iwork[ndimr + i1];
#line 367 "zlalsa.f"
	nlf = ic - nl;
#line 368 "zlalsa.f"
	nrf = ic + 1;

/*        Since B and BX are complex, the following call to DGEMM */
/*        is performed in two steps (real and imaginary parts). */

/*        CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, */
/*     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX ) */

#line 376 "zlalsa.f"
	j = nl * *nrhs << 1;
#line 377 "zlalsa.f"
	i__2 = *nrhs;
#line 377 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 378 "zlalsa.f"
	    i__3 = nlf + nl - 1;
#line 378 "zlalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 379 "zlalsa.f"
		++j;
#line 380 "zlalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 380 "zlalsa.f"
		rwork[j] = b[i__4].r;
#line 381 "zlalsa.f"
/* L10: */
#line 381 "zlalsa.f"
	    }
#line 382 "zlalsa.f"
/* L20: */
#line 382 "zlalsa.f"
	}
#line 383 "zlalsa.f"
	dgemm_("T", "N", &nl, nrhs, &nl, &c_b9, &u[nlf + u_dim1], ldu, &rwork[
		(nl * *nrhs << 1) + 1], &nl, &c_b10, &rwork[1], &nl, (ftnlen)
		1, (ftnlen)1);
#line 385 "zlalsa.f"
	j = nl * *nrhs << 1;
#line 386 "zlalsa.f"
	i__2 = *nrhs;
#line 386 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 387 "zlalsa.f"
	    i__3 = nlf + nl - 1;
#line 387 "zlalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 388 "zlalsa.f"
		++j;
#line 389 "zlalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 390 "zlalsa.f"
/* L30: */
#line 390 "zlalsa.f"
	    }
#line 391 "zlalsa.f"
/* L40: */
#line 391 "zlalsa.f"
	}
#line 392 "zlalsa.f"
	dgemm_("T", "N", &nl, nrhs, &nl, &c_b9, &u[nlf + u_dim1], ldu, &rwork[
		(nl * *nrhs << 1) + 1], &nl, &c_b10, &rwork[nl * *nrhs + 1], &
		nl, (ftnlen)1, (ftnlen)1);
#line 395 "zlalsa.f"
	jreal = 0;
#line 396 "zlalsa.f"
	jimag = nl * *nrhs;
#line 397 "zlalsa.f"
	i__2 = *nrhs;
#line 397 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 398 "zlalsa.f"
	    i__3 = nlf + nl - 1;
#line 398 "zlalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 399 "zlalsa.f"
		++jreal;
#line 400 "zlalsa.f"
		++jimag;
#line 401 "zlalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 401 "zlalsa.f"
		i__5 = jreal;
#line 401 "zlalsa.f"
		i__6 = jimag;
#line 401 "zlalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 401 "zlalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 403 "zlalsa.f"
/* L50: */
#line 403 "zlalsa.f"
	    }
#line 404 "zlalsa.f"
/* L60: */
#line 404 "zlalsa.f"
	}

/*        Since B and BX are complex, the following call to DGEMM */
/*        is performed in two steps (real and imaginary parts). */

/*        CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, */
/*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX ) */

#line 412 "zlalsa.f"
	j = nr * *nrhs << 1;
#line 413 "zlalsa.f"
	i__2 = *nrhs;
#line 413 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 414 "zlalsa.f"
	    i__3 = nrf + nr - 1;
#line 414 "zlalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 415 "zlalsa.f"
		++j;
#line 416 "zlalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 416 "zlalsa.f"
		rwork[j] = b[i__4].r;
#line 417 "zlalsa.f"
/* L70: */
#line 417 "zlalsa.f"
	    }
#line 418 "zlalsa.f"
/* L80: */
#line 418 "zlalsa.f"
	}
#line 419 "zlalsa.f"
	dgemm_("T", "N", &nr, nrhs, &nr, &c_b9, &u[nrf + u_dim1], ldu, &rwork[
		(nr * *nrhs << 1) + 1], &nr, &c_b10, &rwork[1], &nr, (ftnlen)
		1, (ftnlen)1);
#line 421 "zlalsa.f"
	j = nr * *nrhs << 1;
#line 422 "zlalsa.f"
	i__2 = *nrhs;
#line 422 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 423 "zlalsa.f"
	    i__3 = nrf + nr - 1;
#line 423 "zlalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 424 "zlalsa.f"
		++j;
#line 425 "zlalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 426 "zlalsa.f"
/* L90: */
#line 426 "zlalsa.f"
	    }
#line 427 "zlalsa.f"
/* L100: */
#line 427 "zlalsa.f"
	}
#line 428 "zlalsa.f"
	dgemm_("T", "N", &nr, nrhs, &nr, &c_b9, &u[nrf + u_dim1], ldu, &rwork[
		(nr * *nrhs << 1) + 1], &nr, &c_b10, &rwork[nr * *nrhs + 1], &
		nr, (ftnlen)1, (ftnlen)1);
#line 431 "zlalsa.f"
	jreal = 0;
#line 432 "zlalsa.f"
	jimag = nr * *nrhs;
#line 433 "zlalsa.f"
	i__2 = *nrhs;
#line 433 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 434 "zlalsa.f"
	    i__3 = nrf + nr - 1;
#line 434 "zlalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 435 "zlalsa.f"
		++jreal;
#line 436 "zlalsa.f"
		++jimag;
#line 437 "zlalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 437 "zlalsa.f"
		i__5 = jreal;
#line 437 "zlalsa.f"
		i__6 = jimag;
#line 437 "zlalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 437 "zlalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 439 "zlalsa.f"
/* L110: */
#line 439 "zlalsa.f"
	    }
#line 440 "zlalsa.f"
/* L120: */
#line 440 "zlalsa.f"
	}

#line 442 "zlalsa.f"
/* L130: */
#line 442 "zlalsa.f"
    }

/*     Next copy the rows of B that correspond to unchanged rows */
/*     in the bidiagonal matrix to BX. */

#line 447 "zlalsa.f"
    i__1 = nd;
#line 447 "zlalsa.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 448 "zlalsa.f"
	ic = iwork[inode + i__ - 1];
#line 449 "zlalsa.f"
	zcopy_(nrhs, &b[ic + b_dim1], ldb, &bx[ic + bx_dim1], ldbx);
#line 450 "zlalsa.f"
/* L140: */
#line 450 "zlalsa.f"
    }

/*     Finally go through the left singular vector matrices of all */
/*     the other subproblems bottom-up on the tree. */

#line 455 "zlalsa.f"
    j = pow_ii(&c__2, &nlvl);
#line 456 "zlalsa.f"
    sqre = 0;

#line 458 "zlalsa.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {
#line 459 "zlalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        find the first node LF and last node LL on */
/*        the current level LVL */

#line 464 "zlalsa.f"
	if (lvl == 1) {
#line 465 "zlalsa.f"
	    lf = 1;
#line 466 "zlalsa.f"
	    ll = 1;
#line 467 "zlalsa.f"
	} else {
#line 468 "zlalsa.f"
	    i__1 = lvl - 1;
#line 468 "zlalsa.f"
	    lf = pow_ii(&c__2, &i__1);
#line 469 "zlalsa.f"
	    ll = (lf << 1) - 1;
#line 470 "zlalsa.f"
	}
#line 471 "zlalsa.f"
	i__1 = ll;
#line 471 "zlalsa.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 472 "zlalsa.f"
	    im1 = i__ - 1;
#line 473 "zlalsa.f"
	    ic = iwork[inode + im1];
#line 474 "zlalsa.f"
	    nl = iwork[ndiml + im1];
#line 475 "zlalsa.f"
	    nr = iwork[ndimr + im1];
#line 476 "zlalsa.f"
	    nlf = ic - nl;
#line 477 "zlalsa.f"
	    nrf = ic + 1;
#line 478 "zlalsa.f"
	    --j;
#line 479 "zlalsa.f"
	    zlals0_(icompq, &nl, &nr, &sqre, nrhs, &bx[nlf + bx_dim1], ldbx, &
		    b[nlf + b_dim1], ldb, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &rwork[1], info);
#line 486 "zlalsa.f"
/* L150: */
#line 486 "zlalsa.f"
	}
#line 487 "zlalsa.f"
/* L160: */
#line 487 "zlalsa.f"
    }
#line 488 "zlalsa.f"
    goto L330;

/*     ICOMPQ = 1: applying back the right singular vector factors. */

#line 492 "zlalsa.f"
L170:

/*     First now go through the right singular vector matrices of all */
/*     the tree nodes top-down. */

#line 497 "zlalsa.f"
    j = 0;
#line 498 "zlalsa.f"
    i__1 = nlvl;
#line 498 "zlalsa.f"
    for (lvl = 1; lvl <= i__1; ++lvl) {
#line 499 "zlalsa.f"
	lvl2 = (lvl << 1) - 1;

/*        Find the first node LF and last node LL on */
/*        the current level LVL. */

#line 504 "zlalsa.f"
	if (lvl == 1) {
#line 505 "zlalsa.f"
	    lf = 1;
#line 506 "zlalsa.f"
	    ll = 1;
#line 507 "zlalsa.f"
	} else {
#line 508 "zlalsa.f"
	    i__2 = lvl - 1;
#line 508 "zlalsa.f"
	    lf = pow_ii(&c__2, &i__2);
#line 509 "zlalsa.f"
	    ll = (lf << 1) - 1;
#line 510 "zlalsa.f"
	}
#line 511 "zlalsa.f"
	i__2 = lf;
#line 511 "zlalsa.f"
	for (i__ = ll; i__ >= i__2; --i__) {
#line 512 "zlalsa.f"
	    im1 = i__ - 1;
#line 513 "zlalsa.f"
	    ic = iwork[inode + im1];
#line 514 "zlalsa.f"
	    nl = iwork[ndiml + im1];
#line 515 "zlalsa.f"
	    nr = iwork[ndimr + im1];
#line 516 "zlalsa.f"
	    nlf = ic - nl;
#line 517 "zlalsa.f"
	    nrf = ic + 1;
#line 518 "zlalsa.f"
	    if (i__ == ll) {
#line 519 "zlalsa.f"
		sqre = 0;
#line 520 "zlalsa.f"
	    } else {
#line 521 "zlalsa.f"
		sqre = 1;
#line 522 "zlalsa.f"
	    }
#line 523 "zlalsa.f"
	    ++j;
#line 524 "zlalsa.f"
	    zlals0_(icompq, &nl, &nr, &sqre, nrhs, &b[nlf + b_dim1], ldb, &bx[
		    nlf + bx_dim1], ldbx, &perm[nlf + lvl * perm_dim1], &
		    givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, &
		    givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 *
		     poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + 
		    lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[
		    j], &s[j], &rwork[1], info);
#line 531 "zlalsa.f"
/* L180: */
#line 531 "zlalsa.f"
	}
#line 532 "zlalsa.f"
/* L190: */
#line 532 "zlalsa.f"
    }

/*     The nodes on the bottom level of the tree were solved */
/*     by DLASDQ. The corresponding right singular vector */
/*     matrices are in explicit form. Apply them back. */

#line 538 "zlalsa.f"
    ndb1 = (nd + 1) / 2;
#line 539 "zlalsa.f"
    i__1 = nd;
#line 539 "zlalsa.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {
#line 540 "zlalsa.f"
	i1 = i__ - 1;
#line 541 "zlalsa.f"
	ic = iwork[inode + i1];
#line 542 "zlalsa.f"
	nl = iwork[ndiml + i1];
#line 543 "zlalsa.f"
	nr = iwork[ndimr + i1];
#line 544 "zlalsa.f"
	nlp1 = nl + 1;
#line 545 "zlalsa.f"
	if (i__ == nd) {
#line 546 "zlalsa.f"
	    nrp1 = nr;
#line 547 "zlalsa.f"
	} else {
#line 548 "zlalsa.f"
	    nrp1 = nr + 1;
#line 549 "zlalsa.f"
	}
#line 550 "zlalsa.f"
	nlf = ic - nl;
#line 551 "zlalsa.f"
	nrf = ic + 1;

/*        Since B and BX are complex, the following call to DGEMM is */
/*        performed in two steps (real and imaginary parts). */

/*        CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, */
/*    $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX ) */

#line 559 "zlalsa.f"
	j = nlp1 * *nrhs << 1;
#line 560 "zlalsa.f"
	i__2 = *nrhs;
#line 560 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 561 "zlalsa.f"
	    i__3 = nlf + nlp1 - 1;
#line 561 "zlalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 562 "zlalsa.f"
		++j;
#line 563 "zlalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 563 "zlalsa.f"
		rwork[j] = b[i__4].r;
#line 564 "zlalsa.f"
/* L200: */
#line 564 "zlalsa.f"
	    }
#line 565 "zlalsa.f"
/* L210: */
#line 565 "zlalsa.f"
	}
#line 566 "zlalsa.f"
	dgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b9, &vt[nlf + vt_dim1], ldu, &
		rwork[(nlp1 * *nrhs << 1) + 1], &nlp1, &c_b10, &rwork[1], &
		nlp1, (ftnlen)1, (ftnlen)1);
#line 569 "zlalsa.f"
	j = nlp1 * *nrhs << 1;
#line 570 "zlalsa.f"
	i__2 = *nrhs;
#line 570 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 571 "zlalsa.f"
	    i__3 = nlf + nlp1 - 1;
#line 571 "zlalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 572 "zlalsa.f"
		++j;
#line 573 "zlalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 574 "zlalsa.f"
/* L220: */
#line 574 "zlalsa.f"
	    }
#line 575 "zlalsa.f"
/* L230: */
#line 575 "zlalsa.f"
	}
#line 576 "zlalsa.f"
	dgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b9, &vt[nlf + vt_dim1], ldu, &
		rwork[(nlp1 * *nrhs << 1) + 1], &nlp1, &c_b10, &rwork[nlp1 * *
		nrhs + 1], &nlp1, (ftnlen)1, (ftnlen)1);
#line 579 "zlalsa.f"
	jreal = 0;
#line 580 "zlalsa.f"
	jimag = nlp1 * *nrhs;
#line 581 "zlalsa.f"
	i__2 = *nrhs;
#line 581 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 582 "zlalsa.f"
	    i__3 = nlf + nlp1 - 1;
#line 582 "zlalsa.f"
	    for (jrow = nlf; jrow <= i__3; ++jrow) {
#line 583 "zlalsa.f"
		++jreal;
#line 584 "zlalsa.f"
		++jimag;
#line 585 "zlalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 585 "zlalsa.f"
		i__5 = jreal;
#line 585 "zlalsa.f"
		i__6 = jimag;
#line 585 "zlalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 585 "zlalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 587 "zlalsa.f"
/* L240: */
#line 587 "zlalsa.f"
	    }
#line 588 "zlalsa.f"
/* L250: */
#line 588 "zlalsa.f"
	}

/*        Since B and BX are complex, the following call to DGEMM is */
/*        performed in two steps (real and imaginary parts). */

/*        CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, */
/*    $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX ) */

#line 596 "zlalsa.f"
	j = nrp1 * *nrhs << 1;
#line 597 "zlalsa.f"
	i__2 = *nrhs;
#line 597 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 598 "zlalsa.f"
	    i__3 = nrf + nrp1 - 1;
#line 598 "zlalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 599 "zlalsa.f"
		++j;
#line 600 "zlalsa.f"
		i__4 = jrow + jcol * b_dim1;
#line 600 "zlalsa.f"
		rwork[j] = b[i__4].r;
#line 601 "zlalsa.f"
/* L260: */
#line 601 "zlalsa.f"
	    }
#line 602 "zlalsa.f"
/* L270: */
#line 602 "zlalsa.f"
	}
#line 603 "zlalsa.f"
	dgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b9, &vt[nrf + vt_dim1], ldu, &
		rwork[(nrp1 * *nrhs << 1) + 1], &nrp1, &c_b10, &rwork[1], &
		nrp1, (ftnlen)1, (ftnlen)1);
#line 606 "zlalsa.f"
	j = nrp1 * *nrhs << 1;
#line 607 "zlalsa.f"
	i__2 = *nrhs;
#line 607 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 608 "zlalsa.f"
	    i__3 = nrf + nrp1 - 1;
#line 608 "zlalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 609 "zlalsa.f"
		++j;
#line 610 "zlalsa.f"
		rwork[j] = d_imag(&b[jrow + jcol * b_dim1]);
#line 611 "zlalsa.f"
/* L280: */
#line 611 "zlalsa.f"
	    }
#line 612 "zlalsa.f"
/* L290: */
#line 612 "zlalsa.f"
	}
#line 613 "zlalsa.f"
	dgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b9, &vt[nrf + vt_dim1], ldu, &
		rwork[(nrp1 * *nrhs << 1) + 1], &nrp1, &c_b10, &rwork[nrp1 * *
		nrhs + 1], &nrp1, (ftnlen)1, (ftnlen)1);
#line 616 "zlalsa.f"
	jreal = 0;
#line 617 "zlalsa.f"
	jimag = nrp1 * *nrhs;
#line 618 "zlalsa.f"
	i__2 = *nrhs;
#line 618 "zlalsa.f"
	for (jcol = 1; jcol <= i__2; ++jcol) {
#line 619 "zlalsa.f"
	    i__3 = nrf + nrp1 - 1;
#line 619 "zlalsa.f"
	    for (jrow = nrf; jrow <= i__3; ++jrow) {
#line 620 "zlalsa.f"
		++jreal;
#line 621 "zlalsa.f"
		++jimag;
#line 622 "zlalsa.f"
		i__4 = jrow + jcol * bx_dim1;
#line 622 "zlalsa.f"
		i__5 = jreal;
#line 622 "zlalsa.f"
		i__6 = jimag;
#line 622 "zlalsa.f"
		z__1.r = rwork[i__5], z__1.i = rwork[i__6];
#line 622 "zlalsa.f"
		bx[i__4].r = z__1.r, bx[i__4].i = z__1.i;
#line 624 "zlalsa.f"
/* L300: */
#line 624 "zlalsa.f"
	    }
#line 625 "zlalsa.f"
/* L310: */
#line 625 "zlalsa.f"
	}

#line 627 "zlalsa.f"
/* L320: */
#line 627 "zlalsa.f"
    }

#line 629 "zlalsa.f"
L330:

#line 631 "zlalsa.f"
    return 0;

/*     End of ZLALSA */

} /* zlalsa_ */

