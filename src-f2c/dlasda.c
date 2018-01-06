#line 1 "dlasda.f"
/* dlasda.f -- translated by f2c (version 20100827).
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

#line 1 "dlasda.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b11 = 0.;
static doublereal c_b12 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DLASDA computes the singular value decomposition (SVD) of a real upper bidiagonal matrix with d
iagonal d and off-diagonal e. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASDA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasda.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasda.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasda.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, */
/*                          DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, */
/*                          PERM, GIVNUM, C, S, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), */
/*      $                   K( * ), PERM( LDGCOL, * ) */
/*       DOUBLE PRECISION   C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ), */
/*      $                   E( * ), GIVNUM( LDU, * ), POLES( LDU, * ), */
/*      $                   S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), */
/*      $                   Z( LDU, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Using a divide and conquer approach, DLASDA computes the singular */
/* > value decomposition (SVD) of a real upper bidiagonal N-by-M matrix */
/* > B with diagonal D and offdiagonal E, where M = N + SQRE. The */
/* > algorithm computes the singular values in the SVD B = U * S * VT. */
/* > The orthogonal matrices U and VT are optionally computed in */
/* > compact form. */
/* > */
/* > A related subroutine, DLASD0, computes the singular values and */
/* > the singular vectors in explicit form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >         Specifies whether singular vectors are to be computed */
/* >         in compact form, as follows */
/* >         = 0: Compute singular values only. */
/* >         = 1: Compute singular vectors of upper bidiagonal */
/* >              matrix in compact form. */
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
/* >         The row dimension of the upper bidiagonal matrix. This is */
/* >         also the dimension of the main diagonal array D. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* >          SQRE is INTEGER */
/* >         Specifies the column dimension of the bidiagonal matrix. */
/* >         = 0: The bidiagonal matrix has column dimension M = N; */
/* >         = 1: The bidiagonal matrix has column dimension M = N + 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension ( N ) */
/* >         On entry D contains the main diagonal of the bidiagonal */
/* >         matrix. On exit D, if INFO = 0, contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension ( M-1 ) */
/* >         Contains the subdiagonal entries of the bidiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, */
/* >         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced */
/* >         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left */
/* >         singular vector matrices of all subproblems at the bottom */
/* >         level. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER, LDU = > N. */
/* >         The leading dimension of arrays U, VT, DIFL, DIFR, POLES, */
/* >         GIVNUM, and Z. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, */
/* >         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced */
/* >         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT**T contains the right */
/* >         singular vector matrices of all subproblems at the bottom */
/* >         level. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER array, */
/* >         dimension ( N ) if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0. */
/* >         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th */
/* >         secular equation on the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFL */
/* > \verbatim */
/* >          DIFL is DOUBLE PRECISION array, dimension ( LDU, NLVL ), */
/* >         where NLVL = floor(log_2 (N/SMLSIZ))). */
/* > \endverbatim */
/* > */
/* > \param[out] DIFR */
/* > \verbatim */
/* >          DIFR is DOUBLE PRECISION array, */
/* >                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and */
/* >                  dimension ( N ) if ICOMPQ = 0. */
/* >         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1) */
/* >         record distances between singular values on the I-th */
/* >         level and singular values on the (I -1)-th level, and */
/* >         DIFR(1:N, 2 * I ) contains the normalizing factors for */
/* >         the right singular vector matrix. See DLASD8 for details. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, */
/* >                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and */
/* >                  dimension ( N ) if ICOMPQ = 0. */
/* >         The first K elements of Z(1, I) contain the components of */
/* >         the deflation-adjusted updating row vector for subproblems */
/* >         on the I-th level. */
/* > \endverbatim */
/* > */
/* > \param[out] POLES */
/* > \verbatim */
/* >          POLES is DOUBLE PRECISION array, */
/* >         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced */
/* >         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and */
/* >         POLES(1, 2*I) contain  the new and old singular values */
/* >         involved in the secular equations on the I-th level. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER array, */
/* >         dimension ( N ) if ICOMPQ = 1, and not referenced if */
/* >         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records */
/* >         the number of Givens rotations performed on the I-th */
/* >         problem on the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, */
/* >         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not */
/* >         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I, */
/* >         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations */
/* >         of Givens rotations performed on the I-th level on the */
/* >         computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* >          LDGCOL is INTEGER, LDGCOL = > N. */
/* >         The leading dimension of arrays GIVCOL and PERM. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, */
/* >         dimension ( LDGCOL, NLVL ) if ICOMPQ = 1, and not referenced */
/* >         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records */
/* >         permutations done on the I-th level of the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is DOUBLE PRECISION array, */
/* >         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not */
/* >         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I, */
/* >         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S- */
/* >         values of Givens rotations performed on the I-th level on */
/* >         the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, */
/* >         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. */
/* >         If ICOMPQ = 1 and the I-th subproblem is not square, on exit, */
/* >         C( I ) contains the C-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension ( N ) if */
/* >         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1 */
/* >         and the I-th subproblem is not square, on exit, S( I ) */
/* >         contains the S-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension */
/* >         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (7*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = 1, a singular value did not converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasda_(integer *icompq, integer *smlsiz, integer *n, 
	integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer 
	*ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, 
	doublereal *z__, doublereal *poles, integer *givptr, integer *givcol, 
	integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__, 
	doublereal *s, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, 
	    difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, 
	    z_dim1, z_offset, i__1, i__2;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, m, i1, ic, lf, nd, ll, nl, vf, nr, vl, im1, ncc, 
	    nlf, nrf, vfi, iwk, vli, lvl, nru, ndb1, nlp1, lvl2, nrp1;
    static doublereal beta;
    static integer idxq, nlvl;
    static doublereal alpha;
    static integer inode, ndiml, ndimr, idxqi, itemp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer sqrei;
    extern /* Subroutine */ int dlasd6_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    static integer nwork1, nwork2;
    extern /* Subroutine */ int dlasdq_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlasdt_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), dlaset_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer smlszp;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
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

#line 313 "dlasda.f"
    /* Parameter adjustments */
#line 313 "dlasda.f"
    --d__;
#line 313 "dlasda.f"
    --e;
#line 313 "dlasda.f"
    givnum_dim1 = *ldu;
#line 313 "dlasda.f"
    givnum_offset = 1 + givnum_dim1;
#line 313 "dlasda.f"
    givnum -= givnum_offset;
#line 313 "dlasda.f"
    poles_dim1 = *ldu;
#line 313 "dlasda.f"
    poles_offset = 1 + poles_dim1;
#line 313 "dlasda.f"
    poles -= poles_offset;
#line 313 "dlasda.f"
    z_dim1 = *ldu;
#line 313 "dlasda.f"
    z_offset = 1 + z_dim1;
#line 313 "dlasda.f"
    z__ -= z_offset;
#line 313 "dlasda.f"
    difr_dim1 = *ldu;
#line 313 "dlasda.f"
    difr_offset = 1 + difr_dim1;
#line 313 "dlasda.f"
    difr -= difr_offset;
#line 313 "dlasda.f"
    difl_dim1 = *ldu;
#line 313 "dlasda.f"
    difl_offset = 1 + difl_dim1;
#line 313 "dlasda.f"
    difl -= difl_offset;
#line 313 "dlasda.f"
    vt_dim1 = *ldu;
#line 313 "dlasda.f"
    vt_offset = 1 + vt_dim1;
#line 313 "dlasda.f"
    vt -= vt_offset;
#line 313 "dlasda.f"
    u_dim1 = *ldu;
#line 313 "dlasda.f"
    u_offset = 1 + u_dim1;
#line 313 "dlasda.f"
    u -= u_offset;
#line 313 "dlasda.f"
    --k;
#line 313 "dlasda.f"
    --givptr;
#line 313 "dlasda.f"
    perm_dim1 = *ldgcol;
#line 313 "dlasda.f"
    perm_offset = 1 + perm_dim1;
#line 313 "dlasda.f"
    perm -= perm_offset;
#line 313 "dlasda.f"
    givcol_dim1 = *ldgcol;
#line 313 "dlasda.f"
    givcol_offset = 1 + givcol_dim1;
#line 313 "dlasda.f"
    givcol -= givcol_offset;
#line 313 "dlasda.f"
    --c__;
#line 313 "dlasda.f"
    --s;
#line 313 "dlasda.f"
    --work;
#line 313 "dlasda.f"
    --iwork;
#line 313 "dlasda.f"

#line 313 "dlasda.f"
    /* Function Body */
#line 313 "dlasda.f"
    *info = 0;

#line 315 "dlasda.f"
    if (*icompq < 0 || *icompq > 1) {
#line 316 "dlasda.f"
	*info = -1;
#line 317 "dlasda.f"
    } else if (*smlsiz < 3) {
#line 318 "dlasda.f"
	*info = -2;
#line 319 "dlasda.f"
    } else if (*n < 0) {
#line 320 "dlasda.f"
	*info = -3;
#line 321 "dlasda.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 322 "dlasda.f"
	*info = -4;
#line 323 "dlasda.f"
    } else if (*ldu < *n + *sqre) {
#line 324 "dlasda.f"
	*info = -8;
#line 325 "dlasda.f"
    } else if (*ldgcol < *n) {
#line 326 "dlasda.f"
	*info = -17;
#line 327 "dlasda.f"
    }
#line 328 "dlasda.f"
    if (*info != 0) {
#line 329 "dlasda.f"
	i__1 = -(*info);
#line 329 "dlasda.f"
	xerbla_("DLASDA", &i__1, (ftnlen)6);
#line 330 "dlasda.f"
	return 0;
#line 331 "dlasda.f"
    }

#line 333 "dlasda.f"
    m = *n + *sqre;

/*     If the input matrix is too small, call DLASDQ to find the SVD. */

#line 337 "dlasda.f"
    if (*n <= *smlsiz) {
#line 338 "dlasda.f"
	if (*icompq == 0) {
#line 339 "dlasda.f"
	    dlasdq_("U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[
		    vt_offset], ldu, &u[u_offset], ldu, &u[u_offset], ldu, &
		    work[1], info, (ftnlen)1);
#line 341 "dlasda.f"
	} else {
#line 342 "dlasda.f"
	    dlasdq_("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset]
		    , ldu, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], 
		    info, (ftnlen)1);
#line 344 "dlasda.f"
	}
#line 345 "dlasda.f"
	return 0;
#line 346 "dlasda.f"
    }

/*     Book-keeping and  set up the computation tree. */

#line 350 "dlasda.f"
    inode = 1;
#line 351 "dlasda.f"
    ndiml = inode + *n;
#line 352 "dlasda.f"
    ndimr = ndiml + *n;
#line 353 "dlasda.f"
    idxq = ndimr + *n;
#line 354 "dlasda.f"
    iwk = idxq + *n;

#line 356 "dlasda.f"
    ncc = 0;
#line 357 "dlasda.f"
    nru = 0;

#line 359 "dlasda.f"
    smlszp = *smlsiz + 1;
#line 360 "dlasda.f"
    vf = 1;
#line 361 "dlasda.f"
    vl = vf + m;
#line 362 "dlasda.f"
    nwork1 = vl + m;
#line 363 "dlasda.f"
    nwork2 = nwork1 + smlszp * smlszp;

#line 365 "dlasda.f"
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     for the nodes on bottom level of the tree, solve */
/*     their subproblems by DLASDQ. */

#line 371 "dlasda.f"
    ndb1 = (nd + 1) / 2;
#line 372 "dlasda.f"
    i__1 = nd;
#line 372 "dlasda.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*        IC : center row of each node */
/*        NL : number of rows of left  subproblem */
/*        NR : number of rows of right subproblem */
/*        NLF: starting row of the left   subproblem */
/*        NRF: starting row of the right  subproblem */

#line 380 "dlasda.f"
	i1 = i__ - 1;
#line 381 "dlasda.f"
	ic = iwork[inode + i1];
#line 382 "dlasda.f"
	nl = iwork[ndiml + i1];
#line 383 "dlasda.f"
	nlp1 = nl + 1;
#line 384 "dlasda.f"
	nr = iwork[ndimr + i1];
#line 385 "dlasda.f"
	nlf = ic - nl;
#line 386 "dlasda.f"
	nrf = ic + 1;
#line 387 "dlasda.f"
	idxqi = idxq + nlf - 2;
#line 388 "dlasda.f"
	vfi = vf + nlf - 1;
#line 389 "dlasda.f"
	vli = vl + nlf - 1;
#line 390 "dlasda.f"
	sqrei = 1;
#line 391 "dlasda.f"
	if (*icompq == 0) {
#line 392 "dlasda.f"
	    dlaset_("A", &nlp1, &nlp1, &c_b11, &c_b12, &work[nwork1], &smlszp,
		     (ftnlen)1);
#line 394 "dlasda.f"
	    dlasdq_("U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &
		    work[nwork1], &smlszp, &work[nwork2], &nl, &work[nwork2], 
		    &nl, &work[nwork2], info, (ftnlen)1);
#line 398 "dlasda.f"
	    itemp = nwork1 + nl * smlszp;
#line 399 "dlasda.f"
	    dcopy_(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
#line 400 "dlasda.f"
	    dcopy_(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
#line 401 "dlasda.f"
	} else {
#line 402 "dlasda.f"
	    dlaset_("A", &nl, &nl, &c_b11, &c_b12, &u[nlf + u_dim1], ldu, (
		    ftnlen)1);
#line 403 "dlasda.f"
	    dlaset_("A", &nlp1, &nlp1, &c_b11, &c_b12, &vt[nlf + vt_dim1], 
		    ldu, (ftnlen)1);
#line 404 "dlasda.f"
	    dlasdq_("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &
		    vt[nlf + vt_dim1], ldu, &u[nlf + u_dim1], ldu, &u[nlf + 
		    u_dim1], ldu, &work[nwork1], info, (ftnlen)1);
#line 407 "dlasda.f"
	    dcopy_(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
#line 408 "dlasda.f"
	    dcopy_(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
#line 409 "dlasda.f"
	}
#line 410 "dlasda.f"
	if (*info != 0) {
#line 411 "dlasda.f"
	    return 0;
#line 412 "dlasda.f"
	}
#line 413 "dlasda.f"
	i__2 = nl;
#line 413 "dlasda.f"
	for (j = 1; j <= i__2; ++j) {
#line 414 "dlasda.f"
	    iwork[idxqi + j] = j;
#line 415 "dlasda.f"
/* L10: */
#line 415 "dlasda.f"
	}
#line 416 "dlasda.f"
	if (i__ == nd && *sqre == 0) {
#line 417 "dlasda.f"
	    sqrei = 0;
#line 418 "dlasda.f"
	} else {
#line 419 "dlasda.f"
	    sqrei = 1;
#line 420 "dlasda.f"
	}
#line 421 "dlasda.f"
	idxqi += nlp1;
#line 422 "dlasda.f"
	vfi += nlp1;
#line 423 "dlasda.f"
	vli += nlp1;
#line 424 "dlasda.f"
	nrp1 = nr + sqrei;
#line 425 "dlasda.f"
	if (*icompq == 0) {
#line 426 "dlasda.f"
	    dlaset_("A", &nrp1, &nrp1, &c_b11, &c_b12, &work[nwork1], &smlszp,
		     (ftnlen)1);
#line 428 "dlasda.f"
	    dlasdq_("U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &
		    work[nwork1], &smlszp, &work[nwork2], &nr, &work[nwork2], 
		    &nr, &work[nwork2], info, (ftnlen)1);
#line 432 "dlasda.f"
	    itemp = nwork1 + (nrp1 - 1) * smlszp;
#line 433 "dlasda.f"
	    dcopy_(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
#line 434 "dlasda.f"
	    dcopy_(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
#line 435 "dlasda.f"
	} else {
#line 436 "dlasda.f"
	    dlaset_("A", &nr, &nr, &c_b11, &c_b12, &u[nrf + u_dim1], ldu, (
		    ftnlen)1);
#line 437 "dlasda.f"
	    dlaset_("A", &nrp1, &nrp1, &c_b11, &c_b12, &vt[nrf + vt_dim1], 
		    ldu, (ftnlen)1);
#line 438 "dlasda.f"
	    dlasdq_("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &
		    vt[nrf + vt_dim1], ldu, &u[nrf + u_dim1], ldu, &u[nrf + 
		    u_dim1], ldu, &work[nwork1], info, (ftnlen)1);
#line 441 "dlasda.f"
	    dcopy_(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
#line 442 "dlasda.f"
	    dcopy_(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
#line 443 "dlasda.f"
	}
#line 444 "dlasda.f"
	if (*info != 0) {
#line 445 "dlasda.f"
	    return 0;
#line 446 "dlasda.f"
	}
#line 447 "dlasda.f"
	i__2 = nr;
#line 447 "dlasda.f"
	for (j = 1; j <= i__2; ++j) {
#line 448 "dlasda.f"
	    iwork[idxqi + j] = j;
#line 449 "dlasda.f"
/* L20: */
#line 449 "dlasda.f"
	}
#line 450 "dlasda.f"
/* L30: */
#line 450 "dlasda.f"
    }

/*     Now conquer each subproblem bottom-up. */

#line 454 "dlasda.f"
    j = pow_ii(&c__2, &nlvl);
#line 455 "dlasda.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {
#line 456 "dlasda.f"
	lvl2 = (lvl << 1) - 1;

/*        Find the first node LF and last node LL on */
/*        the current level LVL. */

#line 461 "dlasda.f"
	if (lvl == 1) {
#line 462 "dlasda.f"
	    lf = 1;
#line 463 "dlasda.f"
	    ll = 1;
#line 464 "dlasda.f"
	} else {
#line 465 "dlasda.f"
	    i__1 = lvl - 1;
#line 465 "dlasda.f"
	    lf = pow_ii(&c__2, &i__1);
#line 466 "dlasda.f"
	    ll = (lf << 1) - 1;
#line 467 "dlasda.f"
	}
#line 468 "dlasda.f"
	i__1 = ll;
#line 468 "dlasda.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 469 "dlasda.f"
	    im1 = i__ - 1;
#line 470 "dlasda.f"
	    ic = iwork[inode + im1];
#line 471 "dlasda.f"
	    nl = iwork[ndiml + im1];
#line 472 "dlasda.f"
	    nr = iwork[ndimr + im1];
#line 473 "dlasda.f"
	    nlf = ic - nl;
#line 474 "dlasda.f"
	    nrf = ic + 1;
#line 475 "dlasda.f"
	    if (i__ == ll) {
#line 476 "dlasda.f"
		sqrei = *sqre;
#line 477 "dlasda.f"
	    } else {
#line 478 "dlasda.f"
		sqrei = 1;
#line 479 "dlasda.f"
	    }
#line 480 "dlasda.f"
	    vfi = vf + nlf - 1;
#line 481 "dlasda.f"
	    vli = vl + nlf - 1;
#line 482 "dlasda.f"
	    idxqi = idxq + nlf - 1;
#line 483 "dlasda.f"
	    alpha = d__[ic];
#line 484 "dlasda.f"
	    beta = e[ic];
#line 485 "dlasda.f"
	    if (*icompq == 0) {
#line 486 "dlasda.f"
		dlasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[
			perm_offset], &givptr[1], &givcol[givcol_offset], 
			ldgcol, &givnum[givnum_offset], ldu, &poles[
			poles_offset], &difl[difl_offset], &difr[difr_offset],
			 &z__[z_offset], &k[1], &c__[1], &s[1], &work[nwork1],
			 &iwork[iwk], info);
#line 492 "dlasda.f"
	    } else {
#line 493 "dlasda.f"
		--j;
#line 494 "dlasda.f"
		dlasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[nlf + 
			lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * 
			givcol_dim1], ldgcol, &givnum[nlf + lvl2 * 
			givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &
			difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * 
			difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[j], 
			&s[j], &work[nwork1], &iwork[iwk], info);
#line 503 "dlasda.f"
	    }
#line 504 "dlasda.f"
	    if (*info != 0) {
#line 505 "dlasda.f"
		return 0;
#line 506 "dlasda.f"
	    }
#line 507 "dlasda.f"
/* L40: */
#line 507 "dlasda.f"
	}
#line 508 "dlasda.f"
/* L50: */
#line 508 "dlasda.f"
    }

#line 510 "dlasda.f"
    return 0;

/*     End of DLASDA */

} /* dlasda_ */

