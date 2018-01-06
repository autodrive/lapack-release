#line 1 "slasda.f"
/* slasda.f -- translated by f2c (version 20100827).
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

#line 1 "slasda.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b11 = 0.;
static doublereal c_b12 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLASDA computes the singular value decomposition (SVD) of a real upper bidiagonal matrix with d
iagonal d and off-diagonal e. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASDA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasda.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasda.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasda.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, */
/*                          DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, */
/*                          PERM, GIVNUM, C, S, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), */
/*      $                   K( * ), PERM( LDGCOL, * ) */
/*       REAL               C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ), */
/*      $                   E( * ), GIVNUM( LDU, * ), POLES( LDU, * ), */
/*      $                   S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), */
/*      $                   Z( LDU, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Using a divide and conquer approach, SLASDA computes the singular */
/* > value decomposition (SVD) of a real upper bidiagonal N-by-M matrix */
/* > B with diagonal D and offdiagonal E, where M = N + SQRE. The */
/* > algorithm computes the singular values in the SVD B = U * S * VT. */
/* > The orthogonal matrices U and VT are optionally computed in */
/* > compact form. */
/* > */
/* > A related subroutine, SLASD0, computes the singular values and */
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
/* >          D is REAL array, dimension ( N ) */
/* >         On entry D contains the main diagonal of the bidiagonal */
/* >         matrix. On exit D, if INFO = 0, contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension ( M-1 ) */
/* >         Contains the subdiagonal entries of the bidiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, */
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
/* >          VT is REAL array, */
/* >         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced */
/* >         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT**T contains the right */
/* >         singular vector matrices of all subproblems at the bottom */
/* >         level. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER array, dimension ( N ) */
/* >         if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0. */
/* >         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th */
/* >         secular equation on the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFL */
/* > \verbatim */
/* >          DIFL is REAL array, dimension ( LDU, NLVL ), */
/* >         where NLVL = floor(log_2 (N/SMLSIZ))). */
/* > \endverbatim */
/* > */
/* > \param[out] DIFR */
/* > \verbatim */
/* >          DIFR is REAL array, */
/* >                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and */
/* >                  dimension ( N ) if ICOMPQ = 0. */
/* >         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1) */
/* >         record distances between singular values on the I-th */
/* >         level and singular values on the (I -1)-th level, and */
/* >         DIFR(1:N, 2 * I ) contains the normalizing factors for */
/* >         the right singular vector matrix. See SLASD8 for details. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, */
/* >                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and */
/* >                  dimension ( N ) if ICOMPQ = 0. */
/* >         The first K elements of Z(1, I) contain the components of */
/* >         the deflation-adjusted updating row vector for subproblems */
/* >         on the I-th level. */
/* > \endverbatim */
/* > */
/* > \param[out] POLES */
/* > \verbatim */
/* >          POLES is REAL array, */
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
/* >          PERM is INTEGER array, dimension ( LDGCOL, NLVL ) */
/* >         if ICOMPQ = 1, and not referenced */
/* >         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records */
/* >         permutations done on the I-th level of the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is REAL array, */
/* >         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not */
/* >         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I, */
/* >         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S- */
/* >         values of Givens rotations performed on the I-th level on */
/* >         the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL array, */
/* >         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. */
/* >         If ICOMPQ = 1 and the I-th subproblem is not square, on exit, */
/* >         C( I ) contains the C-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL array, dimension ( N ) if */
/* >         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1 */
/* >         and the I-th subproblem is not square, on exit, S( I ) */
/* >         contains the S-value of a Givens rotation related to */
/* >         the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension */
/* >         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (7*N). */
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

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasda_(integer *icompq, integer *smlsiz, integer *n, 
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
    static integer inode, ndiml, ndimr, idxqi, itemp, sqrei;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slasd6_(integer *, integer *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    static integer nwork1, nwork2;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slasdq_(
	    char *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), slasdt_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), slaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static integer smlszp;


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 313 "slasda.f"
    /* Parameter adjustments */
#line 313 "slasda.f"
    --d__;
#line 313 "slasda.f"
    --e;
#line 313 "slasda.f"
    givnum_dim1 = *ldu;
#line 313 "slasda.f"
    givnum_offset = 1 + givnum_dim1;
#line 313 "slasda.f"
    givnum -= givnum_offset;
#line 313 "slasda.f"
    poles_dim1 = *ldu;
#line 313 "slasda.f"
    poles_offset = 1 + poles_dim1;
#line 313 "slasda.f"
    poles -= poles_offset;
#line 313 "slasda.f"
    z_dim1 = *ldu;
#line 313 "slasda.f"
    z_offset = 1 + z_dim1;
#line 313 "slasda.f"
    z__ -= z_offset;
#line 313 "slasda.f"
    difr_dim1 = *ldu;
#line 313 "slasda.f"
    difr_offset = 1 + difr_dim1;
#line 313 "slasda.f"
    difr -= difr_offset;
#line 313 "slasda.f"
    difl_dim1 = *ldu;
#line 313 "slasda.f"
    difl_offset = 1 + difl_dim1;
#line 313 "slasda.f"
    difl -= difl_offset;
#line 313 "slasda.f"
    vt_dim1 = *ldu;
#line 313 "slasda.f"
    vt_offset = 1 + vt_dim1;
#line 313 "slasda.f"
    vt -= vt_offset;
#line 313 "slasda.f"
    u_dim1 = *ldu;
#line 313 "slasda.f"
    u_offset = 1 + u_dim1;
#line 313 "slasda.f"
    u -= u_offset;
#line 313 "slasda.f"
    --k;
#line 313 "slasda.f"
    --givptr;
#line 313 "slasda.f"
    perm_dim1 = *ldgcol;
#line 313 "slasda.f"
    perm_offset = 1 + perm_dim1;
#line 313 "slasda.f"
    perm -= perm_offset;
#line 313 "slasda.f"
    givcol_dim1 = *ldgcol;
#line 313 "slasda.f"
    givcol_offset = 1 + givcol_dim1;
#line 313 "slasda.f"
    givcol -= givcol_offset;
#line 313 "slasda.f"
    --c__;
#line 313 "slasda.f"
    --s;
#line 313 "slasda.f"
    --work;
#line 313 "slasda.f"
    --iwork;
#line 313 "slasda.f"

#line 313 "slasda.f"
    /* Function Body */
#line 313 "slasda.f"
    *info = 0;

#line 315 "slasda.f"
    if (*icompq < 0 || *icompq > 1) {
#line 316 "slasda.f"
	*info = -1;
#line 317 "slasda.f"
    } else if (*smlsiz < 3) {
#line 318 "slasda.f"
	*info = -2;
#line 319 "slasda.f"
    } else if (*n < 0) {
#line 320 "slasda.f"
	*info = -3;
#line 321 "slasda.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 322 "slasda.f"
	*info = -4;
#line 323 "slasda.f"
    } else if (*ldu < *n + *sqre) {
#line 324 "slasda.f"
	*info = -8;
#line 325 "slasda.f"
    } else if (*ldgcol < *n) {
#line 326 "slasda.f"
	*info = -17;
#line 327 "slasda.f"
    }
#line 328 "slasda.f"
    if (*info != 0) {
#line 329 "slasda.f"
	i__1 = -(*info);
#line 329 "slasda.f"
	xerbla_("SLASDA", &i__1, (ftnlen)6);
#line 330 "slasda.f"
	return 0;
#line 331 "slasda.f"
    }

#line 333 "slasda.f"
    m = *n + *sqre;

/*     If the input matrix is too small, call SLASDQ to find the SVD. */

#line 337 "slasda.f"
    if (*n <= *smlsiz) {
#line 338 "slasda.f"
	if (*icompq == 0) {
#line 339 "slasda.f"
	    slasdq_("U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[
		    vt_offset], ldu, &u[u_offset], ldu, &u[u_offset], ldu, &
		    work[1], info, (ftnlen)1);
#line 341 "slasda.f"
	} else {
#line 342 "slasda.f"
	    slasdq_("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset]
		    , ldu, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], 
		    info, (ftnlen)1);
#line 344 "slasda.f"
	}
#line 345 "slasda.f"
	return 0;
#line 346 "slasda.f"
    }

/*     Book-keeping and  set up the computation tree. */

#line 350 "slasda.f"
    inode = 1;
#line 351 "slasda.f"
    ndiml = inode + *n;
#line 352 "slasda.f"
    ndimr = ndiml + *n;
#line 353 "slasda.f"
    idxq = ndimr + *n;
#line 354 "slasda.f"
    iwk = idxq + *n;

#line 356 "slasda.f"
    ncc = 0;
#line 357 "slasda.f"
    nru = 0;

#line 359 "slasda.f"
    smlszp = *smlsiz + 1;
#line 360 "slasda.f"
    vf = 1;
#line 361 "slasda.f"
    vl = vf + m;
#line 362 "slasda.f"
    nwork1 = vl + m;
#line 363 "slasda.f"
    nwork2 = nwork1 + smlszp * smlszp;

#line 365 "slasda.f"
    slasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     for the nodes on bottom level of the tree, solve */
/*     their subproblems by SLASDQ. */

#line 371 "slasda.f"
    ndb1 = (nd + 1) / 2;
#line 372 "slasda.f"
    i__1 = nd;
#line 372 "slasda.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*        IC : center row of each node */
/*        NL : number of rows of left  subproblem */
/*        NR : number of rows of right subproblem */
/*        NLF: starting row of the left   subproblem */
/*        NRF: starting row of the right  subproblem */

#line 380 "slasda.f"
	i1 = i__ - 1;
#line 381 "slasda.f"
	ic = iwork[inode + i1];
#line 382 "slasda.f"
	nl = iwork[ndiml + i1];
#line 383 "slasda.f"
	nlp1 = nl + 1;
#line 384 "slasda.f"
	nr = iwork[ndimr + i1];
#line 385 "slasda.f"
	nlf = ic - nl;
#line 386 "slasda.f"
	nrf = ic + 1;
#line 387 "slasda.f"
	idxqi = idxq + nlf - 2;
#line 388 "slasda.f"
	vfi = vf + nlf - 1;
#line 389 "slasda.f"
	vli = vl + nlf - 1;
#line 390 "slasda.f"
	sqrei = 1;
#line 391 "slasda.f"
	if (*icompq == 0) {
#line 392 "slasda.f"
	    slaset_("A", &nlp1, &nlp1, &c_b11, &c_b12, &work[nwork1], &smlszp,
		     (ftnlen)1);
#line 394 "slasda.f"
	    slasdq_("U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &
		    work[nwork1], &smlszp, &work[nwork2], &nl, &work[nwork2], 
		    &nl, &work[nwork2], info, (ftnlen)1);
#line 398 "slasda.f"
	    itemp = nwork1 + nl * smlszp;
#line 399 "slasda.f"
	    scopy_(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
#line 400 "slasda.f"
	    scopy_(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
#line 401 "slasda.f"
	} else {
#line 402 "slasda.f"
	    slaset_("A", &nl, &nl, &c_b11, &c_b12, &u[nlf + u_dim1], ldu, (
		    ftnlen)1);
#line 403 "slasda.f"
	    slaset_("A", &nlp1, &nlp1, &c_b11, &c_b12, &vt[nlf + vt_dim1], 
		    ldu, (ftnlen)1);
#line 404 "slasda.f"
	    slasdq_("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &
		    vt[nlf + vt_dim1], ldu, &u[nlf + u_dim1], ldu, &u[nlf + 
		    u_dim1], ldu, &work[nwork1], info, (ftnlen)1);
#line 407 "slasda.f"
	    scopy_(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
#line 408 "slasda.f"
	    scopy_(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
#line 409 "slasda.f"
	}
#line 410 "slasda.f"
	if (*info != 0) {
#line 411 "slasda.f"
	    return 0;
#line 412 "slasda.f"
	}
#line 413 "slasda.f"
	i__2 = nl;
#line 413 "slasda.f"
	for (j = 1; j <= i__2; ++j) {
#line 414 "slasda.f"
	    iwork[idxqi + j] = j;
#line 415 "slasda.f"
/* L10: */
#line 415 "slasda.f"
	}
#line 416 "slasda.f"
	if (i__ == nd && *sqre == 0) {
#line 417 "slasda.f"
	    sqrei = 0;
#line 418 "slasda.f"
	} else {
#line 419 "slasda.f"
	    sqrei = 1;
#line 420 "slasda.f"
	}
#line 421 "slasda.f"
	idxqi += nlp1;
#line 422 "slasda.f"
	vfi += nlp1;
#line 423 "slasda.f"
	vli += nlp1;
#line 424 "slasda.f"
	nrp1 = nr + sqrei;
#line 425 "slasda.f"
	if (*icompq == 0) {
#line 426 "slasda.f"
	    slaset_("A", &nrp1, &nrp1, &c_b11, &c_b12, &work[nwork1], &smlszp,
		     (ftnlen)1);
#line 428 "slasda.f"
	    slasdq_("U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &
		    work[nwork1], &smlszp, &work[nwork2], &nr, &work[nwork2], 
		    &nr, &work[nwork2], info, (ftnlen)1);
#line 432 "slasda.f"
	    itemp = nwork1 + (nrp1 - 1) * smlszp;
#line 433 "slasda.f"
	    scopy_(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
#line 434 "slasda.f"
	    scopy_(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
#line 435 "slasda.f"
	} else {
#line 436 "slasda.f"
	    slaset_("A", &nr, &nr, &c_b11, &c_b12, &u[nrf + u_dim1], ldu, (
		    ftnlen)1);
#line 437 "slasda.f"
	    slaset_("A", &nrp1, &nrp1, &c_b11, &c_b12, &vt[nrf + vt_dim1], 
		    ldu, (ftnlen)1);
#line 438 "slasda.f"
	    slasdq_("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &
		    vt[nrf + vt_dim1], ldu, &u[nrf + u_dim1], ldu, &u[nrf + 
		    u_dim1], ldu, &work[nwork1], info, (ftnlen)1);
#line 441 "slasda.f"
	    scopy_(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
#line 442 "slasda.f"
	    scopy_(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
#line 443 "slasda.f"
	}
#line 444 "slasda.f"
	if (*info != 0) {
#line 445 "slasda.f"
	    return 0;
#line 446 "slasda.f"
	}
#line 447 "slasda.f"
	i__2 = nr;
#line 447 "slasda.f"
	for (j = 1; j <= i__2; ++j) {
#line 448 "slasda.f"
	    iwork[idxqi + j] = j;
#line 449 "slasda.f"
/* L20: */
#line 449 "slasda.f"
	}
#line 450 "slasda.f"
/* L30: */
#line 450 "slasda.f"
    }

/*     Now conquer each subproblem bottom-up. */

#line 454 "slasda.f"
    j = pow_ii(&c__2, &nlvl);
#line 455 "slasda.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {
#line 456 "slasda.f"
	lvl2 = (lvl << 1) - 1;

/*        Find the first node LF and last node LL on */
/*        the current level LVL. */

#line 461 "slasda.f"
	if (lvl == 1) {
#line 462 "slasda.f"
	    lf = 1;
#line 463 "slasda.f"
	    ll = 1;
#line 464 "slasda.f"
	} else {
#line 465 "slasda.f"
	    i__1 = lvl - 1;
#line 465 "slasda.f"
	    lf = pow_ii(&c__2, &i__1);
#line 466 "slasda.f"
	    ll = (lf << 1) - 1;
#line 467 "slasda.f"
	}
#line 468 "slasda.f"
	i__1 = ll;
#line 468 "slasda.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 469 "slasda.f"
	    im1 = i__ - 1;
#line 470 "slasda.f"
	    ic = iwork[inode + im1];
#line 471 "slasda.f"
	    nl = iwork[ndiml + im1];
#line 472 "slasda.f"
	    nr = iwork[ndimr + im1];
#line 473 "slasda.f"
	    nlf = ic - nl;
#line 474 "slasda.f"
	    nrf = ic + 1;
#line 475 "slasda.f"
	    if (i__ == ll) {
#line 476 "slasda.f"
		sqrei = *sqre;
#line 477 "slasda.f"
	    } else {
#line 478 "slasda.f"
		sqrei = 1;
#line 479 "slasda.f"
	    }
#line 480 "slasda.f"
	    vfi = vf + nlf - 1;
#line 481 "slasda.f"
	    vli = vl + nlf - 1;
#line 482 "slasda.f"
	    idxqi = idxq + nlf - 1;
#line 483 "slasda.f"
	    alpha = d__[ic];
#line 484 "slasda.f"
	    beta = e[ic];
#line 485 "slasda.f"
	    if (*icompq == 0) {
#line 486 "slasda.f"
		slasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[
			perm_offset], &givptr[1], &givcol[givcol_offset], 
			ldgcol, &givnum[givnum_offset], ldu, &poles[
			poles_offset], &difl[difl_offset], &difr[difr_offset],
			 &z__[z_offset], &k[1], &c__[1], &s[1], &work[nwork1],
			 &iwork[iwk], info);
#line 492 "slasda.f"
	    } else {
#line 493 "slasda.f"
		--j;
#line 494 "slasda.f"
		slasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[nlf + 
			lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * 
			givcol_dim1], ldgcol, &givnum[nlf + lvl2 * 
			givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &
			difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * 
			difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[j], 
			&s[j], &work[nwork1], &iwork[iwk], info);
#line 503 "slasda.f"
	    }
#line 504 "slasda.f"
	    if (*info != 0) {
#line 505 "slasda.f"
		return 0;
#line 506 "slasda.f"
	    }
#line 507 "slasda.f"
/* L40: */
#line 507 "slasda.f"
	}
#line 508 "slasda.f"
/* L50: */
#line 508 "slasda.f"
    }

#line 510 "slasda.f"
    return 0;

/*     End of SLASDA */

} /* slasda_ */

