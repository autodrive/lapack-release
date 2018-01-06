#line 1 "dlasd0.f"
/* dlasd0.f -- translated by f2c (version 20100827).
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

#line 1 "dlasd0.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__2 = 2;

/* > \brief \b DLASD0 computes the singular values of a real upper bidiagonal n-by-m matrix B with diagonal d 
and off-diagonal e. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASD0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDU, LDVT, N, SMLSIZ, SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Using a divide and conquer approach, DLASD0 computes the singular */
/* > value decomposition (SVD) of a real upper bidiagonal N-by-M */
/* > matrix B with diagonal D and offdiagonal E, where M = N + SQRE. */
/* > The algorithm computes orthogonal matrices U and VT such that */
/* > B = U * S * VT. The singular values S are overwritten on D. */
/* > */
/* > A related subroutine, DLASDA, computes only the singular values, */
/* > and optionally, the singular vectors in compact form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         On entry, the row dimension of the upper bidiagonal matrix. */
/* >         This is also the dimension of the main diagonal array D. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* >          SQRE is INTEGER */
/* >         Specifies the column dimension of the bidiagonal matrix. */
/* >         = 0: The bidiagonal matrix has column dimension M = N; */
/* >         = 1: The bidiagonal matrix has column dimension M = N+1; */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry D contains the main diagonal of the bidiagonal */
/* >         matrix. */
/* >         On exit D, if INFO = 0, contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (M-1) */
/* >         Contains the subdiagonal entries of the bidiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension at least (LDQ, N) */
/* >         On exit, U contains the left singular vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >         On entry, leading dimension of U. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension at least (LDVT, M) */
/* >         On exit, VT**T contains the right singular vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >         On entry, leading dimension of VT. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLSIZ */
/* > \verbatim */
/* >          SMLSIZ is INTEGER */
/* >         On entry, maximum size of the subproblems at the */
/* >         bottom of the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER work array. */
/* >         Dimension must be at least (8 * N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION work array. */
/* >         Dimension must be at least (3 * M**2 + 2 * M) */
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
/* Subroutine */ int dlasd0_(integer *n, integer *sqre, doublereal *d__, 
	doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *
	ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *
	info)
{
    /* System generated locals */
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, m, i1, ic, lf, nd, ll, nl, nr, im1, ncc, nlf, nrf, 
	    iwk, lvl, ndb1, nlp1, nrp1;
    static doublereal beta;
    static integer idxq, nlvl;
    static doublereal alpha;
    static integer inode, ndiml, idxqc, ndimr, itemp, sqrei;
    extern /* Subroutine */ int dlasd1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *), dlasdq_(char *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlasdt_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 184 "dlasd0.f"
    /* Parameter adjustments */
#line 184 "dlasd0.f"
    --d__;
#line 184 "dlasd0.f"
    --e;
#line 184 "dlasd0.f"
    u_dim1 = *ldu;
#line 184 "dlasd0.f"
    u_offset = 1 + u_dim1;
#line 184 "dlasd0.f"
    u -= u_offset;
#line 184 "dlasd0.f"
    vt_dim1 = *ldvt;
#line 184 "dlasd0.f"
    vt_offset = 1 + vt_dim1;
#line 184 "dlasd0.f"
    vt -= vt_offset;
#line 184 "dlasd0.f"
    --iwork;
#line 184 "dlasd0.f"
    --work;
#line 184 "dlasd0.f"

#line 184 "dlasd0.f"
    /* Function Body */
#line 184 "dlasd0.f"
    *info = 0;

#line 186 "dlasd0.f"
    if (*n < 0) {
#line 187 "dlasd0.f"
	*info = -1;
#line 188 "dlasd0.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 189 "dlasd0.f"
	*info = -2;
#line 190 "dlasd0.f"
    }

#line 192 "dlasd0.f"
    m = *n + *sqre;

#line 194 "dlasd0.f"
    if (*ldu < *n) {
#line 195 "dlasd0.f"
	*info = -6;
#line 196 "dlasd0.f"
    } else if (*ldvt < m) {
#line 197 "dlasd0.f"
	*info = -8;
#line 198 "dlasd0.f"
    } else if (*smlsiz < 3) {
#line 199 "dlasd0.f"
	*info = -9;
#line 200 "dlasd0.f"
    }
#line 201 "dlasd0.f"
    if (*info != 0) {
#line 202 "dlasd0.f"
	i__1 = -(*info);
#line 202 "dlasd0.f"
	xerbla_("DLASD0", &i__1, (ftnlen)6);
#line 203 "dlasd0.f"
	return 0;
#line 204 "dlasd0.f"
    }

/*     If the input matrix is too small, call DLASDQ to find the SVD. */

#line 208 "dlasd0.f"
    if (*n <= *smlsiz) {
#line 209 "dlasd0.f"
	dlasdq_("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], 
		ldvt, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info, (
		ftnlen)1);
#line 211 "dlasd0.f"
	return 0;
#line 212 "dlasd0.f"
    }

/*     Set up the computation tree. */

#line 216 "dlasd0.f"
    inode = 1;
#line 217 "dlasd0.f"
    ndiml = inode + *n;
#line 218 "dlasd0.f"
    ndimr = ndiml + *n;
#line 219 "dlasd0.f"
    idxq = ndimr + *n;
#line 220 "dlasd0.f"
    iwk = idxq + *n;
#line 221 "dlasd0.f"
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     For the nodes on bottom level of the tree, solve */
/*     their subproblems by DLASDQ. */

#line 227 "dlasd0.f"
    ndb1 = (nd + 1) / 2;
#line 228 "dlasd0.f"
    ncc = 0;
#line 229 "dlasd0.f"
    i__1 = nd;
#line 229 "dlasd0.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*     IC : center row of each node */
/*     NL : number of rows of left  subproblem */
/*     NR : number of rows of right subproblem */
/*     NLF: starting row of the left   subproblem */
/*     NRF: starting row of the right  subproblem */

#line 237 "dlasd0.f"
	i1 = i__ - 1;
#line 238 "dlasd0.f"
	ic = iwork[inode + i1];
#line 239 "dlasd0.f"
	nl = iwork[ndiml + i1];
#line 240 "dlasd0.f"
	nlp1 = nl + 1;
#line 241 "dlasd0.f"
	nr = iwork[ndimr + i1];
#line 242 "dlasd0.f"
	nrp1 = nr + 1;
#line 243 "dlasd0.f"
	nlf = ic - nl;
#line 244 "dlasd0.f"
	nrf = ic + 1;
#line 245 "dlasd0.f"
	sqrei = 1;
#line 246 "dlasd0.f"
	dlasdq_("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[
		nlf + nlf * vt_dim1], ldvt, &u[nlf + nlf * u_dim1], ldu, &u[
		nlf + nlf * u_dim1], ldu, &work[1], info, (ftnlen)1);
#line 249 "dlasd0.f"
	if (*info != 0) {
#line 250 "dlasd0.f"
	    return 0;
#line 251 "dlasd0.f"
	}
#line 252 "dlasd0.f"
	itemp = idxq + nlf - 2;
#line 253 "dlasd0.f"
	i__2 = nl;
#line 253 "dlasd0.f"
	for (j = 1; j <= i__2; ++j) {
#line 254 "dlasd0.f"
	    iwork[itemp + j] = j;
#line 255 "dlasd0.f"
/* L10: */
#line 255 "dlasd0.f"
	}
#line 256 "dlasd0.f"
	if (i__ == nd) {
#line 257 "dlasd0.f"
	    sqrei = *sqre;
#line 258 "dlasd0.f"
	} else {
#line 259 "dlasd0.f"
	    sqrei = 1;
#line 260 "dlasd0.f"
	}
#line 261 "dlasd0.f"
	nrp1 = nr + sqrei;
#line 262 "dlasd0.f"
	dlasdq_("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[
		nrf + nrf * vt_dim1], ldvt, &u[nrf + nrf * u_dim1], ldu, &u[
		nrf + nrf * u_dim1], ldu, &work[1], info, (ftnlen)1);
#line 265 "dlasd0.f"
	if (*info != 0) {
#line 266 "dlasd0.f"
	    return 0;
#line 267 "dlasd0.f"
	}
#line 268 "dlasd0.f"
	itemp = idxq + ic;
#line 269 "dlasd0.f"
	i__2 = nr;
#line 269 "dlasd0.f"
	for (j = 1; j <= i__2; ++j) {
#line 270 "dlasd0.f"
	    iwork[itemp + j - 1] = j;
#line 271 "dlasd0.f"
/* L20: */
#line 271 "dlasd0.f"
	}
#line 272 "dlasd0.f"
/* L30: */
#line 272 "dlasd0.f"
    }

/*     Now conquer each subproblem bottom-up. */

#line 276 "dlasd0.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {

/*        Find the first node LF and last node LL on the */
/*        current level LVL. */

#line 281 "dlasd0.f"
	if (lvl == 1) {
#line 282 "dlasd0.f"
	    lf = 1;
#line 283 "dlasd0.f"
	    ll = 1;
#line 284 "dlasd0.f"
	} else {
#line 285 "dlasd0.f"
	    i__1 = lvl - 1;
#line 285 "dlasd0.f"
	    lf = pow_ii(&c__2, &i__1);
#line 286 "dlasd0.f"
	    ll = (lf << 1) - 1;
#line 287 "dlasd0.f"
	}
#line 288 "dlasd0.f"
	i__1 = ll;
#line 288 "dlasd0.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 289 "dlasd0.f"
	    im1 = i__ - 1;
#line 290 "dlasd0.f"
	    ic = iwork[inode + im1];
#line 291 "dlasd0.f"
	    nl = iwork[ndiml + im1];
#line 292 "dlasd0.f"
	    nr = iwork[ndimr + im1];
#line 293 "dlasd0.f"
	    nlf = ic - nl;
#line 294 "dlasd0.f"
	    if (*sqre == 0 && i__ == ll) {
#line 295 "dlasd0.f"
		sqrei = *sqre;
#line 296 "dlasd0.f"
	    } else {
#line 297 "dlasd0.f"
		sqrei = 1;
#line 298 "dlasd0.f"
	    }
#line 299 "dlasd0.f"
	    idxqc = idxq + nlf - 1;
#line 300 "dlasd0.f"
	    alpha = d__[ic];
#line 301 "dlasd0.f"
	    beta = e[ic];
#line 302 "dlasd0.f"
	    dlasd1_(&nl, &nr, &sqrei, &d__[nlf], &alpha, &beta, &u[nlf + nlf *
		     u_dim1], ldu, &vt[nlf + nlf * vt_dim1], ldvt, &iwork[
		    idxqc], &iwork[iwk], &work[1], info);

/*        Report the possible convergence failure. */

#line 308 "dlasd0.f"
	    if (*info != 0) {
#line 309 "dlasd0.f"
		return 0;
#line 310 "dlasd0.f"
	    }
#line 311 "dlasd0.f"
/* L40: */
#line 311 "dlasd0.f"
	}
#line 312 "dlasd0.f"
/* L50: */
#line 312 "dlasd0.f"
    }

#line 314 "dlasd0.f"
    return 0;

/*     End of DLASD0 */

} /* dlasd0_ */

