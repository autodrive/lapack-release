#line 1 "slasd0.f"
/* slasd0.f -- translated by f2c (version 20100827).
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

#line 1 "slasd0.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__2 = 2;

/* > \brief \b SLASD0 computes the singular values of a real upper bidiagonal n-by-m matrix B with diagonal d 
and off-diagonal e. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASD0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd0.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd0.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd0.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDU, LDVT, N, SMLSIZ, SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Using a divide and conquer approach, SLASD0 computes the singular */
/* > value decomposition (SVD) of a real upper bidiagonal N-by-M */
/* > matrix B with diagonal D and offdiagonal E, where M = N + SQRE. */
/* > The algorithm computes orthogonal matrices U and VT such that */
/* > B = U * S * VT. The singular values S are overwritten on D. */
/* > */
/* > A related subroutine, SLASDA, computes only the singular values, */
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
/* >          D is REAL array, dimension (N) */
/* >         On entry D contains the main diagonal of the bidiagonal */
/* >         matrix. */
/* >         On exit D, if INFO = 0, contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (M-1) */
/* >         Contains the subdiagonal entries of the bidiagonal matrix. */
/* >         On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, dimension at least (LDQ, N) */
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
/* >          VT is REAL array, dimension at least (LDVT, M) */
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
/* >          IWORK is INTEGER array, dimension (8*N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (3*M**2+2*M) */
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

/* > \date November 2015 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasd0_(integer *n, integer *sqre, doublereal *d__, 
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
    extern /* Subroutine */ int slasd1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *), xerbla_(char *, integer *, ftnlen), slasdq_(char *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slasdt_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 182 "slasd0.f"
    /* Parameter adjustments */
#line 182 "slasd0.f"
    --d__;
#line 182 "slasd0.f"
    --e;
#line 182 "slasd0.f"
    u_dim1 = *ldu;
#line 182 "slasd0.f"
    u_offset = 1 + u_dim1;
#line 182 "slasd0.f"
    u -= u_offset;
#line 182 "slasd0.f"
    vt_dim1 = *ldvt;
#line 182 "slasd0.f"
    vt_offset = 1 + vt_dim1;
#line 182 "slasd0.f"
    vt -= vt_offset;
#line 182 "slasd0.f"
    --iwork;
#line 182 "slasd0.f"
    --work;
#line 182 "slasd0.f"

#line 182 "slasd0.f"
    /* Function Body */
#line 182 "slasd0.f"
    *info = 0;

#line 184 "slasd0.f"
    if (*n < 0) {
#line 185 "slasd0.f"
	*info = -1;
#line 186 "slasd0.f"
    } else if (*sqre < 0 || *sqre > 1) {
#line 187 "slasd0.f"
	*info = -2;
#line 188 "slasd0.f"
    }

#line 190 "slasd0.f"
    m = *n + *sqre;

#line 192 "slasd0.f"
    if (*ldu < *n) {
#line 193 "slasd0.f"
	*info = -6;
#line 194 "slasd0.f"
    } else if (*ldvt < m) {
#line 195 "slasd0.f"
	*info = -8;
#line 196 "slasd0.f"
    } else if (*smlsiz < 3) {
#line 197 "slasd0.f"
	*info = -9;
#line 198 "slasd0.f"
    }
#line 199 "slasd0.f"
    if (*info != 0) {
#line 200 "slasd0.f"
	i__1 = -(*info);
#line 200 "slasd0.f"
	xerbla_("SLASD0", &i__1, (ftnlen)6);
#line 201 "slasd0.f"
	return 0;
#line 202 "slasd0.f"
    }

/*     If the input matrix is too small, call SLASDQ to find the SVD. */

#line 206 "slasd0.f"
    if (*n <= *smlsiz) {
#line 207 "slasd0.f"
	slasdq_("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], 
		ldvt, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info, (
		ftnlen)1);
#line 209 "slasd0.f"
	return 0;
#line 210 "slasd0.f"
    }

/*     Set up the computation tree. */

#line 214 "slasd0.f"
    inode = 1;
#line 215 "slasd0.f"
    ndiml = inode + *n;
#line 216 "slasd0.f"
    ndimr = ndiml + *n;
#line 217 "slasd0.f"
    idxq = ndimr + *n;
#line 218 "slasd0.f"
    iwk = idxq + *n;
#line 219 "slasd0.f"
    slasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

/*     For the nodes on bottom level of the tree, solve */
/*     their subproblems by SLASDQ. */

#line 225 "slasd0.f"
    ndb1 = (nd + 1) / 2;
#line 226 "slasd0.f"
    ncc = 0;
#line 227 "slasd0.f"
    i__1 = nd;
#line 227 "slasd0.f"
    for (i__ = ndb1; i__ <= i__1; ++i__) {

/*     IC : center row of each node */
/*     NL : number of rows of left  subproblem */
/*     NR : number of rows of right subproblem */
/*     NLF: starting row of the left   subproblem */
/*     NRF: starting row of the right  subproblem */

#line 235 "slasd0.f"
	i1 = i__ - 1;
#line 236 "slasd0.f"
	ic = iwork[inode + i1];
#line 237 "slasd0.f"
	nl = iwork[ndiml + i1];
#line 238 "slasd0.f"
	nlp1 = nl + 1;
#line 239 "slasd0.f"
	nr = iwork[ndimr + i1];
#line 240 "slasd0.f"
	nrp1 = nr + 1;
#line 241 "slasd0.f"
	nlf = ic - nl;
#line 242 "slasd0.f"
	nrf = ic + 1;
#line 243 "slasd0.f"
	sqrei = 1;
#line 244 "slasd0.f"
	slasdq_("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[
		nlf + nlf * vt_dim1], ldvt, &u[nlf + nlf * u_dim1], ldu, &u[
		nlf + nlf * u_dim1], ldu, &work[1], info, (ftnlen)1);
#line 247 "slasd0.f"
	if (*info != 0) {
#line 248 "slasd0.f"
	    return 0;
#line 249 "slasd0.f"
	}
#line 250 "slasd0.f"
	itemp = idxq + nlf - 2;
#line 251 "slasd0.f"
	i__2 = nl;
#line 251 "slasd0.f"
	for (j = 1; j <= i__2; ++j) {
#line 252 "slasd0.f"
	    iwork[itemp + j] = j;
#line 253 "slasd0.f"
/* L10: */
#line 253 "slasd0.f"
	}
#line 254 "slasd0.f"
	if (i__ == nd) {
#line 255 "slasd0.f"
	    sqrei = *sqre;
#line 256 "slasd0.f"
	} else {
#line 257 "slasd0.f"
	    sqrei = 1;
#line 258 "slasd0.f"
	}
#line 259 "slasd0.f"
	nrp1 = nr + sqrei;
#line 260 "slasd0.f"
	slasdq_("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[
		nrf + nrf * vt_dim1], ldvt, &u[nrf + nrf * u_dim1], ldu, &u[
		nrf + nrf * u_dim1], ldu, &work[1], info, (ftnlen)1);
#line 263 "slasd0.f"
	if (*info != 0) {
#line 264 "slasd0.f"
	    return 0;
#line 265 "slasd0.f"
	}
#line 266 "slasd0.f"
	itemp = idxq + ic;
#line 267 "slasd0.f"
	i__2 = nr;
#line 267 "slasd0.f"
	for (j = 1; j <= i__2; ++j) {
#line 268 "slasd0.f"
	    iwork[itemp + j - 1] = j;
#line 269 "slasd0.f"
/* L20: */
#line 269 "slasd0.f"
	}
#line 270 "slasd0.f"
/* L30: */
#line 270 "slasd0.f"
    }

/*     Now conquer each subproblem bottom-up. */

#line 274 "slasd0.f"
    for (lvl = nlvl; lvl >= 1; --lvl) {

/*        Find the first node LF and last node LL on the */
/*        current level LVL. */

#line 279 "slasd0.f"
	if (lvl == 1) {
#line 280 "slasd0.f"
	    lf = 1;
#line 281 "slasd0.f"
	    ll = 1;
#line 282 "slasd0.f"
	} else {
#line 283 "slasd0.f"
	    i__1 = lvl - 1;
#line 283 "slasd0.f"
	    lf = pow_ii(&c__2, &i__1);
#line 284 "slasd0.f"
	    ll = (lf << 1) - 1;
#line 285 "slasd0.f"
	}
#line 286 "slasd0.f"
	i__1 = ll;
#line 286 "slasd0.f"
	for (i__ = lf; i__ <= i__1; ++i__) {
#line 287 "slasd0.f"
	    im1 = i__ - 1;
#line 288 "slasd0.f"
	    ic = iwork[inode + im1];
#line 289 "slasd0.f"
	    nl = iwork[ndiml + im1];
#line 290 "slasd0.f"
	    nr = iwork[ndimr + im1];
#line 291 "slasd0.f"
	    nlf = ic - nl;
#line 292 "slasd0.f"
	    if (*sqre == 0 && i__ == ll) {
#line 293 "slasd0.f"
		sqrei = *sqre;
#line 294 "slasd0.f"
	    } else {
#line 295 "slasd0.f"
		sqrei = 1;
#line 296 "slasd0.f"
	    }
#line 297 "slasd0.f"
	    idxqc = idxq + nlf - 1;
#line 298 "slasd0.f"
	    alpha = d__[ic];
#line 299 "slasd0.f"
	    beta = e[ic];
#line 300 "slasd0.f"
	    slasd1_(&nl, &nr, &sqrei, &d__[nlf], &alpha, &beta, &u[nlf + nlf *
		     u_dim1], ldu, &vt[nlf + nlf * vt_dim1], ldvt, &iwork[
		    idxqc], &iwork[iwk], &work[1], info);

/*     Report the possible convergence failure. */

#line 306 "slasd0.f"
	    if (*info != 0) {
#line 307 "slasd0.f"
		return 0;
#line 308 "slasd0.f"
	    }
#line 309 "slasd0.f"
/* L40: */
#line 309 "slasd0.f"
	}
#line 310 "slasd0.f"
/* L50: */
#line 310 "slasd0.f"
    }

#line 312 "slasd0.f"
    return 0;

/*     End of SLASD0 */

} /* slasd0_ */

