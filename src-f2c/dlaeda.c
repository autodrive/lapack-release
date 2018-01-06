#line 1 "dlaeda.f"
/* dlaeda.f -- translated by f2c (version 20100827).
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

#line 1 "dlaeda.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b24 = 1.;
static doublereal c_b26 = 0.;

/* > \brief \b DLAEDA used by sstedc. Computes the Z vector determining the rank-one modification of the diago
nal matrix. Used when the original matrix is dense. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAEDA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaeda.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaeda.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaeda.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, */
/*                          GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            CURLVL, CURPBM, INFO, N, TLVLS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ), */
/*      $                   PRMPTR( * ), QPTR( * ) */
/*       DOUBLE PRECISION   GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAEDA computes the Z vector corresponding to the merge step in the */
/* > CURLVLth step of the merge process with TLVLS steps for the CURPBMth */
/* > problem. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The dimension of the symmetric tridiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TLVLS */
/* > \verbatim */
/* >          TLVLS is INTEGER */
/* >         The total number of merging levels in the overall divide and */
/* >         conquer tree. */
/* > \endverbatim */
/* > */
/* > \param[in] CURLVL */
/* > \verbatim */
/* >          CURLVL is INTEGER */
/* >         The current level in the overall merge routine, */
/* >         0 <= curlvl <= tlvls. */
/* > \endverbatim */
/* > */
/* > \param[in] CURPBM */
/* > \verbatim */
/* >          CURPBM is INTEGER */
/* >         The current problem in the current level in the overall */
/* >         merge routine (counting from upper left to lower right). */
/* > \endverbatim */
/* > */
/* > \param[in] PRMPTR */
/* > \verbatim */
/* >          PRMPTR is INTEGER array, dimension (N lg N) */
/* >         Contains a list of pointers which indicate where in PERM a */
/* >         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i) */
/* >         indicates the size of the permutation and incidentally the */
/* >         size of the full, non-deflated problem. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension (N lg N) */
/* >         Contains the permutations (from deflation and sorting) to be */
/* >         applied to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER array, dimension (N lg N) */
/* >         Contains a list of pointers which indicate where in GIVCOL a */
/* >         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i) */
/* >         indicates the number of Givens rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, dimension (2, N lg N) */
/* >         Each pair of numbers indicates a pair of columns to take place */
/* >         in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is DOUBLE PRECISION array, dimension (2, N lg N) */
/* >         Each number indicates the S value to be used in the */
/* >         corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (N**2) */
/* >         Contains the square eigenblocks from previous levels, the */
/* >         starting positions for blocks are given by QPTR. */
/* > \endverbatim */
/* > */
/* > \param[in] QPTR */
/* > \verbatim */
/* >          QPTR is INTEGER array, dimension (N+2) */
/* >         Contains a list of pointers which indicate where in Q an */
/* >         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates */
/* >         the size of the block. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (N) */
/* >         On output this vector contains the updating vector (the last */
/* >         row of the first sub-eigenvector matrix and the first row of */
/* >         the second sub-eigenvector matrix). */
/* > \endverbatim */
/* > */
/* > \param[out] ZTEMP */
/* > \verbatim */
/* >          ZTEMP is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlaeda_(integer *n, integer *tlvls, integer *curlvl, 
	integer *curpbm, integer *prmptr, integer *perm, integer *givptr, 
	integer *givcol, doublereal *givnum, doublereal *q, integer *qptr, 
	doublereal *z__, doublereal *ztemp, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k, mid, ptr;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer curr, bsiz1, bsiz2, psiz1, psiz2, zptr1;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), xerbla_(char *,
	     integer *, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 203 "dlaeda.f"
    /* Parameter adjustments */
#line 203 "dlaeda.f"
    --ztemp;
#line 203 "dlaeda.f"
    --z__;
#line 203 "dlaeda.f"
    --qptr;
#line 203 "dlaeda.f"
    --q;
#line 203 "dlaeda.f"
    givnum -= 3;
#line 203 "dlaeda.f"
    givcol -= 3;
#line 203 "dlaeda.f"
    --givptr;
#line 203 "dlaeda.f"
    --perm;
#line 203 "dlaeda.f"
    --prmptr;
#line 203 "dlaeda.f"

#line 203 "dlaeda.f"
    /* Function Body */
#line 203 "dlaeda.f"
    *info = 0;

#line 205 "dlaeda.f"
    if (*n < 0) {
#line 206 "dlaeda.f"
	*info = -1;
#line 207 "dlaeda.f"
    }
#line 208 "dlaeda.f"
    if (*info != 0) {
#line 209 "dlaeda.f"
	i__1 = -(*info);
#line 209 "dlaeda.f"
	xerbla_("DLAEDA", &i__1, (ftnlen)6);
#line 210 "dlaeda.f"
	return 0;
#line 211 "dlaeda.f"
    }

/*     Quick return if possible */

#line 215 "dlaeda.f"
    if (*n == 0) {
#line 215 "dlaeda.f"
	return 0;
#line 215 "dlaeda.f"
    }

/*     Determine location of first number in second half. */

#line 220 "dlaeda.f"
    mid = *n / 2 + 1;

/*     Gather last/first rows of appropriate eigenblocks into center of Z */

#line 224 "dlaeda.f"
    ptr = 1;

/*     Determine location of lowest level subproblem in the full storage */
/*     scheme */

#line 229 "dlaeda.f"
    i__1 = *curlvl - 1;
#line 229 "dlaeda.f"
    curr = ptr + *curpbm * pow_ii(&c__2, curlvl) + pow_ii(&c__2, &i__1) - 1;

/*     Determine size of these matrices.  We add HALF to the value of */
/*     the SQRT in case the machine underestimates one of these square */
/*     roots. */

#line 235 "dlaeda.f"
    bsiz1 = (integer) (sqrt((doublereal) (qptr[curr + 1] - qptr[curr])) + .5);
#line 236 "dlaeda.f"
    bsiz2 = (integer) (sqrt((doublereal) (qptr[curr + 2] - qptr[curr + 1])) + 
	    .5);
#line 237 "dlaeda.f"
    i__1 = mid - bsiz1 - 1;
#line 237 "dlaeda.f"
    for (k = 1; k <= i__1; ++k) {
#line 238 "dlaeda.f"
	z__[k] = 0.;
#line 239 "dlaeda.f"
/* L10: */
#line 239 "dlaeda.f"
    }
#line 240 "dlaeda.f"
    dcopy_(&bsiz1, &q[qptr[curr] + bsiz1 - 1], &bsiz1, &z__[mid - bsiz1], &
	    c__1);
#line 242 "dlaeda.f"
    dcopy_(&bsiz2, &q[qptr[curr + 1]], &bsiz2, &z__[mid], &c__1);
#line 243 "dlaeda.f"
    i__1 = *n;
#line 243 "dlaeda.f"
    for (k = mid + bsiz2; k <= i__1; ++k) {
#line 244 "dlaeda.f"
	z__[k] = 0.;
#line 245 "dlaeda.f"
/* L20: */
#line 245 "dlaeda.f"
    }

/*     Loop through remaining levels 1 -> CURLVL applying the Givens */
/*     rotations and permutation and then multiplying the center matrices */
/*     against the current Z. */

#line 251 "dlaeda.f"
    ptr = pow_ii(&c__2, tlvls) + 1;
#line 252 "dlaeda.f"
    i__1 = *curlvl - 1;
#line 252 "dlaeda.f"
    for (k = 1; k <= i__1; ++k) {
#line 253 "dlaeda.f"
	i__2 = *curlvl - k;
#line 253 "dlaeda.f"
	i__3 = *curlvl - k - 1;
#line 253 "dlaeda.f"
	curr = ptr + *curpbm * pow_ii(&c__2, &i__2) + pow_ii(&c__2, &i__3) - 
		1;
#line 254 "dlaeda.f"
	psiz1 = prmptr[curr + 1] - prmptr[curr];
#line 255 "dlaeda.f"
	psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
#line 256 "dlaeda.f"
	zptr1 = mid - psiz1;

/*       Apply Givens at CURR and CURR+1 */

#line 260 "dlaeda.f"
	i__2 = givptr[curr + 1] - 1;
#line 260 "dlaeda.f"
	for (i__ = givptr[curr]; i__ <= i__2; ++i__) {
#line 261 "dlaeda.f"
	    drot_(&c__1, &z__[zptr1 + givcol[(i__ << 1) + 1] - 1], &c__1, &
		    z__[zptr1 + givcol[(i__ << 1) + 2] - 1], &c__1, &givnum[(
		    i__ << 1) + 1], &givnum[(i__ << 1) + 2]);
#line 264 "dlaeda.f"
/* L30: */
#line 264 "dlaeda.f"
	}
#line 265 "dlaeda.f"
	i__2 = givptr[curr + 2] - 1;
#line 265 "dlaeda.f"
	for (i__ = givptr[curr + 1]; i__ <= i__2; ++i__) {
#line 266 "dlaeda.f"
	    drot_(&c__1, &z__[mid - 1 + givcol[(i__ << 1) + 1]], &c__1, &z__[
		    mid - 1 + givcol[(i__ << 1) + 2]], &c__1, &givnum[(i__ << 
		    1) + 1], &givnum[(i__ << 1) + 2]);
#line 269 "dlaeda.f"
/* L40: */
#line 269 "dlaeda.f"
	}
#line 270 "dlaeda.f"
	psiz1 = prmptr[curr + 1] - prmptr[curr];
#line 271 "dlaeda.f"
	psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
#line 272 "dlaeda.f"
	i__2 = psiz1 - 1;
#line 272 "dlaeda.f"
	for (i__ = 0; i__ <= i__2; ++i__) {
#line 273 "dlaeda.f"
	    ztemp[i__ + 1] = z__[zptr1 + perm[prmptr[curr] + i__] - 1];
#line 274 "dlaeda.f"
/* L50: */
#line 274 "dlaeda.f"
	}
#line 275 "dlaeda.f"
	i__2 = psiz2 - 1;
#line 275 "dlaeda.f"
	for (i__ = 0; i__ <= i__2; ++i__) {
#line 276 "dlaeda.f"
	    ztemp[psiz1 + i__ + 1] = z__[mid + perm[prmptr[curr + 1] + i__] - 
		    1];
#line 277 "dlaeda.f"
/* L60: */
#line 277 "dlaeda.f"
	}

/*        Multiply Blocks at CURR and CURR+1 */

/*        Determine size of these matrices.  We add HALF to the value of */
/*        the SQRT in case the machine underestimates one of these */
/*        square roots. */

#line 285 "dlaeda.f"
	bsiz1 = (integer) (sqrt((doublereal) (qptr[curr + 1] - qptr[curr])) + 
		.5);
#line 286 "dlaeda.f"
	bsiz2 = (integer) (sqrt((doublereal) (qptr[curr + 2] - qptr[curr + 1])
		) + .5);
#line 288 "dlaeda.f"
	if (bsiz1 > 0) {
#line 289 "dlaeda.f"
	    dgemv_("T", &bsiz1, &bsiz1, &c_b24, &q[qptr[curr]], &bsiz1, &
		    ztemp[1], &c__1, &c_b26, &z__[zptr1], &c__1, (ftnlen)1);
#line 291 "dlaeda.f"
	}
#line 292 "dlaeda.f"
	i__2 = psiz1 - bsiz1;
#line 292 "dlaeda.f"
	dcopy_(&i__2, &ztemp[bsiz1 + 1], &c__1, &z__[zptr1 + bsiz1], &c__1);
#line 294 "dlaeda.f"
	if (bsiz2 > 0) {
#line 295 "dlaeda.f"
	    dgemv_("T", &bsiz2, &bsiz2, &c_b24, &q[qptr[curr + 1]], &bsiz2, &
		    ztemp[psiz1 + 1], &c__1, &c_b26, &z__[mid], &c__1, (
		    ftnlen)1);
#line 297 "dlaeda.f"
	}
#line 298 "dlaeda.f"
	i__2 = psiz2 - bsiz2;
#line 298 "dlaeda.f"
	dcopy_(&i__2, &ztemp[psiz1 + bsiz2 + 1], &c__1, &z__[mid + bsiz2], &
		c__1);

#line 301 "dlaeda.f"
	i__2 = *tlvls - k;
#line 301 "dlaeda.f"
	ptr += pow_ii(&c__2, &i__2);
#line 302 "dlaeda.f"
/* L70: */
#line 302 "dlaeda.f"
    }

#line 304 "dlaeda.f"
    return 0;

/*     End of DLAEDA */

} /* dlaeda_ */

