#line 1 "slasdt.f"
/* slasdt.f -- translated by f2c (version 20100827).
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

#line 1 "slasdt.f"
/* > \brief \b SLASDT creates a tree of subproblems for bidiagonal divide and conquer. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASDT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasdt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasdt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasdt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LVL, MSUB, N, ND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            INODE( * ), NDIML( * ), NDIMR( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASDT creates a tree of subproblems for bidiagonal divide and */
/* > conquer. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          On entry, the number of diagonal elements of the */
/* >          bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] LVL */
/* > \verbatim */
/* >          LVL is INTEGER */
/* >          On exit, the number of levels on the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] ND */
/* > \verbatim */
/* >          ND is INTEGER */
/* >          On exit, the number of nodes on the tree. */
/* > \endverbatim */
/* > */
/* > \param[out] INODE */
/* > \verbatim */
/* >          INODE is INTEGER array, dimension ( N ) */
/* >          On exit, centers of subproblems. */
/* > \endverbatim */
/* > */
/* > \param[out] NDIML */
/* > \verbatim */
/* >          NDIML is INTEGER array, dimension ( N ) */
/* >          On exit, row dimensions of left children. */
/* > \endverbatim */
/* > */
/* > \param[out] NDIMR */
/* > \verbatim */
/* >          NDIMR is INTEGER array, dimension ( N ) */
/* >          On exit, row dimensions of right children. */
/* > \endverbatim */
/* > */
/* > \param[in] MSUB */
/* > \verbatim */
/* >          MSUB is INTEGER */
/* >          On entry, the maximum row dimension each subproblem at the */
/* >          bottom of the tree can be of. */
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
/* Subroutine */ int slasdt_(integer *n, integer *lvl, integer *nd, integer *
	inode, integer *ndiml, integer *ndimr, integer *msub)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, il, ir, maxn;
    static doublereal temp;
    static integer nlvl, llst, ncrnt;


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
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Find the number of levels on the tree. */

#line 137 "slasdt.f"
    /* Parameter adjustments */
#line 137 "slasdt.f"
    --ndimr;
#line 137 "slasdt.f"
    --ndiml;
#line 137 "slasdt.f"
    --inode;
#line 137 "slasdt.f"

#line 137 "slasdt.f"
    /* Function Body */
#line 137 "slasdt.f"
    maxn = max(1,*n);
#line 138 "slasdt.f"
    temp = log((doublereal) maxn / (doublereal) (*msub + 1)) / log(2.);
#line 139 "slasdt.f"
    *lvl = (integer) temp + 1;

#line 141 "slasdt.f"
    i__ = *n / 2;
#line 142 "slasdt.f"
    inode[1] = i__ + 1;
#line 143 "slasdt.f"
    ndiml[1] = i__;
#line 144 "slasdt.f"
    ndimr[1] = *n - i__ - 1;
#line 145 "slasdt.f"
    il = 0;
#line 146 "slasdt.f"
    ir = 1;
#line 147 "slasdt.f"
    llst = 1;
#line 148 "slasdt.f"
    i__1 = *lvl - 1;
#line 148 "slasdt.f"
    for (nlvl = 1; nlvl <= i__1; ++nlvl) {

/*        Constructing the tree at (NLVL+1)-st level. The number of */
/*        nodes created on this level is LLST * 2. */

#line 153 "slasdt.f"
	i__2 = llst - 1;
#line 153 "slasdt.f"
	for (i__ = 0; i__ <= i__2; ++i__) {
#line 154 "slasdt.f"
	    il += 2;
#line 155 "slasdt.f"
	    ir += 2;
#line 156 "slasdt.f"
	    ncrnt = llst + i__;
#line 157 "slasdt.f"
	    ndiml[il] = ndiml[ncrnt] / 2;
#line 158 "slasdt.f"
	    ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
#line 159 "slasdt.f"
	    inode[il] = inode[ncrnt] - ndimr[il] - 1;
#line 160 "slasdt.f"
	    ndiml[ir] = ndimr[ncrnt] / 2;
#line 161 "slasdt.f"
	    ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
#line 162 "slasdt.f"
	    inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
#line 163 "slasdt.f"
/* L10: */
#line 163 "slasdt.f"
	}
#line 164 "slasdt.f"
	llst <<= 1;
#line 165 "slasdt.f"
/* L20: */
#line 165 "slasdt.f"
    }
#line 166 "slasdt.f"
    *nd = (llst << 1) - 1;

#line 168 "slasdt.f"
    return 0;

/*     End of SLASDT */

} /* slasdt_ */

