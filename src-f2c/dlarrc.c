#line 1 "dlarrc.f"
/* dlarrc.f -- translated by f2c (version 20100827).
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

#line 1 "dlarrc.f"
/* > \brief \b DLARRC computes the number of eigenvalues of the symmetric tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRC( JOBT, N, VL, VU, D, E, PIVMIN, */
/*                                   EIGCNT, LCNT, RCNT, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBT */
/*       INTEGER            EIGCNT, INFO, LCNT, N, RCNT */
/*       DOUBLE PRECISION   PIVMIN, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Find the number of eigenvalues of the symmetric tridiagonal matrix T */
/* > that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T */
/* > if JOBT = 'L'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBT */
/* > \verbatim */
/* >          JOBT is CHARACTER*1 */
/* >          = 'T':  Compute Sturm count for matrix T. */
/* >          = 'L':  Compute Sturm count for matrix L D L^T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix. N > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is DOUBLE PRECISION */
/* >          The lower and upper bounds for the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T. */
/* >          JOBT = 'L': The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          JOBT = 'T': The N-1 offdiagonal elements of the matrix T. */
/* >          JOBT = 'L': The N-1 offdiagonal elements of the matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum pivot in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[out] EIGCNT */
/* > \verbatim */
/* >          EIGCNT is INTEGER */
/* >          The number of eigenvalues of the symmetric tridiagonal matrix T */
/* >          that are in the interval (VL,VU] */
/* > \endverbatim */
/* > */
/* > \param[out] LCNT */
/* > \verbatim */
/* >          LCNT is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] RCNT */
/* > \verbatim */
/* >          RCNT is INTEGER */
/* >          The left and right negcounts of the interval. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
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
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlarrc_(char *jobt, integer *n, doublereal *vl, 
	doublereal *vu, doublereal *d__, doublereal *e, doublereal *pivmin, 
	integer *eigcnt, integer *lcnt, integer *rcnt, integer *info, ftnlen 
	jobt_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal sl, su, tmp, tmp2;
    static logical matt;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal lpivot, rpivot;


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
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 171 "dlarrc.f"
    /* Parameter adjustments */
#line 171 "dlarrc.f"
    --e;
#line 171 "dlarrc.f"
    --d__;
#line 171 "dlarrc.f"

#line 171 "dlarrc.f"
    /* Function Body */
#line 171 "dlarrc.f"
    *info = 0;
#line 172 "dlarrc.f"
    *lcnt = 0;
#line 173 "dlarrc.f"
    *rcnt = 0;
#line 174 "dlarrc.f"
    *eigcnt = 0;
#line 175 "dlarrc.f"
    matt = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 178 "dlarrc.f"
    if (matt) {
/*        Sturm sequence count on T */
#line 180 "dlarrc.f"
	lpivot = d__[1] - *vl;
#line 181 "dlarrc.f"
	rpivot = d__[1] - *vu;
#line 182 "dlarrc.f"
	if (lpivot <= 0.) {
#line 183 "dlarrc.f"
	    ++(*lcnt);
#line 184 "dlarrc.f"
	}
#line 185 "dlarrc.f"
	if (rpivot <= 0.) {
#line 186 "dlarrc.f"
	    ++(*rcnt);
#line 187 "dlarrc.f"
	}
#line 188 "dlarrc.f"
	i__1 = *n - 1;
#line 188 "dlarrc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 189 "dlarrc.f"
	    d__1 = e[i__];
#line 189 "dlarrc.f"
	    tmp = d__1 * d__1;
#line 190 "dlarrc.f"
	    lpivot = d__[i__ + 1] - *vl - tmp / lpivot;
#line 191 "dlarrc.f"
	    rpivot = d__[i__ + 1] - *vu - tmp / rpivot;
#line 192 "dlarrc.f"
	    if (lpivot <= 0.) {
#line 193 "dlarrc.f"
		++(*lcnt);
#line 194 "dlarrc.f"
	    }
#line 195 "dlarrc.f"
	    if (rpivot <= 0.) {
#line 196 "dlarrc.f"
		++(*rcnt);
#line 197 "dlarrc.f"
	    }
#line 198 "dlarrc.f"
/* L10: */
#line 198 "dlarrc.f"
	}
#line 199 "dlarrc.f"
    } else {
/*        Sturm sequence count on L D L^T */
#line 201 "dlarrc.f"
	sl = -(*vl);
#line 202 "dlarrc.f"
	su = -(*vu);
#line 203 "dlarrc.f"
	i__1 = *n - 1;
#line 203 "dlarrc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "dlarrc.f"
	    lpivot = d__[i__] + sl;
#line 205 "dlarrc.f"
	    rpivot = d__[i__] + su;
#line 206 "dlarrc.f"
	    if (lpivot <= 0.) {
#line 207 "dlarrc.f"
		++(*lcnt);
#line 208 "dlarrc.f"
	    }
#line 209 "dlarrc.f"
	    if (rpivot <= 0.) {
#line 210 "dlarrc.f"
		++(*rcnt);
#line 211 "dlarrc.f"
	    }
#line 212 "dlarrc.f"
	    tmp = e[i__] * d__[i__] * e[i__];

#line 214 "dlarrc.f"
	    tmp2 = tmp / lpivot;
#line 215 "dlarrc.f"
	    if (tmp2 == 0.) {
#line 216 "dlarrc.f"
		sl = tmp - *vl;
#line 217 "dlarrc.f"
	    } else {
#line 218 "dlarrc.f"
		sl = sl * tmp2 - *vl;
#line 219 "dlarrc.f"
	    }

#line 221 "dlarrc.f"
	    tmp2 = tmp / rpivot;
#line 222 "dlarrc.f"
	    if (tmp2 == 0.) {
#line 223 "dlarrc.f"
		su = tmp - *vu;
#line 224 "dlarrc.f"
	    } else {
#line 225 "dlarrc.f"
		su = su * tmp2 - *vu;
#line 226 "dlarrc.f"
	    }
#line 227 "dlarrc.f"
/* L20: */
#line 227 "dlarrc.f"
	}
#line 228 "dlarrc.f"
	lpivot = d__[*n] + sl;
#line 229 "dlarrc.f"
	rpivot = d__[*n] + su;
#line 230 "dlarrc.f"
	if (lpivot <= 0.) {
#line 231 "dlarrc.f"
	    ++(*lcnt);
#line 232 "dlarrc.f"
	}
#line 233 "dlarrc.f"
	if (rpivot <= 0.) {
#line 234 "dlarrc.f"
	    ++(*rcnt);
#line 235 "dlarrc.f"
	}
#line 236 "dlarrc.f"
    }
#line 237 "dlarrc.f"
    *eigcnt = *rcnt - *lcnt;
#line 239 "dlarrc.f"
    return 0;

/*     end of DLARRC */

} /* dlarrc_ */

