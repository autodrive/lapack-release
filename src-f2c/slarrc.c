#line 1 "slarrc.f"
/* slarrc.f -- translated by f2c (version 20100827).
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

#line 1 "slarrc.f"
/* > \brief \b SLARRC computes the number of eigenvalues of the symmetric tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARRC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARRC( JOBT, N, VL, VU, D, E, PIVMIN, */
/*                                   EIGCNT, LCNT, RCNT, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBT */
/*       INTEGER            EIGCNT, INFO, LCNT, N, RCNT */
/*       REAL               PIVMIN, VL, VU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ) */
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
/* >          VL is REAL */
/* >          The lower bound for the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >          VU is REAL */
/* >          The upper bound for the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T. */
/* >          JOBT = 'L': The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* >          JOBT = 'T': The N-1 offdiagonal elements of the matrix T. */
/* >          JOBT = 'L': The N-1 offdiagonal elements of the matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is REAL */
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

/* > \date June 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int slarrc_(char *jobt, integer *n, doublereal *vl, 
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


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 172 "slarrc.f"
    /* Parameter adjustments */
#line 172 "slarrc.f"
    --e;
#line 172 "slarrc.f"
    --d__;
#line 172 "slarrc.f"

#line 172 "slarrc.f"
    /* Function Body */
#line 172 "slarrc.f"
    *info = 0;

/*     Quick return if possible */

#line 176 "slarrc.f"
    if (*n <= 0) {
#line 177 "slarrc.f"
	return 0;
#line 178 "slarrc.f"
    }

#line 180 "slarrc.f"
    *lcnt = 0;
#line 181 "slarrc.f"
    *rcnt = 0;
#line 182 "slarrc.f"
    *eigcnt = 0;
#line 183 "slarrc.f"
    matt = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 186 "slarrc.f"
    if (matt) {
/*        Sturm sequence count on T */
#line 188 "slarrc.f"
	lpivot = d__[1] - *vl;
#line 189 "slarrc.f"
	rpivot = d__[1] - *vu;
#line 190 "slarrc.f"
	if (lpivot <= 0.) {
#line 191 "slarrc.f"
	    ++(*lcnt);
#line 192 "slarrc.f"
	}
#line 193 "slarrc.f"
	if (rpivot <= 0.) {
#line 194 "slarrc.f"
	    ++(*rcnt);
#line 195 "slarrc.f"
	}
#line 196 "slarrc.f"
	i__1 = *n - 1;
#line 196 "slarrc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 197 "slarrc.f"
	    d__1 = e[i__];
#line 197 "slarrc.f"
	    tmp = d__1 * d__1;
#line 198 "slarrc.f"
	    lpivot = d__[i__ + 1] - *vl - tmp / lpivot;
#line 199 "slarrc.f"
	    rpivot = d__[i__ + 1] - *vu - tmp / rpivot;
#line 200 "slarrc.f"
	    if (lpivot <= 0.) {
#line 201 "slarrc.f"
		++(*lcnt);
#line 202 "slarrc.f"
	    }
#line 203 "slarrc.f"
	    if (rpivot <= 0.) {
#line 204 "slarrc.f"
		++(*rcnt);
#line 205 "slarrc.f"
	    }
#line 206 "slarrc.f"
/* L10: */
#line 206 "slarrc.f"
	}
#line 207 "slarrc.f"
    } else {
/*        Sturm sequence count on L D L^T */
#line 209 "slarrc.f"
	sl = -(*vl);
#line 210 "slarrc.f"
	su = -(*vu);
#line 211 "slarrc.f"
	i__1 = *n - 1;
#line 211 "slarrc.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "slarrc.f"
	    lpivot = d__[i__] + sl;
#line 213 "slarrc.f"
	    rpivot = d__[i__] + su;
#line 214 "slarrc.f"
	    if (lpivot <= 0.) {
#line 215 "slarrc.f"
		++(*lcnt);
#line 216 "slarrc.f"
	    }
#line 217 "slarrc.f"
	    if (rpivot <= 0.) {
#line 218 "slarrc.f"
		++(*rcnt);
#line 219 "slarrc.f"
	    }
#line 220 "slarrc.f"
	    tmp = e[i__] * d__[i__] * e[i__];

#line 222 "slarrc.f"
	    tmp2 = tmp / lpivot;
#line 223 "slarrc.f"
	    if (tmp2 == 0.) {
#line 224 "slarrc.f"
		sl = tmp - *vl;
#line 225 "slarrc.f"
	    } else {
#line 226 "slarrc.f"
		sl = sl * tmp2 - *vl;
#line 227 "slarrc.f"
	    }

#line 229 "slarrc.f"
	    tmp2 = tmp / rpivot;
#line 230 "slarrc.f"
	    if (tmp2 == 0.) {
#line 231 "slarrc.f"
		su = tmp - *vu;
#line 232 "slarrc.f"
	    } else {
#line 233 "slarrc.f"
		su = su * tmp2 - *vu;
#line 234 "slarrc.f"
	    }
#line 235 "slarrc.f"
/* L20: */
#line 235 "slarrc.f"
	}
#line 236 "slarrc.f"
	lpivot = d__[*n] + sl;
#line 237 "slarrc.f"
	rpivot = d__[*n] + su;
#line 238 "slarrc.f"
	if (lpivot <= 0.) {
#line 239 "slarrc.f"
	    ++(*lcnt);
#line 240 "slarrc.f"
	}
#line 241 "slarrc.f"
	if (rpivot <= 0.) {
#line 242 "slarrc.f"
	    ++(*rcnt);
#line 243 "slarrc.f"
	}
#line 244 "slarrc.f"
    }
#line 245 "slarrc.f"
    *eigcnt = *rcnt - *lcnt;
#line 247 "slarrc.f"
    return 0;

/*     end of SLARRC */

} /* slarrc_ */

