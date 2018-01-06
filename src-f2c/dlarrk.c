#line 1 "dlarrk.f"
/* dlarrk.f -- translated by f2c (version 20100827).
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

#line 1 "dlarrk.f"
/* > \brief \b DLARRK computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrk.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrk.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrk.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARRK( N, IW, GL, GU, */
/*                           D, E2, PIVMIN, RELTOL, W, WERR, INFO) */

/*       .. Scalar Arguments .. */
/*       INTEGER   INFO, IW, N */
/*       DOUBLE PRECISION    PIVMIN, RELTOL, GL, GU, W, WERR */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARRK computes one eigenvalue of a symmetric tridiagonal */
/* > matrix T to suitable accuracy. This is an auxiliary code to be */
/* > called from DSTEMR. */
/* > */
/* > To avoid overflow, the matrix must be scaled so that its */
/* > largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
 */
/* > accuracy, it should not be much smaller than that. */
/* > */
/* > See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal */
/* > Matrix", Report CS41, Computer Science Dept., Stanford */
/* > University, July 21, 1966. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the tridiagonal matrix T.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] IW */
/* > \verbatim */
/* >          IW is INTEGER */
/* >          The index of the eigenvalues to be returned. */
/* > \endverbatim */
/* > */
/* > \param[in] GL */
/* > \verbatim */
/* >          GL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] GU */
/* > \verbatim */
/* >          GU is DOUBLE PRECISION */
/* >          An upper and a lower bound on the eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E2 */
/* > \verbatim */
/* >          E2 is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) squared off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* >          PIVMIN is DOUBLE PRECISION */
/* >          The minimum pivot allowed in the Sturm sequence for T. */
/* > \endverbatim */
/* > */
/* > \param[in] RELTOL */
/* > \verbatim */
/* >          RELTOL is DOUBLE PRECISION */
/* >          The minimum relative width of an interval.  When an interval */
/* >          is narrower than RELTOL times the larger (in */
/* >          magnitude) endpoint, then it is considered to be */
/* >          sufficiently small, i.e., converged.  Note: this should */
/* >          always be at least radix*machine epsilon. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] WERR */
/* > \verbatim */
/* >          WERR is DOUBLE PRECISION */
/* >          The error bound on the corresponding eigenvalue approximation */
/* >          in W. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:       Eigenvalue converged */
/* >          = -1:      Eigenvalue did NOT converge */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  FUDGE   DOUBLE PRECISION, default = 2 */
/* >          A "fudge factor" to widen the Gershgorin intervals. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlarrk_(integer *n, integer *iw, doublereal *gl, 
	doublereal *gu, doublereal *d__, doublereal *e2, doublereal *pivmin, 
	doublereal *reltol, doublereal *w, doublereal *werr, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, it;
    static doublereal mid, eps, tmp1, tmp2, left, atoli, right;
    static integer itmax;
    static doublereal rtoli, tnorm;
    extern doublereal dlamch_(char *, ftnlen);
    static integer negcnt;


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 184 "dlarrk.f"
    /* Parameter adjustments */
#line 184 "dlarrk.f"
    --e2;
#line 184 "dlarrk.f"
    --d__;
#line 184 "dlarrk.f"

#line 184 "dlarrk.f"
    /* Function Body */
#line 184 "dlarrk.f"
    if (*n <= 0) {
#line 185 "dlarrk.f"
	*info = 0;
#line 186 "dlarrk.f"
	return 0;
#line 187 "dlarrk.f"
    }

/*     Get machine constants */
#line 190 "dlarrk.f"
    eps = dlamch_("P", (ftnlen)1);
/* Computing MAX */
#line 192 "dlarrk.f"
    d__1 = abs(*gl), d__2 = abs(*gu);
#line 192 "dlarrk.f"
    tnorm = max(d__1,d__2);
#line 193 "dlarrk.f"
    rtoli = *reltol;
#line 194 "dlarrk.f"
    atoli = *pivmin * 4.;
#line 196 "dlarrk.f"
    itmax = (integer) ((log(tnorm + *pivmin) - log(*pivmin)) / log(2.)) + 2;
#line 199 "dlarrk.f"
    *info = -1;
#line 201 "dlarrk.f"
    left = *gl - tnorm * 2. * eps * *n - *pivmin * 4.;
#line 202 "dlarrk.f"
    right = *gu + tnorm * 2. * eps * *n + *pivmin * 4.;
#line 203 "dlarrk.f"
    it = 0;
#line 205 "dlarrk.f"
L10:

/*     Check if interval converged or maximum number of iterations reached */

#line 209 "dlarrk.f"
    tmp1 = (d__1 = right - left, abs(d__1));
/* Computing MAX */
#line 210 "dlarrk.f"
    d__1 = abs(right), d__2 = abs(left);
#line 210 "dlarrk.f"
    tmp2 = max(d__1,d__2);
/* Computing MAX */
#line 211 "dlarrk.f"
    d__1 = max(atoli,*pivmin), d__2 = rtoli * tmp2;
#line 211 "dlarrk.f"
    if (tmp1 < max(d__1,d__2)) {
#line 212 "dlarrk.f"
	*info = 0;
#line 213 "dlarrk.f"
	goto L30;
#line 214 "dlarrk.f"
    }
#line 215 "dlarrk.f"
    if (it > itmax) {
#line 215 "dlarrk.f"
	goto L30;
#line 215 "dlarrk.f"
    }

/*     Count number of negative pivots for mid-point */

#line 221 "dlarrk.f"
    ++it;
#line 222 "dlarrk.f"
    mid = (left + right) * .5;
#line 223 "dlarrk.f"
    negcnt = 0;
#line 224 "dlarrk.f"
    tmp1 = d__[1] - mid;
#line 225 "dlarrk.f"
    if (abs(tmp1) < *pivmin) {
#line 225 "dlarrk.f"
	tmp1 = -(*pivmin);
#line 225 "dlarrk.f"
    }
#line 227 "dlarrk.f"
    if (tmp1 <= 0.) {
#line 227 "dlarrk.f"
	++negcnt;
#line 227 "dlarrk.f"
    }

#line 230 "dlarrk.f"
    i__1 = *n;
#line 230 "dlarrk.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 231 "dlarrk.f"
	tmp1 = d__[i__] - e2[i__ - 1] / tmp1 - mid;
#line 232 "dlarrk.f"
	if (abs(tmp1) < *pivmin) {
#line 232 "dlarrk.f"
	    tmp1 = -(*pivmin);
#line 232 "dlarrk.f"
	}
#line 234 "dlarrk.f"
	if (tmp1 <= 0.) {
#line 234 "dlarrk.f"
	    ++negcnt;
#line 234 "dlarrk.f"
	}
#line 236 "dlarrk.f"
/* L20: */
#line 236 "dlarrk.f"
    }
#line 238 "dlarrk.f"
    if (negcnt >= *iw) {
#line 239 "dlarrk.f"
	right = mid;
#line 240 "dlarrk.f"
    } else {
#line 241 "dlarrk.f"
	left = mid;
#line 242 "dlarrk.f"
    }
#line 243 "dlarrk.f"
    goto L10;
#line 245 "dlarrk.f"
L30:

/*     Converged or maximum number of iterations reached */

#line 249 "dlarrk.f"
    *w = (left + right) * .5;
#line 250 "dlarrk.f"
    *werr = (d__1 = right - left, abs(d__1)) * .5;
#line 252 "dlarrk.f"
    return 0;

/*     End of DLARRK */

} /* dlarrk_ */

