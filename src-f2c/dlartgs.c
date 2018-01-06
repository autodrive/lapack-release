#line 1 "dlartgs.f"
/* dlartgs.f -- translated by f2c (version 20100827).
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

#line 1 "dlartgs.f"
/* > \brief \b DLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for t
he bidiagonal SVD problem. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARTGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartgs
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartgs
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartgs
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARTGS( X, Y, SIGMA, CS, SN ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION        CS, SIGMA, SN, X, Y */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARTGS generates a plane rotation designed to introduce a bulge in */
/* > Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD */
/* > problem. X and Y are the top-row entries, and SIGMA is the shift. */
/* > The computed CS and SN define a plane rotation satisfying */
/* > */
/* >    [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ], */
/* >    [ -SN  CS  ]     [    X * Y    ]     [ 0 ] */
/* > */
/* > with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the */
/* > rotation is by PI/2. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION */
/* >          The (1,1) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION */
/* >          The (1,2) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >          The shift. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* >          CS is DOUBLE PRECISION */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* >          SN is DOUBLE PRECISION */
/* >          The sine of the rotation. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlartgs_(doublereal *x, doublereal *y, doublereal *sigma,
	 doublereal *cs, doublereal *sn)
{
    static doublereal r__, s, w, z__;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal thresh;
    extern /* Subroutine */ int dlartgp_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  =================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

#line 116 "dlartgs.f"
    thresh = dlamch_("E", (ftnlen)1);

/*     Compute the first column of B**T*B - SIGMA^2*I, up to a scale */
/*     factor. */

#line 121 "dlartgs.f"
    if (*sigma == 0. && abs(*x) < thresh || abs(*x) == *sigma && *y == 0.) {
#line 123 "dlartgs.f"
	z__ = 0.;
#line 124 "dlartgs.f"
	w = 0.;
#line 125 "dlartgs.f"
    } else if (*sigma == 0.) {
#line 126 "dlartgs.f"
	if (*x >= 0.) {
#line 127 "dlartgs.f"
	    z__ = *x;
#line 128 "dlartgs.f"
	    w = *y;
#line 129 "dlartgs.f"
	} else {
#line 130 "dlartgs.f"
	    z__ = -(*x);
#line 131 "dlartgs.f"
	    w = -(*y);
#line 132 "dlartgs.f"
	}
#line 133 "dlartgs.f"
    } else if (abs(*x) < thresh) {
#line 134 "dlartgs.f"
	z__ = -(*sigma) * *sigma;
#line 135 "dlartgs.f"
	w = 0.;
#line 136 "dlartgs.f"
    } else {
#line 137 "dlartgs.f"
	if (*x >= 0.) {
#line 138 "dlartgs.f"
	    s = 1.;
#line 139 "dlartgs.f"
	} else {
#line 140 "dlartgs.f"
	    s = -1.;
#line 141 "dlartgs.f"
	}
#line 142 "dlartgs.f"
	z__ = s * (abs(*x) - *sigma) * (s + *sigma / *x);
#line 143 "dlartgs.f"
	w = s * *y;
#line 144 "dlartgs.f"
    }

/*     Generate the rotation. */
/*     CALL DLARTGP( Z, W, CS, SN, R ) might seem more natural; */
/*     reordering the arguments ensures that if Z = 0 then the rotation */
/*     is by PI/2. */

#line 151 "dlartgs.f"
    dlartgp_(&w, &z__, sn, cs, &r__);

#line 153 "dlartgs.f"
    return 0;

/*     End DLARTGS */

} /* dlartgs_ */

