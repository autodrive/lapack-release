#line 1 "slartgs.f"
/* slartgs.f -- translated by f2c (version 20100827).
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

#line 1 "slartgs.f"
/* > \brief \b SLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for t
he bidiagonal SVD problem. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARTGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartgs
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartgs
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartgs
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARTGS( X, Y, SIGMA, CS, SN ) */

/*       .. Scalar Arguments .. */
/*       REAL                    CS, SIGMA, SN, X, Y */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARTGS generates a plane rotation designed to introduce a bulge in */
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
/* >          X is REAL */
/* >          The (1,1) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is REAL */
/* >          The (1,2) entry of an upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is REAL */
/* >          The shift. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* >          CS is REAL */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* >          SN is REAL */
/* >          The sine of the rotation. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int slartgs_(doublereal *x, doublereal *y, doublereal *sigma,
	 doublereal *cs, doublereal *sn)
{
    static doublereal r__, s, w, z__;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal thresh;
    extern /* Subroutine */ int slartgp_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  =================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

#line 119 "slartgs.f"
    thresh = slamch_("E", (ftnlen)1);

/*     Compute the first column of B**T*B - SIGMA^2*I, up to a scale */
/*     factor. */

#line 124 "slartgs.f"
    if (*sigma == 0. && abs(*x) < thresh || abs(*x) == *sigma && *y == 0.) {
#line 126 "slartgs.f"
	z__ = 0.;
#line 127 "slartgs.f"
	w = 0.;
#line 128 "slartgs.f"
    } else if (*sigma == 0.) {
#line 129 "slartgs.f"
	if (*x >= 0.) {
#line 130 "slartgs.f"
	    z__ = *x;
#line 131 "slartgs.f"
	    w = *y;
#line 132 "slartgs.f"
	} else {
#line 133 "slartgs.f"
	    z__ = -(*x);
#line 134 "slartgs.f"
	    w = -(*y);
#line 135 "slartgs.f"
	}
#line 136 "slartgs.f"
    } else if (abs(*x) < thresh) {
#line 137 "slartgs.f"
	z__ = -(*sigma) * *sigma;
#line 138 "slartgs.f"
	w = 0.;
#line 139 "slartgs.f"
    } else {
#line 140 "slartgs.f"
	if (*x >= 0.) {
#line 141 "slartgs.f"
	    s = 1.;
#line 142 "slartgs.f"
	} else {
#line 143 "slartgs.f"
	    s = -1.;
#line 144 "slartgs.f"
	}
#line 145 "slartgs.f"
	z__ = s * (abs(*x) - *sigma) * (s + *sigma / *x);
#line 146 "slartgs.f"
	w = s * *y;
#line 147 "slartgs.f"
    }

/*     Generate the rotation. */
/*     CALL SLARTGP( Z, W, CS, SN, R ) might seem more natural; */
/*     reordering the arguments ensures that if Z = 0 then the rotation */
/*     is by PI/2. */

#line 154 "slartgs.f"
    slartgp_(&w, &z__, sn, cs, &r__);

#line 156 "slartgs.f"
    return 0;

/*     End SLARTGS */

} /* slartgs_ */

