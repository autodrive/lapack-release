#line 1 "slapy3.f"
/* slapy3.f -- translated by f2c (version 20100827).
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

#line 1 "slapy3.f"
/* > \brief \b SLAPY3 returns sqrt(x2+y2+z2). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAPY3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLAPY3( X, Y, Z ) */

/*       .. Scalar Arguments .. */
/*       REAL               X, Y, Z */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause */
/* > unnecessary overflow. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL */
/* >          X, Y and Z specify the values x, y and z. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
doublereal slapy3_(doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal w, xabs, yabs, zabs;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 94 "slapy3.f"
    xabs = abs(*x);
#line 95 "slapy3.f"
    yabs = abs(*y);
#line 96 "slapy3.f"
    zabs = abs(*z__);
/* Computing MAX */
#line 97 "slapy3.f"
    d__1 = max(xabs,yabs);
#line 97 "slapy3.f"
    w = max(d__1,zabs);
#line 98 "slapy3.f"
    if (w == 0.) {
/*     W can be zero for max(0,nan,0) */
/*     adding all three entries together will make sure */
/*     NaN will not disappear. */
#line 102 "slapy3.f"
	ret_val = xabs + yabs + zabs;
#line 103 "slapy3.f"
    } else {
/* Computing 2nd power */
#line 104 "slapy3.f"
	d__1 = xabs / w;
/* Computing 2nd power */
#line 104 "slapy3.f"
	d__2 = yabs / w;
/* Computing 2nd power */
#line 104 "slapy3.f"
	d__3 = zabs / w;
#line 104 "slapy3.f"
	ret_val = w * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
#line 106 "slapy3.f"
    }
#line 107 "slapy3.f"
    return ret_val;

/*     End of SLAPY3 */

} /* slapy3_ */

