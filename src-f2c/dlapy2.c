#line 1 "dlapy2.f"
/* dlapy2.f -- translated by f2c (version 20100827).
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

#line 1 "dlapy2.f"
/* > \brief \b DLAPY2 returns sqrt(x2+y2). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAPY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLAPY2( X, Y ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   X, Y */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary */
/* > overflow. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION */
/* >          X and Y specify the values x and y. */
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
doublereal dlapy2_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static logical x_is_nan__, y_is_nan__;
    static doublereal w, z__, xabs, yabs;
    extern logical disnan_(doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*     .. Scalar Arguments .. */
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

#line 96 "dlapy2.f"
    x_is_nan__ = disnan_(x);
#line 97 "dlapy2.f"
    y_is_nan__ = disnan_(y);
#line 98 "dlapy2.f"
    if (x_is_nan__) {
#line 98 "dlapy2.f"
	ret_val = *x;
#line 98 "dlapy2.f"
    }
#line 99 "dlapy2.f"
    if (y_is_nan__) {
#line 99 "dlapy2.f"
	ret_val = *y;
#line 99 "dlapy2.f"
    }

#line 101 "dlapy2.f"
    if (! (x_is_nan__ || y_is_nan__)) {
#line 102 "dlapy2.f"
	xabs = abs(*x);
#line 103 "dlapy2.f"
	yabs = abs(*y);
#line 104 "dlapy2.f"
	w = max(xabs,yabs);
#line 105 "dlapy2.f"
	z__ = min(xabs,yabs);
#line 106 "dlapy2.f"
	if (z__ == 0.) {
#line 107 "dlapy2.f"
	    ret_val = w;
#line 108 "dlapy2.f"
	} else {
/* Computing 2nd power */
#line 109 "dlapy2.f"
	    d__1 = z__ / w;
#line 109 "dlapy2.f"
	    ret_val = w * sqrt(d__1 * d__1 + 1.);
#line 110 "dlapy2.f"
	}
#line 111 "dlapy2.f"
    }
#line 112 "dlapy2.f"
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */

