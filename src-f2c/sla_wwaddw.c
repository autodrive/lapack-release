#line 1 "sla_wwaddw.f"
/* sla_wwaddw.f -- translated by f2c (version 20100827).
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

#line 1 "sla_wwaddw.f"
/* > \brief \b SLA_WWADDW adds a vector into a doubled-single vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_WWADDW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_wwa
ddw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_wwa
ddw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_wwa
ddw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_WWADDW( N, X, Y, W ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               X( * ), Y( * ), W( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLA_WWADDW adds a vector W into a doubled-single vector (X, Y). */
/* > */
/* >    This works for all extant IBM's hex and binary floating point */
/* >    arithmetics, but not for decimal. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >            The length of vectors X, Y, and W. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (N) */
/* >            The first part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array, dimension (N) */
/* >            The second part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is REAL array, dimension (N) */
/* >            The vector to be added. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sla_wwaddw__(integer *n, doublereal *x, doublereal *y, 
	doublereal *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal s;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 104 "sla_wwaddw.f"
    /* Parameter adjustments */
#line 104 "sla_wwaddw.f"
    --w;
#line 104 "sla_wwaddw.f"
    --y;
#line 104 "sla_wwaddw.f"
    --x;
#line 104 "sla_wwaddw.f"

#line 104 "sla_wwaddw.f"
    /* Function Body */
#line 104 "sla_wwaddw.f"
    i__1 = *n;
#line 104 "sla_wwaddw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 105 "sla_wwaddw.f"
	s = x[i__] + w[i__];
#line 106 "sla_wwaddw.f"
	s = s + s - s;
#line 107 "sla_wwaddw.f"
	y[i__] = x[i__] - s + w[i__] + y[i__];
#line 108 "sla_wwaddw.f"
	x[i__] = s;
#line 109 "sla_wwaddw.f"
/* L10: */
#line 109 "sla_wwaddw.f"
    }
#line 110 "sla_wwaddw.f"
    return 0;
} /* sla_wwaddw__ */

