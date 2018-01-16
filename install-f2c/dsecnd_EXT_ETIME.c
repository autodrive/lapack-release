#line 1 "../INSTALL/dsecnd_EXT_ETIME.f"
/* ../INSTALL/dsecnd_EXT_ETIME.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/dsecnd_EXT_ETIME.f"
/* > \brief \b DSECND  Using ETIME */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*      DOUBLE PRECISION FUNCTION DSECND( ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  DSECND returns the user time for a process in seconds. */
/* >  This version gets the time from the EXTERNAL system function ETIME. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
doublereal dsecnd_(void)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal t1;
    extern doublereal etime_(doublereal *);
    static doublereal tarray[2];


/*  -- LAPACK auxiliary routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */


/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 58 "../INSTALL/dsecnd_EXT_ETIME.f"
    t1 = etime_(tarray);
#line 59 "../INSTALL/dsecnd_EXT_ETIME.f"
    ret_val = tarray[0];
#line 60 "../INSTALL/dsecnd_EXT_ETIME.f"
    return ret_val;

/*     End of DSECND */

} /* dsecnd_ */

