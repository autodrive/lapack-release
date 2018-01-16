#line 1 "../INSTALL/dsecnd_EXT_ETIME_.f"
/* ../INSTALL/dsecnd_EXT_ETIME_.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/dsecnd_EXT_ETIME_.f"
/* > \brief \b DSECND  Using ETIME_ */

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
/* >  This version gets the time from the system function ETIME_. */
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
    extern doublereal etime___(doublereal *);
    static doublereal tarray[2];


/*  -- LAPACK auxiliary routine (version 3.4.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     February 2007 */
/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 55 "../INSTALL/dsecnd_EXT_ETIME_.f"
    t1 = etime___(tarray);
#line 56 "../INSTALL/dsecnd_EXT_ETIME_.f"
    ret_val = tarray[0];
#line 57 "../INSTALL/dsecnd_EXT_ETIME_.f"
    return ret_val;

/*     End of DSECND */

} /* dsecnd_ */

