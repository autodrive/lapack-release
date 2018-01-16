#line 1 "../INSTALL/ilaver.f"
/* ../INSTALL/ilaver.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/ilaver.f"
/* > \brief \b ILAVER returns the LAPACK version. */
/* * */
/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH ) */

/*     INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  This subroutine returns the LAPACK version. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* >  \param[out] VERS_MAJOR */
/* >      return the lapack major version */
/* > */
/* >  \param[out] VERS_MINOR */
/* >      return the lapack minor version from the major version */
/* > */
/* >  \param[out] VERS_PATCH */
/* >      return the lapack patch version from the minor version */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int ilaver_(integer *vers_major__, integer *vers_minor__, 
	integer *vers_patch__)
{

/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*  ===================================================================== */

/*  ===================================================================== */
#line 60 "../INSTALL/ilaver.f"
    *vers_major__ = 3;
#line 61 "../INSTALL/ilaver.f"
    *vers_minor__ = 8;
#line 62 "../INSTALL/ilaver.f"
    *vers_patch__ = 0;
/*  ===================================================================== */

#line 65 "../INSTALL/ilaver.f"
    return 0;
} /* ilaver_ */

