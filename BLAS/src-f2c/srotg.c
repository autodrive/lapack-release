#line 1 "srotg.f"
/* srotg.f -- translated by f2c (version 20100827).
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

#line 1 "srotg.f"
/* Table of constant values */

static doublereal c_b2 = 1.;

/* > \brief \b SROTG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SROTG(SA,SB,C,S) */

/*       .. Scalar Arguments .. */
/*       REAL C,S,SA,SB */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SROTG construct givens plane rotation. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int srotg_(doublereal *sa, doublereal *sb, doublereal *c__, 
	doublereal *s)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal r__, z__, roe, scale;


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 66 "srotg.f"
    roe = *sb;
#line 67 "srotg.f"
    if (abs(*sa) > abs(*sb)) {
#line 67 "srotg.f"
	roe = *sa;
#line 67 "srotg.f"
    }
#line 68 "srotg.f"
    scale = abs(*sa) + abs(*sb);
#line 69 "srotg.f"
    if (scale == 0.) {
#line 70 "srotg.f"
	*c__ = 1.;
#line 71 "srotg.f"
	*s = 0.;
#line 72 "srotg.f"
	r__ = 0.;
#line 73 "srotg.f"
	z__ = 0.;
#line 74 "srotg.f"
    } else {
/* Computing 2nd power */
#line 75 "srotg.f"
	d__1 = *sa / scale;
/* Computing 2nd power */
#line 75 "srotg.f"
	d__2 = *sb / scale;
#line 75 "srotg.f"
	r__ = scale * sqrt(d__1 * d__1 + d__2 * d__2);
#line 76 "srotg.f"
	r__ = d_sign(&c_b2, &roe) * r__;
#line 77 "srotg.f"
	*c__ = *sa / r__;
#line 78 "srotg.f"
	*s = *sb / r__;
#line 79 "srotg.f"
	z__ = 1.;
#line 80 "srotg.f"
	if (abs(*sa) > abs(*sb)) {
#line 80 "srotg.f"
	    z__ = *s;
#line 80 "srotg.f"
	}
#line 81 "srotg.f"
	if (abs(*sb) >= abs(*sa) && *c__ != 0.) {
#line 81 "srotg.f"
	    z__ = 1. / *c__;
#line 81 "srotg.f"
	}
#line 82 "srotg.f"
    }
#line 83 "srotg.f"
    *sa = r__;
#line 84 "srotg.f"
    *sb = z__;
#line 85 "srotg.f"
    return 0;
} /* srotg_ */

