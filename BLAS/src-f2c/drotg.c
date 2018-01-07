#line 1 "drotg.f"
/* drotg.f -- translated by f2c (version 20100827).
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

#line 1 "drotg.f"
/* Table of constant values */

static doublereal c_b2 = 1.;

/* > \brief \b DROTG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DROTG(DA,DB,C,S) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION C,DA,DB,S */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DROTG construct givens plane rotation. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int drotg_(doublereal *da, doublereal *db, doublereal *c__, 
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
#line 66 "drotg.f"
    roe = *db;
#line 67 "drotg.f"
    if (abs(*da) > abs(*db)) {
#line 67 "drotg.f"
	roe = *da;
#line 67 "drotg.f"
    }
#line 68 "drotg.f"
    scale = abs(*da) + abs(*db);
#line 69 "drotg.f"
    if (scale == 0.) {
#line 70 "drotg.f"
	*c__ = 1.;
#line 71 "drotg.f"
	*s = 0.;
#line 72 "drotg.f"
	r__ = 0.;
#line 73 "drotg.f"
	z__ = 0.;
#line 74 "drotg.f"
    } else {
/* Computing 2nd power */
#line 75 "drotg.f"
	d__1 = *da / scale;
/* Computing 2nd power */
#line 75 "drotg.f"
	d__2 = *db / scale;
#line 75 "drotg.f"
	r__ = scale * sqrt(d__1 * d__1 + d__2 * d__2);
#line 76 "drotg.f"
	r__ = d_sign(&c_b2, &roe) * r__;
#line 77 "drotg.f"
	*c__ = *da / r__;
#line 78 "drotg.f"
	*s = *db / r__;
#line 79 "drotg.f"
	z__ = 1.;
#line 80 "drotg.f"
	if (abs(*da) > abs(*db)) {
#line 80 "drotg.f"
	    z__ = *s;
#line 80 "drotg.f"
	}
#line 81 "drotg.f"
	if (abs(*db) >= abs(*da) && *c__ != 0.) {
#line 81 "drotg.f"
	    z__ = 1. / *c__;
#line 81 "drotg.f"
	}
#line 82 "drotg.f"
    }
#line 83 "drotg.f"
    *da = r__;
#line 84 "drotg.f"
    *db = z__;
#line 85 "drotg.f"
    return 0;
} /* drotg_ */

