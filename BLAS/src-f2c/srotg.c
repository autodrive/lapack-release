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

/*  Arguments: */
/*  ========== */

/* > \param[in] SA */
/* > \verbatim */
/* >          SA is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] SB */
/* > \verbatim */
/* >          SB is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

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


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 89 "srotg.f"
    roe = *sb;
#line 90 "srotg.f"
    if (abs(*sa) > abs(*sb)) {
#line 90 "srotg.f"
	roe = *sa;
#line 90 "srotg.f"
    }
#line 91 "srotg.f"
    scale = abs(*sa) + abs(*sb);
#line 92 "srotg.f"
    if (scale == 0.) {
#line 93 "srotg.f"
	*c__ = 1.;
#line 94 "srotg.f"
	*s = 0.;
#line 95 "srotg.f"
	r__ = 0.;
#line 96 "srotg.f"
	z__ = 0.;
#line 97 "srotg.f"
    } else {
/* Computing 2nd power */
#line 98 "srotg.f"
	d__1 = *sa / scale;
/* Computing 2nd power */
#line 98 "srotg.f"
	d__2 = *sb / scale;
#line 98 "srotg.f"
	r__ = scale * sqrt(d__1 * d__1 + d__2 * d__2);
#line 99 "srotg.f"
	r__ = d_sign(&c_b2, &roe) * r__;
#line 100 "srotg.f"
	*c__ = *sa / r__;
#line 101 "srotg.f"
	*s = *sb / r__;
#line 102 "srotg.f"
	z__ = 1.;
#line 103 "srotg.f"
	if (abs(*sa) > abs(*sb)) {
#line 103 "srotg.f"
	    z__ = *s;
#line 103 "srotg.f"
	}
#line 104 "srotg.f"
	if (abs(*sb) >= abs(*sa) && *c__ != 0.) {
#line 104 "srotg.f"
	    z__ = 1. / *c__;
#line 104 "srotg.f"
	}
#line 105 "srotg.f"
    }
#line 106 "srotg.f"
    *sa = r__;
#line 107 "srotg.f"
    *sb = z__;
#line 108 "srotg.f"
    return 0;
} /* srotg_ */

