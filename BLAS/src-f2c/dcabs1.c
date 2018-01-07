#line 1 "dcabs1.f"
/* dcabs1.f -- translated by f2c (version 20100827).
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

#line 1 "dcabs1.f"
/* > \brief \b DCABS1 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DCABS1(Z) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 Z */
/*       .. */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DCABS1 computes absolute value of a double complex number */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup double_blas_level1 */

/*  ===================================================================== */
doublereal dcabs1_(doublecomplex *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. */
/*  ===================================================================== */

/*     .. Intrinsic Functions .. */

#line 56 "dcabs1.f"
    ret_val = (d__1 = z__->r, abs(d__1)) + (d__2 = d_imag(z__), abs(d__2));
#line 57 "dcabs1.f"
    return ret_val;
} /* dcabs1_ */

