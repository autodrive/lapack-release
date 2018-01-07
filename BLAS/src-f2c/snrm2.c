#line 1 "snrm2.f"
/* snrm2.f -- translated by f2c (version 20100827).
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

#line 1 "snrm2.f"
/* > \brief \b SNRM2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SNRM2(N,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SNRM2 returns the euclidean norm of a vector via the function */
/* > name, so that */
/* > */
/* >    SNRM2 := sqrt( x'*x ). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of SX */
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
/* >  -- This version written on 25-October-1982. */
/* >     Modified on 14-October-1993 to inline the call to SLASSQ. */
/* >     Sven Hammarling, Nag Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal snrm2_(integer *n, doublereal *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ix;
    static doublereal ssq, norm, scale, absxi;


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 102 "snrm2.f"
    /* Parameter adjustments */
#line 102 "snrm2.f"
    --x;
#line 102 "snrm2.f"

#line 102 "snrm2.f"
    /* Function Body */
#line 102 "snrm2.f"
    if (*n < 1 || *incx < 1) {
#line 103 "snrm2.f"
	norm = 0.;
#line 104 "snrm2.f"
    } else if (*n == 1) {
#line 105 "snrm2.f"
	norm = abs(x[1]);
#line 106 "snrm2.f"
    } else {
#line 107 "snrm2.f"
	scale = 0.;
#line 108 "snrm2.f"
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK */
/*        auxiliary routine: */
/*        CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

#line 113 "snrm2.f"
	i__1 = (*n - 1) * *incx + 1;
#line 113 "snrm2.f"
	i__2 = *incx;
#line 113 "snrm2.f"
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
#line 114 "snrm2.f"
	    if (x[ix] != 0.) {
#line 115 "snrm2.f"
		absxi = (d__1 = x[ix], abs(d__1));
#line 116 "snrm2.f"
		if (scale < absxi) {
/* Computing 2nd power */
#line 117 "snrm2.f"
		    d__1 = scale / absxi;
#line 117 "snrm2.f"
		    ssq = ssq * (d__1 * d__1) + 1.;
#line 118 "snrm2.f"
		    scale = absxi;
#line 119 "snrm2.f"
		} else {
/* Computing 2nd power */
#line 120 "snrm2.f"
		    d__1 = absxi / scale;
#line 120 "snrm2.f"
		    ssq += d__1 * d__1;
#line 121 "snrm2.f"
		}
#line 122 "snrm2.f"
	    }
#line 123 "snrm2.f"
/* L10: */
#line 123 "snrm2.f"
	}
#line 124 "snrm2.f"
	norm = scale * sqrt(ssq);
#line 125 "snrm2.f"
    }

#line 127 "snrm2.f"
    ret_val = norm;
#line 128 "snrm2.f"
    return ret_val;

/*     End of SNRM2. */

} /* snrm2_ */

