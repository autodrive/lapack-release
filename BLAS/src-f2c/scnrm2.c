#line 1 "scnrm2.f"
/* scnrm2.f -- translated by f2c (version 20100827).
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

#line 1 "scnrm2.f"
/* > \brief \b SCNRM2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SCNRM2(N,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SCNRM2 returns the euclidean norm of a vector via the function */
/* > name, so that */
/* > */
/* >    SCNRM2 := sqrt( x**H*x ) */
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
/* >  -- This version written on 25-October-1982. */
/* >     Modified on 14-October-1993 to inline the call to CLASSQ. */
/* >     Sven Hammarling, Nag Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal scnrm2_(integer *n, doublecomplex *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer ix;
    static doublereal ssq, temp, norm, scale;


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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
#line 82 "scnrm2.f"
    /* Parameter adjustments */
#line 82 "scnrm2.f"
    --x;
#line 82 "scnrm2.f"

#line 82 "scnrm2.f"
    /* Function Body */
#line 82 "scnrm2.f"
    if (*n < 1 || *incx < 1) {
#line 83 "scnrm2.f"
	norm = 0.;
#line 84 "scnrm2.f"
    } else {
#line 85 "scnrm2.f"
	scale = 0.;
#line 86 "scnrm2.f"
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK */
/*        auxiliary routine: */
/*        CALL CLASSQ( N, X, INCX, SCALE, SSQ ) */

#line 91 "scnrm2.f"
	i__1 = (*n - 1) * *incx + 1;
#line 91 "scnrm2.f"
	i__2 = *incx;
#line 91 "scnrm2.f"
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
#line 92 "scnrm2.f"
	    i__3 = ix;
#line 92 "scnrm2.f"
	    if (x[i__3].r != 0.) {
#line 93 "scnrm2.f"
		i__3 = ix;
#line 93 "scnrm2.f"
		temp = (d__1 = x[i__3].r, abs(d__1));
#line 94 "scnrm2.f"
		if (scale < temp) {
/* Computing 2nd power */
#line 95 "scnrm2.f"
		    d__1 = scale / temp;
#line 95 "scnrm2.f"
		    ssq = ssq * (d__1 * d__1) + 1.;
#line 96 "scnrm2.f"
		    scale = temp;
#line 97 "scnrm2.f"
		} else {
/* Computing 2nd power */
#line 98 "scnrm2.f"
		    d__1 = temp / scale;
#line 98 "scnrm2.f"
		    ssq += d__1 * d__1;
#line 99 "scnrm2.f"
		}
#line 100 "scnrm2.f"
	    }
#line 101 "scnrm2.f"
	    if (d_imag(&x[ix]) != 0.) {
#line 102 "scnrm2.f"
		temp = (d__1 = d_imag(&x[ix]), abs(d__1));
#line 103 "scnrm2.f"
		if (scale < temp) {
/* Computing 2nd power */
#line 104 "scnrm2.f"
		    d__1 = scale / temp;
#line 104 "scnrm2.f"
		    ssq = ssq * (d__1 * d__1) + 1.;
#line 105 "scnrm2.f"
		    scale = temp;
#line 106 "scnrm2.f"
		} else {
/* Computing 2nd power */
#line 107 "scnrm2.f"
		    d__1 = temp / scale;
#line 107 "scnrm2.f"
		    ssq += d__1 * d__1;
#line 108 "scnrm2.f"
		}
#line 109 "scnrm2.f"
	    }
#line 110 "scnrm2.f"
/* L10: */
#line 110 "scnrm2.f"
	}
#line 111 "scnrm2.f"
	norm = scale * sqrt(ssq);
#line 112 "scnrm2.f"
    }

#line 114 "scnrm2.f"
    ret_val = norm;
#line 115 "scnrm2.f"
    return ret_val;

/*     End of SCNRM2. */

} /* scnrm2_ */

