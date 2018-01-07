#line 1 "dznrm2.f"
/* dznrm2.f -- translated by f2c (version 20100827).
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

#line 1 "dznrm2.f"
/* > \brief \b DZNRM2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DZNRM2 returns the euclidean norm of a vector via the function */
/* > name, so that */
/* > */
/* >    DZNRM2 := sqrt( x**H*x ) */
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
/* >          X is COMPLEX*16 array, dimension (N) */
/* >         complex vector with N elements */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of X */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  -- This version written on 25-October-1982. */
/* >     Modified on 14-October-1993 to inline the call to ZLASSQ. */
/* >     Sven Hammarling, Nag Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal dznrm2_(integer *n, doublecomplex *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer ix;
    static doublereal ssq, temp, norm, scale;


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
#line 103 "dznrm2.f"
    /* Parameter adjustments */
#line 103 "dznrm2.f"
    --x;
#line 103 "dznrm2.f"

#line 103 "dznrm2.f"
    /* Function Body */
#line 103 "dznrm2.f"
    if (*n < 1 || *incx < 1) {
#line 104 "dznrm2.f"
	norm = 0.;
#line 105 "dznrm2.f"
    } else {
#line 106 "dznrm2.f"
	scale = 0.;
#line 107 "dznrm2.f"
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK */
/*        auxiliary routine: */
/*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ ) */

#line 112 "dznrm2.f"
	i__1 = (*n - 1) * *incx + 1;
#line 112 "dznrm2.f"
	i__2 = *incx;
#line 112 "dznrm2.f"
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
#line 113 "dznrm2.f"
	    i__3 = ix;
#line 113 "dznrm2.f"
	    if (x[i__3].r != 0.) {
#line 114 "dznrm2.f"
		i__3 = ix;
#line 114 "dznrm2.f"
		temp = (d__1 = x[i__3].r, abs(d__1));
#line 115 "dznrm2.f"
		if (scale < temp) {
/* Computing 2nd power */
#line 116 "dznrm2.f"
		    d__1 = scale / temp;
#line 116 "dznrm2.f"
		    ssq = ssq * (d__1 * d__1) + 1.;
#line 117 "dznrm2.f"
		    scale = temp;
#line 118 "dznrm2.f"
		} else {
/* Computing 2nd power */
#line 119 "dznrm2.f"
		    d__1 = temp / scale;
#line 119 "dznrm2.f"
		    ssq += d__1 * d__1;
#line 120 "dznrm2.f"
		}
#line 121 "dznrm2.f"
	    }
#line 122 "dznrm2.f"
	    if (d_imag(&x[ix]) != 0.) {
#line 123 "dznrm2.f"
		temp = (d__1 = d_imag(&x[ix]), abs(d__1));
#line 124 "dznrm2.f"
		if (scale < temp) {
/* Computing 2nd power */
#line 125 "dznrm2.f"
		    d__1 = scale / temp;
#line 125 "dznrm2.f"
		    ssq = ssq * (d__1 * d__1) + 1.;
#line 126 "dznrm2.f"
		    scale = temp;
#line 127 "dznrm2.f"
		} else {
/* Computing 2nd power */
#line 128 "dznrm2.f"
		    d__1 = temp / scale;
#line 128 "dznrm2.f"
		    ssq += d__1 * d__1;
#line 129 "dznrm2.f"
		}
#line 130 "dznrm2.f"
	    }
#line 131 "dznrm2.f"
/* L10: */
#line 131 "dznrm2.f"
	}
#line 132 "dznrm2.f"
	norm = scale * sqrt(ssq);
#line 133 "dznrm2.f"
    }

#line 135 "dznrm2.f"
    ret_val = norm;
#line 136 "dznrm2.f"
    return ret_val;

/*     End of DZNRM2. */

} /* dznrm2_ */

