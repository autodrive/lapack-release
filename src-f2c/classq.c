#line 1 "classq.f"
/* classq.f -- translated by f2c (version 20100827).
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

#line 1 "classq.f"
/* > \brief \b CLASSQ updates a sum of squares represented in scaled form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASSQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/classq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/classq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/classq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       REAL               SCALE, SUMSQ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASSQ returns the values scl and ssq such that */
/* > */
/* >    ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, */
/* > */
/* > where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is */
/* > assumed to be at least unity and the value of ssq will then satisfy */
/* > */
/* >    1.0 .le. ssq .le. ( sumsq + 2*n ). */
/* > */
/* > scale is assumed to be non-negative and scl returns the value */
/* > */
/* >    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ), */
/* >           i */
/* > */
/* > scale and sumsq must be supplied in SCALE and SUMSQ respectively. */
/* > SCALE and SUMSQ are overwritten by scl and ssq respectively. */
/* > */
/* > The routine makes only one pass through the vector X. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of elements to be used from the vector X. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (N) */
/* >          The vector x as described above. */
/* >             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of the vector X. */
/* >          INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          On entry, the value  scale  in the equation above. */
/* >          On exit, SCALE is overwritten with the value  scl . */
/* > \endverbatim */
/* > */
/* > \param[in,out] SUMSQ */
/* > \verbatim */
/* >          SUMSQ is REAL */
/* >          On entry, the value  sumsq  in the equation above. */
/* >          On exit, SUMSQ is overwritten with the value  ssq . */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int classq_(integer *n, doublecomplex *x, integer *incx, 
	doublereal *scale, doublereal *sumsq)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer ix;
    static doublereal temp1;
    extern logical sisnan_(doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 141 "classq.f"
    /* Parameter adjustments */
#line 141 "classq.f"
    --x;
#line 141 "classq.f"

#line 141 "classq.f"
    /* Function Body */
#line 141 "classq.f"
    if (*n > 0) {
#line 142 "classq.f"
	i__1 = (*n - 1) * *incx + 1;
#line 142 "classq.f"
	i__2 = *incx;
#line 142 "classq.f"
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
#line 143 "classq.f"
	    i__3 = ix;
#line 143 "classq.f"
	    temp1 = (d__1 = x[i__3].r, abs(d__1));
#line 144 "classq.f"
	    if (temp1 > 0. || sisnan_(&temp1)) {
#line 145 "classq.f"
		if (*scale < temp1) {
/* Computing 2nd power */
#line 146 "classq.f"
		    d__1 = *scale / temp1;
#line 146 "classq.f"
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
#line 147 "classq.f"
		    *scale = temp1;
#line 148 "classq.f"
		} else {
/* Computing 2nd power */
#line 149 "classq.f"
		    d__1 = temp1 / *scale;
#line 149 "classq.f"
		    *sumsq += d__1 * d__1;
#line 150 "classq.f"
		}
#line 151 "classq.f"
	    }
#line 152 "classq.f"
	    temp1 = (d__1 = d_imag(&x[ix]), abs(d__1));
#line 153 "classq.f"
	    if (temp1 > 0. || sisnan_(&temp1)) {
#line 154 "classq.f"
		if (*scale < temp1 || sisnan_(&temp1)) {
/* Computing 2nd power */
#line 155 "classq.f"
		    d__1 = *scale / temp1;
#line 155 "classq.f"
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
#line 156 "classq.f"
		    *scale = temp1;
#line 157 "classq.f"
		} else {
/* Computing 2nd power */
#line 158 "classq.f"
		    d__1 = temp1 / *scale;
#line 158 "classq.f"
		    *sumsq += d__1 * d__1;
#line 159 "classq.f"
		}
#line 160 "classq.f"
	    }
#line 161 "classq.f"
/* L10: */
#line 161 "classq.f"
	}
#line 162 "classq.f"
    }

#line 164 "classq.f"
    return 0;

/*     End of CLASSQ */

} /* classq_ */

