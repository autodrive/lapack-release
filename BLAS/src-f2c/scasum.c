#line 1 "scasum.f"
/* scasum.f -- translated by f2c (version 20100827).
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

#line 1 "scasum.f"
/* > \brief \b SCASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SCASUM(N,CX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX CX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SCASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and */
/* >    returns a single precision result. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup single_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal scasum_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, nincx;
    static doublereal stemp;


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 76 "scasum.f"
    /* Parameter adjustments */
#line 76 "scasum.f"
    --cx;
#line 76 "scasum.f"

#line 76 "scasum.f"
    /* Function Body */
#line 76 "scasum.f"
    ret_val = 0.;
#line 77 "scasum.f"
    stemp = 0.;
#line 78 "scasum.f"
    if (*n <= 0 || *incx <= 0) {
#line 78 "scasum.f"
	return ret_val;
#line 78 "scasum.f"
    }
#line 79 "scasum.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 83 "scasum.f"
	i__1 = *n;
#line 83 "scasum.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 84 "scasum.f"
	    i__2 = i__;
#line 84 "scasum.f"
	    stemp = stemp + (d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&
		    cx[i__]), abs(d__2));
#line 85 "scasum.f"
	}
#line 86 "scasum.f"
    } else {

/*        code for increment not equal to 1 */

#line 90 "scasum.f"
	nincx = *n * *incx;
#line 91 "scasum.f"
	i__1 = nincx;
#line 91 "scasum.f"
	i__2 = *incx;
#line 91 "scasum.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 92 "scasum.f"
	    i__3 = i__;
#line 92 "scasum.f"
	    stemp = stemp + (d__1 = cx[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		    cx[i__]), abs(d__2));
#line 93 "scasum.f"
	}
#line 94 "scasum.f"
    }
#line 95 "scasum.f"
    ret_val = stemp;
#line 96 "scasum.f"
    return ret_val;
} /* scasum_ */

