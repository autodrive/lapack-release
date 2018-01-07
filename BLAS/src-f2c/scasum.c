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

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* >          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
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


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 96 "scasum.f"
    /* Parameter adjustments */
#line 96 "scasum.f"
    --cx;
#line 96 "scasum.f"

#line 96 "scasum.f"
    /* Function Body */
#line 96 "scasum.f"
    ret_val = 0.;
#line 97 "scasum.f"
    stemp = 0.;
#line 98 "scasum.f"
    if (*n <= 0 || *incx <= 0) {
#line 98 "scasum.f"
	return ret_val;
#line 98 "scasum.f"
    }
#line 99 "scasum.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 103 "scasum.f"
	i__1 = *n;
#line 103 "scasum.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 104 "scasum.f"
	    i__2 = i__;
#line 104 "scasum.f"
	    stemp = stemp + (d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&
		    cx[i__]), abs(d__2));
#line 105 "scasum.f"
	}
#line 106 "scasum.f"
    } else {

/*        code for increment not equal to 1 */

#line 110 "scasum.f"
	nincx = *n * *incx;
#line 111 "scasum.f"
	i__1 = nincx;
#line 111 "scasum.f"
	i__2 = *incx;
#line 111 "scasum.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 112 "scasum.f"
	    i__3 = i__;
#line 112 "scasum.f"
	    stemp = stemp + (d__1 = cx[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		    cx[i__]), abs(d__2));
#line 113 "scasum.f"
	}
#line 114 "scasum.f"
    }
#line 115 "scasum.f"
    ret_val = stemp;
#line 116 "scasum.f"
    return ret_val;
} /* scasum_ */

