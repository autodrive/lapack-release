#line 1 "icmax1.f"
/* icmax1.f -- translated by f2c (version 20100827).
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

#line 1 "icmax1.f"
/* > \brief \b ICMAX1 finds the index of the first vector element of maximum absolute value. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ICMAX1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/icmax1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/icmax1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/icmax1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER          FUNCTION ICMAX1( N, CX, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            CX( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ICMAX1 finds the index of the first vector element of maximum absolute value. */
/* > */
/* > Based on ICAMAX from Level 1 BLAS. */
/* > The change is to use the 'genuine' absolute value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of elements in the vector CX. */
/* > \endverbatim */
/* > */
/* > \param[in] CX */
/* > \verbatim */
/* >          CX is COMPLEX array, dimension (N) */
/* >          The vector CX. The ICMAX1 function returns the index of its first */
/* >          element of maximum absolute value. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The spacing between successive values of CX.  INCX >= 1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date February 2014 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Nick Higham for use with CLACON. */

/*  ===================================================================== */
integer icmax1_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, ix;
    static doublereal smax;


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     February 2014 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 107 "icmax1.f"
    /* Parameter adjustments */
#line 107 "icmax1.f"
    --cx;
#line 107 "icmax1.f"

#line 107 "icmax1.f"
    /* Function Body */
#line 107 "icmax1.f"
    ret_val = 0;
#line 108 "icmax1.f"
    if (*n < 1 || *incx <= 0) {
#line 108 "icmax1.f"
	return ret_val;
#line 108 "icmax1.f"
    }
#line 109 "icmax1.f"
    ret_val = 1;
#line 110 "icmax1.f"
    if (*n == 1) {
#line 110 "icmax1.f"
	return ret_val;
#line 110 "icmax1.f"
    }
#line 111 "icmax1.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 115 "icmax1.f"
	smax = z_abs(&cx[1]);
#line 116 "icmax1.f"
	i__1 = *n;
#line 116 "icmax1.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 117 "icmax1.f"
	    if (z_abs(&cx[i__]) > smax) {
#line 118 "icmax1.f"
		ret_val = i__;
#line 119 "icmax1.f"
		smax = z_abs(&cx[i__]);
#line 120 "icmax1.f"
	    }
#line 121 "icmax1.f"
	}
#line 122 "icmax1.f"
    } else {

/*        code for increment not equal to 1 */

#line 126 "icmax1.f"
	ix = 1;
#line 127 "icmax1.f"
	smax = z_abs(&cx[1]);
#line 128 "icmax1.f"
	ix += *incx;
#line 129 "icmax1.f"
	i__1 = *n;
#line 129 "icmax1.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 130 "icmax1.f"
	    if (z_abs(&cx[ix]) > smax) {
#line 131 "icmax1.f"
		ret_val = i__;
#line 132 "icmax1.f"
		smax = z_abs(&cx[ix]);
#line 133 "icmax1.f"
	    }
#line 134 "icmax1.f"
	    ix += *incx;
#line 135 "icmax1.f"
	}
#line 136 "icmax1.f"
    }
#line 137 "icmax1.f"
    return ret_val;

/*     End of ICMAX1 */

} /* icmax1_ */

