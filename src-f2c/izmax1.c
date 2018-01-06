#line 1 "izmax1.f"
/* izmax1.f -- translated by f2c (version 20100827).
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

#line 1 "izmax1.f"
/* > \brief \b IZMAX1 finds the index of the vector element whose real part has maximum absolute value. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IZMAX1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/izmax1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/izmax1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/izmax1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER          FUNCTION IZMAX1( N, CX, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         CX( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > IZMAX1 finds the index of the element whose real part has maximum */
/* > absolute value. */
/* > */
/* > Based on IZAMAX from Level 1 BLAS. */
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
/* >          CX is COMPLEX*16 array, dimension (N) */
/* >          The vector whose elements will be summed. */
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

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Nick Higham for use with ZLACON. */

/*  ===================================================================== */
integer izmax1_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, ix;
    static doublereal smax;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */

/*     NEXT LINE IS THE ONLY MODIFICATION. */
/*     .. */
/*     .. Executable Statements .. */

#line 116 "izmax1.f"
    /* Parameter adjustments */
#line 116 "izmax1.f"
    --cx;
#line 116 "izmax1.f"

#line 116 "izmax1.f"
    /* Function Body */
#line 116 "izmax1.f"
    ret_val = 0;
#line 117 "izmax1.f"
    if (*n < 1) {
#line 117 "izmax1.f"
	return ret_val;
#line 117 "izmax1.f"
    }
#line 119 "izmax1.f"
    ret_val = 1;
#line 120 "izmax1.f"
    if (*n == 1) {
#line 120 "izmax1.f"
	return ret_val;
#line 120 "izmax1.f"
    }
#line 122 "izmax1.f"
    if (*incx == 1) {
#line 122 "izmax1.f"
	goto L30;
#line 122 "izmax1.f"
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

#line 127 "izmax1.f"
    ix = 1;
#line 128 "izmax1.f"
    smax = z_abs(&cx[1]);
#line 129 "izmax1.f"
    ix += *incx;
#line 130 "izmax1.f"
    i__1 = *n;
#line 130 "izmax1.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 131 "izmax1.f"
	if (z_abs(&cx[ix]) <= smax) {
#line 131 "izmax1.f"
	    goto L10;
#line 131 "izmax1.f"
	}
#line 133 "izmax1.f"
	ret_val = i__;
#line 134 "izmax1.f"
	smax = z_abs(&cx[ix]);
#line 135 "izmax1.f"
L10:
#line 136 "izmax1.f"
	ix += *incx;
#line 137 "izmax1.f"
/* L20: */
#line 137 "izmax1.f"
    }
#line 138 "izmax1.f"
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

#line 142 "izmax1.f"
L30:
#line 143 "izmax1.f"
    smax = z_abs(&cx[1]);
#line 144 "izmax1.f"
    i__1 = *n;
#line 144 "izmax1.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 145 "izmax1.f"
	if (z_abs(&cx[i__]) <= smax) {
#line 145 "izmax1.f"
	    goto L40;
#line 145 "izmax1.f"
	}
#line 147 "izmax1.f"
	ret_val = i__;
#line 148 "izmax1.f"
	smax = z_abs(&cx[i__]);
#line 149 "izmax1.f"
L40:
#line 149 "izmax1.f"
	;
#line 149 "izmax1.f"
    }
#line 150 "izmax1.f"
    return ret_val;

/*     End of IZMAX1 */

} /* izmax1_ */

