#line 1 "scsum1.f"
/* scsum1.f -- translated by f2c (version 20100827).
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

#line 1 "scsum1.f"
/* > \brief \b SCSUM1 forms the 1-norm of the complex vector using the true absolute value. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SCSUM1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/scsum1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/scsum1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/scsum1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SCSUM1( N, CX, INCX ) */

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
/* > SCSUM1 takes the sum of the absolute values of a complex */
/* > vector and returns a single precision result. */
/* > */
/* > Based on SCASUM from the Level 1 BLAS. */
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
/* >          The vector whose elements will be summed. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The spacing between successive values of CX.  INCX > 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Nick Higham for use with CLACON. */

/*  ===================================================================== */
doublereal scsum1_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, nincx;
    static doublereal stemp;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
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
/*     .. Executable Statements .. */

#line 107 "scsum1.f"
    /* Parameter adjustments */
#line 107 "scsum1.f"
    --cx;
#line 107 "scsum1.f"

#line 107 "scsum1.f"
    /* Function Body */
#line 107 "scsum1.f"
    ret_val = 0.;
#line 108 "scsum1.f"
    stemp = 0.;
#line 109 "scsum1.f"
    if (*n <= 0) {
#line 109 "scsum1.f"
	return ret_val;
#line 109 "scsum1.f"
    }
#line 111 "scsum1.f"
    if (*incx == 1) {
#line 111 "scsum1.f"
	goto L20;
#line 111 "scsum1.f"
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

#line 116 "scsum1.f"
    nincx = *n * *incx;
#line 117 "scsum1.f"
    i__1 = nincx;
#line 117 "scsum1.f"
    i__2 = *incx;
#line 117 "scsum1.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*        NEXT LINE MODIFIED. */

#line 121 "scsum1.f"
	stemp += z_abs(&cx[i__]);
#line 122 "scsum1.f"
/* L10: */
#line 122 "scsum1.f"
    }
#line 123 "scsum1.f"
    ret_val = stemp;
#line 124 "scsum1.f"
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

#line 128 "scsum1.f"
L20:
#line 129 "scsum1.f"
    i__2 = *n;
#line 129 "scsum1.f"
    for (i__ = 1; i__ <= i__2; ++i__) {

/*        NEXT LINE MODIFIED. */

#line 133 "scsum1.f"
	stemp += z_abs(&cx[i__]);
#line 134 "scsum1.f"
/* L30: */
#line 134 "scsum1.f"
    }
#line 135 "scsum1.f"
    ret_val = stemp;
#line 136 "scsum1.f"
    return ret_val;

/*     End of SCSUM1 */

} /* scsum1_ */

