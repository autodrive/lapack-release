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
/* > \brief \b IZMAX1 finds the index of the first vector element of maximum absolute value. */

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

/*       INTEGER          FUNCTION IZMAX1( N, ZX, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         ZX( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > IZMAX1 finds the index of the first vector element of maximum absolute value. */
/* > */
/* > Based on IZAMAX from Level 1 BLAS. */
/* > The change is to use the 'genuine' absolute value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of elements in the vector ZX. */
/* > \endverbatim */
/* > */
/* > \param[in] ZX */
/* > \verbatim */
/* >          ZX is COMPLEX*16 array, dimension (N) */
/* >          The vector ZX. The IZMAX1 function returns the index of its first */
/* >          element of maximum absolute value. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The spacing between successive values of ZX.  INCX >= 1. */
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
/* > Nick Higham for use with ZLACON. */

/*  ===================================================================== */
integer izmax1_(integer *n, doublecomplex *zx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, ix;
    static doublereal dmax__;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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

#line 107 "izmax1.f"
    /* Parameter adjustments */
#line 107 "izmax1.f"
    --zx;
#line 107 "izmax1.f"

#line 107 "izmax1.f"
    /* Function Body */
#line 107 "izmax1.f"
    ret_val = 0;
#line 108 "izmax1.f"
    if (*n < 1 || *incx <= 0) {
#line 108 "izmax1.f"
	return ret_val;
#line 108 "izmax1.f"
    }
#line 109 "izmax1.f"
    ret_val = 1;
#line 110 "izmax1.f"
    if (*n == 1) {
#line 110 "izmax1.f"
	return ret_val;
#line 110 "izmax1.f"
    }
#line 111 "izmax1.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 115 "izmax1.f"
	dmax__ = z_abs(&zx[1]);
#line 116 "izmax1.f"
	i__1 = *n;
#line 116 "izmax1.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 117 "izmax1.f"
	    if (z_abs(&zx[i__]) > dmax__) {
#line 118 "izmax1.f"
		ret_val = i__;
#line 119 "izmax1.f"
		dmax__ = z_abs(&zx[i__]);
#line 120 "izmax1.f"
	    }
#line 121 "izmax1.f"
	}
#line 122 "izmax1.f"
    } else {

/*        code for increment not equal to 1 */

#line 126 "izmax1.f"
	ix = 1;
#line 127 "izmax1.f"
	dmax__ = z_abs(&zx[1]);
#line 128 "izmax1.f"
	ix += *incx;
#line 129 "izmax1.f"
	i__1 = *n;
#line 129 "izmax1.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 130 "izmax1.f"
	    if (z_abs(&zx[ix]) > dmax__) {
#line 131 "izmax1.f"
		ret_val = i__;
#line 132 "izmax1.f"
		dmax__ = z_abs(&zx[ix]);
#line 133 "izmax1.f"
	    }
#line 134 "izmax1.f"
	    ix += *incx;
#line 135 "izmax1.f"
	}
#line 136 "izmax1.f"
    }
#line 137 "izmax1.f"
    return ret_val;

/*     End of IZMAX1 */

} /* izmax1_ */

