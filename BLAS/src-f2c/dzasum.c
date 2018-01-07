#line 1 "dzasum.f"
/* dzasum.f -- translated by f2c (version 20100827).
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

#line 1 "dzasum.f"
/* > \brief \b DZASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 ZX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DZASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and */
/* >    returns a single precision result. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal dzasum_(integer *n, doublecomplex *zx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i__, nincx;
    static doublereal stemp;
    extern doublereal dcabs1_(doublecomplex *);


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
/*     .. External Functions .. */
/*     .. */
#line 77 "dzasum.f"
    /* Parameter adjustments */
#line 77 "dzasum.f"
    --zx;
#line 77 "dzasum.f"

#line 77 "dzasum.f"
    /* Function Body */
#line 77 "dzasum.f"
    ret_val = 0.;
#line 78 "dzasum.f"
    stemp = 0.;
#line 79 "dzasum.f"
    if (*n <= 0 || *incx <= 0) {
#line 79 "dzasum.f"
	return ret_val;
#line 79 "dzasum.f"
    }
#line 80 "dzasum.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 84 "dzasum.f"
	i__1 = *n;
#line 84 "dzasum.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 85 "dzasum.f"
	    stemp += dcabs1_(&zx[i__]);
#line 86 "dzasum.f"
	}
#line 87 "dzasum.f"
    } else {

/*        code for increment not equal to 1 */

#line 91 "dzasum.f"
	nincx = *n * *incx;
#line 92 "dzasum.f"
	i__1 = nincx;
#line 92 "dzasum.f"
	i__2 = *incx;
#line 92 "dzasum.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 93 "dzasum.f"
	    stemp += dcabs1_(&zx[i__]);
#line 94 "dzasum.f"
	}
#line 95 "dzasum.f"
    }
#line 96 "dzasum.f"
    ret_val = stemp;
#line 97 "dzasum.f"
    return ret_val;
} /* dzasum_ */

