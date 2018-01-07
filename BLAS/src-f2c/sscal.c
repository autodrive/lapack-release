#line 1 "sscal.f"
/* sscal.f -- translated by f2c (version 20100827).
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

#line 1 "sscal.f"
/* > \brief \b SSCAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSCAL(N,SA,SX,INCX) */

/*       .. Scalar Arguments .. */
/*       REAL SA */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SSCAL scales a vector by a constant. */
/* >    uses unrolled loops for increment equal to 1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] SA */
/* > \verbatim */
/* >          SA is REAL */
/* >           On entry, SA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX */
/* > \verbatim */
/* >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
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
/* Subroutine */ int sscal_(integer *n, doublereal *sa, doublereal *sx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, mp1, nincx;


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
#line 103 "sscal.f"
    /* Parameter adjustments */
#line 103 "sscal.f"
    --sx;
#line 103 "sscal.f"

#line 103 "sscal.f"
    /* Function Body */
#line 103 "sscal.f"
    if (*n <= 0 || *incx <= 0) {
#line 103 "sscal.f"
	return 0;
#line 103 "sscal.f"
    }
#line 104 "sscal.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */


/*        clean-up loop */

#line 111 "sscal.f"
	m = *n % 5;
#line 112 "sscal.f"
	if (m != 0) {
#line 113 "sscal.f"
	    i__1 = m;
#line 113 "sscal.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 114 "sscal.f"
		sx[i__] = *sa * sx[i__];
#line 115 "sscal.f"
	    }
#line 116 "sscal.f"
	    if (*n < 5) {
#line 116 "sscal.f"
		return 0;
#line 116 "sscal.f"
	    }
#line 117 "sscal.f"
	}
#line 118 "sscal.f"
	mp1 = m + 1;
#line 119 "sscal.f"
	i__1 = *n;
#line 119 "sscal.f"
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
#line 120 "sscal.f"
	    sx[i__] = *sa * sx[i__];
#line 121 "sscal.f"
	    sx[i__ + 1] = *sa * sx[i__ + 1];
#line 122 "sscal.f"
	    sx[i__ + 2] = *sa * sx[i__ + 2];
#line 123 "sscal.f"
	    sx[i__ + 3] = *sa * sx[i__ + 3];
#line 124 "sscal.f"
	    sx[i__ + 4] = *sa * sx[i__ + 4];
#line 125 "sscal.f"
	}
#line 126 "sscal.f"
    } else {

/*        code for increment not equal to 1 */

#line 130 "sscal.f"
	nincx = *n * *incx;
#line 131 "sscal.f"
	i__1 = nincx;
#line 131 "sscal.f"
	i__2 = *incx;
#line 131 "sscal.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 132 "sscal.f"
	    sx[i__] = *sa * sx[i__];
#line 133 "sscal.f"
	}
#line 134 "sscal.f"
    }
#line 135 "sscal.f"
    return 0;
} /* sscal_ */

