#line 1 "zlacgv.f"
/* zlacgv.f -- translated by f2c (version 20100827).
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

#line 1 "zlacgv.f"
/* > \brief \b ZLACGV conjugates a complex vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLACGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacgv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacgv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacgv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLACGV( N, X, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACGV conjugates a complex vector of length N. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The length of the vector X.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension */
/* >                         (1+(N-1)*abs(INCX)) */
/* >          On entry, the vector of length N to be conjugated. */
/* >          On exit, X is overwritten with conjg(X). */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The spacing between successive elements of X. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlacgv_(integer *n, doublecomplex *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ioff;


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
/*     .. Executable Statements .. */

#line 99 "zlacgv.f"
    /* Parameter adjustments */
#line 99 "zlacgv.f"
    --x;
#line 99 "zlacgv.f"

#line 99 "zlacgv.f"
    /* Function Body */
#line 99 "zlacgv.f"
    if (*incx == 1) {
#line 100 "zlacgv.f"
	i__1 = *n;
#line 100 "zlacgv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 101 "zlacgv.f"
	    i__2 = i__;
#line 101 "zlacgv.f"
	    d_cnjg(&z__1, &x[i__]);
#line 101 "zlacgv.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 102 "zlacgv.f"
/* L10: */
#line 102 "zlacgv.f"
	}
#line 103 "zlacgv.f"
    } else {
#line 104 "zlacgv.f"
	ioff = 1;
#line 105 "zlacgv.f"
	if (*incx < 0) {
#line 105 "zlacgv.f"
	    ioff = 1 - (*n - 1) * *incx;
#line 105 "zlacgv.f"
	}
#line 107 "zlacgv.f"
	i__1 = *n;
#line 107 "zlacgv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 108 "zlacgv.f"
	    i__2 = ioff;
#line 108 "zlacgv.f"
	    d_cnjg(&z__1, &x[ioff]);
#line 108 "zlacgv.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 109 "zlacgv.f"
	    ioff += *incx;
#line 110 "zlacgv.f"
/* L20: */
#line 110 "zlacgv.f"
	}
#line 111 "zlacgv.f"
    }
#line 112 "zlacgv.f"
    return 0;

/*     End of ZLACGV */

} /* zlacgv_ */

