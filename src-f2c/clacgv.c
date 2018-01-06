#line 1 "clacgv.f"
/* clacgv.f -- translated by f2c (version 20100827).
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

#line 1 "clacgv.f"
/* > \brief \b CLACGV conjugates a complex vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLACGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacgv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacgv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacgv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLACGV( N, X, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACGV conjugates a complex vector of length N. */
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
/* >          X is COMPLEX array, dimension */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clacgv_(integer *n, doublecomplex *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ioff;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 99 "clacgv.f"
    /* Parameter adjustments */
#line 99 "clacgv.f"
    --x;
#line 99 "clacgv.f"

#line 99 "clacgv.f"
    /* Function Body */
#line 99 "clacgv.f"
    if (*incx == 1) {
#line 100 "clacgv.f"
	i__1 = *n;
#line 100 "clacgv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 101 "clacgv.f"
	    i__2 = i__;
#line 101 "clacgv.f"
	    d_cnjg(&z__1, &x[i__]);
#line 101 "clacgv.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 102 "clacgv.f"
/* L10: */
#line 102 "clacgv.f"
	}
#line 103 "clacgv.f"
    } else {
#line 104 "clacgv.f"
	ioff = 1;
#line 105 "clacgv.f"
	if (*incx < 0) {
#line 105 "clacgv.f"
	    ioff = 1 - (*n - 1) * *incx;
#line 105 "clacgv.f"
	}
#line 107 "clacgv.f"
	i__1 = *n;
#line 107 "clacgv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 108 "clacgv.f"
	    i__2 = ioff;
#line 108 "clacgv.f"
	    d_cnjg(&z__1, &x[ioff]);
#line 108 "clacgv.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 109 "clacgv.f"
	    ioff += *incx;
#line 110 "clacgv.f"
/* L20: */
#line 110 "clacgv.f"
	}
#line 111 "clacgv.f"
    }
#line 112 "clacgv.f"
    return 0;

/*     End of CLACGV */

} /* clacgv_ */

