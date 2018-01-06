#line 1 "cla_wwaddw.f"
/* cla_wwaddw.f -- translated by f2c (version 20100827).
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

#line 1 "cla_wwaddw.f"
/* > \brief \b CLA_WWADDW adds a vector into a doubled-single vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_WWADDW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_wwa
ddw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_wwa
ddw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_wwa
ddw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLA_WWADDW( N, X, Y, W ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            X( * ), Y( * ), W( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLA_WWADDW adds a vector W into a doubled-single vector (X, Y). */
/* > */
/* >    This works for all extant IBM's hex and binary floating point */
/* >    arithmetics, but not for decimal. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >            The length of vectors X, Y, and W. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (N) */
/* >            The first part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension (N) */
/* >            The second part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (N) */
/* >            The vector to be added. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cla_wwaddw__(integer *n, doublecomplex *x, doublecomplex 
	*y, doublecomplex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__;
    static doublecomplex s;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 104 "cla_wwaddw.f"
    /* Parameter adjustments */
#line 104 "cla_wwaddw.f"
    --w;
#line 104 "cla_wwaddw.f"
    --y;
#line 104 "cla_wwaddw.f"
    --x;
#line 104 "cla_wwaddw.f"

#line 104 "cla_wwaddw.f"
    /* Function Body */
#line 104 "cla_wwaddw.f"
    i__1 = *n;
#line 104 "cla_wwaddw.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 105 "cla_wwaddw.f"
	i__2 = i__;
#line 105 "cla_wwaddw.f"
	i__3 = i__;
#line 105 "cla_wwaddw.f"
	z__1.r = x[i__2].r + w[i__3].r, z__1.i = x[i__2].i + w[i__3].i;
#line 105 "cla_wwaddw.f"
	s.r = z__1.r, s.i = z__1.i;
#line 106 "cla_wwaddw.f"
	z__2.r = s.r + s.r, z__2.i = s.i + s.i;
#line 106 "cla_wwaddw.f"
	z__1.r = z__2.r - s.r, z__1.i = z__2.i - s.i;
#line 106 "cla_wwaddw.f"
	s.r = z__1.r, s.i = z__1.i;
#line 107 "cla_wwaddw.f"
	i__2 = i__;
#line 107 "cla_wwaddw.f"
	i__3 = i__;
#line 107 "cla_wwaddw.f"
	z__3.r = x[i__3].r - s.r, z__3.i = x[i__3].i - s.i;
#line 107 "cla_wwaddw.f"
	i__4 = i__;
#line 107 "cla_wwaddw.f"
	z__2.r = z__3.r + w[i__4].r, z__2.i = z__3.i + w[i__4].i;
#line 107 "cla_wwaddw.f"
	i__5 = i__;
#line 107 "cla_wwaddw.f"
	z__1.r = z__2.r + y[i__5].r, z__1.i = z__2.i + y[i__5].i;
#line 107 "cla_wwaddw.f"
	y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 108 "cla_wwaddw.f"
	i__2 = i__;
#line 108 "cla_wwaddw.f"
	x[i__2].r = s.r, x[i__2].i = s.i;
#line 109 "cla_wwaddw.f"
/* L10: */
#line 109 "cla_wwaddw.f"
    }
#line 110 "cla_wwaddw.f"
    return 0;
} /* cla_wwaddw__ */

