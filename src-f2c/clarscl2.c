#line 1 "clarscl2.f"
/* clarscl2.f -- translated by f2c (version 20100827).
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

#line 1 "clarscl2.f"
/* > \brief \b CLARSCL2 performs reciprocal diagonal scaling on a vector. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARSCL2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarscl
2.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarscl
2.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarscl
2.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARSCL2 ( M, N, D, X, LDX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            M, N, LDX */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            X( LDX, * ) */
/*       REAL               D( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARSCL2 performs a reciprocal diagonal scaling on an vector: */
/* >   x <-- inv(D) * x */
/* > where the REAL diagonal matrix D is stored as a vector. */
/* > */
/* > Eventually to be replaced by BLAS_cge_diag_scale in the new BLAS */
/* > standard. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >     The number of rows of D and X. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >     The number of columns of X. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, length M */
/* >     Diagonal matrix D, stored as a vector of length M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (LDX,N) */
/* >     On entry, the vector X to be scaled by D. */
/* >     On exit, the scaled vector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >     The leading dimension of the vector X. LDX >= M. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int clarscl2_(integer *m, integer *n, doublereal *d__, 
	doublecomplex *x, integer *ldx)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 114 "clarscl2.f"
    /* Parameter adjustments */
#line 114 "clarscl2.f"
    --d__;
#line 114 "clarscl2.f"
    x_dim1 = *ldx;
#line 114 "clarscl2.f"
    x_offset = 1 + x_dim1;
#line 114 "clarscl2.f"
    x -= x_offset;
#line 114 "clarscl2.f"

#line 114 "clarscl2.f"
    /* Function Body */
#line 114 "clarscl2.f"
    i__1 = *n;
#line 114 "clarscl2.f"
    for (j = 1; j <= i__1; ++j) {
#line 115 "clarscl2.f"
	i__2 = *m;
#line 115 "clarscl2.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 116 "clarscl2.f"
	    i__3 = i__ + j * x_dim1;
#line 116 "clarscl2.f"
	    i__4 = i__ + j * x_dim1;
#line 116 "clarscl2.f"
	    i__5 = i__;
#line 116 "clarscl2.f"
	    z__1.r = x[i__4].r / d__[i__5], z__1.i = x[i__4].i / d__[i__5];
#line 116 "clarscl2.f"
	    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 117 "clarscl2.f"
	}
#line 118 "clarscl2.f"
    }
#line 120 "clarscl2.f"
    return 0;
} /* clarscl2_ */

