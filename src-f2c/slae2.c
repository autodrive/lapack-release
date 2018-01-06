#line 1 "slae2.f"
/* slae2.f -- translated by f2c (version 20100827).
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

#line 1 "slae2.f"
/* > \brief \b SLAE2 computes the eigenvalues of a 2-by-2 symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAE2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slae2.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slae2.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slae2.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAE2( A, B, C, RT1, RT2 ) */

/*       .. Scalar Arguments .. */
/*       REAL               A, B, C, RT1, RT2 */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix */
/* >    [  A   B  ] */
/* >    [  B   C  ]. */
/* > On return, RT1 is the eigenvalue of larger absolute value, and RT2 */
/* > is the eigenvalue of smaller absolute value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL */
/* >          The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL */
/* >          The (1,2) and (2,1) elements of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL */
/* >          The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* >          RT1 is REAL */
/* >          The eigenvalue of larger absolute value. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* >          RT2 is REAL */
/* >          The eigenvalue of smaller absolute value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  RT1 is accurate to a few ulps barring over/underflow. */
/* > */
/* >  RT2 may be inaccurate if there is massive cancellation in the */
/* >  determinant A*C-B*B; higher precision or correctly rounded or */
/* >  correctly truncated arithmetic would be needed to compute RT2 */
/* >  accurately in all cases. */
/* > */
/* >  Overflow is possible only if RT1 is within a factor of 5 of overflow. */
/* >  Underflow is harmless if the input data is 0 or exceeds */
/* >     underflow_threshold / macheps. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slae2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal ab, df, tb, sm, rt, adf, acmn, acmx;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Compute the eigenvalues */

#line 136 "slae2.f"
    sm = *a + *c__;
#line 137 "slae2.f"
    df = *a - *c__;
#line 138 "slae2.f"
    adf = abs(df);
#line 139 "slae2.f"
    tb = *b + *b;
#line 140 "slae2.f"
    ab = abs(tb);
#line 141 "slae2.f"
    if (abs(*a) > abs(*c__)) {
#line 142 "slae2.f"
	acmx = *a;
#line 143 "slae2.f"
	acmn = *c__;
#line 144 "slae2.f"
    } else {
#line 145 "slae2.f"
	acmx = *c__;
#line 146 "slae2.f"
	acmn = *a;
#line 147 "slae2.f"
    }
#line 148 "slae2.f"
    if (adf > ab) {
/* Computing 2nd power */
#line 149 "slae2.f"
	d__1 = ab / adf;
#line 149 "slae2.f"
	rt = adf * sqrt(d__1 * d__1 + 1.);
#line 150 "slae2.f"
    } else if (adf < ab) {
/* Computing 2nd power */
#line 151 "slae2.f"
	d__1 = adf / ab;
#line 151 "slae2.f"
	rt = ab * sqrt(d__1 * d__1 + 1.);
#line 152 "slae2.f"
    } else {

/*        Includes case AB=ADF=0 */

#line 156 "slae2.f"
	rt = ab * sqrt(2.);
#line 157 "slae2.f"
    }
#line 158 "slae2.f"
    if (sm < 0.) {
#line 159 "slae2.f"
	*rt1 = (sm - rt) * .5;

/*        Order of execution important. */
/*        To get fully accurate smaller eigenvalue, */
/*        next line needs to be executed in higher precision. */

#line 165 "slae2.f"
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
#line 166 "slae2.f"
    } else if (sm > 0.) {
#line 167 "slae2.f"
	*rt1 = (sm + rt) * .5;

/*        Order of execution important. */
/*        To get fully accurate smaller eigenvalue, */
/*        next line needs to be executed in higher precision. */

#line 173 "slae2.f"
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
#line 174 "slae2.f"
    } else {

/*        Includes case RT1 = RT2 = 0 */

#line 178 "slae2.f"
	*rt1 = rt * .5;
#line 179 "slae2.f"
	*rt2 = rt * -.5;
#line 180 "slae2.f"
    }
#line 181 "slae2.f"
    return 0;

/*     End of SLAE2 */

} /* slae2_ */

