#line 1 "dlaev2.f"
/* dlaev2.f -- translated by f2c (version 20100827).
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

#line 1 "dlaev2.f"
/* > \brief \b DLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAEV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaev2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaev2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaev2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1 */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix */
/* >    [  A   B  ] */
/* >    [  B   C  ]. */
/* > On return, RT1 is the eigenvalue of larger absolute value, RT2 is the */
/* > eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right */
/* > eigenvector for RT1, giving the decomposition */
/* > */
/* >    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ] */
/* >    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ]. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION */
/* >          The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION */
/* >          The (1,2) element and the conjugate of the (2,1) element of */
/* >          the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >          The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* >          RT1 is DOUBLE PRECISION */
/* >          The eigenvalue of larger absolute value. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* >          RT2 is DOUBLE PRECISION */
/* >          The eigenvalue of smaller absolute value. */
/* > \endverbatim */
/* > */
/* > \param[out] CS1 */
/* > \verbatim */
/* >          CS1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] SN1 */
/* > \verbatim */
/* >          SN1 is DOUBLE PRECISION */
/* >          The vector (CS1, SN1) is a unit right eigenvector for RT1. */
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
/* >  CS1 and SN1 are accurate to a few ulps barring over/underflow. */
/* > */
/* >  Overflow is possible only if RT1 is within a factor of 5 of overflow. */
/* >  Underflow is harmless if the input data is 0 or exceeds */
/* >     underflow_threshold / macheps. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaev2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    static integer sgn1, sgn2;
    static doublereal acmn, acmx;


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

#line 156 "dlaev2.f"
    sm = *a + *c__;
#line 157 "dlaev2.f"
    df = *a - *c__;
#line 158 "dlaev2.f"
    adf = abs(df);
#line 159 "dlaev2.f"
    tb = *b + *b;
#line 160 "dlaev2.f"
    ab = abs(tb);
#line 161 "dlaev2.f"
    if (abs(*a) > abs(*c__)) {
#line 162 "dlaev2.f"
	acmx = *a;
#line 163 "dlaev2.f"
	acmn = *c__;
#line 164 "dlaev2.f"
    } else {
#line 165 "dlaev2.f"
	acmx = *c__;
#line 166 "dlaev2.f"
	acmn = *a;
#line 167 "dlaev2.f"
    }
#line 168 "dlaev2.f"
    if (adf > ab) {
/* Computing 2nd power */
#line 169 "dlaev2.f"
	d__1 = ab / adf;
#line 169 "dlaev2.f"
	rt = adf * sqrt(d__1 * d__1 + 1.);
#line 170 "dlaev2.f"
    } else if (adf < ab) {
/* Computing 2nd power */
#line 171 "dlaev2.f"
	d__1 = adf / ab;
#line 171 "dlaev2.f"
	rt = ab * sqrt(d__1 * d__1 + 1.);
#line 172 "dlaev2.f"
    } else {

/*        Includes case AB=ADF=0 */

#line 176 "dlaev2.f"
	rt = ab * sqrt(2.);
#line 177 "dlaev2.f"
    }
#line 178 "dlaev2.f"
    if (sm < 0.) {
#line 179 "dlaev2.f"
	*rt1 = (sm - rt) * .5;
#line 180 "dlaev2.f"
	sgn1 = -1;

/*        Order of execution important. */
/*        To get fully accurate smaller eigenvalue, */
/*        next line needs to be executed in higher precision. */

#line 186 "dlaev2.f"
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
#line 187 "dlaev2.f"
    } else if (sm > 0.) {
#line 188 "dlaev2.f"
	*rt1 = (sm + rt) * .5;
#line 189 "dlaev2.f"
	sgn1 = 1;

/*        Order of execution important. */
/*        To get fully accurate smaller eigenvalue, */
/*        next line needs to be executed in higher precision. */

#line 195 "dlaev2.f"
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
#line 196 "dlaev2.f"
    } else {

/*        Includes case RT1 = RT2 = 0 */

#line 200 "dlaev2.f"
	*rt1 = rt * .5;
#line 201 "dlaev2.f"
	*rt2 = rt * -.5;
#line 202 "dlaev2.f"
	sgn1 = 1;
#line 203 "dlaev2.f"
    }

/*     Compute the eigenvector */

#line 207 "dlaev2.f"
    if (df >= 0.) {
#line 208 "dlaev2.f"
	cs = df + rt;
#line 209 "dlaev2.f"
	sgn2 = 1;
#line 210 "dlaev2.f"
    } else {
#line 211 "dlaev2.f"
	cs = df - rt;
#line 212 "dlaev2.f"
	sgn2 = -1;
#line 213 "dlaev2.f"
    }
#line 214 "dlaev2.f"
    acs = abs(cs);
#line 215 "dlaev2.f"
    if (acs > ab) {
#line 216 "dlaev2.f"
	ct = -tb / cs;
#line 217 "dlaev2.f"
	*sn1 = 1. / sqrt(ct * ct + 1.);
#line 218 "dlaev2.f"
	*cs1 = ct * *sn1;
#line 219 "dlaev2.f"
    } else {
#line 220 "dlaev2.f"
	if (ab == 0.) {
#line 221 "dlaev2.f"
	    *cs1 = 1.;
#line 222 "dlaev2.f"
	    *sn1 = 0.;
#line 223 "dlaev2.f"
	} else {
#line 224 "dlaev2.f"
	    tn = -cs / tb;
#line 225 "dlaev2.f"
	    *cs1 = 1. / sqrt(tn * tn + 1.);
#line 226 "dlaev2.f"
	    *sn1 = tn * *cs1;
#line 227 "dlaev2.f"
	}
#line 228 "dlaev2.f"
    }
#line 229 "dlaev2.f"
    if (sgn1 == sgn2) {
#line 230 "dlaev2.f"
	tn = *cs1;
#line 231 "dlaev2.f"
	*cs1 = -(*sn1);
#line 232 "dlaev2.f"
	*sn1 = tn;
#line 233 "dlaev2.f"
    }
#line 234 "dlaev2.f"
    return 0;

/*     End of DLAEV2 */

} /* dlaev2_ */

