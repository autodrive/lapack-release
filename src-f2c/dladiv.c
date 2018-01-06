#line 1 "dladiv.f"
/* dladiv.f -- translated by f2c (version 20100827).
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

#line 1 "dladiv.f"
/* > \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLADIV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLADIV( A, B, C, D, P, Q ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   A, B, C, D, P, Q */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLADIV performs complex division in  real arithmetic */
/* > */
/* >                       a + i*b */
/* >            p + i*q = --------- */
/* >                       c + i*d */
/* > */
/* > The algorithm is due to Michael Baudin and Robert L. Smith */
/* > and can be found in the paper */
/* > "A Robust Complex Division in Scilab" */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION */
/* >          The scalars a, b, c, and d in the above expression. */
/* > \endverbatim */
/* > */
/* > \param[out] P */
/* > \verbatim */
/* >          P is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION */
/* >          The scalars p and q in the above expression. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2013 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dladiv_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal s, aa, ab, bb, cc, cd, dd, be, un, ov, eps;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dladiv1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2013 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 128 "dladiv.f"
    aa = *a;
#line 129 "dladiv.f"
    bb = *b;
#line 130 "dladiv.f"
    cc = *c__;
#line 131 "dladiv.f"
    dd = *d__;
/* Computing MAX */
#line 132 "dladiv.f"
    d__1 = abs(*a), d__2 = abs(*b);
#line 132 "dladiv.f"
    ab = max(d__1,d__2);
/* Computing MAX */
#line 133 "dladiv.f"
    d__1 = abs(*c__), d__2 = abs(*d__);
#line 133 "dladiv.f"
    cd = max(d__1,d__2);
#line 134 "dladiv.f"
    s = 1.;
#line 136 "dladiv.f"
    ov = dlamch_("Overflow threshold", (ftnlen)18);
#line 137 "dladiv.f"
    un = dlamch_("Safe minimum", (ftnlen)12);
#line 138 "dladiv.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 139 "dladiv.f"
    be = 2. / (eps * eps);
#line 141 "dladiv.f"
    if (ab >= ov * .5) {
#line 142 "dladiv.f"
	aa *= .5;
#line 143 "dladiv.f"
	bb *= .5;
#line 144 "dladiv.f"
	s *= 2.;
#line 145 "dladiv.f"
    }
#line 146 "dladiv.f"
    if (cd >= ov * .5) {
#line 147 "dladiv.f"
	cc *= .5;
#line 148 "dladiv.f"
	dd *= .5;
#line 149 "dladiv.f"
	s *= .5;
#line 150 "dladiv.f"
    }
#line 151 "dladiv.f"
    if (ab <= un * 2. / eps) {
#line 152 "dladiv.f"
	aa *= be;
#line 153 "dladiv.f"
	bb *= be;
#line 154 "dladiv.f"
	s /= be;
#line 155 "dladiv.f"
    }
#line 156 "dladiv.f"
    if (cd <= un * 2. / eps) {
#line 157 "dladiv.f"
	cc *= be;
#line 158 "dladiv.f"
	dd *= be;
#line 159 "dladiv.f"
	s *= be;
#line 160 "dladiv.f"
    }
#line 161 "dladiv.f"
    if (abs(*d__) <= abs(*c__)) {
#line 162 "dladiv.f"
	dladiv1_(&aa, &bb, &cc, &dd, p, q);
#line 163 "dladiv.f"
    } else {
#line 164 "dladiv.f"
	dladiv1_(&bb, &aa, &dd, &cc, p, q);
#line 165 "dladiv.f"
	*q = -(*q);
#line 166 "dladiv.f"
    }
#line 167 "dladiv.f"
    *p *= s;
#line 168 "dladiv.f"
    *q *= s;

#line 170 "dladiv.f"
    return 0;

/*     End of DLADIV */

} /* dladiv_ */

/* > \ingroup doubleOTHERauxiliary */
/* Subroutine */ int dladiv1_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q)
{
    static doublereal r__, t;
    extern doublereal dladiv2_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2013 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 205 "dladiv.f"
    r__ = *d__ / *c__;
#line 206 "dladiv.f"
    t = 1. / (*c__ + *d__ * r__);
#line 207 "dladiv.f"
    *p = dladiv2_(a, b, c__, d__, &r__, &t);
#line 208 "dladiv.f"
    *a = -(*a);
#line 209 "dladiv.f"
    *q = dladiv2_(b, a, c__, d__, &r__, &t);

#line 211 "dladiv.f"
    return 0;

/*     End of DLADIV1 */

} /* dladiv1_ */

/* > \ingroup doubleOTHERauxiliary */
doublereal dladiv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal 
	*d__, doublereal *r__, doublereal *t)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal br;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2013 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 241 "dladiv.f"
    if (*r__ != 0.) {
#line 242 "dladiv.f"
	br = *b * *r__;
#line 243 "dladiv.f"
	if (br != 0.) {
#line 244 "dladiv.f"
	    ret_val = (*a + br) * *t;
#line 245 "dladiv.f"
	} else {
#line 246 "dladiv.f"
	    ret_val = *a * *t + *b * *t * *r__;
#line 247 "dladiv.f"
	}
#line 248 "dladiv.f"
    } else {
#line 249 "dladiv.f"
	ret_val = (*a + *d__ * (*b / *c__)) * *t;
#line 250 "dladiv.f"
    }

#line 252 "dladiv.f"
    return ret_val;

/*     End of DLADIV12 */

} /* dladiv2_ */

