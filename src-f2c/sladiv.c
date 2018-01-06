#line 1 "sladiv.f"
/* sladiv.f -- translated by f2c (version 20100827).
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

#line 1 "sladiv.f"
/* > \brief \b SLADIV performs complex division in real arithmetic, avoiding unnecessary overflow. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLADIV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sladiv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sladiv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sladiv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLADIV( A, B, C, D, P, Q ) */

/*       .. Scalar Arguments .. */
/*       REAL               A, B, C, D, P, Q */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLADIV performs complex division in  real arithmetic */
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
/* >          A is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL */
/* >          The scalars a, b, c, and d in the above expression. */
/* > \endverbatim */
/* > */
/* > \param[out] P */
/* > \verbatim */
/* >          P is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is REAL */
/* >          The scalars p and q in the above expression. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2013 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int sladiv_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal s, aa, ab, bb, cc, cd, dd, be, un, ov, eps;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int sladiv1_(doublereal *, doublereal *, 
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

#line 128 "sladiv.f"
    aa = *a;
#line 129 "sladiv.f"
    bb = *b;
#line 130 "sladiv.f"
    cc = *c__;
#line 131 "sladiv.f"
    dd = *d__;
/* Computing MAX */
#line 132 "sladiv.f"
    d__1 = abs(*a), d__2 = abs(*b);
#line 132 "sladiv.f"
    ab = max(d__1,d__2);
/* Computing MAX */
#line 133 "sladiv.f"
    d__1 = abs(*c__), d__2 = abs(*d__);
#line 133 "sladiv.f"
    cd = max(d__1,d__2);
#line 134 "sladiv.f"
    s = 1.;
#line 136 "sladiv.f"
    ov = slamch_("Overflow threshold", (ftnlen)18);
#line 137 "sladiv.f"
    un = slamch_("Safe minimum", (ftnlen)12);
#line 138 "sladiv.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 139 "sladiv.f"
    be = 2. / (eps * eps);
#line 141 "sladiv.f"
    if (ab >= ov * .5) {
#line 142 "sladiv.f"
	aa *= .5;
#line 143 "sladiv.f"
	bb *= .5;
#line 144 "sladiv.f"
	s *= 2.;
#line 145 "sladiv.f"
    }
#line 146 "sladiv.f"
    if (cd >= ov * .5) {
#line 147 "sladiv.f"
	cc *= .5;
#line 148 "sladiv.f"
	dd *= .5;
#line 149 "sladiv.f"
	s *= .5;
#line 150 "sladiv.f"
    }
#line 151 "sladiv.f"
    if (ab <= un * 2. / eps) {
#line 152 "sladiv.f"
	aa *= be;
#line 153 "sladiv.f"
	bb *= be;
#line 154 "sladiv.f"
	s /= be;
#line 155 "sladiv.f"
    }
#line 156 "sladiv.f"
    if (cd <= un * 2. / eps) {
#line 157 "sladiv.f"
	cc *= be;
#line 158 "sladiv.f"
	dd *= be;
#line 159 "sladiv.f"
	s *= be;
#line 160 "sladiv.f"
    }
#line 161 "sladiv.f"
    if (abs(*d__) <= abs(*c__)) {
#line 162 "sladiv.f"
	sladiv1_(&aa, &bb, &cc, &dd, p, q);
#line 163 "sladiv.f"
    } else {
#line 164 "sladiv.f"
	sladiv1_(&bb, &aa, &dd, &cc, p, q);
#line 165 "sladiv.f"
	*q = -(*q);
#line 166 "sladiv.f"
    }
#line 167 "sladiv.f"
    *p *= s;
#line 168 "sladiv.f"
    *q *= s;

#line 170 "sladiv.f"
    return 0;

/*     End of SLADIV */

} /* sladiv_ */

/* > \ingroup realOTHERauxiliary */
/* Subroutine */ int sladiv1_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q)
{
    static doublereal r__, t;
    extern doublereal sladiv2_(doublereal *, doublereal *, doublereal *, 
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

#line 205 "sladiv.f"
    r__ = *d__ / *c__;
#line 206 "sladiv.f"
    t = 1. / (*c__ + *d__ * r__);
#line 207 "sladiv.f"
    *p = sladiv2_(a, b, c__, d__, &r__, &t);
#line 208 "sladiv.f"
    *a = -(*a);
#line 209 "sladiv.f"
    *q = sladiv2_(b, a, c__, d__, &r__, &t);

#line 211 "sladiv.f"
    return 0;

/*     End of SLADIV1 */

} /* sladiv1_ */

/* > \ingroup realOTHERauxiliary */
doublereal sladiv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal 
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

#line 241 "sladiv.f"
    if (*r__ != 0.) {
#line 242 "sladiv.f"
	br = *b * *r__;
#line 243 "sladiv.f"
	if (br != 0.) {
#line 244 "sladiv.f"
	    ret_val = (*a + br) * *t;
#line 245 "sladiv.f"
	} else {
#line 246 "sladiv.f"
	    ret_val = *a * *t + *b * *t * *r__;
#line 247 "sladiv.f"
	}
#line 248 "sladiv.f"
    } else {
#line 249 "sladiv.f"
	ret_val = (*a + *d__ * (*b / *c__)) * *t;
#line 250 "sladiv.f"
    }

#line 252 "sladiv.f"
    return ret_val;

/*     End of SLADIV */

} /* sladiv2_ */

