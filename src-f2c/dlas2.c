#line 1 "dlas2.f"
/* dlas2.f -- translated by f2c (version 20100827).
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

#line 1 "dlas2.f"
/* > \brief \b DLAS2 computes singular values of a 2-by-2 triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlas2.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlas2.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlas2.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   F, G, H, SSMAX, SSMIN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAS2  computes the singular values of the 2-by-2 matrix */
/* >    [  F   G  ] */
/* >    [  0   H  ]. */
/* > On return, SSMIN is the smaller singular value and SSMAX is the */
/* > larger singular value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is DOUBLE PRECISION */
/* >          The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is DOUBLE PRECISION */
/* >          The (1,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION */
/* >          The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* >          SSMIN is DOUBLE PRECISION */
/* >          The smaller singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMAX */
/* > \verbatim */
/* >          SSMAX is DOUBLE PRECISION */
/* >          The larger singular value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Barring over/underflow, all output quantities are correct to within */
/* >  a few units in the last place (ulps), even in the absence of a guard */
/* >  digit in addition/subtraction. */
/* > */
/* >  In IEEE arithmetic, the code works correctly if one matrix element is */
/* >  infinite. */
/* > */
/* >  Overflow will not occur unless the largest singular value itself */
/* >  overflows, or is within a few ulps of overflow. (On machines with */
/* >  partial overflow, like the Cray, overflow may occur if the largest */
/* >  singular value is within a factor of 2 of overflow.) */
/* > */
/* >  Underflow is harmless if underflow is gradual. Otherwise, results */
/* >  may correspond to a matrix modified by perturbations of size near */
/* >  the underflow threshold. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlas2_(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, fa, ga, ha, as, at, au, fhmn, fhmx;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ==================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 137 "dlas2.f"
    fa = abs(*f);
#line 138 "dlas2.f"
    ga = abs(*g);
#line 139 "dlas2.f"
    ha = abs(*h__);
#line 140 "dlas2.f"
    fhmn = min(fa,ha);
#line 141 "dlas2.f"
    fhmx = max(fa,ha);
#line 142 "dlas2.f"
    if (fhmn == 0.) {
#line 143 "dlas2.f"
	*ssmin = 0.;
#line 144 "dlas2.f"
	if (fhmx == 0.) {
#line 145 "dlas2.f"
	    *ssmax = ga;
#line 146 "dlas2.f"
	} else {
/* Computing 2nd power */
#line 147 "dlas2.f"
	    d__1 = min(fhmx,ga) / max(fhmx,ga);
#line 147 "dlas2.f"
	    *ssmax = max(fhmx,ga) * sqrt(d__1 * d__1 + 1.);
#line 149 "dlas2.f"
	}
#line 150 "dlas2.f"
    } else {
#line 151 "dlas2.f"
	if (ga < fhmx) {
#line 152 "dlas2.f"
	    as = fhmn / fhmx + 1.;
#line 153 "dlas2.f"
	    at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
#line 154 "dlas2.f"
	    d__1 = ga / fhmx;
#line 154 "dlas2.f"
	    au = d__1 * d__1;
#line 155 "dlas2.f"
	    c__ = 2. / (sqrt(as * as + au) + sqrt(at * at + au));
#line 156 "dlas2.f"
	    *ssmin = fhmn * c__;
#line 157 "dlas2.f"
	    *ssmax = fhmx / c__;
#line 158 "dlas2.f"
	} else {
#line 159 "dlas2.f"
	    au = fhmx / ga;
#line 160 "dlas2.f"
	    if (au == 0.) {

/*              Avoid possible harmful underflow if exponent range */
/*              asymmetric (true SSMIN may not underflow even if */
/*              AU underflows) */

#line 166 "dlas2.f"
		*ssmin = fhmn * fhmx / ga;
#line 167 "dlas2.f"
		*ssmax = ga;
#line 168 "dlas2.f"
	    } else {
#line 169 "dlas2.f"
		as = fhmn / fhmx + 1.;
#line 170 "dlas2.f"
		at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
#line 171 "dlas2.f"
		d__1 = as * au;
/* Computing 2nd power */
#line 171 "dlas2.f"
		d__2 = at * au;
#line 171 "dlas2.f"
		c__ = 1. / (sqrt(d__1 * d__1 + 1.) + sqrt(d__2 * d__2 + 1.));
#line 173 "dlas2.f"
		*ssmin = fhmn * c__ * au;
#line 174 "dlas2.f"
		*ssmin += *ssmin;
#line 175 "dlas2.f"
		*ssmax = ga / (c__ + c__);
#line 176 "dlas2.f"
	    }
#line 177 "dlas2.f"
	}
#line 178 "dlas2.f"
    }
#line 179 "dlas2.f"
    return 0;

/*     End of DLAS2 */

} /* dlas2_ */

