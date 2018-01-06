#line 1 "dlaed5.f"
/* dlaed5.f -- translated by f2c (version 20100827).
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

#line 1 "dlaed5.f"
/* > \brief \b DLAED5 used by sstedc. Solves the 2-by-2 secular equation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAED5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I */
/*       DOUBLE PRECISION   DLAM, RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( 2 ), DELTA( 2 ), Z( 2 ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine computes the I-th eigenvalue of a symmetric rank-one */
/* > modification of a 2-by-2 diagonal matrix */
/* > */
/* >            diag( D )  +  RHO * Z * transpose(Z) . */
/* > */
/* > The diagonal elements in the array D are assumed to satisfy */
/* > */
/* >            D(i) < D(j)  for  i < j . */
/* > */
/* > We also assume RHO > 0 and that the Euclidean norm of the vector */
/* > Z is one. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] I */
/* > \verbatim */
/* >          I is INTEGER */
/* >         The index of the eigenvalue to be computed.  I = 1 or I = 2. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (2) */
/* >         The original eigenvalues.  We assume D(1) < D(2). */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (2) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is DOUBLE PRECISION array, dimension (2) */
/* >         The vector DELTA contains the information necessary */
/* >         to construct the eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAM */
/* > \verbatim */
/* >          DLAM is DOUBLE PRECISION */
/* >         The computed lambda_I, the I-th updated eigenvalue. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaed5_(integer *i__, doublereal *d__, doublereal *z__, 
	doublereal *delta, doublereal *rho, doublereal *dlam)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal b, c__, w, del, tau, temp;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 139 "dlaed5.f"
    /* Parameter adjustments */
#line 139 "dlaed5.f"
    --delta;
#line 139 "dlaed5.f"
    --z__;
#line 139 "dlaed5.f"
    --d__;
#line 139 "dlaed5.f"

#line 139 "dlaed5.f"
    /* Function Body */
#line 139 "dlaed5.f"
    del = d__[2] - d__[1];
#line 140 "dlaed5.f"
    if (*i__ == 1) {
#line 141 "dlaed5.f"
	w = *rho * 2. * (z__[2] * z__[2] - z__[1] * z__[1]) / del + 1.;
#line 142 "dlaed5.f"
	if (w > 0.) {
#line 143 "dlaed5.f"
	    b = del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
#line 144 "dlaed5.f"
	    c__ = *rho * z__[1] * z__[1] * del;

/*           B > ZERO, always */

#line 148 "dlaed5.f"
	    tau = c__ * 2. / (b + sqrt((d__1 = b * b - c__ * 4., abs(d__1))));
#line 149 "dlaed5.f"
	    *dlam = d__[1] + tau;
#line 150 "dlaed5.f"
	    delta[1] = -z__[1] / tau;
#line 151 "dlaed5.f"
	    delta[2] = z__[2] / (del - tau);
#line 152 "dlaed5.f"
	} else {
#line 153 "dlaed5.f"
	    b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
#line 154 "dlaed5.f"
	    c__ = *rho * z__[2] * z__[2] * del;
#line 155 "dlaed5.f"
	    if (b > 0.) {
#line 156 "dlaed5.f"
		tau = c__ * -2. / (b + sqrt(b * b + c__ * 4.));
#line 157 "dlaed5.f"
	    } else {
#line 158 "dlaed5.f"
		tau = (b - sqrt(b * b + c__ * 4.)) / 2.;
#line 159 "dlaed5.f"
	    }
#line 160 "dlaed5.f"
	    *dlam = d__[2] + tau;
#line 161 "dlaed5.f"
	    delta[1] = -z__[1] / (del + tau);
#line 162 "dlaed5.f"
	    delta[2] = -z__[2] / tau;
#line 163 "dlaed5.f"
	}
#line 164 "dlaed5.f"
	temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
#line 165 "dlaed5.f"
	delta[1] /= temp;
#line 166 "dlaed5.f"
	delta[2] /= temp;
#line 167 "dlaed5.f"
    } else {

/*     Now I=2 */

#line 171 "dlaed5.f"
	b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
#line 172 "dlaed5.f"
	c__ = *rho * z__[2] * z__[2] * del;
#line 173 "dlaed5.f"
	if (b > 0.) {
#line 174 "dlaed5.f"
	    tau = (b + sqrt(b * b + c__ * 4.)) / 2.;
#line 175 "dlaed5.f"
	} else {
#line 176 "dlaed5.f"
	    tau = c__ * 2. / (-b + sqrt(b * b + c__ * 4.));
#line 177 "dlaed5.f"
	}
#line 178 "dlaed5.f"
	*dlam = d__[2] + tau;
#line 179 "dlaed5.f"
	delta[1] = -z__[1] / (del + tau);
#line 180 "dlaed5.f"
	delta[2] = -z__[2] / tau;
#line 181 "dlaed5.f"
	temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
#line 182 "dlaed5.f"
	delta[1] /= temp;
#line 183 "dlaed5.f"
	delta[2] /= temp;
#line 184 "dlaed5.f"
    }
#line 185 "dlaed5.f"
    return 0;

/*     End OF DLAED5 */

} /* dlaed5_ */

