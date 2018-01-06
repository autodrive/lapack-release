#line 1 "slaed5.f"
/* slaed5.f -- translated by f2c (version 20100827).
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

#line 1 "slaed5.f"
/* > \brief \b SLAED5 used by sstedc. Solves the 2-by-2 secular equation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAED5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAED5( I, D, Z, DELTA, RHO, DLAM ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I */
/*       REAL               DLAM, RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( 2 ), DELTA( 2 ), Z( 2 ) */
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
/* >          D is REAL array, dimension (2) */
/* >         The original eigenvalues.  We assume D(1) < D(2). */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (2) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is REAL array, dimension (2) */
/* >         The vector DELTA contains the information necessary */
/* >         to construct the eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is REAL */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAM */
/* > \verbatim */
/* >          DLAM is REAL */
/* >         The computed lambda_I, the I-th updated eigenvalue. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaed5_(integer *i__, doublereal *d__, doublereal *z__, 
	doublereal *delta, doublereal *rho, doublereal *dlam)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal b, c__, w, del, tau, temp;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 139 "slaed5.f"
    /* Parameter adjustments */
#line 139 "slaed5.f"
    --delta;
#line 139 "slaed5.f"
    --z__;
#line 139 "slaed5.f"
    --d__;
#line 139 "slaed5.f"

#line 139 "slaed5.f"
    /* Function Body */
#line 139 "slaed5.f"
    del = d__[2] - d__[1];
#line 140 "slaed5.f"
    if (*i__ == 1) {
#line 141 "slaed5.f"
	w = *rho * 2. * (z__[2] * z__[2] - z__[1] * z__[1]) / del + 1.;
#line 142 "slaed5.f"
	if (w > 0.) {
#line 143 "slaed5.f"
	    b = del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
#line 144 "slaed5.f"
	    c__ = *rho * z__[1] * z__[1] * del;

/*           B > ZERO, always */

#line 148 "slaed5.f"
	    tau = c__ * 2. / (b + sqrt((d__1 = b * b - c__ * 4., abs(d__1))));
#line 149 "slaed5.f"
	    *dlam = d__[1] + tau;
#line 150 "slaed5.f"
	    delta[1] = -z__[1] / tau;
#line 151 "slaed5.f"
	    delta[2] = z__[2] / (del - tau);
#line 152 "slaed5.f"
	} else {
#line 153 "slaed5.f"
	    b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
#line 154 "slaed5.f"
	    c__ = *rho * z__[2] * z__[2] * del;
#line 155 "slaed5.f"
	    if (b > 0.) {
#line 156 "slaed5.f"
		tau = c__ * -2. / (b + sqrt(b * b + c__ * 4.));
#line 157 "slaed5.f"
	    } else {
#line 158 "slaed5.f"
		tau = (b - sqrt(b * b + c__ * 4.)) / 2.;
#line 159 "slaed5.f"
	    }
#line 160 "slaed5.f"
	    *dlam = d__[2] + tau;
#line 161 "slaed5.f"
	    delta[1] = -z__[1] / (del + tau);
#line 162 "slaed5.f"
	    delta[2] = -z__[2] / tau;
#line 163 "slaed5.f"
	}
#line 164 "slaed5.f"
	temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
#line 165 "slaed5.f"
	delta[1] /= temp;
#line 166 "slaed5.f"
	delta[2] /= temp;
#line 167 "slaed5.f"
    } else {

/*     Now I=2 */

#line 171 "slaed5.f"
	b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
#line 172 "slaed5.f"
	c__ = *rho * z__[2] * z__[2] * del;
#line 173 "slaed5.f"
	if (b > 0.) {
#line 174 "slaed5.f"
	    tau = (b + sqrt(b * b + c__ * 4.)) / 2.;
#line 175 "slaed5.f"
	} else {
#line 176 "slaed5.f"
	    tau = c__ * 2. / (-b + sqrt(b * b + c__ * 4.));
#line 177 "slaed5.f"
	}
#line 178 "slaed5.f"
	*dlam = d__[2] + tau;
#line 179 "slaed5.f"
	delta[1] = -z__[1] / (del + tau);
#line 180 "slaed5.f"
	delta[2] = -z__[2] / tau;
#line 181 "slaed5.f"
	temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
#line 182 "slaed5.f"
	delta[1] /= temp;
#line 183 "slaed5.f"
	delta[2] /= temp;
#line 184 "slaed5.f"
    }
#line 185 "slaed5.f"
    return 0;

/*     End OF SLAED5 */

} /* slaed5_ */

