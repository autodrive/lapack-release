#line 1 "sla_lin_berr.f"
/* sla_lin_berr.f -- translated by f2c (version 20100827).
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

#line 1 "sla_lin_berr.f"
/* > \brief \b SLA_LIN_BERR computes a component-wise relative backward error. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_LIN_BERR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_lin
_berr.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_lin
_berr.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_lin
_berr.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, NZ, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AYB( N, NRHS ), BERR( NRHS ) */
/*       REAL               RES( N, NRHS ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLA_LIN_BERR computes componentwise relative backward error from */
/* >    the formula */
/* >        max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/* >    where abs(Z) is the componentwise absolute value of the matrix */
/* >    or vector Z. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >     The number of linear equations, i.e., the order of the */
/* >     matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NZ */
/* > \verbatim */
/* >          NZ is INTEGER */
/* >     We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to */
/* >     guard against spuriously zero residuals. Default value is N. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >     The number of right hand sides, i.e., the number of columns */
/* >     of the matrices AYB, RES, and BERR.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] RES */
/* > \verbatim */
/* >          RES is REAL array, dimension (N,NRHS) */
/* >     The residual matrix, i.e., the matrix R in the relative backward */
/* >     error formula above. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* >          AYB is REAL array, dimension (N, NRHS) */
/* >     The denominator in the relative backward error formula above, i.e., */
/* >     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B */
/* >     are from iterative refinement (see sla_gerfsx_extended.f). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is REAL array, dimension (NRHS) */
/* >     The componentwise relative backward error from the formula above. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sla_lin_berr__(integer *n, integer *nz, integer *nrhs, 
	doublereal *res, doublereal *ayb, doublereal *berr)
{
    /* System generated locals */
    integer ayb_dim1, ayb_offset, res_dim1, res_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal tmp, safe1;
    extern doublereal slamch_(char *, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Adding SAFE1 to the numerator guards against spuriously zero */
/*     residuals.  A similar safeguard is in the SLA_yyAMV routine used */
/*     to compute AYB. */

#line 137 "sla_lin_berr.f"
    /* Parameter adjustments */
#line 137 "sla_lin_berr.f"
    --berr;
#line 137 "sla_lin_berr.f"
    ayb_dim1 = *n;
#line 137 "sla_lin_berr.f"
    ayb_offset = 1 + ayb_dim1;
#line 137 "sla_lin_berr.f"
    ayb -= ayb_offset;
#line 137 "sla_lin_berr.f"
    res_dim1 = *n;
#line 137 "sla_lin_berr.f"
    res_offset = 1 + res_dim1;
#line 137 "sla_lin_berr.f"
    res -= res_offset;
#line 137 "sla_lin_berr.f"

#line 137 "sla_lin_berr.f"
    /* Function Body */
#line 137 "sla_lin_berr.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 138 "sla_lin_berr.f"
    safe1 = (*nz + 1) * safe1;
#line 140 "sla_lin_berr.f"
    i__1 = *nrhs;
#line 140 "sla_lin_berr.f"
    for (j = 1; j <= i__1; ++j) {
#line 141 "sla_lin_berr.f"
	berr[j] = 0.;
#line 142 "sla_lin_berr.f"
	i__2 = *n;
#line 142 "sla_lin_berr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 143 "sla_lin_berr.f"
	    if (ayb[i__ + j * ayb_dim1] != 0.) {
#line 144 "sla_lin_berr.f"
		tmp = (safe1 + (d__1 = res[i__ + j * res_dim1], abs(d__1))) / 
			ayb[i__ + j * ayb_dim1];
/* Computing MAX */
#line 145 "sla_lin_berr.f"
		d__1 = berr[j];
#line 145 "sla_lin_berr.f"
		berr[j] = max(d__1,tmp);
#line 146 "sla_lin_berr.f"
	    }

/*     If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know */
/*     the true residual also must be exactly 0.0. */

#line 151 "sla_lin_berr.f"
	}
#line 152 "sla_lin_berr.f"
    }
#line 153 "sla_lin_berr.f"
    return 0;
} /* sla_lin_berr__ */

