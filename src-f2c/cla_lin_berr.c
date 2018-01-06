#line 1 "cla_lin_berr.f"
/* cla_lin_berr.f -- translated by f2c (version 20100827).
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

#line 1 "cla_lin_berr.f"
/* > \brief \b CLA_LIN_BERR computes a component-wise relative backward error. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_LIN_BERR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_lin
_berr.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_lin
_berr.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_lin
_berr.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, NZ, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AYB( N, NRHS ), BERR( NRHS ) */
/*       COMPLEX            RES( N, NRHS ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLA_LIN_BERR computes componentwise relative backward error from */
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
/* >          RES is COMPLEX array, dimension (N,NRHS) */
/* >     The residual matrix, i.e., the matrix R in the relative backward */
/* >     error formula above. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* >          AYB is REAL array, dimension (N, NRHS) */
/* >     The denominator in the relative backward error formula above, i.e., */
/* >     the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B */
/* >     are from iterative refinement (see cla_gerfsx_extended.f). */
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

/* > \date June 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cla_lin_berr__(integer *n, integer *nz, integer *nrhs, 
	doublecomplex *res, doublereal *ayb, doublereal *berr)
{
    /* System generated locals */
    integer ayb_dim1, ayb_offset, res_dim1, res_offset, i__1, i__2, i__3, 
	    i__4;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal tmp, safe1;
    extern doublereal slamch_(char *, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Adding SAFE1 to the numerator guards against spuriously zero */
/*     residuals.  A similar safeguard is in the CLA_yyAMV routine used */
/*     to compute AYB. */

#line 144 "cla_lin_berr.f"
    /* Parameter adjustments */
#line 144 "cla_lin_berr.f"
    --berr;
#line 144 "cla_lin_berr.f"
    ayb_dim1 = *n;
#line 144 "cla_lin_berr.f"
    ayb_offset = 1 + ayb_dim1;
#line 144 "cla_lin_berr.f"
    ayb -= ayb_offset;
#line 144 "cla_lin_berr.f"
    res_dim1 = *n;
#line 144 "cla_lin_berr.f"
    res_offset = 1 + res_dim1;
#line 144 "cla_lin_berr.f"
    res -= res_offset;
#line 144 "cla_lin_berr.f"

#line 144 "cla_lin_berr.f"
    /* Function Body */
#line 144 "cla_lin_berr.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 145 "cla_lin_berr.f"
    safe1 = (*nz + 1) * safe1;
#line 147 "cla_lin_berr.f"
    i__1 = *nrhs;
#line 147 "cla_lin_berr.f"
    for (j = 1; j <= i__1; ++j) {
#line 148 "cla_lin_berr.f"
	berr[j] = 0.;
#line 149 "cla_lin_berr.f"
	i__2 = *n;
#line 149 "cla_lin_berr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 150 "cla_lin_berr.f"
	    if (ayb[i__ + j * ayb_dim1] != 0.) {
#line 151 "cla_lin_berr.f"
		i__3 = i__ + j * res_dim1;
#line 151 "cla_lin_berr.f"
		d__3 = (d__1 = res[i__3].r, abs(d__1)) + (d__2 = d_imag(&res[
			i__ + j * res_dim1]), abs(d__2));
#line 151 "cla_lin_berr.f"
		z__3.r = d__3, z__3.i = 0.;
#line 151 "cla_lin_berr.f"
		z__2.r = safe1 + z__3.r, z__2.i = z__3.i;
#line 151 "cla_lin_berr.f"
		i__4 = i__ + j * ayb_dim1;
#line 151 "cla_lin_berr.f"
		z__1.r = z__2.r / ayb[i__4], z__1.i = z__2.i / ayb[i__4];
#line 151 "cla_lin_berr.f"
		tmp = z__1.r;
/* Computing MAX */
#line 152 "cla_lin_berr.f"
		d__1 = berr[j];
#line 152 "cla_lin_berr.f"
		berr[j] = max(d__1,tmp);
#line 153 "cla_lin_berr.f"
	    }

/*     If AYB is exactly 0.0 (and if computed by CLA_yyAMV), then we know */
/*     the true residual also must be exactly 0.0. */

#line 158 "cla_lin_berr.f"
	}
#line 159 "cla_lin_berr.f"
    }
#line 160 "cla_lin_berr.f"
    return 0;
} /* cla_lin_berr__ */

