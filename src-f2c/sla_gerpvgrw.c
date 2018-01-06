#line 1 "sla_gerpvgrw.f"
/* sla_gerpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gerpvgrw.f"
/* > \brief \b SLA_GERPVGRW */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GERPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_ger
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_ger
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_ger
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, NCOLS, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), AF( LDAF, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_GERPVGRW computes the reciprocal pivot growth factor */
/* > norm(A)/norm(U). The "max absolute element" norm is used. If this is */
/* > much less than 1, the stability of the LU factorization of the */
/* > (equilibrated) matrix A could be poor. This also means that the */
/* > solution X, estimated condition numbers, and error bounds could be */
/* > unreliable. */
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
/* > \param[in] NCOLS */
/* > \verbatim */
/* >          NCOLS is INTEGER */
/* >     The number of columns of the matrix A. NCOLS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >     The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is REAL array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by SGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
doublereal sla_gerpvgrw__(integer *n, integer *ncols, doublereal *a, integer *
	lda, doublereal *af, integer *ldaf)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal amax, umax, rpvgrw;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 123 "sla_gerpvgrw.f"
    /* Parameter adjustments */
#line 123 "sla_gerpvgrw.f"
    a_dim1 = *lda;
#line 123 "sla_gerpvgrw.f"
    a_offset = 1 + a_dim1;
#line 123 "sla_gerpvgrw.f"
    a -= a_offset;
#line 123 "sla_gerpvgrw.f"
    af_dim1 = *ldaf;
#line 123 "sla_gerpvgrw.f"
    af_offset = 1 + af_dim1;
#line 123 "sla_gerpvgrw.f"
    af -= af_offset;
#line 123 "sla_gerpvgrw.f"

#line 123 "sla_gerpvgrw.f"
    /* Function Body */
#line 123 "sla_gerpvgrw.f"
    rpvgrw = 1.;
#line 125 "sla_gerpvgrw.f"
    i__1 = *ncols;
#line 125 "sla_gerpvgrw.f"
    for (j = 1; j <= i__1; ++j) {
#line 126 "sla_gerpvgrw.f"
	amax = 0.;
#line 127 "sla_gerpvgrw.f"
	umax = 0.;
#line 128 "sla_gerpvgrw.f"
	i__2 = *n;
#line 128 "sla_gerpvgrw.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 129 "sla_gerpvgrw.f"
	    d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 129 "sla_gerpvgrw.f"
	    amax = max(d__2,amax);
#line 130 "sla_gerpvgrw.f"
	}
#line 131 "sla_gerpvgrw.f"
	i__2 = j;
#line 131 "sla_gerpvgrw.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 132 "sla_gerpvgrw.f"
	    d__2 = (d__1 = af[i__ + j * af_dim1], abs(d__1));
#line 132 "sla_gerpvgrw.f"
	    umax = max(d__2,umax);
#line 133 "sla_gerpvgrw.f"
	}
#line 134 "sla_gerpvgrw.f"
	if (umax != 0.) {
/* Computing MIN */
#line 135 "sla_gerpvgrw.f"
	    d__1 = amax / umax;
#line 135 "sla_gerpvgrw.f"
	    rpvgrw = min(d__1,rpvgrw);
#line 136 "sla_gerpvgrw.f"
	}
#line 137 "sla_gerpvgrw.f"
    }
#line 138 "sla_gerpvgrw.f"
    ret_val = rpvgrw;
#line 139 "sla_gerpvgrw.f"
    return ret_val;
} /* sla_gerpvgrw__ */

