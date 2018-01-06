#line 1 "sla_gbrpvgrw.f"
/* sla_gbrpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gbrpvgrw.f"
/* > \brief \b SLA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded m
atrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GBRPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbr
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbr
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbr
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB, */
/*                                   LDAFB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), AFB( LDAFB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_GBRPVGRW computes the reciprocal pivot growth factor */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >     The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >     The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCOLS */
/* > \verbatim */
/* >          NCOLS is INTEGER */
/* >     The number of columns of the matrix A.  NCOLS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >     On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >     The j-th column of A is stored in the j-th column of the */
/* >     array AB as follows: */
/* >     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl) */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >     The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* >          AFB is REAL array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by SGBTRF.  U is stored as an upper triangular */
/* >     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* >     and the multipliers used during the factorization are stored */
/* >     in rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
doublereal sla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, j, kd;
    static doublereal amax, umax, rpvgrw;


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
/*     .. Executable Statements .. */

#line 143 "sla_gbrpvgrw.f"
    /* Parameter adjustments */
#line 143 "sla_gbrpvgrw.f"
    ab_dim1 = *ldab;
#line 143 "sla_gbrpvgrw.f"
    ab_offset = 1 + ab_dim1;
#line 143 "sla_gbrpvgrw.f"
    ab -= ab_offset;
#line 143 "sla_gbrpvgrw.f"
    afb_dim1 = *ldafb;
#line 143 "sla_gbrpvgrw.f"
    afb_offset = 1 + afb_dim1;
#line 143 "sla_gbrpvgrw.f"
    afb -= afb_offset;
#line 143 "sla_gbrpvgrw.f"

#line 143 "sla_gbrpvgrw.f"
    /* Function Body */
#line 143 "sla_gbrpvgrw.f"
    rpvgrw = 1.;
#line 145 "sla_gbrpvgrw.f"
    kd = *ku + 1;
#line 146 "sla_gbrpvgrw.f"
    i__1 = *ncols;
#line 146 "sla_gbrpvgrw.f"
    for (j = 1; j <= i__1; ++j) {
#line 147 "sla_gbrpvgrw.f"
	amax = 0.;
#line 148 "sla_gbrpvgrw.f"
	umax = 0.;
/* Computing MAX */
#line 149 "sla_gbrpvgrw.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 149 "sla_gbrpvgrw.f"
	i__4 = j + *kl;
#line 149 "sla_gbrpvgrw.f"
	i__3 = min(i__4,*n);
#line 149 "sla_gbrpvgrw.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 150 "sla_gbrpvgrw.f"
	    d__2 = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(d__1));
#line 150 "sla_gbrpvgrw.f"
	    amax = max(d__2,amax);
#line 151 "sla_gbrpvgrw.f"
	}
/* Computing MAX */
#line 152 "sla_gbrpvgrw.f"
	i__3 = j - *ku;
#line 152 "sla_gbrpvgrw.f"
	i__2 = j;
#line 152 "sla_gbrpvgrw.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 153 "sla_gbrpvgrw.f"
	    d__2 = (d__1 = afb[kd + i__ - j + j * afb_dim1], abs(d__1));
#line 153 "sla_gbrpvgrw.f"
	    umax = max(d__2,umax);
#line 154 "sla_gbrpvgrw.f"
	}
#line 155 "sla_gbrpvgrw.f"
	if (umax != 0.) {
/* Computing MIN */
#line 156 "sla_gbrpvgrw.f"
	    d__1 = amax / umax;
#line 156 "sla_gbrpvgrw.f"
	    rpvgrw = min(d__1,rpvgrw);
#line 157 "sla_gbrpvgrw.f"
	}
#line 158 "sla_gbrpvgrw.f"
    }
#line 159 "sla_gbrpvgrw.f"
    ret_val = rpvgrw;
#line 160 "sla_gbrpvgrw.f"
    return ret_val;
} /* sla_gbrpvgrw__ */

