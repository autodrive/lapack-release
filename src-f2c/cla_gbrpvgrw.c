#line 1 "cla_gbrpvgrw.f"
/* cla_gbrpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gbrpvgrw.f"
/* > \brief \b CLA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded m
atrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GBRPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbr
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbr
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbr
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB, */
/*                                   LDAFB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_GBRPVGRW computes the reciprocal pivot growth factor */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          AFB is COMPLEX array, dimension (LDAFB,N) */
/* >     Details of the LU factorization of the band matrix A, as */
/* >     computed by CGBTRF.  U is stored as an upper triangular */
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

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
doublereal cla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
	ldafb)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 150 "cla_gbrpvgrw.f"
    /* Parameter adjustments */
#line 150 "cla_gbrpvgrw.f"
    ab_dim1 = *ldab;
#line 150 "cla_gbrpvgrw.f"
    ab_offset = 1 + ab_dim1;
#line 150 "cla_gbrpvgrw.f"
    ab -= ab_offset;
#line 150 "cla_gbrpvgrw.f"
    afb_dim1 = *ldafb;
#line 150 "cla_gbrpvgrw.f"
    afb_offset = 1 + afb_dim1;
#line 150 "cla_gbrpvgrw.f"
    afb -= afb_offset;
#line 150 "cla_gbrpvgrw.f"

#line 150 "cla_gbrpvgrw.f"
    /* Function Body */
#line 150 "cla_gbrpvgrw.f"
    rpvgrw = 1.;
#line 152 "cla_gbrpvgrw.f"
    kd = *ku + 1;
#line 153 "cla_gbrpvgrw.f"
    i__1 = *ncols;
#line 153 "cla_gbrpvgrw.f"
    for (j = 1; j <= i__1; ++j) {
#line 154 "cla_gbrpvgrw.f"
	amax = 0.;
#line 155 "cla_gbrpvgrw.f"
	umax = 0.;
/* Computing MAX */
#line 156 "cla_gbrpvgrw.f"
	i__2 = j - *ku;
/* Computing MIN */
#line 156 "cla_gbrpvgrw.f"
	i__4 = j + *kl;
#line 156 "cla_gbrpvgrw.f"
	i__3 = min(i__4,*n);
#line 156 "cla_gbrpvgrw.f"
	for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 157 "cla_gbrpvgrw.f"
	    i__2 = kd + i__ - j + j * ab_dim1;
#line 157 "cla_gbrpvgrw.f"
	    d__3 = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = d_imag(&ab[kd + 
		    i__ - j + j * ab_dim1]), abs(d__2));
#line 157 "cla_gbrpvgrw.f"
	    amax = max(d__3,amax);
#line 158 "cla_gbrpvgrw.f"
	}
/* Computing MAX */
#line 159 "cla_gbrpvgrw.f"
	i__3 = j - *ku;
#line 159 "cla_gbrpvgrw.f"
	i__2 = j;
#line 159 "cla_gbrpvgrw.f"
	for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
/* Computing MAX */
#line 160 "cla_gbrpvgrw.f"
	    i__3 = kd + i__ - j + j * afb_dim1;
#line 160 "cla_gbrpvgrw.f"
	    d__3 = (d__1 = afb[i__3].r, abs(d__1)) + (d__2 = d_imag(&afb[kd + 
		    i__ - j + j * afb_dim1]), abs(d__2));
#line 160 "cla_gbrpvgrw.f"
	    umax = max(d__3,umax);
#line 161 "cla_gbrpvgrw.f"
	}
#line 162 "cla_gbrpvgrw.f"
	if (umax != 0.) {
/* Computing MIN */
#line 163 "cla_gbrpvgrw.f"
	    d__1 = amax / umax;
#line 163 "cla_gbrpvgrw.f"
	    rpvgrw = min(d__1,rpvgrw);
#line 164 "cla_gbrpvgrw.f"
	}
#line 165 "cla_gbrpvgrw.f"
    }
#line 166 "cla_gbrpvgrw.f"
    ret_val = rpvgrw;
#line 167 "cla_gbrpvgrw.f"
    return ret_val;
} /* cla_gbrpvgrw__ */

