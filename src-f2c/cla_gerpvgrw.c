#line 1 "cla_gerpvgrw.f"
/* cla_gerpvgrw.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gerpvgrw.f"
/* > \brief \b CLA_GERPVGRW multiplies a square real matrix by a complex matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GERPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger
pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger
pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger
pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N, NCOLS, LDA, LDAF */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_GERPVGRW computes the reciprocal pivot growth factor */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by CGETRF. */
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

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
doublereal cla_gerpvgrw__(integer *n, integer *ncols, doublecomplex *a, 
	integer *lda, doublecomplex *af, integer *ldaf)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal amax, umax, rpvgrw;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 131 "cla_gerpvgrw.f"
    /* Parameter adjustments */
#line 131 "cla_gerpvgrw.f"
    a_dim1 = *lda;
#line 131 "cla_gerpvgrw.f"
    a_offset = 1 + a_dim1;
#line 131 "cla_gerpvgrw.f"
    a -= a_offset;
#line 131 "cla_gerpvgrw.f"
    af_dim1 = *ldaf;
#line 131 "cla_gerpvgrw.f"
    af_offset = 1 + af_dim1;
#line 131 "cla_gerpvgrw.f"
    af -= af_offset;
#line 131 "cla_gerpvgrw.f"

#line 131 "cla_gerpvgrw.f"
    /* Function Body */
#line 131 "cla_gerpvgrw.f"
    rpvgrw = 1.;
#line 133 "cla_gerpvgrw.f"
    i__1 = *ncols;
#line 133 "cla_gerpvgrw.f"
    for (j = 1; j <= i__1; ++j) {
#line 134 "cla_gerpvgrw.f"
	amax = 0.;
#line 135 "cla_gerpvgrw.f"
	umax = 0.;
#line 136 "cla_gerpvgrw.f"
	i__2 = *n;
#line 136 "cla_gerpvgrw.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 137 "cla_gerpvgrw.f"
	    i__3 = i__ + j * a_dim1;
#line 137 "cla_gerpvgrw.f"
	    d__3 = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j *
		     a_dim1]), abs(d__2));
#line 137 "cla_gerpvgrw.f"
	    amax = max(d__3,amax);
#line 138 "cla_gerpvgrw.f"
	}
#line 139 "cla_gerpvgrw.f"
	i__2 = j;
#line 139 "cla_gerpvgrw.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 140 "cla_gerpvgrw.f"
	    i__3 = i__ + j * af_dim1;
#line 140 "cla_gerpvgrw.f"
	    d__3 = (d__1 = af[i__3].r, abs(d__1)) + (d__2 = d_imag(&af[i__ + 
		    j * af_dim1]), abs(d__2));
#line 140 "cla_gerpvgrw.f"
	    umax = max(d__3,umax);
#line 141 "cla_gerpvgrw.f"
	}
#line 142 "cla_gerpvgrw.f"
	if (umax != 0.) {
/* Computing MIN */
#line 143 "cla_gerpvgrw.f"
	    d__1 = amax / umax;
#line 143 "cla_gerpvgrw.f"
	    rpvgrw = min(d__1,rpvgrw);
#line 144 "cla_gerpvgrw.f"
	}
#line 145 "cla_gerpvgrw.f"
    }
#line 146 "cla_gerpvgrw.f"
    ret_val = rpvgrw;
#line 147 "cla_gerpvgrw.f"
    return ret_val;
} /* cla_gerpvgrw__ */

