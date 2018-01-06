#line 1 "slaqr1.f"
/* slaqr1.f -- translated by f2c (version 20100827).
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

#line 1 "slaqr1.f"
/* > \brief \b SLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H a
nd specified shifts. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQR1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V ) */

/*       .. Scalar Arguments .. */
/*       REAL               SI1, SI2, SR1, SR2 */
/*       INTEGER            LDH, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               H( LDH, * ), V( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >      Given a 2-by-2 or 3-by-3 matrix H, SLAQR1 sets v to a */
/* >      scalar multiple of the first column of the product */
/* > */
/* >      (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I) */
/* > */
/* >      scaling to avoid overflows and most underflows. It */
/* >      is assumed that either */
/* > */
/* >              1) sr1 = sr2 and si1 = -si2 */
/* >          or */
/* >              2) si1 = si2 = 0. */
/* > */
/* >      This is useful for starting double implicit shift bulges */
/* >      in the QR algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is integer */
/* >              Order of the matrix H. N must be either 2 or 3. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is REAL array of dimension (LDH,N) */
/* >              The 2-by-2 or 3-by-3 matrix H in (*). */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is integer */
/* >              The leading dimension of H as declared in */
/* >              the calling procedure.  LDH.GE.N */
/* > \endverbatim */
/* > */
/* > \param[in] SR1 */
/* > \verbatim */
/* >          SR1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] SI1 */
/* > \verbatim */
/* >          SI1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] SR2 */
/* > \verbatim */
/* >          SR2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] SI2 */
/* > \verbatim */
/* >          SI2 is REAL */
/* >              The shifts in (*). */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is REAL array of dimension N */
/* >              A scalar multiple of the first column of the */
/* >              matrix K in (*). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaqr1_(integer *n, doublereal *h__, integer *ldh, 
	doublereal *sr1, doublereal *si1, doublereal *sr2, doublereal *si2, 
	doublereal *v)
{
    /* System generated locals */
    integer h_dim1, h_offset;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal s, h21s, h31s;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ================================================================ */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 150 "slaqr1.f"
    /* Parameter adjustments */
#line 150 "slaqr1.f"
    h_dim1 = *ldh;
#line 150 "slaqr1.f"
    h_offset = 1 + h_dim1;
#line 150 "slaqr1.f"
    h__ -= h_offset;
#line 150 "slaqr1.f"
    --v;
#line 150 "slaqr1.f"

#line 150 "slaqr1.f"
    /* Function Body */
#line 150 "slaqr1.f"
    if (*n == 2) {
#line 151 "slaqr1.f"
	s = (d__1 = h__[h_dim1 + 1] - *sr2, abs(d__1)) + abs(*si2) + (d__2 = 
		h__[h_dim1 + 2], abs(d__2));
#line 152 "slaqr1.f"
	if (s == 0.) {
#line 153 "slaqr1.f"
	    v[1] = 0.;
#line 154 "slaqr1.f"
	    v[2] = 0.;
#line 155 "slaqr1.f"
	} else {
#line 156 "slaqr1.f"
	    h21s = h__[h_dim1 + 2] / s;
#line 157 "slaqr1.f"
	    v[1] = h21s * h__[(h_dim1 << 1) + 1] + (h__[h_dim1 + 1] - *sr1) * 
		    ((h__[h_dim1 + 1] - *sr2) / s) - *si1 * (*si2 / s);
#line 159 "slaqr1.f"
	    v[2] = h21s * (h__[h_dim1 + 1] + h__[(h_dim1 << 1) + 2] - *sr1 - *
		    sr2);
#line 160 "slaqr1.f"
	}
#line 161 "slaqr1.f"
    } else {
#line 162 "slaqr1.f"
	s = (d__1 = h__[h_dim1 + 1] - *sr2, abs(d__1)) + abs(*si2) + (d__2 = 
		h__[h_dim1 + 2], abs(d__2)) + (d__3 = h__[h_dim1 + 3], abs(
		d__3));
#line 164 "slaqr1.f"
	if (s == 0.) {
#line 165 "slaqr1.f"
	    v[1] = 0.;
#line 166 "slaqr1.f"
	    v[2] = 0.;
#line 167 "slaqr1.f"
	    v[3] = 0.;
#line 168 "slaqr1.f"
	} else {
#line 169 "slaqr1.f"
	    h21s = h__[h_dim1 + 2] / s;
#line 170 "slaqr1.f"
	    h31s = h__[h_dim1 + 3] / s;
#line 171 "slaqr1.f"
	    v[1] = (h__[h_dim1 + 1] - *sr1) * ((h__[h_dim1 + 1] - *sr2) / s) 
		    - *si1 * (*si2 / s) + h__[(h_dim1 << 1) + 1] * h21s + h__[
		    h_dim1 * 3 + 1] * h31s;
#line 173 "slaqr1.f"
	    v[2] = h21s * (h__[h_dim1 + 1] + h__[(h_dim1 << 1) + 2] - *sr1 - *
		    sr2) + h__[h_dim1 * 3 + 2] * h31s;
#line 175 "slaqr1.f"
	    v[3] = h31s * (h__[h_dim1 + 1] + h__[h_dim1 * 3 + 3] - *sr1 - *
		    sr2) + h21s * h__[(h_dim1 << 1) + 3];
#line 177 "slaqr1.f"
	}
#line 178 "slaqr1.f"
    }
#line 179 "slaqr1.f"
    return 0;
} /* slaqr1_ */

