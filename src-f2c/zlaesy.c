#line 1 "zlaesy.f"
/* zlaesy.f -- translated by f2c (version 20100827).
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

#line 1 "zlaesy.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__2 = 2;

/* > \brief \b ZLAESY computes the eigenvalues and eigenvectors of a 2-by-2 complex symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAESY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaesy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaesy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaesy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 ) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16         A, B, C, CS1, EVSCAL, RT1, RT2, SN1 */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix */
/* >    ( ( A, B );( B, C ) ) */
/* > provided the norm of the matrix of eigenvectors is larger than */
/* > some threshold value. */
/* > */
/* > RT1 is the eigenvalue of larger absolute value, and RT2 of */
/* > smaller absolute value.  If the eigenvectors are computed, then */
/* > on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence */
/* > */
/* > [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ] */
/* > [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ] */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 */
/* >          The ( 1, 1 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 */
/* >          The ( 1, 2 ) element of input matrix.  The ( 2, 1 ) element */
/* >          is also given by B, since the 2-by-2 matrix is symmetric. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is COMPLEX*16 */
/* >          The ( 2, 2 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* >          RT1 is COMPLEX*16 */
/* >          The eigenvalue of larger modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* >          RT2 is COMPLEX*16 */
/* >          The eigenvalue of smaller modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] EVSCAL */
/* > \verbatim */
/* >          EVSCAL is COMPLEX*16 */
/* >          The complex value by which the eigenvector matrix was scaled */
/* >          to make it orthonormal.  If EVSCAL is zero, the eigenvectors */
/* >          were not computed.  This means one of two things:  the 2-by-2 */
/* >          matrix could not be diagonalized, or the norm of the matrix */
/* >          of eigenvectors before scaling was larger than the threshold */
/* >          value THRESH (set below). */
/* > \endverbatim */
/* > */
/* > \param[out] CS1 */
/* > \verbatim */
/* >          CS1 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[out] SN1 */
/* > \verbatim */
/* >          SN1 is COMPLEX*16 */
/* >          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector */
/* >          for RT1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16SYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaesy_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2, 
	doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void pow_zi(doublecomplex *, doublecomplex *, integer *), z_sqrt(
	    doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex s, t;
    static doublereal z__;
    static doublecomplex tmp;
    static doublereal babs, tabs, evnorm;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */


/*     Special case:  The matrix is actually diagonal. */
/*     To avoid divide by zero later, we treat this case separately. */

#line 154 "zlaesy.f"
    if (z_abs(b) == 0.) {
#line 155 "zlaesy.f"
	rt1->r = a->r, rt1->i = a->i;
#line 156 "zlaesy.f"
	rt2->r = c__->r, rt2->i = c__->i;
#line 157 "zlaesy.f"
	if (z_abs(rt1) < z_abs(rt2)) {
#line 158 "zlaesy.f"
	    tmp.r = rt1->r, tmp.i = rt1->i;
#line 159 "zlaesy.f"
	    rt1->r = rt2->r, rt1->i = rt2->i;
#line 160 "zlaesy.f"
	    rt2->r = tmp.r, rt2->i = tmp.i;
#line 161 "zlaesy.f"
	    cs1->r = 0., cs1->i = 0.;
#line 162 "zlaesy.f"
	    sn1->r = 1., sn1->i = 0.;
#line 163 "zlaesy.f"
	} else {
#line 164 "zlaesy.f"
	    cs1->r = 1., cs1->i = 0.;
#line 165 "zlaesy.f"
	    sn1->r = 0., sn1->i = 0.;
#line 166 "zlaesy.f"
	}
#line 167 "zlaesy.f"
    } else {

/*        Compute the eigenvalues and eigenvectors. */
/*        The characteristic equation is */
/*           lambda **2 - (A+C) lambda + (A*C - B*B) */
/*        and we solve it using the quadratic formula. */

#line 174 "zlaesy.f"
	z__2.r = a->r + c__->r, z__2.i = a->i + c__->i;
#line 174 "zlaesy.f"
	z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
#line 174 "zlaesy.f"
	s.r = z__1.r, s.i = z__1.i;
#line 175 "zlaesy.f"
	z__2.r = a->r - c__->r, z__2.i = a->i - c__->i;
#line 175 "zlaesy.f"
	z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
#line 175 "zlaesy.f"
	t.r = z__1.r, t.i = z__1.i;

/*        Take the square root carefully to avoid over/under flow. */

#line 179 "zlaesy.f"
	babs = z_abs(b);
#line 180 "zlaesy.f"
	tabs = z_abs(&t);
#line 181 "zlaesy.f"
	z__ = max(babs,tabs);
#line 182 "zlaesy.f"
	if (z__ > 0.) {
#line 182 "zlaesy.f"
	    z__5.r = t.r / z__, z__5.i = t.i / z__;
#line 182 "zlaesy.f"
	    pow_zi(&z__4, &z__5, &c__2);
#line 182 "zlaesy.f"
	    z__7.r = b->r / z__, z__7.i = b->i / z__;
#line 182 "zlaesy.f"
	    pow_zi(&z__6, &z__7, &c__2);
#line 182 "zlaesy.f"
	    z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
#line 182 "zlaesy.f"
	    z_sqrt(&z__2, &z__3);
#line 182 "zlaesy.f"
	    z__1.r = z__ * z__2.r, z__1.i = z__ * z__2.i;
#line 182 "zlaesy.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 182 "zlaesy.f"
	}

/*        Compute the two eigenvalues.  RT1 and RT2 are exchanged */
/*        if necessary so that RT1 will have the greater magnitude. */

#line 188 "zlaesy.f"
	z__1.r = s.r + t.r, z__1.i = s.i + t.i;
#line 188 "zlaesy.f"
	rt1->r = z__1.r, rt1->i = z__1.i;
#line 189 "zlaesy.f"
	z__1.r = s.r - t.r, z__1.i = s.i - t.i;
#line 189 "zlaesy.f"
	rt2->r = z__1.r, rt2->i = z__1.i;
#line 190 "zlaesy.f"
	if (z_abs(rt1) < z_abs(rt2)) {
#line 191 "zlaesy.f"
	    tmp.r = rt1->r, tmp.i = rt1->i;
#line 192 "zlaesy.f"
	    rt1->r = rt2->r, rt1->i = rt2->i;
#line 193 "zlaesy.f"
	    rt2->r = tmp.r, rt2->i = tmp.i;
#line 194 "zlaesy.f"
	}

/*        Choose CS1 = 1 and SN1 to satisfy the first equation, then */
/*        scale the components of this eigenvector so that the matrix */
/*        of eigenvectors X satisfies  X * X**T = I .  (No scaling is */
/*        done if the norm of the eigenvalue matrix is less than THRESH.) */

#line 201 "zlaesy.f"
	z__2.r = rt1->r - a->r, z__2.i = rt1->i - a->i;
#line 201 "zlaesy.f"
	z_div(&z__1, &z__2, b);
#line 201 "zlaesy.f"
	sn1->r = z__1.r, sn1->i = z__1.i;
#line 202 "zlaesy.f"
	tabs = z_abs(sn1);
#line 203 "zlaesy.f"
	if (tabs > 1.) {
/* Computing 2nd power */
#line 204 "zlaesy.f"
	    d__2 = 1. / tabs;
#line 204 "zlaesy.f"
	    d__1 = d__2 * d__2;
#line 204 "zlaesy.f"
	    z__5.r = sn1->r / tabs, z__5.i = sn1->i / tabs;
#line 204 "zlaesy.f"
	    pow_zi(&z__4, &z__5, &c__2);
#line 204 "zlaesy.f"
	    z__3.r = d__1 + z__4.r, z__3.i = z__4.i;
#line 204 "zlaesy.f"
	    z_sqrt(&z__2, &z__3);
#line 204 "zlaesy.f"
	    z__1.r = tabs * z__2.r, z__1.i = tabs * z__2.i;
#line 204 "zlaesy.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 205 "zlaesy.f"
	} else {
#line 206 "zlaesy.f"
	    z__3.r = sn1->r * sn1->r - sn1->i * sn1->i, z__3.i = sn1->r * 
		    sn1->i + sn1->i * sn1->r;
#line 206 "zlaesy.f"
	    z__2.r = z__3.r + 1., z__2.i = z__3.i + 0.;
#line 206 "zlaesy.f"
	    z_sqrt(&z__1, &z__2);
#line 206 "zlaesy.f"
	    t.r = z__1.r, t.i = z__1.i;
#line 207 "zlaesy.f"
	}
#line 208 "zlaesy.f"
	evnorm = z_abs(&t);
#line 209 "zlaesy.f"
	if (evnorm >= .1) {
#line 210 "zlaesy.f"
	    z_div(&z__1, &c_b1, &t);
#line 210 "zlaesy.f"
	    evscal->r = z__1.r, evscal->i = z__1.i;
#line 211 "zlaesy.f"
	    cs1->r = evscal->r, cs1->i = evscal->i;
#line 212 "zlaesy.f"
	    z__1.r = sn1->r * evscal->r - sn1->i * evscal->i, z__1.i = sn1->r 
		    * evscal->i + sn1->i * evscal->r;
#line 212 "zlaesy.f"
	    sn1->r = z__1.r, sn1->i = z__1.i;
#line 213 "zlaesy.f"
	} else {
#line 214 "zlaesy.f"
	    evscal->r = 0., evscal->i = 0.;
#line 215 "zlaesy.f"
	}
#line 216 "zlaesy.f"
    }
#line 217 "zlaesy.f"
    return 0;

/*     End of ZLAESY */

} /* zlaesy_ */

