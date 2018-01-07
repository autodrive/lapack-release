#line 1 "crotg.f"
/* crotg.f -- translated by f2c (version 20100827).
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

#line 1 "crotg.f"
/* > \brief \b CROTG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CROTG(CA,CB,C,S) */

/*       .. Scalar Arguments .. */
/*       COMPLEX CA,CB,S */
/*       REAL C */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CROTG determines a complex Givens rotation. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int crotg_(doublecomplex *ca, doublecomplex *cb, doublereal *
	c__, doublecomplex *s)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal norm;
    static doublecomplex alpha;
    static doublereal scale;


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 61 "crotg.f"
    if (z_abs(ca) == 0.) {
#line 62 "crotg.f"
	*c__ = 0.;
#line 63 "crotg.f"
	s->r = 1., s->i = 0.;
#line 64 "crotg.f"
	ca->r = cb->r, ca->i = cb->i;
#line 65 "crotg.f"
    } else {
#line 66 "crotg.f"
	scale = z_abs(ca) + z_abs(cb);
#line 67 "crotg.f"
	z__1.r = ca->r / scale, z__1.i = ca->i / scale;
/* Computing 2nd power */
#line 67 "crotg.f"
	d__1 = z_abs(&z__1);
#line 67 "crotg.f"
	z__2.r = cb->r / scale, z__2.i = cb->i / scale;
/* Computing 2nd power */
#line 67 "crotg.f"
	d__2 = z_abs(&z__2);
#line 67 "crotg.f"
	norm = scale * sqrt(d__1 * d__1 + d__2 * d__2);
#line 68 "crotg.f"
	d__1 = z_abs(ca);
#line 68 "crotg.f"
	z__1.r = ca->r / d__1, z__1.i = ca->i / d__1;
#line 68 "crotg.f"
	alpha.r = z__1.r, alpha.i = z__1.i;
#line 69 "crotg.f"
	*c__ = z_abs(ca) / norm;
#line 70 "crotg.f"
	d_cnjg(&z__3, cb);
#line 70 "crotg.f"
	z__2.r = alpha.r * z__3.r - alpha.i * z__3.i, z__2.i = alpha.r * 
		z__3.i + alpha.i * z__3.r;
#line 70 "crotg.f"
	z__1.r = z__2.r / norm, z__1.i = z__2.i / norm;
#line 70 "crotg.f"
	s->r = z__1.r, s->i = z__1.i;
#line 71 "crotg.f"
	z__1.r = norm * alpha.r, z__1.i = norm * alpha.i;
#line 71 "crotg.f"
	ca->r = z__1.r, ca->i = z__1.i;
#line 72 "crotg.f"
    }
#line 73 "crotg.f"
    return 0;
} /* crotg_ */

