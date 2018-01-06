#line 1 "dlar2v.f"
/* dlar2v.f -- translated by f2c (version 20100827).
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

#line 1 "dlar2v.f"
/* > \brief \b DLAR2V applies a vector of plane rotations with real cosines and real sines from both sides to 
a sequence of 2-by-2 symmetric/Hermitian matrices. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAR2V + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlar2v.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlar2v.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlar2v.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAR2V( N, X, Y, Z, INCX, C, S, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ), S( * ), X( * ), Y( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAR2V applies a vector of real plane rotations from both sides to */
/* > a sequence of 2-by-2 real symmetric matrices, defined by the elements */
/* > of the vectors x, y and z. For i = 1,2,...,n */
/* > */
/* >    ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) ) */
/* >    ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of plane rotations to be applied. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, */
/* >                         dimension (1+(N-1)*INCX) */
/* >          The vector x. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, */
/* >                         dimension (1+(N-1)*INCX) */
/* >          The vector y. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, */
/* >                         dimension (1+(N-1)*INCX) */
/* >          The vector z. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between elements of X, Y and Z. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* >          The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* >          The sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* >          INCC is INTEGER */
/* >          The increment between elements of C and S. INCC > 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlar2v_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *incx, doublereal *c__, doublereal *s, 
	integer *incc)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t1, t2, t3, t4, t5, t6;
    static integer ic;
    static doublereal ci, si;
    static integer ix;
    static doublereal xi, yi, zi;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
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
/*     .. Executable Statements .. */

#line 133 "dlar2v.f"
    /* Parameter adjustments */
#line 133 "dlar2v.f"
    --s;
#line 133 "dlar2v.f"
    --c__;
#line 133 "dlar2v.f"
    --z__;
#line 133 "dlar2v.f"
    --y;
#line 133 "dlar2v.f"
    --x;
#line 133 "dlar2v.f"

#line 133 "dlar2v.f"
    /* Function Body */
#line 133 "dlar2v.f"
    ix = 1;
#line 134 "dlar2v.f"
    ic = 1;
#line 135 "dlar2v.f"
    i__1 = *n;
#line 135 "dlar2v.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 136 "dlar2v.f"
	xi = x[ix];
#line 137 "dlar2v.f"
	yi = y[ix];
#line 138 "dlar2v.f"
	zi = z__[ix];
#line 139 "dlar2v.f"
	ci = c__[ic];
#line 140 "dlar2v.f"
	si = s[ic];
#line 141 "dlar2v.f"
	t1 = si * zi;
#line 142 "dlar2v.f"
	t2 = ci * zi;
#line 143 "dlar2v.f"
	t3 = t2 - si * xi;
#line 144 "dlar2v.f"
	t4 = t2 + si * yi;
#line 145 "dlar2v.f"
	t5 = ci * xi + t1;
#line 146 "dlar2v.f"
	t6 = ci * yi - t1;
#line 147 "dlar2v.f"
	x[ix] = ci * t5 + si * t4;
#line 148 "dlar2v.f"
	y[ix] = ci * t6 - si * t3;
#line 149 "dlar2v.f"
	z__[ix] = ci * t4 - si * t5;
#line 150 "dlar2v.f"
	ix += *incx;
#line 151 "dlar2v.f"
	ic += *incc;
#line 152 "dlar2v.f"
/* L10: */
#line 152 "dlar2v.f"
    }

/*     End of DLAR2V */

#line 156 "dlar2v.f"
    return 0;
} /* dlar2v_ */

