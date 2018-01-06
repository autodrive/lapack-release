#line 1 "cladiv.f"
/* cladiv.f -- translated by f2c (version 20100827).
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

#line 1 "cladiv.f"
/* > \brief \b CLADIV performs complex division in real arithmetic, avoiding unnecessary overflow. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLADIV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cladiv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cladiv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cladiv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       COMPLEX FUNCTION CLADIV( X, Y ) */

/*       .. Scalar Arguments .. */
/*       COMPLEX            X, Y */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLADIV := X / Y, where X and Y are complex.  The computation of X / Y */
/* > will not overflow on an intermediary step unless the results */
/* > overflows. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is COMPLEX */
/* >          The complex scalars X and Y. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Double Complex */ VOID cladiv_(doublecomplex * ret_val, doublecomplex *x, 
	doublecomplex *y)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal zi, zr;
    extern /* Subroutine */ int sladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 89 "cladiv.f"
    d__1 = x->r;
#line 89 "cladiv.f"
    d__2 = d_imag(x);
#line 89 "cladiv.f"
    d__3 = y->r;
#line 89 "cladiv.f"
    d__4 = d_imag(y);
#line 89 "cladiv.f"
    sladiv_(&d__1, &d__2, &d__3, &d__4, &zr, &zi);
#line 91 "cladiv.f"
    z__1.r = zr, z__1.i = zi;
#line 91 "cladiv.f"
     ret_val->r = z__1.r,  ret_val->i = z__1.i;

#line 93 "cladiv.f"
    return ;

/*     End of CLADIV */

} /* cladiv_ */

