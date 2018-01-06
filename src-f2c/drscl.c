#line 1 "drscl.f"
/* drscl.f -- translated by f2c (version 20100827).
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

#line 1 "drscl.f"
/* > \brief \b DRSCL multiplies a vector by the reciprocal of a real scalar. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DRSCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/drscl.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/drscl.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/drscl.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DRSCL( N, SA, SX, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       DOUBLE PRECISION   SA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   SX( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DRSCL multiplies an n-element real vector x by the real scalar 1/a. */
/* > This is done without overflow or underflow as long as */
/* > the final result x/a does not overflow or underflow. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of components of the vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SA */
/* > \verbatim */
/* >          SA is DOUBLE PRECISION */
/* >          The scalar a which is used to divide each component of x. */
/* >          SA must be >= 0, or the subroutine will divide by zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX */
/* > \verbatim */
/* >          SX is DOUBLE PRECISION array, dimension */
/* >                         (1+(N-1)*abs(INCX)) */
/* >          The n-element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of the vector SX. */
/* >          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int drscl_(integer *n, doublereal *sa, doublereal *sx, 
	integer *incx)
{
    static doublereal mul, cden;
    static logical done;
    static doublereal cnum, cden1, cnum1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 124 "drscl.f"
    /* Parameter adjustments */
#line 124 "drscl.f"
    --sx;
#line 124 "drscl.f"

#line 124 "drscl.f"
    /* Function Body */
#line 124 "drscl.f"
    if (*n <= 0) {
#line 124 "drscl.f"
	return 0;
#line 124 "drscl.f"
    }

/*     Get machine parameters */

#line 129 "drscl.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 130 "drscl.f"
    bignum = 1. / smlnum;
#line 131 "drscl.f"
    dlabad_(&smlnum, &bignum);

/*     Initialize the denominator to SA and the numerator to 1. */

#line 135 "drscl.f"
    cden = *sa;
#line 136 "drscl.f"
    cnum = 1.;

#line 138 "drscl.f"
L10:
#line 139 "drscl.f"
    cden1 = cden * smlnum;
#line 140 "drscl.f"
    cnum1 = cnum / bignum;
#line 141 "drscl.f"
    if (abs(cden1) > abs(cnum) && cnum != 0.) {

/*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */

#line 145 "drscl.f"
	mul = smlnum;
#line 146 "drscl.f"
	done = FALSE_;
#line 147 "drscl.f"
	cden = cden1;
#line 148 "drscl.f"
    } else if (abs(cnum1) > abs(cden)) {

/*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

#line 152 "drscl.f"
	mul = bignum;
#line 153 "drscl.f"
	done = FALSE_;
#line 154 "drscl.f"
	cnum = cnum1;
#line 155 "drscl.f"
    } else {

/*        Multiply X by CNUM / CDEN and return. */

#line 159 "drscl.f"
	mul = cnum / cden;
#line 160 "drscl.f"
	done = TRUE_;
#line 161 "drscl.f"
    }

/*     Scale the vector X by MUL */

#line 165 "drscl.f"
    dscal_(n, &mul, &sx[1], incx);

#line 167 "drscl.f"
    if (! done) {
#line 167 "drscl.f"
	goto L10;
#line 167 "drscl.f"
    }

#line 170 "drscl.f"
    return 0;

/*     End of DRSCL */

} /* drscl_ */

