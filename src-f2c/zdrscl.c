#line 1 "zdrscl.f"
/* zdrscl.f -- translated by f2c (version 20100827).
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

#line 1 "zdrscl.f"
/* > \brief \b ZDRSCL multiplies a vector by the reciprocal of a real scalar. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZDRSCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zdrscl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zdrscl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zdrscl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZDRSCL( N, SA, SX, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       DOUBLE PRECISION   SA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         SX( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZDRSCL multiplies an n-element complex vector x by the real scalar */
/* > 1/a.  This is done without overflow or underflow as long as */
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
/* >          SX is COMPLEX*16 array, dimension */
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

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zdrscl_(integer *n, doublereal *sa, doublecomplex *sx, 
	integer *incx)
{
    static doublereal mul, cden;
    static logical done;
    static doublereal cnum, cden1, cnum1;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 124 "zdrscl.f"
    /* Parameter adjustments */
#line 124 "zdrscl.f"
    --sx;
#line 124 "zdrscl.f"

#line 124 "zdrscl.f"
    /* Function Body */
#line 124 "zdrscl.f"
    if (*n <= 0) {
#line 124 "zdrscl.f"
	return 0;
#line 124 "zdrscl.f"
    }

/*     Get machine parameters */

#line 129 "zdrscl.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 130 "zdrscl.f"
    bignum = 1. / smlnum;
#line 131 "zdrscl.f"
    dlabad_(&smlnum, &bignum);

/*     Initialize the denominator to SA and the numerator to 1. */

#line 135 "zdrscl.f"
    cden = *sa;
#line 136 "zdrscl.f"
    cnum = 1.;

#line 138 "zdrscl.f"
L10:
#line 139 "zdrscl.f"
    cden1 = cden * smlnum;
#line 140 "zdrscl.f"
    cnum1 = cnum / bignum;
#line 141 "zdrscl.f"
    if (abs(cden1) > abs(cnum) && cnum != 0.) {

/*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */

#line 145 "zdrscl.f"
	mul = smlnum;
#line 146 "zdrscl.f"
	done = FALSE_;
#line 147 "zdrscl.f"
	cden = cden1;
#line 148 "zdrscl.f"
    } else if (abs(cnum1) > abs(cden)) {

/*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

#line 152 "zdrscl.f"
	mul = bignum;
#line 153 "zdrscl.f"
	done = FALSE_;
#line 154 "zdrscl.f"
	cnum = cnum1;
#line 155 "zdrscl.f"
    } else {

/*        Multiply X by CNUM / CDEN and return. */

#line 159 "zdrscl.f"
	mul = cnum / cden;
#line 160 "zdrscl.f"
	done = TRUE_;
#line 161 "zdrscl.f"
    }

/*     Scale the vector X by MUL */

#line 165 "zdrscl.f"
    zdscal_(n, &mul, &sx[1], incx);

#line 167 "zdrscl.f"
    if (! done) {
#line 167 "zdrscl.f"
	goto L10;
#line 167 "zdrscl.f"
    }

#line 170 "zdrscl.f"
    return 0;

/*     End of ZDRSCL */

} /* zdrscl_ */

