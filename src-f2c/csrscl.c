#line 1 "csrscl.f"
/* csrscl.f -- translated by f2c (version 20100827).
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

#line 1 "csrscl.f"
/* > \brief \b CSRSCL multiplies a vector by the reciprocal of a real scalar. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSRSCL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csrscl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csrscl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csrscl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSRSCL( N, SA, SX, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, N */
/*       REAL               SA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            SX( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSRSCL multiplies an n-element complex vector x by the real scalar */
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
/* >          SA is REAL */
/* >          The scalar a which is used to divide each component of x. */
/* >          SA must be >= 0, or the subroutine will divide by zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX */
/* > \verbatim */
/* >          SX is COMPLEX array, dimension */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int csrscl_(integer *n, doublereal *sa, doublecomplex *sx, 
	integer *incx)
{
    static doublereal mul, cden;
    static logical done;
    static doublereal cnum, cden1, cnum1;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
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

#line 124 "csrscl.f"
    /* Parameter adjustments */
#line 124 "csrscl.f"
    --sx;
#line 124 "csrscl.f"

#line 124 "csrscl.f"
    /* Function Body */
#line 124 "csrscl.f"
    if (*n <= 0) {
#line 124 "csrscl.f"
	return 0;
#line 124 "csrscl.f"
    }

/*     Get machine parameters */

#line 129 "csrscl.f"
    smlnum = slamch_("S", (ftnlen)1);
#line 130 "csrscl.f"
    bignum = 1. / smlnum;
#line 131 "csrscl.f"
    slabad_(&smlnum, &bignum);

/*     Initialize the denominator to SA and the numerator to 1. */

#line 135 "csrscl.f"
    cden = *sa;
#line 136 "csrscl.f"
    cnum = 1.;

#line 138 "csrscl.f"
L10:
#line 139 "csrscl.f"
    cden1 = cden * smlnum;
#line 140 "csrscl.f"
    cnum1 = cnum / bignum;
#line 141 "csrscl.f"
    if (abs(cden1) > abs(cnum) && cnum != 0.) {

/*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */

#line 145 "csrscl.f"
	mul = smlnum;
#line 146 "csrscl.f"
	done = FALSE_;
#line 147 "csrscl.f"
	cden = cden1;
#line 148 "csrscl.f"
    } else if (abs(cnum1) > abs(cden)) {

/*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

#line 152 "csrscl.f"
	mul = bignum;
#line 153 "csrscl.f"
	done = FALSE_;
#line 154 "csrscl.f"
	cnum = cnum1;
#line 155 "csrscl.f"
    } else {

/*        Multiply X by CNUM / CDEN and return. */

#line 159 "csrscl.f"
	mul = cnum / cden;
#line 160 "csrscl.f"
	done = TRUE_;
#line 161 "csrscl.f"
    }

/*     Scale the vector X by MUL */

#line 165 "csrscl.f"
    csscal_(n, &mul, &sx[1], incx);

#line 167 "csrscl.f"
    if (! done) {
#line 167 "csrscl.f"
	goto L10;
#line 167 "csrscl.f"
    }

#line 170 "csrscl.f"
    return 0;

/*     End of CSRSCL */

} /* csrscl_ */

