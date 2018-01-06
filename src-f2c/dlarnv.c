#line 1 "dlarnv.f"
/* dlarnv.f -- translated by f2c (version 20100827).
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

#line 1 "dlarnv.f"
/* > \brief \b DLARNV returns a vector of random numbers from a uniform or normal distribution. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARNV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarnv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarnv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarnv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARNV( IDIST, ISEED, N, X ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IDIST, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISEED( 4 ) */
/*       DOUBLE PRECISION   X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARNV returns a vector of n random real numbers from a uniform or */
/* > normal distribution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IDIST */
/* > \verbatim */
/* >          IDIST is INTEGER */
/* >          Specifies the distribution of the random numbers: */
/* >          = 1:  uniform (0,1) */
/* >          = 2:  uniform (-1,1) */
/* >          = 3:  normal (0,1) */
/* > \endverbatim */
/* > */
/* > \param[in,out] ISEED */
/* > \verbatim */
/* >          ISEED is INTEGER array, dimension (4) */
/* >          On entry, the seed of the random number generator; the array */
/* >          elements must be between 0 and 4095, and ISEED(4) must be */
/* >          odd. */
/* >          On exit, the seed is updated. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of random numbers to be generated. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (N) */
/* >          The generated random numbers. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  This routine calls the auxiliary routine DLARUV to generate random */
/* >  real numbers from a uniform (0,1) distribution, in batches of up to */
/* >  128 using vectorisable code. The Box-Muller method is used to */
/* >  transform numbers from a uniform to a normal distribution. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlarnv_(integer *idist, integer *iseed, integer *n, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal u[128];
    static integer il, iv, il2;
    extern /* Subroutine */ int dlaruv_(integer *, integer *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 137 "dlarnv.f"
    /* Parameter adjustments */
#line 137 "dlarnv.f"
    --x;
#line 137 "dlarnv.f"
    --iseed;
#line 137 "dlarnv.f"

#line 137 "dlarnv.f"
    /* Function Body */
#line 137 "dlarnv.f"
    i__1 = *n;
#line 137 "dlarnv.f"
    for (iv = 1; iv <= i__1; iv += 64) {
/* Computing MIN */
#line 138 "dlarnv.f"
	i__2 = 64, i__3 = *n - iv + 1;
#line 138 "dlarnv.f"
	il = min(i__2,i__3);
#line 139 "dlarnv.f"
	if (*idist == 3) {
#line 140 "dlarnv.f"
	    il2 = il << 1;
#line 141 "dlarnv.f"
	} else {
#line 142 "dlarnv.f"
	    il2 = il;
#line 143 "dlarnv.f"
	}

/*        Call DLARUV to generate IL2 numbers from a uniform (0,1) */
/*        distribution (IL2 <= LV) */

#line 148 "dlarnv.f"
	dlaruv_(&iseed[1], &il2, u);

#line 150 "dlarnv.f"
	if (*idist == 1) {

/*           Copy generated numbers */

#line 154 "dlarnv.f"
	    i__2 = il;
#line 154 "dlarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 155 "dlarnv.f"
		x[iv + i__ - 1] = u[i__ - 1];
#line 156 "dlarnv.f"
/* L10: */
#line 156 "dlarnv.f"
	    }
#line 157 "dlarnv.f"
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribution */

#line 161 "dlarnv.f"
	    i__2 = il;
#line 161 "dlarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 162 "dlarnv.f"
		x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
#line 163 "dlarnv.f"
/* L20: */
#line 163 "dlarnv.f"
	    }
#line 164 "dlarnv.f"
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distribution */

#line 168 "dlarnv.f"
	    i__2 = il;
#line 168 "dlarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "dlarnv.f"
		x[iv + i__ - 1] = sqrt(log(u[(i__ << 1) - 2]) * -2.) * cos(u[(
			i__ << 1) - 1] * 6.2831853071795864769252867663);
#line 171 "dlarnv.f"
/* L30: */
#line 171 "dlarnv.f"
	    }
#line 172 "dlarnv.f"
	}
#line 173 "dlarnv.f"
/* L40: */
#line 173 "dlarnv.f"
    }
#line 174 "dlarnv.f"
    return 0;

/*     End of DLARNV */

} /* dlarnv_ */

