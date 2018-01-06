#line 1 "clarnv.f"
/* clarnv.f -- translated by f2c (version 20100827).
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

#line 1 "clarnv.f"
/* > \brief \b CLARNV returns a vector of random numbers from a uniform or normal distribution. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARNV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarnv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarnv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarnv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARNV( IDIST, ISEED, N, X ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IDIST, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISEED( 4 ) */
/*       COMPLEX            X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARNV returns a vector of n random complex numbers from a uniform or */
/* > normal distribution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IDIST */
/* > \verbatim */
/* >          IDIST is INTEGER */
/* >          Specifies the distribution of the random numbers: */
/* >          = 1:  real and imaginary parts each uniform (0,1) */
/* >          = 2:  real and imaginary parts each uniform (-1,1) */
/* >          = 3:  real and imaginary parts each normal (0,1) */
/* >          = 4:  uniformly distributed on the disc abs(z) < 1 */
/* >          = 5:  uniformly distributed on the circle abs(z) = 1 */
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
/* >          X is COMPLEX array, dimension (N) */
/* >          The generated random numbers. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  This routine calls the auxiliary routine SLARUV to generate random */
/* >  real numbers from a uniform (0,1) distribution, in batches of up to */
/* >  128 using vectorisable code. The Box-Muller method is used to */
/* >  transform numbers from a uniform to a normal distribution. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clarnv_(integer *idist, integer *iseed, integer *n, 
	doublecomplex *x)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);
    void z_exp(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublereal u[128];
    static integer il, iv;
    extern /* Subroutine */ int slaruv_(integer *, integer *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 139 "clarnv.f"
    /* Parameter adjustments */
#line 139 "clarnv.f"
    --x;
#line 139 "clarnv.f"
    --iseed;
#line 139 "clarnv.f"

#line 139 "clarnv.f"
    /* Function Body */
#line 139 "clarnv.f"
    i__1 = *n;
#line 139 "clarnv.f"
    for (iv = 1; iv <= i__1; iv += 64) {
/* Computing MIN */
#line 140 "clarnv.f"
	i__2 = 64, i__3 = *n - iv + 1;
#line 140 "clarnv.f"
	il = min(i__2,i__3);

/*        Call SLARUV to generate 2*IL real numbers from a uniform (0,1) */
/*        distribution (2*IL <= LV) */

#line 145 "clarnv.f"
	i__2 = il << 1;
#line 145 "clarnv.f"
	slaruv_(&iseed[1], &i__2, u);

#line 147 "clarnv.f"
	if (*idist == 1) {

/*           Copy generated numbers */

#line 151 "clarnv.f"
	    i__2 = il;
#line 151 "clarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 152 "clarnv.f"
		i__3 = iv + i__ - 1;
#line 152 "clarnv.f"
		i__4 = (i__ << 1) - 2;
#line 152 "clarnv.f"
		i__5 = (i__ << 1) - 1;
#line 152 "clarnv.f"
		z__1.r = u[i__4], z__1.i = u[i__5];
#line 152 "clarnv.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 153 "clarnv.f"
/* L10: */
#line 153 "clarnv.f"
	    }
#line 154 "clarnv.f"
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribution */

#line 158 "clarnv.f"
	    i__2 = il;
#line 158 "clarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 159 "clarnv.f"
		i__3 = iv + i__ - 1;
#line 159 "clarnv.f"
		d__1 = u[(i__ << 1) - 2] * 2. - 1.;
#line 159 "clarnv.f"
		d__2 = u[(i__ << 1) - 1] * 2. - 1.;
#line 159 "clarnv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 159 "clarnv.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 161 "clarnv.f"
/* L20: */
#line 161 "clarnv.f"
	    }
#line 162 "clarnv.f"
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distribution */

#line 166 "clarnv.f"
	    i__2 = il;
#line 166 "clarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 167 "clarnv.f"
		i__3 = iv + i__ - 1;
#line 167 "clarnv.f"
		d__1 = sqrt(log(u[(i__ << 1) - 2]) * -2.);
#line 167 "clarnv.f"
		d__2 = u[(i__ << 1) - 1] * 6.2831853071795864769252867663;
#line 167 "clarnv.f"
		z__3.r = 0., z__3.i = d__2;
#line 167 "clarnv.f"
		z_exp(&z__2, &z__3);
#line 167 "clarnv.f"
		z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
#line 167 "clarnv.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 169 "clarnv.f"
/* L30: */
#line 169 "clarnv.f"
	    }
#line 170 "clarnv.f"
	} else if (*idist == 4) {

/*           Convert generated numbers to complex numbers uniformly */
/*           distributed on the unit disk */

#line 175 "clarnv.f"
	    i__2 = il;
#line 175 "clarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 176 "clarnv.f"
		i__3 = iv + i__ - 1;
#line 176 "clarnv.f"
		d__1 = sqrt(u[(i__ << 1) - 2]);
#line 176 "clarnv.f"
		d__2 = u[(i__ << 1) - 1] * 6.2831853071795864769252867663;
#line 176 "clarnv.f"
		z__3.r = 0., z__3.i = d__2;
#line 176 "clarnv.f"
		z_exp(&z__2, &z__3);
#line 176 "clarnv.f"
		z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
#line 176 "clarnv.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 178 "clarnv.f"
/* L40: */
#line 178 "clarnv.f"
	    }
#line 179 "clarnv.f"
	} else if (*idist == 5) {

/*           Convert generated numbers to complex numbers uniformly */
/*           distributed on the unit circle */

#line 184 "clarnv.f"
	    i__2 = il;
#line 184 "clarnv.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 185 "clarnv.f"
		i__3 = iv + i__ - 1;
#line 185 "clarnv.f"
		d__1 = u[(i__ << 1) - 1] * 6.2831853071795864769252867663;
#line 185 "clarnv.f"
		z__2.r = 0., z__2.i = d__1;
#line 185 "clarnv.f"
		z_exp(&z__1, &z__2);
#line 185 "clarnv.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 186 "clarnv.f"
/* L50: */
#line 186 "clarnv.f"
	    }
#line 187 "clarnv.f"
	}
#line 188 "clarnv.f"
/* L60: */
#line 188 "clarnv.f"
    }
#line 189 "clarnv.f"
    return 0;

/*     End of CLARNV */

} /* clarnv_ */

