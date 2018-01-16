#line 1 "../INSTALL/dlamchf77.f"
/* ../INSTALL/dlamchf77.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/dlamchf77.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b32 = 0.;

/* > \brief \b DLAMCHF77 deprecated */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*      DOUBLE PRECISION FUNCTION DLAMCH( CMACH ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAMCHF77 determines double precision machine parameters. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] CMACH */
/* > \verbatim */
/* >          Specifies the value to be returned by DLAMCH: */
/* >          = 'E' or 'e',   DLAMCH := eps */
/* >          = 'S' or 's ,   DLAMCH := sfmin */
/* >          = 'B' or 'b',   DLAMCH := base */
/* >          = 'P' or 'p',   DLAMCH := eps*base */
/* >          = 'N' or 'n',   DLAMCH := t */
/* >          = 'R' or 'r',   DLAMCH := rnd */
/* >          = 'M' or 'm',   DLAMCH := emin */
/* >          = 'U' or 'u',   DLAMCH := rmin */
/* >          = 'L' or 'l',   DLAMCH := emax */
/* >          = 'O' or 'o',   DLAMCH := rmax */
/* >          where */
/* >          eps   = relative machine precision */
/* >          sfmin = safe minimum, such that 1/sfmin does not overflow */
/* >          base  = base of the machine */
/* >          prec  = eps*base */
/* >          t     = number of (base) digits in the mantissa */
/* >          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/* >          emin  = minimum exponent before (gradual) underflow */
/* >          rmin  = underflow threshold - base**(emin-1) */
/* >          emax  = largest exponent before overflow */
/* >          rmax  = overflow threshold  - (base**emax)*(1-eps) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
doublereal dlamch_(char *cmach, ftnlen cmach_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal t;
    static integer it;
    static doublereal rnd, eps, base;
    static integer beta;
    static doublereal emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static doublereal rmin, rmax, rmach;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small, sfmin;
    extern /* Subroutine */ int dlamc2_(integer *, integer *, logical *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

#line 100 "../INSTALL/dlamchf77.f"
    if (first) {
#line 101 "../INSTALL/dlamchf77.f"
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
#line 102 "../INSTALL/dlamchf77.f"
	base = (doublereal) beta;
#line 103 "../INSTALL/dlamchf77.f"
	t = (doublereal) it;
#line 104 "../INSTALL/dlamchf77.f"
	if (lrnd) {
#line 105 "../INSTALL/dlamchf77.f"
	    rnd = 1.;
#line 106 "../INSTALL/dlamchf77.f"
	    i__1 = 1 - it;
#line 106 "../INSTALL/dlamchf77.f"
	    eps = pow_di(&base, &i__1) / 2;
#line 107 "../INSTALL/dlamchf77.f"
	} else {
#line 108 "../INSTALL/dlamchf77.f"
	    rnd = 0.;
#line 109 "../INSTALL/dlamchf77.f"
	    i__1 = 1 - it;
#line 109 "../INSTALL/dlamchf77.f"
	    eps = pow_di(&base, &i__1);
#line 110 "../INSTALL/dlamchf77.f"
	}
#line 111 "../INSTALL/dlamchf77.f"
	prec = eps * base;
#line 112 "../INSTALL/dlamchf77.f"
	emin = (doublereal) imin;
#line 113 "../INSTALL/dlamchf77.f"
	emax = (doublereal) imax;
#line 114 "../INSTALL/dlamchf77.f"
	sfmin = rmin;
#line 115 "../INSTALL/dlamchf77.f"
	small = 1. / rmax;
#line 116 "../INSTALL/dlamchf77.f"
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rounding */
/*           causing overflow when computing  1/sfmin. */

#line 121 "../INSTALL/dlamchf77.f"
	    sfmin = small * (eps + 1.);
#line 122 "../INSTALL/dlamchf77.f"
	}
#line 123 "../INSTALL/dlamchf77.f"
    }

#line 125 "../INSTALL/dlamchf77.f"
    if (lsame_(cmach, "E", (ftnlen)1, (ftnlen)1)) {
#line 126 "../INSTALL/dlamchf77.f"
	rmach = eps;
#line 127 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "S", (ftnlen)1, (ftnlen)1)) {
#line 128 "../INSTALL/dlamchf77.f"
	rmach = sfmin;
#line 129 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "B", (ftnlen)1, (ftnlen)1)) {
#line 130 "../INSTALL/dlamchf77.f"
	rmach = base;
#line 131 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "P", (ftnlen)1, (ftnlen)1)) {
#line 132 "../INSTALL/dlamchf77.f"
	rmach = prec;
#line 133 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "N", (ftnlen)1, (ftnlen)1)) {
#line 134 "../INSTALL/dlamchf77.f"
	rmach = t;
#line 135 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "R", (ftnlen)1, (ftnlen)1)) {
#line 136 "../INSTALL/dlamchf77.f"
	rmach = rnd;
#line 137 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "M", (ftnlen)1, (ftnlen)1)) {
#line 138 "../INSTALL/dlamchf77.f"
	rmach = emin;
#line 139 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "U", (ftnlen)1, (ftnlen)1)) {
#line 140 "../INSTALL/dlamchf77.f"
	rmach = rmin;
#line 141 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "L", (ftnlen)1, (ftnlen)1)) {
#line 142 "../INSTALL/dlamchf77.f"
	rmach = emax;
#line 143 "../INSTALL/dlamchf77.f"
    } else if (lsame_(cmach, "O", (ftnlen)1, (ftnlen)1)) {
#line 144 "../INSTALL/dlamchf77.f"
	rmach = rmax;
#line 145 "../INSTALL/dlamchf77.f"
    }

#line 147 "../INSTALL/dlamchf77.f"
    ret_val = rmach;
#line 148 "../INSTALL/dlamchf77.f"
    first = FALSE_;
#line 149 "../INSTALL/dlamchf77.f"
    return ret_val;

/*     End of DLAMCH */

} /* dlamch_ */


/* *********************************************************************** */

/* > \brief \b DLAMC1 */
/* > \details */
/* > \b Purpose: */
/* > \verbatim */
/* > DLAMC1 determines the machine parameters given by BETA, T, RND, and */
/* > IEEE1. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          The base of the machine. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          The number of ( BETA ) digits in the mantissa. */
/* > \endverbatim */
/* > */
/* > \param[out] RND */
/* > \verbatim */
/* >          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/* >          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/* >          be a reliable guide to the way in which the machine performs */
/* >          its arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[out] IEEE1 */
/* > \verbatim */
/* >          Specifies whether rounding appears to be done in the IEEE */
/* >          'round to nearest' style. */
/* > \endverbatim */
/* > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. 
of Colorado Denver and NAG Ltd.. */
/* > \date April 2012 */
/* > \ingroup auxOTHERauxiliary */
/* > */
/* > \details \b Further \b Details */
/* > \verbatim */
/* > */
/* >  The routine is based on the routine  ENVRON  by Malcolm and */
/* >  incorporates suggestions by Gentleman and Marovich. See */
/* > */
/* >     Malcolm M. A. (1972) Algorithms to reveal properties of */
/* >        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */
/* > */
/* >     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
/* >        that reveal properties of floating point arithmetic units. */
/* >        Comms. of the ACM, 17, 276-277. */
/* > \endverbatim */
/* > */
/* Subroutine */ int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a, b, c__, f, t1, t2;
    static integer lt;
    static doublereal one, qtr;
    static logical lrnd;
    static integer lbeta;
    static doublereal savec;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static logical lieee1;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2010 */

/*     .. Scalar Arguments .. */
/*     .. */
/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

#line 235 "../INSTALL/dlamchf77.f"
    if (first) {
#line 236 "../INSTALL/dlamchf77.f"
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
/*        IEEE1, T and RND. */

/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are  stored and not held in registers,  or */
/*        are not affected by optimizers. */

/*        Compute  a = 2.0**m  with the  smallest positive integer m such */
/*        that */

/*           fl( a + 1.0 ) = a. */

#line 250 "../INSTALL/dlamchf77.f"
	a = 1.;
#line 251 "../INSTALL/dlamchf77.f"
	c__ = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
#line 254 "../INSTALL/dlamchf77.f"
L10:
#line 255 "../INSTALL/dlamchf77.f"
	if (c__ == one) {
#line 256 "../INSTALL/dlamchf77.f"
	    a *= 2;
#line 257 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&a, &one);
#line 258 "../INSTALL/dlamchf77.f"
	    d__1 = -a;
#line 258 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&c__, &d__1);
#line 259 "../INSTALL/dlamchf77.f"
	    goto L10;
#line 260 "../INSTALL/dlamchf77.f"
	}
/* +       END WHILE */

/*        Now compute  b = 2.0**m  with the smallest positive integer m */
/*        such that */

/*           fl( a + b ) .gt. a. */

#line 268 "../INSTALL/dlamchf77.f"
	b = 1.;
#line 269 "../INSTALL/dlamchf77.f"
	c__ = dlamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
#line 272 "../INSTALL/dlamchf77.f"
L20:
#line 273 "../INSTALL/dlamchf77.f"
	if (c__ == a) {
#line 274 "../INSTALL/dlamchf77.f"
	    b *= 2;
#line 275 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&a, &b);
#line 276 "../INSTALL/dlamchf77.f"
	    goto L20;
#line 277 "../INSTALL/dlamchf77.f"
	}
/* +       END WHILE */

/*        Now compute the base.  a and c  are neighbouring floating point */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
/*        their difference is beta. Adding 0.25 to c is to ensure that it */
/*        is truncated to beta and not ( beta - 1 ). */

#line 285 "../INSTALL/dlamchf77.f"
	qtr = one / 4;
#line 286 "../INSTALL/dlamchf77.f"
	savec = c__;
#line 287 "../INSTALL/dlamchf77.f"
	d__1 = -a;
#line 287 "../INSTALL/dlamchf77.f"
	c__ = dlamc3_(&c__, &d__1);
#line 288 "../INSTALL/dlamchf77.f"
	lbeta = (integer) (c__ + qtr);

/*        Now determine whether rounding or chopping occurs,  by adding a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

#line 293 "../INSTALL/dlamchf77.f"
	b = (doublereal) lbeta;
#line 294 "../INSTALL/dlamchf77.f"
	d__1 = b / 2;
#line 294 "../INSTALL/dlamchf77.f"
	d__2 = -b / 100;
#line 294 "../INSTALL/dlamchf77.f"
	f = dlamc3_(&d__1, &d__2);
#line 295 "../INSTALL/dlamchf77.f"
	c__ = dlamc3_(&f, &a);
#line 296 "../INSTALL/dlamchf77.f"
	if (c__ == a) {
#line 297 "../INSTALL/dlamchf77.f"
	    lrnd = TRUE_;
#line 298 "../INSTALL/dlamchf77.f"
	} else {
#line 299 "../INSTALL/dlamchf77.f"
	    lrnd = FALSE_;
#line 300 "../INSTALL/dlamchf77.f"
	}
#line 301 "../INSTALL/dlamchf77.f"
	d__1 = b / 2;
#line 301 "../INSTALL/dlamchf77.f"
	d__2 = b / 100;
#line 301 "../INSTALL/dlamchf77.f"
	f = dlamc3_(&d__1, &d__2);
#line 302 "../INSTALL/dlamchf77.f"
	c__ = dlamc3_(&f, &a);
#line 303 "../INSTALL/dlamchf77.f"
	if (lrnd && c__ == a) {
#line 303 "../INSTALL/dlamchf77.f"
	    lrnd = FALSE_;
#line 303 "../INSTALL/dlamchf77.f"
	}

/*        Try and decide whether rounding is done in the  IEEE  'round to */
/*        nearest' style. B/2 is half a unit in the last place of the two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
/*        A, but adding B/2 to SAVEC should change SAVEC. */

#line 312 "../INSTALL/dlamchf77.f"
	d__1 = b / 2;
#line 312 "../INSTALL/dlamchf77.f"
	t1 = dlamc3_(&d__1, &a);
#line 313 "../INSTALL/dlamchf77.f"
	d__1 = b / 2;
#line 313 "../INSTALL/dlamchf77.f"
	t2 = dlamc3_(&d__1, &savec);
#line 314 "../INSTALL/dlamchf77.f"
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part of */
/*        log to the base beta of a,  however it is safer to determine  t */
/*        by powering.  So we find t as the smallest positive integer for */
/*        which */

/*           fl( beta**t + 1.0 ) = 1.0. */

#line 323 "../INSTALL/dlamchf77.f"
	lt = 0;
#line 324 "../INSTALL/dlamchf77.f"
	a = 1.;
#line 325 "../INSTALL/dlamchf77.f"
	c__ = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
#line 328 "../INSTALL/dlamchf77.f"
L30:
#line 329 "../INSTALL/dlamchf77.f"
	if (c__ == one) {
#line 330 "../INSTALL/dlamchf77.f"
	    ++lt;
#line 331 "../INSTALL/dlamchf77.f"
	    a *= lbeta;
#line 332 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&a, &one);
#line 333 "../INSTALL/dlamchf77.f"
	    d__1 = -a;
#line 333 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&c__, &d__1);
#line 334 "../INSTALL/dlamchf77.f"
	    goto L30;
#line 335 "../INSTALL/dlamchf77.f"
	}
/* +       END WHILE */

#line 338 "../INSTALL/dlamchf77.f"
    }

#line 340 "../INSTALL/dlamchf77.f"
    *beta = lbeta;
#line 341 "../INSTALL/dlamchf77.f"
    *t = lt;
#line 342 "../INSTALL/dlamchf77.f"
    *rnd = lrnd;
#line 343 "../INSTALL/dlamchf77.f"
    *ieee1 = lieee1;
#line 344 "../INSTALL/dlamchf77.f"
    first = FALSE_;
#line 345 "../INSTALL/dlamchf77.f"
    return 0;

/*     End of DLAMC1 */

} /* dlamc1_ */


/* *********************************************************************** */

/* > \brief \b DLAMC2 */
/* > \details */
/* > \b Purpose: */
/* > \verbatim */
/* > DLAMC2 determines the machine parameters specified in its argument */
/* > list. */
/* > \endverbatim */
/* > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. 
of Colorado Denver and NAG Ltd.. */
/* > \date April 2012 */
/* > \ingroup auxOTHERauxiliary */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          The base of the machine. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          The number of ( BETA ) digits in the mantissa. */
/* > \endverbatim */
/* > */
/* > \param[out] RND */
/* > \verbatim */
/* >          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/* >          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/* >          be a reliable guide to the way in which the machine performs */
/* >          its arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[out] EPS */
/* > \verbatim */
/* >          The smallest positive number such that */
/* >             fl( 1.0 - EPS ) .LT. 1.0, */
/* >          where fl denotes the computed value. */
/* > \endverbatim */
/* > */
/* > \param[out] EMIN */
/* > \verbatim */
/* >          The minimum exponent before (gradual) underflow occurs. */
/* > \endverbatim */
/* > */
/* > \param[out] RMIN */
/* > \verbatim */
/* >          The smallest normalized number for the machine, given by */
/* >          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
/* >          of BETA. */
/* > \endverbatim */
/* > */
/* > \param[out] EMAX */
/* > \verbatim */
/* >          The maximum exponent before overflow occurs. */
/* > \endverbatim */
/* > */
/* > \param[out] RMAX */
/* > \verbatim */
/* >          The largest positive number for the machine, given by */
/* >          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
/* >          value of BETA. */
/* > \endverbatim */
/* > */
/* > \details \b Further \b Details */
/* > \verbatim */
/* > */
/* >  The computation of  EPS  is based on a routine PARANOIA by */
/* >  W. Kahan of the University of California at Berkeley. */
/* > \endverbatim */
/* Subroutine */ int dlamc2_(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax)
{
    /* Initialized data */

    static logical first = TRUE_;
    static logical iwarn = FALSE_;

    /* Format strings */
    static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre"
	    "ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the va"
	    "lue EMIN looks\002,\002 acceptable please comment out \002,/\002"
	    " the IF block as marked within the code of routine\002,\002 DLAM"
	    "C2,\002,/\002 otherwise supply EMIN explicitly.\002,/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal a, b, c__;
    static integer i__, lt;
    static doublereal one, two;
    static logical ieee;
    static doublereal half;
    static logical lrnd;
    static doublereal leps, zero;
    static integer lbeta;
    static doublereal rbase;
    static integer lemin, lemax, gnmin;
    static doublereal small;
    static integer gpmin;
    static doublereal third, lrmin, lrmax, sixth;
    extern /* Subroutine */ int dlamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static logical lieee1;
    extern /* Subroutine */ int dlamc4_(integer *, doublereal *, integer *), 
	    dlamc5_(integer *, integer *, integer *, logical *, integer *, 
	    doublereal *);
    static integer ngnmin, ngpmin;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2010 */

/*     .. Scalar Arguments .. */
/*     .. */
/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

#line 458 "../INSTALL/dlamchf77.f"
    if (first) {
#line 459 "../INSTALL/dlamchf77.f"
	zero = 0.;
#line 460 "../INSTALL/dlamchf77.f"
	one = 1.;
#line 461 "../INSTALL/dlamchf77.f"
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */

/*        Throughout this routine  we use the function  DLAMC3  to ensure */
/*        that relevant values are stored  and not held in registers,  or */
/*        are not affected by optimizers. */

/*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

#line 472 "../INSTALL/dlamchf77.f"
	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

#line 476 "../INSTALL/dlamchf77.f"
	b = (doublereal) lbeta;
#line 477 "../INSTALL/dlamchf77.f"
	i__1 = -lt;
#line 477 "../INSTALL/dlamchf77.f"
	a = pow_di(&b, &i__1);
#line 478 "../INSTALL/dlamchf77.f"
	leps = a;

/*        Try some tricks to see whether or not this is the correct  EPS. */

#line 482 "../INSTALL/dlamchf77.f"
	b = two / 3;
#line 483 "../INSTALL/dlamchf77.f"
	half = one / 2;
#line 484 "../INSTALL/dlamchf77.f"
	d__1 = -half;
#line 484 "../INSTALL/dlamchf77.f"
	sixth = dlamc3_(&b, &d__1);
#line 485 "../INSTALL/dlamchf77.f"
	third = dlamc3_(&sixth, &sixth);
#line 486 "../INSTALL/dlamchf77.f"
	d__1 = -half;
#line 486 "../INSTALL/dlamchf77.f"
	b = dlamc3_(&third, &d__1);
#line 487 "../INSTALL/dlamchf77.f"
	b = dlamc3_(&b, &sixth);
#line 488 "../INSTALL/dlamchf77.f"
	b = abs(b);
#line 489 "../INSTALL/dlamchf77.f"
	if (b < leps) {
#line 489 "../INSTALL/dlamchf77.f"
	    b = leps;
#line 489 "../INSTALL/dlamchf77.f"
	}

#line 492 "../INSTALL/dlamchf77.f"
	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
#line 495 "../INSTALL/dlamchf77.f"
L10:
#line 496 "../INSTALL/dlamchf77.f"
	if (leps > b && b > zero) {
#line 497 "../INSTALL/dlamchf77.f"
	    leps = b;
#line 498 "../INSTALL/dlamchf77.f"
	    d__1 = half * leps;
/* Computing 5th power */
#line 498 "../INSTALL/dlamchf77.f"
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
#line 498 "../INSTALL/dlamchf77.f"
	    d__5 = leps;
#line 498 "../INSTALL/dlamchf77.f"
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
#line 498 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&d__1, &d__2);
#line 499 "../INSTALL/dlamchf77.f"
	    d__1 = -c__;
#line 499 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&half, &d__1);
#line 500 "../INSTALL/dlamchf77.f"
	    b = dlamc3_(&half, &c__);
#line 501 "../INSTALL/dlamchf77.f"
	    d__1 = -b;
#line 501 "../INSTALL/dlamchf77.f"
	    c__ = dlamc3_(&half, &d__1);
#line 502 "../INSTALL/dlamchf77.f"
	    b = dlamc3_(&half, &c__);
#line 503 "../INSTALL/dlamchf77.f"
	    goto L10;
#line 504 "../INSTALL/dlamchf77.f"
	}
/* +       END WHILE */

#line 507 "../INSTALL/dlamchf77.f"
	if (a < leps) {
#line 507 "../INSTALL/dlamchf77.f"
	    leps = a;
#line 507 "../INSTALL/dlamchf77.f"
	}

/*        Computation of EPS complete. */

/*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). */
/*        Keep dividing  A by BETA until (gradual) underflow occurs. This */
/*        is detected when we cannot recover the previous A. */

#line 516 "../INSTALL/dlamchf77.f"
	rbase = one / lbeta;
#line 517 "../INSTALL/dlamchf77.f"
	small = one;
#line 518 "../INSTALL/dlamchf77.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 519 "../INSTALL/dlamchf77.f"
	    d__1 = small * rbase;
#line 519 "../INSTALL/dlamchf77.f"
	    small = dlamc3_(&d__1, &zero);
#line 520 "../INSTALL/dlamchf77.f"
/* L20: */
#line 520 "../INSTALL/dlamchf77.f"
	}
#line 521 "../INSTALL/dlamchf77.f"
	a = dlamc3_(&one, &small);
#line 522 "../INSTALL/dlamchf77.f"
	dlamc4_(&ngpmin, &one, &lbeta);
#line 523 "../INSTALL/dlamchf77.f"
	d__1 = -one;
#line 523 "../INSTALL/dlamchf77.f"
	dlamc4_(&ngnmin, &d__1, &lbeta);
#line 524 "../INSTALL/dlamchf77.f"
	dlamc4_(&gpmin, &a, &lbeta);
#line 525 "../INSTALL/dlamchf77.f"
	d__1 = -a;
#line 525 "../INSTALL/dlamchf77.f"
	dlamc4_(&gnmin, &d__1, &lbeta);
#line 526 "../INSTALL/dlamchf77.f"
	ieee = FALSE_;

#line 528 "../INSTALL/dlamchf77.f"
	if (ngpmin == ngnmin && gpmin == gnmin) {
#line 529 "../INSTALL/dlamchf77.f"
	    if (ngpmin == gpmin) {
#line 530 "../INSTALL/dlamchf77.f"
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual underflow; */
/*              e.g.,  VAX ) */
#line 533 "../INSTALL/dlamchf77.f"
	    } else if (gpmin - ngpmin == 3) {
#line 534 "../INSTALL/dlamchf77.f"
		lemin = ngpmin - 1 + lt;
#line 535 "../INSTALL/dlamchf77.f"
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual underflow; */
/*              e.g., IEEE standard followers ) */
#line 538 "../INSTALL/dlamchf77.f"
	    } else {
#line 539 "../INSTALL/dlamchf77.f"
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
#line 541 "../INSTALL/dlamchf77.f"
		iwarn = TRUE_;
#line 542 "../INSTALL/dlamchf77.f"
	    }

#line 544 "../INSTALL/dlamchf77.f"
	} else if (ngpmin == gpmin && ngnmin == gnmin) {
#line 545 "../INSTALL/dlamchf77.f"
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
#line 546 "../INSTALL/dlamchf77.f"
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow; */
/*              e.g., CYBER 205 ) */
#line 549 "../INSTALL/dlamchf77.f"
	    } else {
#line 550 "../INSTALL/dlamchf77.f"
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
#line 552 "../INSTALL/dlamchf77.f"
		iwarn = TRUE_;
#line 553 "../INSTALL/dlamchf77.f"
	    }

#line 555 "../INSTALL/dlamchf77.f"
	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
#line 557 "../INSTALL/dlamchf77.f"
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
#line 558 "../INSTALL/dlamchf77.f"
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflow; */
/*              no known machine ) */
#line 561 "../INSTALL/dlamchf77.f"
	    } else {
#line 562 "../INSTALL/dlamchf77.f"
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
#line 564 "../INSTALL/dlamchf77.f"
		iwarn = TRUE_;
#line 565 "../INSTALL/dlamchf77.f"
	    }

#line 567 "../INSTALL/dlamchf77.f"
	} else {
/* Computing MIN */
#line 568 "../INSTALL/dlamchf77.f"
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
#line 568 "../INSTALL/dlamchf77.f"
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
#line 570 "../INSTALL/dlamchf77.f"
	    iwarn = TRUE_;
#line 571 "../INSTALL/dlamchf77.f"
	}
#line 572 "../INSTALL/dlamchf77.f"
	first = FALSE_;
/* ** */
/* Comment out this if block if EMIN is ok */
#line 575 "../INSTALL/dlamchf77.f"
	if (iwarn) {
#line 576 "../INSTALL/dlamchf77.f"
	    first = TRUE_;
#line 577 "../INSTALL/dlamchf77.f"
	    s_wsfe(&io___58);
#line 577 "../INSTALL/dlamchf77.f"
	    do_fio(&c__1, (char *)&lemin, (ftnlen)sizeof(integer));
#line 577 "../INSTALL/dlamchf77.f"
	    e_wsfe();
#line 578 "../INSTALL/dlamchf77.f"
	}
/* ** */

/*        Assume IEEE arithmetic if we found denormalised  numbers above, */
/*        or if arithmetic seems to round in the  IEEE style,  determined */
/*        in routine DLAMC1. A true IEEE machine should have both  things */
/*        true; however, faulty machines may have one or the other. */

#line 586 "../INSTALL/dlamchf77.f"
	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could compute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
/*        this computation. */

#line 592 "../INSTALL/dlamchf77.f"
	lrmin = 1.;
#line 593 "../INSTALL/dlamchf77.f"
	i__1 = 1 - lemin;
#line 593 "../INSTALL/dlamchf77.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 594 "../INSTALL/dlamchf77.f"
	    d__1 = lrmin * rbase;
#line 594 "../INSTALL/dlamchf77.f"
	    lrmin = dlamc3_(&d__1, &zero);
#line 595 "../INSTALL/dlamchf77.f"
/* L30: */
#line 595 "../INSTALL/dlamchf77.f"
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

#line 599 "../INSTALL/dlamchf77.f"
	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
#line 600 "../INSTALL/dlamchf77.f"
    }

#line 602 "../INSTALL/dlamchf77.f"
    *beta = lbeta;
#line 603 "../INSTALL/dlamchf77.f"
    *t = lt;
#line 604 "../INSTALL/dlamchf77.f"
    *rnd = lrnd;
#line 605 "../INSTALL/dlamchf77.f"
    *eps = leps;
#line 606 "../INSTALL/dlamchf77.f"
    *emin = lemin;
#line 607 "../INSTALL/dlamchf77.f"
    *rmin = lrmin;
#line 608 "../INSTALL/dlamchf77.f"
    *emax = lemax;
#line 609 "../INSTALL/dlamchf77.f"
    *rmax = lrmax;

#line 611 "../INSTALL/dlamchf77.f"
    return 0;


/*     End of DLAMC2 */

} /* dlamc2_ */


/* *********************************************************************** */

/* > \brief \b DLAMC3 */
/* > \details */
/* > \b Purpose: */
/* > \verbatim */
/* > DLAMC3  is intended to force  A  and  B  to be stored prior to doing */
/* > the addition of  A  and  B ,  for use in situations where optimizers */
/* > might hold one of these in a register. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          The values A and B. */
/* > \endverbatim */
doublereal dlamc3_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2010 */

/*     .. Scalar Arguments .. */
/*     .. */
/* ===================================================================== */

/*     .. Executable Statements .. */

#line 655 "../INSTALL/dlamchf77.f"
    ret_val = *a + *b;

#line 657 "../INSTALL/dlamchf77.f"
    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */


/* *********************************************************************** */

/* > \brief \b DLAMC4 */
/* > \details */
/* > \b Purpose: */
/* > \verbatim */
/* > DLAMC4 is a service routine for DLAMC2. */
/* > \endverbatim */
/* > */
/* > \param[out] EMIN */
/* > \verbatim */
/* >          The minimum exponent before (gradual) underflow, computed by */
/* >          setting A = START and dividing by BASE until the previous A */
/* >          can not be recovered. */
/* > \endverbatim */
/* > */
/* > \param[in] START */
/* > \verbatim */
/* >          The starting point for determining EMIN. */
/* > \endverbatim */
/* > */
/* > \param[in] BASE */
/* > \verbatim */
/* >          The base of the machine. */
/* > \endverbatim */
/* > */
/* Subroutine */ int dlamc4_(integer *emin, doublereal *start, integer *base)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal a;
    static integer i__;
    static doublereal b1, b2, c1, c2, d1, d2, one, zero, rbase;
    extern doublereal dlamc3_(doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2010 */

/*     .. Scalar Arguments .. */
/*     .. */
/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 711 "../INSTALL/dlamchf77.f"
    a = *start;
#line 712 "../INSTALL/dlamchf77.f"
    one = 1.;
#line 713 "../INSTALL/dlamchf77.f"
    rbase = one / *base;
#line 714 "../INSTALL/dlamchf77.f"
    zero = 0.;
#line 715 "../INSTALL/dlamchf77.f"
    *emin = 1;
#line 716 "../INSTALL/dlamchf77.f"
    d__1 = a * rbase;
#line 716 "../INSTALL/dlamchf77.f"
    b1 = dlamc3_(&d__1, &zero);
#line 717 "../INSTALL/dlamchf77.f"
    c1 = a;
#line 718 "../INSTALL/dlamchf77.f"
    c2 = a;
#line 719 "../INSTALL/dlamchf77.f"
    d1 = a;
#line 720 "../INSTALL/dlamchf77.f"
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
/*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
#line 723 "../INSTALL/dlamchf77.f"
L10:
#line 724 "../INSTALL/dlamchf77.f"
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
#line 726 "../INSTALL/dlamchf77.f"
	--(*emin);
#line 727 "../INSTALL/dlamchf77.f"
	a = b1;
#line 728 "../INSTALL/dlamchf77.f"
	d__1 = a / *base;
#line 728 "../INSTALL/dlamchf77.f"
	b1 = dlamc3_(&d__1, &zero);
#line 729 "../INSTALL/dlamchf77.f"
	d__1 = b1 * *base;
#line 729 "../INSTALL/dlamchf77.f"
	c1 = dlamc3_(&d__1, &zero);
#line 730 "../INSTALL/dlamchf77.f"
	d1 = zero;
#line 731 "../INSTALL/dlamchf77.f"
	i__1 = *base;
#line 731 "../INSTALL/dlamchf77.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 732 "../INSTALL/dlamchf77.f"
	    d1 += b1;
#line 733 "../INSTALL/dlamchf77.f"
/* L20: */
#line 733 "../INSTALL/dlamchf77.f"
	}
#line 734 "../INSTALL/dlamchf77.f"
	d__1 = a * rbase;
#line 734 "../INSTALL/dlamchf77.f"
	b2 = dlamc3_(&d__1, &zero);
#line 735 "../INSTALL/dlamchf77.f"
	d__1 = b2 / rbase;
#line 735 "../INSTALL/dlamchf77.f"
	c2 = dlamc3_(&d__1, &zero);
#line 736 "../INSTALL/dlamchf77.f"
	d2 = zero;
#line 737 "../INSTALL/dlamchf77.f"
	i__1 = *base;
#line 737 "../INSTALL/dlamchf77.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 738 "../INSTALL/dlamchf77.f"
	    d2 += b2;
#line 739 "../INSTALL/dlamchf77.f"
/* L30: */
#line 739 "../INSTALL/dlamchf77.f"
	}
#line 740 "../INSTALL/dlamchf77.f"
	goto L10;
#line 741 "../INSTALL/dlamchf77.f"
    }
/* +    END WHILE */

#line 744 "../INSTALL/dlamchf77.f"
    return 0;

/*     End of DLAMC4 */

} /* dlamc4_ */


/* *********************************************************************** */

/* > \brief \b DLAMC5 */
/* > \details */
/* > \b Purpose: */
/* > \verbatim */
/* > DLAMC5 attempts to compute RMAX, the largest machine floating-point */
/* > number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
/* > approximately to a power of 2.  It will fail on machines where this */
/* > assumption does not hold, for example, the Cyber 205 (EMIN = -28625, */
/* > EMAX = 28718).  It will also fail if the value supplied for EMIN is */
/* > too large (i.e. too close to zero), probably with overflow. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          The base of floating-point arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          The number of base BETA digits in the mantissa of a */
/* >          floating-point value. */
/* > \endverbatim */
/* > */
/* > \param[in] EMIN */
/* > \verbatim */
/* >          The minimum exponent before (gradual) underflow. */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* >          A logical flag specifying whether or not the arithmetic */
/* >          system is thought to comply with the IEEE standard. */
/* > \endverbatim */
/* > */
/* > \param[out] EMAX */
/* > \verbatim */
/* >          The largest exponent before overflow */
/* > \endverbatim */
/* > */
/* > \param[out] RMAX */
/* > \verbatim */
/* >          The largest machine floating-point number. */
/* > \endverbatim */
/* > */
/* Subroutine */ int dlamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal y, z__;
    static integer try__, lexp;
    static doublereal oldy;
    static integer uexp, nbits;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static doublereal recbas;
    static integer exbits, expsum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2010 */

/*     .. Scalar Arguments .. */
/*     .. */
/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     First compute LEXP and UEXP, two powers of 2 that bound */
/*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
/*     approximately to the bound that is closest to abs(EMIN). */
/*     (EMAX is the exponent of the required number RMAX). */

#line 831 "../INSTALL/dlamchf77.f"
    lexp = 1;
#line 832 "../INSTALL/dlamchf77.f"
    exbits = 1;
#line 833 "../INSTALL/dlamchf77.f"
L10:
#line 834 "../INSTALL/dlamchf77.f"
    try__ = lexp << 1;
#line 835 "../INSTALL/dlamchf77.f"
    if (try__ <= -(*emin)) {
#line 836 "../INSTALL/dlamchf77.f"
	lexp = try__;
#line 837 "../INSTALL/dlamchf77.f"
	++exbits;
#line 838 "../INSTALL/dlamchf77.f"
	goto L10;
#line 839 "../INSTALL/dlamchf77.f"
    }
#line 840 "../INSTALL/dlamchf77.f"
    if (lexp == -(*emin)) {
#line 841 "../INSTALL/dlamchf77.f"
	uexp = lexp;
#line 842 "../INSTALL/dlamchf77.f"
    } else {
#line 843 "../INSTALL/dlamchf77.f"
	uexp = try__;
#line 844 "../INSTALL/dlamchf77.f"
	++exbits;
#line 845 "../INSTALL/dlamchf77.f"
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
/*     than or equal to EMIN. EXBITS is the number of bits needed to */
/*     store the exponent. */

#line 851 "../INSTALL/dlamchf77.f"
    if (uexp + *emin > -lexp - *emin) {
#line 852 "../INSTALL/dlamchf77.f"
	expsum = lexp << 1;
#line 853 "../INSTALL/dlamchf77.f"
    } else {
#line 854 "../INSTALL/dlamchf77.f"
	expsum = uexp << 1;
#line 855 "../INSTALL/dlamchf77.f"
    }

/*     EXPSUM is the exponent range, approximately equal to */
/*     EMAX - EMIN + 1 . */

#line 860 "../INSTALL/dlamchf77.f"
    *emax = expsum + *emin - 1;
#line 861 "../INSTALL/dlamchf77.f"
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a */
/*     floating-point number. */

#line 866 "../INSTALL/dlamchf77.f"
    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a */
/*        floating-point number, which is unlikely, or some bits are */
/*        not used in the representation of numbers, which is possible, */
/*        (e.g. Cray machines) or the mantissa has an implicit bit, */
/*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the */
/*        most likely. We have to assume the last alternative. */
/*        If this is true, then we need to reduce EMAX by one because */
/*        there must be some way of representing zero in an implicit-bit */
/*        system. On machines like Cray, we are reducing EMAX by one */
/*        unnecessarily. */

#line 879 "../INSTALL/dlamchf77.f"
	--(*emax);
#line 880 "../INSTALL/dlamchf77.f"
    }

#line 882 "../INSTALL/dlamchf77.f"
    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent */
/*        for infinity and NaN. */

#line 887 "../INSTALL/dlamchf77.f"
	--(*emax);
#line 888 "../INSTALL/dlamchf77.f"
    }

/*     Now create RMAX, the largest machine number, which should */
/*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

/*     First compute 1.0 - BETA**(-P), being careful that the */
/*     result is less than 1.0 . */

#line 896 "../INSTALL/dlamchf77.f"
    recbas = 1. / *beta;
#line 897 "../INSTALL/dlamchf77.f"
    z__ = *beta - 1.;
#line 898 "../INSTALL/dlamchf77.f"
    y = 0.;
#line 899 "../INSTALL/dlamchf77.f"
    i__1 = *p;
#line 899 "../INSTALL/dlamchf77.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 900 "../INSTALL/dlamchf77.f"
	z__ *= recbas;
#line 901 "../INSTALL/dlamchf77.f"
	if (y < 1.) {
#line 901 "../INSTALL/dlamchf77.f"
	    oldy = y;
#line 901 "../INSTALL/dlamchf77.f"
	}
#line 903 "../INSTALL/dlamchf77.f"
	y = dlamc3_(&y, &z__);
#line 904 "../INSTALL/dlamchf77.f"
/* L20: */
#line 904 "../INSTALL/dlamchf77.f"
    }
#line 905 "../INSTALL/dlamchf77.f"
    if (y >= 1.) {
#line 905 "../INSTALL/dlamchf77.f"
	y = oldy;
#line 905 "../INSTALL/dlamchf77.f"
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

#line 910 "../INSTALL/dlamchf77.f"
    i__1 = *emax;
#line 910 "../INSTALL/dlamchf77.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 911 "../INSTALL/dlamchf77.f"
	d__1 = y * *beta;
#line 911 "../INSTALL/dlamchf77.f"
	y = dlamc3_(&d__1, &c_b32);
#line 912 "../INSTALL/dlamchf77.f"
/* L30: */
#line 912 "../INSTALL/dlamchf77.f"
    }

#line 914 "../INSTALL/dlamchf77.f"
    *rmax = y;
#line 915 "../INSTALL/dlamchf77.f"
    return 0;

/*     End of DLAMC5 */

} /* dlamc5_ */

