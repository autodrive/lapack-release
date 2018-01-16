#line 1 "../INSTALL/slamchtst.f"
/* ../INSTALL/slamchtst.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/slamchtst.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;

/* > \brief \b SLAMCHTST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */


/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  =====================================================================      PROGRAM SLAMCHTST */

/*  -- LAPACK test routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/* ===================================================================== */

/*     .. Local Scalars .. */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal t, rnd, eps, base, emin, prec, emax, rmin, rmax, sfmin;
    extern doublereal slamch_(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };


/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 40 "../INSTALL/slamchtst.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 41 "../INSTALL/slamchtst.f"
    sfmin = slamch_("Safe minimum", (ftnlen)12);
#line 42 "../INSTALL/slamchtst.f"
    base = slamch_("Base", (ftnlen)4);
#line 43 "../INSTALL/slamchtst.f"
    prec = slamch_("Precision", (ftnlen)9);
#line 44 "../INSTALL/slamchtst.f"
    t = slamch_("Number of digits in mantissa", (ftnlen)28);
#line 45 "../INSTALL/slamchtst.f"
    rnd = slamch_("Rounding mode", (ftnlen)13);
#line 46 "../INSTALL/slamchtst.f"
    emin = slamch_("Minimum exponent", (ftnlen)16);
#line 47 "../INSTALL/slamchtst.f"
    rmin = slamch_("Underflow threshold", (ftnlen)19);
#line 48 "../INSTALL/slamchtst.f"
    emax = slamch_("Largest exponent", (ftnlen)16);
#line 49 "../INSTALL/slamchtst.f"
    rmax = slamch_("Overflow threshold", (ftnlen)18);

#line 51 "../INSTALL/slamchtst.f"
    s_wsle(&io___11);
#line 51 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Epsilon                      = ", (ftnlen)32);
#line 51 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
#line 51 "../INSTALL/slamchtst.f"
    e_wsle();
#line 52 "../INSTALL/slamchtst.f"
    s_wsle(&io___12);
#line 52 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Safe minimum                 = ", (ftnlen)32);
#line 52 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&sfmin, (ftnlen)sizeof(doublereal));
#line 52 "../INSTALL/slamchtst.f"
    e_wsle();
#line 53 "../INSTALL/slamchtst.f"
    s_wsle(&io___13);
#line 53 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Base                         = ", (ftnlen)32);
#line 53 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&base, (ftnlen)sizeof(doublereal));
#line 53 "../INSTALL/slamchtst.f"
    e_wsle();
#line 54 "../INSTALL/slamchtst.f"
    s_wsle(&io___14);
#line 54 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Precision                    = ", (ftnlen)32);
#line 54 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&prec, (ftnlen)sizeof(doublereal));
#line 54 "../INSTALL/slamchtst.f"
    e_wsle();
#line 55 "../INSTALL/slamchtst.f"
    s_wsle(&io___15);
#line 55 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Number of digits in mantissa = ", (ftnlen)32);
#line 55 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&t, (ftnlen)sizeof(doublereal));
#line 55 "../INSTALL/slamchtst.f"
    e_wsle();
#line 56 "../INSTALL/slamchtst.f"
    s_wsle(&io___16);
#line 56 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Rounding mode                = ", (ftnlen)32);
#line 56 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&rnd, (ftnlen)sizeof(doublereal));
#line 56 "../INSTALL/slamchtst.f"
    e_wsle();
#line 57 "../INSTALL/slamchtst.f"
    s_wsle(&io___17);
#line 57 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Minimum exponent             = ", (ftnlen)32);
#line 57 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&emin, (ftnlen)sizeof(doublereal));
#line 57 "../INSTALL/slamchtst.f"
    e_wsle();
#line 58 "../INSTALL/slamchtst.f"
    s_wsle(&io___18);
#line 58 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Underflow threshold          = ", (ftnlen)32);
#line 58 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&rmin, (ftnlen)sizeof(doublereal));
#line 58 "../INSTALL/slamchtst.f"
    e_wsle();
#line 59 "../INSTALL/slamchtst.f"
    s_wsle(&io___19);
#line 59 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Largest exponent             = ", (ftnlen)32);
#line 59 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&emax, (ftnlen)sizeof(doublereal));
#line 59 "../INSTALL/slamchtst.f"
    e_wsle();
#line 60 "../INSTALL/slamchtst.f"
    s_wsle(&io___20);
#line 60 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Overflow threshold           = ", (ftnlen)32);
#line 60 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&rmax, (ftnlen)sizeof(doublereal));
#line 60 "../INSTALL/slamchtst.f"
    e_wsle();
#line 61 "../INSTALL/slamchtst.f"
    s_wsle(&io___21);
#line 61 "../INSTALL/slamchtst.f"
    do_lio(&c__9, &c__1, " Reciprocal of safe minimum   = ", (ftnlen)32);
#line 61 "../INSTALL/slamchtst.f"
    d__1 = 1 / sfmin;
#line 61 "../INSTALL/slamchtst.f"
    do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
#line 61 "../INSTALL/slamchtst.f"
    e_wsle();

#line 63 "../INSTALL/slamchtst.f"
    return 0;
} /* MAIN__ */

