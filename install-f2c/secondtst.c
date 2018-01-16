#line 1 "../INSTALL/secondtst.f"
/* ../INSTALL/secondtst.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/secondtst.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__1000 = 1000;

/* > \brief \b SECONDTST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */


/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup auxOTHERcomputational */

/*  =====================================================================      PROGRAM SECONDTST */

/*  -- LAPACK test routine (version 3.8.0) -- */

/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/* ===================================================================== */

/*     .. Parameters .. */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 Time for \002,g10.3,\002 SAXPY ops = "
	    "\002,g10.3,\002 seconds\002)";
    static char fmt_9998[] = "(\002 SAXPY performance rate        = \002,g10"
	    ".3,\002 mflops \002)";
    static char fmt_9994[] = "(\002 *** Warning:  Time for operations was le"
	    "ss or equal\002,\002 than zero => timing in TESTING might be dub"
	    "ious\002)";
    static char fmt_9997[] = "(\002 Including SECOND, time        = \002,g10"
	    ".3,\002 seconds\002)";
    static char fmt_9996[] = "(\002 Average time for SECOND       = \002,g10"
	    ".3,\002 milliseconds\002)";
    static char fmt_9995[] = "(\002 Equivalent floating point ops = \002,g10"
	    ".3,\002 ops\002)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal x[1000], y[1000], t1, t2, avg, alpha, total;
    extern /* Subroutine */ int mysub_(integer *, doublereal *, doublereal *);
    extern doublereal second_(void);
    static doublereal tnosec;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, fmt_9999, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_9994, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_9997, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_9996, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_9995, 0 };


/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*    .. Figure TOTAL flops .. */
#line 56 "../INSTALL/secondtst.f"
    total = 1e8;

/*     Initialize X and Y */

#line 60 "../INSTALL/secondtst.f"
    for (i__ = 1; i__ <= 1000; ++i__) {
#line 61 "../INSTALL/secondtst.f"
	x[i__ - 1] = 1. / (doublereal) i__;
#line 62 "../INSTALL/secondtst.f"
	y[i__ - 1] = (doublereal) (1000 - i__) / 1e3;
#line 63 "../INSTALL/secondtst.f"
/* L10: */
#line 63 "../INSTALL/secondtst.f"
    }
#line 64 "../INSTALL/secondtst.f"
    alpha = .315;

/*     Time TOTAL SAXPY operations */

#line 68 "../INSTALL/secondtst.f"
    t1 = second_();
#line 69 "../INSTALL/secondtst.f"
    for (j = 1; j <= 50000; ++j) {
#line 70 "../INSTALL/secondtst.f"
	for (i__ = 1; i__ <= 1000; ++i__) {
#line 71 "../INSTALL/secondtst.f"
	    y[i__ - 1] += alpha * x[i__ - 1];
#line 72 "../INSTALL/secondtst.f"
/* L20: */
#line 72 "../INSTALL/secondtst.f"
	}
#line 73 "../INSTALL/secondtst.f"
	alpha = -alpha;
#line 74 "../INSTALL/secondtst.f"
/* L30: */
#line 74 "../INSTALL/secondtst.f"
    }
#line 75 "../INSTALL/secondtst.f"
    t2 = second_();
#line 76 "../INSTALL/secondtst.f"
    tnosec = t2 - t1;
#line 77 "../INSTALL/secondtst.f"
    s_wsfe(&io___10);
#line 77 "../INSTALL/secondtst.f"
    do_fio(&c__1, (char *)&total, (ftnlen)sizeof(doublereal));
#line 77 "../INSTALL/secondtst.f"
    do_fio(&c__1, (char *)&tnosec, (ftnlen)sizeof(doublereal));
#line 77 "../INSTALL/secondtst.f"
    e_wsfe();
#line 78 "../INSTALL/secondtst.f"
    if (tnosec > 0.) {
#line 79 "../INSTALL/secondtst.f"
	s_wsfe(&io___11);
#line 79 "../INSTALL/secondtst.f"
	d__1 = total / 1e6 / tnosec;
#line 79 "../INSTALL/secondtst.f"
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
#line 79 "../INSTALL/secondtst.f"
	e_wsfe();
#line 80 "../INSTALL/secondtst.f"
    } else {
#line 81 "../INSTALL/secondtst.f"
	s_wsfe(&io___12);
#line 81 "../INSTALL/secondtst.f"
	e_wsfe();
#line 82 "../INSTALL/secondtst.f"
    }

/*     Time TOTAL SAXPY operations with SECOND in the outer loop */

#line 86 "../INSTALL/secondtst.f"
    t1 = second_();
#line 87 "../INSTALL/secondtst.f"
    for (j = 1; j <= 50000; ++j) {
#line 88 "../INSTALL/secondtst.f"
	for (i__ = 1; i__ <= 1000; ++i__) {
#line 89 "../INSTALL/secondtst.f"
	    y[i__ - 1] += alpha * x[i__ - 1];
#line 90 "../INSTALL/secondtst.f"
/* L40: */
#line 90 "../INSTALL/secondtst.f"
	}
#line 91 "../INSTALL/secondtst.f"
	alpha = -alpha;
#line 92 "../INSTALL/secondtst.f"
	t2 = second_();
#line 93 "../INSTALL/secondtst.f"
/* L50: */
#line 93 "../INSTALL/secondtst.f"
    }

/*     Compute the time used in milliseconds used by an average call */
/*     to SECOND. */

#line 98 "../INSTALL/secondtst.f"
    s_wsfe(&io___13);
#line 98 "../INSTALL/secondtst.f"
    d__1 = t2 - t1;
#line 98 "../INSTALL/secondtst.f"
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
#line 98 "../INSTALL/secondtst.f"
    e_wsfe();
#line 99 "../INSTALL/secondtst.f"
    avg = (t2 - t1 - tnosec) * 1e3 / 5e4;
#line 100 "../INSTALL/secondtst.f"
    if (avg > 0.) {
#line 100 "../INSTALL/secondtst.f"
	s_wsfe(&io___15);
#line 100 "../INSTALL/secondtst.f"
	do_fio(&c__1, (char *)&avg, (ftnlen)sizeof(doublereal));
#line 100 "../INSTALL/secondtst.f"
	e_wsfe();
#line 100 "../INSTALL/secondtst.f"
    }

/*     Compute the equivalent number of floating point operations used */
/*     by an average call to SECOND. */

#line 106 "../INSTALL/secondtst.f"
    if (avg > 0. && tnosec > 0.) {
#line 106 "../INSTALL/secondtst.f"
	s_wsfe(&io___16);
#line 106 "../INSTALL/secondtst.f"
	d__1 = avg / 1000 * total / tnosec;
#line 106 "../INSTALL/secondtst.f"
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
#line 106 "../INSTALL/secondtst.f"
	e_wsfe();
#line 106 "../INSTALL/secondtst.f"
    }

#line 117 "../INSTALL/secondtst.f"
    mysub_(&c__1000, x, y);
#line 118 "../INSTALL/secondtst.f"
    return 0;
} /* MAIN__ */

/* Subroutine */ int mysub_(integer *n, doublereal *x, doublereal *y)
{
#line 122 "../INSTALL/secondtst.f"
    /* Parameter adjustments */
#line 122 "../INSTALL/secondtst.f"
    --y;
#line 122 "../INSTALL/secondtst.f"
    --x;
#line 122 "../INSTALL/secondtst.f"

#line 122 "../INSTALL/secondtst.f"
    /* Function Body */
#line 122 "../INSTALL/secondtst.f"
    return 0;
} /* mysub_ */

