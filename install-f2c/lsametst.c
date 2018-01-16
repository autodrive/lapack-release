#line 1 "../INSTALL/lsametst.f"
/* ../INSTALL/lsametst.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/lsametst.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* > \brief \b LSAMETST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*      PROGRAM LSAMETST */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERauxiliary */

/*  =====================================================================      PROGRAM LSAMETST */

/*  -- LAPACK test routine (version 3.7.0) -- */

/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*  ===================================================================== */
/*     .. Local Scalars .. */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 *** Error:  LSAME( \002,a1,\002, \002,"
	    "a1,\002) is .FALSE.\002)";
    static char fmt_9998[] = "(\002 *** Error:  LSAME( \002,a1,\002, \002,"
	    "a1,\002) is .TRUE.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    static integer i1, i2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, fmt_9999, 0 };
    static cilist io___6 = { 0, 6, 0, fmt_9999, 0 };
    static cilist io___7 = { 0, 6, 0, fmt_9999, 0 };
    static cilist io___8 = { 0, 6, 0, fmt_9999, 0 };
    static cilist io___9 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___10 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_9998, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };


/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */


/*     Determine the character set. */

#line 50 "../INSTALL/lsametst.f"
    i1 = 'A';
#line 51 "../INSTALL/lsametst.f"
    i2 = 'a';
#line 52 "../INSTALL/lsametst.f"
    if (i2 - i1 == 32) {
#line 53 "../INSTALL/lsametst.f"
	s_wsle(&io___3);
#line 53 "../INSTALL/lsametst.f"
	do_lio(&c__9, &c__1, " ASCII character set", (ftnlen)20);
#line 53 "../INSTALL/lsametst.f"
	e_wsle();
#line 54 "../INSTALL/lsametst.f"
    } else {
#line 55 "../INSTALL/lsametst.f"
	s_wsle(&io___4);
#line 55 "../INSTALL/lsametst.f"
	do_lio(&c__9, &c__1, " Non-ASCII character set, IOFF should be ", (
		ftnlen)41);
#line 55 "../INSTALL/lsametst.f"
	i__1 = i2 - i1;
#line 55 "../INSTALL/lsametst.f"
	do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(integer));
#line 55 "../INSTALL/lsametst.f"
	e_wsle();
#line 56 "../INSTALL/lsametst.f"
    }

/*     Test LSAME. */

#line 60 "../INSTALL/lsametst.f"
    if (! lsame_("A", "A", (ftnlen)1, (ftnlen)1)) {
#line 60 "../INSTALL/lsametst.f"
	s_wsfe(&io___5);
#line 60 "../INSTALL/lsametst.f"
	do_fio(&c__1, "A", (ftnlen)1);
#line 60 "../INSTALL/lsametst.f"
	do_fio(&c__1, "A", (ftnlen)1);
#line 60 "../INSTALL/lsametst.f"
	e_wsfe();
#line 60 "../INSTALL/lsametst.f"
    }
#line 62 "../INSTALL/lsametst.f"
    if (! lsame_("A", "a", (ftnlen)1, (ftnlen)1)) {
#line 62 "../INSTALL/lsametst.f"
	s_wsfe(&io___6);
#line 62 "../INSTALL/lsametst.f"
	do_fio(&c__1, "A", (ftnlen)1);
#line 62 "../INSTALL/lsametst.f"
	do_fio(&c__1, "a", (ftnlen)1);
#line 62 "../INSTALL/lsametst.f"
	e_wsfe();
#line 62 "../INSTALL/lsametst.f"
    }
#line 64 "../INSTALL/lsametst.f"
    if (! lsame_("a", "A", (ftnlen)1, (ftnlen)1)) {
#line 64 "../INSTALL/lsametst.f"
	s_wsfe(&io___7);
#line 64 "../INSTALL/lsametst.f"
	do_fio(&c__1, "a", (ftnlen)1);
#line 64 "../INSTALL/lsametst.f"
	do_fio(&c__1, "A", (ftnlen)1);
#line 64 "../INSTALL/lsametst.f"
	e_wsfe();
#line 64 "../INSTALL/lsametst.f"
    }
#line 66 "../INSTALL/lsametst.f"
    if (! lsame_("a", "a", (ftnlen)1, (ftnlen)1)) {
#line 66 "../INSTALL/lsametst.f"
	s_wsfe(&io___8);
#line 66 "../INSTALL/lsametst.f"
	do_fio(&c__1, "a", (ftnlen)1);
#line 66 "../INSTALL/lsametst.f"
	do_fio(&c__1, "a", (ftnlen)1);
#line 66 "../INSTALL/lsametst.f"
	e_wsfe();
#line 66 "../INSTALL/lsametst.f"
    }
#line 68 "../INSTALL/lsametst.f"
    if (lsame_("A", "B", (ftnlen)1, (ftnlen)1)) {
#line 68 "../INSTALL/lsametst.f"
	s_wsfe(&io___9);
#line 68 "../INSTALL/lsametst.f"
	do_fio(&c__1, "A", (ftnlen)1);
#line 68 "../INSTALL/lsametst.f"
	do_fio(&c__1, "B", (ftnlen)1);
#line 68 "../INSTALL/lsametst.f"
	e_wsfe();
#line 68 "../INSTALL/lsametst.f"
    }
#line 70 "../INSTALL/lsametst.f"
    if (lsame_("A", "b", (ftnlen)1, (ftnlen)1)) {
#line 70 "../INSTALL/lsametst.f"
	s_wsfe(&io___10);
#line 70 "../INSTALL/lsametst.f"
	do_fio(&c__1, "A", (ftnlen)1);
#line 70 "../INSTALL/lsametst.f"
	do_fio(&c__1, "b", (ftnlen)1);
#line 70 "../INSTALL/lsametst.f"
	e_wsfe();
#line 70 "../INSTALL/lsametst.f"
    }
#line 72 "../INSTALL/lsametst.f"
    if (lsame_("a", "B", (ftnlen)1, (ftnlen)1)) {
#line 72 "../INSTALL/lsametst.f"
	s_wsfe(&io___11);
#line 72 "../INSTALL/lsametst.f"
	do_fio(&c__1, "a", (ftnlen)1);
#line 72 "../INSTALL/lsametst.f"
	do_fio(&c__1, "B", (ftnlen)1);
#line 72 "../INSTALL/lsametst.f"
	e_wsfe();
#line 72 "../INSTALL/lsametst.f"
    }
#line 74 "../INSTALL/lsametst.f"
    if (lsame_("a", "b", (ftnlen)1, (ftnlen)1)) {
#line 74 "../INSTALL/lsametst.f"
	s_wsfe(&io___12);
#line 74 "../INSTALL/lsametst.f"
	do_fio(&c__1, "a", (ftnlen)1);
#line 74 "../INSTALL/lsametst.f"
	do_fio(&c__1, "b", (ftnlen)1);
#line 74 "../INSTALL/lsametst.f"
	e_wsfe();
#line 74 "../INSTALL/lsametst.f"
    }
#line 76 "../INSTALL/lsametst.f"
    if (lsame_("O", "/", (ftnlen)1, (ftnlen)1)) {
#line 76 "../INSTALL/lsametst.f"
	s_wsfe(&io___13);
#line 76 "../INSTALL/lsametst.f"
	do_fio(&c__1, "O", (ftnlen)1);
#line 76 "../INSTALL/lsametst.f"
	do_fio(&c__1, "/", (ftnlen)1);
#line 76 "../INSTALL/lsametst.f"
	e_wsfe();
#line 76 "../INSTALL/lsametst.f"
    }
#line 78 "../INSTALL/lsametst.f"
    if (lsame_("/", "O", (ftnlen)1, (ftnlen)1)) {
#line 78 "../INSTALL/lsametst.f"
	s_wsfe(&io___14);
#line 78 "../INSTALL/lsametst.f"
	do_fio(&c__1, "/", (ftnlen)1);
#line 78 "../INSTALL/lsametst.f"
	do_fio(&c__1, "O", (ftnlen)1);
#line 78 "../INSTALL/lsametst.f"
	e_wsfe();
#line 78 "../INSTALL/lsametst.f"
    }
#line 80 "../INSTALL/lsametst.f"
    if (lsame_("o", "/", (ftnlen)1, (ftnlen)1)) {
#line 80 "../INSTALL/lsametst.f"
	s_wsfe(&io___15);
#line 80 "../INSTALL/lsametst.f"
	do_fio(&c__1, "o", (ftnlen)1);
#line 80 "../INSTALL/lsametst.f"
	do_fio(&c__1, "/", (ftnlen)1);
#line 80 "../INSTALL/lsametst.f"
	e_wsfe();
#line 80 "../INSTALL/lsametst.f"
    }
#line 82 "../INSTALL/lsametst.f"
    if (lsame_("/", "o", (ftnlen)1, (ftnlen)1)) {
#line 82 "../INSTALL/lsametst.f"
	s_wsfe(&io___16);
#line 82 "../INSTALL/lsametst.f"
	do_fio(&c__1, "/", (ftnlen)1);
#line 82 "../INSTALL/lsametst.f"
	do_fio(&c__1, "o", (ftnlen)1);
#line 82 "../INSTALL/lsametst.f"
	e_wsfe();
#line 82 "../INSTALL/lsametst.f"
    }
#line 84 "../INSTALL/lsametst.f"
    s_wsle(&io___17);
#line 84 "../INSTALL/lsametst.f"
    do_lio(&c__9, &c__1, " Tests completed", (ftnlen)16);
#line 84 "../INSTALL/lsametst.f"
    e_wsle();

#line 88 "../INSTALL/lsametst.f"
    return 0;
} /* MAIN__ */

