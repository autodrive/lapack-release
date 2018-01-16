#line 1 "../INSTALL/tstiee.f"
/* ../INSTALL/tstiee.f -- translated by f2c (version 20100827).
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

#line 1 "../INSTALL/tstiee.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__10 = 10;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__11 = 11;
static integer c__0 = 0;
static doublereal c_b227 = 0.;
static doublereal c_b228 = 1.;

/* > \brief \b TSTIEE */

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

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
/* Main program */ int MAIN__(void)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer ieeeok;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___8 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };



/*  -- LAPACK test routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. External Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 36 "../INSTALL/tstiee.f"
    s_wsle(&io___1);
#line 36 "../INSTALL/tstiee.f"
    do_lio(&c__9, &c__1, "We are about to check whether infinity arithmetic", 
	    (ftnlen)49);
#line 36 "../INSTALL/tstiee.f"
    e_wsle();
#line 38 "../INSTALL/tstiee.f"
    s_wsle(&io___2);
#line 38 "../INSTALL/tstiee.f"
    do_lio(&c__9, &c__1, "can be trusted.  If this test hangs, set", (ftnlen)
	    40);
#line 38 "../INSTALL/tstiee.f"
    e_wsle();
#line 39 "../INSTALL/tstiee.f"
    s_wsle(&io___3);
#line 39 "../INSTALL/tstiee.f"
    do_lio(&c__9, &c__1, "ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f", (
	    ftnlen)48);
#line 39 "../INSTALL/tstiee.f"
    e_wsle();

#line 42 "../INSTALL/tstiee.f"
    ieeeok = ilaenv_(&c__10, "ILAENV", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);
#line 43 "../INSTALL/tstiee.f"
    s_wsle(&io___5);
#line 43 "../INSTALL/tstiee.f"
    e_wsle();

#line 45 "../INSTALL/tstiee.f"
    if (ieeeok == 0) {
#line 46 "../INSTALL/tstiee.f"
	s_wsle(&io___6);
#line 46 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "Infinity arithmetic did not perform per the ie"\
		"ee spec", (ftnlen)53);
#line 46 "../INSTALL/tstiee.f"
	e_wsle();
#line 48 "../INSTALL/tstiee.f"
    } else {
#line 49 "../INSTALL/tstiee.f"
	s_wsle(&io___7);
#line 49 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "Infinity arithmetic performed as per the ieee "\
		"spec.", (ftnlen)51);
#line 49 "../INSTALL/tstiee.f"
	e_wsle();
#line 51 "../INSTALL/tstiee.f"
	s_wsle(&io___8);
#line 51 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "However, this is not an exhaustive test and do"\
		"es not", (ftnlen)52);
#line 51 "../INSTALL/tstiee.f"
	e_wsle();
#line 53 "../INSTALL/tstiee.f"
	s_wsle(&io___9);
#line 53 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "guarantee that infinity arithmetic meets the", (
		ftnlen)44);
#line 53 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, " ieee spec.", (ftnlen)11);
#line 53 "../INSTALL/tstiee.f"
	e_wsle();
#line 56 "../INSTALL/tstiee.f"
    }

#line 58 "../INSTALL/tstiee.f"
    s_wsle(&io___10);
#line 58 "../INSTALL/tstiee.f"
    e_wsle();
#line 59 "../INSTALL/tstiee.f"
    s_wsle(&io___11);
#line 59 "../INSTALL/tstiee.f"
    do_lio(&c__9, &c__1, "We are about to check whether NaN arithmetic", (
	    ftnlen)44);
#line 59 "../INSTALL/tstiee.f"
    e_wsle();
#line 61 "../INSTALL/tstiee.f"
    s_wsle(&io___12);
#line 61 "../INSTALL/tstiee.f"
    do_lio(&c__9, &c__1, "can be trusted.  If this test hangs, set", (ftnlen)
	    40);
#line 61 "../INSTALL/tstiee.f"
    e_wsle();
#line 62 "../INSTALL/tstiee.f"
    s_wsle(&io___13);
#line 62 "../INSTALL/tstiee.f"
    do_lio(&c__9, &c__1, "ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f", (
	    ftnlen)48);
#line 62 "../INSTALL/tstiee.f"
    e_wsle();
#line 64 "../INSTALL/tstiee.f"
    ieeeok = ilaenv_(&c__11, "ILAENV", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

#line 66 "../INSTALL/tstiee.f"
    s_wsle(&io___14);
#line 66 "../INSTALL/tstiee.f"
    e_wsle();
#line 67 "../INSTALL/tstiee.f"
    if (ieeeok == 0) {
#line 68 "../INSTALL/tstiee.f"
	s_wsle(&io___15);
#line 68 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "NaN arithmetic did not perform per the ieee sp"\
		"ec", (ftnlen)48);
#line 68 "../INSTALL/tstiee.f"
	e_wsle();
#line 70 "../INSTALL/tstiee.f"
    } else {
#line 71 "../INSTALL/tstiee.f"
	s_wsle(&io___16);
#line 71 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "NaN arithmetic performed as per the ieee", (
		ftnlen)40);
#line 71 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, " spec.", (ftnlen)6);
#line 71 "../INSTALL/tstiee.f"
	e_wsle();
#line 73 "../INSTALL/tstiee.f"
	s_wsle(&io___17);
#line 73 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "However, this is not an exhaustive test and do"\
		"es not", (ftnlen)52);
#line 73 "../INSTALL/tstiee.f"
	e_wsle();
#line 75 "../INSTALL/tstiee.f"
	s_wsle(&io___18);
#line 75 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, "guarantee that NaN arithmetic meets the", (
		ftnlen)39);
#line 75 "../INSTALL/tstiee.f"
	do_lio(&c__9, &c__1, " ieee spec.", (ftnlen)11);
#line 75 "../INSTALL/tstiee.f"
	e_wsle();
#line 77 "../INSTALL/tstiee.f"
    }
#line 78 "../INSTALL/tstiee.f"
    s_wsle(&io___19);
#line 78 "../INSTALL/tstiee.f"
    e_wsle();

#line 80 "../INSTALL/tstiee.f"
    return 0;
} /* MAIN__ */

integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static char c1[1], c2[2], c3[3], c4[2];
    static integer ic, nb, iz, nx;
    static logical cname, sname;
    static integer nbmin;
    extern integer ieeeck_(integer *, doublereal *, doublereal *);
    static char subnam[6];


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV is called from the LAPACK routines to choose problem-dependent */
/*  parameters for the local environment.  See ISPEC for a description of */
/*  the parameters. */

/*  This version provides a set of parameters which should give good, */
/*  but not optimal, performance on many of the currently available */
/*  computers.  Users are encouraged to modify this subroutine to set */
/*  the tuning parameters for their particular machine using the option */
/*  and problem size information in the arguments. */

/*  This routine will not function correctly if it is converted to all */
/*  lower case.  Converting it to all upper case is allowed. */

/*  Arguments: */
/*  ========== */

/*  ISPEC   (input) INTEGER */
/*          Specifies the parameter to be returned as the value of */
/*          ILAENV. */
/*          = 1: the optimal blocksize; if this value is 1, an unblocked */
/*               algorithm will give the best performance. */
/*          = 2: the minimum block size for which the block routine */
/*               should be used; if the usable block size is less than */
/*               this value, an unblocked routine should be used. */
/*          = 3: the crossover point (in a block routine, for N less */
/*               than this value, an unblocked routine should be used) */
/*          = 4: the number of shifts, used in the nonsymmetric */
/*               eigenvalue routines */
/*          = 5: the minimum column dimension for blocking to be used; */
/*               rectangular blocks must have dimension at least k by m, */
/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
/*          = 6: the crossover point for the SVD (when reducing an m by n */
/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
/*               this value, a QR factorization is used first to reduce */
/*               the matrix to a triangular form.) */
/*          = 7: the number of processors */
/*          = 8: the crossover point for the multishift QR and QZ methods */
/*               for nonsymmetric eigenvalue problems. */
/*          = 9: maximum size of the subproblems at the bottom of the */
/*               computation tree in the divide-and-conquer algorithm */
/*               (used by xGELSD and xGESDD) */
/*          =10: ieee NaN arithmetic can be trusted not to trap */
/*          =11: infinity arithmetic can be trusted not to trap */

/*  NAME    (input) CHARACTER*(*) */
/*          The name of the calling subroutine, in either upper case or */
/*          lower case. */

/*  OPTS    (input) CHARACTER*(*) */
/*          The character options to the subroutine NAME, concatenated */
/*          into a single character string.  For example, UPLO = 'U', */
/*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/*          be specified as OPTS = 'UTN'. */

/*  N1      (input) INTEGER */
/*  N2      (input) INTEGER */
/*  N3      (input) INTEGER */
/*  N4      (input) INTEGER */
/*          Problem dimensions for the subroutine NAME; these may not all */
/*          be required. */

/* (ILAENV) (output) INTEGER */
/*          >= 0: the value of the parameter specified by ISPEC */
/*          < 0:  if ILAENV = -k, the k-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  The following conventions have been used when calling ILAENV from the */
/*  LAPACK routines: */
/*  1)  OPTS is a concatenation of all of the character options to */
/*      subroutine NAME, in the same order that they appear in the */
/*      argument list for NAME, even if they are not used in determining */
/*      the value of the parameter specified by ISPEC. */
/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
/*      that they appear in the argument list for NAME.  N1 is used */
/*      first, N2 second, and so on, and unused problem dimensions are */
/*      passed a value of -1. */
/*  3)  The parameter value returned by ILAENV is checked for validity in */
/*      the calling subroutine.  For example, ILAENV is used to retrieve */
/*      the optimal blocksize for STRTRI as follows: */

/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 200 "../INSTALL/tstiee.f"
    switch (*ispec) {
#line 200 "../INSTALL/tstiee.f"
	case 1:  goto L100;
#line 200 "../INSTALL/tstiee.f"
	case 2:  goto L100;
#line 200 "../INSTALL/tstiee.f"
	case 3:  goto L100;
#line 200 "../INSTALL/tstiee.f"
	case 4:  goto L400;
#line 200 "../INSTALL/tstiee.f"
	case 5:  goto L500;
#line 200 "../INSTALL/tstiee.f"
	case 6:  goto L600;
#line 200 "../INSTALL/tstiee.f"
	case 7:  goto L700;
#line 200 "../INSTALL/tstiee.f"
	case 8:  goto L800;
#line 200 "../INSTALL/tstiee.f"
	case 9:  goto L900;
#line 200 "../INSTALL/tstiee.f"
	case 10:  goto L1000;
#line 200 "../INSTALL/tstiee.f"
	case 11:  goto L1100;
#line 200 "../INSTALL/tstiee.f"
    }

/*     Invalid value for ISPEC */

#line 205 "../INSTALL/tstiee.f"
    ret_val = -1;
#line 206 "../INSTALL/tstiee.f"
    return ret_val;

#line 208 "../INSTALL/tstiee.f"
L100:

/*     Convert NAME to upper case if the first character is lower case. */

#line 212 "../INSTALL/tstiee.f"
    ret_val = 1;
#line 213 "../INSTALL/tstiee.f"
    s_copy(subnam, name__, (ftnlen)6, name_len);
#line 214 "../INSTALL/tstiee.f"
    ic = *(unsigned char *)subnam;
#line 215 "../INSTALL/tstiee.f"
    iz = 'Z';
#line 216 "../INSTALL/tstiee.f"
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

#line 220 "../INSTALL/tstiee.f"
	if (ic >= 97 && ic <= 122) {
#line 221 "../INSTALL/tstiee.f"
	    *(unsigned char *)subnam = (char) (ic - 32);
#line 222 "../INSTALL/tstiee.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 223 "../INSTALL/tstiee.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 224 "../INSTALL/tstiee.f"
		if (ic >= 97 && ic <= 122) {
#line 224 "../INSTALL/tstiee.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 224 "../INSTALL/tstiee.f"
		}
#line 226 "../INSTALL/tstiee.f"
/* L10: */
#line 226 "../INSTALL/tstiee.f"
	    }
#line 227 "../INSTALL/tstiee.f"
	}

#line 229 "../INSTALL/tstiee.f"
    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

#line 233 "../INSTALL/tstiee.f"
	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
#line 236 "../INSTALL/tstiee.f"
	    *(unsigned char *)subnam = (char) (ic + 64);
#line 237 "../INSTALL/tstiee.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 238 "../INSTALL/tstiee.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 239 "../INSTALL/tstiee.f"
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
#line 239 "../INSTALL/tstiee.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
#line 239 "../INSTALL/tstiee.f"
		}
#line 243 "../INSTALL/tstiee.f"
/* L20: */
#line 243 "../INSTALL/tstiee.f"
	    }
#line 244 "../INSTALL/tstiee.f"
	}

#line 246 "../INSTALL/tstiee.f"
    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

#line 250 "../INSTALL/tstiee.f"
	if (ic >= 225 && ic <= 250) {
#line 251 "../INSTALL/tstiee.f"
	    *(unsigned char *)subnam = (char) (ic - 32);
#line 252 "../INSTALL/tstiee.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 253 "../INSTALL/tstiee.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 254 "../INSTALL/tstiee.f"
		if (ic >= 225 && ic <= 250) {
#line 254 "../INSTALL/tstiee.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 254 "../INSTALL/tstiee.f"
		}
#line 256 "../INSTALL/tstiee.f"
/* L30: */
#line 256 "../INSTALL/tstiee.f"
	    }
#line 257 "../INSTALL/tstiee.f"
	}
#line 258 "../INSTALL/tstiee.f"
    }

#line 260 "../INSTALL/tstiee.f"
    *(unsigned char *)c1 = *(unsigned char *)subnam;
#line 261 "../INSTALL/tstiee.f"
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
#line 262 "../INSTALL/tstiee.f"
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
#line 263 "../INSTALL/tstiee.f"
    if (! (cname || sname)) {
#line 263 "../INSTALL/tstiee.f"
	return ret_val;
#line 263 "../INSTALL/tstiee.f"
    }
#line 265 "../INSTALL/tstiee.f"
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
#line 266 "../INSTALL/tstiee.f"
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
#line 267 "../INSTALL/tstiee.f"
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

#line 269 "../INSTALL/tstiee.f"
    switch (*ispec) {
#line 269 "../INSTALL/tstiee.f"
	case 1:  goto L110;
#line 269 "../INSTALL/tstiee.f"
	case 2:  goto L200;
#line 269 "../INSTALL/tstiee.f"
	case 3:  goto L300;
#line 269 "../INSTALL/tstiee.f"
    }

#line 271 "../INSTALL/tstiee.f"
L110:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

#line 279 "../INSTALL/tstiee.f"
    nb = 1;

#line 281 "../INSTALL/tstiee.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 282 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 283 "../INSTALL/tstiee.f"
	    if (sname) {
#line 284 "../INSTALL/tstiee.f"
		nb = 64;
#line 285 "../INSTALL/tstiee.f"
	    } else {
#line 286 "../INSTALL/tstiee.f"
		nb = 64;
#line 287 "../INSTALL/tstiee.f"
	    }
#line 288 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
#line 290 "../INSTALL/tstiee.f"
	    if (sname) {
#line 291 "../INSTALL/tstiee.f"
		nb = 32;
#line 292 "../INSTALL/tstiee.f"
	    } else {
#line 293 "../INSTALL/tstiee.f"
		nb = 32;
#line 294 "../INSTALL/tstiee.f"
	    }
#line 295 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 296 "../INSTALL/tstiee.f"
	    if (sname) {
#line 297 "../INSTALL/tstiee.f"
		nb = 32;
#line 298 "../INSTALL/tstiee.f"
	    } else {
#line 299 "../INSTALL/tstiee.f"
		nb = 32;
#line 300 "../INSTALL/tstiee.f"
	    }
#line 301 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 302 "../INSTALL/tstiee.f"
	    if (sname) {
#line 303 "../INSTALL/tstiee.f"
		nb = 32;
#line 304 "../INSTALL/tstiee.f"
	    } else {
#line 305 "../INSTALL/tstiee.f"
		nb = 32;
#line 306 "../INSTALL/tstiee.f"
	    }
#line 307 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 308 "../INSTALL/tstiee.f"
	    if (sname) {
#line 309 "../INSTALL/tstiee.f"
		nb = 64;
#line 310 "../INSTALL/tstiee.f"
	    } else {
#line 311 "../INSTALL/tstiee.f"
		nb = 64;
#line 312 "../INSTALL/tstiee.f"
	    }
#line 313 "../INSTALL/tstiee.f"
	}
#line 314 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
#line 315 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 316 "../INSTALL/tstiee.f"
	    if (sname) {
#line 317 "../INSTALL/tstiee.f"
		nb = 64;
#line 318 "../INSTALL/tstiee.f"
	    } else {
#line 319 "../INSTALL/tstiee.f"
		nb = 64;
#line 320 "../INSTALL/tstiee.f"
	    }
#line 321 "../INSTALL/tstiee.f"
	}
#line 322 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 323 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 324 "../INSTALL/tstiee.f"
	    if (sname) {
#line 325 "../INSTALL/tstiee.f"
		nb = 64;
#line 326 "../INSTALL/tstiee.f"
	    } else {
#line 327 "../INSTALL/tstiee.f"
		nb = 64;
#line 328 "../INSTALL/tstiee.f"
	    }
#line 329 "../INSTALL/tstiee.f"
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 330 "../INSTALL/tstiee.f"
	    nb = 32;
#line 331 "../INSTALL/tstiee.f"
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
#line 332 "../INSTALL/tstiee.f"
	    nb = 64;
#line 333 "../INSTALL/tstiee.f"
	}
#line 334 "../INSTALL/tstiee.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 335 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 336 "../INSTALL/tstiee.f"
	    nb = 64;
#line 337 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 338 "../INSTALL/tstiee.f"
	    nb = 32;
#line 339 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
#line 340 "../INSTALL/tstiee.f"
	    nb = 64;
#line 341 "../INSTALL/tstiee.f"
	}
#line 342 "../INSTALL/tstiee.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 343 "../INSTALL/tstiee.f"
	if (*(unsigned char *)c3 == 'G') {
#line 344 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 347 "../INSTALL/tstiee.f"
		nb = 32;
#line 348 "../INSTALL/tstiee.f"
	    }
#line 349 "../INSTALL/tstiee.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 350 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 353 "../INSTALL/tstiee.f"
		nb = 32;
#line 354 "../INSTALL/tstiee.f"
	    }
#line 355 "../INSTALL/tstiee.f"
	}
#line 356 "../INSTALL/tstiee.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 357 "../INSTALL/tstiee.f"
	if (*(unsigned char *)c3 == 'G') {
#line 358 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 361 "../INSTALL/tstiee.f"
		nb = 32;
#line 362 "../INSTALL/tstiee.f"
	    }
#line 363 "../INSTALL/tstiee.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 364 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 367 "../INSTALL/tstiee.f"
		nb = 32;
#line 368 "../INSTALL/tstiee.f"
	    }
#line 369 "../INSTALL/tstiee.f"
	}
#line 370 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
#line 371 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 372 "../INSTALL/tstiee.f"
	    if (sname) {
#line 373 "../INSTALL/tstiee.f"
		if (*n4 <= 64) {
#line 374 "../INSTALL/tstiee.f"
		    nb = 1;
#line 375 "../INSTALL/tstiee.f"
		} else {
#line 376 "../INSTALL/tstiee.f"
		    nb = 32;
#line 377 "../INSTALL/tstiee.f"
		}
#line 378 "../INSTALL/tstiee.f"
	    } else {
#line 379 "../INSTALL/tstiee.f"
		if (*n4 <= 64) {
#line 380 "../INSTALL/tstiee.f"
		    nb = 1;
#line 381 "../INSTALL/tstiee.f"
		} else {
#line 382 "../INSTALL/tstiee.f"
		    nb = 32;
#line 383 "../INSTALL/tstiee.f"
		}
#line 384 "../INSTALL/tstiee.f"
	    }
#line 385 "../INSTALL/tstiee.f"
	}
#line 386 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
#line 387 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 388 "../INSTALL/tstiee.f"
	    if (sname) {
#line 389 "../INSTALL/tstiee.f"
		if (*n2 <= 64) {
#line 390 "../INSTALL/tstiee.f"
		    nb = 1;
#line 391 "../INSTALL/tstiee.f"
		} else {
#line 392 "../INSTALL/tstiee.f"
		    nb = 32;
#line 393 "../INSTALL/tstiee.f"
		}
#line 394 "../INSTALL/tstiee.f"
	    } else {
#line 395 "../INSTALL/tstiee.f"
		if (*n2 <= 64) {
#line 396 "../INSTALL/tstiee.f"
		    nb = 1;
#line 397 "../INSTALL/tstiee.f"
		} else {
#line 398 "../INSTALL/tstiee.f"
		    nb = 32;
#line 399 "../INSTALL/tstiee.f"
		}
#line 400 "../INSTALL/tstiee.f"
	    }
#line 401 "../INSTALL/tstiee.f"
	}
#line 402 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
#line 403 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 404 "../INSTALL/tstiee.f"
	    if (sname) {
#line 405 "../INSTALL/tstiee.f"
		nb = 64;
#line 406 "../INSTALL/tstiee.f"
	    } else {
#line 407 "../INSTALL/tstiee.f"
		nb = 64;
#line 408 "../INSTALL/tstiee.f"
	    }
#line 409 "../INSTALL/tstiee.f"
	}
#line 410 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
#line 411 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
#line 412 "../INSTALL/tstiee.f"
	    if (sname) {
#line 413 "../INSTALL/tstiee.f"
		nb = 64;
#line 414 "../INSTALL/tstiee.f"
	    } else {
#line 415 "../INSTALL/tstiee.f"
		nb = 64;
#line 416 "../INSTALL/tstiee.f"
	    }
#line 417 "../INSTALL/tstiee.f"
	}
#line 418 "../INSTALL/tstiee.f"
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
#line 419 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
#line 420 "../INSTALL/tstiee.f"
	    nb = 1;
#line 421 "../INSTALL/tstiee.f"
	}
#line 422 "../INSTALL/tstiee.f"
    }
#line 423 "../INSTALL/tstiee.f"
    ret_val = nb;
#line 424 "../INSTALL/tstiee.f"
    return ret_val;

#line 426 "../INSTALL/tstiee.f"
L200:

/*     ISPEC = 2:  minimum block size */

#line 430 "../INSTALL/tstiee.f"
    nbmin = 2;
#line 431 "../INSTALL/tstiee.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 432 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
#line 434 "../INSTALL/tstiee.f"
	    if (sname) {
#line 435 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 436 "../INSTALL/tstiee.f"
	    } else {
#line 437 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 438 "../INSTALL/tstiee.f"
	    }
#line 439 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 440 "../INSTALL/tstiee.f"
	    if (sname) {
#line 441 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 442 "../INSTALL/tstiee.f"
	    } else {
#line 443 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 444 "../INSTALL/tstiee.f"
	    }
#line 445 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 446 "../INSTALL/tstiee.f"
	    if (sname) {
#line 447 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 448 "../INSTALL/tstiee.f"
	    } else {
#line 449 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 450 "../INSTALL/tstiee.f"
	    }
#line 451 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 452 "../INSTALL/tstiee.f"
	    if (sname) {
#line 453 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 454 "../INSTALL/tstiee.f"
	    } else {
#line 455 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 456 "../INSTALL/tstiee.f"
	    }
#line 457 "../INSTALL/tstiee.f"
	}
#line 458 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 459 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 460 "../INSTALL/tstiee.f"
	    if (sname) {
#line 461 "../INSTALL/tstiee.f"
		nbmin = 8;
#line 462 "../INSTALL/tstiee.f"
	    } else {
#line 463 "../INSTALL/tstiee.f"
		nbmin = 8;
#line 464 "../INSTALL/tstiee.f"
	    }
#line 465 "../INSTALL/tstiee.f"
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 466 "../INSTALL/tstiee.f"
	    nbmin = 2;
#line 467 "../INSTALL/tstiee.f"
	}
#line 468 "../INSTALL/tstiee.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 469 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 470 "../INSTALL/tstiee.f"
	    nbmin = 2;
#line 471 "../INSTALL/tstiee.f"
	}
#line 472 "../INSTALL/tstiee.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 473 "../INSTALL/tstiee.f"
	if (*(unsigned char *)c3 == 'G') {
#line 474 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 477 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 478 "../INSTALL/tstiee.f"
	    }
#line 479 "../INSTALL/tstiee.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 480 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 483 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 484 "../INSTALL/tstiee.f"
	    }
#line 485 "../INSTALL/tstiee.f"
	}
#line 486 "../INSTALL/tstiee.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 487 "../INSTALL/tstiee.f"
	if (*(unsigned char *)c3 == 'G') {
#line 488 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 491 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 492 "../INSTALL/tstiee.f"
	    }
#line 493 "../INSTALL/tstiee.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 494 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 497 "../INSTALL/tstiee.f"
		nbmin = 2;
#line 498 "../INSTALL/tstiee.f"
	    }
#line 499 "../INSTALL/tstiee.f"
	}
#line 500 "../INSTALL/tstiee.f"
    }
#line 501 "../INSTALL/tstiee.f"
    ret_val = nbmin;
#line 502 "../INSTALL/tstiee.f"
    return ret_val;

#line 504 "../INSTALL/tstiee.f"
L300:

/*     ISPEC = 3:  crossover point */

#line 508 "../INSTALL/tstiee.f"
    nx = 0;
#line 509 "../INSTALL/tstiee.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 510 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
#line 512 "../INSTALL/tstiee.f"
	    if (sname) {
#line 513 "../INSTALL/tstiee.f"
		nx = 128;
#line 514 "../INSTALL/tstiee.f"
	    } else {
#line 515 "../INSTALL/tstiee.f"
		nx = 128;
#line 516 "../INSTALL/tstiee.f"
	    }
#line 517 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 518 "../INSTALL/tstiee.f"
	    if (sname) {
#line 519 "../INSTALL/tstiee.f"
		nx = 128;
#line 520 "../INSTALL/tstiee.f"
	    } else {
#line 521 "../INSTALL/tstiee.f"
		nx = 128;
#line 522 "../INSTALL/tstiee.f"
	    }
#line 523 "../INSTALL/tstiee.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 524 "../INSTALL/tstiee.f"
	    if (sname) {
#line 525 "../INSTALL/tstiee.f"
		nx = 128;
#line 526 "../INSTALL/tstiee.f"
	    } else {
#line 527 "../INSTALL/tstiee.f"
		nx = 128;
#line 528 "../INSTALL/tstiee.f"
	    }
#line 529 "../INSTALL/tstiee.f"
	}
#line 530 "../INSTALL/tstiee.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 531 "../INSTALL/tstiee.f"
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 532 "../INSTALL/tstiee.f"
	    nx = 32;
#line 533 "../INSTALL/tstiee.f"
	}
#line 534 "../INSTALL/tstiee.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 535 "../INSTALL/tstiee.f"
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 536 "../INSTALL/tstiee.f"
	    nx = 32;
#line 537 "../INSTALL/tstiee.f"
	}
#line 538 "../INSTALL/tstiee.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 539 "../INSTALL/tstiee.f"
	if (*(unsigned char *)c3 == 'G') {
#line 540 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 543 "../INSTALL/tstiee.f"
		nx = 128;
#line 544 "../INSTALL/tstiee.f"
	    }
#line 545 "../INSTALL/tstiee.f"
	}
#line 546 "../INSTALL/tstiee.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 547 "../INSTALL/tstiee.f"
	if (*(unsigned char *)c3 == 'G') {
#line 548 "../INSTALL/tstiee.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 551 "../INSTALL/tstiee.f"
		nx = 128;
#line 552 "../INSTALL/tstiee.f"
	    }
#line 553 "../INSTALL/tstiee.f"
	}
#line 554 "../INSTALL/tstiee.f"
    }
#line 555 "../INSTALL/tstiee.f"
    ret_val = nx;
#line 556 "../INSTALL/tstiee.f"
    return ret_val;

#line 558 "../INSTALL/tstiee.f"
L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

#line 562 "../INSTALL/tstiee.f"
    ret_val = 6;
#line 563 "../INSTALL/tstiee.f"
    return ret_val;

#line 565 "../INSTALL/tstiee.f"
L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

#line 569 "../INSTALL/tstiee.f"
    ret_val = 2;
#line 570 "../INSTALL/tstiee.f"
    return ret_val;

#line 572 "../INSTALL/tstiee.f"
L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

#line 576 "../INSTALL/tstiee.f"
    ret_val = (integer) ((doublereal) min(*n1,*n2) * 1.6);
#line 577 "../INSTALL/tstiee.f"
    return ret_val;

#line 579 "../INSTALL/tstiee.f"
L700:

/*     ISPEC = 7:  number of processors (not used) */

#line 583 "../INSTALL/tstiee.f"
    ret_val = 1;
#line 584 "../INSTALL/tstiee.f"
    return ret_val;

#line 586 "../INSTALL/tstiee.f"
L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

#line 590 "../INSTALL/tstiee.f"
    ret_val = 50;
#line 591 "../INSTALL/tstiee.f"
    return ret_val;

#line 593 "../INSTALL/tstiee.f"
L900:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
/*                 computation tree in the divide-and-conquer algorithm */
/*                 (used by xGELSD and xGESDD) */

#line 599 "../INSTALL/tstiee.f"
    ret_val = 25;
#line 600 "../INSTALL/tstiee.f"
    return ret_val;

#line 602 "../INSTALL/tstiee.f"
L1000:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */

#line 606 "../INSTALL/tstiee.f"
    ret_val = 1;
#line 607 "../INSTALL/tstiee.f"
    if (ret_val == 1) {
#line 608 "../INSTALL/tstiee.f"
	ret_val = ieeeck_(&c__0, &c_b227, &c_b228);
#line 609 "../INSTALL/tstiee.f"
    }
#line 610 "../INSTALL/tstiee.f"
    return ret_val;

#line 612 "../INSTALL/tstiee.f"
L1100:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap */

#line 616 "../INSTALL/tstiee.f"
    ret_val = 1;
#line 617 "../INSTALL/tstiee.f"
    if (ret_val == 1) {
#line 618 "../INSTALL/tstiee.f"
	ret_val = ieeeck_(&c__1, &c_b227, &c_b228);
#line 619 "../INSTALL/tstiee.f"
    }
#line 620 "../INSTALL/tstiee.f"
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

integer ieeeck_(integer *ispec, doublereal *zero, doublereal *one)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static doublereal nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, 
	    negzro, newzro;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  IEEECK is called from the ILAENV to verify that Inifinity and */
/*  possibly NaN arithmetic is safe (i.e. will not trap). */

/*  Arguments: */
/*  ========== */

/*  ISPEC   (input) INTEGER */
/*          Specifies whether to test just for inifinity arithmetic */
/*          or whether to test for infinity and NaN arithmetic. */
/*          = 0: Verify infinity arithmetic only. */
/*          = 1: Verify infinity and NaN arithmetic. */

/*  ZERO    (input) REAL */
/*          Must contain the value 0.0 */
/*          This is passed to prevent the compiler from optimizing */
/*          away this code. */

/*  ONE     (input) REAL */
/*          Must contain the value 1.0 */
/*          This is passed to prevent the compiler from optimizing */
/*          away this code. */

/*  RETURN VALUE:  INTEGER */
/*          = 0:  Arithmetic failed to produce the correct answers */
/*          = 1:  Arithmetic produced the correct answers */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
#line 670 "../INSTALL/tstiee.f"
    ret_val = 1;
#line 672 "../INSTALL/tstiee.f"
    posinf = *one / *zero;
#line 673 "../INSTALL/tstiee.f"
    if (posinf <= *one) {
#line 674 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 675 "../INSTALL/tstiee.f"
	return ret_val;
#line 676 "../INSTALL/tstiee.f"
    }
#line 678 "../INSTALL/tstiee.f"
    neginf = -(*one) / *zero;
#line 679 "../INSTALL/tstiee.f"
    if (neginf >= *zero) {
#line 680 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 681 "../INSTALL/tstiee.f"
	return ret_val;
#line 682 "../INSTALL/tstiee.f"
    }
#line 684 "../INSTALL/tstiee.f"
    negzro = *one / (neginf + *one);
#line 685 "../INSTALL/tstiee.f"
    if (negzro != *zero) {
#line 686 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 687 "../INSTALL/tstiee.f"
	return ret_val;
#line 688 "../INSTALL/tstiee.f"
    }
#line 690 "../INSTALL/tstiee.f"
    neginf = *one / negzro;
#line 691 "../INSTALL/tstiee.f"
    if (neginf >= *zero) {
#line 692 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 693 "../INSTALL/tstiee.f"
	return ret_val;
#line 694 "../INSTALL/tstiee.f"
    }
#line 696 "../INSTALL/tstiee.f"
    newzro = negzro + *zero;
#line 697 "../INSTALL/tstiee.f"
    if (newzro != *zero) {
#line 698 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 699 "../INSTALL/tstiee.f"
	return ret_val;
#line 700 "../INSTALL/tstiee.f"
    }
#line 702 "../INSTALL/tstiee.f"
    posinf = *one / newzro;
#line 703 "../INSTALL/tstiee.f"
    if (posinf <= *one) {
#line 704 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 705 "../INSTALL/tstiee.f"
	return ret_val;
#line 706 "../INSTALL/tstiee.f"
    }
#line 708 "../INSTALL/tstiee.f"
    neginf *= posinf;
#line 709 "../INSTALL/tstiee.f"
    if (neginf >= *zero) {
#line 710 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 711 "../INSTALL/tstiee.f"
	return ret_val;
#line 712 "../INSTALL/tstiee.f"
    }
#line 714 "../INSTALL/tstiee.f"
    posinf *= posinf;
#line 715 "../INSTALL/tstiee.f"
    if (posinf <= *one) {
#line 716 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 717 "../INSTALL/tstiee.f"
	return ret_val;
#line 718 "../INSTALL/tstiee.f"
    }

/*     Return if we were only asked to check infinity arithmetic */

#line 725 "../INSTALL/tstiee.f"
    if (*ispec == 0) {
#line 725 "../INSTALL/tstiee.f"
	return ret_val;
#line 725 "../INSTALL/tstiee.f"
    }
#line 727 "../INSTALL/tstiee.f"
    nan1 = posinf + neginf;
#line 729 "../INSTALL/tstiee.f"
    nan2 = posinf / neginf;
#line 731 "../INSTALL/tstiee.f"
    nan3 = posinf / posinf;
#line 733 "../INSTALL/tstiee.f"
    nan4 = posinf * *zero;
#line 735 "../INSTALL/tstiee.f"
    nan5 = neginf * negzro;
#line 737 "../INSTALL/tstiee.f"
    nan6 = nan5 * 0.;
#line 739 "../INSTALL/tstiee.f"
    if (nan1 == nan1) {
#line 740 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 741 "../INSTALL/tstiee.f"
	return ret_val;
#line 742 "../INSTALL/tstiee.f"
    }
#line 744 "../INSTALL/tstiee.f"
    if (nan2 == nan2) {
#line 745 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 746 "../INSTALL/tstiee.f"
	return ret_val;
#line 747 "../INSTALL/tstiee.f"
    }
#line 749 "../INSTALL/tstiee.f"
    if (nan3 == nan3) {
#line 750 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 751 "../INSTALL/tstiee.f"
	return ret_val;
#line 752 "../INSTALL/tstiee.f"
    }
#line 754 "../INSTALL/tstiee.f"
    if (nan4 == nan4) {
#line 755 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 756 "../INSTALL/tstiee.f"
	return ret_val;
#line 757 "../INSTALL/tstiee.f"
    }
#line 759 "../INSTALL/tstiee.f"
    if (nan5 == nan5) {
#line 760 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 761 "../INSTALL/tstiee.f"
	return ret_val;
#line 762 "../INSTALL/tstiee.f"
    }
#line 764 "../INSTALL/tstiee.f"
    if (nan6 == nan6) {
#line 765 "../INSTALL/tstiee.f"
	ret_val = 0;
#line 766 "../INSTALL/tstiee.f"
	return ret_val;
#line 767 "../INSTALL/tstiee.f"
    }
#line 769 "../INSTALL/tstiee.f"
    return ret_val;
} /* ieeeck_ */

/* Main program alias */ int tstiee_ () { MAIN__ (); return 0; }
