#line 1 "ilaenv.f"
/* ilaenv.f -- translated by f2c (version 20100827).
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

#line 1 "ilaenv.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b163 = 0.;
static doublereal c_b164 = 1.;
static integer c__0 = 0;

/* > \brief \b ILAENV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILAENV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*( * )    NAME, OPTS */
/*       INTEGER            ISPEC, N1, N2, N3, N4 */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ILAENV is called from the LAPACK routines to choose problem-dependent */
/* > parameters for the local environment.  See ISPEC for a description of */
/* > the parameters. */
/* > */
/* > ILAENV returns an INTEGER */
/* > if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC */
/* > if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value. */
/* > */
/* > This version provides a set of parameters which should give good, */
/* > but not optimal, performance on many of the currently available */
/* > computers.  Users are encouraged to modify this subroutine to set */
/* > the tuning parameters for their particular machine using the option */
/* > and problem size information in the arguments. */
/* > */
/* > This routine will not function correctly if it is converted to all */
/* > lower case.  Converting it to all upper case is allowed. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is INTEGER */
/* >          Specifies the parameter to be returned as the value of */
/* >          ILAENV. */
/* >          = 1: the optimal blocksize; if this value is 1, an unblocked */
/* >               algorithm will give the best performance. */
/* >          = 2: the minimum block size for which the block routine */
/* >               should be used; if the usable block size is less than */
/* >               this value, an unblocked routine should be used. */
/* >          = 3: the crossover point (in a block routine, for N less */
/* >               than this value, an unblocked routine should be used) */
/* >          = 4: the number of shifts, used in the nonsymmetric */
/* >               eigenvalue routines (DEPRECATED) */
/* >          = 5: the minimum column dimension for blocking to be used; */
/* >               rectangular blocks must have dimension at least k by m, */
/* >               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
/* >          = 6: the crossover point for the SVD (when reducing an m by n */
/* >               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
/* >               this value, a QR factorization is used first to reduce */
/* >               the matrix to a triangular form.) */
/* >          = 7: the number of processors */
/* >          = 8: the crossover point for the multishift QR method */
/* >               for nonsymmetric eigenvalue problems (DEPRECATED) */
/* >          = 9: maximum size of the subproblems at the bottom of the */
/* >               computation tree in the divide-and-conquer algorithm */
/* >               (used by xGELSD and xGESDD) */
/* >          =10: ieee NaN arithmetic can be trusted not to trap */
/* >          =11: infinity arithmetic can be trusted not to trap */
/* >          12 <= ISPEC <= 16: */
/* >               xHSEQR or one of its subroutines, */
/* >               see IPARMQ for detailed explanation */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* >          NAME is CHARACTER*(*) */
/* >          The name of the calling subroutine, in either upper case or */
/* >          lower case. */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* >          OPTS is CHARACTER*(*) */
/* >          The character options to the subroutine NAME, concatenated */
/* >          into a single character string.  For example, UPLO = 'U', */
/* >          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/* >          be specified as OPTS = 'UTN'. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N2 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N3 */
/* > \verbatim */
/* >          N3 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N4 */
/* > \verbatim */
/* >          N4 is INTEGER */
/* >          Problem dimensions for the subroutine NAME; these may not all */
/* >          be required. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The following conventions have been used when calling ILAENV from the */
/* >  LAPACK routines: */
/* >  1)  OPTS is a concatenation of all of the character options to */
/* >      subroutine NAME, in the same order that they appear in the */
/* >      argument list for NAME, even if they are not used in determining */
/* >      the value of the parameter specified by ISPEC. */
/* >  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
/* >      that they appear in the argument list for NAME.  N1 is used */
/* >      first, N2 second, and so on, and unused problem dimensions are */
/* >      passed a value of -1. */
/* >  3)  The parameter value returned by ILAENV is checked for validity in */
/* >      the calling subroutine.  For example, ILAENV is used to retrieve */
/* >      the optimal blocksize for STRTRI as follows: */
/* > */
/* >      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/* >      IF( NB.LE.1 ) NB = MAX( 1, N ) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
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
    static logical cname;
    static integer nbmin;
    static logical sname;
    extern integer ieeeck_(integer *, doublereal *, doublereal *);
    static char subnam[6];
    extern integer iparmq_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 191 "ilaenv.f"
    switch (*ispec) {
#line 191 "ilaenv.f"
	case 1:  goto L10;
#line 191 "ilaenv.f"
	case 2:  goto L10;
#line 191 "ilaenv.f"
	case 3:  goto L10;
#line 191 "ilaenv.f"
	case 4:  goto L80;
#line 191 "ilaenv.f"
	case 5:  goto L90;
#line 191 "ilaenv.f"
	case 6:  goto L100;
#line 191 "ilaenv.f"
	case 7:  goto L110;
#line 191 "ilaenv.f"
	case 8:  goto L120;
#line 191 "ilaenv.f"
	case 9:  goto L130;
#line 191 "ilaenv.f"
	case 10:  goto L140;
#line 191 "ilaenv.f"
	case 11:  goto L150;
#line 191 "ilaenv.f"
	case 12:  goto L160;
#line 191 "ilaenv.f"
	case 13:  goto L160;
#line 191 "ilaenv.f"
	case 14:  goto L160;
#line 191 "ilaenv.f"
	case 15:  goto L160;
#line 191 "ilaenv.f"
	case 16:  goto L160;
#line 191 "ilaenv.f"
    }

/*     Invalid value for ISPEC */

#line 196 "ilaenv.f"
    ret_val = -1;
#line 197 "ilaenv.f"
    return ret_val;

#line 199 "ilaenv.f"
L10:

/*     Convert NAME to upper case if the first character is lower case. */

#line 203 "ilaenv.f"
    ret_val = 1;
#line 204 "ilaenv.f"
    s_copy(subnam, name__, (ftnlen)6, name_len);
#line 205 "ilaenv.f"
    ic = *(unsigned char *)subnam;
#line 206 "ilaenv.f"
    iz = 'Z';
#line 207 "ilaenv.f"
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

#line 211 "ilaenv.f"
	if (ic >= 97 && ic <= 122) {
#line 212 "ilaenv.f"
	    *(unsigned char *)subnam = (char) (ic - 32);
#line 213 "ilaenv.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 214 "ilaenv.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 215 "ilaenv.f"
		if (ic >= 97 && ic <= 122) {
#line 215 "ilaenv.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 215 "ilaenv.f"
		}
#line 217 "ilaenv.f"
/* L20: */
#line 217 "ilaenv.f"
	    }
#line 218 "ilaenv.f"
	}

#line 220 "ilaenv.f"
    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

#line 224 "ilaenv.f"
	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
#line 227 "ilaenv.f"
	    *(unsigned char *)subnam = (char) (ic + 64);
#line 228 "ilaenv.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 229 "ilaenv.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 230 "ilaenv.f"
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
#line 230 "ilaenv.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
#line 230 "ilaenv.f"
		}
#line 234 "ilaenv.f"
/* L30: */
#line 234 "ilaenv.f"
	    }
#line 235 "ilaenv.f"
	}

#line 237 "ilaenv.f"
    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

#line 241 "ilaenv.f"
	if (ic >= 225 && ic <= 250) {
#line 242 "ilaenv.f"
	    *(unsigned char *)subnam = (char) (ic - 32);
#line 243 "ilaenv.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 244 "ilaenv.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 245 "ilaenv.f"
		if (ic >= 225 && ic <= 250) {
#line 245 "ilaenv.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 245 "ilaenv.f"
		}
#line 247 "ilaenv.f"
/* L40: */
#line 247 "ilaenv.f"
	    }
#line 248 "ilaenv.f"
	}
#line 249 "ilaenv.f"
    }

#line 251 "ilaenv.f"
    *(unsigned char *)c1 = *(unsigned char *)subnam;
#line 252 "ilaenv.f"
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
#line 253 "ilaenv.f"
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
#line 254 "ilaenv.f"
    if (! (cname || sname)) {
#line 254 "ilaenv.f"
	return ret_val;
#line 254 "ilaenv.f"
    }
#line 256 "ilaenv.f"
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
#line 257 "ilaenv.f"
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
#line 258 "ilaenv.f"
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

#line 260 "ilaenv.f"
    switch (*ispec) {
#line 260 "ilaenv.f"
	case 1:  goto L50;
#line 260 "ilaenv.f"
	case 2:  goto L60;
#line 260 "ilaenv.f"
	case 3:  goto L70;
#line 260 "ilaenv.f"
    }

#line 262 "ilaenv.f"
L50:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

#line 270 "ilaenv.f"
    nb = 1;

#line 272 "ilaenv.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 273 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 274 "ilaenv.f"
	    if (sname) {
#line 275 "ilaenv.f"
		nb = 64;
#line 276 "ilaenv.f"
	    } else {
#line 277 "ilaenv.f"
		nb = 64;
#line 278 "ilaenv.f"
	    }
#line 279 "ilaenv.f"
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
#line 281 "ilaenv.f"
	    if (sname) {
#line 282 "ilaenv.f"
		nb = 32;
#line 283 "ilaenv.f"
	    } else {
#line 284 "ilaenv.f"
		nb = 32;
#line 285 "ilaenv.f"
	    }
#line 286 "ilaenv.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 287 "ilaenv.f"
	    if (sname) {
#line 288 "ilaenv.f"
		nb = 32;
#line 289 "ilaenv.f"
	    } else {
#line 290 "ilaenv.f"
		nb = 32;
#line 291 "ilaenv.f"
	    }
#line 292 "ilaenv.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 293 "ilaenv.f"
	    if (sname) {
#line 294 "ilaenv.f"
		nb = 32;
#line 295 "ilaenv.f"
	    } else {
#line 296 "ilaenv.f"
		nb = 32;
#line 297 "ilaenv.f"
	    }
#line 298 "ilaenv.f"
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 299 "ilaenv.f"
	    if (sname) {
#line 300 "ilaenv.f"
		nb = 64;
#line 301 "ilaenv.f"
	    } else {
#line 302 "ilaenv.f"
		nb = 64;
#line 303 "ilaenv.f"
	    }
#line 304 "ilaenv.f"
	}
#line 305 "ilaenv.f"
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
#line 306 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 307 "ilaenv.f"
	    if (sname) {
#line 308 "ilaenv.f"
		nb = 64;
#line 309 "ilaenv.f"
	    } else {
#line 310 "ilaenv.f"
		nb = 64;
#line 311 "ilaenv.f"
	    }
#line 312 "ilaenv.f"
	}
#line 313 "ilaenv.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 314 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 315 "ilaenv.f"
	    if (sname) {
#line 316 "ilaenv.f"
		nb = 64;
#line 317 "ilaenv.f"
	    } else {
#line 318 "ilaenv.f"
		nb = 64;
#line 319 "ilaenv.f"
	    }
#line 320 "ilaenv.f"
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 321 "ilaenv.f"
	    nb = 32;
#line 322 "ilaenv.f"
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
#line 323 "ilaenv.f"
	    nb = 64;
#line 324 "ilaenv.f"
	}
#line 325 "ilaenv.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 326 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 327 "ilaenv.f"
	    nb = 64;
#line 328 "ilaenv.f"
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 329 "ilaenv.f"
	    nb = 32;
#line 330 "ilaenv.f"
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
#line 331 "ilaenv.f"
	    nb = 64;
#line 332 "ilaenv.f"
	}
#line 333 "ilaenv.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 334 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 335 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 338 "ilaenv.f"
		nb = 32;
#line 339 "ilaenv.f"
	    }
#line 340 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 341 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 344 "ilaenv.f"
		nb = 32;
#line 345 "ilaenv.f"
	    }
#line 346 "ilaenv.f"
	}
#line 347 "ilaenv.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 348 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 349 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 352 "ilaenv.f"
		nb = 32;
#line 353 "ilaenv.f"
	    }
#line 354 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 355 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 358 "ilaenv.f"
		nb = 32;
#line 359 "ilaenv.f"
	    }
#line 360 "ilaenv.f"
	}
#line 361 "ilaenv.f"
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
#line 362 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 363 "ilaenv.f"
	    if (sname) {
#line 364 "ilaenv.f"
		if (*n4 <= 64) {
#line 365 "ilaenv.f"
		    nb = 1;
#line 366 "ilaenv.f"
		} else {
#line 367 "ilaenv.f"
		    nb = 32;
#line 368 "ilaenv.f"
		}
#line 369 "ilaenv.f"
	    } else {
#line 370 "ilaenv.f"
		if (*n4 <= 64) {
#line 371 "ilaenv.f"
		    nb = 1;
#line 372 "ilaenv.f"
		} else {
#line 373 "ilaenv.f"
		    nb = 32;
#line 374 "ilaenv.f"
		}
#line 375 "ilaenv.f"
	    }
#line 376 "ilaenv.f"
	}
#line 377 "ilaenv.f"
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
#line 378 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 379 "ilaenv.f"
	    if (sname) {
#line 380 "ilaenv.f"
		if (*n2 <= 64) {
#line 381 "ilaenv.f"
		    nb = 1;
#line 382 "ilaenv.f"
		} else {
#line 383 "ilaenv.f"
		    nb = 32;
#line 384 "ilaenv.f"
		}
#line 385 "ilaenv.f"
	    } else {
#line 386 "ilaenv.f"
		if (*n2 <= 64) {
#line 387 "ilaenv.f"
		    nb = 1;
#line 388 "ilaenv.f"
		} else {
#line 389 "ilaenv.f"
		    nb = 32;
#line 390 "ilaenv.f"
		}
#line 391 "ilaenv.f"
	    }
#line 392 "ilaenv.f"
	}
#line 393 "ilaenv.f"
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
#line 394 "ilaenv.f"
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 395 "ilaenv.f"
	    if (sname) {
#line 396 "ilaenv.f"
		nb = 64;
#line 397 "ilaenv.f"
	    } else {
#line 398 "ilaenv.f"
		nb = 64;
#line 399 "ilaenv.f"
	    }
#line 400 "ilaenv.f"
	}
#line 401 "ilaenv.f"
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
#line 402 "ilaenv.f"
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
#line 403 "ilaenv.f"
	    if (sname) {
#line 404 "ilaenv.f"
		nb = 64;
#line 405 "ilaenv.f"
	    } else {
#line 406 "ilaenv.f"
		nb = 64;
#line 407 "ilaenv.f"
	    }
#line 408 "ilaenv.f"
	}
#line 409 "ilaenv.f"
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
#line 410 "ilaenv.f"
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
#line 411 "ilaenv.f"
	    nb = 1;
#line 412 "ilaenv.f"
	}
#line 413 "ilaenv.f"
    }
#line 414 "ilaenv.f"
    ret_val = nb;
#line 415 "ilaenv.f"
    return ret_val;

#line 417 "ilaenv.f"
L60:

/*     ISPEC = 2:  minimum block size */

#line 421 "ilaenv.f"
    nbmin = 2;
#line 422 "ilaenv.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 423 "ilaenv.f"
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
#line 425 "ilaenv.f"
	    if (sname) {
#line 426 "ilaenv.f"
		nbmin = 2;
#line 427 "ilaenv.f"
	    } else {
#line 428 "ilaenv.f"
		nbmin = 2;
#line 429 "ilaenv.f"
	    }
#line 430 "ilaenv.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 431 "ilaenv.f"
	    if (sname) {
#line 432 "ilaenv.f"
		nbmin = 2;
#line 433 "ilaenv.f"
	    } else {
#line 434 "ilaenv.f"
		nbmin = 2;
#line 435 "ilaenv.f"
	    }
#line 436 "ilaenv.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 437 "ilaenv.f"
	    if (sname) {
#line 438 "ilaenv.f"
		nbmin = 2;
#line 439 "ilaenv.f"
	    } else {
#line 440 "ilaenv.f"
		nbmin = 2;
#line 441 "ilaenv.f"
	    }
#line 442 "ilaenv.f"
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 443 "ilaenv.f"
	    if (sname) {
#line 444 "ilaenv.f"
		nbmin = 2;
#line 445 "ilaenv.f"
	    } else {
#line 446 "ilaenv.f"
		nbmin = 2;
#line 447 "ilaenv.f"
	    }
#line 448 "ilaenv.f"
	}
#line 449 "ilaenv.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 450 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 451 "ilaenv.f"
	    if (sname) {
#line 452 "ilaenv.f"
		nbmin = 8;
#line 453 "ilaenv.f"
	    } else {
#line 454 "ilaenv.f"
		nbmin = 8;
#line 455 "ilaenv.f"
	    }
#line 456 "ilaenv.f"
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 457 "ilaenv.f"
	    nbmin = 2;
#line 458 "ilaenv.f"
	}
#line 459 "ilaenv.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 460 "ilaenv.f"
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 461 "ilaenv.f"
	    nbmin = 2;
#line 462 "ilaenv.f"
	}
#line 463 "ilaenv.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 464 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 465 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 468 "ilaenv.f"
		nbmin = 2;
#line 469 "ilaenv.f"
	    }
#line 470 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 471 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 474 "ilaenv.f"
		nbmin = 2;
#line 475 "ilaenv.f"
	    }
#line 476 "ilaenv.f"
	}
#line 477 "ilaenv.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 478 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 479 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 482 "ilaenv.f"
		nbmin = 2;
#line 483 "ilaenv.f"
	    }
#line 484 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 485 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 488 "ilaenv.f"
		nbmin = 2;
#line 489 "ilaenv.f"
	    }
#line 490 "ilaenv.f"
	}
#line 491 "ilaenv.f"
    }
#line 492 "ilaenv.f"
    ret_val = nbmin;
#line 493 "ilaenv.f"
    return ret_val;

#line 495 "ilaenv.f"
L70:

/*     ISPEC = 3:  crossover point */

#line 499 "ilaenv.f"
    nx = 0;
#line 500 "ilaenv.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 501 "ilaenv.f"
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
#line 503 "ilaenv.f"
	    if (sname) {
#line 504 "ilaenv.f"
		nx = 128;
#line 505 "ilaenv.f"
	    } else {
#line 506 "ilaenv.f"
		nx = 128;
#line 507 "ilaenv.f"
	    }
#line 508 "ilaenv.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 509 "ilaenv.f"
	    if (sname) {
#line 510 "ilaenv.f"
		nx = 128;
#line 511 "ilaenv.f"
	    } else {
#line 512 "ilaenv.f"
		nx = 128;
#line 513 "ilaenv.f"
	    }
#line 514 "ilaenv.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 515 "ilaenv.f"
	    if (sname) {
#line 516 "ilaenv.f"
		nx = 128;
#line 517 "ilaenv.f"
	    } else {
#line 518 "ilaenv.f"
		nx = 128;
#line 519 "ilaenv.f"
	    }
#line 520 "ilaenv.f"
	}
#line 521 "ilaenv.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 522 "ilaenv.f"
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 523 "ilaenv.f"
	    nx = 32;
#line 524 "ilaenv.f"
	}
#line 525 "ilaenv.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 526 "ilaenv.f"
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 527 "ilaenv.f"
	    nx = 32;
#line 528 "ilaenv.f"
	}
#line 529 "ilaenv.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 530 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 531 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 534 "ilaenv.f"
		nx = 128;
#line 535 "ilaenv.f"
	    }
#line 536 "ilaenv.f"
	}
#line 537 "ilaenv.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 538 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 539 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 542 "ilaenv.f"
		nx = 128;
#line 543 "ilaenv.f"
	    }
#line 544 "ilaenv.f"
	}
#line 545 "ilaenv.f"
    }
#line 546 "ilaenv.f"
    ret_val = nx;
#line 547 "ilaenv.f"
    return ret_val;

#line 549 "ilaenv.f"
L80:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

#line 553 "ilaenv.f"
    ret_val = 6;
#line 554 "ilaenv.f"
    return ret_val;

#line 556 "ilaenv.f"
L90:

/*     ISPEC = 5:  minimum column dimension (not used) */

#line 560 "ilaenv.f"
    ret_val = 2;
#line 561 "ilaenv.f"
    return ret_val;

#line 563 "ilaenv.f"
L100:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

#line 567 "ilaenv.f"
    ret_val = (integer) ((doublereal) min(*n1,*n2) * 1.6);
#line 568 "ilaenv.f"
    return ret_val;

#line 570 "ilaenv.f"
L110:

/*     ISPEC = 7:  number of processors (not used) */

#line 574 "ilaenv.f"
    ret_val = 1;
#line 575 "ilaenv.f"
    return ret_val;

#line 577 "ilaenv.f"
L120:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

#line 581 "ilaenv.f"
    ret_val = 50;
#line 582 "ilaenv.f"
    return ret_val;

#line 584 "ilaenv.f"
L130:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
/*                 computation tree in the divide-and-conquer algorithm */
/*                 (used by xGELSD and xGESDD) */

#line 590 "ilaenv.f"
    ret_val = 25;
#line 591 "ilaenv.f"
    return ret_val;

#line 593 "ilaenv.f"
L140:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */

/*     ILAENV = 0 */
#line 598 "ilaenv.f"
    ret_val = 1;
#line 599 "ilaenv.f"
    if (ret_val == 1) {
#line 600 "ilaenv.f"
	ret_val = ieeeck_(&c__1, &c_b163, &c_b164);
#line 601 "ilaenv.f"
    }
#line 602 "ilaenv.f"
    return ret_val;

#line 604 "ilaenv.f"
L150:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap */

/*     ILAENV = 0 */
#line 609 "ilaenv.f"
    ret_val = 1;
#line 610 "ilaenv.f"
    if (ret_val == 1) {
#line 611 "ilaenv.f"
	ret_val = ieeeck_(&c__0, &c_b163, &c_b164);
#line 612 "ilaenv.f"
    }
#line 613 "ilaenv.f"
    return ret_val;

#line 615 "ilaenv.f"
L160:

/*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. */

#line 619 "ilaenv.f"
    ret_val = iparmq_(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len)
	    ;
#line 620 "ilaenv.f"
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

