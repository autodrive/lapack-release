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
static doublereal c_b173 = 0.;
static doublereal c_b174 = 1.;
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
/* >               xHSEQR or related subroutines, */
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

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

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
    extern integer iparam2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen);
    static logical cname;
    static integer nbmin;
    static logical sname;
    extern integer ieeeck_(integer *, doublereal *, doublereal *);
    static char subnam[6];
    extern integer iparmq_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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
	case 17:  goto L170;
#line 191 "ilaenv.f"
	case 18:  goto L170;
#line 191 "ilaenv.f"
	case 19:  goto L170;
#line 191 "ilaenv.f"
	case 20:  goto L170;
#line 191 "ilaenv.f"
	case 21:  goto L170;
#line 191 "ilaenv.f"
    }

/*     Invalid value for ISPEC */

#line 197 "ilaenv.f"
    ret_val = -1;
#line 198 "ilaenv.f"
    return ret_val;

#line 200 "ilaenv.f"
L10:

/*     Convert NAME to upper case if the first character is lower case. */

#line 204 "ilaenv.f"
    ret_val = 1;
#line 205 "ilaenv.f"
    s_copy(subnam, name__, (ftnlen)6, name_len);
#line 206 "ilaenv.f"
    ic = *(unsigned char *)subnam;
#line 207 "ilaenv.f"
    iz = 'Z';
#line 208 "ilaenv.f"
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

#line 212 "ilaenv.f"
	if (ic >= 97 && ic <= 122) {
#line 213 "ilaenv.f"
	    *(unsigned char *)subnam = (char) (ic - 32);
#line 214 "ilaenv.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 215 "ilaenv.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 216 "ilaenv.f"
		if (ic >= 97 && ic <= 122) {
#line 216 "ilaenv.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 216 "ilaenv.f"
		}
#line 218 "ilaenv.f"
/* L20: */
#line 218 "ilaenv.f"
	    }
#line 219 "ilaenv.f"
	}

#line 221 "ilaenv.f"
    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

#line 225 "ilaenv.f"
	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
#line 228 "ilaenv.f"
	    *(unsigned char *)subnam = (char) (ic + 64);
#line 229 "ilaenv.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 230 "ilaenv.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 231 "ilaenv.f"
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
#line 231 "ilaenv.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
#line 231 "ilaenv.f"
		}
#line 235 "ilaenv.f"
/* L30: */
#line 235 "ilaenv.f"
	    }
#line 236 "ilaenv.f"
	}

#line 238 "ilaenv.f"
    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

#line 242 "ilaenv.f"
	if (ic >= 225 && ic <= 250) {
#line 243 "ilaenv.f"
	    *(unsigned char *)subnam = (char) (ic - 32);
#line 244 "ilaenv.f"
	    for (i__ = 2; i__ <= 6; ++i__) {
#line 245 "ilaenv.f"
		ic = *(unsigned char *)&subnam[i__ - 1];
#line 246 "ilaenv.f"
		if (ic >= 225 && ic <= 250) {
#line 246 "ilaenv.f"
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 246 "ilaenv.f"
		}
#line 248 "ilaenv.f"
/* L40: */
#line 248 "ilaenv.f"
	    }
#line 249 "ilaenv.f"
	}
#line 250 "ilaenv.f"
    }

#line 252 "ilaenv.f"
    *(unsigned char *)c1 = *(unsigned char *)subnam;
#line 253 "ilaenv.f"
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
#line 254 "ilaenv.f"
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
#line 255 "ilaenv.f"
    if (! (cname || sname)) {
#line 255 "ilaenv.f"
	return ret_val;
#line 255 "ilaenv.f"
    }
#line 257 "ilaenv.f"
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
#line 258 "ilaenv.f"
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
#line 259 "ilaenv.f"
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

#line 261 "ilaenv.f"
    switch (*ispec) {
#line 261 "ilaenv.f"
	case 1:  goto L50;
#line 261 "ilaenv.f"
	case 2:  goto L60;
#line 261 "ilaenv.f"
	case 3:  goto L70;
#line 261 "ilaenv.f"
    }

#line 263 "ilaenv.f"
L50:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

#line 271 "ilaenv.f"
    nb = 1;

#line 273 "ilaenv.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 274 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 275 "ilaenv.f"
	    if (sname) {
#line 276 "ilaenv.f"
		nb = 64;
#line 277 "ilaenv.f"
	    } else {
#line 278 "ilaenv.f"
		nb = 64;
#line 279 "ilaenv.f"
	    }
#line 280 "ilaenv.f"
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
#line 282 "ilaenv.f"
	    if (sname) {
#line 283 "ilaenv.f"
		nb = 32;
#line 284 "ilaenv.f"
	    } else {
#line 285 "ilaenv.f"
		nb = 32;
#line 286 "ilaenv.f"
	    }
#line 287 "ilaenv.f"
	} else if (s_cmp(c3, "QR ", (ftnlen)3, (ftnlen)3) == 0) {
#line 288 "ilaenv.f"
	    if (*n3 == 1) {
#line 289 "ilaenv.f"
		if (sname) {
/*     M*N */
#line 291 "ilaenv.f"
		    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
#line 292 "ilaenv.f"
			nb = *n1;
#line 293 "ilaenv.f"
		    } else {
#line 294 "ilaenv.f"
			nb = 32768 / *n2;
#line 295 "ilaenv.f"
		    }
#line 296 "ilaenv.f"
		} else {
#line 297 "ilaenv.f"
		    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
#line 298 "ilaenv.f"
			nb = *n1;
#line 299 "ilaenv.f"
		    } else {
#line 300 "ilaenv.f"
			nb = 32768 / *n2;
#line 301 "ilaenv.f"
		    }
#line 302 "ilaenv.f"
		}
#line 303 "ilaenv.f"
	    } else {
#line 304 "ilaenv.f"
		if (sname) {
#line 305 "ilaenv.f"
		    nb = 1;
#line 306 "ilaenv.f"
		} else {
#line 307 "ilaenv.f"
		    nb = 1;
#line 308 "ilaenv.f"
		}
#line 309 "ilaenv.f"
	    }
#line 310 "ilaenv.f"
	} else if (s_cmp(c3, "LQ ", (ftnlen)3, (ftnlen)3) == 0) {
#line 311 "ilaenv.f"
	    if (*n3 == 2) {
#line 312 "ilaenv.f"
		if (sname) {
/*     M*N */
#line 314 "ilaenv.f"
		    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
#line 315 "ilaenv.f"
			nb = *n1;
#line 316 "ilaenv.f"
		    } else {
#line 317 "ilaenv.f"
			nb = 32768 / *n2;
#line 318 "ilaenv.f"
		    }
#line 319 "ilaenv.f"
		} else {
#line 320 "ilaenv.f"
		    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
#line 321 "ilaenv.f"
			nb = *n1;
#line 322 "ilaenv.f"
		    } else {
#line 323 "ilaenv.f"
			nb = 32768 / *n2;
#line 324 "ilaenv.f"
		    }
#line 325 "ilaenv.f"
		}
#line 326 "ilaenv.f"
	    } else {
#line 327 "ilaenv.f"
		if (sname) {
#line 328 "ilaenv.f"
		    nb = 1;
#line 329 "ilaenv.f"
		} else {
#line 330 "ilaenv.f"
		    nb = 1;
#line 331 "ilaenv.f"
		}
#line 332 "ilaenv.f"
	    }
#line 333 "ilaenv.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 334 "ilaenv.f"
	    if (sname) {
#line 335 "ilaenv.f"
		nb = 32;
#line 336 "ilaenv.f"
	    } else {
#line 337 "ilaenv.f"
		nb = 32;
#line 338 "ilaenv.f"
	    }
#line 339 "ilaenv.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 340 "ilaenv.f"
	    if (sname) {
#line 341 "ilaenv.f"
		nb = 32;
#line 342 "ilaenv.f"
	    } else {
#line 343 "ilaenv.f"
		nb = 32;
#line 344 "ilaenv.f"
	    }
#line 345 "ilaenv.f"
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 346 "ilaenv.f"
	    if (sname) {
#line 347 "ilaenv.f"
		nb = 64;
#line 348 "ilaenv.f"
	    } else {
#line 349 "ilaenv.f"
		nb = 64;
#line 350 "ilaenv.f"
	    }
#line 351 "ilaenv.f"
	}
#line 352 "ilaenv.f"
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
#line 353 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 354 "ilaenv.f"
	    if (sname) {
#line 355 "ilaenv.f"
		nb = 64;
#line 356 "ilaenv.f"
	    } else {
#line 357 "ilaenv.f"
		nb = 64;
#line 358 "ilaenv.f"
	    }
#line 359 "ilaenv.f"
	}
#line 360 "ilaenv.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 361 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 362 "ilaenv.f"
	    if (sname) {
#line 363 "ilaenv.f"
		nb = 64;
#line 364 "ilaenv.f"
	    } else {
#line 365 "ilaenv.f"
		nb = 64;
#line 366 "ilaenv.f"
	    }
#line 367 "ilaenv.f"
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 368 "ilaenv.f"
	    nb = 32;
#line 369 "ilaenv.f"
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
#line 370 "ilaenv.f"
	    nb = 64;
#line 371 "ilaenv.f"
	}
#line 372 "ilaenv.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 373 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 374 "ilaenv.f"
	    nb = 64;
#line 375 "ilaenv.f"
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 376 "ilaenv.f"
	    nb = 32;
#line 377 "ilaenv.f"
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
#line 378 "ilaenv.f"
	    nb = 64;
#line 379 "ilaenv.f"
	}
#line 380 "ilaenv.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 381 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 382 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 385 "ilaenv.f"
		nb = 32;
#line 386 "ilaenv.f"
	    }
#line 387 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 388 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 391 "ilaenv.f"
		nb = 32;
#line 392 "ilaenv.f"
	    }
#line 393 "ilaenv.f"
	}
#line 394 "ilaenv.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 395 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 396 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 399 "ilaenv.f"
		nb = 32;
#line 400 "ilaenv.f"
	    }
#line 401 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 402 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 405 "ilaenv.f"
		nb = 32;
#line 406 "ilaenv.f"
	    }
#line 407 "ilaenv.f"
	}
#line 408 "ilaenv.f"
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
#line 409 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 410 "ilaenv.f"
	    if (sname) {
#line 411 "ilaenv.f"
		if (*n4 <= 64) {
#line 412 "ilaenv.f"
		    nb = 1;
#line 413 "ilaenv.f"
		} else {
#line 414 "ilaenv.f"
		    nb = 32;
#line 415 "ilaenv.f"
		}
#line 416 "ilaenv.f"
	    } else {
#line 417 "ilaenv.f"
		if (*n4 <= 64) {
#line 418 "ilaenv.f"
		    nb = 1;
#line 419 "ilaenv.f"
		} else {
#line 420 "ilaenv.f"
		    nb = 32;
#line 421 "ilaenv.f"
		}
#line 422 "ilaenv.f"
	    }
#line 423 "ilaenv.f"
	}
#line 424 "ilaenv.f"
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
#line 425 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 426 "ilaenv.f"
	    if (sname) {
#line 427 "ilaenv.f"
		if (*n2 <= 64) {
#line 428 "ilaenv.f"
		    nb = 1;
#line 429 "ilaenv.f"
		} else {
#line 430 "ilaenv.f"
		    nb = 32;
#line 431 "ilaenv.f"
		}
#line 432 "ilaenv.f"
	    } else {
#line 433 "ilaenv.f"
		if (*n2 <= 64) {
#line 434 "ilaenv.f"
		    nb = 1;
#line 435 "ilaenv.f"
		} else {
#line 436 "ilaenv.f"
		    nb = 32;
#line 437 "ilaenv.f"
		}
#line 438 "ilaenv.f"
	    }
#line 439 "ilaenv.f"
	}
#line 440 "ilaenv.f"
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
#line 441 "ilaenv.f"
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 442 "ilaenv.f"
	    if (sname) {
#line 443 "ilaenv.f"
		nb = 64;
#line 444 "ilaenv.f"
	    } else {
#line 445 "ilaenv.f"
		nb = 64;
#line 446 "ilaenv.f"
	    }
#line 447 "ilaenv.f"
	} else if (s_cmp(c3, "EVC", (ftnlen)3, (ftnlen)3) == 0) {
#line 448 "ilaenv.f"
	    if (sname) {
#line 449 "ilaenv.f"
		nb = 64;
#line 450 "ilaenv.f"
	    } else {
#line 451 "ilaenv.f"
		nb = 64;
#line 452 "ilaenv.f"
	    }
#line 453 "ilaenv.f"
	}
#line 454 "ilaenv.f"
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
#line 455 "ilaenv.f"
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
#line 456 "ilaenv.f"
	    if (sname) {
#line 457 "ilaenv.f"
		nb = 64;
#line 458 "ilaenv.f"
	    } else {
#line 459 "ilaenv.f"
		nb = 64;
#line 460 "ilaenv.f"
	    }
#line 461 "ilaenv.f"
	}
#line 462 "ilaenv.f"
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
#line 463 "ilaenv.f"
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
#line 464 "ilaenv.f"
	    nb = 1;
#line 465 "ilaenv.f"
	}
#line 466 "ilaenv.f"
    } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
#line 467 "ilaenv.f"
	nb = 32;
#line 468 "ilaenv.f"
	if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
#line 469 "ilaenv.f"
	    if (sname) {
#line 470 "ilaenv.f"
		nb = 32;
#line 471 "ilaenv.f"
	    } else {
#line 472 "ilaenv.f"
		nb = 32;
#line 473 "ilaenv.f"
	    }
#line 474 "ilaenv.f"
	}
#line 475 "ilaenv.f"
    }
#line 476 "ilaenv.f"
    ret_val = nb;
#line 477 "ilaenv.f"
    return ret_val;

#line 479 "ilaenv.f"
L60:

/*     ISPEC = 2:  minimum block size */

#line 483 "ilaenv.f"
    nbmin = 2;
#line 484 "ilaenv.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 485 "ilaenv.f"
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
#line 487 "ilaenv.f"
	    if (sname) {
#line 488 "ilaenv.f"
		nbmin = 2;
#line 489 "ilaenv.f"
	    } else {
#line 490 "ilaenv.f"
		nbmin = 2;
#line 491 "ilaenv.f"
	    }
#line 492 "ilaenv.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 493 "ilaenv.f"
	    if (sname) {
#line 494 "ilaenv.f"
		nbmin = 2;
#line 495 "ilaenv.f"
	    } else {
#line 496 "ilaenv.f"
		nbmin = 2;
#line 497 "ilaenv.f"
	    }
#line 498 "ilaenv.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 499 "ilaenv.f"
	    if (sname) {
#line 500 "ilaenv.f"
		nbmin = 2;
#line 501 "ilaenv.f"
	    } else {
#line 502 "ilaenv.f"
		nbmin = 2;
#line 503 "ilaenv.f"
	    }
#line 504 "ilaenv.f"
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
#line 505 "ilaenv.f"
	    if (sname) {
#line 506 "ilaenv.f"
		nbmin = 2;
#line 507 "ilaenv.f"
	    } else {
#line 508 "ilaenv.f"
		nbmin = 2;
#line 509 "ilaenv.f"
	    }
#line 510 "ilaenv.f"
	}
#line 511 "ilaenv.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 512 "ilaenv.f"
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
#line 513 "ilaenv.f"
	    if (sname) {
#line 514 "ilaenv.f"
		nbmin = 8;
#line 515 "ilaenv.f"
	    } else {
#line 516 "ilaenv.f"
		nbmin = 8;
#line 517 "ilaenv.f"
	    }
#line 518 "ilaenv.f"
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 519 "ilaenv.f"
	    nbmin = 2;
#line 520 "ilaenv.f"
	}
#line 521 "ilaenv.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 522 "ilaenv.f"
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 523 "ilaenv.f"
	    nbmin = 2;
#line 524 "ilaenv.f"
	}
#line 525 "ilaenv.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 526 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 527 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 530 "ilaenv.f"
		nbmin = 2;
#line 531 "ilaenv.f"
	    }
#line 532 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 533 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 536 "ilaenv.f"
		nbmin = 2;
#line 537 "ilaenv.f"
	    }
#line 538 "ilaenv.f"
	}
#line 539 "ilaenv.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 540 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 541 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 544 "ilaenv.f"
		nbmin = 2;
#line 545 "ilaenv.f"
	    }
#line 546 "ilaenv.f"
	} else if (*(unsigned char *)c3 == 'M') {
#line 547 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 550 "ilaenv.f"
		nbmin = 2;
#line 551 "ilaenv.f"
	    }
#line 552 "ilaenv.f"
	}
#line 553 "ilaenv.f"
    } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
#line 554 "ilaenv.f"
	nbmin = 2;
#line 555 "ilaenv.f"
	if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
#line 556 "ilaenv.f"
	    nbmin = 2;
#line 557 "ilaenv.f"
	}
#line 558 "ilaenv.f"
    }
#line 559 "ilaenv.f"
    ret_val = nbmin;
#line 560 "ilaenv.f"
    return ret_val;

#line 562 "ilaenv.f"
L70:

/*     ISPEC = 3:  crossover point */

#line 566 "ilaenv.f"
    nx = 0;
#line 567 "ilaenv.f"
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
#line 568 "ilaenv.f"
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
#line 570 "ilaenv.f"
	    if (sname) {
#line 571 "ilaenv.f"
		nx = 128;
#line 572 "ilaenv.f"
	    } else {
#line 573 "ilaenv.f"
		nx = 128;
#line 574 "ilaenv.f"
	    }
#line 575 "ilaenv.f"
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 576 "ilaenv.f"
	    if (sname) {
#line 577 "ilaenv.f"
		nx = 128;
#line 578 "ilaenv.f"
	    } else {
#line 579 "ilaenv.f"
		nx = 128;
#line 580 "ilaenv.f"
	    }
#line 581 "ilaenv.f"
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 582 "ilaenv.f"
	    if (sname) {
#line 583 "ilaenv.f"
		nx = 128;
#line 584 "ilaenv.f"
	    } else {
#line 585 "ilaenv.f"
		nx = 128;
#line 586 "ilaenv.f"
	    }
#line 587 "ilaenv.f"
	}
#line 588 "ilaenv.f"
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
#line 589 "ilaenv.f"
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 590 "ilaenv.f"
	    nx = 32;
#line 591 "ilaenv.f"
	}
#line 592 "ilaenv.f"
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
#line 593 "ilaenv.f"
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
#line 594 "ilaenv.f"
	    nx = 32;
#line 595 "ilaenv.f"
	}
#line 596 "ilaenv.f"
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
#line 597 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 598 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 601 "ilaenv.f"
		nx = 128;
#line 602 "ilaenv.f"
	    }
#line 603 "ilaenv.f"
	}
#line 604 "ilaenv.f"
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
#line 605 "ilaenv.f"
	if (*(unsigned char *)c3 == 'G') {
#line 606 "ilaenv.f"
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
#line 609 "ilaenv.f"
		nx = 128;
#line 610 "ilaenv.f"
	    }
#line 611 "ilaenv.f"
	}
#line 612 "ilaenv.f"
    } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
#line 613 "ilaenv.f"
	nx = 128;
#line 614 "ilaenv.f"
	if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
#line 615 "ilaenv.f"
	    nx = 128;
#line 616 "ilaenv.f"
	}
#line 617 "ilaenv.f"
    }
#line 618 "ilaenv.f"
    ret_val = nx;
#line 619 "ilaenv.f"
    return ret_val;

#line 621 "ilaenv.f"
L80:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

#line 625 "ilaenv.f"
    ret_val = 6;
#line 626 "ilaenv.f"
    return ret_val;

#line 628 "ilaenv.f"
L90:

/*     ISPEC = 5:  minimum column dimension (not used) */

#line 632 "ilaenv.f"
    ret_val = 2;
#line 633 "ilaenv.f"
    return ret_val;

#line 635 "ilaenv.f"
L100:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

#line 639 "ilaenv.f"
    ret_val = (integer) ((doublereal) min(*n1,*n2) * 1.6);
#line 640 "ilaenv.f"
    return ret_val;

#line 642 "ilaenv.f"
L110:

/*     ISPEC = 7:  number of processors (not used) */

#line 646 "ilaenv.f"
    ret_val = 1;
#line 647 "ilaenv.f"
    return ret_val;

#line 649 "ilaenv.f"
L120:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

#line 653 "ilaenv.f"
    ret_val = 50;
#line 654 "ilaenv.f"
    return ret_val;

#line 656 "ilaenv.f"
L130:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
/*                 computation tree in the divide-and-conquer algorithm */
/*                 (used by xGELSD and xGESDD) */

#line 662 "ilaenv.f"
    ret_val = 25;
#line 663 "ilaenv.f"
    return ret_val;

#line 665 "ilaenv.f"
L140:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */

/*     ILAENV = 0 */
#line 670 "ilaenv.f"
    ret_val = 1;
#line 671 "ilaenv.f"
    if (ret_val == 1) {
#line 672 "ilaenv.f"
	ret_val = ieeeck_(&c__1, &c_b173, &c_b174);
#line 673 "ilaenv.f"
    }
#line 674 "ilaenv.f"
    return ret_val;

#line 676 "ilaenv.f"
L150:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap */

/*     ILAENV = 0 */
#line 681 "ilaenv.f"
    ret_val = 1;
#line 682 "ilaenv.f"
    if (ret_val == 1) {
#line 683 "ilaenv.f"
	ret_val = ieeeck_(&c__0, &c_b173, &c_b174);
#line 684 "ilaenv.f"
    }
#line 685 "ilaenv.f"
    return ret_val;

#line 687 "ilaenv.f"
L160:

/*     12 <= ISPEC <= 16: xHSEQR or related subroutines. */

#line 691 "ilaenv.f"
    ret_val = iparmq_(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len)
	    ;
#line 692 "ilaenv.f"
    return ret_val;

#line 694 "ilaenv.f"
L170:

/*     17 <= ISPEC <= 21: 2stage eigenvalues and SVD or related subroutines. */

#line 698 "ilaenv.f"
    ret_val = iparam2stage_(ispec, name__, opts, n1, n2, n3, n4, name_len, 
	    opts_len);
#line 699 "ilaenv.f"
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

