#line 1 "iparmq.f"
/* iparmq.f -- translated by f2c (version 20100827).
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

#line 1 "iparmq.f"
/* > \brief \b IPARMQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IPARMQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, ISPEC, LWORK, N */
/*       CHARACTER          NAME*( * ), OPTS*( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >      This program sets problem and machine dependent parameters */
/* >      useful for xHSEQR and related subroutines for eigenvalue */
/* >      problems. It is called whenever */
/* >      IPARMQ is called with 12 <= ISPEC <= 16 */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is integer scalar */
/* >              ISPEC specifies which tunable parameter IPARMQ should */
/* >              return. */
/* > */
/* >              ISPEC=12: (INMIN)  Matrices of order nmin or less */
/* >                        are sent directly to xLAHQR, the implicit */
/* >                        double shift QR algorithm.  NMIN must be */
/* >                        at least 11. */
/* > */
/* >              ISPEC=13: (INWIN)  Size of the deflation window. */
/* >                        This is best set greater than or equal to */
/* >                        the number of simultaneous shifts NS. */
/* >                        Larger matrices benefit from larger deflation */
/* >                        windows. */
/* > */
/* >              ISPEC=14: (INIBL) Determines when to stop nibbling and */
/* >                        invest in an (expensive) multi-shift QR sweep. */
/* >                        If the aggressive early deflation subroutine */
/* >                        finds LD converged eigenvalues from an order */
/* >                        NW deflation window and LD.GT.(NW*NIBBLE)/100, */
/* >                        then the next QR sweep is skipped and early */
/* >                        deflation is applied immediately to the */
/* >                        remaining active diagonal block.  Setting */
/* >                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a */
/* >                        multi-shift QR sweep whenever early deflation */
/* >                        finds a converged eigenvalue.  Setting */
/* >                        IPARMQ(ISPEC=14) greater than or equal to 100 */
/* >                        prevents TTQRE from skipping a multi-shift */
/* >                        QR sweep. */
/* > */
/* >              ISPEC=15: (NSHFTS) The number of simultaneous shifts in */
/* >                        a multi-shift QR iteration. */
/* > */
/* >              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the */
/* >                        following meanings. */
/* >                        0:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are not */
/* >                            accumulated when updating the */
/* >                            far-from-diagonal matrix entries. */
/* >                        1:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are */
/* >                            accumulated, and matrix-matrix */
/* >                            multiplication is used to update the */
/* >                            far-from-diagonal matrix entries. */
/* >                        2:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are */
/* >                            accumulated, and 2-by-2 block structure */
/* >                            is exploited during matrix-matrix */
/* >                            multiplies. */
/* >                        (If xTRMM is slower than xGEMM, then */
/* >                        IPARMQ(ISPEC=16)=1 may be more efficient than */
/* >                        IPARMQ(ISPEC=16)=2 despite the greater level of */
/* >                        arithmetic work implied by the latter choice.) */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* >          NAME is character string */
/* >               Name of the calling subroutine */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* >          OPTS is character string */
/* >               This is a concatenation of the string arguments to */
/* >               TTQRE. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is integer scalar */
/* >               N is the order of the Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >               It is assumed that H is already upper triangular */
/* >               in rows and columns 1:ILO-1 and IHI+1:N. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is integer scalar */
/* >               The amount of workspace available. */
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
/* >       Little is known about how best to choose these parameters. */
/* >       It is possible to use different values of the parameters */
/* >       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR. */
/* > */
/* >       It is probably best to choose different parameters for */
/* >       different matrices and different parameters at different */
/* >       times during the iteration, but this has not been */
/* >       implemented --- yet. */
/* > */
/* > */
/* >       The best choices of most of the parameters depend */
/* >       in an ill-understood way on the relative execution */
/* >       rate of xLAQR3 and xLAQR5 and on the nature of each */
/* >       particular eigenvalue problem.  Experiment may be the */
/* >       only practical way to determine which choices are most */
/* >       effective. */
/* > */
/* >       Following is a list of default values supplied by IPARMQ. */
/* >       These defaults may be adjusted in order to attain better */
/* >       performance in any particular computational environment. */
/* > */
/* >       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point. */
/* >                        Default: 75. (Must be at least 11.) */
/* > */
/* >       IPARMQ(ISPEC=13) Recommended deflation window size. */
/* >                        This depends on ILO, IHI and NS, the */
/* >                        number of simultaneous shifts returned */
/* >                        by IPARMQ(ISPEC=15).  The default for */
/* >                        (IHI-ILO+1).LE.500 is NS.  The default */
/* >                        for (IHI-ILO+1).GT.500 is 3*NS/2. */
/* > */
/* >       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14. */
/* > */
/* >       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS. */
/* >                        a multi-shift QR iteration. */
/* > */
/* >                        If IHI-ILO+1 is ... */
/* > */
/* >                        greater than      ...but less    ... the */
/* >                        or equal to ...      than        default is */
/* > */
/* >                                0               30       NS =   2+ */
/* >                               30               60       NS =   4+ */
/* >                               60              150       NS =  10 */
/* >                              150              590       NS =  ** */
/* >                              590             3000       NS =  64 */
/* >                             3000             6000       NS = 128 */
/* >                             6000             infinity   NS = 256 */
/* > */
/* >                    (+)  By default matrices of this order are */
/* >                         passed to the implicit double shift routine */
/* >                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These */
/* >                         values of NS are used only in case of a rare */
/* >                         xLAHQR failure. */
/* > */
/* >                    (**) The asterisks (**) indicate an ad-hoc */
/* >                         function increasing from 10 to 64. */
/* > */
/* >       IPARMQ(ISPEC=16) Select structured matrix multiply. */
/* >                        (See ISPEC=16 above for details.) */
/* >                        Default: 3. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
integer iparmq_(integer *ispec, char *name__, char *opts, integer *n, integer 
	*ilo, integer *ihi, integer *lwork, ftnlen name_len, ftnlen opts_len)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);
    integer i_dnnt(doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, ic, nh, ns, iz;
    static char subnam[6];


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */

/*  ================================================================ */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 254 "iparmq.f"
    if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

/*        ==== Set the number simultaneous shifts ==== */

#line 259 "iparmq.f"
	nh = *ihi - *ilo + 1;
#line 260 "iparmq.f"
	ns = 2;
#line 261 "iparmq.f"
	if (nh >= 30) {
#line 261 "iparmq.f"
	    ns = 4;
#line 261 "iparmq.f"
	}
#line 263 "iparmq.f"
	if (nh >= 60) {
#line 263 "iparmq.f"
	    ns = 10;
#line 263 "iparmq.f"
	}
#line 265 "iparmq.f"
	if (nh >= 150) {
/* Computing MAX */
#line 265 "iparmq.f"
	    d__1 = log((doublereal) nh) / log(2.);
#line 265 "iparmq.f"
	    i__1 = 10, i__2 = nh / i_dnnt(&d__1);
#line 265 "iparmq.f"
	    ns = max(i__1,i__2);
#line 265 "iparmq.f"
	}
#line 267 "iparmq.f"
	if (nh >= 590) {
#line 267 "iparmq.f"
	    ns = 64;
#line 267 "iparmq.f"
	}
#line 269 "iparmq.f"
	if (nh >= 3000) {
#line 269 "iparmq.f"
	    ns = 128;
#line 269 "iparmq.f"
	}
#line 271 "iparmq.f"
	if (nh >= 6000) {
#line 271 "iparmq.f"
	    ns = 256;
#line 271 "iparmq.f"
	}
/* Computing MAX */
#line 273 "iparmq.f"
	i__1 = 2, i__2 = ns - ns % 2;
#line 273 "iparmq.f"
	ns = max(i__1,i__2);
#line 274 "iparmq.f"
    }

#line 276 "iparmq.f"
    if (*ispec == 12) {


/*        ===== Matrices of order smaller than NMIN get sent */
/*        .     to xLAHQR, the classic double shift algorithm. */
/*        .     This must be at least 11. ==== */

#line 283 "iparmq.f"
	ret_val = 75;

#line 285 "iparmq.f"
    } else if (*ispec == 14) {

/*        ==== INIBL: skip a multi-shift qr iteration and */
/*        .    whenever aggressive early deflation finds */
/*        .    at least (NIBBLE*(window size)/100) deflations. ==== */

#line 291 "iparmq.f"
	ret_val = 14;

#line 293 "iparmq.f"
    } else if (*ispec == 15) {

/*        ==== NSHFTS: The number of simultaneous shifts ===== */

#line 297 "iparmq.f"
	ret_val = ns;

#line 299 "iparmq.f"
    } else if (*ispec == 13) {

/*        ==== NW: deflation window size.  ==== */

#line 303 "iparmq.f"
	if (nh <= 500) {
#line 304 "iparmq.f"
	    ret_val = ns;
#line 305 "iparmq.f"
	} else {
#line 306 "iparmq.f"
	    ret_val = ns * 3 / 2;
#line 307 "iparmq.f"
	}

#line 309 "iparmq.f"
    } else if (*ispec == 16) {

/*        ==== IACC22: Whether to accumulate reflections */
/*        .     before updating the far-from-diagonal elements */
/*        .     and whether to use 2-by-2 block structure while */
/*        .     doing it.  A small amount of work could be saved */
/*        .     by making this choice dependent also upon the */
/*        .     NH=IHI-ILO+1. */


/*        Convert NAME to upper case if the first character is lower case. */

#line 321 "iparmq.f"
	ret_val = 0;
#line 322 "iparmq.f"
	s_copy(subnam, name__, (ftnlen)6, name_len);
#line 323 "iparmq.f"
	ic = *(unsigned char *)subnam;
#line 324 "iparmq.f"
	iz = 'Z';
#line 325 "iparmq.f"
	if (iz == 90 || iz == 122) {

/*           ASCII character set */

#line 329 "iparmq.f"
	    if (ic >= 97 && ic <= 122) {
#line 330 "iparmq.f"
		*(unsigned char *)subnam = (char) (ic - 32);
#line 331 "iparmq.f"
		for (i__ = 2; i__ <= 6; ++i__) {
#line 332 "iparmq.f"
		    ic = *(unsigned char *)&subnam[i__ - 1];
#line 333 "iparmq.f"
		    if (ic >= 97 && ic <= 122) {
#line 333 "iparmq.f"
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 333 "iparmq.f"
		    }
#line 335 "iparmq.f"
		}
#line 336 "iparmq.f"
	    }

#line 338 "iparmq.f"
	} else if (iz == 233 || iz == 169) {

/*           EBCDIC character set */

#line 342 "iparmq.f"
	    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 
		    && ic <= 169) {
#line 345 "iparmq.f"
		*(unsigned char *)subnam = (char) (ic + 64);
#line 346 "iparmq.f"
		for (i__ = 2; i__ <= 6; ++i__) {
#line 347 "iparmq.f"
		    ic = *(unsigned char *)&subnam[i__ - 1];
#line 348 "iparmq.f"
		    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || 
			    ic >= 162 && ic <= 169) {
#line 348 "iparmq.f"
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
#line 348 "iparmq.f"
		    }
#line 352 "iparmq.f"
		}
#line 353 "iparmq.f"
	    }

#line 355 "iparmq.f"
	} else if (iz == 218 || iz == 250) {

/*           Prime machines:  ASCII+128 */

#line 359 "iparmq.f"
	    if (ic >= 225 && ic <= 250) {
#line 360 "iparmq.f"
		*(unsigned char *)subnam = (char) (ic - 32);
#line 361 "iparmq.f"
		for (i__ = 2; i__ <= 6; ++i__) {
#line 362 "iparmq.f"
		    ic = *(unsigned char *)&subnam[i__ - 1];
#line 363 "iparmq.f"
		    if (ic >= 225 && ic <= 250) {
#line 363 "iparmq.f"
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 363 "iparmq.f"
		    }
#line 365 "iparmq.f"
		}
#line 366 "iparmq.f"
	    }
#line 367 "iparmq.f"
	}

#line 369 "iparmq.f"
	if (s_cmp(subnam + 1, "GGHRD", (ftnlen)5, (ftnlen)5) == 0 || s_cmp(
		subnam + 1, "GGHD3", (ftnlen)5, (ftnlen)5) == 0) {
#line 371 "iparmq.f"
	    ret_val = 1;
#line 372 "iparmq.f"
	    if (nh >= 14) {
#line 372 "iparmq.f"
		ret_val = 2;
#line 372 "iparmq.f"
	    }
#line 374 "iparmq.f"
	} else if (s_cmp(subnam + 3, "EXC", (ftnlen)3, (ftnlen)3) == 0) {
#line 375 "iparmq.f"
	    if (nh >= 14) {
#line 375 "iparmq.f"
		ret_val = 1;
#line 375 "iparmq.f"
	    }
#line 377 "iparmq.f"
	    if (nh >= 14) {
#line 377 "iparmq.f"
		ret_val = 2;
#line 377 "iparmq.f"
	    }
#line 379 "iparmq.f"
	} else if (s_cmp(subnam + 1, "HSEQR", (ftnlen)5, (ftnlen)5) == 0 || 
		s_cmp(subnam + 1, "LAQR", (ftnlen)4, (ftnlen)4) == 0) {
#line 381 "iparmq.f"
	    if (ns >= 14) {
#line 381 "iparmq.f"
		ret_val = 1;
#line 381 "iparmq.f"
	    }
#line 383 "iparmq.f"
	    if (ns >= 14) {
#line 383 "iparmq.f"
		ret_val = 2;
#line 383 "iparmq.f"
	    }
#line 385 "iparmq.f"
	}

#line 387 "iparmq.f"
    } else {
/*        ===== invalid value of ispec ===== */
#line 389 "iparmq.f"
	ret_val = -1;

#line 391 "iparmq.f"
    }

/*     ==== End of IPARMQ ==== */

#line 395 "iparmq.f"
    return ret_val;
} /* iparmq_ */

