#line 1 "clacn2.f"
/* clacn2.f -- translated by f2c (version 20100827).
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

#line 1 "clacn2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLACN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacn2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacn2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacn2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       REAL               EST */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISAVE( 3 ) */
/*       COMPLEX            V( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACN2 estimates the 1-norm of a square, complex matrix A. */
/* > Reverse communication is used for evaluating matrix-vector products. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The order of the matrix.  N >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (N) */
/* >         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
/* >         (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (N) */
/* >         On an intermediate return, X should be overwritten by */
/* >               A * X,   if KASE=1, */
/* >               A**H * X,  if KASE=2, */
/* >         where A**H is the conjugate transpose of A, and CLACN2 must be */
/* >         re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* >          EST is REAL */
/* >         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be */
/* >         unchanged from the previous call to CLACN2. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to CLACN2, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**H * X. */
/* >         On the final return from CLACN2, KASE will again be 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ISAVE */
/* > \verbatim */
/* >          ISAVE is INTEGER array, dimension (3) */
/* >         ISAVE is used to save variables between calls to SLACN2 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Originally named CONEST, dated March 16, 1988. */
/* > */
/* >  Last modified:  April, 1999 */
/* > */
/* >  This is a thread safe version of CLACON, which uses the array ISAVE */
/* >  in place of a SAVE statement, as follows: */
/* > */
/* >     CLACON     CLACN2 */
/* >      JUMP     ISAVE(1) */
/* >      J        ISAVE(2) */
/* >      ITER     ISAVE(3) */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Nick Higham, University of Manchester */

/* > \par References: */
/*  ================ */
/* > */
/* >  N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/* >  a real or complex matrix, with applications to condition estimation", */
/* >  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clacn2_(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase, integer *isave)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublereal temp, absxi;
    static integer jlast;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern integer icmax1_(integer *, doublecomplex *, integer *);
    extern doublereal scsum1_(integer *, doublecomplex *, integer *), slamch_(
	    char *, ftnlen);
    static doublereal safmin, altsgn, estold;


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 178 "clacn2.f"
    /* Parameter adjustments */
#line 178 "clacn2.f"
    --isave;
#line 178 "clacn2.f"
    --x;
#line 178 "clacn2.f"
    --v;
#line 178 "clacn2.f"

#line 178 "clacn2.f"
    /* Function Body */
#line 178 "clacn2.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 179 "clacn2.f"
    if (*kase == 0) {
#line 180 "clacn2.f"
	i__1 = *n;
#line 180 "clacn2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 181 "clacn2.f"
	    i__2 = i__;
#line 181 "clacn2.f"
	    d__1 = 1. / (doublereal) (*n);
#line 181 "clacn2.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 181 "clacn2.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 182 "clacn2.f"
/* L10: */
#line 182 "clacn2.f"
	}
#line 183 "clacn2.f"
	*kase = 1;
#line 184 "clacn2.f"
	isave[1] = 1;
#line 185 "clacn2.f"
	return 0;
#line 186 "clacn2.f"
    }

#line 188 "clacn2.f"
    switch (isave[1]) {
#line 188 "clacn2.f"
	case 1:  goto L20;
#line 188 "clacn2.f"
	case 2:  goto L40;
#line 188 "clacn2.f"
	case 3:  goto L70;
#line 188 "clacn2.f"
	case 4:  goto L90;
#line 188 "clacn2.f"
	case 5:  goto L120;
#line 188 "clacn2.f"
    }

/*     ................ ENTRY   (ISAVE( 1 ) = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 193 "clacn2.f"
L20:
#line 194 "clacn2.f"
    if (*n == 1) {
#line 195 "clacn2.f"
	v[1].r = x[1].r, v[1].i = x[1].i;
#line 196 "clacn2.f"
	*est = z_abs(&v[1]);
/*        ... QUIT */
#line 198 "clacn2.f"
	goto L130;
#line 199 "clacn2.f"
    }
#line 200 "clacn2.f"
    *est = scsum1_(n, &x[1], &c__1);

#line 202 "clacn2.f"
    i__1 = *n;
#line 202 "clacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "clacn2.f"
	absxi = z_abs(&x[i__]);
#line 204 "clacn2.f"
	if (absxi > safmin) {
#line 205 "clacn2.f"
	    i__2 = i__;
#line 205 "clacn2.f"
	    i__3 = i__;
#line 205 "clacn2.f"
	    d__1 = x[i__3].r / absxi;
#line 205 "clacn2.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 205 "clacn2.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 205 "clacn2.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 207 "clacn2.f"
	} else {
#line 208 "clacn2.f"
	    i__2 = i__;
#line 208 "clacn2.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 209 "clacn2.f"
	}
#line 210 "clacn2.f"
/* L30: */
#line 210 "clacn2.f"
    }
#line 211 "clacn2.f"
    *kase = 2;
#line 212 "clacn2.f"
    isave[1] = 2;
#line 213 "clacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 218 "clacn2.f"
L40:
#line 219 "clacn2.f"
    isave[2] = icmax1_(n, &x[1], &c__1);
#line 220 "clacn2.f"
    isave[3] = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 224 "clacn2.f"
L50:
#line 225 "clacn2.f"
    i__1 = *n;
#line 225 "clacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 226 "clacn2.f"
	i__2 = i__;
#line 226 "clacn2.f"
	x[i__2].r = 0., x[i__2].i = 0.;
#line 227 "clacn2.f"
/* L60: */
#line 227 "clacn2.f"
    }
#line 228 "clacn2.f"
    i__1 = isave[2];
#line 228 "clacn2.f"
    x[i__1].r = 1., x[i__1].i = 0.;
#line 229 "clacn2.f"
    *kase = 1;
#line 230 "clacn2.f"
    isave[1] = 3;
#line 231 "clacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 236 "clacn2.f"
L70:
#line 237 "clacn2.f"
    ccopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 238 "clacn2.f"
    estold = *est;
#line 239 "clacn2.f"
    *est = scsum1_(n, &v[1], &c__1);

/*     TEST FOR CYCLING. */
#line 242 "clacn2.f"
    if (*est <= estold) {
#line 242 "clacn2.f"
	goto L100;
#line 242 "clacn2.f"
    }

#line 245 "clacn2.f"
    i__1 = *n;
#line 245 "clacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "clacn2.f"
	absxi = z_abs(&x[i__]);
#line 247 "clacn2.f"
	if (absxi > safmin) {
#line 248 "clacn2.f"
	    i__2 = i__;
#line 248 "clacn2.f"
	    i__3 = i__;
#line 248 "clacn2.f"
	    d__1 = x[i__3].r / absxi;
#line 248 "clacn2.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 248 "clacn2.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 248 "clacn2.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 250 "clacn2.f"
	} else {
#line 251 "clacn2.f"
	    i__2 = i__;
#line 251 "clacn2.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 252 "clacn2.f"
	}
#line 253 "clacn2.f"
/* L80: */
#line 253 "clacn2.f"
    }
#line 254 "clacn2.f"
    *kase = 2;
#line 255 "clacn2.f"
    isave[1] = 4;
#line 256 "clacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 4) */
/*     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 261 "clacn2.f"
L90:
#line 262 "clacn2.f"
    jlast = isave[2];
#line 263 "clacn2.f"
    isave[2] = icmax1_(n, &x[1], &c__1);
#line 264 "clacn2.f"
    if (z_abs(&x[jlast]) != z_abs(&x[isave[2]]) && isave[3] < 5) {
#line 266 "clacn2.f"
	++isave[3];
#line 267 "clacn2.f"
	goto L50;
#line 268 "clacn2.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 272 "clacn2.f"
L100:
#line 273 "clacn2.f"
    altsgn = 1.;
#line 274 "clacn2.f"
    i__1 = *n;
#line 274 "clacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "clacn2.f"
	i__2 = i__;
#line 275 "clacn2.f"
	d__1 = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
#line 275 "clacn2.f"
	z__1.r = d__1, z__1.i = 0.;
#line 275 "clacn2.f"
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 276 "clacn2.f"
	altsgn = -altsgn;
#line 277 "clacn2.f"
/* L110: */
#line 277 "clacn2.f"
    }
#line 278 "clacn2.f"
    *kase = 1;
#line 279 "clacn2.f"
    isave[1] = 5;
#line 280 "clacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 285 "clacn2.f"
L120:
#line 286 "clacn2.f"
    temp = scsum1_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 287 "clacn2.f"
    if (temp > *est) {
#line 288 "clacn2.f"
	ccopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 289 "clacn2.f"
	*est = temp;
#line 290 "clacn2.f"
    }

#line 292 "clacn2.f"
L130:
#line 293 "clacn2.f"
    *kase = 0;
#line 294 "clacn2.f"
    return 0;

/*     End of CLACN2 */

} /* clacn2_ */

