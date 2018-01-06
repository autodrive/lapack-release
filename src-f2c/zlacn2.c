#line 1 "zlacn2.f"
/* zlacn2.f -- translated by f2c (version 20100827).
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

#line 1 "zlacn2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLACN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacn2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacn2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacn2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLACN2( N, V, X, EST, KASE, ISAVE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       DOUBLE PRECISION   EST */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISAVE( 3 ) */
/*       COMPLEX*16         V( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACN2 estimates the 1-norm of a square, complex matrix A. */
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
/* >          V is COMPLEX*16 array, dimension (N) */
/* >         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
/* >         (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (N) */
/* >         On an intermediate return, X should be overwritten by */
/* >               A * X,   if KASE=1, */
/* >               A**H * X,  if KASE=2, */
/* >         where A**H is the conjugate transpose of A, and ZLACN2 must be */
/* >         re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* >          EST is DOUBLE PRECISION */
/* >         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be */
/* >         unchanged from the previous call to ZLACN2. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to ZLACN2, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**H * X. */
/* >         On the final return from ZLACN2, KASE will again be 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ISAVE */
/* > \verbatim */
/* >          ISAVE is INTEGER array, dimension (3) */
/* >         ISAVE is used to save variables between calls to ZLACN2 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Originally named CONEST, dated March 16, 1988. */
/* > */
/* >  Last modified:  April, 1999 */
/* > */
/* >  This is a thread safe version of ZLACON, which uses the array ISAVE */
/* >  in place of a SAVE statement, as follows: */
/* > */
/* >     ZLACON     ZLACN2 */
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
/* Subroutine */ int zlacn2_(integer *n, doublecomplex *v, doublecomplex *x, 
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
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern integer izmax1_(integer *, doublecomplex *, integer *);
    extern doublereal dzsum1_(integer *, doublecomplex *, integer *), dlamch_(
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

#line 178 "zlacn2.f"
    /* Parameter adjustments */
#line 178 "zlacn2.f"
    --isave;
#line 178 "zlacn2.f"
    --x;
#line 178 "zlacn2.f"
    --v;
#line 178 "zlacn2.f"

#line 178 "zlacn2.f"
    /* Function Body */
#line 178 "zlacn2.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 179 "zlacn2.f"
    if (*kase == 0) {
#line 180 "zlacn2.f"
	i__1 = *n;
#line 180 "zlacn2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 181 "zlacn2.f"
	    i__2 = i__;
#line 181 "zlacn2.f"
	    d__1 = 1. / (doublereal) (*n);
#line 181 "zlacn2.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 181 "zlacn2.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 182 "zlacn2.f"
/* L10: */
#line 182 "zlacn2.f"
	}
#line 183 "zlacn2.f"
	*kase = 1;
#line 184 "zlacn2.f"
	isave[1] = 1;
#line 185 "zlacn2.f"
	return 0;
#line 186 "zlacn2.f"
    }

#line 188 "zlacn2.f"
    switch (isave[1]) {
#line 188 "zlacn2.f"
	case 1:  goto L20;
#line 188 "zlacn2.f"
	case 2:  goto L40;
#line 188 "zlacn2.f"
	case 3:  goto L70;
#line 188 "zlacn2.f"
	case 4:  goto L90;
#line 188 "zlacn2.f"
	case 5:  goto L120;
#line 188 "zlacn2.f"
    }

/*     ................ ENTRY   (ISAVE( 1 ) = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 193 "zlacn2.f"
L20:
#line 194 "zlacn2.f"
    if (*n == 1) {
#line 195 "zlacn2.f"
	v[1].r = x[1].r, v[1].i = x[1].i;
#line 196 "zlacn2.f"
	*est = z_abs(&v[1]);
/*        ... QUIT */
#line 198 "zlacn2.f"
	goto L130;
#line 199 "zlacn2.f"
    }
#line 200 "zlacn2.f"
    *est = dzsum1_(n, &x[1], &c__1);

#line 202 "zlacn2.f"
    i__1 = *n;
#line 202 "zlacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "zlacn2.f"
	absxi = z_abs(&x[i__]);
#line 204 "zlacn2.f"
	if (absxi > safmin) {
#line 205 "zlacn2.f"
	    i__2 = i__;
#line 205 "zlacn2.f"
	    i__3 = i__;
#line 205 "zlacn2.f"
	    d__1 = x[i__3].r / absxi;
#line 205 "zlacn2.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 205 "zlacn2.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 205 "zlacn2.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 207 "zlacn2.f"
	} else {
#line 208 "zlacn2.f"
	    i__2 = i__;
#line 208 "zlacn2.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 209 "zlacn2.f"
	}
#line 210 "zlacn2.f"
/* L30: */
#line 210 "zlacn2.f"
    }
#line 211 "zlacn2.f"
    *kase = 2;
#line 212 "zlacn2.f"
    isave[1] = 2;
#line 213 "zlacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 218 "zlacn2.f"
L40:
#line 219 "zlacn2.f"
    isave[2] = izmax1_(n, &x[1], &c__1);
#line 220 "zlacn2.f"
    isave[3] = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 224 "zlacn2.f"
L50:
#line 225 "zlacn2.f"
    i__1 = *n;
#line 225 "zlacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 226 "zlacn2.f"
	i__2 = i__;
#line 226 "zlacn2.f"
	x[i__2].r = 0., x[i__2].i = 0.;
#line 227 "zlacn2.f"
/* L60: */
#line 227 "zlacn2.f"
    }
#line 228 "zlacn2.f"
    i__1 = isave[2];
#line 228 "zlacn2.f"
    x[i__1].r = 1., x[i__1].i = 0.;
#line 229 "zlacn2.f"
    *kase = 1;
#line 230 "zlacn2.f"
    isave[1] = 3;
#line 231 "zlacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 236 "zlacn2.f"
L70:
#line 237 "zlacn2.f"
    zcopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 238 "zlacn2.f"
    estold = *est;
#line 239 "zlacn2.f"
    *est = dzsum1_(n, &v[1], &c__1);

/*     TEST FOR CYCLING. */
#line 242 "zlacn2.f"
    if (*est <= estold) {
#line 242 "zlacn2.f"
	goto L100;
#line 242 "zlacn2.f"
    }

#line 245 "zlacn2.f"
    i__1 = *n;
#line 245 "zlacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "zlacn2.f"
	absxi = z_abs(&x[i__]);
#line 247 "zlacn2.f"
	if (absxi > safmin) {
#line 248 "zlacn2.f"
	    i__2 = i__;
#line 248 "zlacn2.f"
	    i__3 = i__;
#line 248 "zlacn2.f"
	    d__1 = x[i__3].r / absxi;
#line 248 "zlacn2.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 248 "zlacn2.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 248 "zlacn2.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 250 "zlacn2.f"
	} else {
#line 251 "zlacn2.f"
	    i__2 = i__;
#line 251 "zlacn2.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 252 "zlacn2.f"
	}
#line 253 "zlacn2.f"
/* L80: */
#line 253 "zlacn2.f"
    }
#line 254 "zlacn2.f"
    *kase = 2;
#line 255 "zlacn2.f"
    isave[1] = 4;
#line 256 "zlacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 4) */
/*     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 261 "zlacn2.f"
L90:
#line 262 "zlacn2.f"
    jlast = isave[2];
#line 263 "zlacn2.f"
    isave[2] = izmax1_(n, &x[1], &c__1);
#line 264 "zlacn2.f"
    if (z_abs(&x[jlast]) != z_abs(&x[isave[2]]) && isave[3] < 5) {
#line 266 "zlacn2.f"
	++isave[3];
#line 267 "zlacn2.f"
	goto L50;
#line 268 "zlacn2.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 272 "zlacn2.f"
L100:
#line 273 "zlacn2.f"
    altsgn = 1.;
#line 274 "zlacn2.f"
    i__1 = *n;
#line 274 "zlacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "zlacn2.f"
	i__2 = i__;
#line 275 "zlacn2.f"
	d__1 = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
#line 275 "zlacn2.f"
	z__1.r = d__1, z__1.i = 0.;
#line 275 "zlacn2.f"
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 276 "zlacn2.f"
	altsgn = -altsgn;
#line 277 "zlacn2.f"
/* L110: */
#line 277 "zlacn2.f"
    }
#line 278 "zlacn2.f"
    *kase = 1;
#line 279 "zlacn2.f"
    isave[1] = 5;
#line 280 "zlacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 285 "zlacn2.f"
L120:
#line 286 "zlacn2.f"
    temp = dzsum1_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 287 "zlacn2.f"
    if (temp > *est) {
#line 288 "zlacn2.f"
	zcopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 289 "zlacn2.f"
	*est = temp;
#line 290 "zlacn2.f"
    }

#line 292 "zlacn2.f"
L130:
#line 293 "zlacn2.f"
    *kase = 0;
#line 294 "zlacn2.f"
    return 0;

/*     End of ZLACN2 */

} /* zlacn2_ */

