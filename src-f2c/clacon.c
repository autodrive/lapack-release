#line 1 "clacon.f"
/* clacon.f -- translated by f2c (version 20100827).
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

#line 1 "clacon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLACON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLACON( N, V, X, EST, KASE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       REAL               EST */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            V( N ), X( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACON estimates the 1-norm of a square, complex matrix A. */
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
/* >         where A**H is the conjugate transpose of A, and CLACON must be */
/* >         re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* >          EST is REAL */
/* >         On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* >         unchanged from the previous call to CLACON. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to CLACON, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**H * X. */
/* >         On the final return from CLACON, KASE will again be 0. */
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
/* >  Originally named CONEST, dated March 16, 1988. \n */
/* >  Last modified:  April, 1999 */

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
/* Subroutine */ int clacon_(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, iter;
    static doublereal temp;
    static integer jump;
    static doublereal absxi;
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
/*     .. Save statement .. */
/*     .. */
/*     .. Executable Statements .. */

#line 161 "clacon.f"
    /* Parameter adjustments */
#line 161 "clacon.f"
    --x;
#line 161 "clacon.f"
    --v;
#line 161 "clacon.f"

#line 161 "clacon.f"
    /* Function Body */
#line 161 "clacon.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 162 "clacon.f"
    if (*kase == 0) {
#line 163 "clacon.f"
	i__1 = *n;
#line 163 "clacon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 164 "clacon.f"
	    i__2 = i__;
#line 164 "clacon.f"
	    d__1 = 1. / (doublereal) (*n);
#line 164 "clacon.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 164 "clacon.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 165 "clacon.f"
/* L10: */
#line 165 "clacon.f"
	}
#line 166 "clacon.f"
	*kase = 1;
#line 167 "clacon.f"
	jump = 1;
#line 168 "clacon.f"
	return 0;
#line 169 "clacon.f"
    }

#line 171 "clacon.f"
    switch (jump) {
#line 171 "clacon.f"
	case 1:  goto L20;
#line 171 "clacon.f"
	case 2:  goto L40;
#line 171 "clacon.f"
	case 3:  goto L70;
#line 171 "clacon.f"
	case 4:  goto L90;
#line 171 "clacon.f"
	case 5:  goto L120;
#line 171 "clacon.f"
    }

/*     ................ ENTRY   (JUMP = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 176 "clacon.f"
L20:
#line 177 "clacon.f"
    if (*n == 1) {
#line 178 "clacon.f"
	v[1].r = x[1].r, v[1].i = x[1].i;
#line 179 "clacon.f"
	*est = z_abs(&v[1]);
/*        ... QUIT */
#line 181 "clacon.f"
	goto L130;
#line 182 "clacon.f"
    }
#line 183 "clacon.f"
    *est = scsum1_(n, &x[1], &c__1);

#line 185 "clacon.f"
    i__1 = *n;
#line 185 "clacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 186 "clacon.f"
	absxi = z_abs(&x[i__]);
#line 187 "clacon.f"
	if (absxi > safmin) {
#line 188 "clacon.f"
	    i__2 = i__;
#line 188 "clacon.f"
	    i__3 = i__;
#line 188 "clacon.f"
	    d__1 = x[i__3].r / absxi;
#line 188 "clacon.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 188 "clacon.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 188 "clacon.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 190 "clacon.f"
	} else {
#line 191 "clacon.f"
	    i__2 = i__;
#line 191 "clacon.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 192 "clacon.f"
	}
#line 193 "clacon.f"
/* L30: */
#line 193 "clacon.f"
    }
#line 194 "clacon.f"
    *kase = 2;
#line 195 "clacon.f"
    jump = 2;
#line 196 "clacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 201 "clacon.f"
L40:
#line 202 "clacon.f"
    j = icmax1_(n, &x[1], &c__1);
#line 203 "clacon.f"
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 207 "clacon.f"
L50:
#line 208 "clacon.f"
    i__1 = *n;
#line 208 "clacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "clacon.f"
	i__2 = i__;
#line 209 "clacon.f"
	x[i__2].r = 0., x[i__2].i = 0.;
#line 210 "clacon.f"
/* L60: */
#line 210 "clacon.f"
    }
#line 211 "clacon.f"
    i__1 = j;
#line 211 "clacon.f"
    x[i__1].r = 1., x[i__1].i = 0.;
#line 212 "clacon.f"
    *kase = 1;
#line 213 "clacon.f"
    jump = 3;
#line 214 "clacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 219 "clacon.f"
L70:
#line 220 "clacon.f"
    ccopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 221 "clacon.f"
    estold = *est;
#line 222 "clacon.f"
    *est = scsum1_(n, &v[1], &c__1);

/*     TEST FOR CYCLING. */
#line 225 "clacon.f"
    if (*est <= estold) {
#line 225 "clacon.f"
	goto L100;
#line 225 "clacon.f"
    }

#line 228 "clacon.f"
    i__1 = *n;
#line 228 "clacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "clacon.f"
	absxi = z_abs(&x[i__]);
#line 230 "clacon.f"
	if (absxi > safmin) {
#line 231 "clacon.f"
	    i__2 = i__;
#line 231 "clacon.f"
	    i__3 = i__;
#line 231 "clacon.f"
	    d__1 = x[i__3].r / absxi;
#line 231 "clacon.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 231 "clacon.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 231 "clacon.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 233 "clacon.f"
	} else {
#line 234 "clacon.f"
	    i__2 = i__;
#line 234 "clacon.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 235 "clacon.f"
	}
#line 236 "clacon.f"
/* L80: */
#line 236 "clacon.f"
    }
#line 237 "clacon.f"
    *kase = 2;
#line 238 "clacon.f"
    jump = 4;
#line 239 "clacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 4) */
/*     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 244 "clacon.f"
L90:
#line 245 "clacon.f"
    jlast = j;
#line 246 "clacon.f"
    j = icmax1_(n, &x[1], &c__1);
#line 247 "clacon.f"
    if (z_abs(&x[jlast]) != z_abs(&x[j]) && iter < 5) {
#line 249 "clacon.f"
	++iter;
#line 250 "clacon.f"
	goto L50;
#line 251 "clacon.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 255 "clacon.f"
L100:
#line 256 "clacon.f"
    altsgn = 1.;
#line 257 "clacon.f"
    i__1 = *n;
#line 257 "clacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "clacon.f"
	i__2 = i__;
#line 258 "clacon.f"
	d__1 = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
#line 258 "clacon.f"
	z__1.r = d__1, z__1.i = 0.;
#line 258 "clacon.f"
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 259 "clacon.f"
	altsgn = -altsgn;
#line 260 "clacon.f"
/* L110: */
#line 260 "clacon.f"
    }
#line 261 "clacon.f"
    *kase = 1;
#line 262 "clacon.f"
    jump = 5;
#line 263 "clacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 268 "clacon.f"
L120:
#line 269 "clacon.f"
    temp = scsum1_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 270 "clacon.f"
    if (temp > *est) {
#line 271 "clacon.f"
	ccopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 272 "clacon.f"
	*est = temp;
#line 273 "clacon.f"
    }

#line 275 "clacon.f"
L130:
#line 276 "clacon.f"
    *kase = 0;
#line 277 "clacon.f"
    return 0;

/*     End of CLACON */

} /* clacon_ */

