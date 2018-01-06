#line 1 "zlacon.f"
/* zlacon.f -- translated by f2c (version 20100827).
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

#line 1 "zlacon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLACON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLACON( N, V, X, EST, KASE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       DOUBLE PRECISION   EST */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         V( N ), X( N ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACON estimates the 1-norm of a square, complex matrix A. */
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
/* >         where A**H is the conjugate transpose of A, and ZLACON must be */
/* >         re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* >          EST is DOUBLE PRECISION */
/* >         On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* >         unchanged from the previous call to ZLACON. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to ZLACON, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**H * X. */
/* >         On the final return from ZLACON, KASE will again be 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

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
/* Subroutine */ int zlacon_(integer *n, doublecomplex *v, doublecomplex *x, 
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
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern integer izmax1_(integer *, doublecomplex *, integer *);
    extern doublereal dzsum1_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    static doublereal safmin, altsgn, estold;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 161 "zlacon.f"
    /* Parameter adjustments */
#line 161 "zlacon.f"
    --x;
#line 161 "zlacon.f"
    --v;
#line 161 "zlacon.f"

#line 161 "zlacon.f"
    /* Function Body */
#line 161 "zlacon.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 162 "zlacon.f"
    if (*kase == 0) {
#line 163 "zlacon.f"
	i__1 = *n;
#line 163 "zlacon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 164 "zlacon.f"
	    i__2 = i__;
#line 164 "zlacon.f"
	    d__1 = 1. / (doublereal) (*n);
#line 164 "zlacon.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 164 "zlacon.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 165 "zlacon.f"
/* L10: */
#line 165 "zlacon.f"
	}
#line 166 "zlacon.f"
	*kase = 1;
#line 167 "zlacon.f"
	jump = 1;
#line 168 "zlacon.f"
	return 0;
#line 169 "zlacon.f"
    }

#line 171 "zlacon.f"
    switch (jump) {
#line 171 "zlacon.f"
	case 1:  goto L20;
#line 171 "zlacon.f"
	case 2:  goto L40;
#line 171 "zlacon.f"
	case 3:  goto L70;
#line 171 "zlacon.f"
	case 4:  goto L90;
#line 171 "zlacon.f"
	case 5:  goto L120;
#line 171 "zlacon.f"
    }

/*     ................ ENTRY   (JUMP = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 176 "zlacon.f"
L20:
#line 177 "zlacon.f"
    if (*n == 1) {
#line 178 "zlacon.f"
	v[1].r = x[1].r, v[1].i = x[1].i;
#line 179 "zlacon.f"
	*est = z_abs(&v[1]);
/*        ... QUIT */
#line 181 "zlacon.f"
	goto L130;
#line 182 "zlacon.f"
    }
#line 183 "zlacon.f"
    *est = dzsum1_(n, &x[1], &c__1);

#line 185 "zlacon.f"
    i__1 = *n;
#line 185 "zlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 186 "zlacon.f"
	absxi = z_abs(&x[i__]);
#line 187 "zlacon.f"
	if (absxi > safmin) {
#line 188 "zlacon.f"
	    i__2 = i__;
#line 188 "zlacon.f"
	    i__3 = i__;
#line 188 "zlacon.f"
	    d__1 = x[i__3].r / absxi;
#line 188 "zlacon.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 188 "zlacon.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 188 "zlacon.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 190 "zlacon.f"
	} else {
#line 191 "zlacon.f"
	    i__2 = i__;
#line 191 "zlacon.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 192 "zlacon.f"
	}
#line 193 "zlacon.f"
/* L30: */
#line 193 "zlacon.f"
    }
#line 194 "zlacon.f"
    *kase = 2;
#line 195 "zlacon.f"
    jump = 2;
#line 196 "zlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 201 "zlacon.f"
L40:
#line 202 "zlacon.f"
    j = izmax1_(n, &x[1], &c__1);
#line 203 "zlacon.f"
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 207 "zlacon.f"
L50:
#line 208 "zlacon.f"
    i__1 = *n;
#line 208 "zlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "zlacon.f"
	i__2 = i__;
#line 209 "zlacon.f"
	x[i__2].r = 0., x[i__2].i = 0.;
#line 210 "zlacon.f"
/* L60: */
#line 210 "zlacon.f"
    }
#line 211 "zlacon.f"
    i__1 = j;
#line 211 "zlacon.f"
    x[i__1].r = 1., x[i__1].i = 0.;
#line 212 "zlacon.f"
    *kase = 1;
#line 213 "zlacon.f"
    jump = 3;
#line 214 "zlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 219 "zlacon.f"
L70:
#line 220 "zlacon.f"
    zcopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 221 "zlacon.f"
    estold = *est;
#line 222 "zlacon.f"
    *est = dzsum1_(n, &v[1], &c__1);

/*     TEST FOR CYCLING. */
#line 225 "zlacon.f"
    if (*est <= estold) {
#line 225 "zlacon.f"
	goto L100;
#line 225 "zlacon.f"
    }

#line 228 "zlacon.f"
    i__1 = *n;
#line 228 "zlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "zlacon.f"
	absxi = z_abs(&x[i__]);
#line 230 "zlacon.f"
	if (absxi > safmin) {
#line 231 "zlacon.f"
	    i__2 = i__;
#line 231 "zlacon.f"
	    i__3 = i__;
#line 231 "zlacon.f"
	    d__1 = x[i__3].r / absxi;
#line 231 "zlacon.f"
	    d__2 = d_imag(&x[i__]) / absxi;
#line 231 "zlacon.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 231 "zlacon.f"
	    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 233 "zlacon.f"
	} else {
#line 234 "zlacon.f"
	    i__2 = i__;
#line 234 "zlacon.f"
	    x[i__2].r = 1., x[i__2].i = 0.;
#line 235 "zlacon.f"
	}
#line 236 "zlacon.f"
/* L80: */
#line 236 "zlacon.f"
    }
#line 237 "zlacon.f"
    *kase = 2;
#line 238 "zlacon.f"
    jump = 4;
#line 239 "zlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 4) */
/*     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */

#line 244 "zlacon.f"
L90:
#line 245 "zlacon.f"
    jlast = j;
#line 246 "zlacon.f"
    j = izmax1_(n, &x[1], &c__1);
#line 247 "zlacon.f"
    if (z_abs(&x[jlast]) != z_abs(&x[j]) && iter < 5) {
#line 249 "zlacon.f"
	++iter;
#line 250 "zlacon.f"
	goto L50;
#line 251 "zlacon.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 255 "zlacon.f"
L100:
#line 256 "zlacon.f"
    altsgn = 1.;
#line 257 "zlacon.f"
    i__1 = *n;
#line 257 "zlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "zlacon.f"
	i__2 = i__;
#line 258 "zlacon.f"
	d__1 = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
#line 258 "zlacon.f"
	z__1.r = d__1, z__1.i = 0.;
#line 258 "zlacon.f"
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 259 "zlacon.f"
	altsgn = -altsgn;
#line 260 "zlacon.f"
/* L110: */
#line 260 "zlacon.f"
    }
#line 261 "zlacon.f"
    *kase = 1;
#line 262 "zlacon.f"
    jump = 5;
#line 263 "zlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 268 "zlacon.f"
L120:
#line 269 "zlacon.f"
    temp = dzsum1_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 270 "zlacon.f"
    if (temp > *est) {
#line 271 "zlacon.f"
	zcopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 272 "zlacon.f"
	*est = temp;
#line 273 "zlacon.f"
    }

#line 275 "zlacon.f"
L130:
#line 276 "zlacon.f"
    *kase = 0;
#line 277 "zlacon.f"
    return 0;

/*     End of ZLACON */

} /* zlacon_ */

