#line 1 "slacon.f"
/* slacon.f -- translated by f2c (version 20100827).
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

#line 1 "slacon.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 1.;

/* > \brief \b SLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLACON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLACON( N, V, X, ISGN, EST, KASE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       REAL               EST */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISGN( * ) */
/*       REAL               V( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLACON estimates the 1-norm of a square, real matrix A. */
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
/* >          V is REAL array, dimension (N) */
/* >         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
/* >         (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (N) */
/* >         On an intermediate return, X should be overwritten by */
/* >               A * X,   if KASE=1, */
/* >               A**T * X,  if KASE=2, */
/* >         and SLACON must be re-called with all the other parameters */
/* >         unchanged. */
/* > \endverbatim */
/* > */
/* > \param[out] ISGN */
/* > \verbatim */
/* >          ISGN is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* >          EST is REAL */
/* >         On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* >         unchanged from the previous call to SLACON. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to SLACON, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**T * X. */
/* >         On the final return from SLACON, KASE will again be 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >  Nick Higham, University of Manchester. \n */
/* >  Originally named SONEST, dated March 16, 1988. */

/* > \par References: */
/*  ================ */
/* > */
/* >  N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/* >  a real or complex matrix, with applications to condition estimation", */
/* >  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slacon_(integer *n, doublereal *v, doublereal *x, 
	integer *isgn, doublereal *est, integer *kase)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, iter;
    static doublereal temp;
    static integer jump, jlast;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal altsgn, estold;


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

#line 160 "slacon.f"
    /* Parameter adjustments */
#line 160 "slacon.f"
    --isgn;
#line 160 "slacon.f"
    --x;
#line 160 "slacon.f"
    --v;
#line 160 "slacon.f"

#line 160 "slacon.f"
    /* Function Body */
#line 160 "slacon.f"
    if (*kase == 0) {
#line 161 "slacon.f"
	i__1 = *n;
#line 161 "slacon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 162 "slacon.f"
	    x[i__] = 1. / (doublereal) (*n);
#line 163 "slacon.f"
/* L10: */
#line 163 "slacon.f"
	}
#line 164 "slacon.f"
	*kase = 1;
#line 165 "slacon.f"
	jump = 1;
#line 166 "slacon.f"
	return 0;
#line 167 "slacon.f"
    }

#line 169 "slacon.f"
    switch (jump) {
#line 169 "slacon.f"
	case 1:  goto L20;
#line 169 "slacon.f"
	case 2:  goto L40;
#line 169 "slacon.f"
	case 3:  goto L70;
#line 169 "slacon.f"
	case 4:  goto L110;
#line 169 "slacon.f"
	case 5:  goto L140;
#line 169 "slacon.f"
    }

/*     ................ ENTRY   (JUMP = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 174 "slacon.f"
L20:
#line 175 "slacon.f"
    if (*n == 1) {
#line 176 "slacon.f"
	v[1] = x[1];
#line 177 "slacon.f"
	*est = abs(v[1]);
/*        ... QUIT */
#line 179 "slacon.f"
	goto L150;
#line 180 "slacon.f"
    }
#line 181 "slacon.f"
    *est = sasum_(n, &x[1], &c__1);

#line 183 "slacon.f"
    i__1 = *n;
#line 183 "slacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 184 "slacon.f"
	x[i__] = d_sign(&c_b11, &x[i__]);
#line 185 "slacon.f"
	isgn[i__] = i_dnnt(&x[i__]);
#line 186 "slacon.f"
/* L30: */
#line 186 "slacon.f"
    }
#line 187 "slacon.f"
    *kase = 2;
#line 188 "slacon.f"
    jump = 2;
#line 189 "slacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

#line 194 "slacon.f"
L40:
#line 195 "slacon.f"
    j = isamax_(n, &x[1], &c__1);
#line 196 "slacon.f"
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 200 "slacon.f"
L50:
#line 201 "slacon.f"
    i__1 = *n;
#line 201 "slacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 202 "slacon.f"
	x[i__] = 0.;
#line 203 "slacon.f"
/* L60: */
#line 203 "slacon.f"
    }
#line 204 "slacon.f"
    x[j] = 1.;
#line 205 "slacon.f"
    *kase = 1;
#line 206 "slacon.f"
    jump = 3;
#line 207 "slacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 212 "slacon.f"
L70:
#line 213 "slacon.f"
    scopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 214 "slacon.f"
    estold = *est;
#line 215 "slacon.f"
    *est = sasum_(n, &v[1], &c__1);
#line 216 "slacon.f"
    i__1 = *n;
#line 216 "slacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 217 "slacon.f"
	d__1 = d_sign(&c_b11, &x[i__]);
#line 217 "slacon.f"
	if (i_dnnt(&d__1) != isgn[i__]) {
#line 217 "slacon.f"
	    goto L90;
#line 217 "slacon.f"
	}
#line 219 "slacon.f"
/* L80: */
#line 219 "slacon.f"
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
#line 221 "slacon.f"
    goto L120;

#line 223 "slacon.f"
L90:
/*     TEST FOR CYCLING. */
#line 225 "slacon.f"
    if (*est <= estold) {
#line 225 "slacon.f"
	goto L120;
#line 225 "slacon.f"
    }

#line 228 "slacon.f"
    i__1 = *n;
#line 228 "slacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "slacon.f"
	x[i__] = d_sign(&c_b11, &x[i__]);
#line 230 "slacon.f"
	isgn[i__] = i_dnnt(&x[i__]);
#line 231 "slacon.f"
/* L100: */
#line 231 "slacon.f"
    }
#line 232 "slacon.f"
    *kase = 2;
#line 233 "slacon.f"
    jump = 4;
#line 234 "slacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 4) */
/*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

#line 239 "slacon.f"
L110:
#line 240 "slacon.f"
    jlast = j;
#line 241 "slacon.f"
    j = isamax_(n, &x[1], &c__1);
#line 242 "slacon.f"
    if (x[jlast] != (d__1 = x[j], abs(d__1)) && iter < 5) {
#line 243 "slacon.f"
	++iter;
#line 244 "slacon.f"
	goto L50;
#line 245 "slacon.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 249 "slacon.f"
L120:
#line 250 "slacon.f"
    altsgn = 1.;
#line 251 "slacon.f"
    i__1 = *n;
#line 251 "slacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 252 "slacon.f"
	x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 
		1.);
#line 253 "slacon.f"
	altsgn = -altsgn;
#line 254 "slacon.f"
/* L130: */
#line 254 "slacon.f"
    }
#line 255 "slacon.f"
    *kase = 1;
#line 256 "slacon.f"
    jump = 5;
#line 257 "slacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 262 "slacon.f"
L140:
#line 263 "slacon.f"
    temp = sasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 264 "slacon.f"
    if (temp > *est) {
#line 265 "slacon.f"
	scopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 266 "slacon.f"
	*est = temp;
#line 267 "slacon.f"
    }

#line 269 "slacon.f"
L150:
#line 270 "slacon.f"
    *kase = 0;
#line 271 "slacon.f"
    return 0;

/*     End of SLACON */

} /* slacon_ */

