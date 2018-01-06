#line 1 "slacn2.f"
/* slacn2.f -- translated by f2c (version 20100827).
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

#line 1 "slacn2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 1.;

/* > \brief \b SLACN2 estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLACN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacn2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacn2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacn2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLACN2( N, V, X, ISGN, EST, KASE, ISAVE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       REAL               EST */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISGN( * ), ISAVE( 3 ) */
/*       REAL               V( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLACN2 estimates the 1-norm of a square, real matrix A. */
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
/* >         and SLACN2 must be re-called with all the other parameters */
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
/* >         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be */
/* >         unchanged from the previous call to SLACN2. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to SLACN2, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**T * X. */
/* >         On the final return from SLACN2, KASE will again be 0. */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Originally named SONEST, dated March 16, 1988. */
/* > */
/* >  This is a thread safe version of SLACON, which uses the array ISAVE */
/* >  in place of a SAVE statement, as follows: */
/* > */
/* >     SLACON     SLACN2 */
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
/* Subroutine */ int slacn2_(integer *n, doublereal *v, doublereal *x, 
	integer *isgn, doublereal *est, integer *kase, integer *isave)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal temp;
    static integer jlast;
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
/*     .. Executable Statements .. */

#line 178 "slacn2.f"
    /* Parameter adjustments */
#line 178 "slacn2.f"
    --isave;
#line 178 "slacn2.f"
    --isgn;
#line 178 "slacn2.f"
    --x;
#line 178 "slacn2.f"
    --v;
#line 178 "slacn2.f"

#line 178 "slacn2.f"
    /* Function Body */
#line 178 "slacn2.f"
    if (*kase == 0) {
#line 179 "slacn2.f"
	i__1 = *n;
#line 179 "slacn2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 180 "slacn2.f"
	    x[i__] = 1. / (doublereal) (*n);
#line 181 "slacn2.f"
/* L10: */
#line 181 "slacn2.f"
	}
#line 182 "slacn2.f"
	*kase = 1;
#line 183 "slacn2.f"
	isave[1] = 1;
#line 184 "slacn2.f"
	return 0;
#line 185 "slacn2.f"
    }

#line 187 "slacn2.f"
    switch (isave[1]) {
#line 187 "slacn2.f"
	case 1:  goto L20;
#line 187 "slacn2.f"
	case 2:  goto L40;
#line 187 "slacn2.f"
	case 3:  goto L70;
#line 187 "slacn2.f"
	case 4:  goto L110;
#line 187 "slacn2.f"
	case 5:  goto L140;
#line 187 "slacn2.f"
    }

/*     ................ ENTRY   (ISAVE( 1 ) = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 192 "slacn2.f"
L20:
#line 193 "slacn2.f"
    if (*n == 1) {
#line 194 "slacn2.f"
	v[1] = x[1];
#line 195 "slacn2.f"
	*est = abs(v[1]);
/*        ... QUIT */
#line 197 "slacn2.f"
	goto L150;
#line 198 "slacn2.f"
    }
#line 199 "slacn2.f"
    *est = sasum_(n, &x[1], &c__1);

#line 201 "slacn2.f"
    i__1 = *n;
#line 201 "slacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 202 "slacn2.f"
	x[i__] = d_sign(&c_b11, &x[i__]);
#line 203 "slacn2.f"
	isgn[i__] = i_dnnt(&x[i__]);
#line 204 "slacn2.f"
/* L30: */
#line 204 "slacn2.f"
    }
#line 205 "slacn2.f"
    *kase = 2;
#line 206 "slacn2.f"
    isave[1] = 2;
#line 207 "slacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

#line 212 "slacn2.f"
L40:
#line 213 "slacn2.f"
    isave[2] = isamax_(n, &x[1], &c__1);
#line 214 "slacn2.f"
    isave[3] = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 218 "slacn2.f"
L50:
#line 219 "slacn2.f"
    i__1 = *n;
#line 219 "slacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "slacn2.f"
	x[i__] = 0.;
#line 221 "slacn2.f"
/* L60: */
#line 221 "slacn2.f"
    }
#line 222 "slacn2.f"
    x[isave[2]] = 1.;
#line 223 "slacn2.f"
    *kase = 1;
#line 224 "slacn2.f"
    isave[1] = 3;
#line 225 "slacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 230 "slacn2.f"
L70:
#line 231 "slacn2.f"
    scopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 232 "slacn2.f"
    estold = *est;
#line 233 "slacn2.f"
    *est = sasum_(n, &v[1], &c__1);
#line 234 "slacn2.f"
    i__1 = *n;
#line 234 "slacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "slacn2.f"
	d__1 = d_sign(&c_b11, &x[i__]);
#line 235 "slacn2.f"
	if (i_dnnt(&d__1) != isgn[i__]) {
#line 235 "slacn2.f"
	    goto L90;
#line 235 "slacn2.f"
	}
#line 237 "slacn2.f"
/* L80: */
#line 237 "slacn2.f"
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
#line 239 "slacn2.f"
    goto L120;

#line 241 "slacn2.f"
L90:
/*     TEST FOR CYCLING. */
#line 243 "slacn2.f"
    if (*est <= estold) {
#line 243 "slacn2.f"
	goto L120;
#line 243 "slacn2.f"
    }

#line 246 "slacn2.f"
    i__1 = *n;
#line 246 "slacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 247 "slacn2.f"
	x[i__] = d_sign(&c_b11, &x[i__]);
#line 248 "slacn2.f"
	isgn[i__] = i_dnnt(&x[i__]);
#line 249 "slacn2.f"
/* L100: */
#line 249 "slacn2.f"
    }
#line 250 "slacn2.f"
    *kase = 2;
#line 251 "slacn2.f"
    isave[1] = 4;
#line 252 "slacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 4) */
/*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

#line 257 "slacn2.f"
L110:
#line 258 "slacn2.f"
    jlast = isave[2];
#line 259 "slacn2.f"
    isave[2] = isamax_(n, &x[1], &c__1);
#line 260 "slacn2.f"
    if (x[jlast] != (d__1 = x[isave[2]], abs(d__1)) && isave[3] < 5) {
#line 262 "slacn2.f"
	++isave[3];
#line 263 "slacn2.f"
	goto L50;
#line 264 "slacn2.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 268 "slacn2.f"
L120:
#line 269 "slacn2.f"
    altsgn = 1.;
#line 270 "slacn2.f"
    i__1 = *n;
#line 270 "slacn2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "slacn2.f"
	x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 
		1.);
#line 272 "slacn2.f"
	altsgn = -altsgn;
#line 273 "slacn2.f"
/* L130: */
#line 273 "slacn2.f"
    }
#line 274 "slacn2.f"
    *kase = 1;
#line 275 "slacn2.f"
    isave[1] = 5;
#line 276 "slacn2.f"
    return 0;

/*     ................ ENTRY   (ISAVE( 1 ) = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 281 "slacn2.f"
L140:
#line 282 "slacn2.f"
    temp = sasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 283 "slacn2.f"
    if (temp > *est) {
#line 284 "slacn2.f"
	scopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 285 "slacn2.f"
	*est = temp;
#line 286 "slacn2.f"
    }

#line 288 "slacn2.f"
L150:
#line 289 "slacn2.f"
    *kase = 0;
#line 290 "slacn2.f"
    return 0;

/*     End of SLACN2 */

} /* slacn2_ */

