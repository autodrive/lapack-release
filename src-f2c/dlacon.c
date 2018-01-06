#line 1 "dlacon.f"
/* dlacon.f -- translated by f2c (version 20100827).
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

#line 1 "dlacon.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 1.;

/* > \brief \b DLACON estimates the 1-norm of a square matrix, using reverse communication for evaluating matr
ix-vector products. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLACON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            KASE, N */
/*       DOUBLE PRECISION   EST */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            ISGN( * ) */
/*       DOUBLE PRECISION   V( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLACON estimates the 1-norm of a square, real matrix A. */
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
/* >          V is DOUBLE PRECISION array, dimension (N) */
/* >         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
/* >         (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (N) */
/* >         On an intermediate return, X should be overwritten by */
/* >               A * X,   if KASE=1, */
/* >               A**T * X,  if KASE=2, */
/* >         and DLACON must be re-called with all the other parameters */
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
/* >          EST is DOUBLE PRECISION */
/* >         On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* >         unchanged from the previous call to DLACON. */
/* >         On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* >          KASE is INTEGER */
/* >         On the initial call to DLACON, KASE should be 0. */
/* >         On an intermediate return, KASE will be 1 or 2, indicating */
/* >         whether X should be overwritten by A * X  or A**T * X. */
/* >         On the final return from DLACON, KASE will again be 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

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
/* Subroutine */ int dlacon_(integer *n, doublereal *v, doublereal *x, 
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
    static integer jump;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer jlast;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
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

#line 160 "dlacon.f"
    /* Parameter adjustments */
#line 160 "dlacon.f"
    --isgn;
#line 160 "dlacon.f"
    --x;
#line 160 "dlacon.f"
    --v;
#line 160 "dlacon.f"

#line 160 "dlacon.f"
    /* Function Body */
#line 160 "dlacon.f"
    if (*kase == 0) {
#line 161 "dlacon.f"
	i__1 = *n;
#line 161 "dlacon.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 162 "dlacon.f"
	    x[i__] = 1. / (doublereal) (*n);
#line 163 "dlacon.f"
/* L10: */
#line 163 "dlacon.f"
	}
#line 164 "dlacon.f"
	*kase = 1;
#line 165 "dlacon.f"
	jump = 1;
#line 166 "dlacon.f"
	return 0;
#line 167 "dlacon.f"
    }

#line 169 "dlacon.f"
    switch (jump) {
#line 169 "dlacon.f"
	case 1:  goto L20;
#line 169 "dlacon.f"
	case 2:  goto L40;
#line 169 "dlacon.f"
	case 3:  goto L70;
#line 169 "dlacon.f"
	case 4:  goto L110;
#line 169 "dlacon.f"
	case 5:  goto L140;
#line 169 "dlacon.f"
    }

/*     ................ ENTRY   (JUMP = 1) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

#line 174 "dlacon.f"
L20:
#line 175 "dlacon.f"
    if (*n == 1) {
#line 176 "dlacon.f"
	v[1] = x[1];
#line 177 "dlacon.f"
	*est = abs(v[1]);
/*        ... QUIT */
#line 179 "dlacon.f"
	goto L150;
#line 180 "dlacon.f"
    }
#line 181 "dlacon.f"
    *est = dasum_(n, &x[1], &c__1);

#line 183 "dlacon.f"
    i__1 = *n;
#line 183 "dlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 184 "dlacon.f"
	x[i__] = d_sign(&c_b11, &x[i__]);
#line 185 "dlacon.f"
	isgn[i__] = i_dnnt(&x[i__]);
#line 186 "dlacon.f"
/* L30: */
#line 186 "dlacon.f"
    }
#line 187 "dlacon.f"
    *kase = 2;
#line 188 "dlacon.f"
    jump = 2;
#line 189 "dlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 2) */
/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

#line 194 "dlacon.f"
L40:
#line 195 "dlacon.f"
    j = idamax_(n, &x[1], &c__1);
#line 196 "dlacon.f"
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

#line 200 "dlacon.f"
L50:
#line 201 "dlacon.f"
    i__1 = *n;
#line 201 "dlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 202 "dlacon.f"
	x[i__] = 0.;
#line 203 "dlacon.f"
/* L60: */
#line 203 "dlacon.f"
    }
#line 204 "dlacon.f"
    x[j] = 1.;
#line 205 "dlacon.f"
    *kase = 1;
#line 206 "dlacon.f"
    jump = 3;
#line 207 "dlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 3) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 212 "dlacon.f"
L70:
#line 213 "dlacon.f"
    dcopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 214 "dlacon.f"
    estold = *est;
#line 215 "dlacon.f"
    *est = dasum_(n, &v[1], &c__1);
#line 216 "dlacon.f"
    i__1 = *n;
#line 216 "dlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 217 "dlacon.f"
	d__1 = d_sign(&c_b11, &x[i__]);
#line 217 "dlacon.f"
	if (i_dnnt(&d__1) != isgn[i__]) {
#line 217 "dlacon.f"
	    goto L90;
#line 217 "dlacon.f"
	}
#line 219 "dlacon.f"
/* L80: */
#line 219 "dlacon.f"
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
#line 221 "dlacon.f"
    goto L120;

#line 223 "dlacon.f"
L90:
/*     TEST FOR CYCLING. */
#line 225 "dlacon.f"
    if (*est <= estold) {
#line 225 "dlacon.f"
	goto L120;
#line 225 "dlacon.f"
    }

#line 228 "dlacon.f"
    i__1 = *n;
#line 228 "dlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "dlacon.f"
	x[i__] = d_sign(&c_b11, &x[i__]);
#line 230 "dlacon.f"
	isgn[i__] = i_dnnt(&x[i__]);
#line 231 "dlacon.f"
/* L100: */
#line 231 "dlacon.f"
    }
#line 232 "dlacon.f"
    *kase = 2;
#line 233 "dlacon.f"
    jump = 4;
#line 234 "dlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 4) */
/*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

#line 239 "dlacon.f"
L110:
#line 240 "dlacon.f"
    jlast = j;
#line 241 "dlacon.f"
    j = idamax_(n, &x[1], &c__1);
#line 242 "dlacon.f"
    if (x[jlast] != (d__1 = x[j], abs(d__1)) && iter < 5) {
#line 243 "dlacon.f"
	++iter;
#line 244 "dlacon.f"
	goto L50;
#line 245 "dlacon.f"
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

#line 249 "dlacon.f"
L120:
#line 250 "dlacon.f"
    altsgn = 1.;
#line 251 "dlacon.f"
    i__1 = *n;
#line 251 "dlacon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 252 "dlacon.f"
	x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 
		1.);
#line 253 "dlacon.f"
	altsgn = -altsgn;
#line 254 "dlacon.f"
/* L130: */
#line 254 "dlacon.f"
    }
#line 255 "dlacon.f"
    *kase = 1;
#line 256 "dlacon.f"
    jump = 5;
#line 257 "dlacon.f"
    return 0;

/*     ................ ENTRY   (JUMP = 5) */
/*     X HAS BEEN OVERWRITTEN BY A*X. */

#line 262 "dlacon.f"
L140:
#line 263 "dlacon.f"
    temp = dasum_(n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
#line 264 "dlacon.f"
    if (temp > *est) {
#line 265 "dlacon.f"
	dcopy_(n, &x[1], &c__1, &v[1], &c__1);
#line 266 "dlacon.f"
	*est = temp;
#line 267 "dlacon.f"
    }

#line 269 "dlacon.f"
L150:
#line 270 "dlacon.f"
    *kase = 0;
#line 271 "dlacon.f"
    return 0;

/*     End of DLACON */

} /* dlacon_ */

