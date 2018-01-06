#line 1 "dlaed9.f"
/* dlaed9.f -- translated by f2c (version 20100827).
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

#line 1 "dlaed9.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLAED9 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Us
ed when the original matrix is dense. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAED9 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed9.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed9.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed9.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W, */
/*                          S, LDS, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, KSTART, KSTOP, LDQ, LDS, N */
/*       DOUBLE PRECISION   RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ), */
/*      $                   W( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED9 finds the roots of the secular equation, as defined by the */
/* > values in D, Z, and RHO, between KSTART and KSTOP.  It makes the */
/* > appropriate calls to DLAED4 and then stores the new matrix of */
/* > eigenvectors for use in calculating the next level of Z vectors. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of terms in the rational function to be solved by */
/* >          DLAED4.  K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KSTART */
/* > \verbatim */
/* >          KSTART is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] KSTOP */
/* > \verbatim */
/* >          KSTOP is INTEGER */
/* >          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP */
/* >          are to be computed.  1 <= KSTART <= KSTOP <= K. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of rows and columns in the Q matrix. */
/* >          N >= K (delation may result in N > K). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          D(I) contains the updated eigenvalues */
/* >          for KSTART <= I <= KSTOP. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >          The value of the parameter in the rank one update equation. */
/* >          RHO >= 0 required. */
/* > \endverbatim */
/* > */
/* > \param[in] DLAMDA */
/* > \verbatim */
/* >          DLAMDA is DOUBLE PRECISION array, dimension (K) */
/* >          The first K elements of this array contain the old roots */
/* >          of the deflated updating problem.  These are the poles */
/* >          of the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (K) */
/* >          The first K elements of this array contain the components */
/* >          of the deflation-adjusted updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (LDS, K) */
/* >          Will contain the eigenvectors of the repaired matrix which */
/* >          will be stored for subsequent Z vector calculation and */
/* >          multiplied by the previously accumulated eigenvectors */
/* >          to update the system. */
/* > \endverbatim */
/* > */
/* > \param[in] LDS */
/* > \verbatim */
/* >          LDS is INTEGER */
/* >          The leading dimension of S.  LDS >= max( 1, K ). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = 1, an eigenvalue did not converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlaed9_(integer *k, integer *kstart, integer *kstop, 
	integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *
	rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, 
	integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, s_dim1, s_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j;
    static doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaed4_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 193 "dlaed9.f"
    /* Parameter adjustments */
#line 193 "dlaed9.f"
    --d__;
#line 193 "dlaed9.f"
    q_dim1 = *ldq;
#line 193 "dlaed9.f"
    q_offset = 1 + q_dim1;
#line 193 "dlaed9.f"
    q -= q_offset;
#line 193 "dlaed9.f"
    --dlamda;
#line 193 "dlaed9.f"
    --w;
#line 193 "dlaed9.f"
    s_dim1 = *lds;
#line 193 "dlaed9.f"
    s_offset = 1 + s_dim1;
#line 193 "dlaed9.f"
    s -= s_offset;
#line 193 "dlaed9.f"

#line 193 "dlaed9.f"
    /* Function Body */
#line 193 "dlaed9.f"
    *info = 0;

#line 195 "dlaed9.f"
    if (*k < 0) {
#line 196 "dlaed9.f"
	*info = -1;
#line 197 "dlaed9.f"
    } else if (*kstart < 1 || *kstart > max(1,*k)) {
#line 198 "dlaed9.f"
	*info = -2;
#line 199 "dlaed9.f"
    } else if (max(1,*kstop) < *kstart || *kstop > max(1,*k)) {
#line 201 "dlaed9.f"
	*info = -3;
#line 202 "dlaed9.f"
    } else if (*n < *k) {
#line 203 "dlaed9.f"
	*info = -4;
#line 204 "dlaed9.f"
    } else if (*ldq < max(1,*k)) {
#line 205 "dlaed9.f"
	*info = -7;
#line 206 "dlaed9.f"
    } else if (*lds < max(1,*k)) {
#line 207 "dlaed9.f"
	*info = -12;
#line 208 "dlaed9.f"
    }
#line 209 "dlaed9.f"
    if (*info != 0) {
#line 210 "dlaed9.f"
	i__1 = -(*info);
#line 210 "dlaed9.f"
	xerbla_("DLAED9", &i__1, (ftnlen)6);
#line 211 "dlaed9.f"
	return 0;
#line 212 "dlaed9.f"
    }

/*     Quick return if possible */

#line 216 "dlaed9.f"
    if (*k == 0) {
#line 216 "dlaed9.f"
	return 0;
#line 216 "dlaed9.f"
    }

/*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can */
/*     be computed with high relative accuracy (barring over/underflow). */
/*     This is a problem on machines without a guard digit in */
/*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2). */
/*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I), */
/*     which on any of these machines zeros out the bottommost */
/*     bit of DLAMDA(I) if it is 1; this makes the subsequent */
/*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation */
/*     occurs. On binary machines with a guard digit (almost all */
/*     machines) it does not change DLAMDA(I) at all. On hexadecimal */
/*     and decimal machines with a guard digit, it slightly */
/*     changes the bottommost bits of DLAMDA(I). It does not account */
/*     for hexadecimal or decimal machines without guard digits */
/*     (we know of none). We use a subroutine call to compute */
/*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating */
/*     this code. */

#line 236 "dlaed9.f"
    i__1 = *n;
#line 236 "dlaed9.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 237 "dlaed9.f"
	dlamda[i__] = dlamc3_(&dlamda[i__], &dlamda[i__]) - dlamda[i__];
#line 238 "dlaed9.f"
/* L10: */
#line 238 "dlaed9.f"
    }

#line 240 "dlaed9.f"
    i__1 = *kstop;
#line 240 "dlaed9.f"
    for (j = *kstart; j <= i__1; ++j) {
#line 241 "dlaed9.f"
	dlaed4_(k, &j, &dlamda[1], &w[1], &q[j * q_dim1 + 1], rho, &d__[j], 
		info);

/*        If the zero finder fails, the computation is terminated. */

#line 245 "dlaed9.f"
	if (*info != 0) {
#line 245 "dlaed9.f"
	    goto L120;
#line 245 "dlaed9.f"
	}
#line 247 "dlaed9.f"
/* L20: */
#line 247 "dlaed9.f"
    }

#line 249 "dlaed9.f"
    if (*k == 1 || *k == 2) {
#line 250 "dlaed9.f"
	i__1 = *k;
#line 250 "dlaed9.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "dlaed9.f"
	    i__2 = *k;
#line 251 "dlaed9.f"
	    for (j = 1; j <= i__2; ++j) {
#line 252 "dlaed9.f"
		s[j + i__ * s_dim1] = q[j + i__ * q_dim1];
#line 253 "dlaed9.f"
/* L30: */
#line 253 "dlaed9.f"
	    }
#line 254 "dlaed9.f"
/* L40: */
#line 254 "dlaed9.f"
	}
#line 255 "dlaed9.f"
	goto L120;
#line 256 "dlaed9.f"
    }

/*     Compute updated W. */

#line 260 "dlaed9.f"
    dcopy_(k, &w[1], &c__1, &s[s_offset], &c__1);

/*     Initialize W(I) = Q(I,I) */

#line 264 "dlaed9.f"
    i__1 = *ldq + 1;
#line 264 "dlaed9.f"
    dcopy_(k, &q[q_offset], &i__1, &w[1], &c__1);
#line 265 "dlaed9.f"
    i__1 = *k;
#line 265 "dlaed9.f"
    for (j = 1; j <= i__1; ++j) {
#line 266 "dlaed9.f"
	i__2 = j - 1;
#line 266 "dlaed9.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 267 "dlaed9.f"
	    w[i__] *= q[i__ + j * q_dim1] / (dlamda[i__] - dlamda[j]);
#line 268 "dlaed9.f"
/* L50: */
#line 268 "dlaed9.f"
	}
#line 269 "dlaed9.f"
	i__2 = *k;
#line 269 "dlaed9.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 270 "dlaed9.f"
	    w[i__] *= q[i__ + j * q_dim1] / (dlamda[i__] - dlamda[j]);
#line 271 "dlaed9.f"
/* L60: */
#line 271 "dlaed9.f"
	}
#line 272 "dlaed9.f"
/* L70: */
#line 272 "dlaed9.f"
    }
#line 273 "dlaed9.f"
    i__1 = *k;
#line 273 "dlaed9.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "dlaed9.f"
	d__1 = sqrt(-w[i__]);
#line 274 "dlaed9.f"
	w[i__] = d_sign(&d__1, &s[i__ + s_dim1]);
#line 275 "dlaed9.f"
/* L80: */
#line 275 "dlaed9.f"
    }

/*     Compute eigenvectors of the modified rank-1 modification. */

#line 279 "dlaed9.f"
    i__1 = *k;
#line 279 "dlaed9.f"
    for (j = 1; j <= i__1; ++j) {
#line 280 "dlaed9.f"
	i__2 = *k;
#line 280 "dlaed9.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 281 "dlaed9.f"
	    q[i__ + j * q_dim1] = w[i__] / q[i__ + j * q_dim1];
#line 282 "dlaed9.f"
/* L90: */
#line 282 "dlaed9.f"
	}
#line 283 "dlaed9.f"
	temp = dnrm2_(k, &q[j * q_dim1 + 1], &c__1);
#line 284 "dlaed9.f"
	i__2 = *k;
#line 284 "dlaed9.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 285 "dlaed9.f"
	    s[i__ + j * s_dim1] = q[i__ + j * q_dim1] / temp;
#line 286 "dlaed9.f"
/* L100: */
#line 286 "dlaed9.f"
	}
#line 287 "dlaed9.f"
/* L110: */
#line 287 "dlaed9.f"
    }

#line 289 "dlaed9.f"
L120:
#line 290 "dlaed9.f"
    return 0;

/*     End of DLAED9 */

} /* dlaed9_ */

