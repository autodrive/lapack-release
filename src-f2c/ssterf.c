#line 1 "ssterf.f"
/* ssterf.f -- translated by f2c (version 20100827).
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

#line 1 "ssterf.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b32 = 1.;

/* > \brief \b SSTERF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSTERF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssterf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssterf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssterf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSTERF( N, D, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTERF computes all eigenvalues of a symmetric tridiagonal matrix */
/* > using the Pal-Walker-Kahan variant of the QL or QR algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  the algorithm failed to find all of the eigenvalues in */
/* >                a total of 30*N iterations; if INFO = i, then i */
/* >                elements of E have not converged to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssterf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal c__;
    static integer i__, l, m;
    static doublereal p, r__, s;
    static integer l1;
    static doublereal bb, rt1, rt2, eps, rte;
    static integer lsv;
    static doublereal eps2, oldc;
    static integer lend, jtot;
    extern /* Subroutine */ int slae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal gamma, alpha, sigma, anorm;
    extern doublereal slapy2_(doublereal *, doublereal *);
    static integer iscale;
    static doublereal oldgam;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer lendsv;
    static doublereal ssfmin;
    static integer nmaxit;
    static doublereal ssfmax;
    extern doublereal slanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int slasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

/*     Test the input parameters. */

#line 131 "ssterf.f"
    /* Parameter adjustments */
#line 131 "ssterf.f"
    --e;
#line 131 "ssterf.f"
    --d__;
#line 131 "ssterf.f"

#line 131 "ssterf.f"
    /* Function Body */
#line 131 "ssterf.f"
    *info = 0;

/*     Quick return if possible */

#line 135 "ssterf.f"
    if (*n < 0) {
#line 136 "ssterf.f"
	*info = -1;
#line 137 "ssterf.f"
	i__1 = -(*info);
#line 137 "ssterf.f"
	xerbla_("SSTERF", &i__1, (ftnlen)6);
#line 138 "ssterf.f"
	return 0;
#line 139 "ssterf.f"
    }
#line 140 "ssterf.f"
    if (*n <= 1) {
#line 140 "ssterf.f"
	return 0;
#line 140 "ssterf.f"
    }

/*     Determine the unit roundoff for this environment. */

#line 145 "ssterf.f"
    eps = slamch_("E", (ftnlen)1);
/* Computing 2nd power */
#line 146 "ssterf.f"
    d__1 = eps;
#line 146 "ssterf.f"
    eps2 = d__1 * d__1;
#line 147 "ssterf.f"
    safmin = slamch_("S", (ftnlen)1);
#line 148 "ssterf.f"
    safmax = 1. / safmin;
#line 149 "ssterf.f"
    ssfmax = sqrt(safmax) / 3.;
#line 150 "ssterf.f"
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues of the tridiagonal matrix. */

#line 154 "ssterf.f"
    nmaxit = *n * 30;
#line 155 "ssterf.f"
    sigma = 0.;
#line 156 "ssterf.f"
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

#line 162 "ssterf.f"
    l1 = 1;

#line 164 "ssterf.f"
L10:
#line 165 "ssterf.f"
    if (l1 > *n) {
#line 165 "ssterf.f"
	goto L170;
#line 165 "ssterf.f"
    }
#line 167 "ssterf.f"
    if (l1 > 1) {
#line 167 "ssterf.f"
	e[l1 - 1] = 0.;
#line 167 "ssterf.f"
    }
#line 169 "ssterf.f"
    i__1 = *n - 1;
#line 169 "ssterf.f"
    for (m = l1; m <= i__1; ++m) {
#line 170 "ssterf.f"
	if ((d__3 = e[m], abs(d__3)) <= sqrt((d__1 = d__[m], abs(d__1))) * 
		sqrt((d__2 = d__[m + 1], abs(d__2))) * eps) {
#line 172 "ssterf.f"
	    e[m] = 0.;
#line 173 "ssterf.f"
	    goto L30;
#line 174 "ssterf.f"
	}
#line 175 "ssterf.f"
/* L20: */
#line 175 "ssterf.f"
    }
#line 176 "ssterf.f"
    m = *n;

#line 178 "ssterf.f"
L30:
#line 179 "ssterf.f"
    l = l1;
#line 180 "ssterf.f"
    lsv = l;
#line 181 "ssterf.f"
    lend = m;
#line 182 "ssterf.f"
    lendsv = lend;
#line 183 "ssterf.f"
    l1 = m + 1;
#line 184 "ssterf.f"
    if (lend == l) {
#line 184 "ssterf.f"
	goto L10;
#line 184 "ssterf.f"
    }

/*     Scale submatrix in rows and columns L to LEND */

#line 189 "ssterf.f"
    i__1 = lend - l + 1;
#line 189 "ssterf.f"
    anorm = slanst_("M", &i__1, &d__[l], &e[l], (ftnlen)1);
#line 190 "ssterf.f"
    iscale = 0;
#line 191 "ssterf.f"
    if (anorm == 0.) {
#line 191 "ssterf.f"
	goto L10;
#line 191 "ssterf.f"
    }
#line 193 "ssterf.f"
    if (anorm > ssfmax) {
#line 194 "ssterf.f"
	iscale = 1;
#line 195 "ssterf.f"
	i__1 = lend - l + 1;
#line 195 "ssterf.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 197 "ssterf.f"
	i__1 = lend - l;
#line 197 "ssterf.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 199 "ssterf.f"
    } else if (anorm < ssfmin) {
#line 200 "ssterf.f"
	iscale = 2;
#line 201 "ssterf.f"
	i__1 = lend - l + 1;
#line 201 "ssterf.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 203 "ssterf.f"
	i__1 = lend - l;
#line 203 "ssterf.f"
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 205 "ssterf.f"
    }

#line 207 "ssterf.f"
    i__1 = lend - 1;
#line 207 "ssterf.f"
    for (i__ = l; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 208 "ssterf.f"
	d__1 = e[i__];
#line 208 "ssterf.f"
	e[i__] = d__1 * d__1;
#line 209 "ssterf.f"
/* L40: */
#line 209 "ssterf.f"
    }

/*     Choose between QL and QR iteration */

#line 213 "ssterf.f"
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
#line 214 "ssterf.f"
	lend = lsv;
#line 215 "ssterf.f"
	l = lendsv;
#line 216 "ssterf.f"
    }

#line 218 "ssterf.f"
    if (lend >= l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

#line 224 "ssterf.f"
L50:
#line 225 "ssterf.f"
	if (l != lend) {
#line 226 "ssterf.f"
	    i__1 = lend - 1;
#line 226 "ssterf.f"
	    for (m = l; m <= i__1; ++m) {
#line 227 "ssterf.f"
		if ((d__2 = e[m], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
			+ 1], abs(d__1))) {
#line 227 "ssterf.f"
		    goto L70;
#line 227 "ssterf.f"
		}
#line 229 "ssterf.f"
/* L60: */
#line 229 "ssterf.f"
	    }
#line 230 "ssterf.f"
	}
#line 231 "ssterf.f"
	m = lend;

#line 233 "ssterf.f"
L70:
#line 234 "ssterf.f"
	if (m < lend) {
#line 234 "ssterf.f"
	    e[m] = 0.;
#line 234 "ssterf.f"
	}
#line 236 "ssterf.f"
	p = d__[l];
#line 237 "ssterf.f"
	if (m == l) {
#line 237 "ssterf.f"
	    goto L90;
#line 237 "ssterf.f"
	}

/*        If remaining matrix is 2 by 2, use SLAE2 to compute its */
/*        eigenvalues. */

#line 243 "ssterf.f"
	if (m == l + 1) {
#line 244 "ssterf.f"
	    rte = sqrt(e[l]);
#line 245 "ssterf.f"
	    slae2_(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
#line 246 "ssterf.f"
	    d__[l] = rt1;
#line 247 "ssterf.f"
	    d__[l + 1] = rt2;
#line 248 "ssterf.f"
	    e[l] = 0.;
#line 249 "ssterf.f"
	    l += 2;
#line 250 "ssterf.f"
	    if (l <= lend) {
#line 250 "ssterf.f"
		goto L50;
#line 250 "ssterf.f"
	    }
#line 252 "ssterf.f"
	    goto L150;
#line 253 "ssterf.f"
	}

#line 255 "ssterf.f"
	if (jtot == nmaxit) {
#line 255 "ssterf.f"
	    goto L150;
#line 255 "ssterf.f"
	}
#line 257 "ssterf.f"
	++jtot;

/*        Form shift. */

#line 261 "ssterf.f"
	rte = sqrt(e[l]);
#line 262 "ssterf.f"
	sigma = (d__[l + 1] - p) / (rte * 2.);
#line 263 "ssterf.f"
	r__ = slapy2_(&sigma, &c_b32);
#line 264 "ssterf.f"
	sigma = p - rte / (sigma + d_sign(&r__, &sigma));

#line 266 "ssterf.f"
	c__ = 1.;
#line 267 "ssterf.f"
	s = 0.;
#line 268 "ssterf.f"
	gamma = d__[m] - sigma;
#line 269 "ssterf.f"
	p = gamma * gamma;

/*        Inner loop */

#line 273 "ssterf.f"
	i__1 = l;
#line 273 "ssterf.f"
	for (i__ = m - 1; i__ >= i__1; --i__) {
#line 274 "ssterf.f"
	    bb = e[i__];
#line 275 "ssterf.f"
	    r__ = p + bb;
#line 276 "ssterf.f"
	    if (i__ != m - 1) {
#line 276 "ssterf.f"
		e[i__ + 1] = s * r__;
#line 276 "ssterf.f"
	    }
#line 278 "ssterf.f"
	    oldc = c__;
#line 279 "ssterf.f"
	    c__ = p / r__;
#line 280 "ssterf.f"
	    s = bb / r__;
#line 281 "ssterf.f"
	    oldgam = gamma;
#line 282 "ssterf.f"
	    alpha = d__[i__];
#line 283 "ssterf.f"
	    gamma = c__ * (alpha - sigma) - s * oldgam;
#line 284 "ssterf.f"
	    d__[i__ + 1] = oldgam + (alpha - gamma);
#line 285 "ssterf.f"
	    if (c__ != 0.) {
#line 286 "ssterf.f"
		p = gamma * gamma / c__;
#line 287 "ssterf.f"
	    } else {
#line 288 "ssterf.f"
		p = oldc * bb;
#line 289 "ssterf.f"
	    }
#line 290 "ssterf.f"
/* L80: */
#line 290 "ssterf.f"
	}

#line 292 "ssterf.f"
	e[l] = s * p;
#line 293 "ssterf.f"
	d__[l] = sigma + gamma;
#line 294 "ssterf.f"
	goto L50;

/*        Eigenvalue found. */

#line 298 "ssterf.f"
L90:
#line 299 "ssterf.f"
	d__[l] = p;

#line 301 "ssterf.f"
	++l;
#line 302 "ssterf.f"
	if (l <= lend) {
#line 302 "ssterf.f"
	    goto L50;
#line 302 "ssterf.f"
	}
#line 304 "ssterf.f"
	goto L150;

#line 306 "ssterf.f"
    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

#line 312 "ssterf.f"
L100:
#line 313 "ssterf.f"
	i__1 = lend + 1;
#line 313 "ssterf.f"
	for (m = l; m >= i__1; --m) {
#line 314 "ssterf.f"
	    if ((d__2 = e[m - 1], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
		    - 1], abs(d__1))) {
#line 314 "ssterf.f"
		goto L120;
#line 314 "ssterf.f"
	    }
#line 316 "ssterf.f"
/* L110: */
#line 316 "ssterf.f"
	}
#line 317 "ssterf.f"
	m = lend;

#line 319 "ssterf.f"
L120:
#line 320 "ssterf.f"
	if (m > lend) {
#line 320 "ssterf.f"
	    e[m - 1] = 0.;
#line 320 "ssterf.f"
	}
#line 322 "ssterf.f"
	p = d__[l];
#line 323 "ssterf.f"
	if (m == l) {
#line 323 "ssterf.f"
	    goto L140;
#line 323 "ssterf.f"
	}

/*        If remaining matrix is 2 by 2, use SLAE2 to compute its */
/*        eigenvalues. */

#line 329 "ssterf.f"
	if (m == l - 1) {
#line 330 "ssterf.f"
	    rte = sqrt(e[l - 1]);
#line 331 "ssterf.f"
	    slae2_(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
#line 332 "ssterf.f"
	    d__[l] = rt1;
#line 333 "ssterf.f"
	    d__[l - 1] = rt2;
#line 334 "ssterf.f"
	    e[l - 1] = 0.;
#line 335 "ssterf.f"
	    l += -2;
#line 336 "ssterf.f"
	    if (l >= lend) {
#line 336 "ssterf.f"
		goto L100;
#line 336 "ssterf.f"
	    }
#line 338 "ssterf.f"
	    goto L150;
#line 339 "ssterf.f"
	}

#line 341 "ssterf.f"
	if (jtot == nmaxit) {
#line 341 "ssterf.f"
	    goto L150;
#line 341 "ssterf.f"
	}
#line 343 "ssterf.f"
	++jtot;

/*        Form shift. */

#line 347 "ssterf.f"
	rte = sqrt(e[l - 1]);
#line 348 "ssterf.f"
	sigma = (d__[l - 1] - p) / (rte * 2.);
#line 349 "ssterf.f"
	r__ = slapy2_(&sigma, &c_b32);
#line 350 "ssterf.f"
	sigma = p - rte / (sigma + d_sign(&r__, &sigma));

#line 352 "ssterf.f"
	c__ = 1.;
#line 353 "ssterf.f"
	s = 0.;
#line 354 "ssterf.f"
	gamma = d__[m] - sigma;
#line 355 "ssterf.f"
	p = gamma * gamma;

/*        Inner loop */

#line 359 "ssterf.f"
	i__1 = l - 1;
#line 359 "ssterf.f"
	for (i__ = m; i__ <= i__1; ++i__) {
#line 360 "ssterf.f"
	    bb = e[i__];
#line 361 "ssterf.f"
	    r__ = p + bb;
#line 362 "ssterf.f"
	    if (i__ != m) {
#line 362 "ssterf.f"
		e[i__ - 1] = s * r__;
#line 362 "ssterf.f"
	    }
#line 364 "ssterf.f"
	    oldc = c__;
#line 365 "ssterf.f"
	    c__ = p / r__;
#line 366 "ssterf.f"
	    s = bb / r__;
#line 367 "ssterf.f"
	    oldgam = gamma;
#line 368 "ssterf.f"
	    alpha = d__[i__ + 1];
#line 369 "ssterf.f"
	    gamma = c__ * (alpha - sigma) - s * oldgam;
#line 370 "ssterf.f"
	    d__[i__] = oldgam + (alpha - gamma);
#line 371 "ssterf.f"
	    if (c__ != 0.) {
#line 372 "ssterf.f"
		p = gamma * gamma / c__;
#line 373 "ssterf.f"
	    } else {
#line 374 "ssterf.f"
		p = oldc * bb;
#line 375 "ssterf.f"
	    }
#line 376 "ssterf.f"
/* L130: */
#line 376 "ssterf.f"
	}

#line 378 "ssterf.f"
	e[l - 1] = s * p;
#line 379 "ssterf.f"
	d__[l] = sigma + gamma;
#line 380 "ssterf.f"
	goto L100;

/*        Eigenvalue found. */

#line 384 "ssterf.f"
L140:
#line 385 "ssterf.f"
	d__[l] = p;

#line 387 "ssterf.f"
	--l;
#line 388 "ssterf.f"
	if (l >= lend) {
#line 388 "ssterf.f"
	    goto L100;
#line 388 "ssterf.f"
	}
#line 390 "ssterf.f"
	goto L150;

#line 392 "ssterf.f"
    }

/*     Undo scaling if necessary */

#line 396 "ssterf.f"
L150:
#line 397 "ssterf.f"
    if (iscale == 1) {
#line 397 "ssterf.f"
	i__1 = lendsv - lsv + 1;
#line 397 "ssterf.f"
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 397 "ssterf.f"
    }
#line 400 "ssterf.f"
    if (iscale == 2) {
#line 400 "ssterf.f"
	i__1 = lendsv - lsv + 1;
#line 400 "ssterf.f"
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 400 "ssterf.f"
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

#line 407 "ssterf.f"
    if (jtot < nmaxit) {
#line 407 "ssterf.f"
	goto L10;
#line 407 "ssterf.f"
    }
#line 409 "ssterf.f"
    i__1 = *n - 1;
#line 409 "ssterf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 410 "ssterf.f"
	if (e[i__] != 0.) {
#line 410 "ssterf.f"
	    ++(*info);
#line 410 "ssterf.f"
	}
#line 412 "ssterf.f"
/* L160: */
#line 412 "ssterf.f"
    }
#line 413 "ssterf.f"
    goto L180;

/*     Sort eigenvalues in increasing order. */

#line 417 "ssterf.f"
L170:
#line 418 "ssterf.f"
    slasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 420 "ssterf.f"
L180:
#line 421 "ssterf.f"
    return 0;

/*     End of SSTERF */

} /* ssterf_ */

