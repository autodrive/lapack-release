#line 1 "dsterf.f"
/* dsterf.f -- translated by f2c (version 20100827).
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

#line 1 "dsterf.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b33 = 1.;

/* > \brief \b DSTERF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTERF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsterf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsterf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsterf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTERF( N, D, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTERF computes all eigenvalues of a symmetric tridiagonal matrix */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
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

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsterf_(integer *n, doublereal *d__, doublereal *e, 
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
    static integer lend;
    static doublereal rmax;
    static integer jtot;
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal gamma, alpha, sigma, anorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    static integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal oldgam, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal safmax;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static integer lendsv;
    static doublereal ssfmin;
    static integer nmaxit;
    static doublereal ssfmax;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

/*     Test the input parameters. */

#line 131 "dsterf.f"
    /* Parameter adjustments */
#line 131 "dsterf.f"
    --e;
#line 131 "dsterf.f"
    --d__;
#line 131 "dsterf.f"

#line 131 "dsterf.f"
    /* Function Body */
#line 131 "dsterf.f"
    *info = 0;

/*     Quick return if possible */

#line 135 "dsterf.f"
    if (*n < 0) {
#line 136 "dsterf.f"
	*info = -1;
#line 137 "dsterf.f"
	i__1 = -(*info);
#line 137 "dsterf.f"
	xerbla_("DSTERF", &i__1, (ftnlen)6);
#line 138 "dsterf.f"
	return 0;
#line 139 "dsterf.f"
    }
#line 140 "dsterf.f"
    if (*n <= 1) {
#line 140 "dsterf.f"
	return 0;
#line 140 "dsterf.f"
    }

/*     Determine the unit roundoff for this environment. */

#line 145 "dsterf.f"
    eps = dlamch_("E", (ftnlen)1);
/* Computing 2nd power */
#line 146 "dsterf.f"
    d__1 = eps;
#line 146 "dsterf.f"
    eps2 = d__1 * d__1;
#line 147 "dsterf.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 148 "dsterf.f"
    safmax = 1. / safmin;
#line 149 "dsterf.f"
    ssfmax = sqrt(safmax) / 3.;
#line 150 "dsterf.f"
    ssfmin = sqrt(safmin) / eps2;
#line 151 "dsterf.f"
    rmax = dlamch_("O", (ftnlen)1);

/*     Compute the eigenvalues of the tridiagonal matrix. */

#line 155 "dsterf.f"
    nmaxit = *n * 30;
#line 156 "dsterf.f"
    sigma = 0.;
#line 157 "dsterf.f"
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

#line 163 "dsterf.f"
    l1 = 1;

#line 165 "dsterf.f"
L10:
#line 166 "dsterf.f"
    if (l1 > *n) {
#line 166 "dsterf.f"
	goto L170;
#line 166 "dsterf.f"
    }
#line 168 "dsterf.f"
    if (l1 > 1) {
#line 168 "dsterf.f"
	e[l1 - 1] = 0.;
#line 168 "dsterf.f"
    }
#line 170 "dsterf.f"
    i__1 = *n - 1;
#line 170 "dsterf.f"
    for (m = l1; m <= i__1; ++m) {
#line 171 "dsterf.f"
	if ((d__3 = e[m], abs(d__3)) <= sqrt((d__1 = d__[m], abs(d__1))) * 
		sqrt((d__2 = d__[m + 1], abs(d__2))) * eps) {
#line 173 "dsterf.f"
	    e[m] = 0.;
#line 174 "dsterf.f"
	    goto L30;
#line 175 "dsterf.f"
	}
#line 176 "dsterf.f"
/* L20: */
#line 176 "dsterf.f"
    }
#line 177 "dsterf.f"
    m = *n;

#line 179 "dsterf.f"
L30:
#line 180 "dsterf.f"
    l = l1;
#line 181 "dsterf.f"
    lsv = l;
#line 182 "dsterf.f"
    lend = m;
#line 183 "dsterf.f"
    lendsv = lend;
#line 184 "dsterf.f"
    l1 = m + 1;
#line 185 "dsterf.f"
    if (lend == l) {
#line 185 "dsterf.f"
	goto L10;
#line 185 "dsterf.f"
    }

/*     Scale submatrix in rows and columns L to LEND */

#line 190 "dsterf.f"
    i__1 = lend - l + 1;
#line 190 "dsterf.f"
    anorm = dlanst_("M", &i__1, &d__[l], &e[l], (ftnlen)1);
#line 191 "dsterf.f"
    iscale = 0;
#line 192 "dsterf.f"
    if (anorm == 0.) {
#line 192 "dsterf.f"
	goto L10;
#line 192 "dsterf.f"
    }
#line 194 "dsterf.f"
    if (anorm > ssfmax) {
#line 195 "dsterf.f"
	iscale = 1;
#line 196 "dsterf.f"
	i__1 = lend - l + 1;
#line 196 "dsterf.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 198 "dsterf.f"
	i__1 = lend - l;
#line 198 "dsterf.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 200 "dsterf.f"
    } else if (anorm < ssfmin) {
#line 201 "dsterf.f"
	iscale = 2;
#line 202 "dsterf.f"
	i__1 = lend - l + 1;
#line 202 "dsterf.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info, (ftnlen)1);
#line 204 "dsterf.f"
	i__1 = lend - l;
#line 204 "dsterf.f"
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info, (ftnlen)1);
#line 206 "dsterf.f"
    }

#line 208 "dsterf.f"
    i__1 = lend - 1;
#line 208 "dsterf.f"
    for (i__ = l; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 209 "dsterf.f"
	d__1 = e[i__];
#line 209 "dsterf.f"
	e[i__] = d__1 * d__1;
#line 210 "dsterf.f"
/* L40: */
#line 210 "dsterf.f"
    }

/*     Choose between QL and QR iteration */

#line 214 "dsterf.f"
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
#line 215 "dsterf.f"
	lend = lsv;
#line 216 "dsterf.f"
	l = lendsv;
#line 217 "dsterf.f"
    }

#line 219 "dsterf.f"
    if (lend >= l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

#line 225 "dsterf.f"
L50:
#line 226 "dsterf.f"
	if (l != lend) {
#line 227 "dsterf.f"
	    i__1 = lend - 1;
#line 227 "dsterf.f"
	    for (m = l; m <= i__1; ++m) {
#line 228 "dsterf.f"
		if ((d__2 = e[m], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
			+ 1], abs(d__1))) {
#line 228 "dsterf.f"
		    goto L70;
#line 228 "dsterf.f"
		}
#line 230 "dsterf.f"
/* L60: */
#line 230 "dsterf.f"
	    }
#line 231 "dsterf.f"
	}
#line 232 "dsterf.f"
	m = lend;

#line 234 "dsterf.f"
L70:
#line 235 "dsterf.f"
	if (m < lend) {
#line 235 "dsterf.f"
	    e[m] = 0.;
#line 235 "dsterf.f"
	}
#line 237 "dsterf.f"
	p = d__[l];
#line 238 "dsterf.f"
	if (m == l) {
#line 238 "dsterf.f"
	    goto L90;
#line 238 "dsterf.f"
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its */
/*        eigenvalues. */

#line 244 "dsterf.f"
	if (m == l + 1) {
#line 245 "dsterf.f"
	    rte = sqrt(e[l]);
#line 246 "dsterf.f"
	    dlae2_(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
#line 247 "dsterf.f"
	    d__[l] = rt1;
#line 248 "dsterf.f"
	    d__[l + 1] = rt2;
#line 249 "dsterf.f"
	    e[l] = 0.;
#line 250 "dsterf.f"
	    l += 2;
#line 251 "dsterf.f"
	    if (l <= lend) {
#line 251 "dsterf.f"
		goto L50;
#line 251 "dsterf.f"
	    }
#line 253 "dsterf.f"
	    goto L150;
#line 254 "dsterf.f"
	}

#line 256 "dsterf.f"
	if (jtot == nmaxit) {
#line 256 "dsterf.f"
	    goto L150;
#line 256 "dsterf.f"
	}
#line 258 "dsterf.f"
	++jtot;

/*        Form shift. */

#line 262 "dsterf.f"
	rte = sqrt(e[l]);
#line 263 "dsterf.f"
	sigma = (d__[l + 1] - p) / (rte * 2.);
#line 264 "dsterf.f"
	r__ = dlapy2_(&sigma, &c_b33);
#line 265 "dsterf.f"
	sigma = p - rte / (sigma + d_sign(&r__, &sigma));

#line 267 "dsterf.f"
	c__ = 1.;
#line 268 "dsterf.f"
	s = 0.;
#line 269 "dsterf.f"
	gamma = d__[m] - sigma;
#line 270 "dsterf.f"
	p = gamma * gamma;

/*        Inner loop */

#line 274 "dsterf.f"
	i__1 = l;
#line 274 "dsterf.f"
	for (i__ = m - 1; i__ >= i__1; --i__) {
#line 275 "dsterf.f"
	    bb = e[i__];
#line 276 "dsterf.f"
	    r__ = p + bb;
#line 277 "dsterf.f"
	    if (i__ != m - 1) {
#line 277 "dsterf.f"
		e[i__ + 1] = s * r__;
#line 277 "dsterf.f"
	    }
#line 279 "dsterf.f"
	    oldc = c__;
#line 280 "dsterf.f"
	    c__ = p / r__;
#line 281 "dsterf.f"
	    s = bb / r__;
#line 282 "dsterf.f"
	    oldgam = gamma;
#line 283 "dsterf.f"
	    alpha = d__[i__];
#line 284 "dsterf.f"
	    gamma = c__ * (alpha - sigma) - s * oldgam;
#line 285 "dsterf.f"
	    d__[i__ + 1] = oldgam + (alpha - gamma);
#line 286 "dsterf.f"
	    if (c__ != 0.) {
#line 287 "dsterf.f"
		p = gamma * gamma / c__;
#line 288 "dsterf.f"
	    } else {
#line 289 "dsterf.f"
		p = oldc * bb;
#line 290 "dsterf.f"
	    }
#line 291 "dsterf.f"
/* L80: */
#line 291 "dsterf.f"
	}

#line 293 "dsterf.f"
	e[l] = s * p;
#line 294 "dsterf.f"
	d__[l] = sigma + gamma;
#line 295 "dsterf.f"
	goto L50;

/*        Eigenvalue found. */

#line 299 "dsterf.f"
L90:
#line 300 "dsterf.f"
	d__[l] = p;

#line 302 "dsterf.f"
	++l;
#line 303 "dsterf.f"
	if (l <= lend) {
#line 303 "dsterf.f"
	    goto L50;
#line 303 "dsterf.f"
	}
#line 305 "dsterf.f"
	goto L150;

#line 307 "dsterf.f"
    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

#line 313 "dsterf.f"
L100:
#line 314 "dsterf.f"
	i__1 = lend + 1;
#line 314 "dsterf.f"
	for (m = l; m >= i__1; --m) {
#line 315 "dsterf.f"
	    if ((d__2 = e[m - 1], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
		    - 1], abs(d__1))) {
#line 315 "dsterf.f"
		goto L120;
#line 315 "dsterf.f"
	    }
#line 317 "dsterf.f"
/* L110: */
#line 317 "dsterf.f"
	}
#line 318 "dsterf.f"
	m = lend;

#line 320 "dsterf.f"
L120:
#line 321 "dsterf.f"
	if (m > lend) {
#line 321 "dsterf.f"
	    e[m - 1] = 0.;
#line 321 "dsterf.f"
	}
#line 323 "dsterf.f"
	p = d__[l];
#line 324 "dsterf.f"
	if (m == l) {
#line 324 "dsterf.f"
	    goto L140;
#line 324 "dsterf.f"
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its */
/*        eigenvalues. */

#line 330 "dsterf.f"
	if (m == l - 1) {
#line 331 "dsterf.f"
	    rte = sqrt(e[l - 1]);
#line 332 "dsterf.f"
	    dlae2_(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
#line 333 "dsterf.f"
	    d__[l] = rt1;
#line 334 "dsterf.f"
	    d__[l - 1] = rt2;
#line 335 "dsterf.f"
	    e[l - 1] = 0.;
#line 336 "dsterf.f"
	    l += -2;
#line 337 "dsterf.f"
	    if (l >= lend) {
#line 337 "dsterf.f"
		goto L100;
#line 337 "dsterf.f"
	    }
#line 339 "dsterf.f"
	    goto L150;
#line 340 "dsterf.f"
	}

#line 342 "dsterf.f"
	if (jtot == nmaxit) {
#line 342 "dsterf.f"
	    goto L150;
#line 342 "dsterf.f"
	}
#line 344 "dsterf.f"
	++jtot;

/*        Form shift. */

#line 348 "dsterf.f"
	rte = sqrt(e[l - 1]);
#line 349 "dsterf.f"
	sigma = (d__[l - 1] - p) / (rte * 2.);
#line 350 "dsterf.f"
	r__ = dlapy2_(&sigma, &c_b33);
#line 351 "dsterf.f"
	sigma = p - rte / (sigma + d_sign(&r__, &sigma));

#line 353 "dsterf.f"
	c__ = 1.;
#line 354 "dsterf.f"
	s = 0.;
#line 355 "dsterf.f"
	gamma = d__[m] - sigma;
#line 356 "dsterf.f"
	p = gamma * gamma;

/*        Inner loop */

#line 360 "dsterf.f"
	i__1 = l - 1;
#line 360 "dsterf.f"
	for (i__ = m; i__ <= i__1; ++i__) {
#line 361 "dsterf.f"
	    bb = e[i__];
#line 362 "dsterf.f"
	    r__ = p + bb;
#line 363 "dsterf.f"
	    if (i__ != m) {
#line 363 "dsterf.f"
		e[i__ - 1] = s * r__;
#line 363 "dsterf.f"
	    }
#line 365 "dsterf.f"
	    oldc = c__;
#line 366 "dsterf.f"
	    c__ = p / r__;
#line 367 "dsterf.f"
	    s = bb / r__;
#line 368 "dsterf.f"
	    oldgam = gamma;
#line 369 "dsterf.f"
	    alpha = d__[i__ + 1];
#line 370 "dsterf.f"
	    gamma = c__ * (alpha - sigma) - s * oldgam;
#line 371 "dsterf.f"
	    d__[i__] = oldgam + (alpha - gamma);
#line 372 "dsterf.f"
	    if (c__ != 0.) {
#line 373 "dsterf.f"
		p = gamma * gamma / c__;
#line 374 "dsterf.f"
	    } else {
#line 375 "dsterf.f"
		p = oldc * bb;
#line 376 "dsterf.f"
	    }
#line 377 "dsterf.f"
/* L130: */
#line 377 "dsterf.f"
	}

#line 379 "dsterf.f"
	e[l - 1] = s * p;
#line 380 "dsterf.f"
	d__[l] = sigma + gamma;
#line 381 "dsterf.f"
	goto L100;

/*        Eigenvalue found. */

#line 385 "dsterf.f"
L140:
#line 386 "dsterf.f"
	d__[l] = p;

#line 388 "dsterf.f"
	--l;
#line 389 "dsterf.f"
	if (l >= lend) {
#line 389 "dsterf.f"
	    goto L100;
#line 389 "dsterf.f"
	}
#line 391 "dsterf.f"
	goto L150;

#line 393 "dsterf.f"
    }

/*     Undo scaling if necessary */

#line 397 "dsterf.f"
L150:
#line 398 "dsterf.f"
    if (iscale == 1) {
#line 398 "dsterf.f"
	i__1 = lendsv - lsv + 1;
#line 398 "dsterf.f"
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 398 "dsterf.f"
    }
#line 401 "dsterf.f"
    if (iscale == 2) {
#line 401 "dsterf.f"
	i__1 = lendsv - lsv + 1;
#line 401 "dsterf.f"
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info, (ftnlen)1);
#line 401 "dsterf.f"
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

#line 408 "dsterf.f"
    if (jtot < nmaxit) {
#line 408 "dsterf.f"
	goto L10;
#line 408 "dsterf.f"
    }
#line 410 "dsterf.f"
    i__1 = *n - 1;
#line 410 "dsterf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 411 "dsterf.f"
	if (e[i__] != 0.) {
#line 411 "dsterf.f"
	    ++(*info);
#line 411 "dsterf.f"
	}
#line 413 "dsterf.f"
/* L160: */
#line 413 "dsterf.f"
    }
#line 414 "dsterf.f"
    goto L180;

/*     Sort eigenvalues in increasing order. */

#line 418 "dsterf.f"
L170:
#line 419 "dsterf.f"
    dlasrt_("I", n, &d__[1], info, (ftnlen)1);

#line 421 "dsterf.f"
L180:
#line 422 "dsterf.f"
    return 0;

/*     End of DSTERF */

} /* dsterf_ */

